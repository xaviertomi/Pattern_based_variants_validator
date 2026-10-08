#!/usr/bin/env python3
"""Owned input preparation, deterministic partitions and checksum provenance."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

# GT="alt" in bcftools 1.23.1 omits partial-missing ALT calls; retain every positive allele index.
ALT_FILTER_EXPR = 'GT="alt" || GT~"[1-9]"'


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False).encode('utf-8')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_bytes(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_bytes() == data:
        return
    fd, tmp = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def atomic_json(path, value):
    atomic_bytes(path, canonical(value) + b'\n')


def read_tsv(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent)
    os.close(fd)
    try:
        if str(path).endswith('.gz'):
            raw = open(tmp, 'wb')
            zipped = gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0)
            import io
            handle = io.TextIOWrapper(zipped, encoding='utf-8', newline='')
        else:
            handle = open(tmp, 'w', encoding='utf-8', newline='')
            raw = None
        with handle:
            writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)
        if raw:
            raw.close()
        if path.exists() and sha256(path) == sha256(tmp):
            os.unlink(tmp)
        else:
            os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def dataset_id(role, name):
    if role not in ('pacbio', 'illumina', 'other_illumina'):
        raise ValueError('Unknown dataset role: ' + role)
    return role + '_' + hashlib.sha256(name.encode('utf-8')).hexdigest()[:16]


def _integer(value, key, minimum=0, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise ValueError(f'{key} must be an integer in [{minimum},{maximum if maximum is not None else "Inf"}]')
    return value


def _positive(value, key, upper=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0 or (upper is not None and value > upper):
        raise ValueError(f'{key} must be positive' + (f' and <= {upper}' if upper else ''))
    return value


def _names(values, key):
    if not isinstance(values, list) or not values or any(not isinstance(x, str) or not x or '\n' in x or '\r' in x or '\t' in x for x in values) or len(set(values)) != len(values):
        raise ValueError(key + ' must contain nonempty unique exact sample names')
    return values


def validate_config(config, root, scope='selection'):
    c = json.loads(json.dumps(config))
    root = Path(root).resolve()
    reference_checkout = (root.parent / 'clone_vcf_filter').resolve()
    for key in ('reference_fasta', 'pacbio_vcf', 'illumina_vcf', 'output_dir', 'scratch_dir'):
        if not isinstance(c.get(key), str) or not c[key]:
            raise ValueError('Missing path setting: ' + key)
    for key in ('reference_fasta', 'pacbio_vcf', 'illumina_vcf', 'other_illumina', 'output_dir', 'scratch_dir'):
        if key in c:
            if not isinstance(c[key], str) or not c[key]:
                raise ValueError('Invalid path setting: ' + key)
            path = Path(c[key])
            if key == 'other_illumina':
                # Selection identity never stats even an unavailable external path.
                c[key] = os.path.abspath(root / path)
            else:
                c[key] = str((root / path).resolve() if not path.is_absolute() else path.resolve())
    execution = c.setdefault('execution', {})
    execution.setdefault('use_containers', True)
    execution.setdefault('scan_windows_per_shard', 50000)
    execution.setdefault('extra_bind_paths', [])
    execution.setdefault('image_cache', str(root / 'cache/images'))
    execution['image_cache'] = str((root / execution['image_cache']).resolve())
    writable = {k: c[k] for k in ('output_dir', 'scratch_dir')}
    writable.update(image_cache=execution['image_cache'], work=str(root / '.snakemake'), logs=str(root / 'logs'))
    for key, value in writable.items():
        resolved = Path(value).resolve()
        if resolved == reference_checkout or reference_checkout in resolved.parents:
            raise ValueError(f'Writable path inside read-only infrastructure reference: {key}: {resolved}')
    pb = _names(c.get('pacbio_samples'), 'pacbio_samples')
    mapping = c.get('illumina_corresponding')
    if not isinstance(mapping, dict) or set(mapping) != set(pb):
        raise ValueError('illumina_corresponding keys must equal pacbio_samples')
    il = []
    for name in pb:
        il.extend(_names(mapping[name], 'illumina_corresponding.' + name))
    if len(set(il)) != len(il):
        raise ValueError('Corresponding Illumina samples cannot be assigned twice')
    ids = [dataset_id(role, name) for role, names in [('pacbio', pb), ('illumina', il)] for name in names]
    if len(ids) != len(set(ids)):
        raise ValueError('Dataset identity collision')
    vcf_outputs = c.setdefault('vcf_outputs', {})
    if not isinstance(vcf_outputs, dict):
        raise ValueError('vcf_outputs must be a mapping')
    vcf_outputs.setdefault('consensus_fraction', 0.75)
    vcf_outputs.setdefault('groups', {})
    _positive(vcf_outputs['consensus_fraction'], 'vcf_outputs.consensus_fraction', 1)
    groups = vcf_outputs['groups']
    if not isinstance(groups, dict):
        raise ValueError('vcf_outputs.groups must be a mapping')
    illumina_names = set(il)
    for group_id, group in groups.items():
        if not isinstance(group_id, str) or not re.fullmatch(r'[a-z0-9][a-z0-9_-]*', group_id):
            raise ValueError(f'Unsafe VCF consensus group ID: {group_id!r}')
        if not isinstance(group, dict) or group.get('role') not in ('illumina', 'other_illumina'):
            raise ValueError(f'vcf_outputs.groups.{group_id}.role must be illumina or other_illumina')
        members = _names(group.get('samples'), f'vcf_outputs.groups.{group_id}.samples')
        if group['role'] == 'illumina' and not set(members) <= illumina_names:
            unknown = sorted(set(members) - illumina_names)
            raise ValueError(f'Unknown corresponding sample(s) in VCF group {group_id}: {unknown}')
    metadata = c.setdefault('sample_metadata', {})
    for role, names in [('pacbio', pb), ('illumina', il)]:
        values = metadata.setdefault(role, {})
        if not isinstance(values, dict) or set(values) - set(names) or any(not isinstance(v, dict) for v in values.values()):
            raise ValueError('Invalid sample_metadata.' + role)
    metadata.setdefault('other_illumina', {})
    if not isinstance(metadata['other_illumina'], dict) or any(not isinstance(value, dict) for value in metadata['other_illumina'].values()):
        raise ValueError('Invalid sample_metadata.other_illumina')
    for name in pb:
        clone = metadata['pacbio'].get(name, {}).get('clone_id')
        for other in mapping[name]:
            il_clone = metadata['illumina'].get(other, {}).get('clone_id')
            if clone is not None and il_clone is not None and clone != il_clone:
                raise ValueError('Clone metadata mismatch: ' + name + ' / ' + other)
    if 'fold_beds' in c or 'fold_beds' in c.get('folds', {}):
        raise ValueError('fold_beds is obsolete: folds are automatically generated from the reference index')
    g = c.setdefault('folds', {}).setdefault('grouping', {})
    g.setdefault('enabled', False)
    g.setdefault('regex', None)
    g.setdefault('unmatched_policy', 'error')
    if type(g['enabled']) is not bool:
        raise ValueError('folds.grouping.enabled must be boolean')
    if g['unmatched_policy'] not in ('error', 'exclude'):
        raise ValueError('Invalid folds.grouping.unmatched_policy')
    if g['enabled']:
        try:
            expression = re.compile(g['regex']) if isinstance(g['regex'], str) and g['regex'] else None
        except re.error as e:
            raise ValueError('Invalid grouping regex: ' + str(e)) from e
        if expression is None or expression.groups < 1:
            raise ValueError('Enabled grouping requires a valid regex with at least one capture')
    a = c.setdefault('analysis', {})
    defaults = dict(flank_bp=50, center_radius_bp=10, rotations=[1, 2, 3, 4, 5], shuffle_kmer=3, base_seed=42, seed_overrides=[], streme_min_width=6, streme_max_width=15, streme_total_length=10000000, streme_e_max=.05, fimo_p_threshold=.0001)
    for k, v in defaults.items():
        a.setdefault(k, v)
    flank = _integer(a['flank_bp'], 'analysis.flank_bp')
    _integer(a['center_radius_bp'], 'analysis.center_radius_bp', 0, flank)
    if a['rotations'] != [1, 2, 3, 4, 5] or any(type(x) is not int for x in a['rotations']):
        raise ValueError('analysis.rotations must be exactly [1,2,3,4,5]')
    _integer(a['shuffle_kmer'], 'analysis.shuffle_kmer', 1, 6)
    _integer(a['base_seed'], 'analysis.base_seed')
    low = _integer(a['streme_min_width'], 'analysis.streme_min_width', 1, 2 * flank + 1)
    _integer(a['streme_max_width'], 'analysis.streme_max_width', low, 2 * flank + 1)
    _integer(a['streme_total_length'], 'analysis.streme_total_length', 1)
    for k in ('streme_e_max', 'fimo_p_threshold'):
        _positive(a[k], 'analysis.' + k, 1)
    s = c.setdefault('selection', {})
    for k, v in dict(pacbio_coverage_ratio_min=1., illumina_coverage_ratio_min=1., require_all_pacbio_clones=True, corresponding_dataset_support_fraction=1., min_supported_rotations=3, require_central_enrichment=False).items():
        s.setdefault(k, v)
    for k in ('pacbio_coverage_ratio_min', 'illumina_coverage_ratio_min'):
        _positive(s[k], 'selection.' + k)
    _positive(s['corresponding_dataset_support_fraction'], 'selection.corresponding_dataset_support_fraction', 1)
    _integer(s['min_supported_rotations'], 'selection.min_supported_rotations', 1, 5)
    for obj, key in [(s, 'require_all_pacbio_clones'), (s, 'require_central_enrichment'), (execution, 'use_containers')]:
        if type(obj[key]) is not bool:
            raise ValueError(key + ' must be boolean')
    _integer(execution['scan_windows_per_shard'], 'execution.scan_windows_per_shard', 1)
    fm = s.setdefault('family_matching', {})
    fm.setdefault('method', 'tomtom')
    fm.setdefault('q_max', .05)
    fm.setdefault('min_overlap_bp', 6)
    if fm['method'] != 'tomtom':
        raise ValueError('family_matching.method must be tomtom')
    _positive(fm['q_max'], 'family_matching.q_max', 1)
    _integer(fm['min_overlap_bp'], 'family_matching.min_overlap_bp', 1, low)
    if 'pilot' in c:
        pilot = c['pilot']
        if not isinstance(pilot, dict) or pilot.get('pacbio_sample') not in pb:
            raise ValueError('pilot requires an exact selected PacBio sample')
        _integer(pilot.get('rotation'), 'pilot.rotation', 1, 5)
    if not isinstance(a['seed_overrides'], list):
        raise ValueError('analysis.seed_overrides must be a list')
    seed_rows(c)
    return c


def selection_config(c):
    return {k: c[k] for k in ('reference_fasta', 'pacbio_vcf', 'illumina_vcf', 'pacbio_samples', 'illumina_corresponding', 'folds', 'analysis', 'selection', 'container') if k in c} | {'sample_metadata': {k: c['sample_metadata'][k] for k in ('pacbio', 'illumina')}}


def configure_scratch(config, rule, identity='global', rotation=0, disk_mb=1000, scope='selection'):
    if scope == 'application':
        scientific = {'freeze_sha256': sha256(Path(config['output_dir']) / 'library/freeze.json'),
                      'external_metadata': config.get('sample_metadata', {}).get('other_illumina', {})}
    else:
        scientific = selection_config(config)
    signature = hashlib.sha256(canonical(scientific)).hexdigest()[:16]
    parent = Path(config['scratch_dir']) / signature / rule / identity / str(rotation)
    parent.mkdir(parents=True, exist_ok=True)
    required = _integer(disk_mb, 'disk_mb', 1) * 1024 * 1024
    if shutil.disk_usage(parent).free < required:
        raise ValueError(f'Insufficient scratch free space: {parent}; need {disk_mb} MiB; site quota must also permit this allocation')
    attempt = Path(tempfile.mkdtemp(prefix=os.environ.get('SLURM_JOB_ID', 'local') + '-' + str(os.getpid()) + '-', dir=parent))
    os.environ['TMPDIR'] = str(attempt)
    tempfile.tempdir = str(attempt)
    return attempt


def seed_rows(c):
    tuples = []
    for platform, names in [('pacbio', c['pacbio_samples']), ('illumina', [x for names in c['illumina_corresponding'].values() for x in names])]:
        for name in names:
            for r in range(1, 6):
                for role in (('train_shuffle', 'streme', 'test_shuffle') if platform == 'pacbio' else ('test_shuffle',)):
                    tuples.append((platform, name, r, role))
    overrides = {}
    for entry in c['analysis'].get('seed_overrides', []):
        if not isinstance(entry, dict):
            raise ValueError('Seed override must be an object')
        key = (entry.get('platform'), entry.get('sample'), entry.get('rotation'), entry.get('role'))
        if key not in tuples or key in overrides or type(entry.get('rotation')) is not int:
            raise ValueError('Unknown or duplicate seed override tuple: ' + repr(key))
        overrides[key] = _integer(entry.get('seed'), 'seed override', 1, 2147483646)
    by_seed = {}
    for key, seed in overrides.items():
        by_seed.setdefault(seed, []).append(key)
    for seed, keys in by_seed.items():
        if len(keys) > 1 and not (len(keys) == 2 and keys[0][:3] == keys[1][:3] and keys[0][0] == 'pacbio' and {k[3] for k in keys} == {'train_shuffle', 'streme'}):
            raise ValueError('Conflicting explicit seeds: ' + str(seed))
    used = set(by_seed)
    rows = []
    for key in sorted(tuples, key=lambda item: (item[3], item[0], item[1], item[2])):
        platform, name, r, role = key
        seed = overrides.get(key)
        if seed is None:
            seed = int.from_bytes(hashlib.sha256(canonical([c['analysis']['base_seed'], role, platform, name, r])).digest()[:4], 'big') % 2147483646 + 1
            while seed in used:
                seed = seed % 2147483646 + 1
        used.add(seed)
        rows.append(dict(platform=platform, source_role=platform, dataset_id=dataset_id(platform, name), sample_name=name, rotation=r, role=role, seed=seed, origin='overridden' if key in overrides else 'generated'))
    return rows


def parse_gt(gt, alt_count):
    if not isinstance(gt, str) or not re.fullmatch(r'(?:\.|[0-9]+)(?:[/|](?:\.|[0-9]+))*', gt):
        raise ValueError('Malformed GT: ' + repr(gt))
    tokens = re.split(r'[/|]', gt)
    values = [None if x == '.' else int(x) for x in tokens]
    if any(x is not None and x > alt_count for x in values):
        raise ValueError(f'Out-of-range GT {gt}; ALT count {alt_count}')
    return {'alleles': values, 'ploidy': len(values), 'partial_missing': any(x is None for x in values) and any(x is not None for x in values), 'gt_status': 'alt' if any(x is not None and x > 0 for x in values) else ('missing' if all(x is None for x in values) else 'reference_only')}


def read_fai(path):
    rows = {}
    with open(path, encoding='utf-8') as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip('\n\r').split('\t')
            name = fields[0]
            if len(fields) < 5 or not name or name in rows or not re.fullmatch('[0-9]+', fields[1]) or int(fields[1]) <= 0:
                raise ValueError(f'Invalid reference index {path}: line {line_number}, contig {name!r} (malformed, duplicate or nonpositive length)')
            if any(not re.fullmatch('[0-9]+', x) for x in fields[2:5]) or int(fields[3]) <= 0 or int(fields[4]) < int(fields[3]):
                raise ValueError(f'Invalid reference index {path}: line {line_number}, contig {name!r}')
            rows[name] = int(fields[1])
    if not rows:
        raise ValueError('Empty reference index: ' + str(path))
    return rows


def _assign(lengths, grouping):
    enabled = grouping.get('enabled', False)
    expression = re.compile(grouping['regex']) if enabled else None
    groups, excluded, assignments = {}, [], []
    for contig, length in sorted(lengths.items()):
        match = expression.match(contig) if enabled else None
        if enabled and not match:
            excluded.append(dict(contig=contig, length_bp=length, group_id='', fold_id='', inclusion_status='excluded', reason='unmatched_group_regex'))
            continue
        group = match.group(1) if enabled else contig
        if not group:
            raise ValueError('Empty grouping capture: ' + contig)
        groups.setdefault(group, []).append((contig, length))
    totals = [0] * 5
    # ponytail: indivisible groups and greedy assignment cannot guarantee equal loads; optimize only if measured imbalance makes folds unusable.
    for group in sorted(groups, key=lambda k: (-sum(n for _, n in groups[k]), k)):
        fold = min(range(5), key=lambda r: (totals[r], r))
        for contig, length in groups[group]:
            assignments.append(dict(contig=contig, length_bp=length, group_id=group, fold_id=fold + 1, inclusion_status='included', reason=''))
        totals[fold] += sum(n for _, n in groups[group])
    return sorted(assignments + excluded, key=lambda row: row['contig']), groups, excluded, totals


def build_folds(reference_fai, reference_sha256, output_dir, grouping):
    output_dir = Path(output_dir)
    grouping = {'enabled': False, 'regex': None, 'unmatched_policy': 'error'} | grouping
    lengths = read_fai(reference_fai)
    assignments, groups, excluded, totals = _assign(lengths, grouping)
    report = output_dir / 'reports/unmatched_contigs.tsv'
    write_tsv(report, [dict(contig=x['contig'], length=x['length_bp'], policy=grouping['unmatched_policy'], action='error' if grouping['unmatched_policy'] == 'error' else 'exclude', reason=x['reason']) for x in excluded], ['contig', 'length', 'policy', 'action', 'reason'])
    errors = []
    if excluded and grouping['unmatched_policy'] == 'error':
        errors.append(f'Unmatched contigs: regex={grouping["regex"]!r}; count={len(excluded)}; names={",".join(x["contig"] for x in excluded)}; report={report}')
    if len(groups) < 5:
        errors.append(f'Cannot generate five folds: only {len(groups)} included groups are available; need at least 5 (grouping_enabled={grouping["enabled"]}, unmatched_policy={grouping["unmatched_policy"]}, excluded_contigs={len(excluded)}).')
    if errors:
        raise ValueError('\n'.join(errors))
    fold_rows = []
    for fold in range(1, 6):
        members = [x for x in assignments if x['fold_id'] == fold]
        path = output_dir / f'folds/fold{fold}.bed'
        atomic_bytes(path, ''.join(f'{x["contig"]}\t0\t{x["length_bp"]}\n' for x in members).encode('utf-8'))
        fold_rows.append(dict(fold_id=fold, bed=f'folds/fold{fold}.bed', sha256=sha256(path), total_bp=totals[fold - 1], group_count=len({x['group_id'] for x in members}), contig_count=len(members)))
    write_tsv(output_dir / 'manifests/fold_assignments.tsv', assignments, ['contig', 'length_bp', 'group_id', 'fold_id', 'inclusion_status', 'reason'])
    write_tsv(output_dir / 'reports/fold_balance.tsv', [x | {'largest_bp': max(totals), 'smallest_bp': min(totals), 'largest_smallest_ratio': max(totals) / min(totals)} for x in fold_rows], ['fold_id', 'total_bp', 'group_count', 'contig_count', 'largest_bp', 'smallest_bp', 'largest_smallest_ratio'])
    manifest = dict(schema_version=1, algorithm='lpt_group_bp_v1', reference_sha256=reference_sha256, index_sha256=sha256(reference_fai), grouping=grouping, assignments=assignments, exclusions=excluded, folds=fold_rows, assignment_sha256=hashlib.sha256(canonical(assignments)).hexdigest())
    atomic_json(output_dir / 'manifests/folds.json', manifest)
    verify_folds(reference_fai, output_dir, reference_sha256, grouping)
    return manifest


def verify_folds(fai, output_dir, reference_sha=None, grouping=None):
    output_dir = Path(output_dir)
    manifest = json.loads((output_dir / 'manifests/folds.json').read_text(encoding='utf-8'))
    lengths = read_fai(fai)
    g = grouping if grouping is not None else manifest['grouping']
    expected, groups, excluded, totals = _assign(lengths, g)
    if len(groups) < 5 or (excluded and g['unmatched_policy'] == 'error') or manifest['algorithm'] != 'lpt_group_bp_v1' or manifest['grouping'] != g or manifest['assignments'] != expected or manifest['assignment_sha256'] != hashlib.sha256(canonical(expected)).hexdigest() or manifest['index_sha256'] != sha256(fai) or (reference_sha and manifest['reference_sha256'] != reference_sha):
        raise ValueError('Invalid automatic fold manifest: ' + str(output_dir / 'manifests/folds.json'))
    assignment_rows = list(read_tsv(output_dir / 'manifests/fold_assignments.tsv'))
    if assignment_rows != [{key: str(value) for key, value in row.items()} for row in expected] or manifest['exclusions'] != excluded:
        raise ValueError('Invalid fold assignment table or exclusions')
    if len(manifest['folds']) != 5:
        raise ValueError('Invalid five-fold manifest')
    for r, entry in enumerate(manifest['folds'], 1):
        path = output_dir / f'folds/fold{r}.bed'
        expected_bytes = ''.join(f'{x["contig"]}\t0\t{x["length_bp"]}\n' for x in expected if x['fold_id'] == r).encode('utf-8')
        if path.read_bytes() != expected_bytes or entry['fold_id'] != r or entry['sha256'] != sha256(path) or entry['total_bp'] != totals[r - 1]:
            raise ValueError('Invalid automatic fold partition: ' + str(path))
    return manifest


def input_signature(input_paths, signature_params):
    return hashlib.sha256(canonical({'inputs': {str(p): sha256(p) for p in sorted(map(str, input_paths))}, 'params': signature_params})).hexdigest()


def fingerprint_inputs(paths, params, destination):
    result = {'schema_version': 1, 'inputs': {key: {'path': str(path), 'sha256': sha256(path)} for key, path in sorted(paths.items())}, 'params': params}
    result['signature'] = hashlib.sha256(canonical(result)).hexdigest()
    atomic_json(destination, result)
    return result


def publish_provenance(primary, inputs, params, outputs, scope='selection', producer_target=None):
    params = dict(params)
    code_paths = sorted(set(map(str, params.get('code_paths', []))) | {str(Path(__file__).resolve())})
    params['code_paths'] = code_paths
    params['code_hashes'] = {path: sha256(path) for path in code_paths}
    paths = sorted(set(map(str, inputs)))
    record = dict(schema_version=1, scope=scope, producer_target=producer_target or str(primary), input_paths=paths, signature_params=params, input_signature=input_signature(paths, params), required_outputs={str(p): sha256(p) for p in outputs})
    atomic_json(str(primary) + '.provenance.json', record)
    return record


def check_reuse(provenance_path, expected_signature, expected_outputs):
    try:
        record = json.loads(Path(provenance_path).read_text(encoding='utf-8'))
        expected = set(map(str, expected_outputs))
        return record['schema_version'] == 1 and record['scope'] in ('selection', 'application') and isinstance(record['producer_target'], str) and record['input_signature'] == expected_signature and set(record['required_outputs']) == expected and all(sha256(p) == record['required_outputs'][p] for p in expected)
    except (OSError, ValueError, KeyError, TypeError):
        return False


def run_command(args, stdout=None, log=None):
    if log is not None:
        log = Path(log)
        log.parent.mkdir(parents=True, exist_ok=True)
        err = open(log, 'ab')
    else:
        err = subprocess.PIPE
    try:
        result = subprocess.run(list(map(str, args)), stdout=stdout if stdout is not None else subprocess.PIPE, stderr=err, check=False)
        if result.returncode:
            detail = result.stderr.decode('utf-8', errors='replace') if result.stderr else f'see {log}'
            raise RuntimeError(f'Command failed ({result.returncode}): {args!r}: {detail}')
        return result.stdout.decode('utf-8') if result.stdout is not None else ''
    finally:
        if log is not None:
            err.close()


def stage_reference(source, output_dir):
    source = Path(source)
    target = Path(output_dir) / 'reference/reference.fa'
    target.parent.mkdir(parents=True, exist_ok=True)
    if not source.is_file():
        raise ValueError('Unreadable reference: ' + str(source))
    with open(source, 'rb') as f:
        compressed = f.read(2) == b'\x1f\x8b'
    if compressed:
        with gzip.open(source, 'rb') as src, open(str(target) + '.tmp', 'wb') as dst:
            shutil.copyfileobj(src, dst)
        os.replace(str(target) + '.tmp', target)
    else:
        # Copying keeps indexes owned and freeze portable without symlink-target mount assumptions.
        shutil.copyfile(source, str(target) + '.tmp')
        os.replace(str(target) + '.tmp', target)
    return target


def index_reference(reference):
    names = set()
    with open(reference, encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('>'):
                tokens = line[1:].split()
                if not tokens or tokens[0] in names:
                    raise ValueError('Duplicate/empty reference contig: ' + line.strip())
                names.add(tokens[0])
    run_command(['samtools', 'faidx', reference])
    read_fai(str(reference) + '.fai')
    return Path(str(reference) + '.fai')


def stage_vcf(source, output_dir, role):
    target = Path(output_dir) / f'inputs/{role}.vcf.gz'
    target.parent.mkdir(parents=True, exist_ok=True)
    if not Path(source).is_file():
        raise ValueError('Unreadable VCF: ' + str(source))
    tmp = Path(str(target) + '.tmp.gz')
    run_command(['bcftools', 'view', '-Oz', '-o', tmp, source], log=Path(output_dir) / f'logs/stage_vcf/{role}.stderr.log')
    run_command(['bcftools', 'index', '--csi', '--force', tmp], log=Path(output_dir) / f'logs/stage_vcf/{role}.stderr.log')
    os.replace(tmp, target)
    os.replace(str(tmp) + '.csi', str(target) + '.csi')
    return target


def header_samples(vcf):
    names = run_command(['bcftools', 'query', '-l', vcf]).splitlines()
    if len(names) != len(set(names)):
        raise ValueError('Duplicate VCF header samples: ' + str(vcf))
    return names


def validate_header(vcf, reference, selected):
    return validate_header_text(vcf, reference, selected, header_samples(vcf), run_command(['bcftools', 'view', '-h', vcf]))


def validate_header_text(vcf, reference, selected, names, text):
    if len(names) != len(set(names)):
        raise ValueError('Duplicate VCF header samples: ' + str(vcf))
    missing = [name for name in selected if name not in names]
    if missing:
        raise ValueError(f'Missing samples in {vcf}: {missing}; available: {names}')
    definitions = [line for line in text.splitlines() if re.match(r'##FORMAT=<ID=GT(?:,|>)', line)]
    if len(definitions) != 1 or not re.search(r'(?:<|,)Number=1(?:,|>)', definitions[0]) or not re.search(r'(?:<|,)Type=String(?:,|>)', definitions[0]):
        raise ValueError('Missing/invalid FORMAT/GT definition: ' + str(vcf))
    lengths = read_fai(str(reference) + '.fai')
    for line in text.splitlines():
        if line.startswith('##contig=<'):
            match = re.search(r'(?:<|,)ID=([^,>]+)', line)
            lm = re.search(r'(?:<|,)length=([0-9]+)', line, re.IGNORECASE)
            name = match.group(1).strip('"') if match else ''
            if name not in lengths or (lm and int(lm.group(1)) != lengths[name]):
                raise ValueError('VCF/reference contig mismatch: ' + str(vcf) + ': ' + name)
    return names


def variant_type(ref, alt):
    alleles = alt.split(',')
    if any(not re.fullmatch('[ACGTNacgtn]+', x) for x in alleles):
        return 'symbolic'
    kinds = {'SNP' if len(ref) == len(x) == 1 else ('MNV' if len(ref) == len(x) else 'indel') for x in alleles}
    return next(iter(kinds)) if len(kinds) == 1 else 'mixed'


WINDOW_FIELDS = ['window_id', 'CHROM', 'POS', 'interval_start', 'interval_end', 'fold_id', 'fold_status', 'extraction_status', 'extraction_reason', 'ambiguity_count']
GENOTYPE_FIELDS = ['source_vcf', 'source_record_index', 'record_id', 'dataset_id', 'sample_name', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'GT', 'variant_type', 'gt_status', 'partial_missing', 'ploidy'] + [field for field in WINDOW_FIELDS if field not in ('CHROM', 'POS')]


def _query_rows(vcf):
    process = subprocess.Popen(['bcftools', 'query', '-f', '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n', str(vcf)], stdout=subprocess.PIPE, stderr=subprocess.TemporaryFile(), text=True)
    try:
        for line in process.stdout:
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 6:
                raise ValueError('Malformed single-sample query row: ' + line.strip())
            yield dict(zip(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'GT'], fields))
        if process.wait():
            process.stderr.seek(0)
            raise RuntimeError('bcftools query failed: ' + process.stderr.read())
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()
        process.stdout.close()
        process.stderr.close()

def query_file_rows(path):
    with open(path, encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            values = line.rstrip('\n').split('\t')
            if len(values) != 6:
                raise ValueError(f'Malformed query stream {path}: line {number}')
            yield dict(zip(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'GT'], values))



def prepare_sample(vcf, reference, source_vcf, source_hash, role, name, output_dir, flank, fold_manifest=None, phase='all', header_validated=False):
    output_dir = Path(output_dir)
    did = dataset_id(role, name)
    directory = output_dir / 'prepared' / did
    directory.mkdir(parents=True, exist_ok=True)
    if phase in ('all', 'records', 'records-from-query'):
        if not header_validated:
            validate_header(vcf, reference, [name])
        sample_file = directory / 'exact_sample.txt'
        atomic_bytes(sample_file, (name + '\n').encode('utf-8'))
        sample_vcf = directory / 'sample.vcf.gz'
        alt_vcf = directory / 'sample.alt.vcf.gz'
        log = output_dir / 'logs/prepare_sample' / did / 'records.stderr.log'
        if phase != 'records-from-query':
            run_command(['bcftools', 'view', '-S', sample_file, '-Oz', '-o', sample_vcf, vcf], log=log)
            run_command(['bcftools', 'index', '--csi', '--force', sample_vcf], log=log)
        lengths = read_fai(str(reference) + '.fai')
        def genotypes():
            stream = query_file_rows(directory / 'raw.records.tsv') if phase == 'records-from-query' else _query_rows(sample_vcf)
            for index, row in enumerate(stream, 1):
                try:
                    gt = parse_gt(row['GT'], 0 if row['ALT'] == '.' else len(row['ALT'].split(',')))
                    pos = int(row['POS'])
                    if row['CHROM'] not in lengths or pos < 1 or pos > lengths[row['CHROM']]:
                        raise ValueError('Unknown contig or out-of-range POS')
                except ValueError as e:
                    raise ValueError(f'{source_vcf}: sample {name}, record {index}: {e}') from e
                yield row | {'source_vcf': str(source_vcf), 'source_record_index': index, 'record_id': f'{role}_{source_hash}_{index}', 'dataset_id': did, 'sample_name': name, 'variant_type': variant_type(row['REF'], row['ALT']), 'gt_status': gt['gt_status'], 'partial_missing': int(gt['partial_missing']), 'ploidy': gt['ploidy']}
        write_tsv(directory / 'records.tsv.gz', genotypes(), GENOTYPE_FIELDS)
        if phase != 'records-from-query':
            run_command(['bcftools', 'view', '-i', ALT_FILTER_EXPR, '-Oz', '-o', alt_vcf, sample_vcf], log=log)
            run_command(['bcftools', 'index', '--csi', '--force', alt_vcf], log=log)
            with open(os.devnull, 'wb') as discarded:
                run_command(['bcftools', 'norm', '-f', reference, '-c', 'e', '-Ou', alt_vcf], stdout=discarded, log=log)
        original_alt = (row for row in read_tsv(directory / 'records.tsv.gz') if row['gt_status'] == 'alt')
        stream = query_file_rows(directory / 'raw.alt.tsv') if phase == 'records-from-query' else _query_rows(alt_vcf)
        for filtered in stream:
            original = next(original_alt, None)
            if original is None or any(filtered[k] != original[k] for k in ('CHROM', 'POS', 'REF', 'ALT', 'GT')):
                raise ValueError('ALT occurrence merge disagrees with validated GT stream: ' + name)
        if next(original_alt, None) is not None:
            raise ValueError('ALT occurrence merge omitted original records: ' + name)
        atomic_json(directory / 'records.complete.json', {'dataset_id': did, 'source_sha256': source_hash, 'sample_name': name})
        if phase in ('records', 'records-from-query'):
            return {'dataset_id': did, 'directory': str(directory)}
    lengths = read_fai(str(reference) + '.fai')
    ref_sha = sha256(reference)
    if isinstance(fold_manifest, (str, Path)):
        fold_manifest = json.loads(Path(fold_manifest).read_text(encoding='utf-8'))
    assignments = {row['contig']: row for row in fold_manifest['assignments']} if fold_manifest else {}
    import sqlite3
    fd, temporary_database = tempfile.mkstemp(prefix=did + '.windows.', suffix='.sqlite', dir=os.environ.get('TMPDIR'))
    os.close(fd)
    db_path = Path(temporary_database)
    db = sqlite3.connect(db_path)
    db.execute('PRAGMA cache_size=-65536')
    db.execute('PRAGMA temp_store=FILE')
    db.execute('CREATE TABLE windows (id TEXT PRIMARY KEY, chrom TEXT, pos INTEGER, start INTEGER, stop INTEGER, fold TEXT, fold_status TEXT, status TEXT, reason TEXT, sequence TEXT, ambiguity INTEGER)')
    def attach_windows():
        for row in read_tsv(directory / 'records.tsv.gz'):
            if row['gt_status'] != 'alt':
                yield row
                continue
            pos = int(row['POS'])
            start, stop = pos - flank, pos + flank
            wid = 'w_' + hashlib.sha256(canonical([ref_sha, row['CHROM'], start, stop])).hexdigest()[:16]
            assignment = assignments.get(row['CHROM'])
            fold = assignment['fold_id'] if assignment else ''
            fold_status = assignment['inclusion_status'] if assignment else 'application'
            status, reason = ('not_scanned', 'boundary_window') if start < 1 or stop > lengths[row['CHROM']] else ('eligible', '')
            existing = db.execute('SELECT chrom,start,stop FROM windows WHERE id=?', (wid,)).fetchone()
            if existing and existing != (row['CHROM'], start, stop):
                raise ValueError('Window identity collision: ' + wid)
            db.execute('INSERT OR IGNORE INTO windows VALUES (?,?,?,?,?,?,?,?,?,?,?)', (wid, row['CHROM'], pos, start, stop, str(fold), fold_status, status, reason, '', 0))
            yield row | dict(window_id=wid, interval_start=start, interval_end=stop, fold_id=fold, fold_status=fold_status, extraction_status=status, extraction_reason=reason)
        db.commit()
    write_tsv(directory / 'genotypes.tsv.gz', attach_windows(), GENOTYPE_FIELDS)
    regions = directory / 'regions.txt'
    if phase != 'contexts-from-fasta':
        with open(regions, 'w', encoding='utf-8', newline='\n') as handle:
            for chrom, start, stop in db.execute("SELECT chrom,start,stop FROM windows WHERE status='eligible' ORDER BY id"):
                handle.write(f'{chrom}:{start}-{stop}\n')
    if phase == 'contexts-plan':
        db.close()
        db_path.unlink()
        atomic_json(directory / 'contexts.plan.json', {'dataset_id': did, 'regions_sha256': sha256(regions)})
        return {'dataset_id': did, 'directory': str(directory)}
    fasta = directory / 'extracted.fasta'
    if regions.stat().st_size:
        if phase != 'contexts-from-fasta':
            with open(fasta, 'wb') as handle:
                run_command(['samtools', 'faidx', '-r', regions, reference], stdout=handle, log=output_dir / 'logs/prepare_sample' / did / 'contexts.stderr.log')
        expected = iter(db.execute("SELECT id,chrom,start,stop FROM windows WHERE status='eligible' ORDER BY id"))
        def save_sequence(header, sequence):
            row = next(expected, None)
            if row is None or header != f'{row[1]}:{row[2]}-{row[3]}' or len(sequence) != 2 * flank + 1 or not re.fullmatch('[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+', sequence):
                raise ValueError('Failed/invalid reference extraction: ' + str(header))
            db.execute('UPDATE windows SET sequence=?,ambiguity=? WHERE id=?', (sequence.upper(), sum(ch.upper() not in 'ACGT' for ch in sequence), row[0]))
        header, parts = None, []
        with open(fasta, encoding='utf-8') as handle:
            for line in handle:
                if line.startswith('>'):
                    if header is not None:
                        save_sequence(header, ''.join(parts))
                    header, parts = line[1:].strip(), []
                else:
                    parts.append(line.strip())
            if header is not None:
                save_sequence(header, ''.join(parts))
        if next(expected, None) is not None:
            raise ValueError('Reference extraction omitted windows: ' + did)
    db.commit()
    def add_ambiguity_counts():
        for row in read_tsv(directory / 'genotypes.tsv.gz'):
            if row['window_id']:
                row['ambiguity_count'] = db.execute('SELECT ambiguity FROM windows WHERE id=?', (row['window_id'],)).fetchone()[0]
            yield row
    write_tsv(directory / 'genotypes.tsv.gz', add_ambiguity_counts(), GENOTYPE_FIELDS)
    def windows(fold=None):
        query = 'SELECT id,chrom,pos,start,stop,fold,fold_status,status,reason,ambiguity FROM windows'
        args = ()
        if fold is not None:
            query += ' WHERE fold=?'
            args = (str(fold),)
        for values in db.execute(query + ' ORDER BY id', args):
            yield dict(zip(WINDOW_FIELDS, values))
    write_tsv(directory / 'windows.tsv', windows(), WINDOW_FIELDS)
    counts = []
    for fold in [None] + list(range(1, 6)):
        prefix = 'all' if fold is None else f'fold{fold}'
        path = directory / (prefix + '.fasta')
        with open(str(path) + '.tmp', 'w', encoding='utf-8', newline='\n') as handle:
            query = "SELECT id,sequence FROM windows WHERE status='eligible'" + ('' if fold is None else ' AND fold=?') + ' ORDER BY id'
            count = 0
            for wid, seq in db.execute(query, () if fold is None else (str(fold),)):
                handle.write(f'>{wid}\n{seq}\n')
                count += 1
        os.replace(str(path) + '.tmp', path)
        if fold is not None:
            write_tsv(directory / (prefix + '.windows.tsv'), windows(fold), WINDOW_FIELDS)
            counts.append(dict(dataset_id=did, sample_name=name, role=role, fold_id=fold, eligible_windows=count))
    write_tsv(directory / 'fold_counts.tsv', counts, ['dataset_id', 'sample_name', 'role', 'fold_id', 'eligible_windows'])
    db.close()
    db_path.unlink()
    result = dict(dataset_id=did, sample_name=name, directory=str(directory), reference_sha256=ref_sha)
    atomic_json(directory / 'sample.complete.json', result)
    return result


def datasets(c):
    result = {}
    for role, names in [('pacbio', c['pacbio_samples']), ('illumina', [x for names in c['illumina_corresponding'].values() for x in names])]:
        for name in names:
            result[dataset_id(role, name)] = dict(dataset_id=dataset_id(role, name), role=role, sample_name=name, source_vcf=c[role + '_vcf'], corresponding_pacbio=name if role == 'pacbio' else next(pb for pb, il in c['illumina_corresponding'].items() if name in il), metadata=c['sample_metadata'][role].get(name, {}))
    return result


def preflight(tool, directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    expected = {'bcftools': '1.23.1', 'samtools': '1.23.1', 'meme': '5.5.7', 'python': '3.12.10', 'r': '4.5.1'}
    if tool == 'meme':
        from evidence import preflight_meme
        result = preflight_meme(directory)
        result['expected_version'] = expected[tool]
        return result
    commands = {'bcftools': ['bcftools', '--version'], 'samtools': ['samtools', '--version'], 'python': ['python3', '--version'], 'r': ['Rscript', '--version']}
    response = run_command(commands[tool])
    if expected[tool] not in response.splitlines()[0]:
        raise ValueError(f'Pinned {tool} version mismatch: expected {expected[tool]}, got {response}')
    atomic_bytes(directory / 'version.txt', response.encode('utf-8'))
    if tool == 'samtools':
        fasta = directory / 'tiny.fa'
        atomic_bytes(fasta, b'>tiny\nACGTACGTACGT\n')
        index_reference(fasta)
        observed = run_command(['samtools', 'faidx', fasta, 'tiny:2-5'])
        if ''.join(observed.splitlines()[1:]) != 'CGTA':
            raise ValueError('samtools faidx coordinate preflight failed')
    if tool == 'bcftools':
        fasta = directory / 'tiny.fa'
        atomic_bytes(fasta, b'>tiny\nACGTACGTACGT\n')
        atomic_bytes(str(fasta) + '.fai', b'tiny\t12\t6\t12\t13\n')
        vcf = directory / 'tiny.vcf'
        atomic_bytes(vcf, b'##fileformat=VCFv4.2\n##contig=<ID=tiny,length=12>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\ntiny\t2\t.\tC\tT\t.\tPASS\t.\tGT\t1/.\n')
        staged = stage_vcf(vcf, directory, 'pacbio')
        if header_samples(staged) != ['S'] or list(_query_rows(staged))[0]['GT'] != '1/.':
            raise ValueError('bcftools sample/query preflight failed')
        with open(os.devnull, 'wb') as handle:
            run_command(['bcftools', 'norm', '-f', fasta, '-c', 'e', '-Ou', staged], stdout=handle)
    if tool == 'r':
        result = run_command(['Rscript', '-e', 'cat(sum(c(1,2,3)))'])
        if result.strip() != '6':
            raise ValueError('R worker execution failed')
    probe = directory / '.write_probe'
    atomic_bytes(probe, b'worker write check\n')
    probe.unlink()
    return dict(tool=tool, expected_version=expected[tool], version=response.strip(), commands=commands[tool])


def reusable_action(primary, params, expected_outputs):
    try:
        record = json.loads(Path(str(primary) + '.provenance.json').read_text(encoding='utf-8'))
        current = dict(record['signature_params'])
        current.update(params)
        current['code_hashes'] = {path: sha256(path) for path in current['code_paths']}
        signature = input_signature(record['input_paths'], current)
        return check_reuse(str(primary) + '.provenance.json', signature, expected_outputs)
    except (OSError, ValueError, KeyError, TypeError):
        return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('action')
    parser.add_argument('--config', required=True)
    parser.add_argument('--root', default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument('--role')
    parser.add_argument('--dataset-id')
    parser.add_argument('--phase', default='all')
    parser.add_argument('--tool')
    parser.add_argument('--disk-mb', type=int, default=1000)
    args = parser.parse_args()
    c = validate_config(json.loads(Path(args.config).read_text(encoding='utf-8')), args.root)
    if args.action != 'fingerprint-selection':
        configure_scratch(c, args.action, args.dataset_id or args.tool or args.role or 'global', disk_mb=args.disk_mb)
    out = Path(c['output_dir'])
    reference = out / 'reference/reference.fa'
    ds = datasets(c)
    action = args.action
    relevant = {
        'stage-reference': ['reference_fasta', 'container.python'],
        'index-reference': ['container.samtools'],
        'stage-vcf': [f'{args.role}_vcf', 'container.bcftools'],
        'headers': ['pacbio_samples', 'illumina_corresponding', 'container.bcftools'],
        'folds': ['folds.grouping', 'container.python'],
        'validate-folds': ['folds.grouping', 'container.python'],
        'seeds': ['analysis.base_seed', 'analysis.seed_overrides', 'pacbio_samples', 'illumina_corresponding', 'container.python'],
        'prepare-sample': ['analysis.flank_bp', 'pacbio_samples', 'illumina_corresponding', 'sample_metadata.pacbio', 'sample_metadata.illumina', 'container.python'],
        'sample-name': ['pacbio_samples', 'illumina_corresponding', 'container.python'],
        'capture-sample-tools': ['container.bcftools', 'container.python'],
        'capture-extraction': ['container.samtools', 'container.python'],
        'preflight': ['container.' + str(args.tool)],
        'preflight-gate': ['container.python'],
        'gate': ['pacbio_samples', 'illumina_corresponding', 'container.python'],
    }.get(action, [])
    def config_value(key):
        value = c
        for part in key.split('.'):
            value = value.get(part) if isinstance(value, dict) else None
        return value
    inputs, outputs, params = [], [], {'action': action, 'config_paths': relevant, 'config_values': {key: config_value(key) for key in relevant}}
    if action == 'fingerprint-selection':
        scientific = selection_config(c)
        atomic_json(out / 'manifests/selection_config.json', scientific)
        fingerprint_inputs({k: Path(c[k]) for k in ('reference_fasta', 'pacbio_vcf', 'illumina_vcf')} | {'prepare_code': Path(__file__)}, scientific, out / 'manifests/selection_fingerprints.json')
        write_tsv(out / 'manifests/datasets.tsv', [row | {'metadata': json.dumps(row['metadata'], sort_keys=True)} for row in ds.values()], ['dataset_id', 'role', 'sample_name', 'source_vcf', 'corresponding_pacbio', 'metadata'])
        return
    if action == 'preflight-fixtures':
        outputs = [out / 'preflight/fixtures.complete.json']
        for tool in ('bcftools', 'samtools', 'r'):
            directory = out / 'preflight' / tool
            atomic_bytes(directory / 'tiny.fa', b'>tiny\nACGTACGTACGT\n')
            atomic_bytes(directory / 'tiny.fa.fai', b'tiny\t12\t6\t12\t13\n')
            atomic_bytes(directory / 'tiny.vcf', b'##fileformat=VCFv4.2\n##contig=<ID=tiny,length=12>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\ntiny\t2\t.\tC\tT\t.\tPASS\t.\tGT\t1/.\n')
            outputs += [directory / x for x in ('tiny.fa', 'tiny.fa.fai', 'tiny.vcf')]
        atomic_json(outputs[0], {'schema_version': 1})
        publish_provenance(outputs[0], [], params, outputs)
        return
    if action in ('sample-name', 'capture-sample-tools', 'capture-extraction'):
        row = ds[args.dataset_id]
        directory = out / 'prepared' / args.dataset_id
        params.update(sample=row, phase=args.phase)
        if action == 'sample-name':
            outputs = [directory / 'exact_sample.txt']
            atomic_bytes(outputs[0], (row['sample_name'] + '\n').encode('utf-8'))
            inputs = []
        elif action == 'capture-sample-tools':
            inputs = [out / ('inputs/' + row['role'] + '.vcf.gz'), directory / 'exact_sample.txt', reference, str(reference) + '.fai']
            outputs = [directory / x for x in ('raw.records.tsv', 'raw.alt.tsv', 'sample.vcf.gz', 'sample.vcf.gz.csi', 'sample.alt.vcf.gz', 'sample.alt.vcf.gz.csi', 'records.tool.success')]
        else:
            inputs = [directory / 'regions.txt', reference, str(reference) + '.fai']
            outputs = [directory / 'extracted.fasta', directory / 'extraction.tool.success']
        publish_provenance(outputs[0], inputs, params, outputs)
        return
    if action == 'stage-reference':
        if reusable_action(reference, params, [reference]):
            return
        stage_reference(c['reference_fasta'], out)
        inputs, outputs = [c['reference_fasta']], [reference]
    elif action == 'index-reference':
        if args.phase == 'capture':
            index = Path(str(reference) + '.fai')
            read_fai(index)
            names = set()
            with open(reference, encoding='utf-8') as handle:
                for line in handle:
                    if line.startswith('>'):
                        name = line[1:].split()[0]
                        if name in names:
                            raise ValueError('Duplicate reference contig: ' + name)
                        names.add(name)
            outputs = [index]
        else:
            outputs = [index_reference(reference)]
        inputs = [reference]
    elif action == 'stage-vcf':
        staged = out / f'inputs/{args.role}.vcf.gz' if args.phase == 'capture' else stage_vcf(c[args.role + '_vcf'], out, args.role)
        inputs, outputs = [c[args.role + '_vcf']], [staged, Path(str(staged) + '.csi')]
        if args.phase == 'capture':
            outputs += [out / f'inputs/{args.role}.samples.txt', out / f'inputs/{args.role}.header.txt']
    elif action == 'headers':
        all_headers = {}
        for role in ('pacbio', 'illumina'):
            staged = out / f'inputs/{role}.vcf.gz'
            names = [row['sample_name'] for row in ds.values() if row['role'] == role]
            if args.phase == 'capture':
                all_headers[role] = validate_header_text(staged, reference, names, (out / f'inputs/{role}.samples.txt').read_text().splitlines(), (out / f'inputs/{role}.header.txt').read_text())
            else:
                all_headers[role] = validate_header(staged, reference, names)
            inputs += [staged]
            if args.phase == 'capture':
                inputs += [out / f'inputs/{role}.samples.txt', out / f'inputs/{role}.header.txt']
            write_tsv(out / f'manifests/{role}_header_samples.tsv', ({'sample_name': name} for name in all_headers[role]), ['sample_name'])
        inputs += [str(reference) + '.fai']
        atomic_json(out / 'manifests/headers.ok.json', all_headers)
        outputs = [out / 'manifests/headers.ok.json'] + [out / f'manifests/{role}_header_samples.tsv' for role in ('pacbio', 'illumina')]
    elif action in ('folds', 'validate-folds'):
        inputs = [reference, str(reference) + '.fai']
        params['grouping'] = c['folds']['grouping']
        if action == 'folds':
            build_folds(str(reference) + '.fai', sha256(reference), out, c['folds']['grouping'])
            outputs = [out / 'manifests/folds.json', out / 'manifests/fold_assignments.tsv', out / 'reports/fold_balance.tsv', out / 'reports/unmatched_contigs.tsv'] + [out / f'folds/fold{r}.bed' for r in range(1, 6)]
        else:
            verify_folds(str(reference) + '.fai', out, sha256(reference), c['folds']['grouping'])
            inputs += [out / 'manifests/folds.json'] + [out / f'folds/fold{r}.bed' for r in range(1, 6)]
            outputs = [out / 'manifests/folds.validated.json']
            atomic_json(outputs[0], {'assignment_sha256': json.loads((out / 'manifests/folds.json').read_text())['assignment_sha256']})
    elif action == 'seeds':
        rows = seed_rows(c)
        outputs = [out / 'manifests/seeds.tsv']
        write_tsv(outputs[0], rows, list(rows[0]))
        params['analysis'] = c['analysis']
        params['identities'] = c['illumina_corresponding']
    elif action == 'prepare-sample':
        row = ds[args.dataset_id]
        staged = out / f'inputs/{row["role"]}.vcf.gz'
        fingerprints = json.loads((out / 'manifests/selection_fingerprints.json').read_text(encoding='utf-8'))
        source_hash = fingerprints['inputs'][row['role'] + '_vcf']['sha256']
        prepare_sample(staged, reference, row['source_vcf'], source_hash, row['role'], row['sample_name'], out, c['analysis']['flank_bp'], out / 'manifests/folds.json', args.phase, header_validated=True)
        directory = out / 'prepared' / args.dataset_id
        if args.phase in ('records', 'records-from-query'):
            inputs = [staged, reference, str(reference) + '.fai', out / 'manifests/headers.ok.json']
            if args.phase == 'records-from-query':
                inputs += [directory / 'raw.records.tsv', directory / 'raw.alt.tsv', directory / 'records.tool.success']
            outputs = [directory / x for x in ('records.complete.json', 'records.tsv.gz')]
        elif args.phase == 'contexts-plan':
            inputs = [directory / 'records.tsv.gz', reference, str(reference) + '.fai', out / 'manifests/folds.json']
            outputs = [directory / 'contexts.plan.json', directory / 'regions.txt']
        else:
            inputs = [directory / 'records.tsv.gz', reference, str(reference) + '.fai', out / 'manifests/folds.json']
            if args.phase == 'contexts-from-fasta':
                inputs.append(directory / 'extracted.fasta')
                inputs.append(directory / 'extraction.tool.success')
            outputs = [directory / x for x in ('sample.complete.json', 'genotypes.tsv.gz', 'windows.tsv', 'all.fasta', 'fold_counts.tsv')] + [directory / f'fold{r}.{suffix}' for r in range(1, 6) for suffix in ('fasta', 'windows.tsv')]
        params.update(sample=row, flank=c['analysis']['flank_bp'], phase=args.phase)
    elif action == 'gate':
        verify_folds(str(reference) + '.fai', out, sha256(reference), c['folds']['grouping'])
        inputs = [out / 'manifests/headers.ok.json', out / 'manifests/folds.validated.json', out / 'manifests/preflight.ok.json'] + [out / 'prepared' / did / 'sample.complete.json' for did in ds]
        for path in inputs:
            json.loads(path.read_text(encoding='utf-8'))
        counts = [row for did in ds for row in read_tsv(out / 'prepared' / did / 'fold_counts.tsv')]
        write_tsv(out / 'reports/dataset_fold_counts.tsv', counts, ['dataset_id', 'sample_name', 'role', 'fold_id', 'eligible_windows'])
        def excluded_rows():
            for did in ds:
                for row in read_tsv(out / 'prepared' / did / 'genotypes.tsv.gz'):
                    if row['gt_status'] == 'alt' and (row['fold_status'] == 'excluded' or row['extraction_status'] == 'not_scanned'):
                        yield row | {'reason': row['extraction_reason'] or 'unmatched_group_regex'}
        write_tsv(out / 'reports/exclusions.tsv', excluded_rows(), ['dataset_id', 'sample_name', 'record_id', 'CHROM', 'POS', 'fold_status', 'extraction_status', 'reason'])
        outputs = [out / 'manifests/preparation.ok.json', out / 'reports/dataset_fold_counts.tsv', out / 'reports/exclusions.tsv']
        atomic_json(outputs[0], {'schema_version': 1, 'datasets': sorted(ds), 'inputs': {str(path): sha256(path) for path in inputs}})
    elif action == 'preflight':
        outputs = [out / f'manifests/preflight/{args.tool}.ok.json']
        directory = out / 'preflight' / args.tool
        if args.phase == 'capture':
            versions = {'bcftools': '1.23.1', 'samtools': '1.23.1', 'r': '4.5.1'}
            version = (directory / 'version.txt').read_text()
            if versions[args.tool] not in version.splitlines()[0]:
                raise ValueError('Pinned version mismatch: ' + args.tool + ': ' + version)
            if args.tool == 'samtools' and ''.join((directory / 'extracted.fasta').read_text().splitlines()[1:]) != 'CGTA':
                raise ValueError('samtools coordinate preflight failed')
            if args.tool == 'bcftools' and list(query_file_rows(directory / 'query.tsv'))[0]['GT'] != '1/.':
                raise ValueError('bcftools GT/query preflight failed')
            if args.tool == 'r' and (directory / 'result.txt').read_text().strip() != '6':
                raise ValueError('R worker preflight failed')
            inputs = [directory / 'version.txt', directory / 'tool.success']
            inputs.append(directory / {'bcftools': 'query.tsv', 'samtools': 'extracted.fasta', 'r': 'result.txt'}[args.tool])
            publish_provenance(directory / 'tool.success', [out / 'preflight/fixtures.complete.json'], params, inputs)
            result = {'tool': args.tool, 'expected_version': versions[args.tool], 'version': version.strip()}
        else:
            result = preflight(args.tool, directory)
        atomic_json(outputs[0], result)
        params['tool'] = args.tool
    elif action == 'preflight-gate':
        inputs = [out / f'manifests/preflight/{tool}.ok.json' for tool in ('bcftools', 'samtools', 'meme', 'python', 'r')]
        outputs = [out / 'manifests/preflight.ok.json']
        atomic_json(outputs[0], {path.stem: json.loads(path.read_text()) for path in inputs})
    else:
        raise ValueError('Unknown preparation action: ' + action)
    publish_provenance(outputs[0], inputs, params, outputs)


if __name__ == '__main__':
    main()
