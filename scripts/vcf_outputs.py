"""Original-record sample exports and a separate, traceable allele identity stream."""
import argparse
import gzip
import json
import os
import re
import sqlite3
import tempfile
from urllib.parse import quote
from fractions import Fraction
from pathlib import Path

from annotate import FIELDS, annotation_paths, verify_freeze
from prepare import (atomic_json, check_reuse, configure_scratch, dataset_id,
                     fingerprint_inputs, input_signature, parse_gt, publish_provenance,
                     read_fai, read_tsv, sha256, write_tsv)

ALLELE_FIELDS = FIELDS + ['source_sha256', 'allele_id', 'original_alt_index',
                         'original_alt', 'original_format', 'original_sample',
                         'original_record_fields', 'annotation_path', 'annotation_sha256',
                         'freeze_sha256', 'canonical_CHROM', 'canonical_POS',
                         'canonical_REF', 'canonical_ALT', 'normalization_status',
                         'exclusion_reason']


def load_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def verify_artifact(primary, required=()):
    """Require current inputs and all declared output checksums, not just a marker."""
    sidecar = Path(str(primary) + '.provenance.json')
    try:
        record = load_json(sidecar)
        params = record['signature_params']
        if record['scope'] != 'application':
            raise ValueError('not application provenance')
        if any(sha256(path) != digest for path, digest in params.get('code_hashes', {}).items()):
            raise ValueError('producer code changed')
        signature = input_signature(record['input_paths'], params)
        if not set(map(str, required)) <= set(record['required_outputs']):
            raise ValueError('unchecked required outputs')
        if not check_reuse(sidecar, signature, record['required_outputs']):
            raise ValueError('stale input signature or output checksum')
        return record
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise ValueError(f'Incomplete or stale artifact {primary}: {error}') from error


def provenance(primary, inputs, params, outputs):
    params = dict(params, code_paths=[str(Path(__file__).resolve()),
                                    str(Path(__file__).with_name('annotate.py').resolve()),
                                    str(Path(__file__).with_name('prepare.py').resolve())])
    return publish_provenance(primary, inputs, params, outputs, scope='application')


def publish_complete(directory, name, inputs, params, outputs, details):
    for path in outputs:
        if path.suffix != '.csi':
            provenance(path, inputs, params, outputs)
    marker = directory / name
    atomic_json(marker, dict(schema_version=1, outputs={str(p): sha256(p) for p in outputs}, **details))
    provenance(marker, inputs, params, outputs + [marker])
    return marker


def manifest(config, output_dir):
    output_dir = Path(output_dir)
    freeze = verify_freeze(output_dir)
    freeze_path = output_dir / 'library/freeze.json'
    freeze_hash = sha256(freeze_path)
    datasets, inputs = [], [freeze_path, Path(freeze['_reference']), Path(freeze['_reference'] + '.fai')]
    expected_corresponding = [name for names in config['illumina_corresponding'].values() for name in names]
    frozen_corresponding = [name for names in freeze['selection_config']['illumina_corresponding'].values() for name in names]
    if config['illumina_corresponding'] != freeze['selection_config']['illumina_corresponding']:
        raise ValueError('Corresponding membership differs from frozen selection')
    for role, source_key in [('illumina', 'illumina_vcf'), ('other_illumina', 'other_illumina')]:
        app, _ = annotation_paths(output_dir, role)
        application = load_json(app / 'manifest.json')
        annotation = app / 'genotypes.tsv.gz'
        complete = app / 'annotation.complete.json'
        verify_artifact(complete, [complete, annotation])
        completion = load_json(complete)
        if completion.get('role') != role or completion.get('freeze_sha256') != freeze_hash:
            raise ValueError('Annotation freeze/role mismatch: ' + str(complete))
        source = Path(config[source_key])
        source_hash = sha256(source)
        if str(source) != application['source_vcf'] or source_hash != application['source_sha256']:
            raise ValueError('Source differs from annotated input: ' + str(source))
        names = [row['sample_name'] for row in application['datasets']]
        if not names or len(names) != len(set(names)):
            raise ValueError('Application manifest requires unique nonempty samples: ' + role)
        if role == 'illumina' and names != expected_corresponding:
            raise ValueError('Corresponding manifest membership differs from frozen selection: ' + repr(frozen_corresponding))
        annotation_hash = sha256(annotation)
        for row in application['datasets']:
            name = row['sample_name']
            did = dataset_id(role, name)
            if row['dataset_id'] != did or row['role'] != role or row['source_vcf'] != str(source):
                raise ValueError('Application dataset identity mismatch: ' + name)
            context = 'selection_reuse' if role == 'illumina' else (
                'external_application_metadata_checked' if application['identity_check_status'] == 'explicit_clone_metadata'
                else 'external_application_identity_unverified')
            datasets.append(dict(dataset_id=did, role=role, sample_name=name, source_vcf=str(source),
                                 source_sha256=source_hash, annotation_path=str(annotation),
                                 annotation_sha256=annotation_hash, annotation_complete=str(complete),
                                 freeze_sha256=freeze_hash, selection_participant=role == 'illumina',
                                 evidence_context=context, identity_check_status=application['identity_check_status']))
        inputs += [source, app / 'manifest.json', annotation, complete, Path(str(complete) + '.provenance.json')]
    if len({row['dataset_id'] for row in datasets}) != len(datasets):
        raise ValueError('Dataset ID collision')
    settings = config.get('vcf_outputs', {'consensus_fraction': 0.75, 'groups': {}})
    value = settings.get('consensus_fraction', 0.75)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError('consensus_fraction must be a finite number in (0,1]')
    try:
        fraction = Fraction(str(value))
    except ValueError as error:
        raise ValueError('consensus_fraction must be a finite number in (0,1]') from error
    if not 0 < fraction <= 1:
        raise ValueError('consensus_fraction must be in (0,1]')
    membership = {(row['role'], row['sample_name']): row['dataset_id'] for row in datasets}
    configured = settings.get('groups', {})
    if not isinstance(configured, dict):
        raise ValueError('vcf_outputs.groups must be a mapping')
    groups = []
    for group_id, group in configured.items():
        if not isinstance(group_id, str) or not re.fullmatch(r'[a-z0-9][a-z0-9_-]*', group_id):
            raise ValueError('Unsafe VCF group ID: ' + repr(group_id))
        if not isinstance(group, dict) or group.get('role') not in ('illumina', 'other_illumina'):
            raise ValueError('Invalid role for VCF group: ' + group_id)
        names = group.get('samples')
        if not isinstance(names, list) or not names or any(not isinstance(name, str) or not name for name in names):
            raise ValueError('VCF group requires nonempty exact sample names: ' + group_id)
        if len(names) != len(set(names)):
            raise ValueError('Duplicate samples in VCF group: ' + group_id)
        unknown = [name for name in names if (group['role'], name) not in membership]
        if unknown:
            raise ValueError(f'Unknown {group["role"]} members in group {group_id}: {unknown!r}')
        n = len(names)
        groups.append(dict(group_id=group_id, role=group['role'], samples=names, n=n,
                           required_support=(fraction.numerator * n + fraction.denominator - 1) // fraction.denominator,
                           member_dataset_ids=[membership[group['role'], name] for name in names]))
    directory = output_dir / 'vcf_outputs'
    directory.mkdir(parents=True, exist_ok=True)
    # Dataset planning has no group parameters: regrouping cannot invalidate raw exports.
    fingerprint = fingerprint_inputs({str(p): p for p in inputs}, {'datasets': datasets}, directory / 'input_fingerprints.json')
    atomic_json(directory / 'datasets.json', {'schema_version': 1, 'datasets': datasets})
    provenance(directory / 'datasets.json', inputs, {'fingerprint_signature': fingerprint['signature']}, [directory / 'datasets.json'])
    atomic_json(directory / 'groups.json', {'schema_version': 1, 'consensus_fraction': str(fraction), 'groups': groups})
    provenance(directory / 'groups.json', inputs + [directory / 'datasets.json'],
               {'config_paths': ['vcf_outputs'], 'config_values': {'vcf_outputs': settings}}, [directory / 'groups.json'])
    return datasets, groups


def sample_context(output_dir, identifier):
    output_dir = Path(output_dir)
    plan = output_dir / 'vcf_outputs/datasets.json'
    verify_artifact(plan, [plan])
    rows = [row for row in load_json(plan)['datasets'] if row['dataset_id'] == identifier]
    if len(rows) != 1:
        raise ValueError('Unknown or duplicate dataset: ' + str(identifier))
    row = rows[0]
    freeze = verify_freeze(output_dir)
    freeze_path = output_dir / 'library/freeze.json'
    annotation = Path(row['annotation_path'])
    complete = Path(row['annotation_complete'])
    verify_artifact(complete, [complete, annotation])
    for path, digest in [(row['source_vcf'], row['source_sha256']), (annotation, row['annotation_sha256']),
                         (freeze_path, row['freeze_sha256'])]:
        if sha256(path) != digest:
            raise ValueError('Changed export input: ' + str(path))
    directory = output_dir / 'vcf_outputs' / identifier
    directory.mkdir(parents=True, exist_ok=True)
    inputs = [plan, Path(str(plan) + '.provenance.json'), Path(row['source_vcf']), annotation,
              complete, Path(str(complete) + '.provenance.json'), freeze_path,
              Path(freeze['_reference']), Path(freeze['_reference'] + '.fai')]
    return row, freeze, directory, inputs


def open_vcf(path):
    with open(path, 'rb') as handle:
        compressed = handle.read(2) == b'\x1f\x8b'
    return (gzip.open if compressed else open)(path, 'rt', encoding='utf-8', newline='')


def verify_bgzf_csi(vcf, index, raw=None):
    with Path(vcf).open('rb') as handle:
        header = handle.read(12)
        if len(header) != 12 or header[:3] != b'\x1f\x8b\x08' or not header[3] & 4:
            raise ValueError('VCF is not BGZF: ' + str(vcf))
        extra = handle.read(int.from_bytes(header[10:12], 'little'))
        offset, bgzf = 0, False
        while offset + 4 <= len(extra):
            size = int.from_bytes(extra[offset + 2:offset + 4], 'little')
            if offset + 4 + size > len(extra):
                raise ValueError('Malformed BGZF extra fields: ' + str(vcf))
            bgzf |= extra[offset:offset + 2] == b'BC' and size == 2
            offset += 4 + size
        if not bgzf or offset != len(extra):
            raise ValueError('VCF lacks the BGZF block header: ' + str(vcf))
    with Path(index).open('rb') as handle:
        compressed_index = handle.read(2) == b'\x1f\x8b'
    with (gzip.open if compressed_index else open)(index, 'rb') as handle:
        if handle.read(4) != b'CSI\x01':
            raise ValueError('Missing or malformed CSI: ' + str(index))
    if raw is not None:
        with Path(raw).open('rb') as original, gzip.open(vcf, 'rb') as compressed:
            while True:
                block = original.read(1024 * 1024)
                if compressed.read(len(block) if block else 1) != block:
                    raise ValueError('Compressed VCF differs from raw export: ' + str(vcf))
                if not block:
                    break


def exclusion_reason(ref, alt):
    if not re.fullmatch('[ACGTNacgtn]+', ref):
        raise ValueError('Malformed REF: ' + repr(ref))
    if alt in ('<NON_REF>', '<*>'):
        return 'unspecified_alt'
    if alt == '*':
        return 'spanning_deletion'
    if re.fullmatch(r'<[A-Za-z0-9_.:-]+>', alt):
        return 'symbolic_alt'
    if '[' in alt or ']' in alt:
        match = re.fullmatch(r'([ACGTNacgtn]*)([\[\]])([^\s\[\]]+):([1-9][0-9]*)\2([ACGTNacgtn]*)', alt)
        if not match or bool(match[1]) == bool(match[5]):
            raise ValueError('Malformed breakend ALT: ' + repr(alt))
        return 'breakend_alt'
    if re.fullmatch(r'(?:\.[ACGTNacgtn]+|[ACGTNacgtn]+\.)', alt):
        return 'breakend_alt'
    if not re.fullmatch('[ACGTNacgtn]+', alt):
        raise ValueError('Malformed ALT: ' + repr(alt))
    return 'ambiguous_sequence' if 'N' in (ref + alt).upper() else ''


def sample_inputs(config, output_dir, identifier):
    row, freeze, directory, inputs = sample_context(output_dir, identifier)
    contigs = read_fai(freeze['_reference'] + '.fai')
    with tempfile.TemporaryDirectory(prefix='.sample-', dir=directory) as scratch:
        scratch = Path(scratch)
        db = sqlite3.connect(scratch / 'annotations.sqlite')
        try:
            db.execute('CREATE TABLE annotations (ordinal INTEGER PRIMARY KEY, row TEXT NOT NULL)')
            previous = 0
            for annotation in read_tsv(row['annotation_path']):
                if annotation['dataset_id'] != identifier:
                    continue
                ordinal_text = annotation['source_record_index']
                if not re.fullmatch('[1-9][0-9]*', ordinal_text) or int(ordinal_text) <= previous:
                    raise ValueError('Annotations are not in unique source ordinal order: ' + identifier)
                previous = int(ordinal_text)
                db.execute('INSERT INTO annotations VALUES (?,?)', (previous, json.dumps(annotation)))
            db.commit()
            raw = scratch / 'promising.raw.vcf'
            sites = scratch / 'source_alleles.vcf'
            sidecar = scratch / 'source_alleles.tsv.gz'
            counts = {'source_record_count': 0, 'promising_record_count': 0, 'allele_count': 0, 'supported_allele_count': 0}
            with open_vcf(row['source_vcf']) as source, raw.open('w', encoding='utf-8', newline='') as export, sites.open('w', encoding='utf-8', newline='') as site:
                header = None
                for line in source:
                    if line.startswith('##'):
                        export.write(line)
                    elif line.startswith('#CHROM\t'):
                        header = line.rstrip('\r\n').split('\t')
                        break
                    else:
                        raise ValueError('Malformed VCF metadata/header: ' + row['source_vcf'])
                if header is None or header[:9] != ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']:
                    raise ValueError('Missing or malformed sample VCF header')
                names = header[9:]
                if len(names) != len(set(names)) or any(not name for name in names) or names.count(row['sample_name']) != 1:
                    raise ValueError('VCF header lacks unique exact sample: ' + row['sample_name'])
                sample_column = header.index(row['sample_name'], 9)
                export.write('##motif_label=prioritization_only\n')
                export.write('##motif_freeze_sha256=' + row['freeze_sha256'] + '\n')
                export.write('##motif_selection_participant=' + str(row['selection_participant']).lower() + '\n')
                export.write('##motif_evidence_context=' + row['evidence_context'] + '\n')
                export.write('\t'.join(header[:9] + [row['sample_name']]) + '\n')
                site.write('##fileformat=VCFv4.2\n')
                for chrom, length in contigs.items():
                    site.write(f'##contig=<ID={chrom},length={length}>\n')
                site.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

                def occurrences():
                    for ordinal, line in enumerate(source, 1):
                        fields = line.rstrip('\r\n').split('\t')
                        if len(fields) != len(header) or any(not value for value in fields):
                            raise ValueError(f'Malformed VCF row {ordinal}')
                        chrom, pos, record_id, ref, alts = fields[:5]
                        if chrom not in contigs or not re.fullmatch('[1-9][0-9]*', pos) or int(pos) + len(ref) - 1 > contigs[chrom]:
                            raise ValueError(f'Invalid VCF coordinates at source row {ordinal}')
                        if not re.fullmatch('[ACGTNacgtn]+', ref):
                            raise ValueError('Malformed REF at source row ' + str(ordinal))
                        alt_list = [] if alts == '.' else alts.split(',')
                        reasons = [exclusion_reason(ref, alt) for alt in alt_list]
                        keys = fields[8].split(':')
                        if fields[8] != '.' and (len(keys) != len(set(keys)) or any(not re.fullmatch('[A-Za-z_][A-Za-z0-9_.]*', key) for key in keys)):
                            raise ValueError('Malformed FORMAT at source row ' + str(ordinal))
                        values = fields[sample_column].split(':')
                        if len(values) > len(keys):
                            raise ValueError('Too many sample FORMAT values at source row ' + str(ordinal))
                        gt_index = keys.index('GT') if 'GT' in keys else None
                        gt = values[gt_index] if gt_index is not None and gt_index < len(values) and values[gt_index] else '.'
                        parsed = parse_gt(gt, len(alt_list))
                        found = db.execute('SELECT row FROM annotations WHERE ordinal=?', (ordinal,)).fetchone()
                        if found is None:
                            raise ValueError(f'Missing annotation for {identifier} source row {ordinal}')
                        annotation = json.loads(found[0])
                        expected = dict(source_vcf=row['source_vcf'], source_record_index=str(ordinal),
                                        record_id=f'{row["role"]}_{row["source_sha256"]}_{ordinal}',
                                        dataset_id=identifier, sample_name=row['sample_name'], role=row['role'],
                                        selection_participant=str(row['selection_participant']).lower(), evidence_context=row['evidence_context'],
                                        CHROM=chrom, POS=pos, ID=record_id, REF=ref, ALT=alts, GT=gt)
                        if any(annotation.get(key) != value for key, value in expected.items()):
                            bad = [key for key, value in expected.items() if annotation.get(key) != value]
                            raise ValueError(f'Source/annotation identity mismatch at row {ordinal}: {bad}')
                        counts['source_record_count'] += 1
                        if annotation['annotation_status'] == 'promising' and parsed['gt_status'] == 'alt':
                            export.write('\t'.join(fields[:9] + [fields[sample_column]]) + '\n')
                            counts['promising_record_count'] += 1
                        for alt_index, (alt, reason) in enumerate(zip(alt_list, reasons), 1):
                            allele_id = f'a{ordinal}_{alt_index}'
                            occurrence = dict(annotation, source_sha256=row['source_sha256'], allele_id=allele_id,
                                              original_alt_index=alt_index, original_alt=alt,
                                              original_format=fields[8], original_sample=fields[sample_column],
                                              original_record_fields=json.dumps(fields[:9] + [fields[sample_column]], separators=(',', ':')),
                                              annotation_path=row['annotation_path'], annotation_sha256=row['annotation_sha256'],
                                              freeze_sha256=row['freeze_sha256'], canonical_CHROM='', canonical_POS='',
                                              canonical_REF='', canonical_ALT='', normalization_status='excluded' if reason else 'pending',
                                              exclusion_reason=reason)
                            counts['allele_count'] += 1
                            if not reason:
                                site.write('\t'.join([chrom, pos, allele_id, ref.upper(), alt.upper(), '.', '.', '.']) + '\n')
                                counts['supported_allele_count'] += 1
                            yield occurrence
                    if counts['source_record_count'] != previous or db.execute('SELECT count(*) FROM annotations').fetchone()[0] != counts['source_record_count']:
                        raise ValueError('Annotation/source record counts disagree: ' + identifier)
                write_tsv(sidecar, occurrences(), ALLELE_FIELDS)
            if sha256(row['source_vcf']) != row['source_sha256'] or sha256(row['annotation_path']) != row['annotation_sha256']:
                raise ValueError('Source or annotation changed during export')
            outputs = [directory / path.name for path in (raw, sites, sidecar)]
            for temporary, destination in zip((raw, sites, sidecar), outputs):
                os.replace(temporary, destination)
            return publish_complete(directory, 'inputs.complete.json', inputs, {'dataset': row}, outputs,
                                    dict(dataset_id=identifier, source_sha256=row['source_sha256'], freeze_sha256=row['freeze_sha256'], **counts))
        finally:
            db.close()


def capture_alleles(config, output_dir, identifier, normalized_vcf=None):
    row, freeze, directory, inputs = sample_context(output_dir, identifier)
    marker = directory / 'inputs.complete.json'
    source_map = directory / 'source_alleles.tsv.gz'
    source_sites = directory / 'source_alleles.vcf'
    verify_artifact(marker, [marker, source_map, source_sites])
    normalized = Path(normalized_vcf) if normalized_vcf else directory / 'normalized_alleles.vcf'
    contigs = read_fai(freeze['_reference'] + '.fai')
    with tempfile.TemporaryDirectory(prefix='.alleles-', dir=directory) as scratch:
        scratch = Path(scratch)
        db = sqlite3.connect(scratch / 'alleles.sqlite')
        try:
            db.execute('CREATE TABLE alleles (id TEXT PRIMARY KEY, row TEXT NOT NULL, canonical TEXT)')
            expected_count = 0
            for occurrence in read_tsv(source_map):
                if occurrence['normalization_status'] not in ('pending', 'excluded'):
                    raise ValueError('Invalid source normalization status')
                db.execute('INSERT INTO alleles VALUES (?,?,NULL)', (occurrence['allele_id'], json.dumps(occurrence)))
                expected_count += occurrence['normalization_status'] == 'pending'
            db.commit()
            observed = 0
            with open_vcf(normalized) as handle:
                header_seen = False
                for line in handle:
                    if line.startswith('##') and not header_seen:
                        continue
                    if line.startswith('#CHROM\t') and not header_seen:
                        if line.rstrip('\r\n').split('\t') != ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']:
                            raise ValueError('Normalization must return a site-only VCF')
                        header_seen = True
                        continue
                    fields = line.rstrip('\r\n').split('\t')
                    if not header_seen or len(fields) != 8:
                        raise ValueError('Malformed normalized site VCF')
                    chrom, pos, allele_id, ref, alt = fields[:5]
                    if chrom not in contigs or not re.fullmatch('[1-9][0-9]*', pos) or not re.fullmatch('[ACGT]+', ref) or not re.fullmatch('[ACGT]+', alt) or int(pos) + len(ref) - 1 > contigs[chrom]:
                        raise ValueError('Invalid normalized canonical allele: ' + allele_id)
                    found = db.execute('SELECT row,canonical FROM alleles WHERE id=?', (allele_id,)).fetchone()
                    if found is None or json.loads(found[0])['normalization_status'] != 'pending' or found[1] is not None:
                        raise ValueError('Unknown, excluded or duplicate normalization ID: ' + allele_id)
                    original = json.loads(found[0])
                    if original['CHROM'] != chrom:
                        raise ValueError('Normalization changed contig: ' + allele_id)
                    db.execute('UPDATE alleles SET canonical=? WHERE id=?', (json.dumps([chrom, pos, ref, alt]), allele_id))
                    observed += 1
                if not header_seen or observed != expected_count:
                    raise ValueError(f'Lost normalization IDs: expected {expected_count}, observed {observed}')
            db.commit()

            def normalized_rows():
                for occurrence in read_tsv(source_map):
                    if occurrence['normalization_status'] == 'pending':
                        canonical = db.execute('SELECT canonical FROM alleles WHERE id=?', (occurrence['allele_id'],)).fetchone()[0]
                        if canonical is None:
                            raise ValueError('Lost normalization ID: ' + occurrence['allele_id'])
                        occurrence.update(zip(['canonical_CHROM', 'canonical_POS', 'canonical_REF', 'canonical_ALT'], json.loads(canonical)))
                        occurrence['normalization_status'] = 'normalized'
                    yield occurrence
            temporary = scratch / 'alleles.tsv.gz'
            write_tsv(temporary, normalized_rows(), ALLELE_FIELDS)
            output = directory / 'alleles.tsv.gz'
            os.replace(temporary, output)
            inputs += [marker, Path(str(marker) + '.provenance.json'), source_map, source_sites, normalized]
            return publish_complete(directory, 'alleles.complete.json', inputs,
                                    {'dataset': row, 'normalization': 'bcftools norm --no-version -f frozen_reference -c e -Ov'}, [output],
                                    {'dataset_id': identifier, 'normalized_allele_count': observed,
                                     'source_sha256': row['source_sha256'], 'freeze_sha256': row['freeze_sha256']})
        finally:
            db.close()


def capture_export(config, output_dir, identifier, bgzip_version=None):
    row, _, directory, inputs = sample_context(output_dir, identifier)
    marker = directory / 'inputs.complete.json'
    raw = directory / 'promising.raw.vcf'
    verify_artifact(marker, [marker, raw])
    version_path = Path(bgzip_version) if bgzip_version else directory / 'bgzip.version.txt'
    version = version_path.read_text(encoding='utf-8').strip()
    if not version or not re.search(r'\bbgzip\b', version, re.IGNORECASE):
        raise ValueError('Missing observed bgzip --version output: ' + str(version_path))
    vcf, index = directory / 'promising.vcf.gz', directory / 'promising.vcf.gz.csi'
    verify_bgzf_csi(vcf, index, raw)
    inputs += [marker, Path(str(marker) + '.provenance.json'), raw, version_path]
    return publish_complete(directory, 'export.complete.json', inputs, {'dataset': row, 'bgzip_version': version},
                            [vcf, index], {'dataset_id': identifier, 'source_sha256': row['source_sha256'],
                                           'freeze_sha256': row['freeze_sha256'], 'bgzip_version': version})


def project_gt(gt, target_indexes, alt_count):
    alleles = parse_gt(gt, alt_count)['alleles']
    targets = set(target_indexes)
    if any(type(index) is not int or not 1 <= index <= alt_count for index in targets):
        raise ValueError('Invalid target ALT index')
    parts = re.split(r'([/|])', gt)
    projected = []
    for part in parts:
        if part in ('/', '|'):
            projected.append(part)
        elif part == '.':
            projected.append('.')
        else:
            index = int(part)
            projected.append('0' if index == 0 else ('1' if index in targets else '.'))
    if len([part for part in projected if part not in ('/', '|')]) != len(alleles):
        raise ValueError('GT projection changed ploidy')
    return ''.join(projected)


def resolve_member(occurrences):
    """Reduce duplicate source records by allele presence, then resolve GT separately."""
    observations = []
    for occurrence in sorted(occurrences, key=lambda row: int(row['source_record_index'])):
        gt = occurrence.get('GT', '.')
        parsed = parse_gt(gt, occurrence['alt_count'])
        targets = set(occurrence['target_indexes'])
        called_targets = sum(value in targets for value in parsed['alleles'] if value is not None)
        presence = 'present' if called_targets else ('unknown' if None in parsed['alleles'] else 'absent')
        reference_count = parsed['alleles'].count(0)
        missing_count = parsed['alleles'].count(None)
        other = []
        for value in parsed['alleles']:
            if value is not None and value > 0 and value not in targets:
                other.append(occurrence['allele_identities'][value])
        signature = (parsed['ploidy'], reference_count, called_targets, missing_count,
                     tuple(sorted(json.dumps(value, separators=(',', ':'), ensure_ascii=False) for value in other)))
        projected = project_gt(gt, targets, occurrence['alt_count'])
        item = dict(occurrence, parsed=parsed, target_presence=presence, signature=signature,
                    projected_GT=projected, target_copies=called_targets,
                    promising=occurrence.get('annotation_status') == 'promising')
        observations.append(item)

    if not observations:
        return dict(call_state='absent', presence_state='unknown', presence_conflict=False,
                    projected_GT='.', support=False, duplicate_count=0, dosage_disagreement=False,
                    ploidy_disagreement=False, duplicate_genotype_disagreement=False,
                    missingness=True, annotation_disagreement=False, source_occurrences=[],
                    comparison_signature=None, has_called_allele=False)

    presences = {item['target_presence'] for item in observations}
    present = 'present' in presences
    absent = 'absent' in presences
    conflict = present and absent
    if conflict:
        presence_state = 'conflict'
    elif present:
        presence_state = 'present'
    elif absent:
        presence_state = 'absent'
    else:
        presence_state = 'unknown'

    complete = [item for item in observations if None not in item['parsed']['alleles']]
    target_dosages = {item['target_copies'] for item in complete}
    dosage_disagreement = len(target_dosages) > 1
    ploidies = {item['parsed']['ploidy'] for item in observations if item['GT'] != '.'}
    ploidy_disagreement = len(ploidies) > 1
    complete_signatures = {item['signature'] for item in complete}
    duplicate_genotype_disagreement = len(complete_signatures) > 1
    signatures = {item['signature'] for item in observations}
    missingness = any(None in item['parsed']['alleles'] for item in observations)
    annotation_disagreement = len({item.get('annotation_status', '') for item in observations}) > 1
    if conflict:
        call_state = 'duplicate_conflict'
    elif present:
        call_state = 'target_alt'
    elif any(value is not None and value > 0 for item in observations for value in item['parsed']['alleles']):
        call_state = 'other_alt_only'
    elif any(value == 0 for item in observations for value in item['parsed']['alleles']):
        call_state = 'reference_only'
    else:
        call_state = 'missing'

    support = present and not absent and any(
        item['target_presence'] == 'present' and item['promising'] for item in observations)
    if len(signatures) == 1 and not conflict:
        projected_gt = observations[0]['projected_GT']
        comparison_signature = observations[0]['signature']
    else:
        projected_gt = '/'.join(['.'] * next(iter(ploidies))) if len(ploidies) == 1 else '.'
        comparison_signature = None
    source_occurrences = []
    for item in observations:
        source = dict(item.get('source_occurrence', {}))
        source['target_presence'] = item['target_presence']
        source['target_alt_indexes'] = sorted(item['target_indexes'])
        source_occurrences.append(source)
    return dict(call_state=call_state, presence_state=presence_state, presence_conflict=conflict,
                projected_GT=projected_gt, support=support, duplicate_count=len(observations),
                dosage_disagreement=dosage_disagreement, ploidy_disagreement=ploidy_disagreement,
                duplicate_genotype_disagreement=duplicate_genotype_disagreement,
                missingness=missingness, annotation_disagreement=annotation_disagreement,
                source_occurrences=source_occurrences, comparison_signature=comparison_signature,
                has_called_allele=any(value is not None for item in observations for value in item['parsed']['alleles']))


def group_context(output_dir, group_id):
    output_dir = Path(output_dir)
    datasets_path = output_dir / 'vcf_outputs/datasets.json'
    groups_path = output_dir / 'vcf_outputs/groups.json'
    verify_artifact(datasets_path, [datasets_path])
    verify_artifact(groups_path, [groups_path])
    datasets = load_json(datasets_path)['datasets']
    matches = [group for group in load_json(groups_path)['groups'] if group['group_id'] == group_id]
    if len(matches) != 1:
        raise ValueError('Unknown or duplicate VCF group: ' + group_id)
    group = matches[0]
    freeze = verify_freeze(output_dir)
    freeze_path = output_dir / 'library/freeze.json'
    member_rows = {row['dataset_id']: row for row in datasets}
    if len(group['member_dataset_ids']) != group['n'] or group['n'] < 1:
        raise ValueError('Invalid fixed group denominator: ' + group_id)
    if any(did not in member_rows or member_rows[did]['role'] != group['role'] for did in group['member_dataset_ids']):
        raise ValueError('Group members differ from validated dataset manifest: ' + group_id)
    directory = output_dir / 'vcf_outputs/groups' / group_id
    directory.mkdir(parents=True, exist_ok=True)
    inputs = [datasets_path, Path(str(datasets_path) + '.provenance.json'),
              groups_path, Path(str(groups_path) + '.provenance.json'), freeze_path,
              Path(freeze['_reference']), Path(freeze['_reference'] + '.fai')]
    members = []
    for did in group['member_dataset_ids']:
        row = member_rows[did]
        sample_dir = output_dir / 'vcf_outputs' / did
        export = sample_dir / 'export.complete.json'
        export_vcf = sample_dir / 'promising.vcf.gz'
        export_index = Path(str(export_vcf) + '.csi')
        alleles = sample_dir / 'alleles.complete.json'
        allele_table = sample_dir / 'alleles.tsv.gz'
        verify_artifact(export, [export, export_vcf, export_index])
        verify_artifact(alleles, [alleles, allele_table])
        export_record, allele_record = load_json(export), load_json(alleles)
        if (export_record.get('dataset_id') != did or allele_record.get('dataset_id') != did or
                export_record.get('source_sha256') != row['source_sha256'] or
                allele_record.get('source_sha256') != row['source_sha256'] or
                export_record.get('freeze_sha256') != row['freeze_sha256'] or
                allele_record.get('freeze_sha256') != row['freeze_sha256']):
            raise ValueError('Member completion identity mismatch: ' + did)
        annotation_complete = Path(row['annotation_complete'])
        verify_artifact(annotation_complete, [annotation_complete, Path(row['annotation_path'])])
        members.append(row)
        inputs += [Path(row['source_vcf']), Path(row['annotation_path']), annotation_complete,
                   Path(str(annotation_complete) + '.provenance.json'), export, Path(str(export) + '.provenance.json'),
                   export_vcf, export_index, alleles, Path(str(alleles) + '.provenance.json'), allele_table]
    return group, members, freeze, directory, inputs


def _canonical_key(row):
    return json.dumps([row['canonical_CHROM'], row['canonical_POS'], row['canonical_REF'], row['canonical_ALT']],
                      separators=(',', ':'))


def _record_observation(rows, target_key):
    first = rows[0]
    raw_fields = json.loads(first['original_record_fields'])
    alts = [] if raw_fields[4] == '.' else raw_fields[4].split(',')
    alt_rows = {int(row['original_alt_index']): row for row in rows}
    if len(alt_rows) != len(rows) or set(alt_rows) != set(range(1, len(alts) + 1)):
        raise ValueError('Missing or duplicate original ALT index in allele map')
    if any(row['original_record_fields'] != first['original_record_fields'] or row['GT'] != first['GT']
           for row in rows[1:]):
        raise ValueError('Inconsistent source-record identity in allele map')
    target_indexes = {index for index, row in alt_rows.items()
                      if row['normalization_status'] == 'normalized' and _canonical_key(row) == target_key}
    if not target_indexes:
        raise ValueError('Canonical candidate lacks source ALT index')
    identities, alt_alleles = {}, []
    for index, alt in enumerate(alts, 1):
        row = alt_rows[index]
        if row['normalization_status'] == 'normalized':
            identities[index] = ('canonical', row['canonical_CHROM'], row['canonical_POS'],
                                 row['canonical_REF'], row['canonical_ALT'])
            canonical = [row['canonical_CHROM'], row['canonical_POS'], row['canonical_REF'], row['canonical_ALT']]
        else:
            identities[index] = ('raw', row['CHROM'], row['POS'], row['REF'], row['original_alt'])
            canonical = None
        alt_alleles.append(dict(alt_index=index, allele=alt, canonical=canonical,
                                normalization_status=row['normalization_status'],
                                exclusion_reason=row['exclusion_reason']))
    source = dict(source_vcf=first['source_vcf'], source_sha256=first['source_sha256'],
                  source_record_index=first['source_record_index'], source_record_id=first['record_id'],
                  original_CHROM=first['CHROM'], original_POS=first['POS'], original_ID=first['ID'],
                  original_REF=first['REF'], original_ALT=first['ALT'], original_GT=first['GT'],
                  original_format=first['original_format'], original_sample=first['original_sample'],
                  original_ALT_list=alts, original_alt_indexes=list(range(1, len(alts) + 1)),
                  original_alt_alleles=alt_alleles,
                  window_id=first['window_id'], interval_start=first['interval_start'],
                  interval_end=first['interval_end'], annotation_status=first['annotation_status'],
                  annotation_reason=first['annotation_reason'], best_final_motif_id=first['best_final_motif_id'],
                  best_family_id=first['best_family_id'], reported_p_value=first['reported_p_value'],
                  local_start=first['local_start'], local_stop=first['local_stop'],
                  annotation_path=first['annotation_path'], annotation_sha256=first['annotation_sha256'],
                  freeze_sha256=first['freeze_sha256'])
    return dict(source_record_index=first['source_record_index'], GT=first['GT'], alt_count=len(alts),
                target_indexes=target_indexes, allele_identities=identities,
                annotation_status=first['annotation_status'], source_occurrence=source)


def _member_evidence_row(group, candidate, member, result):
    return dict(group_id=group['group_id'], role=group['role'], dataset_id=member['dataset_id'],
                sample_name=member['sample_name'], CHROM=candidate[0], POS=candidate[1],
                REF=candidate[2], ALT=candidate[3], selection_participant=str(member['selection_participant']).lower(),
                evidence_context=member['evidence_context'], call_state=result['call_state'],
                presence_state=result['presence_state'], presence_conflict=str(result['presence_conflict']).lower(),
                projected_GT=result['projected_GT'], support=str(result['support']).lower(),
                duplicate_count=result['duplicate_count'], dosage_disagreement=str(result['dosage_disagreement']).lower(),
                ploidy_disagreement=str(result['ploidy_disagreement']).lower(),
                duplicate_genotype_disagreement=str(result['duplicate_genotype_disagreement']).lower(),
                missingness=str(result['missingness']).lower(),
                annotation_disagreement=str(result['annotation_disagreement']).lower(),
                source_occurrences=json.dumps(result['source_occurrences'], ensure_ascii=False, separators=(',', ':')))


def _support_row(group, candidate, results, members, genotype_disagreement):
    supporters = [member['sample_name'] for member, result in zip(members, results) if result['support']]
    conflicts = [member['sample_name'] for member, result in zip(members, results) if result['presence_conflict']]
    count, n = len(supporters), group['n']
    return dict(group_id=group['group_id'], role=group['role'], CHROM=candidate[0], POS=candidate[1],
                REF=candidate[2], ALT=candidate[3], group_n=n, required_support=group['required_support'],
                support_count=count, support_percent=f'{100 * count / n:.6f}',
                supporting_samples=json.dumps(supporters, ensure_ascii=False, separators=(',', ':')),
                genotype_disagreement=str(genotype_disagreement).lower(),
                conflicting_samples=json.dumps(conflicts, ensure_ascii=False, separators=(',', ':')),
                retained=str(count >= group['required_support']).lower())


def write_consensus(config, output_dir, group_id):
    output_dir = Path(output_dir)
    group, members, freeze, directory, _ = group_context(output_dir, group_id)
    reference_order = {chrom: index for index, chrom in enumerate(read_fai(freeze['_reference'] + '.fai'))}
    settings = config.get('vcf_outputs', {'consensus_fraction': 0.75, 'groups': {}})
    fraction = Fraction(str(settings['consensus_fraction']))
    required = (fraction.numerator * group['n'] + fraction.denominator - 1) // fraction.denominator
    if required != group['required_support']:
        raise ValueError('Group threshold differs from validated configuration: ' + group_id)
    scratch_parent = Path(tempfile.gettempdir())
    with tempfile.TemporaryDirectory(prefix='.consensus-', dir=scratch_parent) as scratch:
        db = sqlite3.connect(Path(scratch) / 'alleles.sqlite')
        try:
            db.execute('CREATE TABLE alleles (dataset_id TEXT NOT NULL, ordinal INTEGER NOT NULL, alt_index INTEGER NOT NULL, canonical_key TEXT, canonical_chrom TEXT, canonical_pos INTEGER, canonical_ref TEXT, canonical_alt TEXT, status TEXT NOT NULL, row_json TEXT NOT NULL, PRIMARY KEY(dataset_id, ordinal, alt_index))')
            db.execute('CREATE INDEX allele_candidate ON alleles(canonical_key, dataset_id, ordinal)')
            db.execute('CREATE TABLE candidates (canonical_key TEXT PRIMARY KEY, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, contig_order INTEGER)')
            db.execute('CREATE TABLE evidence_rows (row_id INTEGER PRIMARY KEY, row_json TEXT NOT NULL)')
            for member in members:
                table = Path(output_dir) / 'vcf_outputs' / member['dataset_id'] / 'alleles.tsv.gz'
                for row in read_tsv(table):
                    ordinal, alt_index = int(row['source_record_index']), int(row['original_alt_index'])
                    if row['dataset_id'] != member['dataset_id'] or row['source_sha256'] != member['source_sha256']:
                        raise ValueError('Allele map member/source mismatch: ' + member['dataset_id'])
                    status = row['normalization_status']
                    if status not in ('normalized', 'excluded'):
                        raise ValueError('Incomplete normalization map: ' + member['dataset_id'])
                    if status == 'normalized':
                        key = _canonical_key(row)
                        chrom, pos, ref, alt = (row['canonical_CHROM'], row['canonical_POS'],
                                               row['canonical_REF'], row['canonical_ALT'])
                        if chrom not in reference_order:
                            raise ValueError('Unknown canonical contig: ' + chrom)
                        raw_fields = json.loads(row['original_record_fields'])
                        gt = row['GT']
                        parsed = parse_gt(gt, len(raw_fields[4].split(',')) if raw_fields[4] != '.' else 0)
                        if alt_index in parsed['alleles']:
                            db.execute('INSERT OR IGNORE INTO candidates VALUES (?,?,?,?,?,?)',
                                       (key, chrom, int(pos), ref, alt, reference_order[chrom]))
                    else:
                        key = chrom = pos = ref = alt = None
                    db.execute('INSERT INTO alleles VALUES (?,?,?,?,?,?,?,?,?,?)',
                               (member['dataset_id'], ordinal, alt_index, key, chrom, int(pos) if pos else None,
                                ref, alt, status, json.dumps(row, ensure_ascii=False, separators=(',', ':'))))
            db.commit()
            support_path = directory / 'support.tsv.gz'
            evidence_path = directory / 'member_evidence.tsv.gz'
            excluded_path = directory / 'excluded_alleles.tsv.gz'
            support_fields = ['group_id', 'role', 'CHROM', 'POS', 'REF', 'ALT', 'group_n', 'required_support',
                              'support_count', 'support_percent', 'supporting_samples', 'genotype_disagreement',
                              'conflicting_samples', 'retained']
            evidence_fields = ['group_id', 'role', 'dataset_id', 'sample_name', 'CHROM', 'POS', 'REF', 'ALT',
                               'selection_participant', 'evidence_context', 'call_state', 'presence_state',
                               'presence_conflict', 'projected_GT', 'support', 'duplicate_count',
                               'dosage_disagreement', 'ploidy_disagreement', 'duplicate_genotype_disagreement',
                               'missingness', 'annotation_disagreement', 'source_occurrences']
            excluded_fields = ['group_id', 'role', 'dataset_id', 'sample_name', 'source_vcf', 'source_sha256',
                               'source_record_index', 'record_id', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'GT',
                               'original_alt_index', 'original_alt', 'annotation_status', 'window_id',
                               'interval_start', 'interval_end', 'selection_participant', 'evidence_context',
                               'exclusion_reason']
            raw_vcf = directory / 'consensus.raw.vcf'
            with raw_vcf.open('w', encoding='utf-8', newline='') as vcf:
                vcf.write('##fileformat=VCFv4.2\n')
                vcf.write('##motif_label=prioritization_only\n')
                vcf.write('##motif_freeze_sha256=' + sha256(Path(output_dir) / 'library/freeze.json') + '\n')
                vcf.write('##motif_consensus_group=' + group_id + '\n')
                vcf.write(f'##motif_consensus_fraction={fraction.numerator}/{fraction.denominator}\n')
                vcf.write('##motif_role=' + group['role'] + '\n')
                participant = group['role'] == 'illumina'
                vcf.write('##motif_selection_participant=' + str(participant).lower() + '\n')
                contexts = sorted({member['evidence_context'] for member in members})
                vcf.write('##motif_evidence_context=' + (contexts[0] if len(contexts) == 1 else 'mixed') + '\n')
                for chrom, length in read_fai(freeze['_reference'] + '.fai').items():
                    vcf.write(f'##contig=<ID={chrom},length={length}>\n')
                info_definitions = [
                    ('GROUP_N', '1', 'Integer', 'Configured group member count'),
                    ('REQUIRED_SUPPORT', '1', 'Integer', 'Minimum promising member support'),
                    ('PROMISING_SUPPORT', '1', 'Integer', 'Members with a promising target-ALT occurrence'),
                    ('PROMISING_PERCENT', '1', 'Float', 'Promising support as percent of all group members'),
                    ('PROMISING_SAMPLES', '.', 'String', 'Percent-escaped exact sample names supporting this allele'),
                    ('CONFLICTING_SAMPLES', '.', 'String', 'Percent-escaped samples with target-presence conflicts'),
                ]
                for name, number, kind, description in info_definitions:
                    vcf.write(f'##INFO=<ID={name},Number={number},Type={kind},Description="{description}">\n')
                vcf.write('##INFO=<ID=GT_DISAGREEMENT,Number=0,Type=Flag,Description="Genotype or ploidy disagreement; does not alter support">\n')
                vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Projected individual genotype for this canonical ALT">\n')
                vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +
                          '\t'.join(member['sample_name'] for member in members) + '\n')

                def support_rows():
                    for key, chrom, pos, ref, alt, _rank in db.execute(
                            'SELECT canonical_key, chrom, pos, ref, alt, contig_order FROM candidates ORDER BY contig_order,pos,ref,alt'):
                        candidate = (chrom, str(pos), ref, alt)
                        results = []
                        for member in members:
                            target_rows = [json.loads(found[0]) for found in db.execute(
                                'SELECT row_json FROM alleles WHERE dataset_id=? AND canonical_key=? ORDER BY ordinal,alt_index',
                                (member['dataset_id'], key))]
                            by_ordinal = {}
                            for row in target_rows:
                                by_ordinal.setdefault(int(row['source_record_index']), []).append(row)
                            occurrences = []
                            for ordinal in sorted(by_ordinal):
                                record_rows = [json.loads(found[0]) for found in db.execute(
                                    'SELECT row_json FROM alleles WHERE dataset_id=? AND ordinal=? ORDER BY alt_index',
                                    (member['dataset_id'], ordinal))]
                                occurrences.append(_record_observation(record_rows, key))
                            results.append(resolve_member(occurrences))
                        signatures = {result['comparison_signature'] for member, result in zip(members, results)
                                      if not result['presence_conflict'] and result['has_called_allele']
                                      and result['comparison_signature'] is not None}
                        genotype_disagreement = (len(signatures) > 1 or any(
                            result['duplicate_genotype_disagreement'] or result['ploidy_disagreement']
                            for result in results))
                        support_row = _support_row(group, candidate, results, members, genotype_disagreement)
                        for member, result in zip(members, results):
                            evidence = _member_evidence_row(group, candidate, member, result)
                            db.execute('INSERT INTO evidence_rows(row_json) VALUES (?)',
                                       (json.dumps(evidence, ensure_ascii=False, separators=(',', ':')),))
                        yield support_row
                        supporters = [quote(member['sample_name'], safe='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-')
                                      for member, result in zip(members, results) if result['support']]
                        conflicts = [quote(member['sample_name'], safe='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-')
                                     for member, result in zip(members, results) if result['presence_conflict']]
                        if int(support_row['support_count']) >= group['required_support']:
                            info = [f"GROUP_N={group['n']}", f"REQUIRED_SUPPORT={group['required_support']}",
                                    f"PROMISING_SUPPORT={support_row['support_count']}",
                                    f"PROMISING_PERCENT={support_row['support_percent']}",
                                    'PROMISING_SAMPLES=' + ','.join(supporters),
                                    'CONFLICTING_SAMPLES=' + (','.join(conflicts) if conflicts else '.')]
                            if genotype_disagreement:
                                info.append('GT_DISAGREEMENT')
                            projected = [result['projected_GT'] for result in results]
                            vcf.write('\t'.join([chrom, str(pos), '.', ref, alt, '.', '.', ';'.join(info),
                                                 'GT', *projected]) + '\n')

                write_tsv(support_path, support_rows(), support_fields)
                write_tsv(evidence_path,
                          (json.loads(row[0]) for row in db.execute('SELECT row_json FROM evidence_rows ORDER BY row_id')),
                          evidence_fields)

            def excluded_rows():
                for member in members:
                    for found in db.execute(
                            "SELECT row_json FROM alleles WHERE dataset_id=? AND status='excluded' ORDER BY ordinal,alt_index",
                            (member['dataset_id'],)):
                        row = json.loads(found[0])
                        raw_fields = json.loads(row['original_record_fields'])
                        alts = [] if raw_fields[4] == '.' else raw_fields[4].split(',')
                        parsed = parse_gt(row['GT'], len(alts))
                        index = int(row['original_alt_index'])
                        if index not in parsed['alleles']:
                            continue
                        yield dict(group_id=group_id, role=group['role'], dataset_id=member['dataset_id'],
                                   sample_name=member['sample_name'], source_vcf=row['source_vcf'],
                                   source_sha256=row['source_sha256'], source_record_index=row['source_record_index'],
                                   record_id=row['record_id'], CHROM=row['CHROM'], POS=row['POS'], ID=row['ID'],
                                   REF=row['REF'], ALT=row['ALT'], GT=row['GT'], original_alt_index=index,
                                   original_alt=row['original_alt'], annotation_status=row['annotation_status'],
                                   window_id=row['window_id'], interval_start=row['interval_start'],
                                   interval_end=row['interval_end'], selection_participant=row['selection_participant'],
                                   evidence_context=row['evidence_context'], exclusion_reason=row['exclusion_reason'])
            write_tsv(excluded_path, excluded_rows(), excluded_fields)
            return dict(group_id=group_id, candidate_count=db.execute('SELECT count(*) FROM candidates').fetchone()[0],
                        member_count=group['n'], required_support=group['required_support'],
                        freeze_sha256=sha256(Path(output_dir) / 'library/freeze.json'),
                        consensus_fraction=f'{fraction.numerator}/{fraction.denominator}')
        finally:
            db.close()


def capture_group(config, output_dir, group_id, raw_vcf=None):
    group, _, _, directory, inputs = group_context(output_dir, group_id)
    raw = Path(raw_vcf) if raw_vcf else directory / 'consensus.raw.vcf'
    vcf = directory / 'consensus.vcf.gz'
    index = Path(str(vcf) + '.csi')
    support = directory / 'support.tsv.gz'
    evidence = directory / 'member_evidence.tsv.gz'
    excluded = directory / 'excluded_alleles.tsv.gz'
    verify_bgzf_csi(vcf, index, raw)
    detail = {
        'group_id': group_id, 'role': group['role'], 'member_dataset_ids': group['member_dataset_ids'],
        'member_count': group['n'], 'required_support': group['required_support'],
        'freeze_sha256': sha256(Path(output_dir) / 'library/freeze.json'),
    }
    settings = config.get('vcf_outputs', {'consensus_fraction': 0.75, 'groups': {}})
    params = {'group_id': group_id, 'config_paths': ['vcf_outputs'],
              'config_values': {'vcf_outputs': settings}}
    outputs = [vcf, support, evidence, excluded]
    return publish_complete(directory, 'complete.json', inputs, params, outputs + [index], detail)


def complete_outputs(config, output_dir):
    output_dir = Path(output_dir)
    datasets_path = output_dir / 'vcf_outputs/datasets.json'
    groups_path = output_dir / 'vcf_outputs/groups.json'
    verify_artifact(datasets_path, [datasets_path])
    verify_artifact(groups_path, [groups_path])
    datasets = load_json(datasets_path)['datasets']
    groups = load_json(groups_path)['groups']
    freeze = verify_freeze(output_dir)
    inputs = [datasets_path, Path(str(datasets_path) + '.provenance.json'),
              groups_path, Path(str(groups_path) + '.provenance.json'),
              output_dir / 'library/freeze.json', Path(freeze['_reference']), Path(freeze['_reference'] + '.fai')]
    for row in datasets:
        directory = output_dir / 'vcf_outputs' / row['dataset_id']
        marker = directory / 'export.complete.json'
        vcf = directory / 'promising.vcf.gz'
        index = Path(str(vcf) + '.csi')
        verify_artifact(marker, [marker, vcf, index])
        record = load_json(marker)
        if record.get('dataset_id') != row['dataset_id'] or record.get('freeze_sha256') != row['freeze_sha256']:
            raise ValueError('Sample export completion differs from manifest: ' + row['dataset_id'])
        inputs += [marker, Path(str(marker) + '.provenance.json'), vcf, index]
    for group in groups:
        group_context(output_dir, group['group_id'])
        directory = output_dir / 'vcf_outputs/groups' / group['group_id']
        marker = directory / 'complete.json'
        vcf = directory / 'consensus.vcf.gz'
        index = Path(str(vcf) + '.csi')
        support, evidence, excluded = (directory / name for name in
                                       ('support.tsv.gz', 'member_evidence.tsv.gz', 'excluded_alleles.tsv.gz'))
        verify_artifact(marker, [marker, vcf, index, support, evidence, excluded])
        record = load_json(marker)
        if (record.get('group_id') != group['group_id'] or
                record.get('member_dataset_ids') != group['member_dataset_ids'] or
                record.get('member_count') != group['n'] or
                record.get('required_support') != group['required_support'] or
                record.get('freeze_sha256') != sha256(output_dir / 'library/freeze.json')):
            raise ValueError('Group completion differs from validated membership: ' + group['group_id'])
        inputs += [marker, Path(str(marker) + '.provenance.json'), vcf, index, support, evidence, excluded]
    settings = config.get('vcf_outputs', {'consensus_fraction': 0.75, 'groups': {}})
    return publish_complete(output_dir / 'vcf_outputs', 'outputs.complete.json', inputs,
                            {'config_paths': ['vcf_outputs'], 'config_values': {'vcf_outputs': settings}},
                            [], {'dataset_count': len(datasets), 'group_count': len(groups),
                                 'freeze_sha256': sha256(output_dir / 'library/freeze.json')})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['manifest', 'sample-inputs', 'capture-alleles', 'capture-export',
                                          'consensus', 'capture-group', 'complete'])
    parser.add_argument('--config', required=True)
    parser.add_argument('--dataset-id')
    parser.add_argument('--group-id')
    parser.add_argument('--normalized-vcf')
    parser.add_argument('--raw-vcf')
    parser.add_argument('--bgzip-version')
    parser.add_argument('--disk-mb', type=int)
    args = parser.parse_args()
    if args.action in ('sample-inputs', 'capture-alleles', 'capture-export') and not args.dataset_id:
        parser.error(args.action + ' requires --dataset-id')
    if args.action in ('consensus', 'capture-group') and not args.group_id:
        parser.error(args.action + ' requires --group-id')
    config = load_json(args.config)
    identity = args.dataset_id or args.group_id or 'global'
    default_disk = {'sample-inputs': 30000, 'capture-alleles': 1000, 'capture-export': 1000,
                    'consensus': 60000, 'capture-group': 1000}.get(args.action, 1000)
    configure_scratch(config, args.action, identity, disk_mb=args.disk_mb or default_disk, scope='application')
    output_dir = Path(config['output_dir'])
    if args.action == 'manifest':
        manifest(config, output_dir)
    elif args.action == 'sample-inputs':
        sample_inputs(config, output_dir, args.dataset_id)
    elif args.action == 'capture-alleles':
        capture_alleles(config, output_dir, args.dataset_id, args.normalized_vcf)
    elif args.action == 'capture-export':
        capture_export(config, output_dir, args.dataset_id, args.bgzip_version)
    elif args.action == 'consensus':
        write_consensus(config, output_dir, args.group_id)
    elif args.action == 'capture-group':
        capture_group(config, output_dir, args.group_id, args.raw_vcf)
    else:
        complete_outputs(config, output_dir)


if __name__ == '__main__':
    main()
