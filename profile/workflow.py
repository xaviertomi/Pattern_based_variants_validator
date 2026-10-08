"""Host-side configuration, deployment and checksum inventory for the launcher."""
import argparse
import copy
import hashlib
import importlib.metadata
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]
class UniqueLoader(yaml.SafeLoader):
    pass

def mapping(loader, node, deep=False):
    result = {}
    for keynode, valnode in node.value:
        key = loader.construct_object(keynode, deep=deep)
        if key in result:
            raise ValueError(f'Duplicate YAML key {key!r} at line {keynode.start_mark.line + 1}')
        result[key] = loader.construct_object(valnode, deep=deep)
    return result
UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, mapping)

def load_yaml(path):
    with open(path, encoding='utf-8') as handle:
        value = yaml.load(handle, Loader=UniqueLoader)
    if not isinstance(value, dict):
        raise ValueError(f'Configuration must be a mapping: {path}')
    return value

def merge(a, b):
    for key, value in b.items():
        if isinstance(value, dict) and isinstance(a.get(key), dict):
            merge(a[key], value)
        else:
            a[key] = copy.deepcopy(value)
    return a

def split_args(args):
    files, values, rest = [], [], []
    i = 0
    while i < len(args):
        token = args[i]
        if token in ('--configfile', '--configfiles', '--config'):
            kind = token
            i += 1
            collected = []
            while i < len(args) and not args[i].startswith('-'):
                collected.append(args[i]); i += 1
            if not collected:
                raise ValueError(f'{kind} requires values')
            (values if kind == '--config' else files).extend(collected)
        elif token.startswith('--configfile='):
            files.append(token.split('=', 1)[1]); i += 1
        else:
            rest.append(token); i += 1
    return files, values, rest

def effective(args):
    files, values, rest = split_args(args)
    cfg = load_yaml(ROOT / 'config.yaml')
    files = [str((ROOT / name).resolve()) for name in files]
    for path in files:
        merge(cfg, load_yaml(path))
    for value in values:
        if '=' not in value:
            raise ValueError(f'Expected --config key=value: {value}')
        key, text = value.split('=', 1)
        cfg[key] = yaml.safe_load(text)
    sys.path.insert(0, str(ROOT / 'scripts'))
    from prepare import validate_config
    cfg = validate_config(cfg, ROOT)
    forbidden = (ROOT.parent / 'clone_vcf_filter').resolve()
    writable = {'logs': ROOT / 'logs', 'work': ROOT / '.snakemake',
                'output_dir': cfg['output_dir'], 'scratch_dir': cfg['scratch_dir'],
                'image_cache': cfg['execution']['image_cache']}
    for index, item in enumerate(cfg['execution'].get('extra_bind_paths', [])):
        if isinstance(item, dict) and item.get('mode', 'ro') == 'rw':
            writable[f'extra_bind_paths[{index}]'] = item['path']
    for key, path in writable.items():
        resolved = Path(path).resolve()
        if resolved == forbidden or forbidden in resolved.parents:
            raise ValueError(f'Writable path inside read-only infrastructure reference: {key}: {resolved}')
    overrides = []
    if files:
        overrides += ['--configfile', *files]
    if values:
        overrides += ['--config', *values]
    return cfg, overrides, rest

def execute(args, capture=False):
    p = subprocess.run(args, cwd=ROOT, text=True, capture_output=capture)
    if p.returncode:
        if capture and p.stderr:
            print(p.stderr, file=sys.stderr)
        raise RuntimeError(f'{args[0]} exited {p.returncode}')
    return p.stdout if capture else None

def host_check(cluster=False):
    version = execute(['snakemake', '--version'], True).strip()
    if version != '9.4.0':
        raise RuntimeError(f'Required Snakemake 9.4.0; found {version}')
    if cluster and importlib.metadata.version('snakemake-executor-plugin-cluster-generic') != '1.0.9':
        raise RuntimeError('Required cluster-generic plugin 1.0.9')

def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def inventory(base, target, scope, cfg):
    """Read only declared outputs, then invalidate current signatures and descendants."""
    summary = execute(['snakemake', target, *base, '--profile', 'none', '--executor', 'local', '--cores', '1', '--summary'], True)
    outputs = []
    for line in summary.splitlines():
        parts = line.split('\t')
        if len(parts) >= 5 and parts[0] != 'output_file':
            outputs.append(Path(parts[0]))
    invalid = set()
    if scope == 'selection':
        out = Path(cfg['output_dir'])
        fai = out / 'reference' / 'reference.fa.fai'
        folds = out / 'manifests' / 'folds.json'
        if fai.exists() and folds.exists():
            from prepare import verify_folds
            try:
                verify_folds(fai, out, digest(out / 'reference' / 'reference.fa'), cfg['folds']['grouping'])
            except (OSError, ValueError, RuntimeError, KeyError):
                invalid.add(folds)
    dependencies = {}
    artifacts = {}
    for output in outputs:
        if str(output).endswith('.provenance.json'):
            primary = Path(str(output)[:-len('.provenance.json')])
            try:
                record = json.loads(output.read_text())
                if record.get('scope') != scope:
                    continue
                paths = record['input_paths']
                dependencies[primary] = {Path(p) for p in paths}
                params = copy.deepcopy(record['signature_params'])
                if 'code_paths' in params:
                    params['code_hashes'] = {p: digest(p) for p in params['code_paths']}
                if 'config_paths' in params:
                    params['config_values'] = {}
                    for key in params['config_paths']:
                        value = cfg
                        for component in key.split('.'):
                            value = value[component]
                        params['config_values'][key] = value
                signature = hashlib.sha256(canonical({'inputs': {p: digest(p) for p in paths}, 'params': params})).hexdigest()
                required = record['required_outputs']
                for path in required:
                    artifacts[Path(path)] = primary
                from prepare import check_reuse
                if not required or not check_reuse(output, signature, [Path(p) for p in required]):
                    invalid.add(primary)
                    invalid.update(Path(p) for p, h in required.items()
                                   if not Path(p).is_file() or digest(p) != h)
            except (OSError, KeyError, ValueError, TypeError):
                invalid.add(primary)
    # Missing declared scientific outputs cannot be accepted based on newer mtimes.
    invalid.update(p for p in outputs if not p.exists() and not str(p).endswith('.provenance.json'))
    changed = True
    while changed:
        changed = False
        for primary, inputs in dependencies.items():
            if primary not in invalid and any(p in invalid or artifacts.get(p) in invalid for p in inputs):
                invalid.add(primary); changed = True
    return sorted(str(p) for p in invalid)

def binds(cfg, scope, *, source_keys=()):
    mounts = {}
    def add(path, mode):
        p = str(Path(path).absolute())
        if p in mounts and mounts[p] != mode:
            raise ValueError(f'Conflicting container bind modes: {p}')
        mounts[p] = mode
    add(ROOT, 'ro')
    for key in ('output_dir', 'scratch_dir'):
        add(cfg[key], 'rw')
    add(ROOT / 'logs', 'rw'); add(ROOT / '.snakemake', 'rw')
    add(cfg['execution']['image_cache'], 'rw')
    if scope == 'selection':
        for key in ('reference_fasta', 'pacbio_vcf', 'illumina_vcf'):
            path = Path(cfg[key])
            add(path.absolute(), 'ro'); add(path.resolve(), 'ro')
    elif scope == 'application':
        allowed = {'illumina_vcf', 'other_illumina'}
        if not source_keys:
            raise ValueError('Application scope requires target-selected source keys')
        unknown = set(source_keys) - allowed
        if unknown:
            raise ValueError('Unsupported application source bind: ' + repr(sorted(unknown)))
        for key in source_keys:
            add(Path(cfg[key]).absolute(), 'ro')
            add(Path(cfg[key]).resolve(), 'ro')
    for item in cfg['execution'].get('extra_bind_paths', []):
        if isinstance(item, dict):
            add(item['path'], item.get('mode', 'ro'))
        else:
            add(item, 'ro')
    return ','.join(f'{p}:{p}:{mode}' for p, mode in sorted(mounts.items()))

def deployment(cfg, scope, provision=False, *, source_keys=()):
    if not cfg['execution']['use_containers']:
        return []
    cache = Path(cfg['execution']['image_cache'])
    scratch = Path(cfg['scratch_dir'])
    for p in (cache, cache / 'oci', scratch / 'apptainer'):
        p.mkdir(parents=True, exist_ok=True)
    os.environ['APPTAINER_CACHEDIR'] = str(cache / 'oci')
    os.environ['APPTAINER_TMPDIR'] = str(scratch / 'apptainer')
    if not shutil.which('apptainer'):
        raise RuntimeError('Containers require Apptainer on PATH')
    if provision:
        # Provision immutable images serially before any scientific workers are submitted.
        for image in sorted(set(cfg['container'].values())):
            if image.startswith('docker://'):
                name = hashlib.md5(image.encode()).hexdigest() + '.simg'
                destination = cache / name
                if not destination.exists():
                    temporary = cache / (name + '.partial')
                    execute(['apptainer', 'pull', '--force', str(temporary), image])
                    temporary.replace(destination)
            elif not Path(image).is_file():
                raise RuntimeError(f'Pinned image is unavailable: {image}')
    return ['--sdm', 'apptainer', '--apptainer-prefix', str(cache),
            '--apptainer-args', '--bind ' + binds(cfg, scope, source_keys=source_keys)]

def frozen_complete(out):
    from annotate import verify_freeze
    try:
        verify_freeze(out)
    except ValueError as error:
        raise RuntimeError(str(error)) from error

def run(mode, args):
    from snakemake.cli import get_argument_parser
    cfg, overrides, rest = effective(args)
    targets = get_argument_parser().parse_args(rest).targets
    if not targets:
        rest.insert(0, 'pilot' if mode == 'pilot' else 'all')
        targets = [rest[0]]
    base = ['--snakefile', str(ROOT / 'Snakefile'), '--directory', str(ROOT), *overrides]
    host_check(mode not in ('dry', 'dag', 'unlock'))
    application_targets = {'annotate_other', 'annotate_corresponding', 'export_promising',
                           'consensus_groups', 'vcf_deliverables'}
    if mode in ('dry', 'dag', 'unlock'):
        if targets and all(target in application_targets for target in targets):
            if '--config' not in base:
                base.append('--config')
            base.append('workflow_scope=application')
        return
    if not os.environ.get('SLURM_JOB_ID') and mode != 'controller':
        c = load_yaml(ROOT / 'profile' / 'controller.yaml')
        (ROOT / 'logs' / 'controller').mkdir(parents=True, exist_ok=True)
        cmd = ['sbatch', '--parsable', '--export=ALL', '--account', str(c['account']), '--partition', str(c['partition']),
               '--cpus-per-task', str(c['cpus']), '--mem', str(c['mem_mb']), '--time', str(c['time']),
               '--output', str(ROOT / 'logs/controller/motif-%j.out'), '--error', str(ROOT / 'logs/controller/motif-%j.err')]
        if c.get('mail_user'):
            cmd += ['--mail-user', str(c['mail_user']), '--mail-type', ','.join(c['mail_type'])]
        cmd += [str(ROOT / 'snakemake.sh'), 'controller', mode, *args]
        response = execute(cmd, True).strip().split(';')
        if not response[0].isdigit() or len(response) > 2:
            raise RuntimeError('Invalid controller sbatch response')
        from slurm import save
        save(response[0], {'cluster': response[1] if len(response) == 2 else None})
        print(response[0])
        return
    if mode == 'controller':
        raise ValueError('Internal controller mode must receive original mode')
    (ROOT / 'logs' / 'jobs').mkdir(parents=True, exist_ok=True)
    if targets == ['preflight']:
        execute(['snakemake', *rest, *base, '--profile', str(ROOT / 'profile'), *deployment(cfg, 'preflight', True)])
        return
    selection = any(t not in application_targets for t in targets)
    combined = any(t in ('all', 'reports') for t in targets)
    def production(scope, target, requested):
        scope_base = list(base)
        if scope == 'application':
            if '--config' not in scope_base:
                scope_base.append('--config')
            scope_base.append('workflow_scope=application')
        fingerprints = ['fingerprint_selection'] if scope == 'selection' else []
        source_keys = set()
        if scope == 'application':
            external = {'annotate_other', 'reports', 'all', 'export_promising',
                        'consensus_groups', 'vcf_deliverables'}
            corresponding = {'annotate_corresponding', 'all', 'export_promising',
                             'consensus_groups', 'vcf_deliverables'}
            requested_targets = set(targets)
            if requested_targets & external:
                fingerprints.append('fingerprint_application')
                source_keys.add('other_illumina')
            if requested_targets & corresponding:
                fingerprints.append('fingerprint_corresponding')
                source_keys.add('illumina_vcf')
        for fingerprint in fingerprints:
            execute(['snakemake', fingerprint, *scope_base, '--profile', 'none', '--executor', 'local',
                     '--cores', '1', '--forcerun', fingerprint])
        force = inventory(scope_base, target, scope, cfg)
        # Checksums drive scientific invalidation; metadata timestamps do not.
        cmd = ['snakemake', *requested, *scope_base, '--profile', str(ROOT / 'profile'),
               '--rerun-triggers', 'input', *deployment(cfg, scope, True, source_keys=source_keys)]
        if mode == 'rerun':
            cmd += ['--forceall']
        elif force:
            cmd += ['--forcerun', *force]
        execute(cmd)
    if selection:
        production('selection', 'freeze_library' if combined else targets[0], ['freeze_library'] if combined else rest)
    if combined or not selection:
        frozen_complete(Path(cfg['output_dir']))
        production('application', targets[0], rest)

def main():
    if len(sys.argv) < 2:
        raise ValueError('Expected launcher mode')
    mode, args = sys.argv[1], sys.argv[2:]
    if mode == 'controller':
        if not os.environ.get('SLURM_JOB_ID'):
            raise ValueError('Internal controller mode requires an existing SLURM allocation')
        if not args:
            raise ValueError('Missing controller execution mode')
        mode, args = args[0], args[1:]
    run(mode, args)

if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, RuntimeError, importlib.metadata.PackageNotFoundError) as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)
