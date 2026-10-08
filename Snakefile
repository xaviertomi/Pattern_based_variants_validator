"""Standalone reference-context motif selection and frozen-library application."""
import sys
from pathlib import Path

ROOT = Path(workflow.basedir).resolve()
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'profile'))
from prepare import validate_config, dataset_id, ALT_FILTER_EXPR
from workflow import load_yaml

# No input inspection or writes occur while constructing the rule graph.
DEFAULT_CONFIG = load_yaml(ROOT / 'config.yaml')
for config_path in workflow.configfiles:
    load_yaml(config_path)
def merge_config(a, b):
    for key, value in b.items():
        if isinstance(value, dict) and isinstance(a.get(key), dict):
            merge_config(a[key], value)
        else:
            a[key] = value
    return a
config = validate_config(merge_config(DEFAULT_CONFIG, config), ROOT)
OUT = Path(config['output_dir'])
CFG = str(OUT / 'manifests' / 'effective_config.json')
CONFIG_JSON = CFG
CFG_JSON = CFG
ROTATIONS = [1, 2, 3, 4, 5]
PB_IDS = [dataset_id('pacbio', name) for name in config['pacbio_samples']]
IL_IDS = [dataset_id('illumina', name) for pb in config['pacbio_samples'] for name in config['illumina_corresponding'][pb]]
DATASETS = {}
for name in config['pacbio_samples']:
    DATASETS[dataset_id('pacbio', name)] = {'role': 'pacbio', 'sample_name': name,
        'source_vcf': config['pacbio_vcf'], 'corresponding_pacbio': name}
    for il in config['illumina_corresponding'][name]:
        DATASETS[dataset_id('illumina', il)] = {'role': 'illumina', 'sample_name': il,
            'source_vcf': config['illumina_vcf'], 'corresponding_pacbio': name}

def IMAGE(tool):
    return config['container'][tool] if config['execution']['use_containers'] else None

shell.executable('/bin/bash')
shell.prefix('set -euo pipefail; ')

rule all:
    input:
        str(OUT / 'library' / 'freeze.json'),
        str(OUT / 'application' / 'genotypes.tsv.gz'),
        str(OUT / 'reports' / 'report.complete.json'),
        str(OUT / 'vcf_outputs' / 'outputs.complete.json')
rule prepare:
    input: str(OUT / 'manifests' / 'preparation.ok.json')

rule preflight:
    input: str(OUT / 'manifests' / 'preflight.ok.json')

rule discover:
    input: [str(OUT / 'evidence' / pb / f'r{r:02d}' / 'candidates.json') for pb in PB_IDS for r in ROTATIONS]

rule select_corresponding:
    input: [str(OUT / 'evidence' / pb / f'r{r:02d}' / 'transfer.complete.json') for pb in PB_IDS for r in ROTATIONS]

rule freeze_library:
    input: str(OUT / 'library' / 'freeze.json')

rule annotate_other:
    input: str(OUT / 'application' / 'genotypes.tsv.gz')
rule annotate_corresponding:
    input: str(OUT / 'application' / 'corresponding' / 'genotypes.tsv.gz')

rule reports:
    input: str(OUT / 'reports' / 'report.complete.json')

def pilot_input(wildcards):
    pilot = config.get('pilot')
    if not pilot:
        raise ValueError('Pilot requires explicit pilot.pacbio_sample and pilot.rotation settings.')
    return str(OUT / 'evidence' / dataset_id('pacbio', pilot['pacbio_sample']) / f"r{pilot['rotation']:02d}" / 'transfer.complete.json')

rule pilot:
    input: pilot_input

include: 'rules/prepare.smk'
if config.get('workflow_scope', 'selection') != 'application':
    include: 'rules/discovery.smk'
    include: 'rules/transfer.smk'
    include: 'rules/selection.smk'
include: 'rules/application.smk'


rule export_promising:
    input:
        datasets=lambda wc: checkpoints.plan_vcf_outputs.get().output.datasets,
        exports=vcf_export_completions


rule consensus_groups:
    input:
        groups=lambda wc: checkpoints.plan_vcf_outputs.get().output.groups,
        groups_provenance=lambda wc: checkpoints.plan_vcf_outputs.get().output.groups_provenance,
        completions=vcf_group_completions


rule vcf_deliverables:
    input: str(OUT / 'vcf_outputs' / 'outputs.complete.json')
