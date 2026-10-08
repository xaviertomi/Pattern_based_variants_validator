"""External application with real single-image tool boundaries and dynamic shards."""
import json

APP_CFG = str(OUT / 'application' / 'effective_config.json')
CORR_CFG = str(OUT / 'application/corresponding/effective_config.json')

localrules: fingerprint_application, fingerprint_corresponding

rule fingerprint_application:
    input: freeze=str(OUT / 'library/freeze.json')
    output:
        fingerprint=str(OUT / 'application/input_fingerprints.json'), config=APP_CFG
    run:
        import sys
        sys.path.insert(0, str(ROOT / 'scripts'))
        from annotate import application_fingerprint
        from prepare import atomic_json
        atomic_json(output.config, dict(config))
        application_fingerprint(dict(config), OUT)



def other_manifest(wildcards):
    with open(checkpoints.prepare_other_manifest.get().output.manifest, encoding='utf-8') as handle:
        return json.load(handle)


def other_prepared(wildcards):
    return [str(OUT / 'prepared' / row['dataset_id'] / 'sample.complete.json') for row in other_manifest(wildcards)['datasets']]



def corresponding_inputs(wildcards=None):
    paths = [OUT / 'library/freeze.json', OUT / 'manifests/selection_fingerprints.json',
             OUT / 'manifests/selection_config.json', Path(config['illumina_vcf']),
             OUT / 'reference/reference.fa', OUT / 'reference/reference.fa.fai']
    for name in (sample for members in config['illumina_corresponding'].values() for sample in members):
        did = dataset_id('illumina', name)
        directory = OUT / 'prepared' / did
        paths.extend(directory / file for file in (
            'records.complete.json', 'records.complete.json.provenance.json', 'records.tsv.gz',
            'sample.complete.json', 'sample.complete.json.provenance.json',
            'genotypes.tsv.gz', 'all.fasta', 'windows.tsv'))
    return list(map(str, paths))


def corresponding_manifest_rows(wildcards):
    with open(checkpoints.prepare_corresponding_manifest.get().output.manifest, encoding='utf-8') as handle:
        return json.load(handle)


def corresponding_prepared(wildcards):
    return [str(OUT / 'prepared' / row['dataset_id'] / 'sample.complete.json')
            for row in corresponding_manifest_rows(wildcards)['datasets']]


def corresponding_shard_manifest(wildcards):
    with open(checkpoints.shard_corresponding_windows.get().output.manifest, encoding='utf-8') as handle:
        return json.load(handle)


def corresponding_scans(wildcards):
    return [str(OUT / 'application/corresponding/scans' / (row['shard_id'] + '.tsv'))
            for row in corresponding_shard_manifest(wildcards)['shards']]


def corresponding_success(wildcards):
    return [path + '.success.json' for path in corresponding_scans(wildcards)]


def corresponding_fasta(wildcards):
    matches = [row['fasta'] for row in corresponding_shard_manifest(wildcards)['shards']
               if row['shard_id'] == wildcards.shard]
    if len(matches) != 1:
        raise ValueError('Unknown corresponding scan shard: ' + wildcards.shard)
    return matches[0]


def corresponding_fastas(wildcards):
    return [str(OUT / 'prepared' / row['dataset_id'] / 'all.fasta')
            for row in corresponding_manifest_rows(wildcards)['datasets']]

rule fingerprint_corresponding:
    input: artifacts=corresponding_inputs
    output:
        fingerprint=str(OUT / 'application/corresponding/input_fingerprints.json'),
        config=CORR_CFG
    run:
        import sys
        sys.path.insert(0, str(ROOT / 'scripts'))
        from annotate import corresponding_fingerprint
        from prepare import atomic_json
        atomic_json(output.config, dict(config))
        corresponding_fingerprint(dict(config), OUT)

def other_shard_manifest(wildcards):
    with open(checkpoints.shard_other_windows.get().output.manifest, encoding='utf-8') as handle:
        return json.load(handle)


def other_scans(wildcards):
    return [str(OUT / 'application/scans' / (row['shard_id'] + '.tsv')) for row in other_shard_manifest(wildcards)['shards']]


def other_success(wildcards):
    return [path + '.success.json' for path in other_scans(wildcards)]


def other_fasta(wildcards):
    matches = [row['fasta'] for row in other_shard_manifest(wildcards)['shards'] if row['shard_id'] == wildcards.shard]
    if len(matches) != 1:
        raise ValueError('Unknown external scan shard: ' + wildcards.shard)
    return matches[0]


rule stage_other_vcf:
    input:
        preflight=str(OUT / 'manifests/preflight.ok.json'),
        fingerprint=str(OUT / 'application/input_fingerprints.json'),
        freeze=str(OUT / 'library/freeze.json'),
        source=lambda wc: config['other_illumina']
    output:
        vcf=str(OUT / 'inputs/other_illumina.vcf.gz'),
        index=str(OUT / 'inputs/other_illumina.vcf.gz.csi'),
        samples=str(OUT / 'application/header.samples.txt'),
        header=str(OUT / 'application/header.vcf.txt')
    log: str(OUT / 'logs/stage_other_vcf.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='02:00:00', disk_mb=30000
    container: IMAGE('bcftools')
    shell:
        'bcftools view -Oz -o {output.vcf:q} {input.source:q} 2> {log:q}; '
        'bcftools index --csi --force {output.vcf:q} 2>> {log:q}; '
        'bcftools query -l {output.vcf:q} > {output.samples:q} 2>> {log:q}; '
        'bcftools view -h {output.vcf:q} > {output.header:q} 2>> {log:q}'


checkpoint prepare_other_manifest:
    input:
        vcf=str(OUT / 'inputs/other_illumina.vcf.gz'),
        index=str(OUT / 'inputs/other_illumina.vcf.gz.csi'),
        samples=str(OUT / 'application/header.samples.txt'),
        header=str(OUT / 'application/header.vcf.txt'),
        freeze=str(OUT / 'library/freeze.json')
    output:
        manifest=str(OUT / 'application/manifest.json'),
        datasets=str(OUT / 'application/datasets.tsv'),
        configuration=str(OUT / 'application/application_config.json'),
        provenance=str(OUT / 'application/manifest.json.provenance.json')
    params: script=str(ROOT / 'scripts/annotate.py'), config=APP_CFG
    log: str(OUT / 'logs/prepare_other_manifest.stderr.log')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    container: IMAGE('python')
    shell:
        'python3 {params.script:q} manifest --prepared-headers --config {params.config:q} --disk-mb {resources.disk_mb} 2> {log:q}'


rule plan_other_records:
    input: manifest=str(OUT / 'application/manifest.json')
    output:
        plan=str(OUT / 'prepared/{other_dataset}/records.plan.json'),
        sample=str(OUT / 'prepared/{other_dataset}/exact_sample.txt'),
        provenance=str(OUT / 'prepared/{other_dataset}/records.plan.json.provenance.json')
    wildcard_constraints: other_dataset='other_illumina_[a-f0-9]{16}'
    params: script=str(ROOT / 'scripts/annotate.py'), config=APP_CFG
    log: str(OUT / 'logs/plan_other_records/{other_dataset}.stderr.log')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    container: IMAGE('python')
    shell:
        'python3 {params.script:q} prepare-sample --phase records-plan --config {params.config:q} --dataset-id {wildcards.other_dataset:q} --disk-mb {resources.disk_mb} 2> {log:q}'


rule query_other_records:
    input:
        plan=str(OUT / 'prepared/{other_dataset}/records.plan.json'),
        sample=str(OUT / 'prepared/{other_dataset}/exact_sample.txt'),
        vcf=str(OUT / 'inputs/other_illumina.vcf.gz'),
        reference=str(OUT / 'reference/reference.fa')
    output:
        sample=str(OUT / 'prepared/{other_dataset}/sample.vcf.gz'),
        sample_index=str(OUT / 'prepared/{other_dataset}/sample.vcf.gz.csi'),
        alt=str(OUT / 'prepared/{other_dataset}/sample.alt.vcf.gz'),
        alt_index=str(OUT / 'prepared/{other_dataset}/sample.alt.vcf.gz.csi'),
        records=str(OUT / 'prepared/{other_dataset}/raw.records.tsv'),
        alt_records=str(OUT / 'prepared/{other_dataset}/raw.alt.tsv')
    wildcard_constraints: other_dataset='other_illumina_[a-f0-9]{16}'
    params: fmt='%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n', alt_filter=ALT_FILTER_EXPR
    log: str(OUT / 'logs/query_other_records/{other_dataset}.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='02:00:00', disk_mb=30000
    container: IMAGE('bcftools')
    shell:
        'bcftools view -S {input.sample:q} -Oz -o {output.sample:q} {input.vcf:q} 2> {log:q}; '
        'bcftools index --csi --force {output.sample:q} 2>> {log:q}; '
        'bcftools query -f {params.fmt:q} {output.sample:q} > {output.records:q} 2>> {log:q}; '
        'bcftools view -i {params.alt_filter:q} -Oz -o {output.alt:q} {output.sample:q} 2>> {log:q}; '
        'bcftools index --csi --force {output.alt:q} 2>> {log:q}; '
        'bcftools norm -f {input.reference:q} -c e -Ou {output.alt:q} > /dev/null 2>> {log:q}; '
        'bcftools query -f {params.fmt:q} {output.alt:q} > {output.alt_records:q} 2>> {log:q}'


rule prepare_other_records:
    input:
        manifest=str(OUT / 'application/manifest.json'),
        records=str(OUT / 'prepared/{other_dataset}/raw.records.tsv'),
        alt=str(OUT / 'prepared/{other_dataset}/raw.alt.tsv')
    output:
        complete=str(OUT / 'prepared/{other_dataset}/records.complete.json'),
        records=str(OUT / 'prepared/{other_dataset}/records.tsv.gz'),
        provenance=str(OUT / 'prepared/{other_dataset}/records.complete.json.provenance.json')
    wildcard_constraints: other_dataset='other_illumina_[a-f0-9]{16}'
    params: script=str(ROOT / 'scripts/annotate.py'), config=APP_CFG
    log: str(OUT / 'logs/prepare_other_records/{other_dataset}.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='02:00:00', disk_mb=30000
    container: IMAGE('python')
    shell:
        'python3 {params.script:q} prepare-sample --phase records-from-query --config {params.config:q} --dataset-id {wildcards.other_dataset:q} --disk-mb {resources.disk_mb} 2> {log:q}'


rule plan_other_contexts:
    input: complete=str(OUT / 'prepared/{other_dataset}/records.complete.json')
    output:
        regions=str(OUT / 'prepared/{other_dataset}/regions.txt'),
        plan=str(OUT / 'prepared/{other_dataset}/contexts.plan.json'),
        provenance=str(OUT / 'prepared/{other_dataset}/contexts.plan.json.provenance.json')
    wildcard_constraints: other_dataset='other_illumina_[a-f0-9]{16}'
    params: script=str(ROOT / 'scripts/annotate.py'), config=APP_CFG
    log: str(OUT / 'logs/plan_other_contexts/{other_dataset}.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=10000
    container: IMAGE('python')
    shell:
        'python3 {params.script:q} prepare-sample --phase contexts-plan --config {params.config:q} --dataset-id {wildcards.other_dataset:q} --disk-mb {resources.disk_mb} 2> {log:q}'


rule extract_other_contexts:
    input:
        regions=str(OUT / 'prepared/{other_dataset}/regions.txt'),
        reference=str(OUT / 'reference/reference.fa'),
        index=str(OUT / 'reference/reference.fa.fai')
    output: fasta=str(OUT / 'prepared/{other_dataset}/extracted.fasta')
    wildcard_constraints: other_dataset='other_illumina_[a-f0-9]{16}'
    log: str(OUT / 'logs/extract_other_contexts/{other_dataset}.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=10000
    container: IMAGE('samtools')
    shell:
        'if test -s {input.regions:q}; then samtools faidx -r {input.regions:q} {input.reference:q} > {output.fasta:q} 2> {log:q}; else : > {output.fasta:q}; : > {log:q}; fi'


rule prepare_other_sample:
    input:
        manifest=str(OUT / 'application/manifest.json'),
        plan=str(OUT / 'prepared/{other_dataset}/contexts.plan.json'),
        extracted=str(OUT / 'prepared/{other_dataset}/extracted.fasta')
    output:
        complete=str(OUT / 'prepared/{other_dataset}/sample.complete.json'),
        genotypes=str(OUT / 'prepared/{other_dataset}/genotypes.tsv.gz'),
        windows=str(OUT / 'prepared/{other_dataset}/windows.tsv'),
        fasta=str(OUT / 'prepared/{other_dataset}/all.fasta'),
        provenance=str(OUT / 'prepared/{other_dataset}/sample.complete.json.provenance.json')
    wildcard_constraints: other_dataset='other_illumina_[a-f0-9]{16}'
    params: script=str(ROOT / 'scripts/annotate.py'), config=APP_CFG
    log: str(OUT / 'logs/prepare_other_sample/{other_dataset}.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=10000
    container: IMAGE('python')
    shell:
        'python3 {params.script:q} prepare-sample --phase contexts-from-fasta --config {params.config:q} --dataset-id {wildcards.other_dataset:q} --disk-mb {resources.disk_mb} 2> {log:q}'


checkpoint shard_other_windows:
    input: samples=other_prepared, manifest=str(OUT / 'application/manifest.json')
    output:
        manifest=str(OUT / 'application/shards/manifest.json'),
        windows=str(OUT / 'application/shards/windows.sqlite'),
        provenance=str(OUT / 'application/shards/manifest.json.provenance.json')
    params: script=str(ROOT / 'scripts/annotate.py'), config=APP_CFG
    log: str(OUT / 'logs/shard_other_windows.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=10000
    container: IMAGE('python')
    shell:
        'python3 {params.script:q} shard --role other_illumina --config {params.config:q} --disk-mb {resources.disk_mb} 2> {log:q}'


rule fimo_other:
    input:
        fasta=other_fasta,
        models=str(OUT / 'library/final_motifs.meme'),
        background=str(OUT / 'library/scoring_background.bfile'),
        freeze=str(OUT / 'library/freeze.json')
    output:
        scan=str(OUT / 'application/scans/{shard}.tsv'),
        success=str(OUT / 'application/scans/{shard}.tsv.success.json'),
        provenance=str(OUT / 'application/scans/{shard}.tsv.provenance.json')
    wildcard_constraints: shard='s[0-9]{6}'
    params: script=str(ROOT / 'scripts/annotate.py'), config=APP_CFG
    log: str(OUT / 'logs/fimo_other/{shard}.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='04:00:00', disk_mb=30000
    container: IMAGE('meme')
    shell:
        'python3 {params.script:q} scan --config {params.config:q} --disk-mb {resources.disk_mb} --fasta {input.fasta:q} --output {output.scan:q} --log {log:q}'


rule annotate_genotypes:
    input:
        scans=other_scans, success=other_success, samples=other_prepared,
        freeze=str(OUT / 'library/freeze.json'), shards=str(OUT / 'application/shards/manifest.json')
    output:
        complete=str(OUT / 'application/annotation.complete.json'),
        genotypes=str(OUT / 'application/genotypes.tsv.gz'),
        alt=str(OUT / 'application/alt_annotations.tsv.gz'),
        matches=str(OUT / 'application/all_matches.tsv.gz'),
        counts=str(OUT / 'reports/annotation_counts.tsv'),
        positions=str(OUT / 'reports/application_positions.tsv'),
        provenance=str(OUT / 'application/annotation.complete.json.provenance.json')
    params: script=str(ROOT / 'scripts/annotate.py'), config=APP_CFG
    log: str(OUT / 'logs/annotate_genotypes.stderr.log')
    threads: 1
    resources: mem_mb=8000, time='04:00:00', disk_mb=60000
    container: IMAGE('python')
    shell:
        'python3 {params.script:q} annotate --role other_illumina --config {params.config:q} --disk-mb {resources.disk_mb} --scans {input.scans:q} 2> {log:q}'


checkpoint prepare_corresponding_manifest:
    input:
        fingerprint=str(OUT / 'application/corresponding/input_fingerprints.json'),
        artifacts=corresponding_inputs
    output:
        manifest=str(OUT / 'application/corresponding/manifest.json'),
        datasets=str(OUT / 'application/corresponding/datasets.tsv'),
        provenance=str(OUT / 'application/corresponding/manifest.json.provenance.json')
    params: script=str(ROOT / 'scripts/annotate.py'), config=CORR_CFG
    log: str(OUT / 'logs/prepare_corresponding_manifest.stderr.log')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    container: IMAGE('python')
    shell:
        'python3 {params.script:q} corresponding-manifest --config {params.config:q} --disk-mb {resources.disk_mb} 2> {log:q}'


checkpoint shard_corresponding_windows:
    input:
        samples=corresponding_prepared, fasta=corresponding_fastas,
        manifest=str(OUT / 'application/corresponding/manifest.json')
    output:
        manifest=str(OUT / 'application/corresponding/shards/manifest.json'),
        windows=str(OUT / 'application/corresponding/shards/windows.sqlite'),
        provenance=str(OUT / 'application/corresponding/shards/manifest.json.provenance.json')
    params: script=str(ROOT / 'scripts/annotate.py'), config=CORR_CFG
    log: str(OUT / 'logs/shard_corresponding_windows.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=10000
    container: IMAGE('python')
    shell:
        'python3 {params.script:q} shard --role illumina --config {params.config:q} --disk-mb {resources.disk_mb} 2> {log:q}'


rule fimo_corresponding_final:
    input:
        fasta=corresponding_fasta,
        models=str(OUT / 'library/final_motifs.meme'),
        background=str(OUT / 'library/scoring_background.bfile'),
        freeze=str(OUT / 'library/freeze.json')
    output:
        scan=str(OUT / 'application/corresponding/scans/{shard}.tsv'),
        success=str(OUT / 'application/corresponding/scans/{shard}.tsv.success.json'),
        provenance=str(OUT / 'application/corresponding/scans/{shard}.tsv.provenance.json')
    wildcard_constraints: shard='s[0-9]{6}'
    params: script=str(ROOT / 'scripts/annotate.py'), config=CORR_CFG
    log: str(OUT / 'logs/fimo_corresponding_final/{shard}.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='04:00:00', disk_mb=30000
    container: IMAGE('meme')
    shell:
        'python3 {params.script:q} scan --config {params.config:q} --disk-mb {resources.disk_mb} --fasta {input.fasta:q} --output {output.scan:q} --log {log:q}'


rule annotate_corresponding_genotypes:
    input:
        scans=corresponding_scans, success=corresponding_success,
        samples=corresponding_prepared,
        freeze=str(OUT / 'library/freeze.json'),
        manifest=str(OUT / 'application/corresponding/manifest.json'),
        shards=str(OUT / 'application/corresponding/shards/manifest.json')
    output:
        complete=str(OUT / 'application/corresponding/annotation.complete.json'),
        genotypes=str(OUT / 'application/corresponding/genotypes.tsv.gz'),
        alt=str(OUT / 'application/corresponding/alt_annotations.tsv.gz'),
        matches=str(OUT / 'application/corresponding/all_matches.tsv.gz'),
        counts=str(OUT / 'reports/corresponding/annotation_counts.tsv'),
        positions=str(OUT / 'reports/corresponding/application_positions.tsv'),
        provenance=str(OUT / 'application/corresponding/annotation.complete.json.provenance.json')
    params: script=str(ROOT / 'scripts/annotate.py'), config=CORR_CFG
    log: str(OUT / 'logs/annotate_corresponding_genotypes.stderr.log')
    threads: 1
    resources: mem_mb=8000, time='04:00:00', disk_mb=60000
    container: IMAGE('python')
    shell:
        'python3 {params.script:q} annotate --role illumina --config {params.config:q} --disk-mb {resources.disk_mb} --scans {input.scans:q} 2> {log:q}'


rule render_reports:
    input:
        freeze=str(OUT / 'library/freeze.json'),
        counts=str(OUT / 'reports/annotation_counts.tsv'),
        coverage=str(OUT / 'reports/motif_coverage.tsv'),
        decisions=str(OUT / 'reports/selection_decisions.tsv'),
        representatives=str(OUT / 'reports/representative_validation.tsv'),
        bins=str(OUT / 'reports/conditional_position_bins.tsv')
    output:
        selection=str(OUT / 'reports/motif_selection.pdf'),
        positions=str(OUT / 'reports/conditional_positions.pdf'),
        readme=str(OUT / 'reports/README.md')
    params: script=str(ROOT / 'scripts/report.R'), out=str(OUT)
    log: str(OUT / 'logs/render_reports.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=5000
    container: IMAGE('r')
    shell:
        'Rscript {params.script:q} {params.out:q} 2> {log:q}'


rule report_complete:
    input:
        selection=str(OUT / 'reports/motif_selection.pdf'),
        positions=str(OUT / 'reports/conditional_positions.pdf'),
        readme=str(OUT / 'reports/README.md')
    output:
        complete=str(OUT / 'reports/report.complete.json'),
        provenance=str(OUT / 'reports/report.complete.json.provenance.json')
    params: script=str(ROOT / 'scripts/annotate.py'), config=APP_CFG
    log: str(OUT / 'logs/report_complete.stderr.log')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    container: IMAGE('python')
    shell:
        'python3 {params.script:q} report-complete --config {params.config:q} --disk-mb {resources.disk_mb} 2> {log:q}'


def vcf_output_datasets(wildcards=None):
    with open(checkpoints.plan_vcf_outputs.get().output.datasets, encoding='utf-8') as handle:
        return json.load(handle)


def vcf_output_groups(wildcards=None):
    with open(checkpoints.plan_vcf_outputs.get().output.groups, encoding='utf-8') as handle:
        return json.load(handle)


def vcf_export_completions(wildcards=None):
    return [str(OUT / 'vcf_outputs' / row['dataset_id'] / 'export.complete.json')
            for row in vcf_output_datasets(wildcards)['datasets']]


def vcf_group_member_ids(wildcards=None):
    return list(dict.fromkeys(did for group in vcf_output_groups(wildcards)['groups']
                              for did in group['member_dataset_ids']))


def vcf_sample_row(wildcards):
    matches = [row for row in vcf_output_datasets(wildcards)['datasets']
               if row['dataset_id'] == wildcards.vcf_dataset]
    if len(matches) != 1:
        raise ValueError('Unknown VCF export dataset: ' + wildcards.vcf_dataset)
    return matches[0]


def vcf_sample_artifacts(wildcards):
    row = vcf_sample_row(wildcards)
    app = OUT / ('application/corresponding' if row['role'] == 'illumina' else 'application')
    prepared = OUT / 'prepared' / row['dataset_id']
    paths = [Path(row['source_vcf']), Path(row['annotation_path']),
             Path(row['annotation_complete']), Path(row['annotation_complete'] + '.provenance.json')]
    paths.extend(app / name for name in (
        'manifest.json', 'manifest.json.provenance.json', 'alt_annotations.tsv.gz'))
    paths.extend(prepared / name for name in (
        'records.tsv.gz', 'records.complete.json', 'records.complete.json.provenance.json',
        'genotypes.tsv.gz', 'windows.tsv', 'sample.complete.json',
        'sample.complete.json.provenance.json'))
    return list(map(str, paths))


def vcf_group_source_alleles(wildcards):
    vcf_sample_row(wildcards)
    if wildcards.vcf_dataset not in vcf_group_member_ids(wildcards):
        raise ValueError('Normalization requires a configured group member: ' + wildcards.vcf_dataset)
    return str(OUT / 'vcf_outputs' / wildcards.vcf_dataset / 'source_alleles.vcf')

def vcf_group_row(wildcards):
    matches = [row for row in vcf_output_groups(wildcards)['groups']
               if row['group_id'] == wildcards.vcf_group]
    if len(matches) != 1:
        raise ValueError('Unknown VCF consensus group: ' + wildcards.vcf_group)
    return matches[0]


def vcf_group_artifacts(wildcards):
    datasets = {row['dataset_id']: row for row in vcf_output_datasets(wildcards)['datasets']}
    group = vcf_group_row(wildcards)
    paths = []
    for dataset_id in group['member_dataset_ids']:
        row = datasets[dataset_id]
        sample = OUT / 'vcf_outputs' / dataset_id
        app = OUT / ('application/corresponding' if row['role'] == 'illumina' else 'application')
        prepared = OUT / 'prepared' / dataset_id
        complete = Path(row['annotation_complete'])
        paths.extend([Path(row['source_vcf']), Path(row['annotation_path']), complete,
                      Path(str(complete) + '.provenance.json')])
        paths.extend(app / name for name in (
            'manifest.json', 'manifest.json.provenance.json', 'genotypes.tsv.gz',
            'alt_annotations.tsv.gz', 'all_matches.tsv.gz'))
        paths.extend(prepared / name for name in (
            'records.tsv.gz', 'records.complete.json', 'records.complete.json.provenance.json',
            'genotypes.tsv.gz', 'all.fasta', 'windows.tsv', 'sample.complete.json',
            'sample.complete.json.provenance.json'))
        paths.extend(sample / name for name in (
            'export.complete.json', 'export.complete.json.provenance.json', 'promising.vcf.gz',
            'promising.vcf.gz.csi', 'promising.vcf.gz.provenance.json', 'bgzip.version.txt',
            'promising.raw.vcf', 'source_alleles.vcf', 'source_alleles.tsv.gz',
            'source_alleles.tsv.gz.provenance.json', 'inputs.complete.json',
            'inputs.complete.json.provenance.json', 'normalized_alleles.vcf',
            'alleles.tsv.gz', 'alleles.tsv.gz.provenance.json', 'alleles.complete.json',
            'alleles.complete.json.provenance.json'))
    return list(map(str, dict.fromkeys(paths)))


def vcf_group_completions(wildcards=None):
    return [str(OUT / 'vcf_outputs/groups' / group['group_id'] / 'complete.json')
            for group in vcf_output_groups(wildcards)['groups']]


def vcf_output_sample_artifacts(wildcards=None):
    paths = []
    for row in vcf_output_datasets(wildcards)['datasets']:
        sample = OUT / 'vcf_outputs' / row['dataset_id']
        paths.extend(sample / name for name in (
            'export.complete.json', 'export.complete.json.provenance.json', 'promising.vcf.gz',
            'promising.vcf.gz.csi', 'promising.vcf.gz.provenance.json', 'bgzip.version.txt'))
    return list(map(str, paths))


def vcf_output_group_artifacts(wildcards=None):
    paths = []
    for group in vcf_output_groups(wildcards)['groups']:
        directory = OUT / 'vcf_outputs/groups' / group['group_id']
        paths.extend(directory / name for name in (
            'complete.json', 'complete.json.provenance.json', 'consensus.vcf.gz',
            'consensus.vcf.gz.csi', 'consensus.vcf.gz.provenance.json', 'support.tsv.gz',
            'support.tsv.gz.provenance.json', 'member_evidence.tsv.gz',
            'member_evidence.tsv.gz.provenance.json', 'excluded_alleles.tsv.gz',
            'excluded_alleles.tsv.gz.provenance.json'))
    return list(map(str, paths))

checkpoint plan_vcf_outputs:
    input:
        external=str(OUT / 'application/annotation.complete.json'),
        external_provenance=str(OUT / 'application/annotation.complete.json.provenance.json'),
        external_manifest=str(OUT / 'application/manifest.json'),
        external_annotations=str(OUT / 'application/genotypes.tsv.gz'),
        corresponding=str(OUT / 'application/corresponding/annotation.complete.json'),
        corresponding_provenance=str(OUT / 'application/corresponding/annotation.complete.json.provenance.json'),
        corresponding_manifest=str(OUT / 'application/corresponding/manifest.json'),
        corresponding_annotations=str(OUT / 'application/corresponding/genotypes.tsv.gz'),
        sources=lambda wc: [config['illumina_vcf'], config['other_illumina']],
        freeze=str(OUT / 'library/freeze.json'),
        reference=str(OUT / 'reference/reference.fa'),
        index=str(OUT / 'reference/reference.fa.fai'),
        config=APP_CFG, corresponding_config=CORR_CFG,
        script=str(ROOT / 'scripts/vcf_outputs.py')
    output:
        datasets=str(OUT / 'vcf_outputs/datasets.json'),
        groups=str(OUT / 'vcf_outputs/groups.json'),
        fingerprints=str(OUT / 'vcf_outputs/input_fingerprints.json'),
        datasets_provenance=str(OUT / 'vcf_outputs/datasets.json.provenance.json'),
        groups_provenance=str(OUT / 'vcf_outputs/groups.json.provenance.json')
    log: str(OUT / 'logs/plan_vcf_outputs.stderr.log')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    container: IMAGE('python')
    shell:
        'python3 {input.script:q} manifest --config {input.config:q} --disk-mb {resources.disk_mb} 2> {log:q}'


rule plan_vcf_sample_inputs:
    input:
        datasets=lambda wc: checkpoints.plan_vcf_outputs.get().output.datasets,
        datasets_provenance=lambda wc: checkpoints.plan_vcf_outputs.get().output.datasets_provenance,
        artifacts=vcf_sample_artifacts,
        freeze=str(OUT / 'library/freeze.json'),
        reference=str(OUT / 'reference/reference.fa'),
        index=str(OUT / 'reference/reference.fa.fai'),
        script=str(ROOT / 'scripts/vcf_outputs.py')
    output:
        promising=str(OUT / 'vcf_outputs/{vcf_dataset}/promising.raw.vcf'),
        sites=str(OUT / 'vcf_outputs/{vcf_dataset}/source_alleles.vcf'),
        source_map=str(OUT / 'vcf_outputs/{vcf_dataset}/source_alleles.tsv.gz'),
        complete=str(OUT / 'vcf_outputs/{vcf_dataset}/inputs.complete.json'),
        provenance=str(OUT / 'vcf_outputs/{vcf_dataset}/inputs.complete.json.provenance.json'),
        raw_provenance=str(OUT / 'vcf_outputs/{vcf_dataset}/promising.raw.vcf.provenance.json'),
        sites_provenance=str(OUT / 'vcf_outputs/{vcf_dataset}/source_alleles.vcf.provenance.json'),
        map_provenance=str(OUT / 'vcf_outputs/{vcf_dataset}/source_alleles.tsv.gz.provenance.json')
    wildcard_constraints: vcf_dataset='(?:illumina|other_illumina)_[a-f0-9]{16}'
    params: config=APP_CFG
    log: str(OUT / 'logs/plan_vcf_sample_inputs/{vcf_dataset}.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='02:00:00', disk_mb=30000
    container: IMAGE('python')
    shell:
        'python3 {input.script:q} sample-inputs --config {params.config:q} --dataset-id {wildcards.vcf_dataset:q} --disk-mb {resources.disk_mb} 2> {log:q}'


rule compress_promising_vcf:
    input:
        raw=str(OUT / 'vcf_outputs/{vcf_dataset}/promising.raw.vcf'),
        complete=str(OUT / 'vcf_outputs/{vcf_dataset}/inputs.complete.json'),
        provenance=str(OUT / 'vcf_outputs/{vcf_dataset}/inputs.complete.json.provenance.json')
    output:
        vcf=str(OUT / 'vcf_outputs/{vcf_dataset}/promising.vcf.gz'),
        index=str(OUT / 'vcf_outputs/{vcf_dataset}/promising.vcf.gz.csi'),
        version=str(OUT / 'vcf_outputs/{vcf_dataset}/bgzip.version.txt')
    wildcard_constraints: vcf_dataset='(?:illumina|other_illumina)_[a-f0-9]{16}'
    log: str(OUT / 'logs/compress_promising_vcf/{vcf_dataset}.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='02:00:00', disk_mb=30000
    container: IMAGE('bcftools')
    shell:
        'bgzip --version > {output.version:q}.tmp 2> {log:q}; '
        'bgzip -c {input.raw:q} > {output.vcf:q}.tmp 2>> {log:q}; '
        'bcftools index --csi --force -o {output.index:q}.tmp {output.vcf:q}.tmp 2>> {log:q}; '
        'mv {output.vcf:q}.tmp {output.vcf:q}; '
        'mv {output.index:q}.tmp {output.index:q}; '
        'mv {output.version:q}.tmp {output.version:q}'


rule capture_promising_export:
    input:
        datasets=lambda wc: checkpoints.plan_vcf_outputs.get().output.datasets,
        artifacts=vcf_sample_artifacts,
        vcf=str(OUT / 'vcf_outputs/{vcf_dataset}/promising.vcf.gz'),
        index=str(OUT / 'vcf_outputs/{vcf_dataset}/promising.vcf.gz.csi'),
        version=str(OUT / 'vcf_outputs/{vcf_dataset}/bgzip.version.txt'),
        raw=str(OUT / 'vcf_outputs/{vcf_dataset}/promising.raw.vcf'),
        source_map=str(OUT / 'vcf_outputs/{vcf_dataset}/source_alleles.tsv.gz'),
        sites=str(OUT / 'vcf_outputs/{vcf_dataset}/source_alleles.vcf'),
        complete=str(OUT / 'vcf_outputs/{vcf_dataset}/inputs.complete.json'),
        input_provenance=str(OUT / 'vcf_outputs/{vcf_dataset}/inputs.complete.json.provenance.json'),
        freeze=str(OUT / 'library/freeze.json'),
        reference=str(OUT / 'reference/reference.fa'),
        reference_index=str(OUT / 'reference/reference.fa.fai'),
        script=str(ROOT / 'scripts/vcf_outputs.py')
    output:
        complete=str(OUT / 'vcf_outputs/{vcf_dataset}/export.complete.json'),
        provenance=str(OUT / 'vcf_outputs/{vcf_dataset}/export.complete.json.provenance.json'),
        vcf_provenance=str(OUT / 'vcf_outputs/{vcf_dataset}/promising.vcf.gz.provenance.json')
    wildcard_constraints: vcf_dataset='(?:illumina|other_illumina)_[a-f0-9]{16}'
    params: config=APP_CFG
    log: str(OUT / 'logs/capture_promising_export/{vcf_dataset}.stderr.log')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    container: IMAGE('python')
    shell:
        'python3 {input.script:q} capture-export --config {params.config:q} --dataset-id {wildcards.vcf_dataset:q} --bgzip-version {input.version:q} --disk-mb {resources.disk_mb} 2> {log:q}'


rule normalize_vcf_alleles:
    input:
        sites=vcf_group_source_alleles,
        groups=lambda wc: checkpoints.plan_vcf_outputs.get().output.groups,
        source_map=str(OUT / 'vcf_outputs/{vcf_dataset}/source_alleles.tsv.gz'),
        complete=str(OUT / 'vcf_outputs/{vcf_dataset}/inputs.complete.json'),
        provenance=str(OUT / 'vcf_outputs/{vcf_dataset}/inputs.complete.json.provenance.json'),
        freeze=str(OUT / 'library/freeze.json'),
        reference=str(OUT / 'reference/reference.fa'),
        index=str(OUT / 'reference/reference.fa.fai')
    output: normalized=str(OUT / 'vcf_outputs/{vcf_dataset}/normalized_alleles.vcf')
    wildcard_constraints: vcf_dataset='(?:illumina|other_illumina)_[a-f0-9]{16}'
    log: str(OUT / 'logs/normalize_vcf_alleles/{vcf_dataset}.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='02:00:00', disk_mb=30000
    container: IMAGE('bcftools')
    shell:
        'bcftools norm --no-version -f {input.reference:q} -c e -Ov '
        '-o {output.normalized:q}.tmp {input.sites:q} 2> {log:q}; '
        'mv {output.normalized:q}.tmp {output.normalized:q}'


rule capture_vcf_alleles:
    input:
        datasets=lambda wc: checkpoints.plan_vcf_outputs.get().output.datasets,
        groups=lambda wc: checkpoints.plan_vcf_outputs.get().output.groups,
        sites=vcf_group_source_alleles,
        artifacts=vcf_sample_artifacts,
        normalized=str(OUT / 'vcf_outputs/{vcf_dataset}/normalized_alleles.vcf'),
        source_map=str(OUT / 'vcf_outputs/{vcf_dataset}/source_alleles.tsv.gz'),
        complete=str(OUT / 'vcf_outputs/{vcf_dataset}/inputs.complete.json'),
        input_provenance=str(OUT / 'vcf_outputs/{vcf_dataset}/inputs.complete.json.provenance.json'),
        freeze=str(OUT / 'library/freeze.json'),
        reference=str(OUT / 'reference/reference.fa'),
        index=str(OUT / 'reference/reference.fa.fai'),
        script=str(ROOT / 'scripts/vcf_outputs.py')
    output:
        alleles=str(OUT / 'vcf_outputs/{vcf_dataset}/alleles.tsv.gz'),
        complete=str(OUT / 'vcf_outputs/{vcf_dataset}/alleles.complete.json'),
        provenance=str(OUT / 'vcf_outputs/{vcf_dataset}/alleles.complete.json.provenance.json'),
        map_provenance=str(OUT / 'vcf_outputs/{vcf_dataset}/alleles.tsv.gz.provenance.json')
    wildcard_constraints: vcf_dataset='(?:illumina|other_illumina)_[a-f0-9]{16}'
    params: config=APP_CFG
    log: str(OUT / 'logs/capture_vcf_alleles/{vcf_dataset}.stderr.log')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    container: IMAGE('python')
    shell:
        'python3 {input.script:q} capture-alleles --config {params.config:q} --dataset-id {wildcards.vcf_dataset:q} --normalized-vcf {input.normalized:q} --disk-mb {resources.disk_mb} 2> {log:q}'


rule group_consensus:
    input:
        datasets=lambda wc: checkpoints.plan_vcf_outputs.get().output.datasets,
        datasets_provenance=lambda wc: checkpoints.plan_vcf_outputs.get().output.datasets_provenance,
        groups=lambda wc: checkpoints.plan_vcf_outputs.get().output.groups,
        groups_provenance=lambda wc: checkpoints.plan_vcf_outputs.get().output.groups_provenance,
        artifacts=vcf_group_artifacts,
        freeze=str(OUT / 'library/freeze.json'),
        reference=str(OUT / 'reference/reference.fa'),
        reference_index=str(OUT / 'reference/reference.fa.fai'),
        script=str(ROOT / 'scripts/vcf_outputs.py')
    output:
        raw=temp(str(OUT / 'vcf_outputs/groups/{vcf_group}/consensus.raw.vcf')),
        support=str(OUT / 'vcf_outputs/groups/{vcf_group}/support.tsv.gz'),
        evidence=str(OUT / 'vcf_outputs/groups/{vcf_group}/member_evidence.tsv.gz'),
        excluded=str(OUT / 'vcf_outputs/groups/{vcf_group}/excluded_alleles.tsv.gz')
    wildcard_constraints: vcf_group='[a-z0-9][a-z0-9_-]*'
    params: config=APP_CFG
    log: str(OUT / 'logs/group_consensus/{vcf_group}.stderr.log')
    threads: 1
    resources: mem_mb=8000, time='04:00:00', disk_mb=60000
    container: IMAGE('python')
    shell:
        'python3 {input.script:q} consensus --config {params.config:q} --group-id {wildcards.vcf_group:q} --disk-mb {resources.disk_mb} 2> {log:q}'


rule compress_consensus_vcf:
    input: raw=str(OUT / 'vcf_outputs/groups/{vcf_group}/consensus.raw.vcf')
    output:
        vcf=str(OUT / 'vcf_outputs/groups/{vcf_group}/consensus.vcf.gz'),
        index=str(OUT / 'vcf_outputs/groups/{vcf_group}/consensus.vcf.gz.csi')
    wildcard_constraints: vcf_group='[a-z0-9][a-z0-9_-]*'
    log: str(OUT / 'logs/compress_consensus_vcf/{vcf_group}.stderr.log')
    threads: 1
    resources: mem_mb=4000, time='02:00:00', disk_mb=30000
    container: IMAGE('bcftools')
    shell:
        'bgzip -c {input.raw:q} > {output.vcf:q}.tmp 2> {log:q}; '
        'bcftools index --csi --force -o {output.index:q}.tmp {output.vcf:q}.tmp 2>> {log:q}; '
        'mv {output.vcf:q}.tmp {output.vcf:q}; '
        'mv {output.index:q}.tmp {output.index:q}'


rule capture_group_consensus:
    input:
        datasets=lambda wc: checkpoints.plan_vcf_outputs.get().output.datasets,
        datasets_provenance=lambda wc: checkpoints.plan_vcf_outputs.get().output.datasets_provenance,
        groups=lambda wc: checkpoints.plan_vcf_outputs.get().output.groups,
        groups_provenance=lambda wc: checkpoints.plan_vcf_outputs.get().output.groups_provenance,
        artifacts=vcf_group_artifacts,
        raw=str(OUT / 'vcf_outputs/groups/{vcf_group}/consensus.raw.vcf'),
        vcf=str(OUT / 'vcf_outputs/groups/{vcf_group}/consensus.vcf.gz'),
        index=str(OUT / 'vcf_outputs/groups/{vcf_group}/consensus.vcf.gz.csi'),
        support=str(OUT / 'vcf_outputs/groups/{vcf_group}/support.tsv.gz'),
        evidence=str(OUT / 'vcf_outputs/groups/{vcf_group}/member_evidence.tsv.gz'),
        excluded=str(OUT / 'vcf_outputs/groups/{vcf_group}/excluded_alleles.tsv.gz'),
        freeze=str(OUT / 'library/freeze.json'),
        reference=str(OUT / 'reference/reference.fa'),
        reference_index=str(OUT / 'reference/reference.fa.fai'),
        script=str(ROOT / 'scripts/vcf_outputs.py')
    output:
        complete=str(OUT / 'vcf_outputs/groups/{vcf_group}/complete.json'),
        provenance=str(OUT / 'vcf_outputs/groups/{vcf_group}/complete.json.provenance.json'),
        vcf_provenance=str(OUT / 'vcf_outputs/groups/{vcf_group}/consensus.vcf.gz.provenance.json'),
        support_provenance=str(OUT / 'vcf_outputs/groups/{vcf_group}/support.tsv.gz.provenance.json'),
        evidence_provenance=str(OUT / 'vcf_outputs/groups/{vcf_group}/member_evidence.tsv.gz.provenance.json'),
        excluded_provenance=str(OUT / 'vcf_outputs/groups/{vcf_group}/excluded_alleles.tsv.gz.provenance.json')
    wildcard_constraints: vcf_group='[a-z0-9][a-z0-9_-]*'
    params: config=APP_CFG
    log: str(OUT / 'logs/capture_group_consensus/{vcf_group}.stderr.log')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    container: IMAGE('python')
    shell:
        'python3 {input.script:q} capture-group --config {params.config:q} --group-id {wildcards.vcf_group:q} --raw-vcf {input.raw:q} --disk-mb {resources.disk_mb} 2> {log:q}'




rule complete_vcf_outputs:
    input:
        datasets=lambda wc: checkpoints.plan_vcf_outputs.get().output.datasets,
        datasets_provenance=lambda wc: checkpoints.plan_vcf_outputs.get().output.datasets_provenance,
        groups=lambda wc: checkpoints.plan_vcf_outputs.get().output.groups,
        groups_provenance=lambda wc: checkpoints.plan_vcf_outputs.get().output.groups_provenance,
        sample_exports=vcf_export_completions,
        sample_artifacts=vcf_output_sample_artifacts,
        group_completions=vcf_group_completions,
        group_artifacts=vcf_output_group_artifacts,
        freeze=str(OUT / 'library/freeze.json'),
        reference=str(OUT / 'reference/reference.fa'),
        reference_index=str(OUT / 'reference/reference.fa.fai'),
        script=str(ROOT / 'scripts/vcf_outputs.py')
    output:
        complete=str(OUT / 'vcf_outputs/outputs.complete.json'),
        provenance=str(OUT / 'vcf_outputs/outputs.complete.json.provenance.json')
    params: config=APP_CFG
    log: str(OUT / 'logs/complete_vcf_outputs.stderr.log')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    container: IMAGE('python')
    shell:
        'python3 {input.script:q} complete --config {params.config:q} --disk-mb {resources.disk_mb} 2> {log:q}'


