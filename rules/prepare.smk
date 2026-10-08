# Each executable runs only in its own image; Python consumes declared captures.
PREP_SCRIPT = str(ROOT / 'scripts/prepare.py')
PREP_IDS = PB_IDS + IL_IDS
PREP_REF = str(OUT / 'reference/reference.fa')
PREP_FAI = PREP_REF + '.fai'
PREP_GATE = str(OUT / 'manifests/preparation.ok.json')
PREP_FOLDS = str(OUT / 'manifests/folds.json')
PREP_PREFLIGHT = str(OUT / 'manifests/preflight.ok.json')

localrules: write_effective_config, fingerprint_selection

rule write_effective_config:
    output: CONFIG_JSON
    params: effective={key: value for key, value in config.items() if key != 'workflow_scope'}
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    run:
        import sys
        sys.path.insert(0, str(ROOT / 'scripts'))
        from prepare import atomic_json, validate_config
        atomic_json(output[0], validate_config(params.effective, ROOT))

rule fingerprint_selection:
    input: config=CONFIG_JSON
    output:
        fingerprint=str(OUT / 'manifests/selection_fingerprints.json'),
        scientific=str(OUT / 'manifests/selection_config.json'),
        datasets=str(OUT / 'manifests/datasets.tsv')
    params: effective=config
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    run:
        import sys, json
        sys.path.insert(0, str(ROOT / 'scripts'))
        from prepare import atomic_json, validate_config, selection_config, fingerprint_inputs, datasets, write_tsv
        effective = validate_config(params.effective, ROOT)
        scientific = selection_config(effective)
        atomic_json(output.scientific, scientific)
        paths = {key: effective[key] for key in ('reference_fasta', 'pacbio_vcf', 'illumina_vcf')}
        paths.update({str(ROOT / 'scripts' / name): ROOT / 'scripts' / name for name in ('prepare.py', 'evidence.py')})
        paths.update({str(path): path for path in sorted((ROOT / 'rules').glob('*.smk')) if path.name != 'application.smk'})
        paths['Snakefile'] = ROOT / 'Snakefile'
        fingerprint_inputs(paths, scientific, output.fingerprint)
        write_tsv(output.datasets, (row | {'metadata': json.dumps(row['metadata'], sort_keys=True)} for row in datasets(effective).values()), ['dataset_id', 'role', 'sample_name', 'source_vcf', 'corresponding_pacbio', 'metadata'])

rule stage_reference:
    input: config=CONFIG_JSON, source=config['reference_fasta']
    output: reference=PREP_REF, provenance=PREP_REF + '.provenance.json'
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=4000, time='02:00:00', disk_mb=30000
    shell: 'python3 {PREP_SCRIPT:q} stage-reference --config {input.config:q} --root {ROOT:q} --disk-mb {resources.disk_mb}'

rule index_reference:
    input: reference=PREP_REF
    output: index=PREP_FAI
    log: str(OUT / 'logs/index_reference/reference.stderr.log')
    container: IMAGE('samtools')
    threads: 1
    resources: mem_mb=4000, time='02:00:00', disk_mb=30000
    shell: 'mkdir -p $(dirname {log:q}); samtools faidx {input.reference:q} 2> {log:q}'

rule capture_reference_index:
    input: config=CONFIG_JSON, reference=PREP_REF, index=PREP_FAI
    output: provenance=PREP_FAI + '.provenance.json'
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=10000
    shell: 'python3 {PREP_SCRIPT:q} index-reference --phase capture --config {input.config:q} --root {ROOT:q} --disk-mb {resources.disk_mb}'

rule make_folds:
    input: config=CONFIG_JSON, reference=PREP_REF, index=PREP_FAI, provenance=PREP_FAI + '.provenance.json'
    output:
        manifest=PREP_FOLDS, provenance=PREP_FOLDS + '.provenance.json',
        assignments=str(OUT / 'manifests/fold_assignments.tsv'), balance=str(OUT / 'reports/fold_balance.tsv'), unmatched=str(OUT / 'reports/unmatched_contigs.tsv'),
        beds=[str(OUT / f'folds/fold{r}.bed') for r in ROTATIONS]
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=10000
    shell: 'python3 {PREP_SCRIPT:q} folds --config {input.config:q} --root {ROOT:q} --disk-mb {resources.disk_mb}'

rule validate_fold_partition:
    input: config=CONFIG_JSON, reference=PREP_REF, index=PREP_FAI, manifest=PREP_FOLDS, beds=[str(OUT / f'folds/fold{r}.bed') for r in ROTATIONS]
    output: gate=str(OUT / 'manifests/folds.validated.json'), provenance=str(OUT / 'manifests/folds.validated.json.provenance.json')
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=10000
    shell: 'python3 {PREP_SCRIPT:q} validate-folds --config {input.config:q} --root {ROOT:q} --disk-mb {resources.disk_mb}'

rule stage_shared_vcf:
    input: fingerprint=str(OUT / 'manifests/selection_fingerprints.json')
    output:
        vcf=str(OUT / 'inputs/{role}.vcf.gz'), index=str(OUT / 'inputs/{role}.vcf.gz.csi'),
        samples=str(OUT / 'inputs/{role}.samples.txt'), header=str(OUT / 'inputs/{role}.header.txt')
    params: source=lambda w: config[w.role + '_vcf']
    wildcard_constraints: role='pacbio|illumina'
    log: str(OUT / 'logs/stage_vcf/{role}.stderr.log')
    container: IMAGE('bcftools')
    threads: 1
    resources: mem_mb=4000, time='02:00:00', disk_mb=30000
    shell:
        'mkdir -p $(dirname {output.vcf:q}) $(dirname {log:q}); '
        'bcftools view -Oz -o {output.vcf:q}.tmp.gz {params.source:q} 2> {log:q}; '
        'bcftools index --csi --force {output.vcf:q}.tmp.gz 2>> {log:q}; '
        'bcftools query -l {output.vcf:q}.tmp.gz > {output.samples:q}.tmp 2>> {log:q}; '
        'bcftools view -h {output.vcf:q}.tmp.gz > {output.header:q}.tmp 2>> {log:q}; '
        'mv {output.vcf:q}.tmp.gz {output.vcf:q}; mv {output.vcf:q}.tmp.gz.csi {output.index:q}; '
        'mv {output.samples:q}.tmp {output.samples:q}; mv {output.header:q}.tmp {output.header:q}'

rule capture_shared_vcf:
    input: config=CONFIG_JSON, vcf=str(OUT / 'inputs/{role}.vcf.gz'), index=str(OUT / 'inputs/{role}.vcf.gz.csi'), samples=str(OUT / 'inputs/{role}.samples.txt'), header=str(OUT / 'inputs/{role}.header.txt')
    output: provenance=str(OUT / 'inputs/{role}.vcf.gz.provenance.json')
    wildcard_constraints: role='pacbio|illumina'
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    shell: 'python3 {PREP_SCRIPT:q} stage-vcf --phase capture --config {input.config:q} --root {ROOT:q} --role {wildcards.role:q} --disk-mb {resources.disk_mb}'

rule validate_shared_headers:
    input:
        config=CONFIG_JSON, reference=PREP_REF, fai=PREP_FAI,
        vcfs=[str(OUT / f'inputs/{role}.vcf.gz') for role in ('pacbio', 'illumina')],
        samples=[str(OUT / f'inputs/{role}.samples.txt') for role in ('pacbio', 'illumina')],
        headers=[str(OUT / f'inputs/{role}.header.txt') for role in ('pacbio', 'illumina')],
        provenance=[str(OUT / f'inputs/{role}.vcf.gz.provenance.json') for role in ('pacbio', 'illumina')]
    output: gate=str(OUT / 'manifests/headers.ok.json'), provenance=str(OUT / 'manifests/headers.ok.json.provenance.json'), lists=[str(OUT / f'manifests/{role}_header_samples.tsv') for role in ('pacbio', 'illumina')]
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    shell: 'python3 {PREP_SCRIPT:q} headers --phase capture --config {input.config:q} --root {ROOT:q} --disk-mb {resources.disk_mb}'

rule selected_exact_sample:
    input: config=CONFIG_JSON
    output: sample=str(OUT / 'prepared/{dataset_id}/exact_sample.txt'), provenance=str(OUT / 'prepared/{dataset_id}/exact_sample.txt.provenance.json')
    wildcard_constraints: dataset_id='(?:pacbio|illumina)_[0-9a-f]{16}'
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    shell: 'python3 {PREP_SCRIPT:q} sample-name --config {input.config:q} --dataset-id {wildcards.dataset_id:q} --disk-mb {resources.disk_mb}'

rule query_selected_records:
    input:
        sample=str(OUT / 'prepared/{dataset_id}/exact_sample.txt'), gate=str(OUT / 'manifests/headers.ok.json'),
        vcf=lambda w: str(OUT / f'inputs/{DATASETS[w.dataset_id]["role"]}.vcf.gz'), reference=PREP_REF, fai=PREP_FAI
    output:
        records=str(OUT / 'prepared/{dataset_id}/raw.records.tsv'), alt_records=str(OUT / 'prepared/{dataset_id}/raw.alt.tsv'),
        vcf=str(OUT / 'prepared/{dataset_id}/sample.vcf.gz'), index=str(OUT / 'prepared/{dataset_id}/sample.vcf.gz.csi'),
        alt=str(OUT / 'prepared/{dataset_id}/sample.alt.vcf.gz'), alt_index=str(OUT / 'prepared/{dataset_id}/sample.alt.vcf.gz.csi'), success=str(OUT / 'prepared/{dataset_id}/records.tool.success')
    wildcard_constraints: dataset_id='(?:pacbio|illumina)_[0-9a-f]{16}'
    log: str(OUT / 'logs/prepare_sample/{dataset_id}/records.stderr.log')
    params: alt_filter=ALT_FILTER_EXPR
    container: IMAGE('bcftools')
    threads: 1
    resources: mem_mb=4000, time='02:00:00', disk_mb=30000
    shell:
        "mkdir -p $(dirname {log:q}); "
        "bcftools view -S {input.sample:q} -Oz -o {output.vcf:q} {input.vcf:q} 2> {log:q}; "
        "bcftools index --csi --force {output.vcf:q} 2>> {log:q}; "
        "bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT[\\t%GT]\\n' {output.vcf:q} > {output.records:q}.tmp 2>> {log:q}; "
        "bcftools view -i {params.alt_filter:q} -Oz -o {output.alt:q} {output.vcf:q} 2>> {log:q}; "
        "bcftools index --csi --force {output.alt:q} 2>> {log:q}; "
        "bcftools norm -f {input.reference:q} -c e -Ou {output.alt:q} > /dev/null 2>> {log:q}; "
        "bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT[\\t%GT]\\n' {output.alt:q} > {output.alt_records:q}.tmp 2>> {log:q}; "
        "mv {output.records:q}.tmp {output.records:q}; mv {output.alt_records:q}.tmp {output.alt_records:q}; printf 'success\\n' > {output.success:q}"

rule prepare_selected_records:
    input:
        config=CONFIG_JSON, raw=str(OUT / 'prepared/{dataset_id}/raw.records.tsv'), alt=str(OUT / 'prepared/{dataset_id}/raw.alt.tsv'), success=str(OUT / 'prepared/{dataset_id}/records.tool.success')
    output:
        complete=str(OUT / 'prepared/{dataset_id}/records.complete.json'), provenance=str(OUT / 'prepared/{dataset_id}/records.complete.json.provenance.json'),
        records=str(OUT / 'prepared/{dataset_id}/records.tsv.gz'), capture=str(OUT / 'prepared/{dataset_id}/raw.records.tsv.provenance.json')
    wildcard_constraints: dataset_id='(?:pacbio|illumina)_[0-9a-f]{16}'
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=4000, time='02:00:00', disk_mb=30000
    shell:
        'python3 {PREP_SCRIPT:q} capture-sample-tools --config {input.config:q} --dataset-id {wildcards.dataset_id:q} --disk-mb {resources.disk_mb}; '
        'python3 {PREP_SCRIPT:q} prepare-sample --config {input.config:q} --dataset-id {wildcards.dataset_id:q} --phase records-from-query --disk-mb {resources.disk_mb}'

rule plan_selected_contexts:
    input: config=CONFIG_JSON, records=str(OUT / 'prepared/{dataset_id}/records.tsv.gz'), complete=str(OUT / 'prepared/{dataset_id}/records.complete.json'), folds=PREP_FOLDS, partition=str(OUT / 'manifests/folds.validated.json')
    output: plan=str(OUT / 'prepared/{dataset_id}/contexts.plan.json'), provenance=protected(str(OUT / 'prepared/{dataset_id}/contexts.plan.json.provenance.json')), regions=protected(str(OUT / 'prepared/{dataset_id}/regions.txt'))
    wildcard_constraints: dataset_id='(?:pacbio|illumina)_[0-9a-f]{16}'
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=10000
    shell: 'python3 {PREP_SCRIPT:q} prepare-sample --config {input.config:q} --dataset-id {wildcards.dataset_id:q} --phase contexts-plan --disk-mb {resources.disk_mb}'

rule extract_selected_contexts:
    input: regions=str(OUT / 'prepared/{dataset_id}/regions.txt'), plan=str(OUT / 'prepared/{dataset_id}/contexts.plan.json'), reference=PREP_REF, fai=PREP_FAI
    output: fasta=str(OUT / 'prepared/{dataset_id}/extracted.fasta'), success=str(OUT / 'prepared/{dataset_id}/extraction.tool.success')
    wildcard_constraints: dataset_id='(?:pacbio|illumina)_[0-9a-f]{16}'
    log: str(OUT / 'logs/prepare_sample/{dataset_id}/contexts.stderr.log')
    container: IMAGE('samtools')
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=10000
    shell:
        'mkdir -p $(dirname {log:q}); if test -s {input.regions:q}; then samtools faidx -r {input.regions:q} {input.reference:q} > {output.fasta:q}.tmp 2> {log:q}; else : > {output.fasta:q}.tmp; : > {log:q}; fi; '
        "mv {output.fasta:q}.tmp {output.fasta:q}; printf 'success\\n' > {output.success:q}"

rule prepare_selected_contexts:
    input: config=CONFIG_JSON, fasta=str(OUT / 'prepared/{dataset_id}/extracted.fasta'), success=str(OUT / 'prepared/{dataset_id}/extraction.tool.success'), records=str(OUT / 'prepared/{dataset_id}/records.tsv.gz')
    output:
        complete=str(OUT / 'prepared/{dataset_id}/sample.complete.json'), provenance=str(OUT / 'prepared/{dataset_id}/sample.complete.json.provenance.json'), capture=str(OUT / 'prepared/{dataset_id}/extracted.fasta.provenance.json'),
        genotypes=str(OUT / 'prepared/{dataset_id}/genotypes.tsv.gz'), windows=str(OUT / 'prepared/{dataset_id}/windows.tsv'), fasta=str(OUT / 'prepared/{dataset_id}/all.fasta'), counts=str(OUT / 'prepared/{dataset_id}/fold_counts.tsv'),
        fold_fasta=[str(OUT / ('prepared/{dataset_id}/fold' + str(r) + '.fasta')) for r in ROTATIONS], fold_windows=[str(OUT / ('prepared/{dataset_id}/fold' + str(r) + '.windows.tsv')) for r in ROTATIONS]
    wildcard_constraints: dataset_id='(?:pacbio|illumina)_[0-9a-f]{16}'
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=10000
    shell:
        'python3 {PREP_SCRIPT:q} capture-extraction --config {input.config:q} --dataset-id {wildcards.dataset_id:q} --disk-mb {resources.disk_mb}; '
        'python3 {PREP_SCRIPT:q} prepare-sample --config {input.config:q} --dataset-id {wildcards.dataset_id:q} --phase contexts-from-fasta --disk-mb {resources.disk_mb}'

rule make_seeds:
    input: config=CONFIG_JSON
    output: seeds=str(OUT / 'manifests/seeds.tsv'), provenance=str(OUT / 'manifests/seeds.tsv.provenance.json')
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    shell: 'python3 {PREP_SCRIPT:q} seeds --config {input.config:q} --root {ROOT:q} --disk-mb {resources.disk_mb}'

rule preflight_fixtures:
    input: config=CONFIG_JSON
    output:
        complete=str(OUT / 'preflight/fixtures.complete.json'), provenance=str(OUT / 'preflight/fixtures.complete.json.provenance.json'),
        fixtures=[str(OUT / 'preflight' / tool / name) for tool in ('bcftools', 'samtools', 'r') for name in ('tiny.fa', 'tiny.fa.fai', 'tiny.vcf')]
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    shell: 'python3 {PREP_SCRIPT:q} preflight-fixtures --config {input.config:q} --disk-mb {resources.disk_mb}'

rule preflight_bcftools_tools:
    input: fixtures=str(OUT / 'preflight/fixtures.complete.json')
    output: success=str(OUT / 'preflight/bcftools/tool.success'), version=str(OUT / 'preflight/bcftools/version.txt'), query=str(OUT / 'preflight/bcftools/query.tsv')
    params: directory=str(OUT / 'preflight/bcftools')
    log: str(OUT / 'logs/preflight/bcftools.stderr.log')
    container: IMAGE('bcftools')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    shell:
        'mkdir -p $(dirname {log:q}); bcftools --version > {output.version:q} 2> {log:q}; '
        'bcftools view -s S -Oz -o {params.directory:q}/tiny.vcf.gz {params.directory:q}/tiny.vcf 2>> {log:q}; '
        'bcftools index --csi --force {params.directory:q}/tiny.vcf.gz 2>> {log:q}; '
        "bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT[\\t%GT]\\n' {params.directory:q}/tiny.vcf.gz > {output.query:q} 2>> {log:q}; "
        "bcftools norm -f {params.directory:q}/tiny.fa -c e -Ou {params.directory:q}/tiny.vcf.gz > /dev/null 2>> {log:q}; printf 'success\\n' > {output.success:q}"

rule preflight_samtools_tools:
    input: fixtures=str(OUT / 'preflight/fixtures.complete.json')
    output: success=str(OUT / 'preflight/samtools/tool.success'), version=str(OUT / 'preflight/samtools/version.txt'), fasta=str(OUT / 'preflight/samtools/extracted.fasta')
    params: directory=str(OUT / 'preflight/samtools')
    log: str(OUT / 'logs/preflight/samtools.stderr.log')
    container: IMAGE('samtools')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    shell:
        'mkdir -p $(dirname {log:q}); samtools --version > {output.version:q} 2> {log:q}; '
        'samtools faidx {params.directory:q}/tiny.fa 2>> {log:q}; '
        "samtools faidx {params.directory:q}/tiny.fa tiny:2-5 > {output.fasta:q} 2>> {log:q}; printf 'success\\n' > {output.success:q}"

rule preflight_r_tools:
    input: fixtures=str(OUT / 'preflight/fixtures.complete.json')
    output: success=str(OUT / 'preflight/r/tool.success'), version=str(OUT / 'preflight/r/version.txt'), result=str(OUT / 'preflight/r/result.txt')
    log: str(OUT / 'logs/preflight/r.stderr.log')
    container: IMAGE('r')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    shell:
        "mkdir -p $(dirname {log:q}); Rscript --version > {output.version:q} 2>&1; Rscript -e 'cat(sum(c(1,2,3)))' > {output.result:q} 2> {log:q}; printf 'success\\n' > {output.success:q}"

rule capture_worker_preflight:
    input:
        config=CONFIG_JSON, success=str(OUT / 'preflight/{tool}/tool.success'), version=str(OUT / 'preflight/{tool}/version.txt'),
        result=lambda w: str(OUT / 'preflight' / w.tool / {'bcftools': 'query.tsv', 'samtools': 'extracted.fasta', 'r': 'result.txt'}[w.tool])
    output:
        gate=str(OUT / 'manifests/preflight/{tool}.ok.json'), provenance=str(OUT / 'manifests/preflight/{tool}.ok.json.provenance.json'),
        capture=str(OUT / 'preflight/{tool}/tool.success.provenance.json')
    wildcard_constraints: tool='bcftools|samtools|r'
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    shell: 'python3 {PREP_SCRIPT:q} preflight --phase capture --config {input.config:q} --tool {wildcards.tool:q} --disk-mb {resources.disk_mb}'

rule preflight_meme:
    input: config=CONFIG_JSON
    output: gate=str(OUT / 'manifests/preflight/meme.ok.json'), provenance=str(OUT / 'manifests/preflight/meme.ok.json.provenance.json')
    container: IMAGE('meme')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    shell: 'python3 {PREP_SCRIPT:q} preflight --config {input.config:q} --root {ROOT:q} --tool meme --disk-mb {resources.disk_mb}'

rule preflight_python:
    input: config=CONFIG_JSON
    output: gate=str(OUT / 'manifests/preflight/python.ok.json'), provenance=str(OUT / 'manifests/preflight/python.ok.json.provenance.json')
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    shell: 'python3 {PREP_SCRIPT:q} preflight --config {input.config:q} --root {ROOT:q} --tool python --disk-mb {resources.disk_mb}'

rule tool_preflight_gate:
    input: config=CONFIG_JSON, tools=[str(OUT / f'manifests/preflight/{tool}.ok.json') for tool in ('bcftools', 'samtools', 'meme', 'python', 'r')]
    output: gate=PREP_PREFLIGHT, provenance=PREP_PREFLIGHT + '.provenance.json'
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=2000, time='00:15:00', disk_mb=1000
    shell: 'python3 {PREP_SCRIPT:q} preflight-gate --config {input.config:q} --root {ROOT:q} --disk-mb {resources.disk_mb}'

rule global_preparation_gate:
    input: config=CONFIG_JSON, headers=str(OUT / 'manifests/headers.ok.json'), folds=str(OUT / 'manifests/folds.validated.json'), preflight=PREP_PREFLIGHT, seeds=str(OUT / 'manifests/seeds.tsv'), samples=[str(OUT / 'prepared' / did / 'sample.complete.json') for did in PREP_IDS]
    output: gate=PREP_GATE, provenance=PREP_GATE + '.provenance.json', counts=str(OUT / 'reports/dataset_fold_counts.tsv'), exclusions=str(OUT / 'reports/exclusions.tsv')
    container: IMAGE('python')
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=10000
    shell: 'python3 {PREP_SCRIPT:q} gate --config {input.config:q} --root {ROOT:q} --disk-mb {resources.disk_mb}'
