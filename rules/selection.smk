SEL = str(OUT / 'selection')
REP = SEL+'/representatives/{member}/{dataset}/r{rotation}'

def representative_manifest(wildcards):
    return checkpoints.prepare_representatives.get().output.main

def representative_model(wildcards):
    manifest=representative_manifest(wildcards)
    with open(manifest,encoding='utf-8') as handle:
        candidates=json.load(handle)['candidates']
    matches=[m for m in candidates if m['id']==wildcards.member]
    if len(matches)!=1 or wildcards.dataset not in matches[0]['evaluated_dataset_ids'] or int(wildcards.rotation)!=matches[0]['source_rotation']:
        raise ValueError('Representative scan not in predetermined held-out evaluation manifest')
    return str(OUT/'selection/representative_models'/wildcards.member/'model.meme')

def representative_summaries(wildcards):
    manifest=representative_manifest(wildcards)
    with open(manifest,encoding='utf-8') as handle: data=json.load(handle)
    return [str(OUT/'selection/representatives'/m['id']/ds/f"r{m['source_rotation']:02d}.metrics.tsv") for m in data['candidates'] for ds in m['evaluated_dataset_ids']]

rule pool_models:
    input: expand(str(OUT/'evidence/{dataset}/r{rotation}/transfer.complete.json'),dataset=PB_IDS,rotation=[f'{r:02d}' for r in ROTATIONS]),BG
    output: main=SEL+'/pool.json',meme=SEL+'/transferable_models.meme',provenance=SEL+'/pool.json.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=8000,time='04:00:00',disk_mb=10000
    container: IMAGE('python')
    log: SEL+'/pool.rule.log'
    shell: 'python3 {params.script:q} pool --config {params.cfg:q} --disk-mb {resources.disk_mb} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.meme:q} > {log:q} 2>&1'

rule tomtom:
    input: pool=SEL+'/pool.json',models=SEL+'/transferable_models.meme',preflight=str(OUT/'manifests/preflight/meme.ok.json')
    output: main=SEL+'/tomtom.tsv',provenance=SEL+'/tomtom.tsv.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=8000,time='04:00:00',disk_mb=10000
    container: IMAGE('meme')
    log: SEL+'/tomtom.rule.log'
    shell: 'python3 {params.script:q} tomtom --config {params.cfg:q} --disk-mb {resources.disk_mb} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} > {log:q} 2>&1'

rule make_families:
    input: pool=SEL+'/pool.json',tomtom=SEL+'/tomtom.tsv'
    output: main=SEL+'/families.json',table=str(OUT/'reports/family_membership.tsv'),provenance=SEL+'/families.json.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=8000,time='04:00:00',disk_mb=10000
    container: IMAGE('python')
    log: SEL+'/families.rule.log'
    shell: 'python3 {params.script:q} families --config {params.cfg:q} --disk-mb {resources.disk_mb} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.table:q} > {log:q} 2>&1'

rule summarize_family_unions:
    input:
        families=SEL+'/families.json',
        pacbio=expand(str(OUT/'evidence/{dataset}/r{rotation}/pacbio.metrics.tsv.sqlite'),dataset=PB_IDS,rotation=[f'{r:02d}' for r in ROTATIONS]),
        corresponding=[str(OUT/'evidence'/pb/f'r{r:02d}'/f'{il}.metrics.tsv.sqlite') for pb in PB_IDS for il in IL_IDS if DATASETS[il]['corresponding_pacbio']==DATASETS[pb]['sample_name'] for r in ROTATIONS],
        windows=[str(OUT/'prepared'/ds/f'fold{r}.fasta') for ds in PB_IDS+IL_IDS for r in ROTATIONS]
    output: main=SEL+'/stability.json',coverage=str(OUT/'reports/family_coverage.tsv'),provenance=SEL+'/stability.json.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=8000,time='04:00:00',disk_mb=60000
    container: IMAGE('python')
    log: SEL+'/stability.rule.log'
    shell: 'python3 {params.script:q} stability --config {params.cfg:q} --disk-mb {resources.disk_mb} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.coverage:q} > {log:q} 2>&1'

checkpoint prepare_representatives:
    input: stability=SEL+'/stability.json',background=BG
    output: main=SEL+'/representatives/manifest.json',models=directory(SEL+'/representative_models'),provenance=SEL+'/representatives/manifest.json.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=4000,time='01:00:00',disk_mb=10000
    container: IMAGE('python')
    log: SEL+'/representatives/prepare.rule.log'
    shell: 'python3 {params.script:q} representatives --config {params.cfg:q} --disk-mb {resources.disk_mb} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} > {log:q} 2>&1'

rule fimo_representative:
    input:
        manifest=representative_manifest,models=representative_model,background=BG,gate=GATE,
        fasta=lambda w:str(OUT / ('prepared' if w.kind=='real' else 'controls') / w.dataset / f"fold{int(w.rotation) if w.kind=='real' else w.rotation}.fasta"),
        preflight=str(OUT/'manifests/preflight/meme.ok.json')
    output: main=REP+'.{kind}.tsv',success=REP+'.{kind}.tsv.success.json',provenance=REP+'.{kind}.tsv.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG,stderr=lambda w:str(OUT/'logs/fimo_representative'/w.member/w.dataset/f'r{w.rotation}'/f'{w.kind}.stderr.log')
    threads: 1
    resources: mem_mb=4000,time='04:00:00',disk_mb=30000
    container: IMAGE('meme')
    log: REP+'.{kind}.rule.log'
    shell: 'python3 {params.script:q} scan --config {params.cfg:q} --disk-mb {resources.disk_mb} --models {input.models:q} --fasta {input.fasta:q} --background {input.background:q} --output {output.main:q} --log {params.stderr:q} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.success:q} > {log:q} 2>&1'

rule summarize_representative:
    input: real=REP+'.real.tsv',control=REP+'.control.tsv',real_ok=REP+'.real.tsv.success.json',control_ok=REP+'.control.tsv.success.json',models=representative_model,manifest=representative_manifest,fasta=lambda w:str(OUT/'prepared'/w.dataset/f'fold{int(w.rotation)}.fasta')
    output: main=REP+'.metrics.tsv',database=REP+'.metrics.tsv.sqlite',real_audit=REP+'.metrics.tsv.real.occurrences.tsv.gz',control_audit=REP+'.metrics.tsv.control.occurrences.tsv.gz',provenance=REP+'.metrics.tsv.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=8000,time='02:00:00',disk_mb=60000
    container: IMAGE('python')
    log: REP+'.summary.rule.log'
    shell: 'python3 {params.script:q} summarize --config {params.cfg:q} --disk-mb {resources.disk_mb} --dataset {wildcards.dataset:q} --member {wildcards.member:q} --rotation {wildcards.rotation:q} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.database:q} {output.real_audit:q} {output.control_audit:q} > {log:q} 2>&1'

rule validate_representatives:
    input: manifest=representative_manifest,summaries=representative_summaries,background=BG
    output: main=SEL+'/representatives/validation.json',table=str(OUT/'reports/representative_validation.tsv'),provenance=SEL+'/representatives/validation.json.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=8000,time='04:00:00',disk_mb=60000
    container: IMAGE('python')
    log: SEL+'/representatives/validation.rule.log'
    shell: 'python3 {params.script:q} validate-representatives --config {params.cfg:q} --disk-mb {resources.disk_mb} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.table:q} > {log:q} 2>&1'

rule rank_families:
    input: stability=SEL+'/stability.json',validation=SEL+'/representatives/validation.json',coverage=str(OUT/'reports/family_coverage.tsv')
    output: main=SEL+'/ranking.json',table=str(OUT/'reports/selection_decisions.tsv'),provenance=SEL+'/ranking.json.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=8000,time='04:00:00',disk_mb=10000
    container: IMAGE('python')
    log: SEL+'/ranking.rule.log'
    shell: 'python3 {params.script:q} rank --config {params.cfg:q} --disk-mb {resources.disk_mb} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.table:q} > {log:q} 2>&1'

rule freeze:
    input:
        ranking=SEL+'/ranking.json',background=BG,selection_config=str(OUT/'manifests/selection_config.json'),
        genotypes=[str(OUT/'prepared'/ds/'genotypes.tsv.gz') for ds in PB_IDS+IL_IDS],
        evidence=[str(OUT/p) for p in ['reports/representative_validation.tsv','selection/representatives/validation.json','selection/stability.json','reports/selection_decisions.tsv','manifests/folds.json','manifests/seeds.tsv','reference/reference.fa','reference/reference.fa.fai','manifests/selection_fingerprints.json','manifests/datasets.tsv','manifests/preflight.ok.json','manifests/preflight/meme.ok.json']],
        summaries=[str(OUT/'evidence'/pb/f'r{r:02d}'/name) for pb in PB_IDS for r in ROTATIONS for name in ['pacbio.metrics.tsv','pacbio.metrics.tsv.sqlite','transfer.tsv','streme.complete.json','candidates.json']]
    output:
        main=str(OUT/'library/freeze.json'),models=str(OUT/'library/final_motifs.meme'),manifest=str(OUT/'library/final_motif_manifest.tsv'),members=str(OUT/'library/final_motif_members.tsv'),background=str(OUT/'library/scoring_background.bfile'),loci=str(OUT/'library/selection_loci.sqlite'),reports=[str(OUT/'reports'/name) for name in ['motif_coverage.tsv','transfer_evidence.tsv','streme_usage.tsv','conditional_position_bins.tsv']],provenance=str(OUT/'library/freeze.json.provenance.json')
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=8000,time='04:00:00',disk_mb=10000
    container: IMAGE('python')
    log: str(OUT/'library/freeze.rule.log')
    shell: 'python3 {params.script:q} freeze --config {params.cfg:q} --disk-mb {resources.disk_mb} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.models:q} {output.manifest:q} {output.members:q} {output.background:q} {output.loci:q} {output.reports:q} > {log:q} 2>&1'
