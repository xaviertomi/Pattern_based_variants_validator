TRANSFER_EV = str(OUT / 'evidence/{source}/r{rotation}')

rule fimo_corresponding:
    input:
        models=TRANSFER_EV+'/candidates.meme', background=BG, gate=GATE,
        fasta=lambda w:str(OUT / ('prepared' if w.kind=='real' else 'controls') / w.dataset / f"fold{int(w.rotation) if w.kind=='real' else w.rotation}.fasta"),
        preflight=str(OUT / 'manifests/preflight/meme.ok.json')
    output: main=TRANSFER_EV+'/{dataset}.{kind}.tsv',success=TRANSFER_EV+'/{dataset}.{kind}.tsv.success.json',provenance=TRANSFER_EV+'/{dataset}.{kind}.tsv.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG,stderr=lambda w:str(OUT/'logs/fimo_corresponding'/w.source/f'r{w.rotation}'/f'{w.dataset}.{w.kind}.stderr.log')
    threads: 1
    resources: mem_mb=4000,time='04:00:00',disk_mb=30000
    container: IMAGE('meme')
    log: TRANSFER_EV+'/{dataset}.{kind}.rule.log'
    shell: 'python3 {params.script:q} scan --config {params.cfg:q} --disk-mb {resources.disk_mb} --models {input.models:q} --fasta {input.fasta:q} --background {input.background:q} --output {output.main:q} --log {params.stderr:q} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.success:q} > {log:q} 2>&1'

rule summarize_transfer:
    input: real=TRANSFER_EV+'/{dataset}.real.tsv',control=TRANSFER_EV+'/{dataset}.control.tsv',real_ok=TRANSFER_EV+'/{dataset}.real.tsv.success.json',control_ok=TRANSFER_EV+'/{dataset}.control.tsv.success.json',models=TRANSFER_EV+'/candidates.meme',fasta=lambda w:str(OUT/'prepared'/w.dataset/f'fold{int(w.rotation)}.fasta'),windows=str(OUT/'prepared/{dataset}/windows.tsv'),source_genotypes=str(OUT/'prepared/{source}/genotypes.tsv.gz')
    output: main=TRANSFER_EV+'/{dataset}.metrics.tsv',database=TRANSFER_EV+'/{dataset}.metrics.tsv.sqlite',subsets=TRANSFER_EV+'/{dataset}.metrics.tsv.subsets.tsv',real_audit=TRANSFER_EV+'/{dataset}.metrics.tsv.real.occurrences.tsv.gz',control_audit=TRANSFER_EV+'/{dataset}.metrics.tsv.control.occurrences.tsv.gz',provenance=TRANSFER_EV+'/{dataset}.metrics.tsv.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=8000,time='02:00:00',disk_mb=60000
    container: IMAGE('python')
    log: TRANSFER_EV+'/{dataset}.summary.rule.log'
    shell: 'python3 {params.script:q} summarize --config {params.cfg:q} --disk-mb {resources.disk_mb} --dataset {wildcards.dataset:q} --source {wildcards.source:q} --rotation {wildcards.rotation:q} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.database:q} {output.subsets:q} {output.real_audit:q} {output.control_audit:q} > {log:q} 2>&1'

rule select_transferable:
    input:
        candidates=EV+'/candidates.json',
        metrics=lambda w:[str(OUT/'evidence'/w.dataset/f'r{w.rotation}'/f'{i}.metrics.tsv') for i in IL_IDS if DATASETS[i]['corresponding_pacbio']==DATASETS[w.dataset]['sample_name']]
    output: main=EV+'/transfer.complete.json',table=EV+'/transfer.tsv',provenance=EV+'/transfer.complete.json.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=2000,time='00:15:00',disk_mb=1000
    container: IMAGE('python')
    log: EV+'/transfer.rule.log'
    shell: 'python3 {params.script:q} transfer --config {params.cfg:q} --disk-mb {resources.disk_mb} --dataset {wildcards.dataset:q} --rotation {wildcards.rotation:q} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.table:q} > {log:q} 2>&1'
