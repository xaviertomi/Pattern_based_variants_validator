EVIDENCE_SCRIPT = str(ROOT / 'scripts/evidence.py')
BG = str(OUT / 'reference/scoring_background.bfile')
GATE = str(OUT / 'manifests/preparation.ok.json')
EV = str(OUT / 'evidence/{dataset}/r{rotation}')

wildcard_constraints:
    rotation='0[1-5]',
    dataset='(?:pacbio|illumina)_[0-9a-f]{16}',
    kind='real|control'

rule make_scoring_background:
    input: reference=str(OUT / 'reference/reference.fa'), gate=GATE
    output: main=BG, provenance=BG+'.provenance.json'
    params: script=EVIDENCE_SCRIPT, cfg=CFG
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=10000
    container: IMAGE('meme')
    log: str(OUT / 'logs/background.rule.log')
    shell: 'python3 {params.script:q} background --config {params.cfg:q} --disk-mb {resources.disk_mb} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} > {log:q} 2>&1'

rule make_train_test:
    input:
        gate=GATE,
        folds=lambda w: [str(OUT / 'prepared' / w.dataset / f'fold{r}.fasta') for r in ROTATIONS],
        assignment=str(OUT / 'manifests/folds.json')
    output: main=EV+'/split.json', train=EV+'/train.fasta', test=EV+'/test.fasta', provenance=EV+'/split.json.provenance.json'
    params: script=EVIDENCE_SCRIPT, cfg=CFG
    threads: 1
    resources: mem_mb=4000, time='01:00:00', disk_mb=10000
    container: IMAGE('python')
    log: EV+'/split.rule.log'
    shell: 'python3 {params.script:q} split --config {params.cfg:q} --disk-mb {resources.disk_mb} --dataset {wildcards.dataset:q} --rotation {wildcards.rotation:q} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.train:q} {output.test:q} > {log:q} 2>&1'

rule shuffle_train:
    input: fasta=EV+'/train.fasta', split=EV+'/split.json', seeds=str(OUT / 'manifests/seeds.tsv'), gate=GATE
    output: main=EV+'/train.control.fasta', mapping=EV+'/train.control.fasta.map.json', provenance=EV+'/train.control.fasta.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=4000,time='01:00:00',disk_mb=10000
    container: IMAGE('meme')
    log: EV+'/shuffle_train.rule.log'
    shell: 'python3 {params.script:q} shuffle --config {params.cfg:q} --disk-mb {resources.disk_mb} --dataset {wildcards.dataset:q} --rotation {wildcards.rotation:q} --kind train --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.mapping:q} > {log:q} 2>&1'

rule shuffle_test:
    input: fasta=lambda w:str(OUT / 'prepared' / w.dataset / f'fold{int(w.rotation)}.fasta'), seeds=str(OUT / 'manifests/seeds.tsv'), gate=GATE
    output: main=str(OUT / 'controls/{dataset}/fold{rotation}.fasta'), mapping=str(OUT / 'controls/{dataset}/fold{rotation}.fasta.map.json'), provenance=str(OUT / 'controls/{dataset}/fold{rotation}.fasta.provenance.json')
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=4000,time='01:00:00',disk_mb=10000
    container: IMAGE('meme')
    log: str(OUT / 'logs/shuffle_test/{dataset}/r{rotation}.log')
    shell: 'python3 {params.script:q} shuffle --config {params.cfg:q} --disk-mb {resources.disk_mb} --dataset {wildcards.dataset:q} --rotation {wildcards.rotation:q} --kind test --output {output.main:q} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.mapping:q} > {log:q} 2>&1'

rule streme:
    input: train=EV+'/train.fasta', control=EV+'/train.control.fasta', seeds=str(OUT / 'manifests/seeds.tsv'), gate=GATE
    output: main=EV+'/streme.complete.json', text=EV+'/streme.txt', xml=EV+'/streme.xml', provenance=EV+'/streme.complete.json.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=64000,time='12:00:00',disk_mb=20000
    container: IMAGE('meme')
    log: EV+'/streme.rule.log'
    shell: 'python3 {params.script:q} streme --config {params.cfg:q} --disk-mb {resources.disk_mb} --dataset {wildcards.dataset:q} --rotation {wildcards.rotation:q} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.text:q} {output.xml:q} > {log:q} 2>&1'

rule parse_models:
    input: completed=EV+'/streme.complete.json', text=EV+'/streme.txt', xml=EV+'/streme.xml', background=BG
    output: main=EV+'/models.json', meme=EV+'/models.meme', provenance=EV+'/models.json.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=4000,time='01:00:00',disk_mb=10000
    container: IMAGE('python')
    log: EV+'/parse_models.rule.log'
    shell: 'python3 {params.script:q} parse-models --config {params.cfg:q} --disk-mb {resources.disk_mb} --dataset {wildcards.dataset:q} --rotation {wildcards.rotation:q} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.meme:q} > {log:q} 2>&1'

rule fimo_pacbio:
    input:
        models=EV+'/models.meme', background=BG, gate=GATE,
        fasta=lambda w:str(OUT / ('prepared' if w.kind=='real' else 'controls') / w.dataset / f"fold{int(w.rotation) if w.kind=='real' else w.rotation}.fasta"),
        preflight=str(OUT / 'manifests/preflight/meme.ok.json')
    output: main=EV+'/pacbio.{kind}.tsv', success=EV+'/pacbio.{kind}.tsv.success.json', provenance=EV+'/pacbio.{kind}.tsv.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG,stderr=lambda w:str(OUT/'logs/fimo_pacbio'/w.dataset/f'r{w.rotation}'/f'{w.kind}.stderr.log')
    threads: 1
    resources: mem_mb=4000,time='04:00:00',disk_mb=30000
    container: IMAGE('meme')
    log: EV+'/pacbio.{kind}.rule.log'
    shell: 'python3 {params.script:q} scan --config {params.cfg:q} --disk-mb {resources.disk_mb} --models {input.models:q} --fasta {input.fasta:q} --background {input.background:q} --output {output.main:q} --log {params.stderr:q} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.success:q} > {log:q} 2>&1'

rule summarize_pacbio:
    input: real=EV+'/pacbio.real.tsv', control=EV+'/pacbio.control.tsv', real_ok=EV+'/pacbio.real.tsv.success.json', control_ok=EV+'/pacbio.control.tsv.success.json', models=EV+'/models.meme', fasta=lambda w:str(OUT/'prepared'/w.dataset/f'fold{int(w.rotation)}.fasta')
    output: main=EV+'/pacbio.metrics.tsv', database=EV+'/pacbio.metrics.tsv.sqlite', real_audit=EV+'/pacbio.metrics.tsv.real.occurrences.tsv.gz', control_audit=EV+'/pacbio.metrics.tsv.control.occurrences.tsv.gz', provenance=EV+'/pacbio.metrics.tsv.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=8000,time='02:00:00',disk_mb=60000
    container: IMAGE('python')
    log: EV+'/summarize_pacbio.rule.log'
    shell: 'python3 {params.script:q} summarize --config {params.cfg:q} --disk-mb {resources.disk_mb} --dataset {wildcards.dataset:q} --rotation {wildcards.rotation:q} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.database:q} {output.real_audit:q} {output.control_audit:q} > {log:q} 2>&1'

rule select_candidates:
    input: models=EV+'/models.json', metrics=EV+'/pacbio.metrics.tsv', background=BG
    output: main=EV+'/candidates.json', meme=EV+'/candidates.meme', provenance=EV+'/candidates.json.provenance.json'
    params: script=EVIDENCE_SCRIPT,cfg=CFG
    threads: 1
    resources: mem_mb=2000,time='00:15:00',disk_mb=1000
    container: IMAGE('python')
    log: EV+'/candidates.rule.log'
    shell: 'python3 {params.script:q} candidates --config {params.cfg:q} --disk-mb {resources.disk_mb} --dataset {wildcards.dataset:q} --rotation {wildcards.rotation:q} --primary {output.main:q} --inputs {input:q} --outputs {output.main:q} {output.meme:q} > {log:q} 2>&1'
