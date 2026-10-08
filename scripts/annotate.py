"""Immutable-library application: shared interval scans, record-preserving annotation."""
import argparse
import csv
import gzip
import json
import math
import os
import re
import sqlite3
import sys
import tempfile
from collections import Counter
from decimal import Decimal, InvalidOperation
from pathlib import Path

from prepare import atomic_json, check_reuse, configure_scratch, dataset_id, fingerprint_inputs, input_signature, prepare_sample, publish_provenance, read_tsv, sha256, stage_vcf, validate_header, write_tsv

FIELDS = ['source_vcf','source_record_index','record_id','dataset_id','sample_name','role','selection_participant','evidence_context','CHROM','POS','ID','REF','ALT','GT','variant_type','partial_missing','ambiguity_count','window_id','interval_start','interval_end','extraction_status','extraction_reason','annotation_status','annotation_reason','identity_check_status','best_final_motif_id','best_family_id','distinct_motif_count','score','reported_p_value','matched_sequence','strand','local_start','local_stop','midpoint_distance','seen_in_selection_loci']
MATCH_FIELDS = ['motif_id','motif_alt_id','sequence_name','start','stop','strand','score','p-value','matched_sequence']


def load_json(path):
    with open(path, encoding='utf-8') as handle:
        return json.load(handle)


def provenance(primary, inputs, params, outputs):
    params = dict(params, code_paths=sorted(set(params.get('code_paths',[])) | {str(Path(__file__).resolve())}))
    return publish_provenance(primary,inputs,params,outputs,scope='application',producer_target=str(primary))

def application_fingerprint(config, output_dir):
    output_dir = Path(output_dir)
    verify_freeze(output_dir)
    paths = {'freeze':output_dir/'library/freeze.json','external':Path(config['other_illumina'])}
    params = {'external_metadata':config.get('sample_metadata',{}).get('other_illumina',{}),'code_sha256':sha256(Path(__file__)),'scan_windows_per_shard':config['execution']['scan_windows_per_shard']}
    return fingerprint_inputs(paths,params,output_dir/'application/input_fingerprints.json')

def annotation_paths(output_dir, role):
    if role == 'other_illumina':
        app = Path(output_dir) / 'application'
        reports = Path(output_dir) / 'reports'
    elif role == 'illumina':
        app = Path(output_dir) / 'application/corresponding'
        reports = Path(output_dir) / 'reports/corresponding'
    else:
        raise ValueError('Unknown application role: ' + role)
    return app, reports


def verify_prepared_artifact(primary, required):
    primary = Path(primary)
    sidecar = Path(str(primary) + '.provenance.json')
    try:
        record = load_json(sidecar)
        paths = record['input_paths']
        if record['scope'] != 'selection' or any(not Path(path).is_file() for path in paths):
            raise ValueError('invalid selection provenance inputs')
        if input_signature(paths, record['signature_params']) != record['input_signature']:
            raise ValueError('prepared input signature mismatch')
        expected = set(map(str, required))
        if not expected <= set(record['required_outputs']) or not check_reuse(sidecar, record['input_signature'], record['required_outputs']):
            raise ValueError('prepared output checksum mismatch')
    except (OSError, KeyError, TypeError, ValueError) as error:
        raise ValueError('Prepared selection artifact is incomplete: ' + str(primary) + ': ' + str(error)) from error
    return record


def corresponding_members(config, freeze):
    selected = freeze['selection_config']['illumina_corresponding']
    if config['illumina_corresponding'] != selected:
        raise ValueError('Configured corresponding samples disagree with the frozen selection; rerun selection and freeze_library')
    return [(pb, name) for pb, names in selected.items() for name in names]


def corresponding_fingerprint(config, output_dir):
    output_dir = Path(output_dir)
    freeze = verify_freeze(output_dir)
    members = corresponding_members(config, freeze)
    source = Path(config['illumina_vcf'])
    frozen_path = output_dir / 'manifests/selection_fingerprints.json'
    fingerprints = load_json(frozen_path)
    if fingerprints.get('params') != freeze['selection_config']:
        raise ValueError('Frozen selection fingerprint does not match freeze.json')
    source_digest = sha256(source)
    frozen_source = fingerprints['inputs']['illumina_vcf']['sha256']
    if source_digest != frozen_source:
        raise ValueError('Corresponding source VCF differs from the frozen selection input')
    if fingerprints['inputs']['reference_fasta']['sha256'] != freeze['reference_sha256']:
        raise ValueError('Frozen selection reference fingerprint mismatch')
    paths = {'freeze': output_dir / 'library/freeze.json',
             'selection_fingerprints': frozen_path,
             'selection_config': output_dir / 'manifests/selection_config.json',
             'source_vcf': source,
             'reference': Path(freeze['_reference']),
             'reference_index': Path(freeze['_reference'] + '.fai')}
    for _, name in members:
        did = dataset_id('illumina', name)
        directory = output_dir / 'prepared' / did
        paths.update({f'{did}/{file}': directory / file for file in (
            'records.complete.json', 'records.complete.json.provenance.json', 'records.tsv.gz',
            'sample.complete.json', 'sample.complete.json.provenance.json',
            'genotypes.tsv.gz', 'all.fasta', 'windows.tsv')})
        records = load_json(directory / 'records.complete.json')
        sample = load_json(directory / 'sample.complete.json')
        if records.get('dataset_id') != did or records.get('sample_name') != name or records.get('source_sha256') != frozen_source:
            raise ValueError(f'Frozen records completion mismatch: {name}')
        if sample.get('dataset_id') != did or sample.get('sample_name') != name or sample.get('reference_sha256') != freeze['reference_sha256']:
            raise ValueError(f'Frozen sample completion mismatch: {name}')
        verify_prepared_artifact(directory / 'records.complete.json',
                                 [directory / file for file in ('records.complete.json', 'records.tsv.gz')])
        verify_prepared_artifact(directory / 'sample.complete.json',
                                 [directory / file for file in ('sample.complete.json', 'genotypes.tsv.gz', 'all.fasta', 'windows.tsv')])
    params = {'freeze_sha256': sha256(output_dir / 'library/freeze.json'),
              'selection_fingerprint_sha256': sha256(frozen_path),
              'source_sha256': source_digest,
              'members': [{'corresponding_pacbio': pb, 'sample_name': name, 'dataset_id': dataset_id('illumina', name)} for pb, name in members],
              'code_sha256': sha256(Path(__file__)),
              'prepare_code_sha256': sha256(Path(__file__).with_name('prepare.py'))}
    return fingerprint_inputs(paths, params, output_dir / 'application/corresponding/input_fingerprints.json')


def corresponding_manifest(config, output_dir):
    output_dir = Path(output_dir)
    freeze = verify_freeze(output_dir)
    members = corresponding_members(config, freeze)
    fingerprint = corresponding_fingerprint(config, output_dir)
    source = Path(config['illumina_vcf'])
    source_hash = fingerprint['inputs']['source_vcf']['sha256']
    app, _ = annotation_paths(output_dir, 'illumina')
    rows = [{'dataset_id': dataset_id('illumina', name), 'role': 'illumina', 'sample_name': name,
             'source_vcf': str(source), 'corresponding_pacbio': pb, 'identity_check_status': 'frozen_selection_manifest',
             'selection_participant': 'true', 'evidence_context': 'selection_reuse'} for pb, name in members]
    fields = ['dataset_id', 'role', 'sample_name', 'source_vcf', 'corresponding_pacbio',
              'selection_participant', 'evidence_context']
    write_tsv(app / 'datasets.tsv', rows, fields)
    manifest = {'schema_version': 1, 'datasets': rows, 'source_vcf': str(source),
                'source_sha256': source_hash, 'reference': freeze['_reference'],
                'reference_sha256': freeze['reference_sha256'], 'flank_bp': freeze['flank_bp'],
                'identity_check_status': 'frozen_selection_manifest',
                'freeze_sha256': sha256(output_dir / 'library/freeze.json'),
                'selection_fingerprints_sha256': sha256(output_dir / 'manifests/selection_fingerprints.json'),
                'fingerprint_signature': fingerprint['signature']}
    atomic_json(app / 'manifest.json', manifest)
    return manifest


def scan_external(config, output_dir, fasta, output, log):
    from evidence import scan_fimo
    output_dir = Path(output_dir)
    freeze = verify_freeze(output_dir)
    flag = load_json(output_dir/'manifests/preflight/meme.ok.json')['fimo_background_flag']
    if flag not in ['--bgfile','--bfile']:
        raise ValueError('Pinned FIMO preflight did not resolve an explicit background flag')
    return scan_fimo(output_dir/'library/final_motifs.meme',Path(fasta),output_dir/'library/scoring_background.bfile',Path(output),freeze['fimo_p_threshold'],Path(log),flag)

def verify_freeze(output_dir):
    library = Path(output_dir) / 'library'
    try:
        freeze = load_json(library / 'freeze.json')
        outputs = freeze['outputs']
        if not outputs or not freeze.get('representative_validation') or not freeze.get('family_stability') or not any(Path(k).name=='representative_validation.tsv' for k in outputs):
            raise ValueError('missing exact-representative validation evidence')
        for name, expected in outputs.items():
            path = Path(name)
            if not path.is_absolute():
                path = Path(output_dir) / path
            if not path.is_file() or sha256(path) != expected:
                raise ValueError('frozen artifact checksum mismatch: ' + str(path))
        for name in ['final_motifs.meme','final_motif_manifest.tsv','final_motif_members.tsv','scoring_background.bfile','selection_loci.sqlite']:
            path = library / name
            if not path.is_file() or not any(Path(k).name == name for k in outputs):
                raise ValueError('unchecked frozen artifact: ' + name)
        ref = Path(freeze['reference_fasta'])
        if not ref.is_absolute():
            ref = Path(output_dir) / ref
        expected = freeze['reference_sha256']
        if not expected or not ref.is_file() or sha256(ref) != expected or not Path(str(ref)+'.fai').is_file():
            raise ValueError('missing or modified frozen reference')
        freeze['_reference'] = str(ref)
        return freeze
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise ValueError('Frozen library is missing or incomplete; run freeze_library first. ' + str(error)) from error


def identity_manifest(config, output_dir, prepared_headers=False):
    freeze = verify_freeze(output_dir)
    source = Path(config['other_illumina'])
    if prepared_headers:
        from prepare import validate_header_text
        staged = Path(output_dir)/'inputs/other_illumina.vcf.gz'
        names = validate_header_text(staged,Path(freeze['_reference']),[],(Path(output_dir)/'application/header.samples.txt').read_text().splitlines(),(Path(output_dir)/'application/header.vcf.txt').read_text())
    else:
        staged = stage_vcf(source, Path(output_dir), 'other_illumina')
        names = validate_header(staged,Path(freeze['_reference']),[])
    if not names or len(names) != len(set(names)):
        raise ValueError('External VCF must contain unique nonempty exact sample names')
    selection = freeze['selection_config']
    selected = set(selection['pacbio_samples'])
    selected.update(x for values in selection['illumina_corresponding'].values() for x in values)
    overlap = sorted(selected.intersection(names))
    if overlap:
        raise ValueError('External header overlaps selected samples: ' + ', '.join(overlap))
    metadata = selection['sample_metadata']
    external = config.get('sample_metadata', {}).get('other_illumina', {})
    extra = sorted(set(external) - set(names))
    if extra:
        raise ValueError('External metadata names absent from header: ' + ', '.join(extra))
    clones = {entry['clone_id'] for role in ['pacbio','illumina'] for entry in metadata.get(role, {}).values() if entry.get('clone_id')}
    for name in names:
        clone = external.get(name, {}).get('clone_id')
        if clone and clone in clones:
            raise ValueError('External biological clone overlaps selection: ' + name + ': ' + clone)
    verified = all(external.get(name, {}).get('clone_id') for name in names) and all(metadata.get(role, {}).get(name, {}).get('clone_id') for role, members in [('pacbio',selection['pacbio_samples']),('illumina',[x for v in selection['illumina_corresponding'].values() for x in v])] for name in members)
    status = 'explicit_clone_metadata' if verified else 'header_names_only'
    if not verified:
        print('WARNING: External identity checked by header names only; alias/biological independence is unverified.', file=sys.stderr)
    rows = [{'dataset_id':dataset_id('other_illumina', name),'role':'other_illumina','sample_name':name,'source_vcf':str(source),'corresponding_pacbio':'','identity_check_status':status,'selection_participant':'false','evidence_context':'external_application_metadata_checked' if verified else 'external_application_identity_unverified'} for name in names]
    if len({r['dataset_id'] for r in rows}) != len(rows):
        raise ValueError('External dataset ID collision')
    app = Path(output_dir) / 'application'
    app.mkdir(parents=True, exist_ok=True)
    write_tsv(app/'datasets.tsv',rows,list(rows[0]))
    manifest = {'schema_version':1,'datasets':rows,'source_vcf':str(source),'source_sha256':sha256(source),'staged_vcf':str(staged),'reference':freeze['_reference'],'flank_bp':freeze['flank_bp'],'identity_check_status':status}
    atomic_json(app/'manifest.json',manifest)
    atomic_json(app/'application_config.json',{'freeze_sha256':sha256(Path(output_dir)/'library/freeze.json'),'source_sha256':manifest['source_sha256'],'samples':names,'metadata':external,'analysis':{'flank_bp':manifest['flank_bp'],'fimo_p_threshold':freeze['fimo_p_threshold']}})
    return manifest


def prepare_external(config, output_dir, identifier, phase='all'):
    manifest = load_json(Path(output_dir)/'application/manifest.json')
    row = next((r for r in manifest['datasets'] if r['dataset_id']==identifier), None)
    if row is None:
        raise ValueError('Unknown external dataset: '+identifier)
    directory = Path(output_dir)/'prepared'/identifier
    directory.mkdir(parents=True,exist_ok=True)
    if phase=='records-plan':
        (directory/'exact_sample.txt').write_text(row['sample_name']+'\n',encoding='utf-8')
        atomic_json(directory/'records.plan.json',{'dataset_id':identifier,'sample_name':row['sample_name']})
        provenance(directory/'records.plan.json',[Path(output_dir)/'application/manifest.json'],{'sample_name':row['sample_name']},[directory/'records.plan.json',directory/'exact_sample.txt'])
        return
    result = prepare_sample(Path(manifest['staged_vcf']),Path(manifest['reference']),manifest['source_vcf'],manifest['source_sha256'],'other_illumina',row['sample_name'],Path(output_dir),manifest['flank_bp'],None,phase=phase,header_validated=True)
    if phase=='contexts-plan':
        atomic_json(directory/'contexts.plan.json',{'dataset_id':identifier,'flank_bp':manifest['flank_bp']})
        provenance(directory/'contexts.plan.json',[directory/'records.tsv.gz',Path(manifest['reference']),Path(manifest['reference']+'.fai')],{'flank_bp':manifest['flank_bp']},[directory/'contexts.plan.json',directory/'regions.txt'])
        return result
    records = phase=='records-from-query'
    names = ['records.tsv.gz','records.complete.json','sample.vcf.gz','sample.vcf.gz.csi','sample.alt.vcf.gz','sample.alt.vcf.gz.csi','raw.records.tsv','raw.alt.tsv'] if records else ['genotypes.tsv.gz','windows.tsv','all.fasta','sample.complete.json','extracted.fasta']
    inputs = [Path(output_dir)/'application/manifest.json',Path(manifest['reference']),Path(manifest['reference']+'.fai'),Path(manifest['staged_vcf']),directory/'records.plan.json'] if records else [directory/'records.tsv.gz',Path(manifest['reference']),Path(manifest['reference']+'.fai'),directory/'contexts.plan.json']
    provenance(directory/('records.complete.json' if records else 'sample.complete.json'),inputs,{'sample':row['sample_name'],'phase':phase,'flank_bp':manifest['flank_bp']},[directory/name for name in names])
    return result


def fasta_records(path):
    name, parts = None, []
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(parts)
                name, parts = line[1:].strip().split()[0], []
            elif line.strip():
                if name is None:
                    raise ValueError('Sequence before FASTA header: '+str(path))
                parts.append(line.strip().upper())
        if name is not None:
            yield name, ''.join(parts)


def shard_windows(output_dir, maximum, *, role):
    if isinstance(maximum,bool) or not isinstance(maximum,int) or maximum<=0:
        raise ValueError('scan_windows_per_shard must be a positive integer')
    output_dir = Path(output_dir)
    app, _ = annotation_paths(output_dir, role)
    manifest = load_json(app/'manifest.json')
    destination = app/'shards'
    destination.mkdir(parents=True,exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='windows-',suffix='.sqlite',dir=destination)
    os.close(fd)
    db = sqlite3.connect(temporary)
    db.execute('PRAGMA cache_size=-65536')
    db.execute('PRAGMA temp_store=FILE')
    db.execute('CREATE TABLE windows(window_id TEXT PRIMARY KEY, sequence TEXT NOT NULL, scoreable INTEGER NOT NULL DEFAULT 1)')
    try:
        for row in manifest['datasets']:
            for name, sequence in fasta_records(output_dir/'prepared'/row['dataset_id']/'all.fasta'):
                old = db.execute('SELECT sequence FROM windows WHERE window_id=?',(name,)).fetchone()
                if old and old[0] != sequence:
                    raise ValueError('Shared-window sequence collision: '+name)
                db.execute('INSERT OR IGNORE INTO windows(window_id,sequence) VALUES(?,?)',(name,sequence))
            db.commit()
        shards, count, handle = [], 0, None
        for name, sequence in db.execute('SELECT window_id,sequence FROM windows ORDER BY window_id'):
            if count % maximum == 0:
                if handle:
                    handle.close()
                identifier = 's%06d' % (len(shards)+1)
                path = destination/(identifier+'.fasta')
                handle = open(path,'w',encoding='utf-8',newline='\n')
                shards.append({'shard_id':identifier,'fasta':str(path),'window_count':0})
            handle.write('>'+name+'\n'+sequence+'\n')
            shards[-1]['window_count'] += 1
            count += 1
        if handle:
            handle.close()
        db.commit()
        db.close()
        os.replace(temporary,destination/'windows.sqlite')
        atomic_json(destination/'manifest.json',{'schema_version':1,'window_count':count,'shards':shards})
    finally:
        db.close()
        if os.path.exists(temporary):
            os.unlink(temporary)


def winner_key(row, rank):
    try:
        p, score = Decimal(str(row['p-value'])), Decimal(str(row['score']))
        start, stop = int(row['start']), int(row['stop'])
        strand = row['strand']
        if not p.is_finite() or not score.is_finite() or p<0 or p>1 or start<1 or stop<start or strand not in ['+','-']:
            raise ValueError('invalid FIMO occurrence')
        return p,int(rank),-score,start,stop,0 if strand=='+' else 1
    except (KeyError,ValueError,InvalidOperation) as error:
        raise ValueError('Malformed FIMO occurrence: '+str(row)) from error


def annotations(output_dir, scans, *, role):
    from evidence import iter_fimo
    output_dir = Path(output_dir)
    freeze = verify_freeze(output_dir)
    app, reports = annotation_paths(output_dir, role)
    manifest = load_json(app/'manifest.json')
    models = {}
    for row in read_tsv(output_dir/'library/final_motif_manifest.tsv'):
        identifier = row['final_motif_id']
        if not identifier or identifier in models:
            raise ValueError('Invalid frozen motif manifest')
        models[identifier] = {'family_id':row['family_id'],'rank':int(row['global_rank']),'width':int(row['width'])}
        if models[identifier]['rank']<1 or models[identifier]['width']<1:
            raise ValueError('Frozen motif rank/width is invalid: '+identifier)
    if len({m['rank'] for m in models.values()})!=len(models) or len(models)!=freeze['motif_count']:
        raise ValueError('Frozen motif ranks or count do not match manifest')
    shards = load_json(app/'shards/manifest.json')
    if len(scans)!=len(shards['shards']) or {Path(p).stem for p in scans}!={r['shard_id'] for r in shards['shards']}:
        raise ValueError(f'Incomplete {role} shard scans')
    fd, scratch_database = tempfile.mkstemp(prefix='annotation-', suffix='.sqlite')
    os.close(fd)
    db = sqlite3.connect(scratch_database)
    db.execute('PRAGMA cache_size=-65536')
    db.execute('PRAGMA temp_store=FILE')
    db.executescript('DROP TABLE IF EXISTS hits; DROP TABLE IF EXISTS positions; CREATE TABLE hits(window_id TEXT,motif_id TEXT,p TEXT,rank INTEGER,score TEXT,start INTEGER,stop INTEGER,strand TEXT,matched TEXT); CREATE INDEX by_window ON hits(window_id); CREATE TABLE positions(dataset_id TEXT,window_id TEXT,motif_id TEXT,bin INTEGER, PRIMARY KEY(dataset_id,window_id));')
    windows = sqlite3.connect('file:'+str((app/'shards/windows.sqlite').resolve()).replace('\\','/')+'?mode=ro',uri=True)
    loci = sqlite3.connect('file:'+str((output_dir/'library/selection_loci.sqlite').resolve()).replace('\\','/')+'?mode=ro',uri=True)
    try:
        temporary = app/'all_matches.tsv.gz.tmp'
        with gzip.open(temporary,'wt',encoding='utf-8',newline='') as audit:
            writer=csv.DictWriter(audit,fieldnames=MATCH_FIELDS,delimiter='\t',extrasaction='ignore')
            writer.writeheader()
            for path in scans:
                success=Path(str(path)+'.success.json')
                if not success.is_file():
                    raise ValueError('Missing successful FIMO completion: '+str(path))
                completion = load_json(success)
                expected_shard = next(r for r in shards['shards'] if r['shard_id']==Path(path).stem)
                if completion.get('sha256')!=sha256(path) or completion.get('models_sha256')!=sha256(output_dir/'library/final_motifs.meme') or completion.get('background_sha256')!=freeze['background_sha256'] or completion.get('fasta_sha256')!=sha256(expected_shard['fasta']) or completion.get('p_threshold')!=freeze['fimo_p_threshold'] or completion.get('status') not in ['complete','empty_models','empty_windows']:
                    raise ValueError('Invalid or stale FIMO completion: '+str(path))
                for row in iter_fimo(path):
                    motif=row['motif_id']; window=row['sequence_name']
                    if motif not in models:
                        raise ValueError('Unknown frozen motif: '+motif)
                    seq=windows.execute('SELECT sequence FROM windows WHERE window_id=?',(window,)).fetchone()
                    if seq is None:
                        raise ValueError('Unknown '+role+' window: '+window)
                    key=winner_key(row,models[motif]['rank'])
                    if key[4]>len(seq[0]):
                        raise ValueError('FIMO coordinates outside window: '+window)
                    writer.writerow(row)
                    db.execute('INSERT INTO hits VALUES(?,?,?,?,?,?,?,?,?)',(window,motif,str(key[0]),key[1],row['score'],key[3],key[4],row['strand'],row.get('matched_sequence','')))
            db.commit()
        os.replace(temporary,app/'all_matches.tsv.gz')
        counts=Counter()
        with gzip.open(app/'genotypes.tsv.gz.tmp','wt',encoding='utf-8',newline='') as all_handle, gzip.open(app/'alt_annotations.tsv.gz.tmp','wt',encoding='utf-8',newline='') as alt_handle:
            all_writer=csv.DictWriter(all_handle,fieldnames=FIELDS,delimiter='\t',extrasaction='ignore'); all_writer.writeheader()
            alt_writer=csv.DictWriter(alt_handle,fieldnames=FIELDS,delimiter='\t',extrasaction='ignore'); alt_writer.writeheader()
            for dataset in manifest['datasets']:
                for original in read_tsv(output_dir/'prepared'/dataset['dataset_id']/'genotypes.tsv.gz'):
                    row={key:original.get(key,'') for key in FIELDS}
                    row.update(role=dataset['role'],selection_participant=dataset['selection_participant'],
                               evidence_context=dataset['evidence_context'],
                               identity_check_status=dataset['identity_check_status'],distinct_motif_count=0)
                    row['seen_in_selection_loci']='true' if loci.execute('SELECT 1 FROM loci WHERE CHROM=? AND POS=? LIMIT 1',(row['CHROM'],int(row['POS']))).fetchone() else 'false'
                    gt=row['GT']; alleles=re.split(r'[/|]',gt)
                    alt=any(x!='.' and int(x)>0 for x in alleles)
                    if not alt:
                        row.update(annotation_status='not_applicable',annotation_reason='missing_gt' if all(x=='.' for x in alleles) else 'reference_only')
                    elif row['extraction_status']!='eligible':
                        row.update(annotation_status='not_scanned',annotation_reason=row['extraction_reason'])
                    elif not models:
                        row.update(annotation_status='unclassified',annotation_reason='no_selected_motifs')
                    else:
                        seq=windows.execute('SELECT sequence FROM windows WHERE window_id=?',(row['window_id'],)).fetchone()
                        if seq is None:
                            raise ValueError('Eligible genotype missing extracted window: '+row['record_id'])
                        if not re.search('[ACGT]{%d,}'%min(m['width'] for m in models.values()),seq[0]):
                            row.update(annotation_status='not_scanned',annotation_reason='no_scoreable_segment')
                        else:
                            best=None; best_key=None; distinct=set()
                            for motif,p,rank,score,start,stop,strand,matched in db.execute('SELECT motif_id,p,rank,score,start,stop,strand,matched FROM hits WHERE window_id=?',(row['window_id'],)):
                                distinct.add(motif)
                                hit={'motif_id':motif,'p-value':p,'score':score,'start':start,'stop':stop,'strand':strand,'matched_sequence':matched}
                                key=winner_key(hit,rank)
                                if best_key is None or key<best_key:
                                    best,best_key=hit,key
                            row['distinct_motif_count']=len(distinct)
                            row.update(annotation_status='promising' if best else 'unclassified',annotation_reason='')
                            if best:
                                midpoint=(best['start']+best['stop'])/2-(manifest['flank_bp']+1)
                                row.update(best_final_motif_id=best['motif_id'],best_family_id=models[best['motif_id']]['family_id'],score=best['score'],reported_p_value=best['p-value'],matched_sequence=best['matched_sequence'],strand=best['strand'],local_start=best['start'],local_stop=best['stop'],midpoint_distance=midpoint)
                                db.execute('INSERT OR IGNORE INTO positions VALUES(?,?,?,?)',(row['dataset_id'],row['window_id'],best['motif_id'],math.floor(midpoint)))
                    all_writer.writerow(row)
                    if alt: alt_writer.writerow(row)
                    counts[(row['dataset_id'],row['sample_name'],row['annotation_status'],row['annotation_reason'])]+=1
        for name in ['genotypes.tsv.gz','alt_annotations.tsv.gz']:
            os.replace(app/(name+'.tmp'),app/name)
        reports.mkdir(parents=True,exist_ok=True)
        write_tsv(reports/'annotation_counts.tsv',({'dataset_id':d,'sample_name':n,'annotation_status':s,'annotation_reason':r,'record_count':v} for (d,n,s,r),v in sorted(counts.items())),['dataset_id','sample_name','annotation_status','annotation_reason','record_count'])
        write_tsv(reports/'application_positions.tsv',({'dataset_id':d,'motif_id':m,'midpoint_bin':b,'window_count':v} for d,m,b,v in db.execute('SELECT dataset_id,motif_id,bin,count(*) FROM positions GROUP BY dataset_id,motif_id,bin ORDER BY dataset_id,motif_id,bin')),['dataset_id','motif_id','midpoint_bin','window_count'])
        atomic_json(app/'annotation.complete.json',{'schema_version':1,'role':role,'freeze_sha256':sha256(output_dir/'library/freeze.json'),'record_count':sum(counts.values()),'identity_check_status':manifest['identity_check_status']})
    finally:
        db.close(); windows.close(); loci.close()
        Path(scratch_database).unlink(missing_ok=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['manifest','corresponding-manifest','prepare-sample','shard','annotate','verify-freeze','fingerprint','fingerprint-corresponding','scan','report-complete'])
    parser.add_argument('--config',required=True)
    parser.add_argument('--role',choices=['illumina','other_illumina'],default='other_illumina')
    parser.add_argument('--dataset-id')
    parser.add_argument('--phase',choices=['records-plan','records-from-query','contexts-plan','contexts-from-fasta','all'],default='all')
    parser.add_argument('--scans',nargs='*',default=[])
    parser.add_argument('--fasta')
    parser.add_argument('--output')
    parser.add_argument('--log')
    parser.add_argument('--prepared-headers',action='store_true')
    parser.add_argument('--disk-mb', type=int)
    args=parser.parse_args()
    config=load_json(args.config); out=Path(config['output_dir'])
    if args.action=='fingerprint': return application_fingerprint(config,out)
    if args.action=='fingerprint-corresponding': return corresponding_fingerprint(config,out)
    disk_mb = args.disk_mb or {'scan':30000, 'annotate':60000, 'shard':30000, 'prepare-sample':30000}.get(args.action,10000)
    configure_scratch(config, args.action, args.dataset_id or Path(args.fasta or 'global').stem, disk_mb=disk_mb, scope='application')
    if args.action=='scan':
        result=scan_external(config,out,args.fasta,args.output,args.log)
        provenance(Path(args.output),[Path(args.fasta),out/'library/final_motifs.meme',out/'library/scoring_background.bfile',out/'library/freeze.json'],{'p_threshold':load_json(out/'library/freeze.json')['fimo_p_threshold'],'container':config.get('container',{}).get('meme'),'code_paths':[str(Path(__file__).resolve()),str(Path(__file__).with_name('evidence.py').resolve())]},[Path(args.output),Path(args.output+'.success.json')])
        return result
    if args.action=='manifest':
        identity_manifest(config,out,args.prepared_headers)
        app=out/'application'
        inputs=[out/'library/freeze.json',app/'input_fingerprints.json',Path(config['other_illumina'])]
        outputs=[app/'manifest.json',app/'datasets.tsv',app/'application_config.json',out/'inputs/other_illumina.vcf.gz',out/'inputs/other_illumina.vcf.gz.csi',app/'header.samples.txt',app/'header.vcf.txt']
        provenance(app/'manifest.json',inputs,{'metadata':config.get('sample_metadata',{}).get('other_illumina',{}),'container':config.get('container',{}).get('bcftools')},outputs)
    elif args.action=='corresponding-manifest':
        manifest=corresponding_manifest(config,out)
        app,_=annotation_paths(out,'illumina')
        fingerprint=load_json(app/'input_fingerprints.json')
        inputs=[out/'library/freeze.json',out/'manifests/selection_fingerprints.json']+[Path(value['path']) for value in fingerprint['inputs'].values()]
        outputs=[app/'manifest.json',app/'datasets.tsv']
        provenance(app/'manifest.json',inputs+[app/'input_fingerprints.json'],{'role':'illumina','selection_participant':True,'evidence_context':'selection_reuse','fingerprint_signature':fingerprint['signature'],'code_paths':[str(Path(__file__).with_name('prepare.py').resolve())]},outputs)
    elif args.action=='prepare-sample': prepare_external(config,out,args.dataset_id,args.phase)
    elif args.action=='shard':
        shard_windows(out,config['execution']['scan_windows_per_shard'],role=args.role)
        app,_=annotation_paths(out,args.role); manifest=load_json(app/'manifest.json'); shards=load_json(app/'shards/manifest.json')
        provenance(app/'shards/manifest.json',[app/'manifest.json']+[out/'prepared'/r['dataset_id']/'all.fasta' for r in manifest['datasets']],{'role':args.role,'maximum':config['execution']['scan_windows_per_shard'],'config_paths':['execution.scan_windows_per_shard'],'config_values':{'execution.scan_windows_per_shard':config['execution']['scan_windows_per_shard']}},[app/'shards/manifest.json',app/'shards/windows.sqlite']+[Path(r['fasta']) for r in shards['shards']])
    elif args.action=='annotate':
        annotations(out,[Path(x) for x in args.scans],role=args.role)
        app,reports=annotation_paths(out,args.role); manifest=load_json(app/'manifest.json')
        inputs=[out/'library/freeze.json',out/'library/final_motif_manifest.tsv',out/'library/selection_loci.sqlite',app/'manifest.json',app/'shards/windows.sqlite']+[Path(x) for x in args.scans]+[Path(x+'.success.json') for x in args.scans]+[out/'prepared'/r['dataset_id']/'genotypes.tsv.gz' for r in manifest['datasets']]
        outputs=[app/name for name in ['annotation.complete.json','genotypes.tsv.gz','alt_annotations.tsv.gz','all_matches.tsv.gz']]+[reports/'annotation_counts.tsv',reports/'application_positions.tsv']
        provenance(app/'annotation.complete.json',inputs,{'role':args.role,'policy':'decimal_p_global_rank_score_v1'},outputs)
    elif args.action=='report-complete':
        reports=out/'reports'
        outputs=[reports/name for name in ['motif_selection.pdf','conditional_positions.pdf','README.md']]
        inputs=[reports/(name+'.tsv') for name in ['annotation_counts','motif_coverage','selection_decisions','representative_validation','conditional_position_bins']]+[out/'library/freeze.json']
        atomic_json(reports/'report.complete.json',{'schema_version':1,'outputs':{str(p):sha256(p) for p in outputs}})
        provenance(reports/'report.complete.json',inputs,{'code_paths':[str(Path(__file__).with_name('report.R').resolve())],'container':config.get('container',{}).get('r')},outputs+[reports/'report.complete.json'])
    else: verify_freeze(out)

if __name__=='__main__':
    main()
