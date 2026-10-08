"""Discovery, disk-backed held-out evidence and immutable exact-matrix selection."""
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
import sqlite3
import statistics
import subprocess
import tempfile
from decimal import Decimal, InvalidOperation
import xml.etree.ElementTree as ET

from prepare import canonical, sha256, atomic_json as _atomic_json, read_tsv, write_tsv, dataset_id, configure_scratch

FIMO_FIELDS = ['motif_id','motif_alt_id','sequence_name','start','stop','strand','score','p-value','q-value','matched_sequence']
METRIC_FIELDS = ['dataset_id','rotation','motif_id','n_windows','real_hits','control_hits','real_coverage','control_coverage','coverage_difference','coverage_ratio','real_central_coverage','control_central_coverage','real_exact_center_coverage','control_exact_center_coverage','real_central_proportion','control_central_proportion','central_difference','status']

def atomic_json(path,value):
    def safe(v):
        if isinstance(v,float) and not math.isfinite(v):return 'Inf' if v>0 else '-Inf'
        if isinstance(v,dict):return {k:safe(x) for k,x in v.items()}
        if isinstance(v,(tuple,list)):return [safe(x) for x in v]
        return v
    return _atomic_json(path,safe(value))


def read_fasta(path):
    name, chunks = None, []
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                if name is None:
                    raise ValueError(f'Invalid FASTA: {path}')
                chunks.append(line.upper())
    if name is not None:
        yield name, ''.join(chunks)


def write_fasta(path, records):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    with open(tmp, 'w', encoding='utf-8', newline='\n') as handle:
        for name, seq in records:
            handle.write(f'>{name}\n{seq}\n')
    os.replace(tmp, path)


def parse_meme(path):
    """Preserve original model order and exact probability text/nsites."""
    lines = Path(path).read_text(encoding='utf-8').splitlines()
    models = []
    for i, line in enumerate(lines):
        if not line.startswith('MOTIF '):
            continue
        parts = line.split()
        model = {'id': parts[1], 'original_id': parts[1], 'consensus': parts[2] if len(parts)>2 else '', 'order': len(models)+1}
        j = i+1
        while j < len(lines) and not lines[j].startswith('letter-probability matrix:'):
            if lines[j].startswith('MOTIF '):
                raise ValueError(f'Missing matrix for {parts[1]} in {path}')
            j += 1
        if j == len(lines):
            raise ValueError(f'Missing matrix for {parts[1]} in {path}')
        attrs = dict(re.findall(r'(\w+)\s*=\s*([^\s]+)', lines[j]))
        width = int(attrs['w']); matrix = []
        for row in lines[j+1:j+1+width]:
            values = row.split()
            if len(values)!=4 or any(not Decimal(v).is_finite() or Decimal(v)<0 for v in values) or abs(sum(Decimal(v) for v in values)-1)>Decimal('0.002'):
                raise ValueError(f'Invalid matrix {parts[1]} in {path}')
            matrix.append(values)
        if len(matrix)!=width or width<1:
            raise ValueError(f'Truncated matrix {parts[1]} in {path}')
        e = attrs.get('E')
        try:
            e = float(e)
            if not math.isfinite(e) or e<0: e=None
        except (TypeError, ValueError): e=None
        model.update(width=width, nsites=attrs.get('nsites','1'), evalue=e, matrix=matrix)
        model['matrix_sha256'] = hashlib.sha256(canonical([width,model['nsites'],matrix])).hexdigest()
        models.append(model)
    if len({m['id'] for m in models}) != len(models):
        raise ValueError(f'Duplicate motif IDs in {path}')
    return models


def background_frequencies(path):
    text = Path(path).read_text(encoding='utf-8')
    vals = dict(re.findall(r'\b([ACGT])\s+([0-9.eE+-]+)', text))
    if set(vals)!=set('ACGT') or any(not math.isfinite(float(v)) or float(v)<=0 for v in vals.values()) or abs(sum(float(v) for v in vals.values())-1)>0.001:
        raise ValueError(f'Invalid zero-order background: {path}')
    return vals

def make_background(reference,output,log):
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    tmp=output.with_name(output.name+'.tmp')
    run_tool(['fasta-get-markov','-m','0','-pseudo','0.1',reference,tmp],log)
    background_frequencies(tmp)
    os.replace(tmp,output)


def write_meme(path, models, background):
    bg = background_frequencies(background)
    path = Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp = path.with_name(path.name+'.tmp')
    with open(tmp,'w',encoding='utf-8',newline='\n') as h:
        h.write('MEME version 5\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\n')
        h.write(' '.join(f'{b} {bg[b]}' for b in 'ACGT')+'\n')
        for m in models:
            efield=f" E= {m['evalue']}" if m.get('evalue') is not None else ''
            h.write(f"\nMOTIF {m['id']} {m.get('consensus','')}\nletter-probability matrix: alength= 4 w= {m['width']} nsites= {m['nsites']}{efield}\n")
            for row in m['matrix']: h.write(' '.join(row)+'\n')
    os.replace(tmp,path)


def run_tool(args, log, stdout=None):
    Path(log).parent.mkdir(parents=True,exist_ok=True)
    with open(log,'wb') as err:
        if stdout is None:
            result = subprocess.run([str(x) for x in args],stderr=err,stdout=subprocess.PIPE,check=False)
        else:
            with open(stdout,'wb') as out:
                result = subprocess.run([str(x) for x in args],stderr=err,stdout=out,check=False)
    if result.returncode:
        raise RuntimeError(f'Tool exited {result.returncode}: {args}; stderr: {log}')
    return result.stdout


def shuffle_fasta(source, output, kmer, seed, log):
    output=Path(output); output.parent.mkdir(parents=True,exist_ok=True)
    original=list(read_fasta(source))
    if not original:
        write_fasta(output,[]); atomic_json(str(output)+'.map.json',[]); return
    tmp=output.with_name(output.name+'.shuffle.tmp')
    run_tool(['fasta-shuffle-letters','-dna','-kmer',kmer,'-seed',seed,source,tmp],log)
    shuffled=list(read_fasta(tmp))
    if len(shuffled)!=len(original) or any(len(a[1])!=len(b[1]) for a,b in zip(original,shuffled)):
        raise ValueError('Shuffle sequence count/length mismatch')
    # MEME shuffler preserves order but changes headers; map explicitly, never infer suffixes.
    write_fasta(output,((a[0],b[1]) for a,b in zip(original,shuffled)))
    atomic_json(str(output)+'.map.json',[{'window_id':a[0],'shuffled_id':b[0]} for a,b in zip(original,shuffled)])
    tmp.unlink()


def scan_fimo(models, fasta, background, output, p_threshold, log, flag='--bgfile'):
    output=Path(output); output.parent.mkdir(parents=True,exist_ok=True)
    model_set=parse_meme(models); windows={n:len(s) for n,s in read_fasta(fasta)}
    tmp=output.with_name(output.name+'.scan.tmp')
    if not model_set or not windows:
        tmp.write_text('\t'.join(FIMO_FIELDS)+'\n',encoding='utf-8')
        Path(log).parent.mkdir(parents=True,exist_ok=True); Path(log).write_text('Empty models or FASTA; FIMO not invoked.\n')
        status='empty_models' if not model_set else 'empty_windows'
    else:
        run_tool(['fimo','--text','--no-pgc',flag,background,'--motif-pseudo','0.1','--thresh',p_threshold,models,fasta],log,tmp)
        # MEME 5.5.7 emits no header when a successful text scan has zero matches.
        if tmp.stat().st_size == 0:
            tmp.write_text('\t'.join(FIMO_FIELDS)+'\n',encoding='utf-8')
        for row in iter_fimo(tmp,{m['id'] for m in model_set},windows):
            if row['stop']>windows[row['sequence_name']]: raise ValueError('FIMO coordinates exceed window')
        status='complete'
    os.replace(tmp,output)
    atomic_json(str(output)+'.success.json',{'status':status,'sha256':sha256(output),'models_sha256':sha256(models),'fasta_sha256':sha256(fasta),'background_sha256':sha256(background),'p_threshold':p_threshold,'background_flag':flag,'pseudocount':0.1,'strands':'+ -'})


def iter_fimo(path,motif_ids=None,window_ids=None):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt',encoding='utf-8',newline='') as h:
        reader=csv.DictReader((l for l in h if l.strip() and not l.startswith('#')),delimiter='\t')
        required={'motif_id','sequence_name','start','stop','strand','score','p-value','matched_sequence'}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f'Malformed FIMO header: {path}')
        for number,row in enumerate(reader,2):
            try:
                if None in row or any(row[k] is None for k in required): raise ValueError('truncated row')
                row['start'],row['stop']=int(row['start']),int(row['stop'])
                score,p=Decimal(row['score']),Decimal(row['p-value'])
                if row['start']<1 or row['stop']<row['start'] or row['strand'] not in ['+','-'] or not score.is_finite() or not p.is_finite() or not 0<=p<=1: raise ValueError('invalid occurrence')
                if motif_ids is not None and row['motif_id'] not in motif_ids: raise ValueError('unknown motif')
                if window_ids is not None and row['sequence_name'] not in window_ids: raise ValueError('unknown window')
            except (ValueError,InvalidOperation,TypeError) as exc:
                raise ValueError(f'Malformed FIMO {path}:{number}: {exc}') from exc
            yield row


def metric(real, control, n, flank, radius):
    def count(rows):
        hits=central=exact=0
        for row in rows:
            hits+=1; start,stop=int(row['start']),int(row['stop'])
            central+=abs((start+stop)/2-(flank+1))<=radius
            exact+=start<=flank+1<=stop
        return hits,central,exact
    rh,rc,re=count(real); ch,cc,ce=count(control)
    cov=lambda v: v/n if n else None
    prop=lambda a,b:a/b if b else None
    r,c=cov(rh),cov(ch); rp,cp=prop(rc,rh),prop(cc,ch)
    return dict(n_windows=n,real_hits=rh,control_hits=ch,real_coverage=r,control_coverage=c,coverage_difference=r-c if n else None,coverage_ratio=(r/c if c else math.inf if r else None) if n else None,real_central_coverage=cov(rc),control_central_coverage=cov(cc),real_exact_center_coverage=cov(re),control_exact_center_coverage=cov(ce),real_central_proportion=rp,control_central_proportion=cp,central_difference=rp-cp if rp is not None and cp is not None else None,status='complete' if n else 'insufficient_evaluation_windows')


def number(value):
    return None if value in ('',None,'NA','None') else float(value)


def support(metrics,ratio_min,central=False,evalue=None,e_max=None):
    ratio,diff,cd=(number(metrics.get(k)) for k in ['coverage_ratio','coverage_difference','central_difference'])
    return bool(ratio is not None and ratio>ratio_min and diff is not None and diff>0 and (not central or cd is not None and cd>0) and (e_max is None or evalue is not None and evalue<e_max))


def display_rows(rows):
    for row in rows:
        yield {k:('NA' if v is None or isinstance(v,float) and v==-math.inf else 'Inf' if isinstance(v,float) and v==math.inf else v) for k,v in row.items()}


def open_evidence_db(path):
    db=sqlite3.connect(path)
    db.execute('PRAGMA cache_size=-65536'); db.execute('PRAGMA temp_store=FILE')
    db.execute('CREATE TABLE IF NOT EXISTS hits (kind TEXT,motif TEXT,window TEXT,start INTEGER,stop INTEGER,strand TEXT,score TEXT,p TEXT,matched TEXT, PRIMARY KEY(kind,motif,window))')
    return db


def positional_key(row):
    return (-Decimal(str(row['score'])),Decimal(str(row.get('p-value',row.get('p')))),int(row['start']),int(row['stop']),0 if row['strand']=='+' else 1)


def summarize_scans(real,control,models,fasta,output,flank,radius,dataset,rotation,db_path=None):
    mids={m['id'] for m in parse_meme(models)}; windows={n for n,_ in read_fasta(fasta)}
    output=Path(output); output.parent.mkdir(parents=True,exist_ok=True)
    db_path=Path(db_path or str(output)+'.sqlite.tmp')
    if db_path.exists(): db_path.unlink()
    db=open_evidence_db(db_path)
    for kind,path in [('real',real),('control',control)]:
        success=Path(str(path)+'.success.json')
        if not success.exists() or json.loads(success.read_text())['sha256']!=sha256(path): raise ValueError(f'Missing valid FIMO producer completion: {path}')
        audit=Path(str(output)+f'.{kind}.occurrences.tsv.gz')
        with gzip.open(str(audit)+'.tmp','wt',encoding='utf-8',newline='') as out:
            writer=csv.DictWriter(out,fieldnames=FIMO_FIELDS,delimiter='\t',extrasaction='ignore'); writer.writeheader()
            for i,row in enumerate(iter_fimo(path,mids,windows),1):
                writer.writerow(row)
                old=db.execute('SELECT start,stop,strand,score,p FROM hits WHERE kind=? AND motif=? AND window=?',(kind,row['motif_id'],row['sequence_name'])).fetchone()
                if old is None or positional_key(row)<positional_key(dict(zip(['start','stop','strand','score','p'],old))):
                    db.execute('INSERT OR REPLACE INTO hits VALUES (?,?,?,?,?,?,?,?,?)',(kind,row['motif_id'],row['sequence_name'],row['start'],row['stop'],row['strand'],row['score'],row['p-value'],row['matched_sequence']))
                if i%10000==0: db.commit()
        os.replace(str(audit)+'.tmp',audit)
    db.commit()
    rows=[]
    for mid in sorted(mids):
        def hits(kind):
            return ({'start':a,'stop':b} for a,b in db.execute('SELECT start,stop FROM hits WHERE kind=? AND motif=?',(kind,mid)))
        rows.append(dict(dataset_id=dataset,rotation=rotation,motif_id=mid,**metric(hits('real'),hits('control'),len(windows),flank,radius)))
    write_tsv(output,display_rows(rows),METRIC_FIELDS)
    db.close()
    if db_path.parent==output.parent:
        os.replace(db_path,str(output)+'.sqlite')
    else:
        staged=str(output)+'.sqlite.tmp'
        shutil.copyfile(db_path,staged);os.replace(staged,str(output)+'.sqlite');db_path.unlink()
    return rows


def tomtom_edges(path,ids,q_max,min_overlap):
    directed={}
    for row in read_tsv(path):
        if row.get('Query_ID','').startswith('#') or not row.get('Target_ID'): continue
        a,b=row['Query_ID'],row['Target_ID']
        if a not in ids or b not in ids: raise ValueError('Tomtom reports unknown model')
        q=Decimal(row['q-value']); overlap=int(row['Overlap'])
        if not q.is_finite() or not 0<=q<=1 or overlap<0: raise ValueError('Invalid Tomtom match')
        if q<=Decimal(str(q_max)) and overlap>=min_overlap: directed[a,b]=row
    return {frozenset((a,b)) for a,b in directed if a!=b and (b,a) in directed}


def complete_link_families(ids,edges):
    families=[]
    # ponytail: deterministic complete-link is quadratic; index adjacency if model pools make it expensive.
    for mid in sorted(ids):
        for members in families:
            if all(frozenset((mid,other)) in edges for other in members):
                members.append(mid); break
        else: families.append([mid])
    result=[]; seen=set()
    for members in families:
        fid='fam_'+hashlib.sha256(canonical(members)).hexdigest()[:16]
        if fid in seen: raise ValueError('Family ID collision')
        seen.add(fid); result.append({'family_id':fid,'members':members})
    return result


def corresponding_support(rows,count,fraction,ratio_min,central=False):
    required=math.ceil(fraction*count)
    successes=sum(support(row,ratio_min,central) for row in rows)
    return {'supporting_count':successes,'required_count':required,'configured_count':count,'support':len(rows)==count and successes>=required}


def family_stability(rows,pacbio_ids,min_rotations,require_all=True):
    seen=set()
    for row in rows:
        key=(row['dataset_id'],int(row['rotation']))
        if key in seen or key[0] not in pacbio_ids or key[1] not in range(1,6):raise ValueError('Invalid or duplicate family clone/rotation row')
        seen.add(key)
    counts={pb:sum(bool(r['support']) for r in rows if r['dataset_id']==pb) for pb in pacbio_ids}
    stable=[pb for pb in pacbio_ids if counts[pb]>=min_rotations]
    passed=len(stable)==len(pacbio_ids) if require_all else bool(stable)
    return {'pass':passed,'counts':counts,'stable_supported':stable,'required_pacbio_ids':list(pacbio_ids) if require_all else stable}


def validate_candidate(rows,required_pacbio,correspondences,selection):
    by={r['dataset_id']:r for r in rows}; groups={}
    if len(by)!=len(rows):raise ValueError('Duplicate representative evaluated-dataset rows')
    central=selection['require_central_enrichment']
    for pb in required_pacbio:
        names=correspondences[pb]
        transfer=corresponding_support([by[i] for i in names if i in by],len(names),selection['corresponding_dataset_support_fraction'],selection['illumina_coverage_ratio_min'],central)
        groups[pb]={**transfer,'pacbio_support':pb in by and support(by[pb],selection['pacbio_coverage_ratio_min'],central)}
        groups[pb]['pass']=groups[pb]['pacbio_support'] and transfer['support']
    return {'pass':bool(required_pacbio) and all(g['pass'] for g in groups.values()),'groups':groups}


def median(values):
    v=sorted(x for x in values if x is not None)
    if not v:return None
    n=len(v)
    return v[n//2] if n%2 or math.isinf(v[n//2]) else (v[n//2-1]+v[n//2])/2


def representative_rank(model,rows,required,correspondences,selection):
    used=[r for r in rows if r['dataset_id'] in required or any(r['dataset_id'] in correspondences[p] for p in required) and support(r,selection['illumina_coverage_ratio_min'],selection['require_central_enrichment'])]
    ratios=[number(r['coverage_ratio']) for r in used]; diffs=[number(r['coverage_difference']) for r in used]
    source=next(r for r in rows if r['dataset_id']==model['dataset_id'])
    return (-min(ratios),-median(ratios),-min(diffs),-number(source['real_coverage']),model['evalue'],model['id'])


def publish_provenance(primary,inputs,params,outputs,scope='selection'):
    from prepare import publish_provenance as publish
    selection_actions={'pool','tomtom','families','stability','representatives','validate-representatives','rank','freeze'}
    rule_file='selection.smk' if params.get('action') in selection_actions or params.get('member') else 'transfer.smk' if params.get('action')=='transfer' or params.get('source') else 'discovery.smk'
    params={**params,'code_paths':sorted(set(params.get('code_paths',[])) | {str(Path(__file__).resolve()),str(Path(__file__).resolve().parents[1]/'rules'/rule_file)})}
    return publish(primary,list(inputs.values()) if isinstance(inputs,dict) else inputs,params,outputs,scope)


def context(config):
    out=Path(config['output_dir']); analysis=config['analysis']; selection=config['selection']
    pb={dataset_id('pacbio',n):n for n in config['pacbio_samples']}
    corr={p:[dataset_id('illumina',n) for n in config['illumina_corresponding'][name]] for p,name in pb.items()}
    return out,analysis,selection,pb,corr


def seed(config,dataset,rotation,role):
    out,_,_,pb,corr=context(config)
    for row in read_tsv(out/'manifests/seeds.tsv'):
        name=pb.get(dataset)
        if name is None:
            name=next(n for values in config['illumina_corresponding'].values() for n in values if dataset_id('illumina',n)==dataset)
        if row.get('sample_name',row.get('sample'))==name and int(row['rotation'])==rotation and row['role']==role and row.get('platform',row.get('source_role','pacbio' if dataset in pb else 'illumina'))==('pacbio' if dataset in pb else 'illumina'):
            return int(row['seed'])
    raise ValueError(f'Missing seed: {dataset}/{rotation}/{role}')


def evdir(out,dataset,rotation):return out/'evidence'/dataset/f'r{rotation:02d}'


def make_train_test(config,dataset,rotation):
    out,*_=context(config); dest=evdir(out,dataset,rotation)
    prepared=out/'prepared'/dataset
    for kind,folds in [('train',[r for r in range(1,6) if r!=rotation]),('test',[rotation])]:
        write_fasta(dest/f'{kind}.fasta',(record for r in folds for record in sorted(read_fasta(prepared/f'fold{r}.fasta'))))
    atomic_json(dest/'split.json',{'test_fold':rotation,'train_folds':[r for r in range(1,6) if r!=rotation],'dataset_id':dataset})


def discover(config,dataset,rotation):
    out,a,_,_,_=context(config); dest=evdir(out,dataset,rotation)
    records=list(read_fasta(dest/'train.fasta'))
    dest.mkdir(parents=True,exist_ok=True)
    if len(records)<2:
        (dest/'streme.txt').write_text('MEME version 5\n',encoding='utf-8')
        (dest/'streme.xml').write_text('<streme status="insufficient_discovery_sequences"/>',encoding='utf-8')
        atomic_json(dest/'streme.complete.json',{'status':'insufficient_discovery_sequences','models':[]});return
    scratch=Path(os.environ.get('TMPDIR', config['scratch_dir'])); scratch.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='streme-',dir=scratch) as temp:
        target=Path(temp)/'streme_out'
        run_tool(['streme','--dna','--p',dest/'train.fasta','--n',dest/'train.control.fasta','--minw',a['streme_min_width'],'--maxw',a['streme_max_width'],'--seed',seed(config,dataset,rotation,'streme'),'--totallength',a['streme_total_length'],'--o',target],dest/'streme.stderr.log')
        for name in ['streme.txt','streme.xml']:
            if not (target/name).exists():raise ValueError(f'STREME missing {name}')
            shutil.copyfile(target/name,dest/name)
    models=parse_meme(dest/'streme.txt')
    for m in models:
        m.update(dataset_id=dataset,rotation=rotation,id=f"{dataset}__r{rotation:02d}__m{m['order']:04d}")
    xml=ET.parse(dest/'streme.xml')
    metadata={node.tag:node.attrib for node in xml.iter() if node.tag in ['train_positives','train_negatives','test_positives','test_negatives','stop','reason_for_stopping']}
    atomic_json(dest/'streme.complete.json',{'status':'complete','models':models,'usage':metadata,'input_sequences':len(records),'input_bases':sum(len(s) for _,s in records),'cap':a['streme_total_length']})


def parse_models_job(config,dataset,rotation):
    out,a,_,_,_=context(config); dest=evdir(out,dataset,rotation)
    data=json.loads((dest/'streme.complete.json').read_text())
    models=data['models']; write_meme(dest/'models.meme',models,out/'reference/scoring_background.bfile')
    atomic_json(dest/'models.json',data)


def summarize_job(config,dataset,rotation,source=None,member=None):
    out,a,_,_,_=context(config)
    if member:
        dest=out/'selection/representatives'/member/dataset
        base=dest/f'r{rotation:02d}'; models=out/'selection/representative_models'/member/'model.meme'
    else:
        dest=evdir(out,source or dataset,rotation)
        base=dest/('pacbio' if source is None else dataset); models=dest/('models.meme' if source is None else 'candidates.meme')
    fasta=out/'prepared'/dataset/f'fold{rotation}.fasta'
    scratch=Path(os.environ.get('TMPDIR', config['scratch_dir']))
    scratch.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='attempt-',dir=scratch) as attempt:
        rows=summarize_scans(str(base)+'.real.tsv',str(base)+'.control.tsv',models,fasta,str(base)+'.metrics.tsv',a['flank_bp'],a['center_radius_bp'],dataset,rotation,Path(attempt)/'hits.sqlite')
    if source:
        db=sqlite3.connect(str(base)+'.metrics.tsv.sqlite')
        db.execute('PRAGMA cache_size=-65536');db.execute('PRAGMA temp_store=FILE')
        db.execute('CREATE TABLE loci (chrom TEXT,pos INTEGER,PRIMARY KEY(chrom,pos))')
        for row in read_tsv(out/'prepared'/source/'genotypes.tsv.gz'):
            if any(t.isdigit() and int(t)>0 for t in re.split(r'[/|]',row['GT'])):
                db.execute('INSERT OR IGNORE INTO loci VALUES (?,?)',(row['CHROM'],int(row['POS'])))
        db.execute('CREATE TABLE subsets (window TEXT PRIMARY KEY, subset TEXT)')
        eligible={n for n,_ in read_fasta(fasta)}
        for row in read_tsv(out/'prepared'/dataset/'windows.tsv'):
            if row['window_id'] in eligible:
                shared=db.execute('SELECT 1 FROM loci WHERE chrom=? AND pos=?',(row['CHROM'],int(row['POS']))).fetchone()
                db.execute('INSERT OR IGNORE INTO subsets VALUES (?,?)',(row['window_id'],'shared' if shared else 'additional'))
        result=[]
        for subset in ['shared','additional']:
            n=db.execute('SELECT COUNT(*) FROM subsets WHERE subset=?',(subset,)).fetchone()[0]
            for m in parse_meme(models):
                def hits(kind):
                    return ({'start':start,'stop':stop} for start,stop in db.execute('SELECT h.start,h.stop FROM hits h JOIN subsets s ON h.window=s.window WHERE h.kind=? AND h.motif=? AND s.subset=?',(kind,m['id'],subset)))
                result.append(dict(dataset_id=dataset,rotation=rotation,motif_id=m['id'],subset=subset,**metric(hits('real'),hits('control'),n,a['flank_bp'],a['center_radius_bp'])))
        write_tsv(str(base)+'.metrics.tsv.subsets.tsv',display_rows(result),METRIC_FIELDS+['subset'])
        db.commit();db.close()
    return rows


def candidates_job(config,dataset,rotation):
    out,a,s,_,_=context(config); dest=evdir(out,dataset,rotation)
    models=json.loads((dest/'models.json').read_text())['models']; rows={r['motif_id']:r for r in read_tsv(dest/'pacbio.metrics.tsv')}
    selected=[]; decisions=[]
    for m in models:
        ok=support(rows[m['id']],s['pacbio_coverage_ratio_min'],s['require_central_enrichment'],m['evalue'],a['streme_e_max'])
        decisions.append({'motif_id':m['id'],'eligible':ok,'reason':'accepted' if ok else 'unavailable_evalue' if m['evalue'] is None else 'heldout_support_or_evalue_failed'})
        if ok:selected.append(m)
    write_meme(dest/'candidates.meme',selected,out/'reference/scoring_background.bfile')
    atomic_json(dest/'candidates.json',{'models':selected,'decisions':decisions,'dataset_id':dataset,'rotation':rotation})


def transfer_job(config,dataset,rotation):
    out,_,s,_,corr=context(config); dest=evdir(out,dataset,rotation)
    models=json.loads((dest/'candidates.json').read_text())['models']; metrics={i:{r['motif_id']:r for r in read_tsv(dest/f'{i}.metrics.tsv')} for i in corr[dataset]}
    chosen=[]; decisions=[]; rows=[]
    for m in models:
        evidence=[metrics[i][m['id']] for i in corr[dataset]]
        group=corresponding_support(evidence,len(corr[dataset]),s['corresponding_dataset_support_fraction'],s['illumina_coverage_ratio_min'],s['require_central_enrichment'])
        decisions.append({'motif_id':m['id'],**group})
        for row in evidence:rows.append({**row,'support':support(row,s['illumina_coverage_ratio_min'],s['require_central_enrichment']),**group,'pacbio_dataset_id':dataset})
        if group['support']:chosen.append(m)
    write_tsv(dest/'transfer.tsv',rows,METRIC_FIELDS+['support','supporting_count','required_count','configured_count','pacbio_dataset_id'])
    atomic_json(dest/'transfer.complete.json',{'models':chosen,'decisions':decisions,'dataset_id':dataset,'rotation':rotation})


def pool_job(config):
    out,_,_,pb,_=context(config); models=[]
    for p in pb:
        for r in range(1,6):models.extend(json.loads((evdir(out,p,r)/'transfer.complete.json').read_text())['models'])
    atomic_json(out/'selection/pool.json',{'models':models})
    write_meme(out/'selection/transferable_models.meme',models,out/'reference/scoring_background.bfile')


def tomtom_job(config):
    out,_,s,_,_=context(config); dest=out/'selection'; models=json.loads((dest/'pool.json').read_text())['models']
    target=dest/'tomtom.tsv'
    if not models:
        write_tsv(target,[],['Query_ID','Target_ID','q-value','Overlap','Orientation','Optimal_offset']);return
    scratch=Path(os.environ.get('TMPDIR', config['scratch_dir']));scratch.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='tomtom-',dir=scratch) as temp:
        result=Path(temp)/'results'
        run_tool(['tomtom','-oc',result,'-dist','ed','-min-overlap',s['family_matching']['min_overlap_bp'],'-thresh',s['family_matching']['q_max'],dest/'transferable_models.meme',dest/'transferable_models.meme'],dest/'tomtom.stderr.log')
        shutil.copyfile(result/'tomtom.tsv',target)


def families_job(config):
    out,_,s,_,_=context(config); models=json.loads((out/'selection/pool.json').read_text())['models']; ids={m['id'] for m in models}
    edges=tomtom_edges(out/'selection/tomtom.tsv',ids,s['family_matching']['q_max'],s['family_matching']['min_overlap_bp'])
    families=complete_link_families(ids,edges)
    atomic_json(out/'selection/families.json',{'families':families,'models':models,'edges':[sorted(e) for e in sorted(edges,key=lambda e:sorted(e))]})
    write_tsv(out/'reports/family_membership.tsv',({'family_id':f['family_id'],'motif_id':m} for f in families for m in f['members']),['family_id','motif_id'])


def union_rows(config,family,pb,rotation,models):
    out,a,s,_,corr=context(config); dest=evdir(out,pb,rotation)
    members=[m for m in models if m['id'] in family['members'] and m['dataset_id']==pb and m['rotation']==rotation]
    rows=[]
    for ds in [pb]+corr[pb]:
        windows={n for n,_ in read_fasta(out/'prepared'/ds/f'fold{rotation}.fasta')}
        unions={'real':{},'control':{}}
        base=dest/('pacbio' if ds==pb else ds)
        if members:
            db=sqlite3.connect(str(base)+'.metrics.tsv.sqlite')
            for kind in unions:
                for m in members:
                    for window,start,stop,strand,score,p,matched in db.execute('SELECT window,start,stop,strand,score,p,matched FROM hits WHERE kind=? AND motif=?',(kind,m['id'])):
                        value={'start':start,'stop':stop,'strand':strand,'score':score,'p':p,'motif_id':m['id']}
                        key=(Decimal(p),m['id'],-Decimal(score),start,stop,0 if strand=='+' else 1)
                        if window not in unions[kind] or key<unions[kind][window][0]:unions[kind][window]=(key,value)
            db.close()
        row=dict(dataset_id=ds,pacbio_dataset_id=pb,rotation=rotation,motif_id=family['family_id'],family_id=family['family_id'],evalue=None,**metric((v[1] for v in unions['real'].values()),(v[1] for v in unions['control'].values()),len(windows),a['flank_bp'],a['center_radius_bp']))
        row['support']=bool(members) and support(row,s['pacbio_coverage_ratio_min'] if ds==pb else s['illumina_coverage_ratio_min'],s['require_central_enrichment'])
        rows.append(row)
    transfer=corresponding_support(rows[1:],len(corr[pb]),s['corresponding_dataset_support_fraction'],s['illumina_coverage_ratio_min'],s['require_central_enrichment'])
    return rows,{'dataset_id':pb,'rotation':rotation,'support':rows[0]['support'] and transfer['support'],'members':[m['id'] for m in members],**{k:v for k,v in transfer.items() if k!='support'}}


def stability_job(config):
    out,_,s,pb,_=context(config); data=json.loads((out/'selection/families.json').read_text()); all_rows=[]; decisions=[]
    for f in data['families']:
        groups=[]
        for p in pb:
            for r in range(1,6):
                rows,group=union_rows(config,f,p,r,data['models']); all_rows.extend(rows);groups.append(group)
        stable=family_stability(groups,list(pb),s['min_supported_rotations'],s['require_all_pacbio_clones'])
        decisions.append({**f,**stable,'groups':groups})
    write_tsv(out/'reports/family_coverage.tsv',display_rows(all_rows),METRIC_FIELDS+['family_id','pacbio_dataset_id','support','evalue'])
    atomic_json(out/'selection/stability.json',{'families':decisions,'models':data['models']})


def representatives_job(config):
    out,_,_,_,corr=context(config); data=json.loads((out/'selection/stability.json').read_text()); models={m['id']:m for m in data['models']}; candidates=[]
    model_root=out/'selection/representative_models';model_root.mkdir(parents=True,exist_ok=True)
    for f in data['families']:
        if not f['pass']:continue
        contributing={mid for g in f['groups'] if g['support'] for mid in g['members']}
        for mid in sorted(contributing):
            m=models[mid]
            if m['dataset_id'] not in f['required_pacbio_ids']:continue
            evaluated=[ds for pb in f['required_pacbio_ids'] for ds in [pb]+corr[pb]]
            entry={**m,'family_id':f['family_id'],'member_id':mid,'required_pacbio_ids':f['required_pacbio_ids'],'evaluated_dataset_ids':evaluated,'source_rotation':m['rotation']}
            candidates.append(entry)
            write_meme(model_root/mid/'model.meme',[m],out/'reference/scoring_background.bfile')
    atomic_json(out/'selection/representatives/manifest.json',{'candidates':candidates})


def validate_representatives_job(config):
    out,_,s,pb,corr=context(config); manifest=json.loads((out/'selection/representatives/manifest.json').read_text()); results=[]; audit=[]
    background_hash=sha256(out/'reference/scoring_background.bfile')
    for m in manifest['candidates']:
        rows=[]
        for ds in m['evaluated_dataset_ids']:
            base=out/'selection/representatives'/m['id']/ds/f"r{m['rotation']:02d}"
            row=next(iter(read_tsv(str(base)+'.metrics.tsv')))
            row={**row,'family_id':m['family_id'],'member_id':m['id'],'matrix_sha256':m['matrix_sha256'],'source_sample':pb[m['dataset_id']],'source_rotation':m['rotation'],'heldout_fold':m['rotation'],'evaluated_role':'pacbio' if ds in pb else 'illumina','input_sha256':sha256(out/'prepared'/ds/f"fold{m['rotation']}.fasta"),'control_sha256':sha256(out/'controls'/ds/f"fold{m['rotation']:02d}.fasta"),'background_sha256':background_hash}
            row['dataset_support']=support(row,s['pacbio_coverage_ratio_min'] if ds in pb else s['illumina_coverage_ratio_min'],s['require_central_enrichment']);rows.append(row)
        decision=validate_candidate(rows,m['required_pacbio_ids'],corr,s)
        rank=representative_rank(m,rows,m['required_pacbio_ids'],corr,s) if decision['pass'] else None
        results.append({'model':m,'decision':decision,'rank':rank,'rows':rows})
        for row in rows:
            group=next(g for p,g in decision['groups'].items() if row['dataset_id']==p or row['dataset_id'] in corr[p])
            audit.append({**row,'group_supporting_count':group['supporting_count'],'group_required_count':group['required_count'],'candidate_pass':decision['pass'],'reason':'validated_exact_matrix' if decision['pass'] else 'required_group_failed'})
    fields=METRIC_FIELDS+['family_id','member_id','matrix_sha256','source_sample','source_rotation','heldout_fold','evaluated_role','input_sha256','control_sha256','background_sha256','dataset_support','group_supporting_count','group_required_count','candidate_pass','reason']
    write_tsv(out/'reports/representative_validation.tsv',audit,fields)
    atomic_json(out/'selection/representatives/validation.json',{'candidates':results})


def rank_job(config):
    out,_,_,_,corr=context(config); stable=json.loads((out/'selection/stability.json').read_text())['families']; validated=json.loads((out/'selection/representatives/validation.json').read_text())['candidates']; coverage=list(read_tsv(out/'reports/family_coverage.tsv')); decisions=[]
    for f in stable:
        candidates=[c for c in validated if c['model']['family_id']==f['family_id'] and c['decision']['pass']]
        chosen=min(candidates,key=lambda c:tuple(float(v) for v in c['rank'][:-1])+(c['rank'][-1],)) if candidates else None
        required=[ds for pb in f['required_pacbio_ids'] for ds in [pb]+corr[pb]]
        dataset_stats=[]; allrat=[]; alldiff=[]; allcov=[]
        for ds in required:
            rows=[r for r in coverage if r['family_id']==f['family_id'] and r['dataset_id']==ds]
            ratios=[number(r['coverage_ratio']) for r in rows]; diffs=[number(r['coverage_difference']) for r in rows]; covs=[number(r['real_coverage']) for r in rows]
            dataset_stats.append({'dataset_id':ds,'median_ratio':median(ratios),'median_difference':median(diffs),'missing_ratio_rows':5-sum(v is not None for v in ratios),'missing_difference_rows':5-sum(v is not None for v in diffs)})
            allrat.extend(v for v in ratios if v is not None);alldiff.extend(v for v in diffs if v is not None);allcov.extend(v for v in covs if v is not None)
        weakest_ratio=min((d['median_ratio'] if d['median_ratio'] is not None else -math.inf for d in dataset_stats),default=-math.inf)
        weakest_diff=min((d['median_difference'] if d['median_difference'] is not None else -math.inf for d in dataset_stats),default=-math.inf)
        supporting=[g for g in f['groups'] if g['support']]
        components=[len(f['stable_supported']),min((f['counts'][p] for p in f['required_pacbio_ids']),default=0),min((g['supporting_count']/g['configured_count'] for g in supporting),default=0),weakest_ratio,median(allrat),weakest_diff,statistics.mean(allcov) if allcov else None]
        decisions.append({'family_id':f['family_id'],'family_stability_pass':f['pass'],'representative_validation_pass':bool(chosen),'included':f['pass'] and bool(chosen),'reason':'selected' if chosen else 'no_validated_representative' if f['pass'] else 'family_stability_failed','representative':chosen['model'] if chosen else None,'required_pacbio_ids':f['required_pacbio_ids'],'components':components,'dataset_statistics':dataset_stats,'stability':f})
    selected=sorted((d for d in decisions if d['included']),key=lambda d:tuple(-(v if v is not None else -math.inf) for v in d['components'])+(d['family_id'],))
    for rank,d in enumerate(selected,1):d['rank']=rank
    atomic_json(out/'selection/ranking.json',{'decisions':decisions,'selected':selected})
    fields=['family_id','family_stability_pass','representative_validation_pass','included','reason','selected_representative','source_rotation','rank','supported_clone_count','minimum_supported_rotations','minimum_corresponding_fraction','weakest_dataset_median_ratio','pooled_median_ratio','weakest_dataset_median_difference','mean_real_coverage','dataset_statistics']
    rows=[]
    for d in decisions:
        m=d['representative']; comp=d['components']
        rows.append(dict(zip(fields,[d['family_id'],d['family_stability_pass'],d['representative_validation_pass'],d['included'],d['reason'],m['id'] if m else '',m['rotation'] if m else '',d.get('rank',''),*comp,json.dumps(d['dataset_statistics'])])))
    write_tsv(out/'reports/selection_decisions.tsv',display_rows(rows),fields)


def selection_reports(config):
    out,a,_,pb,corr=context(config); coverage=[];transfer=[];usage=[];bins=[]
    for p in pb:
        for rotation in range(1,6):
            dest=evdir(out,p,rotation)
            coverage.extend(read_tsv(dest/'pacbio.metrics.tsv'))
            transfer.extend(read_tsv(dest/'transfer.tsv'))
            data=json.loads((dest/'streme.complete.json').read_text())
            usage.append({'dataset_id':p,'rotation':rotation,'status':data['status'],'input_sequences':data.get('input_sequences',0),'input_bases':data.get('input_bases',0),'length_cap':data.get('cap',a['streme_total_length']),'tool_usage_metadata':json.dumps(data.get('usage',{}))})
            shortlist=json.loads((dest/'candidates.json').read_text())['models']
            db=sqlite3.connect(str(dest/'pacbio.metrics.tsv.sqlite'))
            for m in shortlist:
                counts={'real':{},'control':{}};total={}
                for kind in counts:
                    for start,stop in db.execute('SELECT start,stop FROM hits WHERE kind=? AND motif=?',(kind,m['id'])):
                        midpoint=math.floor((start+stop)/2-(a['flank_bp']+1))
                        counts[kind][midpoint]=counts[kind].get(midpoint,0)+1
                    total[kind]=sum(counts[kind].values())
                for midpoint in range(-a['flank_bp'],a['flank_bp']+1):
                    real=counts['real'].get(midpoint,0)/total['real'] if total['real'] else None
                    control=counts['control'].get(midpoint,0)/total['control'] if total['control'] else None
                    bins.append({'dataset_id':p,'rotation':rotation,'motif_id':m['id'],'midpoint_bin':midpoint,'real_conditional_proportion':real,'control_conditional_proportion':control,'conditional_difference':real-control if real is not None and control is not None else None,'real_hit_windows':total['real'],'control_hit_windows':total['control']})
            db.close()
    write_tsv(out/'reports/motif_coverage.tsv',coverage,METRIC_FIELDS)
    write_tsv(out/'reports/transfer_evidence.tsv',transfer,METRIC_FIELDS+['support','supporting_count','required_count','configured_count','pacbio_dataset_id'])
    write_tsv(out/'reports/streme_usage.tsv',usage,['dataset_id','rotation','status','input_sequences','input_bases','length_cap','tool_usage_metadata'])
    write_tsv(out/'reports/conditional_position_bins.tsv',display_rows(bins),['dataset_id','rotation','motif_id','midpoint_bin','real_conditional_proportion','control_conditional_proportion','conditional_difference','real_hit_windows','control_hit_windows'])


def freeze_job(config):
    out,a,_,pb,corr=context(config); library=out/'library';library.mkdir(parents=True,exist_ok=True)
    data=json.loads((out/'selection/ranking.json').read_text());models=[];rows=[];members=[]
    selection_reports(config)
    for d in data['selected']:
        original=d['representative']; final={**original,'id':'final_'+d['family_id']};models.append(final)
        rows.append({'final_motif_id':final['id'],'family_id':d['family_id'],'global_rank':d['rank'],'rank':d['rank'],'source_member_id':original['id'],'source_dataset_id':original['dataset_id'],'source_rotation':original['rotation'],'original_id':original['original_id'],'matrix_sha256':original['matrix_sha256'],'width':original['width'],'nsites':original['nsites'],'evalue':original['evalue'],'required_pacbio_ids':json.dumps(d['required_pacbio_ids']),'family_stability_evidence_sha256':sha256(out/'selection/stability.json'),'representative_validation_sha256':sha256(out/'reports/representative_validation.tsv')})
        for mid in d['stability']['members']:members.append({'final_motif_id':final['id'],'family_id':d['family_id'],'member_id':mid})
    shutil.copyfile(out/'reference/scoring_background.bfile',library/'scoring_background.bfile')
    write_meme(library/'final_motifs.meme',models,library/'scoring_background.bfile')
    fields=['final_motif_id','family_id','global_rank','rank','source_member_id','source_dataset_id','source_rotation','original_id','matrix_sha256','width','nsites','evalue','required_pacbio_ids','family_stability_evidence_sha256','representative_validation_sha256']
    write_tsv(library/'final_motif_manifest.tsv',rows,fields);write_tsv(library/'final_motif_members.tsv',members,['final_motif_id','family_id','member_id'])
    tmp=library/'selection_loci.sqlite.tmp'
    if tmp.exists():tmp.unlink()
    db=sqlite3.connect(tmp);db.execute('CREATE TABLE loci (CHROM TEXT, POS INTEGER, PRIMARY KEY(CHROM,POS))')
    for ds in list(pb)+[i for values in corr.values() for i in values]:
        for row in read_tsv(out/'prepared'/ds/'genotypes.tsv.gz'):
            if row.get('gt_status')=='alt' or any(t.isdigit() and int(t)>0 for t in re.split(r'[/|]',row['GT'])):
                db.execute('INSERT OR IGNORE INTO loci VALUES (?,?)',(row['CHROM'],int(row['POS'])))
    db.commit();db.close();os.replace(tmp,library/'selection_loci.sqlite')
    outputs=[library/n for n in ['final_motifs.meme','final_motif_manifest.tsv','final_motif_members.tsv','scoring_background.bfile','selection_loci.sqlite']]
    evidence=[out/'reports/representative_validation.tsv',out/'selection/representatives/validation.json',out/'selection/stability.json',out/'reports/selection_decisions.tsv',out/'manifests/folds.json',out/'manifests/seeds.tsv',out/'reference/reference.fa',out/'reference/reference.fa.fai']
    evidence.extend(out/'reports'/name for name in ['motif_coverage.tsv','transfer_evidence.tsv','streme_usage.tsv','conditional_position_bins.tsv'])
    evidence.extend(out/'manifests'/name for name in ['selection_config.json','selection_fingerprints.json','datasets.tsv','preflight.ok.json'])
    freeze={'schema_version':1,'motif_count':len(models),'outputs':{str(p.resolve()):sha256(p) for p in outputs+evidence},'reference_fasta':str((out/'reference/reference.fa').resolve()),'reference_sha256':sha256(out/'reference/reference.fa'),'flank_bp':a['flank_bp'],'fimo_p_threshold':a['fimo_p_threshold'],'selection_loci':str((library/'selection_loci.sqlite').resolve()),'background_sha256':sha256(library/'scoring_background.bfile'),'selection_config_sha256':sha256(out/'manifests/selection_config.json'),'family_stability':{'path':str(out/'selection/stability.json'),'sha256':sha256(out/'selection/stability.json'),'rotations':[1,2,3,4,5]},'representative_validation':{'path':str(out/'reports/representative_validation.tsv'),'sha256':sha256(out/'reports/representative_validation.tsv'),'scope':'exact matrix on its source-rotation held-out fold across predetermined required groups'},'selected_samples':list(pb.values())+[n for vals in config['illumina_corresponding'].values() for n in vals],'sample_metadata':{k:v for k,v in config.get('sample_metadata',{}).items() if k!='other_illumina'},'tool_containers':config.get('container',{}),'code_sha256':sha256(__file__)}
    freeze['required_outputs']=freeze['outputs']
    freeze['representative_validation_sha256']=freeze['representative_validation']['sha256']
    freeze['selection_config']=json.loads((out/'manifests/selection_config.json').read_text())
    freeze['producer_code_sha256']={str(Path(__file__).resolve()):sha256(__file__),str(Path(__file__).with_name('prepare.py').resolve()):sha256(Path(__file__).with_name('prepare.py'))}
    freeze['fimo_background_flag']=json.loads((out/'manifests/preflight/meme.ok.json').read_text())['fimo_background_flag']
    atomic_json(library/'freeze.json',freeze)


def preflight_meme(workdir):
    workdir=Path(workdir);workdir.mkdir(parents=True,exist_ok=True)
    workdir=Path(tempfile.mkdtemp(prefix='attempt-',dir=workdir))
    versions={}
    for tool in ['streme','fimo','tomtom']:
        result=subprocess.run([tool,'--version'],capture_output=True,text=True)
        text=result.stdout+result.stderr
        if result.returncode or '5.5.7' not in text:raise RuntimeError(f'Expected {tool} 5.5.7: {text}')
        versions[tool]=text.strip()
    fasta=workdir/'tiny.fasta';write_fasta(fasta,[(f'w{i}',('ACGTTGCACTGA' if i%2 else 'GCTAGACCTTGC')*10) for i in range(100)])
    background=workdir/'background.bfile';run_tool(['fasta-get-markov','-m','0','-pseudo','0.1',fasta,background],workdir/'background.log');background.write_text(background.read_text(encoding='utf-8').replace(str(fasta),'preflight.fasta'),encoding='utf-8')
    shuffled=workdir/'shuffle.fasta';shuffle_fasta(fasta,shuffled,3,42,workdir/'shuffle.log')
    output=workdir/'streme';run_tool(['streme','--dna','--p',fasta,'--n',shuffled,'--minw','6','--maxw','12','--seed','42','--totallength','5000','--o',output],workdir/'streme.log')
    parse_meme(output/'streme.txt')
    matrix={'id':'known','width':12,'nsites':'20','evalue':0.001,'matrix':[['0.97' if b==x else '0.01' for b in 'ACGT'] for x in 'ACGTTGCACTGA']}
    models=workdir/'known.meme';write_meme(models,[matrix],background)
    flag=None
    helptext=subprocess.run(['fimo'],capture_output=True,text=True)
    for candidate in ['--bgfile','--bfile']:
        if candidate in helptext.stdout+helptext.stderr:flag=candidate;break
    if flag is None:raise RuntimeError('FIMO lacks an explicit background flag')
    raw=workdir/'fimo.tsv';scan_fimo(models,fasta,background,raw,0.0001,workdir/'fimo.log',flag)
    occurrences=list(iter_fimo(raw,{'known'},{f'w{i}' for i in range(100)}))
    if not occurrences or any(r['stop']-r['start']+1!=12 for r in occurrences):raise RuntimeError('FIMO coordinate/background smoke failed')
    tom=workdir/'tomtom';run_tool(['tomtom','-oc',tom,'-dist','ed','-min-overlap','6','-thresh','0.05',models,models],workdir/'tomtom.log')
    list(read_tsv(tom/'tomtom.tsv'))
    import sys
    return {'versions':versions,'bundled_python_version':sys.version,'fimo_background_flag':flag,'background_sha256':sha256(background),'commands':'real shuffle/background/capped STREME/text FIMO/Tomtom; tool-only tiny fixtures'}


def main():
    p=argparse.ArgumentParser();p.add_argument('action');p.add_argument('--config',required=True)
    for key in ['dataset','source','member','kind','models','fasta','background','output','log']:p.add_argument('--'+key)
    p.add_argument('--primary');p.add_argument('--inputs',nargs='*',default=[]);p.add_argument('--outputs',nargs='*',default=[])
    p.add_argument('--rotation',type=int);p.add_argument('--disk-mb',type=int);args=p.parse_args();config=json.loads(Path(args.config).read_text());out,a,s,pb,corr=context(config)
    action=args.action;ds=args.dataset;r=args.rotation
    disk_mb = args.disk_mb or {'streme':20000, 'scan':30000, 'summarize':60000, 'stability':60000, 'validate-representatives':60000}.get(action,10000)
    configure_scratch(config, action, args.member or args.source or ds or 'global', r or 0, disk_mb)
    if action=='split':make_train_test(config,ds,r)
    elif action=='shuffle':
        if args.kind=='train':source=evdir(out,ds,r)/'train.fasta';output=evdir(out,ds,r)/'train.control.fasta';role='train_shuffle'
        else:source=out/'prepared'/ds/f'fold{r}.fasta';output=out/'controls'/ds/f'fold{r:02d}.fasta';role='test_shuffle'
        shuffle_fasta(source,output,a['shuffle_kmer'],seed(config,ds,r,role),str(output)+'.stderr.log')
    elif action=='streme':discover(config,ds,r)
    elif action=='parse-models':parse_models_job(config,ds,r)
    elif action=='background':make_background(out/'reference/reference.fa',out/'reference/scoring_background.bfile',out/'logs/background.stderr.log')
    elif action=='scan':
        flag=json.loads((out/'manifests/preflight/meme.ok.json').read_text())['fimo_background_flag']
        scan_fimo(args.models,args.fasta,args.background,args.output,a['fimo_p_threshold'],args.log,flag)
    elif action=='summarize':summarize_job(config,ds,r,args.source,args.member)
    elif action=='candidates':candidates_job(config,ds,r)
    elif action=='transfer':transfer_job(config,ds,r)
    elif action=='pool':pool_job(config)
    elif action=='tomtom':tomtom_job(config)
    elif action=='families':families_job(config)
    elif action=='stability':stability_job(config)
    elif action=='representatives':representatives_job(config)
    elif action=='validate-representatives':validate_representatives_job(config)
    elif action=='rank':rank_job(config)
    elif action=='freeze':freeze_job(config)
    else:raise ValueError(f'Unknown action: {action}')
    if args.primary:
        metric_keys=['analysis.flank_bp','analysis.center_radius_bp']
        paths={
            'split':[], 'background':['container.meme'],
            'shuffle':['analysis.shuffle_kmer','analysis.base_seed','analysis.seed_overrides','container.meme'],
            'streme':['analysis.streme_min_width','analysis.streme_max_width','analysis.streme_total_length','analysis.base_seed','analysis.seed_overrides','container.meme'],
            'parse-models':[],
            'scan':['analysis.fimo_p_threshold','container.meme'],
            'summarize':metric_keys,
            'candidates':['analysis.streme_e_max','selection.pacbio_coverage_ratio_min','selection.require_central_enrichment'],
            'transfer':['illumina_corresponding','selection.illumina_coverage_ratio_min','selection.corresponding_dataset_support_fraction','selection.require_central_enrichment'],
            'pool':[], 'tomtom':['selection.family_matching','container.meme'],
            'families':['selection.family_matching'],
            'stability':metric_keys+['selection','pacbio_samples','illumina_corresponding'],
            'representatives':['selection.require_all_pacbio_clones'],
            'validate-representatives':['selection','illumina_corresponding'],
            'rank':['selection','pacbio_samples','illumina_corresponding'],
            'freeze':['analysis.flank_bp','analysis.fimo_p_threshold','selection','container','sample_metadata.pacbio','sample_metadata.illumina']
        }[action]
        def lookup(key):
            value=config
            for component in key.split('.'):value=value.get(component,{}) if isinstance(value,dict) else {}
            return value
        if action=='representatives':
            args.outputs.extend(str(out/'selection/representative_models'/m['id']/'model.meme') for m in json.loads((out/'selection/representatives/manifest.json').read_text())['candidates'])
        publish_provenance(args.primary,args.inputs,{'action':action,'dataset':ds,'source':args.source,'rotation':r,'member':args.member,'kind':args.kind,'config_paths':paths,'config_values':{key:lookup(key) for key in paths}},args.outputs)


if __name__=='__main__':main()
