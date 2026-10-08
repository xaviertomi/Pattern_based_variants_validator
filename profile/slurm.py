#!/usr/bin/env python3
"""cluster-generic SLURM submission, accounting status and cancellation."""
import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STATUS = ROOT / 'logs' / 'status'
RUNNING = set('PENDING RUNNING CONFIGURING COMPLETING SUSPENDED REQUEUED REQUEUE_FED REQUEUE_HOLD RESIZING STAGE_OUT SIGNALING'.split())
FAILED = set('BOOT_FAIL CANCELLED DEADLINE FAILED NODE_FAIL OUT_OF_MEMORY PREEMPTED REVOKED TIMEOUT'.split())

def command(args):
    p = subprocess.run(args, text=True, capture_output=True)
    if p.stderr:
        print(p.stderr.rstrip(), file=sys.stderr)
    if p.returncode:
        raise RuntimeError(f'{args[0]} exited {p.returncode}')
    return p.stdout.strip()

def job_id(value):
    if not re.fullmatch(r'[0-9]+', value):
        raise ValueError(f'Invalid numeric SLURM job ID: {value}')
    return value

def record_path(jid):
    return STATUS / (job_id(jid) + '.json')

def load(jid):
    p = record_path(jid)
    return json.loads(p.read_text()) if p.exists() else {}

def save(jid, data):
    STATUS.mkdir(parents=True, exist_ok=True)
    p = record_path(jid)
    tmp = p.with_suffix('.tmp')
    tmp.write_text(json.dumps(data, sort_keys=True) + '\n')
    tmp.replace(p)

def cluster_args(data):
    return ['--clusters', data['cluster']] if data.get('cluster') else []

def submit(a):
    for p in (a.output, a.error):
        Path(p).parent.mkdir(parents=True, exist_ok=True)
    result = command(['sbatch', '--parsable', '--export=ALL', '--partition', a.partition,
        '--account', a.account, '--cpus-per-task', str(a.cpus_per_task), '--mem', str(a.mem),
        '--time', a.time, '--job-name', a.job_name, '--output', a.output, '--error', a.error, a.jobscript])
    if not re.fullmatch(r'[0-9]+(?:;[^;\s]+)?', result):
        raise RuntimeError(f'Invalid sbatch parsable response: {result!r}')
    parts = result.split(';')
    jid = parts[0]
    save(jid, {'cluster': parts[1] if len(parts) == 2 else None})
    print(jid)

def status(jid):
    jid = job_id(jid)
    data = load(jid)
    text = command(['sacct', *cluster_args(data), '-X', '-n', '-P', '-j', jid,
                    '-o', 'JobIDRaw,State%40,ExitCode'])
    rows = [line.split('|') for line in text.splitlines() if line.split('|')[0] == jid]
    if rows:
        if len(rows) != 1 or len(rows[0]) < 3:
            raise RuntimeError(f'Ambiguous accounting allocation row for {jid}')
        state = rows[0][1].strip().split()[0].rstrip('+')
        exitcode = rows[0][2].strip()
        if state in RUNNING:
            result = 'running'
        elif state in FAILED:
            result = 'failed'
        elif state == 'COMPLETED':
            result = 'success' if exitcode == '0:0' else 'failed'
        else:
            raise RuntimeError(f'Unknown SLURM state for {jid}: {state}')
        if 'first_absence' in data:
            data.pop('first_absence')
            save(jid, data)
        print(result)
        return
    queue = command(['squeue', *cluster_args(data), '-h', '-j', jid, '-o', '%i|%T'])
    if any(line.split('|')[0] == jid for line in queue.splitlines()):
        print('running')
        return
    now = time.time()
    if 'first_absence' not in data:
        data['first_absence'] = now
        save(jid, data)
    if now - data['first_absence'] < 120:
        print('running')
    else:
        print(f'Job {jid} absent from both sacct and squeue for 120 seconds', file=sys.stderr)
        print('failed')

def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='action', required=True)
    s = sub.add_parser('submit')
    for key in ('partition', 'account', 'cpus-per-task', 'mem', 'time', 'job-name', 'output', 'error'):
        s.add_argument('--' + key, required=True)
    s.add_argument('jobscript')
    sub.add_parser('status').add_argument('jobid')
    sub.add_parser('cancel').add_argument('jobids', nargs='+')
    a = p.parse_args()
    if a.action == 'submit':
        submit(a)
    elif a.action == 'status':
        status(a.jobid)
    else:
        for jid in a.jobids:
            command(['scancel', *cluster_args(load(jid)), job_id(jid)])

if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)
