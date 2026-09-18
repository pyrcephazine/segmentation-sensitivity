#!/usr/bin/env python3
"""Submit one experiment job or the complete workflow to Slurm."""
import argparse
from pathlib import Path
import re
import shlex
import subprocess
import sys

from experiments.jobs import JOBS, PIPELINE

ROOT = Path(__file__).resolve().parent


def batch_script(name, root=ROOT):
    """Render a fail-fast script; preserve all scientific command arguments."""
    lines = ['#!/usr/bin/env bash', 'set -euo pipefail']
    if name == 'check':
        lines.append('export STUDY_WORK=' + shlex.quote(str(root / '.work/refactor')))
    lines.append('source ' + shlex.quote(str(root / 'environment.sh')))
    lines.append(JOBS[name]['commands'])
    return '\n'.join(lines) + '\n'


def submit_job(name, args, dependencies=(), root=ROOT):
    job = JOBS[name]
    logs = root / '.work/logs'
    command = ['sbatch', '--parsable', '--chdir=' + str(root),
               '--job-name=seg-' + name, '--account=' + args.account,
               '--partition=' + args.partition, '--qos=' + args.qos,
               '--nodes=1', '--ntasks=' + str(job['tasks']), '--cpus-per-task=1',
               '--mem=' + job['memory'], '--time=' + job['time'],
               '--output=' + str(logs / (name + '-%j.log'))]
    if dependencies:
        command.append('--dependency=afterok:' + ':'.join(dependencies))
    script = batch_script(name, root)
    if args.dry_run:
        print('# ' + shlex.join(command))
        print(script)
        return '<' + name + '-job-id>'
    logs.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(command, input=script, text=True, capture_output=True, check=True)
    job_id = result.stdout.strip().split(';', 1)[0]
    if not job_id.isdigit():
        raise RuntimeError('Unexpected sbatch output; submission may have succeeded: ' + result.stdout.strip())
    print(name + ': submitted job ' + job_id, flush=True)
    return job_id


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('job', choices=['all', *JOBS])
    parser.add_argument('--dry-run', action='store_true', help='Print commands and scripts without submitting or creating files')
    parser.add_argument('--after', nargs='+', default=[], metavar='JOB_ID', help='Wait for these jobs to succeed before starting')
    parser.add_argument('--account', default='commons')
    parser.add_argument('--partition', default='commons')
    parser.add_argument('--qos', default='nots_commons')
    args = parser.parse_args(argv)
    if any(not re.fullmatch(r'[0-9]+', value) for value in args.after):
        parser.error('--after requires numeric Slurm job IDs')
    if args.job == 'all':
        jobs = [submit_job(name, args, args.after) for name in PIPELINE]
        submit_job('reports', args, jobs)
    else:
        submit_job(args.job, args, args.after)


if __name__ == '__main__':
    try:
        main()
    except subprocess.CalledProcessError as error:
        print('Submission failed: ' + (error.stderr or str(error)).strip(), file=sys.stderr)
        print('Earlier job IDs, if any, are printed above; those jobs remain submitted.', file=sys.stderr)
        sys.exit(1)
    except (OSError, RuntimeError) as error:
        print(str(error), file=sys.stderr)
        sys.exit(1)
