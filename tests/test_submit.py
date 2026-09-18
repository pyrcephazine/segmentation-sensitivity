"""Submission safety and dependency wiring, without contacting Slurm."""
import argparse
import contextlib
import io
from pathlib import Path
import subprocess
import unittest
from unittest.mock import patch

import submit


class SubmissionTests(unittest.TestCase):
    def args(self, **changes):
        options=dict(account='commons',partition='commons',qos='nots_commons',dry_run=False)
        options.update(changes)
        return argparse.Namespace(**options)

    def test_full_workflow_reports_wait_for_every_job(self):
        responses=[subprocess.CompletedProcess([],0,str(i)+';cluster\n','') for i in range(101,107)]
        with patch('submit.subprocess.run',side_effect=responses) as run, patch('submit.Path.mkdir'), contextlib.redirect_stdout(io.StringIO()):
            submit.main(['all','--after','90'])
        self.assertEqual(run.call_count,6)
        for call in run.call_args_list[:5]:
            self.assertIn('--dependency=afterok:90',call.args[0])
        self.assertIn('--dependency=afterok:101:102:103:104:105',run.call_args_list[-1].args[0])
        self.assertNotIn('package.py',run.call_args_list[-1].kwargs['input'])

    def test_dry_run_neither_submits_nor_creates_files(self):
        with patch('submit.subprocess.run') as run, patch('submit.Path.mkdir') as mkdir, contextlib.redirect_stdout(io.StringIO()) as output:
            submit.main(['all','--dry-run'])
        run.assert_not_called();mkdir.assert_not_called()
        self.assertIn('--dependency=afterok:',output.getvalue())
        self.assertIn('experiments.pt002_inflow.run --prepare',output.getvalue())

    def test_submission_failure_stops_pipeline(self):
        success=subprocess.CompletedProcess([],0,'101\n','')
        error=subprocess.CalledProcessError(1,['sbatch'],stderr='quota exceeded')
        with patch('submit.subprocess.run',side_effect=[success,error]) as run, patch('submit.Path.mkdir'), contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(subprocess.CalledProcessError):submit.main(['all'])
        self.assertEqual(run.call_count,2)

    def test_unsafe_dependency_rejected_before_submission(self):
        with patch('submit.subprocess.run') as run, contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):submit.main(['reports','--after','1;touch /tmp/x'])
        run.assert_not_called()

    def test_paths_with_spaces_and_resource_profile(self):
        root=Path('/tmp/a study')
        response=subprocess.CompletedProcess([],0,'123\n','')
        with patch('submit.subprocess.run',return_value=response) as run, patch('submit.Path.mkdir'), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(submit.submit_job('supplied',self.args(),root=root),'123')
        call=run.call_args
        self.assertIn('--chdir=/tmp/a study',call.args[0])
        self.assertIn('--ntasks=4',call.args[0])
        self.assertIn("source '/tmp/a study/environment.sh'",call.kwargs['input'])
        self.assertIn('set -euo pipefail',call.kwargs['input'])

    def test_each_profile_generates_valid_bash(self):
        for name in submit.JOBS:
            with self.subTest(name=name):
                subprocess.run(['bash','-n'],input=submit.batch_script(name),text=True,check=True)


if __name__=='__main__':unittest.main()
