"""Record cluster access and a small, temporary sequential direct-I/O sample."""
from experiments.synthetic_breaks.config import WORK
from pathlib import Path
import datetime
import json
import mmap
import os
import socket
import subprocess
import tempfile
import time


def command(argv):
    try:
        result=subprocess.run(argv,text=True,capture_output=True,timeout=30)
        return dict(command=argv,returncode=result.returncode,
                    stdout=result.stdout,stderr=result.stderr)
    except (OSError,subprocess.TimeoutExpired) as error:
        return dict(command=argv,error=str(error))


def main():
    root=WORK
    paths=['/home/pzz1','/projects/br1/pzz1','/scratch/pzz1','/tmp',
           '/projects/comp416','/projects/comp646','/storage/hpc/work',
           '/rhf/allocations/ea14']
    result=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                host=socket.gethostname(),slurm_job=os.environ.get('SLURM_JOB_ID'),
                paths=[dict(path=p,exists=os.path.exists(p),readable=os.access(p,os.R_OK),
                            writable=os.access(p,os.W_OK),searchable=os.access(p,os.X_OK)) for p in paths],
                commands=[command(c) for c in [
                    ['id'],['df','-hT',*paths[:4]],['quota','-s'],['quota','-gs'],
                    ['sacctmgr','show','assoc','user=pzz1','format=User,Account,Partition,QOS','-n','-P'],
                    ['git','--version'],
                    ['git','ls-remote','https://github.com/pyrcephazine/multihepatic','HEAD'],
                    ['git','ls-remote','https://github.com/scientificcomputing/networks_fenicsx','HEAD']]],
                benchmark_method='One 256 MiB incompressible sequential file per filesystem; 4 MiB operations; O_DIRECT write + fsync, then O_DIRECT read. Temporary files removed. Client cache bypass does not bypass storage-server caches. A single-node sample, not a parallel filesystem benchmark.',
                benchmark=[])
    size=256*2**20; block=4*2**20
    with mmap.mmap(-1,size) as payload:
        payload[:]=os.urandom(size)
        for parent in paths[:4]:
            fd,path=tempfile.mkstemp(prefix='.segmentation-io-check-',dir=parent)
            os.close(fd)
            record=dict(path=parent,bytes=size)
            try:
                fd=os.open(path,os.O_WRONLY|os.O_TRUNC|os.O_DIRECT)
                try:
                    start=time.perf_counter()
                    for offset in range(0,size,block):
                        with memoryview(payload)[offset:offset+block] as view:
                            if os.write(fd,view)!=block:raise RuntimeError('Incomplete benchmark write')
                    os.fsync(fd)
                    record['write_MiB_per_second']=size/2**20/(time.perf_counter()-start)
                finally:os.close(fd)
                fd=os.open(path,os.O_RDONLY|os.O_DIRECT)
                try:
                    start=time.perf_counter()
                    for offset in range(0,size,block):
                        with memoryview(payload)[offset:offset+block] as view:
                            if os.readv(fd,[view])!=block:raise RuntimeError('Incomplete benchmark read')
                    record['read_MiB_per_second']=size/2**20/(time.perf_counter()-start)
                finally:os.close(fd)
            except OSError as error:
                record['error']=str(error)
            finally:Path(path).unlink(missing_ok=True)
            result['benchmark'].append(record)
    (root/'cluster-access.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['benchmark'],indent=2))


if __name__=='__main__':main()
