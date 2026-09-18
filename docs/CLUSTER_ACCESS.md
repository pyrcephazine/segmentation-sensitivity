# NOTS access and storage check

Keep the approximately 2 GB input ZIP in `/projects/br1/pzz1/segmentation-sensitivity/`.
That is where `oilblooddata.zip` was found, and it is suitable for retention beyond
two weeks. Use `/scratch/pzz1/segmentation-sensitivity/` for active simulation I/O,
and copy selected completed results to home or project storage.

## Storage available to this account

* `/home/pzz1`: writable, NFS, 10 GiB user quota; about 5.22 GiB used before packaging this study.
* `/projects/br1/pzz1`: writable, NFS; the br1 group used 36,033 MiB of its shared 100,000 MiB quota at the check.
* `/scratch/pzz1`: writable, VAST over NFSv4, shared across nodes.
* `/tmp`: writable local XFS on the current compute node; not shared across nodes.
* `/projects/comp416`: the group directory is writable by this account; no study files were placed there.
* `/projects/comp646`: the group directory can be read/searched but not written by this account.
* The configured `/storage/hpc/work` path is absent on this node. `/rhf/allocations/ea14` appears in quota output but is inaccessible to this account.

Rice documents no purge interval for home/projects, a 14-day purge interval for
shared scratch, and job-lifetime retention for local temporary storage. These
retention categories make projects the appropriate place for the persistent ZIP.
See [Rice's storage table](https://kb.rice.edu/147975) and
[quota and I/O guidance](https://kb.rice.edu/147979).

## Measured sequential I/O sample

Measured on `bb6u11g1` at 2026-09-11T19:36:24.429868+00:00 in Slurm allocation 1222926.

| Location | Write (MiB/s) | Read (MiB/s) |
| --- | --- | --- |
| `/home/pzz1` | 1248 | 1887 |
| `/projects/br1/pzz1` | 632 | 870 |
| `/scratch/pzz1` | 289 | 905 |
| `/tmp` | 1685 | 4424 |

Each sample wrote and read one 256 MiB incompressible temporary file with 4 MiB
operations. Writes used O_DIRECT and fsync; reads used O_DIRECT. The files were
removed. This bypasses the client page cache but does not bypass server caches.
These are short single-process samples, not sustained multi-node throughput or
metadata benchmarks. In this sample local temporary storage was fastest; home
also performed well. Scratch is the cluster's designated shared location for job
I/O even though this small write sample was slower. Filesystem-wide free space
from `df` is not a user's quota. Raw outputs are in `execution/cluster-access.json`.

## Compute and network access

The session already ran on NOTS compute node bb6u11g1 as pzz1, with groups br1,
comp646, comp416, hpcusers, and nots. Slurm associations permit commons, debug,
long, and scavenge partitions through the commons account. The study and
validation jobs ran in commons. MPI solved the same coupled problem with one,
two, and eight ranks; numerical comparisons are in the study report.

An SSH attempt to `pzz1@nots.rice.edu` from the existing compute session failed
public-key authentication. No usable private key or SSH agent was present in
that session, although an authorized_keys file existed. This did not prevent
running here or submitting Slurm jobs. Nothing was changed in SSH authorization.

Web search successfully retrieved the official Rice storage pages above. Git
2.50.1 successfully queried both upstream repositories over HTTPS. Git is a
module-provided command in this environment; load its compiler dependency first:

```bash
module load GCC/14.3.0 git/2.50.1
git ls-remote https://github.com/pyrcephazine/multihepatic HEAD
git ls-remote https://github.com/scientificcomputing/networks_fenicsx HEAD
```

The returned revisions matched the checked-out starting sources:

* multihepatic: `2283b8eda631fe3d62931abeed35c0c8e2817f74`
* networks_fenicsx: `4d0dd39d6a788185c589c9f4c5017701cd8d27db`

The FEniCSx/PETSc/OpenMPI stack is installed through native Rice modules and a
project-local Python virtual environment. Run `source scripts/environment.sh`
from the study source checkout to activate it. The delivered results are in the
NOTS home directory; the transfer command in START_HERE.md copies the compressed
bundle to a separate workstation.
