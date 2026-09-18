"""Resource profiles and experiment commands used by the root submitter."""

JOBS = {'synthetic': {'tasks': 8,
               'memory': '32G',
               'time': '02:00:00',
               'commands': 'mpirun -n 8 python -u -m experiments.synthetic_breaks.run --geometry '
                           'toys --n 8 --cases gap_at_0.5,gap_trunk,gap_upper --output '
                           '${STUDY_WORK}/synthetic_breaks/validation/mpi8\n'
                           'mpirun -n 8 python -u -m experiments.synthetic_breaks.run --geometry '
                           'toys --n 24\n'
                           'mpirun -n 8 python -u -m experiments.synthetic_breaks.run --geometry '
                           'tree --n 40'},
 'validate': {'tasks': 4,
              'memory': '24G',
              'time': '01:00:00',
              'commands': 'mpirun -n 4 python -u -m experiments.synthetic_breaks.validate\n'
                          'mpirun -n 4 python -u -m experiments.synthetic_breaks.run --geometry '
                          'toys --n 12 --cases gap_at_0.5,gap_trunk,gap_upper --output '
                          '${STUDY_WORK}/synthetic_breaks/validation/mesh\n'
                          'mpirun -n 4 python -u -m experiments.synthetic_breaks.run --geometry '
                          'toys --n 36 --cases gap_at_0.5,gap_trunk,gap_upper --output '
                          '${STUDY_WORK}/synthetic_breaks/validation/mesh\n'
                          'mpirun -n 4 python -u -m experiments.synthetic_breaks.run --geometry '
                          'toys --n 24 --h .005 --cases gap_at_0.5,gap_trunk,gap_upper --output '
                          '${STUDY_WORK}/synthetic_breaks/validation/network_refinement\n'
                          'mpirun -n 4 python -u -m experiments.synthetic_breaks.run --geometry '
                          'toys --n 24 --circle-degree 40 --cases gap_at_0.5,gap_trunk,gap_upper '
                          '--output "${STUDY_WORK}/synthetic_breaks/validation/quadrature"'},
 'mpi-check': {'tasks': 2,
               'memory': '16G',
               'time': '00:30:00',
               'commands': 'mpirun -n 1 python -u -m experiments.synthetic_breaks.run --geometry '
                           'toys --n 8 --cases gap_at_0.5,gap_trunk,gap_upper --output '
                           '${STUDY_WORK}/synthetic_breaks/validation/mpi1\n'
                           'mpirun -n 2 python -u -m experiments.synthetic_breaks.run --geometry '
                           'toys --n 8 --output "${STUDY_WORK}/synthetic_breaks/validation/mpi2"\n'
                           'mpirun -n 1 python -u -m experiments.synthetic_breaks.run --geometry '
                           'toys --n 8 --inlet 2 --cases gap_at_0.5,gap_trunk,gap_upper --output '
                           '${STUDY_WORK}/synthetic_breaks/validation/amplitude2'},
 'mesh-check': {'tasks': 4,
                'memory': '32G',
                'time': '01:00:00',
                'commands': 'mpirun -n 4 python -u -m experiments.synthetic_breaks.run --geometry '
                            'tree --n 32 --cases seed_11_breaks_8,seed_11_breaks_32 --output '
                            '${STUDY_WORK}/synthetic_breaks/validation/large_mesh\n'
                            'mpirun -n 4 python -u -m experiments.synthetic_breaks.run --geometry '
                            'tree --n 48 --cases seed_11_breaks_8,seed_11_breaks_32 --output '
                            '${STUDY_WORK}/synthetic_breaks/validation/large_mesh'},
 'supplied': {'tasks': 4,
              'memory': '16G',
              'time': '00:45:00',
              'commands': 'python -m experiments.pt002_inflow.run --prepare\n'
                          'mpirun -n 1 python -u -m experiments.synthetic_breaks.run --geometry '
                          'line --n 8 --cases intact --output '
                          '${STUDY_WORK}/synthetic_breaks/validation/import_regression\n'
                          'mpirun -n 4 python -u -m experiments.pt002_inflow.run --h .005 --name '
                          'intact\n'
                          'mpirun -n 4 python -u -m experiments.pt002_inflow.run --h .0025 --name '
                          'network_refinement'},
 'reports': {'tasks': 1,
             'memory': '8G',
             'time': '00:20:00',
             'commands': 'python -m experiments.synthetic_breaks.analyze\n'
                         'python -m experiments.synthetic_breaks.report\n'
                         'python -m experiments.pt002_inflow.report\n'
                         'python -m experiments.collect'},
 'check': {'tasks': 2,
           'memory': '8G',
           'time': '00:15:00',
           'commands': 'mpirun -n 2 python -u -m experiments.synthetic_breaks.regression '
                       '--geometry toys --n 8 --cases gap_at_0.5,gap_trunk,gap_upper'}}

PIPELINE = ("synthetic", "validate", "mpi-check", "mesh-check", "supplied")
