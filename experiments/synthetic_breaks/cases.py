"""Fixed synthetic-study cases and nested random-break ensembles."""
import numpy as np
from segmentation_sensitivity.geometry import Break


def cases_for(name, seeds=(11,29,47,83,101)):
    cases = [('intact',[],None)]
    if name == 'line':
        for fraction in [.2,.5,.8]:
            cases.append((f'gap_at_{fraction:.1f}',[Break(0,fraction,.02)],None))
        for width in [.01,.04]:
            cases.append((f'gap_width_{width:.2f}',[Break(0,.5,width)],None))
    elif name == 'y':
        for bid,label in enumerate(['trunk','upper','lower']):
            cases.append((f'gap_{label}',[Break(bid)],None))
        cases.append(('gap_both_daughters',[Break(1),Break(2)],None))
    elif name == 'tree':
        for seed in seeds:
            order = np.random.default_rng(seed).permutation(63)
            for count in [1,2,4,8,16,32]:
                cases.append((f'seed_{seed}_breaks_{count}',
                              [Break(int(b)) for b in order[:count]],seed))
    else:
        raise ValueError(name)
    return cases
