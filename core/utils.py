import numpy as np

def estimate_fdr_and_power_with_arrays(rejections, hypos):
    """
    This function estimate the fdr: Proportion of reject|H0 True,
    and power: Proportion of H1's that are rejected.

    Input: 

    A boolean array of rejections.
    A boolean array of underlying truth.

    Output:

    The proportion of H0's that are rejected.
    The proportion of H1's that are rejected.
    """  

    # rejected H0.
    false_discoveries = np.logical_and(np.logical_not(hypos),rejections)

    # rejected H1.
    true_discoveries = np.logical_and(hypos, rejections)

    num_H0 = float(sum(np.logical_not(hypos)))
    num_H1 = float(sum(hypos)) 
    fdr = sum(false_discoveries) / float(sum(rejections)) if float(sum(rejections)) != 0 else 0

    power = sum(true_discoveries) / num_H1 if num_H1 != 0 else 0

    # return fdr and power.
    return  fdr, power