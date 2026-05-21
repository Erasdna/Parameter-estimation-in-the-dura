from typing import Callable
from scipy.stats import qmc
import numpy as np
from multiprocessing import Pool
import pandas as pd


def sample_parameters(parameter_bounds: dict, n_qmc_samples: int, sampler):
    sample_unit = sampler.random(n_qmc_samples)

    tasks = []
    for i in range(n_qmc_samples):
        p_vals = []
        for j, parameter_bound in enumerate(parameter_bounds.values()):
            lower = float(parameter_bound[0])
            upper = float(parameter_bound[1])

            # sampling
            if lower > 0 and (upper / lower) >= 100.0:
                # log sampling
                log_lower = np.log10(lower)
                log_upper = np.log10(upper)
                val = 10 ** (log_lower + sample_unit[i, j] * (log_upper - log_lower))
            else:
                # linear sampling
                val = lower + sample_unit[i, j] * (upper - lower)

            p_vals.append(val)
        tasks.append((i, p_vals))
    return tasks


def run_qmc(
    simulator: Callable,
    n_procs: int,
    parameter_bounds: dict,
    n_qmc_samples: int,
):
    parameter_names = list(parameter_bounds.keys())
    sampler = qmc.Halton(d=len(parameter_names), scramble=True, seed=42)
    tasks = sample_parameters(parameter_bounds, n_qmc_samples, sampler)

    with Pool(processes=n_procs) as pool:
        results = [r for r in pool.map(simulator, tasks) if r is not None]

    df = pd.DataFrame(
        results,
        columns=["Idx"] + parameter_names + ["Weighted_Loss", "Pts_in_IQR"],
    )
    return df.sort_values("Weighted_Loss", ascending=True)
