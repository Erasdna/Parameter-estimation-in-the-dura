import numpy as np
import scipy.sparse as sp
from scipy.optimize import curve_fit
from functools import partial


# CSF fitting
def gamma_func(t, A, alpha, beta):
    t_safe = np.maximum(t, 1e-6)
    return A * (t_safe**alpha) * np.exp(-t_safe / beta)


def double_exp_func(t, A, alpha, beta):
    return A * (np.exp(-alpha * t) - np.exp(-beta * t))


def setup_matrices(df_edges, df_nodes):
    # Setup of matrices
    M_G_dense = np.zeros((df_nodes.shape[0], df_nodes.shape[0]))
    M_lap_dok = sp.dok_matrix((M_G_dense.shape[0], M_G_dense.shape[1]))

    # Extract previously computed quantities
    Vol_vec = df_nodes["Volume_dm3"].values
    area_over_distance = df_edges["Contact_Area_dm2"] / df_edges["Distance_dm"]
    source_ids = df_edges["Source_ROI"].astype(int) - 1
    target_ids = df_edges["Target_ROI"].astype(int) - 1

    # Construct connectivity matrices
    M_lap_dok[source_ids, target_ids] = area_over_distance / Vol_vec[source_ids]
    M_lap_dok[target_ids, source_ids] = area_over_distance / Vol_vec[target_ids]
    M_G_dense[source_ids, target_ids] = M_G_dense[
        target_ids,
        source_ids,
    ] = area_over_distance
    M_lap_dok[np.arange(M_G_dense.shape[0]), np.arange(M_G_dense.shape[0])] = -np.sum(
        M_lap_dok,
        axis=1,
    )

    return M_G_dense, M_lap_dok.tocsr()


def fit_CSF_data(csf_data, t_data):
    # FITTING CSF DATA (Gamma vs Double Exp)
    print("\n>>> Fitting CSF functions for each patch..")

    best_cT_funcs = []
    fit_stats = {"gamma": 0, "double_exp": 0, "fallback_interp": 0}

    csf_data = np.where(csf_data < 0, 0, csf_data)
    for i in range(csf_data.shape[0]):
        y_data = csf_data[i, :]

        # If no tracer, =0
        if np.max(y_data) <= 1e-6:
            best_cT_funcs.append(lambda t: 0.0)
            continue

        err_g, err_d = np.inf, np.inf
        popt_g, popt_d = None, None

        # fit GAMMA
        try:
            popt_g, _ = curve_fit(
                gamma_func,
                t_data,
                y_data,
                p0=[np.max(y_data), 2.0, 10.0],
                bounds=(0, np.inf),
                maxfev=2000,
            )
            err_g = np.sum((gamma_func(t_data, *popt_g) - y_data) ** 2)
        except RuntimeError:
            continue  # Fall back on double exponential
        except Exception as e:
            print(f"Error during CSF fitting {e}")

        # fit DOUBLE EXPONENTIAL
        try:
            popt_d, _ = curve_fit(
                double_exp_func,
                t_data,
                y_data,
                p0=[np.max(y_data), 0.01, 0.1],
                bounds=(0, np.inf),
                maxfev=2000,
            )
            err_d = np.sum((double_exp_func(t_data, *popt_d) - y_data) ** 2)
        except RuntimeError:
            continue  # Fall back on linear function
        except Exception as e:
            print(f"Error during CSF fitting {e}")

        # Best fit
        if err_g < err_d and err_g != np.inf:
            best_cT_funcs.append(
                partial(gamma_func, A=popt_g[0], alpha=popt_g[1], beta=popt_g[2]),
            )
            fit_stats["gamma"] += 1
        elif err_d != np.inf:
            best_cT_funcs.append(
                partial(double_exp_func, A=popt_d[0], alpha=popt_d[1], beta=popt_d[2]),
            )
            fit_stats["double_exp"] += 1
        else:
            # If both fails, linear interpolation
            best_cT_funcs.append(lambda t, y=y_data: np.interp(t, t_data, y))
            fit_stats["fallback_interp"] += 1

    print(
        f"    [+] Fitting results: {fit_stats['gamma']} Gamma, \
            {fit_stats['double_exp']} DoubleExp, {fit_stats['fallback_interp']} Linear.",
    )

    return best_cT_funcs
