from scipy.integrate import solve_ivp
import numpy as np
import src.constants as constants
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def ode_system(t, y, params, M_adv, q_S_pos, CSF_at_t, M_laplacian, DcLV_mask):
    phi, P_bot, D_L, P_DcLV, phi_L = params
    CD, CL = y[: M_laplacian.shape[0]], y[M_laplacian.shape[0] :]

    eta = (phi * phi_L + phi_L) / constants.phi_CSF
    D_D = 3.5e-5
    beta = 0.228 / (phi * phi_L)
    P_vec = P_bot
    exchange = beta * (CL * phi - CD)

    dCDdt = (
        D_D * M_laplacian.dot(CD)
        + constants.INV_L_DURA * P_vec * (eta * CSF_at_t - CD)
        + exchange
        - q_S_pos * CD
    )
    clearance_vec = P_DcLV * DcLV_mask
    dCLdt = (
        D_L * M_laplacian.dot(CL)
        + (1.0 / phi_L) * M_adv.dot(CL)
        - exchange
        + q_S_pos * CD
        - clearance_vec * CL
    )
    return np.concatenate([dCDdt, dCLdt])


# Objective function
def compute_metrics(solution, dura_data: dict):
    if solution is None:
        return np.inf, 0

    # IQR objective function
    iqr_widths = dura_data["PC75"][:, 1:] - dura_data["PC25"][:, 1:]

    iqr_widths = np.maximum(iqr_widths, 1e-2)

    # Error = (Model - Median)/(IQR widths)
    residuals = (solution[:, 1:] - dura_data["median"][:, 1:]) / iqr_widths
    weighted_loss = np.mean(residuals**2)

    # Diagnostic
    pts_in_iqr = np.sum(
        (solution[:, 1:] >= dura_data["PC25"][:, 1:])
        & (solution[:, 1:] <= dura_data["PC75"][:, 1:]),
    )

    return weighted_loss, pts_in_iqr


def simulate(
    matrices: tuple,
    dura_data: dict,
    fitted_csf_functions,
    Z_coords: np.ndarray,
    args,
):
    # Connectivity matrices
    M_laplacian, M_G_dense = matrices

    # Mask for deep cervical lymphatic vessels (DcLV)
    z_mid = 200
    DcLV_mask = np.where(Z_coords < z_mid, 1.0, 0.0)
    _, node_bot_idx = np.argmax(Z_coords), np.argmin(Z_coords)
    idx, params = args

    # Load parameters
    phi, _, _, _, phi_L = params
    if not (0.1 <= phi * phi_L <= 0.35):
        return None

    N_NODES = M_laplacian.shape[0]

    try:
        A = (
            constants.K_DARCY_FIXED * M_laplacian
            - sp.eye(N_NODES) * constants.GAMMA_FIXED
        )
        b = -constants.GAMMA_FIXED * constants.H_D_FIXED * np.ones(N_NODES)
        # No-Flux on Top:
        A_lil = A.tolil()
        # pressure 0 bottom.
        A_lil[node_bot_idx, :] = 0
        A_lil[node_bot_idx, node_bot_idx] = 1
        b[node_bot_idx] = 0
        HL = spla.spsolve(A_lil.tocsr(), b)

        Delta_H = HL.reshape(1, -1) - HL.reshape(-1, 1)
        Q_mat = constants.K_DARCY_FIXED * M_G_dense * Delta_H
        M_adv = sp.csr_matrix(
            np.maximum(0, Q_mat) - np.diag(np.sum(np.maximum(0, -Q_mat), axis=1)),
        )
        q_S_pos = constants.GAMMA_FIXED * np.maximum(0, constants.H_D_FIXED - HL)

        # Solve
        sol = solve_ivp(
            lambda t, y: ode_system(
                t,
                y,
                params,
                M_adv,
                q_S_pos,
                np.array([f(t) for f in fitted_csf_functions]),
                M_laplacian,
                DcLV_mask,
            ),
            (0, 72),
            np.zeros(N_NODES * 2),
            t_eval=constants.TIME_POINTS,
            method="Radau",
            rtol=1e-3,
        )
        solution = sol.y[:N_NODES, :] + sol.y[N_NODES:, :] if sol.success else None
        loss, pts = compute_metrics(solution, dura_data)
        return [idx] + list(params) + [loss, pts]

    except Exception as e:
        print(f"Errore nella simulazione: {e}")
        return None
