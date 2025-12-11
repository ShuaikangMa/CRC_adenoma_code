import math
import pandas as pd
import torch
import warnings
from tqdm import tqdm
import numpy as np
import json
from sklearn.linear_model import LinearRegression
from scipy.optimize import minimize
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import ScalarFormatter

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def equation(x, power_par):
    if isinstance(power_par, torch.Tensor):
        if power_par.dim() == 1:
            a, b = power_par[0], power_par[1]
            return a * torch.pow(x, b)
        else:
            a = power_par[:, 0].unsqueeze(1)
            b = power_par[:, 1].unsqueeze(1)
            return a * torch.pow(x, b)
    else:  # numpy ndarray
        if power_par.ndim == 1:
            a, b = power_par[0], power_par[1]
            return a * np.power(x, b)
        else:
            a = power_par[:, 0][:, None]
            b = power_par[:, 1][:, None]
            return a * np.power(x, b)

def power_equation_base(x, y):
    try:
        x = np.array(x, dtype=float).reshape(-1, 1)
        y = np.array(y, dtype=float)
        log_x = np.log(x + 1e-6)
        log_y = np.log(y + 1e-6)
        model = LinearRegression().fit(log_x, log_y)
        a = np.exp(model.intercept_)
        b = model.coef_[0]
        return a, b
    except:
        return None

def zero_inflated_power_law(x, y, max_retries=10):
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)

    mask_zero = (y == 0)
    mask_nonzero = ~mask_zero
    x0 = x[mask_zero]
    x1 = x[mask_nonzero]
    y1 = y[mask_nonzero]

    def joint_neg_log_likelihood(par):
        beta1 = par[:2]   # logit
        beta2 = par[2:4]  # power law
        sigma = abs(par[4])
        x0_clip = np.clip(x0, 1e-5, None)
        x1_clip = np.clip(x1, 1e-5, None)

        logit0 = beta1[0] * x0_clip**beta1[1]
        logit1 = beta1[0] * x1_clip**beta1[1]

        p0 = np.exp(logit0) / (1 + np.exp(logit0))
        p1 = 1 / (1 + np.exp(logit1))

        loglik_zeros = np.sum(np.log(np.clip(p0, 1e-10, 1)))
        y1_hat = beta2[0] * x1_clip**beta2[1]
        loglik_nonzeros = np.sum(np.log(np.clip(p1, 1e-10, 1))) + \
                          np.sum(-0.5 * ((y1 - y1_hat) / sigma)**2 - np.log(sigma) - 0.5 * np.log(2 * np.pi))
        return - (loglik_zeros + loglik_nonzeros)

    for attempt in range(max_retries):
        init = np.random.uniform(0.05, 0.5, size=5)
        try:
            result = minimize(joint_neg_log_likelihood, init, method='Nelder-Mead',
                              options={'maxiter': 10000, 'disp': False})
            if result.success:
                beta2 = result.x[2:4]
                return beta2[0], beta2[1]
        except:
            continue
    return None


def equation_all(x, y, max_retry=100):
    if np.sum(np.array(y) == 0) > 10:
        result = zero_inflated_power_law(x, y, max_retries=max_retry)
        if result is not None:
            # print("Fitted zero-inflated power law parameters:", result)
            return result  # (a, b)
        else:
            y_nonzero = np.array(y)[np.array(y) > 0]
            x_nonzero = np.array(x)[np.array(y) > 0]
            result = power_equation_base(x_nonzero, y_nonzero)
            if result is not None:
                return result
    else:
        # print("Fitting power law equation without zero inflation.")
        result = power_equation_base(x, y)
        if result is not None:
            return result
    return None


def logsumexp(v):
    vm = torch.max(v)
    return torch.log(torch.sum(torch.exp(v - vm))) + vm

def make_positive_definite(cov, max_trials=10, base_eps=1e-5):
    eps = base_eps
    identity = torch.eye(cov.shape[0], device=cov.device)
    for _ in range(max_trials):
        try:
            L = torch.linalg.cholesky(cov + eps * identity)
            return cov + eps * identity
        except RuntimeError:
            eps *= 10
    raise RuntimeError("Failed to make covariance matrix positive definite.")

def get_SAD1_covmatrix(par, n):
    eps = 1e-3
    phi = par[0].clamp(-0.99, 0.99)
    gamma = par[1].clamp(min=1e-3)

    indices = torch.arange(1, n + 1, dtype=torch.float32).to(device)
    diag_elements = (1 - torch.pow(phi, 2 * indices)) / (1 - torch.pow(phi, 2))
    sigma = torch.diag(diag_elements) * (gamma ** 2)

    return sigma + eps * torch.eye(n, dtype=torch.float32, device=device)

def robust_kmeans(X, k, min_size=2, max_trials=1000):

    if isinstance(X, pd.DataFrame):
        data_values = X.values
        mean_values = X.mean(axis=1).to_frame(name='mean').values
    else:
        data_values = X
        mean_values = np.mean(X, axis=1).reshape(-1, 1)

    fallback_km = None
    fallback_seed = None

    for i in range(max_trials):

        km = KMeans(n_clusters=k, random_state=np.random.randint(10000))
        labels = km.fit_predict(mean_values)
        bincounts = np.bincount(labels, minlength=k)

        print(f"try {i + 1} ，seed {seed}")

        if np.all(bincounts >= min_size):
            print(f"Yes")
            break

        fallback_km = km
        fallback_seed = seed
    else:
        print(f"No")
        km = fallback_km
        seed = fallback_seed
        labels = km.labels_

    mean_vectors = np.zeros((k, data_values.shape[1]))
    for c in range(k):
        mean_vectors[c, :] = data_values[labels == c, :].mean(axis=0)

    return mean_vectors, labels, seed


def biget_par_int(X, k, times1, times2, n1, n2):
    n, d = X.shape
    X1 = X[:, :n1]
    X2 = X[:, n1:n1 + n2]
  
    cov_int = [0.5, 0.03, 0.6, 0.03]
  
    if k <= 10:
        min_size = 10
    if k > 10 and k <= 15:
        min_size = 5
    if k > 15:
        min_size = 3
    mean_vectors, labels, seed = robust_kmeans(X.cpu().numpy(), k, min_size=min_size)
    times1_np = times1.cpu().numpy()
    times2_np = times2.cpu().numpy()
    fit1_list = []
    fit2_list = []

    for c in range(k):
        X_c = X[labels == c, :].cpu().numpy()
        zero_ratio = np.sum(X_c == 0) / X_c.size

        if zero_ratio <= 0.1:
            print(f"Cluster {c+1} has low zero ratio: {zero_ratio:.2f}, using power law fit.")
            fit1_c = equation_all(times1_np, mean_vectors[c, :n1])
            fit2_c = equation_all(times2_np, mean_vectors[c, n1:n1 + n2])
        else:
            print(f"Cluster {c+1} has high zero ratio: {zero_ratio:.2f}, using robust fitting.")
            X1_c = X_c[:, :n1]
            X2_c = X_c[:, n1:n1 + n2]

            fit1_all = np.array([equation_all(times1_np, X1_c[i, :]) for i in range(X1_c.shape[0])])
            fit2_all = np.array([equation_all(times2_np, X2_c[i, :]) for i in range(X2_c.shape[0])])

            fit1_vals = np.array([equation(times1_np, fit1_all[i, :]) for i in range(fit1_all.shape[0])]).T
            fit2_vals = np.array([equation(times2_np, fit2_all[i, :]) for i in range(fit2_all.shape[0])]).T

            mean_fit1_vals = fit1_vals.mean(axis=1)
            mean_fit2_vals = fit2_vals.mean(axis=1)

            fit1_c = equation_all(times1_np, mean_fit1_vals)
            fit2_c = equation_all(times2_np, mean_fit2_vals)

        fit1_list.append(fit1_c)
        fit2_list.append(fit2_c)

    fit1 = np.array(fit1_list).T
    fit2 = np.array(fit2_list).T

    print("initial params：", fit1.shape, fit2.shape)

    return {
        'initial_cov_params': cov_int,
        'initial_mu_params': np.hstack((fit1, fit2)).flatten(),
        'initial_probibality': np.bincount(labels) / n
    }


def sigmoid(x):
    return 1 / (1 + torch.exp(-x))

def get_SAD1_covmatrix(par, n):
    eps = 1e-3
    phi = par[0].clamp(-0.99, 0.99)
    gamma = par[1].clamp(min=1e-3)

    indices = torch.arange(1, n + 1, dtype=torch.float32).to(device)
    diag_elements = (1 - torch.pow(phi, 2 * indices)) / (1 - torch.pow(phi, 2))
    sigma = torch.diag(diag_elements) * (gamma ** 2)

    return sigma + eps * torch.eye(n, dtype=torch.float32, device=device)

def normal_logpdf(x, mu, sigma):
    var = sigma ** 2
    log_pdf = -0.5 * torch.log(2 * math.pi * var) - 0.5 * ((x - mu) ** 2) / var
    return log_pdf


def zero_inflated_log_likelihood(x, mu_params, sigma_params, zero_params, times):

    times_clip = torch.clamp(times, min=1e-5)
    zero_logits = zero_params[0] + zero_params[1] * torch.log(times_clip)
    p0 = torch.sigmoid(zero_logits)  # [n_features]

    mu = equation(times, mu_params)  # [n_features]

    if sigma_params.ndim == 0:
        sigma = sigma_params * torch.ones_like(mu)
    else:
        sigma = sigma_params

    mask_zero = (x == 0).float()
    mask_nonzero = 1 - mask_zero

    log_prob_zero = torch.log(torch.clamp(p0, min=1e-10))

    log_prob_nonzero = torch.log(torch.clamp(1 - p0, min=1e-10)) + normal_logpdf(x, mu,
                                                                                 sigma)

    loglik = torch.sum(mask_zero * log_prob_zero + mask_nonzero * log_prob_nonzero, dim=1)
    return loglik


def biQ_function(par, prob_log, omega_log, X, k, n1, n2, times1, times2):
    n = X.shape[0]
    n1, n2, k = map(int, [n1, n2, k])
    X1 = X[:, :n1]
    X2 = X[:, n1:n1 + n2]

    par_mu = par[8:]
    zero_params1 = par[4:6]
    zero_params2 = par[6:8]

    cov1 = get_SAD1_covmatrix(par[:2], n1)
    cov2 = get_SAD1_covmatrix(par[2:4], n2)

    sigma1 = torch.sqrt(torch.diag(cov1))
    sigma2 = torch.sqrt(torch.diag(cov2))

    mu1_mat = torch.column_stack((par_mu[:k], par_mu[2 * k:3 * k]))  # [k,2]
    mu2_mat = torch.column_stack((par_mu[k:2 * k], par_mu[3 * k:4 * k]))  # [k,2]

    mvn_log1 = torch.stack([
        zero_inflated_log_likelihood(X1, mu1_mat[i], sigma1, zero_params1, times1)
        for i in range(k)
    ], dim=1)

    mvn_log2 = torch.stack([
        zero_inflated_log_likelihood(X2, mu2_mat[i], sigma2, zero_params2, times2)
        for i in range(k)
    ], dim=1)

    mvn_log = mvn_log1 + mvn_log2
    tmp = prob_log.T + mvn_log - omega_log
    Q = -torch.sum(tmp * torch.exp(omega_log))
    return Q


def bifun_clu(data1, data2, k, Time1=None, Time2=None, trans=torch.log10, initial_pars=None,
              iter_max=500, init_lr=None):
    data1 = data1[:, torch.argsort(torch.sum(data1, axis=0))]
    data2 = data2[:, torch.argsort(torch.sum(data2, axis=0))]

    data = torch.hstack((data1, data2))

    n, d = data.shape
    n1, n2 = data1.shape[1], data2.shape[1]
    epsilon = 1000
    iter = 0

    data = trans(data + 1)
    X1 = data[:, :n1]
    X2 = data[:, n1:n1 + n2]


    times1 = trans(torch.sum(data1, axis=0) + 1).to(device)
    times2 = trans(torch.sum(data2, axis=0) + 1).to(device)
    if initial_pars is None:
        par_int = biget_par_int(data, k, times1, times2, n1, n2)
        prob_int = np.ones(k) / k
        cov_int = par_int['initial_cov_params']
        mu_int = par_int['initial_mu_params']

        zero_inflation_init = np.array([-0.1, -0.1, -0.1, -0.1])
        initial_pars = np.hstack((cov_int, zero_inflation_init, mu_int))
        par = torch.tensor(initial_pars, dtype=torch.float32).to(device)
    else:
        cov_par = initial_pars["cov_par"]
        zero_inflation_params = initial_pars["zero_inflation_params"]
        mu_par = initial_pars["mu_par"]
        prob_int = initial_pars['probility']
        par = torch.tensor(
            cov_par +  # [0:4]
            zero_inflation_params['zero_params1'] +  # [4:6]
            zero_inflation_params['zero_params2'] +  # [6:8]
            sum([sum(row, []) for row in mu_par], [])  # [8:]
            , dtype=torch.float32)
        # print(par)


    par.requires_grad = True
    prob_logi = torch.log(torch.tensor(prob_int, dtype=torch.float32).to(device))

    initial_lr = init_lr
    optimizer = torch.optim.AdamW([par], lr=initial_lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor=0.3, patience=10)

    while abs(epsilon) > 0.1 and iter < iter_max:
        par_mui = par[8:]
        cov1i = get_SAD1_covmatrix(par[:2], n1).to(device).detach()
        cov2i = get_SAD1_covmatrix(par[2:4], n2).to(device).detach()

        zero_params1 = par[4:6].to(device).detach()
        zero_params2 = par[6:8].to(device).detach()

        mu1_mati = torch.column_stack((par_mui[:k], par_mui[2 * k:3 * k])).to(device).detach()
        mu2_mati = torch.column_stack((par_mui[k:2 * k], par_mui[3 * k:4 * k])).to(device).detach()

        times1_ = times1.to(device).detach()
        times2_ = times2.to(device).detach()

        sigma1 = torch.sqrt(torch.diag(cov1i))  # [n_features]
        sigma2 = torch.sqrt(torch.diag(cov2i))  # [n_features]

        mvn_log1i = torch.stack([
            zero_inflated_log_likelihood(X1, mu1_mati[i], sigma1, zero_params1, times1_) for i in range(k)
        ], dim=1).to(device).detach()

        mvn_log2i = torch.stack([
            zero_inflated_log_likelihood(X2, mu2_mati[i], sigma2, zero_params2, times2_) for i in range(k)
        ], dim=1).to(device).detach()
        mvn_logi = mvn_log1i + mvn_log2i
        mvni = mvn_logi + prob_logi

        omega_logi = mvni - torch.logsumexp(mvni, dim=1, keepdim=True)
        omegai = torch.exp(omega_logi)

        LL_mem = biQ_function(par, prob_logi, omega_logi, data, k, n1, n2, times1, times2)

        nk = omegai.sum(dim=0)
        alpha = 1 / 2 * k
        prob = (nk / nk.sum()) + alpha
        prob = prob / prob.sum()
        prob_logi = torch.log(prob)

        for _ in tqdm(range(25), desc="Training Progress"):
            optimizer.zero_grad()
            Q = biQ_function(par, prob_logi, omega_logi, data, k, n1, n2, times1, times2)
            Q.backward()
            torch.nn.utils.clip_grad_norm_(par, max_norm=1.0)
            optimizer.step()

        par_hat = par
        par = par_hat

        LL_next = biQ_function(par, prob_logi, omega_logi, data, k, n1, n2, times1, times2)
        epsilon = LL_next - LL_mem
        LL_mem = LL_next
        scheduler.step(LL_next)
        iter += 1

        print(f"Iter: {iter}, Log-Likelihood: {LL_next}")
        final_assignments = torch.argmax(omegai, dim=1).cpu().numpy()

    assignments = torch.argmax(omegai, dim=1).cpu().numpy()
    omega_np = omegai.cpu().numpy()
    times1_np = times1.cpu().numpy()
    times2_np = times2.cpu().numpy()
    X_np = data.cpu().numpy()
    prob_np = torch.exp(prob_logi).cpu().numpy()

    def restructure_mu(mu_par, k):
        mu_par = np.asarray(mu_par)
        mu_re = np.zeros((k, 2, 2))
        mu_re[:, 0, 0] = mu_par[:k]
        mu_re[:, 0, 1] = mu_par[k:2 * k]
        mu_re[:, 1, 0] = mu_par[2 * k:3 * k]
        mu_re[:, 1, 1] = mu_par[3 * k:4 * k]
        return mu_re

    mu_param_raw = par_hat[8:].detach().cpu().numpy()
    mu_param = restructure_mu(mu_param_raw, k)

    import pandas as pd
    X_df = pd.DataFrame(X_np)
    X_df['cluster'] = assignments
    X_df['name'] = [f"sample_{i}" for i in range(X_df.shape[0])]  # optional
    X_df = X_df.reset_index(drop=True)

    cov_param = par_hat[:4].detach().cpu().numpy()
    par = par.detach().cpu().numpy()
    LL_next = LL_next.detach().cpu().numpy()
    AIC = -2 * LL_next + 2 * (len(par_hat) + k - 1)
    BIC = -2 * LL_next + np.log(n) * (len(par) + k - 1)
    zero_inflation_params = {
        'zero_params1': par_hat[4:6].detach().cpu().numpy(),
        'zero_params2': par_hat[6:8].detach().cpu().numpy()
    }

    cluster_counts = np.bincount(assignments, minlength=k)
    return {
        'cluster_counts': cluster_counts,
        'cluster_number': k,
        'log-likelihood': LL_next.item(),
        'AIC': AIC.item(),
        'BIC': BIC.item(),
        'cov_par': cov_param,
        'mu_par': mu_param,
        'probability': prob_np,
        'omega': omega_np,
        'cluster': X_df,
        'cluster2': {'data': X_np, 'cluster': assignments},
        'Time1': times1_np,
        'Time2': times2_np,
        'original_data': X_np,
        'zero_inflation_params': zero_inflation_params
    }


for k in range(2,16):

    result = bifun_clu(
        data1, data2, k=k, init_lr=1e-3
    )

    def convert_ndarray_to_list(obj):
        if isinstance(obj, dict):
            return {k: convert_ndarray_to_list(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_ndarray_to_list(i) for i in obj]
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return obj


    result_dict = {
        'cluster_counts': result['cluster_counts'].tolist(),
        'cluster_number': int(result['cluster_number']),
        'log-likelihood': float(result['log-likelihood']),
        'BIC': float(result['BIC']),
        'AIC': float(result['AIC']),
        'cov_par': np.array(result['cov_par']),
        'mu_par': np.array(result['mu_par']),
        'probability': np.array(result['probability']),
        'omega': np.array(result['omega']),
        'cluster': np.array(result['cluster']),
        'Time1': np.array(result['Time1']),
        'Time2': np.array(result['Time2']),
        'original_data': np.array(result['original_data']),
        'zero_inflation_params': result['zero_inflation_params']
    }

    result_dict = convert_ndarray_to_list(result_dict)
    save_path = f'k{k}.json'
    with open(save_path, 'w') as f:
        json.dump(result_dict, f, indent=2)

    def bifun_clu_plot_py(result, degree=0.6, show_legend=False,
                          color_L="#FF6B6B", color_R="#4ECDC4",
                          figsize=(14, 8), save_path=None):
        mu_par = result['mu_par']
        X = result['original_data']
        Time1, Time2 = result['Time1'], result['Time2']
        cluster = result['cluster']['cluster'].values
        k = result['cluster_number']

        n_samples = X.shape[0]
        n1, n2 = len(Time1), len(Time2)

        times1_new = np.linspace(min(Time1), max(Time1), 100)
        times2_new = np.linspace(min(Time2), max(Time2), 100)

        def power_fit(x, params):
            a, b = params[:, 0:1], params[:, 1:2]
            return a * (x[np.newaxis, :] ** b)

        mu_fit1 = power_fit(times1_new, mu_par[:, :, 0])
        mu_fit2 = power_fit(times2_new, mu_par[:, :, 1])

        def build_line_df(fits, times, side):
            df = pd.DataFrame(fits, columns=times)
            df['cluster'] = np.arange(k)
            df = df.melt(id_vars='cluster', var_name='x', value_name='y')
            df['type'] = side
            return df

        df_line = pd.concat([
            build_line_df(mu_fit1, times1_new, 'L'),
            build_line_df(mu_fit2, times2_new, 'R')
        ], ignore_index=True)

        def build_point_df(X_part, times, prefix, side):
            df = pd.DataFrame(X_part, columns=times)
            df['cluster'] = cluster
            df['sample'] = [f"{prefix}_{i}" for i in range(n_samples)]
            df = df.melt(id_vars=['cluster', 'sample'], var_name='x', value_name='y')
            df['x'] = df['x'].astype(float)
            df['type'] = side
            return df

        df_points = pd.concat([
            build_point_df(X[:, :n1], Time1, 'a', 'L'),
            build_point_df(X[:, n1:], Time2, 'b', 'R')
        ], ignore_index=True)

        cluster_counts = pd.Series(cluster).value_counts().sort_index()
        alpha_map = (cluster_counts / cluster_counts.max() * degree).to_dict()
        df_points['alpha'] = df_points['cluster'].map(alpha_map)

        sns.set_theme(style="white") 
        ncols = min(4, k)
        nrows = int(np.ceil(k / ncols))
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
        axes = axes.flatten()

        for i in range(k):
            ax = axes[i]
            df_p = df_points[df_points['cluster'] == i]
            df_l = df_line[df_line['cluster'] == i]

            ax.set_facecolor("white")
            ax.grid(False)

            sns.scatterplot(
                data=df_p, x='x', y='y', hue='type',
                alpha=df_p['alpha'], palette={'L': color_L, 'R': color_R},
                ax=ax, legend=False, s=12, edgecolor='none'
            )

            sns.lineplot(
                data=df_l, x='x', y='y', hue='type',
                palette={'L': color_L, 'R': color_R},
                ax=ax, legend=False, lw=2.5
            )

            ax.set_title(f"Cluster {i + 1} (n={cluster_counts[i]})", fontsize=11, fontweight='bold')
            ax.set_xlabel("Time", fontsize=10)
            ax.set_ylabel("Abundance", fontsize=10)
            ax.tick_params(axis='both', labelsize=9)

            ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.xaxis.get_major_formatter().set_scientific(True)
            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.yaxis.get_major_formatter().set_scientific(True)

        for j in range(k, len(axes)):
            fig.delaxes(axes[j])

        if show_legend:
            handles = [
                plt.Line2D([], [], marker='o', color='w', markerfacecolor=color_L, label='Left', markersize=6),
                plt.Line2D([], [], marker='o', color='w', markerfacecolor=color_R, label='Right', markersize=6)
            ]
            fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 1.02),
                       ncol=2, frameon=False, fontsize=10)

        plt.tight_layout()
        plt.subplots_adjust(top=0.93)
        if save_path:
            plt.savefig(save_path, bbox_inches='tight', dpi=300)

        # plt.show()
    save_path = f'bifun_cluster_{k}.png'

    if np.all(result['cluster_counts'] > 0):
        bifun_clu_plot_py(result, save_path=save_path, figsize=(16, 10),
                          degree=1)
    else:
        print(f"Cluster {k} contains empty clusters, skipping plot.")
