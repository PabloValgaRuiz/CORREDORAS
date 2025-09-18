import numpy as np
import pandas as pd
from statsmodels.tsa.api import VAR
from statsmodels.tsa.stattools import grangercausalitytests


# results_dir = 'results_GG'
# dirdatain = '../data_GG/GAM'
# fit_type = 'fit_301_0.01'

results_dir = 'results_basa'
dirdatain = 'GAM_species'
fit_type = 'fit_140_0.01'


# --------------------------
# 1. Example: create synthetic data
# --------------------------
np.random.seed(42)
n_obs = 200

# True VAR(1) system with causality X -> Y
e = np.random.normal(size=(n_obs, 2))
X = np.zeros((n_obs, 2))
for t in range(1, n_obs):
    X[t, 0] = 0.7 * X[t-1, 0] + e[t, 0]          # X depends on past X
    X[t, 1] = 0.3 * X[t-1, 0] + 0.6 * X[t-1, 1] + e[t, 1]  # Y depends on past X and Y

df = pd.DataFrame(X, columns=["X", "Y"])

# --------------------------
# 2. Fit initial VAR model
# --------------------------
model = VAR(df)
lag_order = model.select_order(maxlags=5).aic
res = model.fit(lag_order)

# Store fitted parameters
coef_matrices = res.coefs
intercept = res.intercept
residuals = res.resid
initial_values = df.values[:lag_order]
print("Residuals shape:", residuals.shape)
print(residuals)

# --------------------------
# 3. Bootstrap function
# --------------------------
def bootstrap_var_once():
    # Resample residuals with replacement, keep cross-series alignment
    boot_resids = residuals[np.random.randint(0, len(residuals), size=len(residuals))]
    
    # Create array for synthetic data
    synthetic = np.zeros_like(df.values)
    synthetic[:lag_order] = initial_values  # seed with initial lags
    
    # Generate synthetic series
    for t in range(lag_order, len(df)):
        y_pred = intercept.copy()
        for L in range(lag_order):
            y_pred += coef_matrices[L] @ synthetic[t - L - 1]
        synthetic[t] = y_pred + boot_resids[t - lag_order]
    
    return pd.DataFrame(synthetic, columns=df.columns)

# --------------------------
# 4. Run bootstrap and Granger
# --------------------------
B = 500
edge_X_to_Y_count = 0
alpha = 0.05

for _ in range(B):
    boot_df = bootstrap_var_once()
    # Run Granger test for X -> Y
    gc_res = grangercausalitytests(boot_df[["Y", "X"]], maxlag=lag_order, verbose=False)
    # Here we just check smallest p-value across lags
    min_pval = min(gc_res[lag][0]["ssr_ftest"][1] for lag in gc_res)
    if min_pval < alpha:
        edge_X_to_Y_count += 1

stability = edge_X_to_Y_count / B
print(f"Edge X -> Y appears in {stability*100:.1f}% of bootstraps")
