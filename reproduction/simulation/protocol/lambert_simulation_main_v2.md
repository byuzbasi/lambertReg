# Lambert main simulation v2

```json
{
  "protocol": "lambert_simulation_main_v2",
  "n_train": 250,
  "p": [125, 500],
  "rho": [0, 0.3, 0.6, 0.9],
  "snr": [0.25, 1, 4],
  "profiles": ["equal", "moderate"],
  "support_size": 10,
  "support_indices": [1,2,3,4,5,6,7,8,9,10],
  "signs": [1,1,1,1,1,-1,-1,-1,-1,-1],
  "magnitudes_equal": [1,1,1,1,1,1,1,1,1,1],
  "magnitudes_moderate": [3,3,2,2,2,2,1,1,1,1],
  "magnitudes_extreme": [4,3,2,1.5,1,0.7,0.5,0.3,0.2,0.1],
  "beta_l2_norm": 1,
  "n_test": 2000,
  "repetitions": 300,
  "folds": 5,
  "lambda_points": 60,
  "lambda_fraction_log10_limits": [0,-3],
  "ridge_factor_log10_limits": [3,-4],
  "elastic_net_alpha": [0,0.25,0.5,0.75,1],
  "lambert_shapes": 1,
  "mcp_gamma": 3,
  "scad_a": 3.7,
  "methods": ["lambert_c1","lasso","ridge","elastic_net","mcp","scad"],
  "master_seed": 2026100401,
  "stream_offsets": {"main":0,"pilot":20000,"smoke":30000},
  "stress": {"p":[125,500],"rho":[0,0.9],"snr":1,"profile":"extreme"},
  "pilot_cases": [
    {"p":125,"rho":0,"snr":0.25,"profile":"equal"},
    {"p":125,"rho":0.3,"snr":1,"profile":"moderate"},
    {"p":500,"rho":0.6,"snr":4,"profile":"moderate"},
    {"p":500,"rho":0.9,"snr":0.25,"profile":"equal"},
    {"p":125,"rho":0,"snr":1,"profile":"extreme"},
    {"p":500,"rho":0.9,"snr":1,"profile":"extreme"}
  ],
  "max_coordinate_sweeps_per_start": 10000,
  "normalized_kkt_tolerance": 1e-7,
  "support_reporting_tolerance": 1e-8,
  "implementation": {
    "ncvreg_eps":1e-12,"ncvreg_max_iter":10000,
    "glmnet_thresh":1e-16,"glmnet_maxit":1000000,
    "progress_seconds":30,"tie_tol":1e-10,"objective_increase_tol":1e-9,
    "curvature_join_tol":1e-8,"curvature_negative_tol":1e-8,
    "standardization_tol":1e-12,"execution_order":"replicate_then_scenario",
    "version":"lambert_simulation_main_v2"
  }
}
```
