# Seven-method BIC comparison v1

```json
{
  "protocol": "lambert_bic_comparison_v1",
  "baseline": "outputs/analysis/lambert_main_truba_v2_a01",
  "baseline_manifest_sha256": "a051d00456f0f2565e12974efd971bafb4c296001acb9ce600971b674071e363",
  "repetitions": 300,
  "p": [125,500], "rho": [0,0.3,0.6,0.9], "snr": [1,4],
  "profile": "moderate",
  "methods": ["lambert_c1","lasso","ridge","elastic_net","mcp","scad"],
  "workers": 6,
  "reproduction_tolerance": 0.000001,
  "smoke_stream_offset": 50000,
  "smoke_tasks": 7,
  "smoke_n": 40, "smoke_p": [20,60], "smoke_n_test": 40,
  "smoke_lambda_points": 4,
  "progress_seconds": 30
}
```
