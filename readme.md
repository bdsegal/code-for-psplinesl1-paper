## Code for reproducing simulations and analyses in "P-splines with an l1 penalty for repeated measures" by Segal, et al. (2017).

The accompanying R package is available at [https://github.com/bdsegal/psplinesl1](https://github.com/bdsegal/psplinesl1).

### Contents

1. `simulation`: Code for the simulations.
    1. `generate_data.R`: Generate data. The simulated dataset presented in the paper is provided in the `psplinesl1` package as `simData`.
    2. `plot_data.R`: Plot the simulated dataset.
    3. `l1_estimate.R`: Fit the l1-penalized model.
    4. `l2_estimate.R`: Fit the l2-penalized model.
    5. `bayes_esimate.R`: Fit the Bayesian models. The stan code for fitting models with a laplace, normal, and diffuse prior on the finite order differences in coefficients is in `bayes_lap.stan`, `bayes_norm.stan`, and `bayes_noPen.stan` respectively.
    6. `change_point_simulation_batch.R`: Runs 100 simulations and saves the results in the `batch` folder. To conduct 1,000 simulations, run from the command line and pass in arguments 1-10.
    7. `change_point_simulation_assess.R` Plot the results of the change point simulation.
2. `application`: Code for the analysis of electrodermal activity (EDA) data collected as part of a stress study.
    1. `process_data.R`: pre-process raw data.
    2. `analyze_EDA_l1.R`: Fit the l1-penalized model.
    3. `analyze_EDA_l2.R`: Fit the l2-penalized model.
    4. `analyze_EDA_l2_alt.R`: Fit the l2-penalized model with an alternative correlation structure.
    5. `bayes_norm`: Folder containing files for fitting Bayesian model with normal prior on the "unpenalized" random effect coefficients. `EDA_bayes_sqrt_lap_norm.R` sets up and runs the model. The stan code is in `EDA_lap_sqrt_norm.stan`.
    6. `bayes_cauchy`: Folder containing files for fitting Bayesian model with Cauchy prior on the "unpenalized" random effect coefficients. `EDA_bayes_sqrt_lap_cauchy.R` sets up and runs the model. The stan code is in `EDA_lap_sqrt_cauchy.stan`.
    7. `EDA_bayes_plots.R`: Plot results from Bayesian models.
3. `inference_knots_plots`: Folder containing scripts to check the proof of observation 1 and make plots demonstrating our approach for approximate inference.
    1. `observation1_check.R`: Check proof of observation 1.
    2. `approx_inference_plots.R`: Plots demonstrating our approximate inference approach.

### References
Segal, B. D., Elliott, M. R., Braun, T., and Jiang, H. (2017). P-splines with an l1 penalty for repeated measures. Available at [https://arxiv.org/abs/1707.08933](https://arxiv.org/abs/1707.08933).