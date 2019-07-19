## Code for reproducing simulations and analyses in "P-splines with an l1 penalty for repeated measures" by Segal, et al. (2018).

The accompanying R package is available at [https://github.com/bdsegal/psplinesl1](https://github.com/bdsegal/psplinesl1).

### Contents

1. `simulation`: Code for the simulations.
    1. `generate_data.R`: Generate data for single group. The simulated dataset presented in the paper is provided in the `psplinesl1` package as `simData`.
    2. `generate_data2.R`: Generate data for two groups. The simulated dataset presented in the paper is provided in the `psplinesl1` package as `simData2groups`.
    3. `plot_data.R`: Plot the simulated datasets.
    4. `l1_estimate.R`: Fit the l1-penalized model to simData.
    5. `l2_estimate.R`: Fit the l2-penalized model to simData.
    6. `bayes_esimate.R`: Fit the Bayesian models. The stan code for fitting models with a laplace, normal, and diffuse prior on the finite order differences in coefficients is in `bayes_lap.stan`, `bayes_norm.stan`, and `bayes_noPen.stan` respectively.
    7. `change_point_simulation_batch.R`: Runs 100 simulations and saves the results in the `batch` folder. To conduct 1,000 simulations, run from the command line and pass in arguments 1-10.
    8. `change_point_simulation_assess.R` Plot the results of the change point simulation.
    9. `l1_estimate_2groups.R`: Fit the l1-penalized model to simData2groups.
    10. `l2_estimate_2groups.R`: Fit the l2-penalized model to simData2groups.
    11. `coverage_prob_simulation.R`: Simulate coverage probability and confidence interval width.
2. `application`: Code for the analysis of electrodermal activity (EDA) data collected as part of a stress study.
    1. `process_data.R`: pre-process raw data.
    2. `analyze_EDA_l1.R`: Fit the l1-penalized model.
    3. `analyze_EDA_l2.R`: Fit the l2-penalized model.
    4. `analyze_EDA_l1_alt.R`: Fit the l1-penalized model with an alternative correlation structure.
    5. `analyze_EDA_l2_alt.R`: Fit the l2-penalized model with an alternative correlation structure.
    6. `bayes_norm`: Folder containing files for fitting Bayesian model with normal prior on the "unpenalized" random effect coefficients. `EDA_bayes_sqrt_lap_norm.R` sets up and runs the model. The stan code is in `EDA_lap_sqrt_norm.stan`.
    7. `bayes_cauchy`: Folder containing files for fitting Bayesian model with Cauchy prior on the "unpenalized" random effect coefficients. `EDA_bayes_sqrt_lap_cauchy.R` sets up and runs the model. The stan code is in `EDA_lap_sqrt_cauchy.stan`.
    8. `EDA_bayes_plots.R`: Plot results from Bayesian models.
3. `inference_knots_plots`: Folder containing scripts to check the proof of observation 1 and make plots demonstrating our approach for approximate inference.
    1. `observation1_check.R`: Check proof of observation 1.
    2. `approx_inference_plots.R`: Plots demonstrating our approximate inference approach.

### References
Segal, B. D., Elliott M., Braun T., and Jiang, H. (2018).  P-splines with an l1 penalty for repeated measures. Electronic Journal of Statistics. 12(2), 3554-3600. [doi.org/10.1214/18-EJS1487](https://doi.org/10.1214/18-EJS1487)