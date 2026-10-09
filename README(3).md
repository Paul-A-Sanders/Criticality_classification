# Stability Analysis of Pitchfork and Swift--Hohenberg Models

MATLAB code for simulating stochastic pitchfork and Swift--Hohenberg
systems, estimating model coefficients from simulated data, and testing
whether the estimated cubic coefficient distinguishes two dynamical
regimes.

The repository is organised into simulation, analysis, and
combined-workflow folders for each model:

``` text
.
├── Pitchfork/
│   ├── Simulation/
│   ├── Analysis/
│   └── Combine/
└── Swift_Hohenberg/
    ├── Simulation/
    ├── Analysis/
    └── Combine/
```

## Scientific objective

The project tests whether a nonlinear coefficient estimated from noisy
data can distinguish the two signs of the cubic term. The pitchfork
model provides a spatially homogeneous test case; the Swift--Hohenberg
model extends the workflow to spatially distributed data.

The fitted cubic coefficient, (`\hat{s}`{=tex}), is used as a
classification score to distinguish simulations with true (s=+1) from
those with (s=-1). ROC curves and the area under the ROC curve (AUC)
measure classification performance across parameter values and
observation horizons. The Swift--Hohenberg workflow also analyses
dominant wavenumbers, growth rates, and dispersion relations.

## Requirements

The current scripts are written for **MATLAB R2026b** and require:

-   **Parallel Computing Toolbox**, for `parfor` and parallel simulation
    batches.
-   **Statistics and Machine Learning Toolbox**, for ROC/AUC analysis
    with `perfcurve` and coefficient-density plots using `ksdensity`.

The scripts also use standard MATLAB numerical and graphics
functionality, including FFT operations, cell arrays, tiled layouts, and
figure export functions.

No external dataset is required for the standard workflows. The combined
scripts generate the stochastic data when a compatible processed
MAT-file is not available.

## Repository structure

### `Pitchfork/`

-   `Simulation/` contains the stochastic pitchfork simulator.
-   `Analysis/` contains the coefficient-estimation function.
-   `Combine/` contains the script that coordinates simulation, fitting,
    classification analysis, saving, and plotting.

### `Swift_Hohenberg/`

-   `Simulation/` contains the periodic FFT-based Swift--Hohenberg
    simulator.
-   `Analysis/` contains the coefficient-estimation and
    dispersion-analysis function.
-   `Combine/` contains the script that coordinates simulation, fitting,
    classification analysis, saving, and plotting.

The simulation and analysis functions are reusable components; the
combined scripts define the experiment settings and orchestrate the full
workflow. The script names and function signatures below refer to the
current versions of the code.

## Pitchfork workflow

The combined script is `Pitchfork/Combine/test_code_pitchfork_loops.m`.
It simulates stochastic pitchfork trajectories for a sweep of negative
(`\mu`{=tex}) values and both cubic-coefficient signs.

### Main functions

**`simulate_pitchfork`** generates one stochastic trajectory:

``` matlab
[t_vec, u_vec] = simulate_pitchfork( ...
    mu, s_true, tmax, tsteps, noiseType, sigma, InitialScale)
```

**`analyze_pitchfork`** estimates coefficients from a trajectory or a
prefix of one:

``` matlab
[theta_sub, output2, output3] = analyze_pitchfork(u_sub, t_sub);
```

The combined script stores the first two fitted coefficients, with the
second coefficient (`\hat{s}`{=tex}) used as the classification score
(`kc_shat = 2`).

### Workflow and outputs

1.  Define the (`\mu`{=tex}) sweep, cubic-sign classes, noise settings,
    time step, and observation horizons.
2.  Simulate independent seeded trajectories in parallel and reject
    failed, non-finite, or excessively large trajectories.
3.  Fit the model coefficients at the requested observation horizons. If
    too few realisations remain valid at the longest horizon, adaptively
    reduce the active maximum horizon.
4.  Store fitted coefficients and seed IDs, then evaluate classification
    performance using AUC, ROC curves, and optimal ROC thresholds.
5.  Generate AUC, coefficient-density, ROC, and threshold figures, and
    save the processed data in a settings-specific output folder.

### Current example settings

``` matlab
noiseType = "White";
InitialScale = 0;
sigma = 1;
mu_values = -logspace(-1,2,50);
n_realizations = 30;
tmax_values = [1,2,5,10,20,50,100,200,500,1000,2000,5000];
dt_target = 0.01;
s_values = [1,-1];
kc_shat = 2;
nThetaComp = 2;
```

These settings are defined near the beginning of the combined script and
can be changed for different experiments.

## Swift--Hohenberg workflow

The combined script is `Swift_Hohenberg/Combine/test_code_chatgpt.m`. It
applies a related workflow to a spatial stochastic Swift--Hohenberg
system on a periodic domain.

### Main functions

**`simulate_scaled_SH_Fourier`** generates the spatiotemporal field
using an FFT-based solver:

``` matlab
[tvec_full, xvec_full, u_full] = simulate_scaled_SH_Fourier( ...
    mu, s_true, gamma_value, tmax, tsteps, L, xsteps, ...
    noiseType, correlation_length, NoiseScale, InitialScale)
```

**`analyze_SH_scaled`** estimates model coefficients and dominant-mode
quantities:

``` matlab
[theta_sub, kmax_sub, lambda_sub] = ...
    analyze_SH_scaled(u_sub, t_sub, xvec_full);
```

The combined script stores four fitted coefficients. The fourth
coefficient, (`\hat{s}`{=tex}), is used as the classification score
(`kc_shat = 4`). The analysis also returns the estimated wavenumber of
the most unstable mode (`kmax_sub`) and its growth rate (`lambda_sub`).

### Workflow and outputs

1.  Define the (`\mu`{=tex}) sweep, domain length, spatial and temporal
    resolution, cubic-sign classes, noise settings, and observation
    horizons.
2.  Simulate independent seeded spatiotemporal fields with the periodic
    FFT solver, rejecting failed, non-finite, or excessively large
    realisations.
3.  Fit the coefficients and estimate dominant wavenumbers and growth
    rates at the requested observation horizons. If too few realisations
    remain valid at the longest horizon, adaptively retry with shorter
    horizons.
4.  Store the fitted coefficients, seed IDs, wavenumbers, and growth
    rates. Evaluate classification performance using the fitted cubic
    coefficient and compare estimated and theoretical dispersion
    relations.
5.  Generate AUC, coefficient-density, ROC, threshold, and
    dispersion-summary figures, and save the processed data in a
    settings-specific output folder.

The periodic spatial grid uses unique points on (\[0,L)); the endpoint
(x=L) is not duplicated. The domain length controls the available
discrete Fourier modes.

### Current example settings

``` matlab
noiseType = "White";
correlation_length = 0;
InitialScale = 0;
NoiseScale = 1;
mu_values = -logspace(-1,2,20);
n_realizations = 100;
L_values = 20*pi;
dx_target = 0.2*pi;
dt_target = 2.5e-4;
tmax_values = [0.1,0.2,0.5,1,2,5,10,20];
gamma_value = 1;
s_values = [1,-1];
kc_shat = 4;
nThetaComp = 4;
```

These settings are defined near the beginning of the combined script and
can be changed for different experiments.

## Running the workflows

Run MATLAB from the repository root, then add the repository folders to
the MATLAB path:

``` matlab
addpath(genpath(pwd));
```

Run either combined script:

``` matlab
run(fullfile('Pitchfork','Combine','test_code_pitchfork_loops.m'));
run(fullfile('Swift_Hohenberg','Combine','test_code_chatgpt.m'));
```

Run the scripts separately. Each script defines its own settings,
creates a settings-specific output directory, and either loads
compatible saved results or performs the simulation and analysis.

## Saved data and figures

Both combined scripts save processed data in a `saved_data` subfolder
and place generated figures in subfolders under a settings-specific
output directory. The directory names are constructed from key numerical
settings; consult the output-path section of each script for the exact
names.

The saved MAT-files include fitted coefficient estimates (`per_theta`)
and their associated seed IDs (`per_theta_seeds`). The Swift--Hohenberg
results additionally include estimated dominant wavenumbers (`per_kmax`)
and growth rates (`per_lambda_max`). Parameter values and other settings
needed to identify the experiment are saved alongside the results.

AUC values are stored in `auc_all_all`. Other ROC quantities and figures
are calculated from the processed coefficient estimates when the script
runs. If a compatible processed MAT-file is found, the script can skip
the computationally expensive simulation and fitting stage and continue
with analysis and plotting.

## Numerical and reproducibility notes

-   The two classes are defined by the true cubic-coefficient signs
    (s=+1) and (s=-1); the fitted cubic coefficient is the
    classification score.
-   ROC/AUC analysis is carried out for parameter and
    observation-horizon combinations that satisfy the validity and
    sample-count requirements in the script.
-   Observation horizons are analysed using prefixes of trajectories or
    fields simulated up to the active maximum horizon.
-   Random seed IDs are stored with fitted estimates to make
    realisations traceable and support reproducibility.
-   For reproducible comparisons, record the MATLAB release, repository
    commit, all parameter settings, random-number generator and seed
    convention, and parallel execution settings.

## Citation

If you use this code in academic work, cite the associated paper once
its bibliographic details are final. Add the full reference and DOI here
when available.

## Contact

For questions, bug reports, or collaboration, open a GitHub Issue or
contact Paul Sanders at `p.a.sanders@uu.nl`.
