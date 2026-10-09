# Stability Analysis of Pitchfork and Swift--Hohenberg Models

MATLAB code for simulating stochastic pitchfork and Swift--Hohenberg
systems, estimating model coefficients from simulated data, and testing
whether the estimated cubic coefficient distinguishes the two dynamical
regimes.

The repository separates each model into simulation, analysis, and
combined-workflow code:

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

The project investigates whether coefficients estimated from noisy data
can distinguish the signs of the cubic nonlinearity. The pitchfork model
provides a spatially homogeneous test case; the Swift--Hohenberg model
extends the workflow to spatially distributed data.

The estimated cubic coefficient, (\\hat{s}), is used as a classification
score. ROC curves and the area under the ROC curve (AUC) quantify how
well estimates separate simulations with the two true coefficient signs.
Performance is compared across parameter values and observation
horizons.

## Repository structure

### `Pitchfork/`

-   `Simulation/` --- stochastic pitchfork trajectory simulator.
-   `Analysis/` --- coefficient estimation from trajectory data.
-   `Combine/` --- scripts that run batches of simulations, analyse and
    cache results, and generate figures.

### `Swift_Hohenberg/`

-   `Simulation/` --- stochastic Swift--Hohenberg solver on a periodic
    domain using Fourier/FFT methods.
-   `Analysis/` --- coefficient estimation and analysis of dominant
    wavenumbers, growth rates, and dispersion relations.
-   `Combine/` --- scripts that combine simulation, analysis, result
    caching, and plotting.

The folders separate reusable model functions from scripts that
orchestrate full experiments. Check the relevant folder for the current
MATLAB filenames and settings.

## Requirements

The code is written for MATLAB and uses scripts, functions, cell arrays,
structures, plotting, and (in the combined workflows) `parfor`.

-   **Parallel Computing Toolbox** is required for `parfor`. For
    debugging, `parfor` can be changed to `for`, with longer runtimes
    expected.
-   **Statistics and Machine Learning Toolbox** is required for
    `perfcurve` and for any scripts that use `ksdensity` or `prctile`.
-   MATLAB graphics functionality is used for `exportgraphics`,
    `savefig`, tiled layouts, and figure annotations.

No external dataset is required for the standard workflows. When
compatible cached results are not available, the combined scripts
generate the stochastic data.

## Pitchfork workflow

### `simulate_pitchfork`

Generates one stochastic pitchfork trajectory:

``` matlab
[t_vec, u_vec] = simulate_pitchfork( ...
    mu, s_true, tmax, tsteps, noiseType, sigma, InitialScale)
```

  -----------------------------------------------------------------------
  Argument                            Meaning
  ----------------------------------- -----------------------------------
  `mu`                                Linear coefficient. The example
                                      sweep uses `-logspace(-1,2,50)`.

  `s_true`                            True cubic-coefficient sign
                                      defining the class; the example
                                      uses `s_values = [1, -1]`.

  `tmax`                              Maximum simulated time.

  `tsteps`                            Number of time samples, determined
                                      from `tmax` and the target time
                                      step.

  `noiseType`                         Noise specification; the example
                                      uses `"White"`.

  `sigma`                             Noise amplitude.

  `InitialScale`                      Initial-condition scaling
                                      parameter; the example uses `0`.
  -----------------------------------------------------------------------

The function returns the time vector `t_vec` and simulated state data
`u_vec`. The combined script checks the data dimensions and converts the
trajectory to the orientation expected by the analysis code.

### `analyze_pitchfork`

Estimates model coefficients from a trajectory or trajectory prefix:

``` matlab
[theta_sub, output2, output3] = analyze_pitchfork(u_sub, t_sub);
```

The combined workflow stores the first two components of `theta_sub` and
uses the second as its classification score:

``` matlab
kc_shat = 2;
```

### Combined workflow

The pitchfork combined workflow is represented by
`test_code_pitchfork_loops.m` in `Pitchfork/Combine/`.

It defines the parameter sweep and numerical settings, simulates seeded
realisations (using `parfor` where configured), rejects invalid
trajectories, estimates coefficients at multiple observation horizons,
and adaptively falls back to shorter horizons when too few realisations
remain valid at the longest horizon. It then computes AUC and optimal
ROC operating points, generates AUC, coefficient-distribution, ROC, and
threshold figures, and saves fitted coefficients with their metadata.

Example settings from the current workflow:

``` matlab
noiseType = "White";
InitialScale = 0;
sigma = 1;
mu_values = -logspace(-1,2,50);
n_realizations = 20;
tmax_values = [1, 2, 5];
dt_target = 0.01;
s_values = [1, -1];
kc_shat = 2;
```

These are script-specific example settings, not universal defaults.

## Swift--Hohenberg workflow

### `simulate_scaled_SH_Fourier`

Generates stochastic Swift--Hohenberg data on a periodic spatial domain
using an FFT-based solver:

``` matlab
[tvec_full, xvec_full, u_full] = simulate_scaled_SH_Fourier( ...
    mu, s_true, gamma_value, tmax, tsteps, L, xsteps, ...
    noiseType, correlation_length, NoiseScale, InitialScale)
```

  -----------------------------------------------------------------------
  Argument                            Meaning
  ----------------------------------- -----------------------------------
  `mu`                                Linear coefficient; the example
                                      sweep uses `-logspace(-1,2,20)`.

  `s_true`                            True cubic-coefficient sign
                                      defining the class; the example
                                      uses `s_values = [1, -1]`.

  `gamma_value`                       Swift--Hohenberg parameter; the
                                      example uses `1`.

  `tmax`                              Maximum simulated time.

  `tsteps`                            Number of time samples, determined
                                      from `tmax` and the target time
                                      step.

  `L`                                 Periodic-domain length; the example
                                      uses `20*pi`.

  `xsteps`                            Number of unique spatial grid
                                      points, determined by the domain
                                      length and target spacing.

  `noiseType`                         Noise specification; the example
                                      uses `"White"`.

  `correlation_length`                Spatial noise-correlation length;
                                      the example uses `0`.

  `NoiseScale`                        Noise-strength scaling parameter;
                                      the example uses `1`.

  `InitialScale`                      Initial-condition scaling
                                      parameter; the example uses `0`.
  -----------------------------------------------------------------------

The outputs are the time vector `tvec_full`, spatial grid `xvec_full`,
and spatiotemporal field `u_full`. The periodic grid uses unique points
on `[0,L)`, without duplicating the endpoint `x = L`.

### `analyze_SH_scaled`

Estimates model coefficients and dispersion-related quantities from a
spatiotemporal field:

``` matlab
[theta_sub, kmax_sub, lambda_sub] = ...
    analyze_SH_scaled(u_sub, t_sub, xvec_full);
```

  -----------------------------------------------------------------------
  Output                              Meaning
  ----------------------------------- -----------------------------------
  `theta_sub`                         Estimated model coefficients; the
                                      first four components are stored.

  `kmax_sub`                          Estimated wavenumber associated
                                      with the most unstable mode.

  `lambda_sub`                        Estimated growth rate at that
                                      wavenumber.
  -----------------------------------------------------------------------

The fourth coefficient estimate is used as the classification score:

``` matlab
kc_shat = 4;
nThetaComp = 4;
```

### Combined workflow

The Swift--Hohenberg combined workflow is represented by
`test_code_chatgpt.m` in `Swift_Hohenberg/Combine/` in the current
script naming.

It defines the parameter sweep and spatial grid, simulates seeded
realisations, analyses trajectory prefixes, rejects invalid fields, and
stores coefficient estimates, dominant wavenumbers, growth rates, and
seed IDs. It can use shorter observation horizons when too few
realisations survive at the longest horizon. It then evaluates
classification performance, generates AUC, coefficient-distribution,
ROC, threshold, and dispersion-summary figures, and saves processed
results.

Example settings from the current workflow:

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
tmax_values = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20];
gamma_value = 1;
s_values = [1, -1];
kc_shat = 4;
nThetaComp = 4;
```

If these values differ from the current script, use the script itself as
the source of truth.

## Running the workflows

### 1. Clone the repository

``` bash
git clone https://github.com/<your-username>/<your-repository>.git
cd <your-repository>
```

Replace the placeholder URL with the actual repository URL.

### 2. Add the repository to the MATLAB path

From the repository root:

``` matlab
repoRoot = pwd;
addpath(genpath(repoRoot));
```

### 3. Run a combined script

Pitchfork:

``` matlab
run(fullfile('Pitchfork', 'Combine', 'test_code_pitchfork_loops.m'));
```

Swift--Hohenberg:

``` matlab
run(fullfile('Swift_Hohenberg', 'Combine', 'test_code_chatgpt.m'));
```

Adjust the script names if they have changed in your checkout. Before a
full run, inspect the settings near the start of the script and check
that the required functions resolve to the intended files.

### 4. Test with a smaller run first

For debugging, reduce the number of parameter values and realisations
and use shorter observation horizons. Keep the time step and spatial
grid numerically appropriate; increasing the time step simply to speed
up a run can change the results or destabilise the simulation.

## Saved results and figures

The combined scripts save processed results and figures in
settings-specific output folders. Exact folder names and figure
subdirectories depend on the script version and settings; consult the
output-path definitions in the relevant combined script.

MAT-files may include:

-   `per_theta` --- fitted coefficient estimates indexed by model
    component, parameter value, class, and observation horizon;
    Swift--Hohenberg results may also be indexed by domain length.
-   `per_theta_seeds` --- seed IDs associated with the estimates.
-   `auc_all_all` --- AUC values for the selected classification score.
-   `per_kmax` --- estimated dominant wavenumbers for Swift--Hohenberg.
-   `per_lambda_max` --- estimated growth rates at the dominant
    wavenumbers.
-   Parameter and numerical settings used for the run.

Derived ROC curves, optimal operating points, and thresholds may be
recomputed from cached coefficient estimates. They do not necessarily
need to be saved separately. When a compatible cache is found, the
simulation stage can be skipped and downstream analysis and figures
regenerated.

## Numerical and methodological notes

-   Seed IDs are assigned to realisations so simulations can be tracked
    and, subject to MATLAB and solver details, reproduced.
-   Independent simulations may run in parallel using `parfor`.
-   A realisation can be rejected if simulation or analysis fails, data
    are empty or have unexpected dimensions, values are non-finite, or
    the configured magnitude threshold is exceeded.
-   Observation horizons are evaluated using prefixes of a trajectory or
    field generated up to the active maximum horizon.
-   The classification treats `s > 0` as the positive class and `s < 0`
    as the negative class.
-   ROC/AUC analysis is restricted to parameter/horizon combinations
    that meet the script's validity and sample-count criteria for both
    classes.
-   Swift--Hohenberg simulations use a periodic spatial grid with unique
    points on `[0,L)`. Dispersion summaries compare estimated
    growth-rate information with the theoretical relation implemented in
    the analysis script.

## Reproducibility

For each reported run, record the MATLAB version, operating system,
repository commit, parameter and numerical settings, random-number
generator and seed convention, parallelisation settings, and versions of
the simulation and analysis functions.

Parallel execution, MATLAB versions, numerical libraries, and hardware
can affect execution order and floating-point results. Preserve seed IDs
and settings alongside processed data when comparing runs.

## Troubleshooting

### MATLAB cannot find a function

Check that the required folders are on the path and that each function
resolves to the intended file:

``` matlab
which simulate_pitchfork -all
which analyze_pitchfork -all
which simulate_scaled_SH_Fourier -all
which analyze_SH_scaled -all
```

### `parfor` is unavailable

Check that the Parallel Computing Toolbox is installed and licensed. For
debugging, replace `parfor` with `for` in the relevant loop; runtime may
increase substantially.

### Too few realisations are accepted

Check the time step, noise settings, configured magnitude threshold,
observation horizons, simulator output dimensions, and data format
expected by the analysis function. The scripts may fall back to shorter
horizons when too few realisations survive at the longest horizon.

### A cached file is loaded unexpectedly

Inspect the cache path and the settings used to construct it. To force a
fresh simulation, move or remove the relevant cached result or change
the output location. Preserve results needed for comparison before
deleting anything.

## Citation

If you use this code in academic work, cite the associated paper once
the bibliographic details are final. Add the full reference and DOI here
when available.

## Contact

For questions, bug reports, or collaboration, open a GitHub Issue or
contact Paul Sanders at `p.a.sanders@uu.nl`.
