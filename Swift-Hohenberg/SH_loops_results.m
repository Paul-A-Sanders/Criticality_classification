% SH_loops_results.m
% =========================================================================
% PURPOSE
%   Generate stochastic Swift-Hohenberg data on a periodic domain with an
%   FFT-based solver, estimate the model's linear and cubic coefficients
%   together with the wavenumber and growth rate of the most unstable
%   mode (analyze_SH_scaled), and test how well the fitted cubic
%   coefficient s_hat separates the two signs s = +1 and s = -1:
%
%       dq = [mu*q - gamma*(1+d_xx)^2*q + s*q^3] dtau ,
%
%   with mu < 0 and a fixed gamma (see gamma_value).
%
% WORKFLOW
%   1. If a compatible processed MAT-file exists, load it and skip the
%      simulations.
%   2. Otherwise, for every (mu, s) sweep: simulate independent noise
%      realisations in parallel (parfor), reject diverging or failed runs
%      and fit each usable run for every requested observation horizon
%      tau_max. If too few realisations reach the longest horizon, retry at
%      shorter horizons.
%   3. Compute AUC of s_hat for every (mu, tau_max) and save the processed
%      estimates.
%   4. Draw figures: AUC vs mu, distributions of the estimates, ROC curves,
%      and estimated vs theoretical dispersion relations.
%
% REQUIREMENTS
%   * simulate_scaled_SH_Fourier.m and analyze_SH_scaled.m on the MATLAB path.
%   * Parallel Computing Toolbox (parfor) and Statistics and Machine
%     Learning Toolbox (perfcurve, ksdensity, prctile).
%   * MATLAB R2020a or newer (exportgraphics, tiledlayout); developed with
%     R2026b.
%
% OUTPUTS
%   A folder named after the key settings is created in the current
%   directory:
%       SH_ROC_scaled_FFT_mu_s_real=..._dx=...pi_dt=.../
%           saved_data/                processed MAT-file (see DATA LAYOUT)
%           auc_plots/                 AUC vs mu, one line per tau_max
%           ridge_plots_all/           distributions of the fitted coefficients
%           roc_plots_per_time/        ROC curves, one figure per tau_max
%           dispersion_summary_panels/ estimated vs theoretical dispersion
%
% DATA LAYOUT (saved in the MAT-file)
%   per_theta{comp, iMu, iS, iL, it}   row vector of fitted values of
%                                      theta = [a0; a2; a4; s]_hat for
%                                      mu_values(iMu), s_values(iS),
%                                      L_values(iL), tmax_values(it)
%   per_theta_seeds{...}               RNG seeds of those realisations
%   per_kmax{iMu, iS, iL, it}          fitted most-unstable wavenumber
%   per_lambda_max{iMu, iS, iL, it}    fitted largest growth rate
%   auc_all_all(iMu, it, iL)           AUC of s_hat
%
% REUSING SAVED DATA
%   If the MAT-file for these settings already exists, it is loaded and the
%   simulations are skipped. The values stored in that file (mu_values,
%   tmax_values, L_values, s_values, gamma_value, n_realizations, ...) then
%   REPLACE the ones set below. To re-simulate after changing a setting that
%   is not part of the folder name (for example mu_values), set
%   rerunSimulations = true or delete the saved file.
%
% Review the user settings below before running.
% =========================================================================

clearvars; clc;
% set(0,'DefaultFigureVisible','off');   % Uncomment to run without showing figures.

numCores = feature('numCores');   % Used to size the parallel batches.

%% ========================================================================
%% 1. USER SETTINGS
%% ========================================================================

% --- Noise ---
noiseType          = "White";   % 'White', 'Correlated', 'Uniform' or 'None'
correlation_length = 0;         % Smoothing length; only used for 'Correlated'
InitialScale       = 0;         % Std. dev. of the random initial condition; 0 starts at q = 0
NoiseScale         = 1;         % Noise amplitude

% --- Parameter sweep ---
mu_values = -logspace(-1,2,20);   % mu from -0.1 to -100 (log-spaced). mu < 0 is below the
                                  % linear instability; values near 0 are close to onset.
gamma_value = 1;                  % Fixed gamma in this experiment
s_values  = [1, -1];              % Sign of the cubic term: the two classes to distinguish
nS        = numel(s_values);

% --- Ensemble and observation horizons ---
n_realizations        = 100;      % Target number of realisations per (mu, s, tau_max)
u_explosion_threshold = 1e10;     % Reject trajectories whose |q| exceeds this value
tmax_values = [0.1,0.2,0.5,1,2,5,10,20];   % Observation horizons tau_max
nT = numel(tmax_values);

% --- Grid ---
L_values  = 20*pi;       % Periodic domain length(s)
dx_target = 0.2*pi;      % Grid spacing: Nx = L/dx_target = 100 points
dt_target = 2.5e-4;      % Explicit-Euler time step; must be small enough for stability
                         % (simulate_scaled_SH_Fourier checks dt*max|lambda_k| < 2).
                         % tmax and L must be integer multiples of dt_target and dx_target.

% --- Classification score ---
kc_shat    = 4;   % theta = [a0; a2; a4; s]_hat; component 4 (s_hat) is the ROC score
nThetaComp = 4;   % Number of fitted components that are kept

if kc_shat > nThetaComp
    error('kc_shat out of range');
end

% --- Run control ---
rerunSimulations = false;   % true: ignore any saved MAT-file and simulate again

nMu = numel(mu_values);
nL  = numel(L_values);

%% ========================================================================
%% 2. OUTPUT LOCATIONS
%% ========================================================================
% The key settings are part of the folder name, so different settings never
% overwrite each other.
base_dir_root = fullfile(pwd, sprintf( ...
    'SH_ROC_scaled_FFT_mu_s_real=%d_dx=%.3gpi_dt=%.3g', ...
    n_realizations, dx_target/pi, dt_target));
if ~exist(base_dir_root,'dir'), mkdir(base_dir_root); end

outDataDir = fullfile(base_dir_root, 'saved_data');
saveName   = sprintf('SH_ROC_data_FFT_mu_s_real%d_dx%.3gpi_dt%.3g.mat', ...
    n_realizations, dx_target/pi, dt_target);
savePath   = fullfile(outDataDir, saveName);

%% ========================================================================
%% 3. LOAD PROCESSED DATA (IF AVAILABLE)
%% ========================================================================
% A saved file lets the script skip the expensive simulation and fitting.
loadedFromFile = false;
if exist(savePath, 'file') && ~rerunSimulations
    fprintf('Found saved data at %s. Loading and skipping simulations...\n', savePath);
    S = load(savePath);

    requiredFields = {'per_theta','per_theta_seeds','per_kmax','per_lambda_max', ...
        'mu_values','tmax_values','L_values','s_values','kc_shat','nThetaComp','n_realizations'};
    missingFields = requiredFields(~isfield(S, requiredFields));
    if ~isempty(missingFields)
        error('Saved file is missing required fields: %s', strjoin(missingFields, ', '));
    end

    per_theta       = S.per_theta;
    per_theta_seeds = S.per_theta_seeds;
    per_kmax        = S.per_kmax;
    per_lambda_max  = S.per_lambda_max;
    mu_values       = S.mu_values;
    tmax_values     = S.tmax_values;
    L_values        = S.L_values;
    s_values        = S.s_values;
    kc_shat         = S.kc_shat;
    nThetaComp      = S.nThetaComp;
    n_realizations  = S.n_realizations;
    if isfield(S,'gamma_value'), gamma_value = S.gamma_value; end

    nMu = numel(mu_values);
    nT  = numel(tmax_values);
    nL  = numel(L_values);
    nS  = numel(s_values);
    loadedFromFile = true;
end

%% ========================================================================
%% 4. SIMULATION AND PARAMETER ESTIMATION
%% ========================================================================
% Skipped when processed data were loaded above.
if ~loadedFromFile

    per_theta       = cell(nThetaComp, nMu, nS, nL, nT);   % See DATA LAYOUT in the header.
    per_theta_seeds = cell(nThetaComp, nMu, nS, nL, nT);
    per_kmax        = cell(nMu, nS, nL, nT);
    per_lambda_max  = cell(nMu, nS, nL, nT);

    for iMu = 1:nMu
        mu = mu_values(iMu);
        fprintf('mu (%d/%d) = %g\n', iMu, nMu, mu);

        for iL = 1:nL
            % L sets the available spatial wavelengths and therefore the
            % discrete Fourier modes represented by the simulation.
            L = L_values(iL);

            % Periodic FFT grid: Nx unique points on [0,L), so x = L is not duplicated.
            xsteps = round(L/dx_target);
            dx_actual = L/xsteps;
            if abs(dx_actual-dx_target) > 1e-12*max(1,abs(dx_target))
                error('L must be an integer multiple of dx_target.');
            end

            fprintf('L=%.6g, Nx=%d, dx=%.12g\n', L, xsteps, dx_actual);

            for iS = 1:nS
                % Both classes are run at the same mu and L so that their
                % fitted coefficient distributions can be compared directly.
                s_true = s_values(iS);
                fprintf('\ns = %g\n', s_true);

                % Start every (mu, s, L) at the longest requested horizon. If a
                % setting (typically the subcritical case) needs a shorter
                % usable horizon, this does not carry over to the next s.
                % activeIt indexes tmax_values.
                activeIt = nT;
                tmax_max = tmax_values(activeIt);
                [tsteps_max, tidx_ref, dt_actual] = horizonGrid(tmax_values, tmax_max, dt_target);

                fprintf('Starting this setting at tau_max=%.3g (Nt=%d, dt=%.12g)\n', ...
                    tmax_max, tsteps_max, dt_actual);

                batchSize       = n_realizations;       % Upper limit on seeds per batch
                maxBatches      = 20;
                maxSeedAttempts = batchSize*maxBatches; % Give up after this many seeds; the
                                                        % setting then keeps what it has collected.
                seedStart = 0;
                batch     = 0;

                % 'successful' pools the accepted realisations. Each entry has
                %   seed        RNG seed of the realisation
                %   maxIt       index of the longest tau_max with a valid fit
                %   theta       cell(1,nT) of fitted coefficient vectors
                %   kmax        1-by-nT fitted most-unstable wavenumbers
                %   lambdaMax   1-by-nT fitted largest growth rates
                %   validTimes  logical(1,nT) marking horizons with a valid fit
                successful = struct('valid',false,'seed',[],'maxIt',0,'theta',{{}}, ...
                    'kmax',[],'lambdaMax',[],'validTimes',false(1,nT));
                successful = successful([]);   % Start with an empty struct array

                while seedStart < maxSeedAttempts

                    collected = sum([successful.maxIt] >= activeIt);
                    if collected >= n_realizations
                        break
                    end

                    batch = batch + 1;
                    fprintf('\nBatch %d\n', batch);

                    numSeedsThisBatch = nextBatchSize(n_realizations - collected, numCores, batchSize);
                    seeds     = seedStart + (1:numSeedsThisBatch);
                    seedStart = seedStart + numSeedsThisBatch;

                    results = struct('valid',false,'seed',0,'maxIt',0,'theta',{{}}, ...
                        'kmax',[],'lambdaMax',[],'validTimes',false(1,nT));
                    results = repmat(results, numel(seeds), 1);
                    fprintf('Trying seeds %d -- %d\n', seeds(1), seeds(end));

                    % ---- Parallel simulation and fitting, one seed per iteration ----
                    parfor ib = 1:numel(seeds)

                        seed = seeds(ib);
                        rng(seed,'combRecursive');   % Same seed = same independent realisation

                        theta_store  = cell(1,nT);
                        kmax_store   = nan(1,nT);
                        lambda_store = nan(1,nT);
                        validTimes   = false(1,nT);
                        local_maxIt  = 0;

                        try
                            % One stochastic space-time realisation up to the active
                            % horizon; it is truncated below for every tau_max.
                            [tvec_full,xvec_full,u_full] = simulate_scaled_SH_Fourier( ...
                                mu, s_true, gamma_value, tmax_max, tsteps_max, L, xsteps, ...
                                noiseType, correlation_length, NoiseScale, InitialScale);
                        catch
                            % Any error skips this seed silently. This includes the
                            % explicit-Euler stability check in simulate_scaled_SH_Fourier:
                            % if dt_target is too large for the largest |mu|, EVERY seed
                            % fails and the horizon-reduction logic below ends in
                            % "Cannot reduce tmax further". Test one call without try/catch
                            % to see the real error message.
                            continue
                        end

                        if isempty(u_full)
                            continue
                        end

                        % Make sure the time dimension runs along the rows.
                        [n1,n2] = size(u_full);
                        if n1 == numel(tvec_full)
                            u_sim = u_full;
                        elseif n2 == numel(tvec_full)
                            u_sim = u_full';
                        else
                            continue
                        end

                        % Fit the first tau_max time units for every horizon.
                        for it = 1:nT

                            idx_end = tidx_ref(it);
                            if isnan(idx_end) || idx_end < 3
                                continue
                            end

                            t_sub = tvec_full(1:idx_end);
                            u_sub = u_sim(1:idx_end,:);

                            % Reject non-finite or exploding trajectory segments.
                            if any(isnan(u_sub(:))) || any(isinf(u_sub(:)))
                                continue
                            end
                            if max(abs(u_sub(:))) > u_explosion_threshold
                                continue
                            end

                            try
                                % theta = [a0; a2; a4; s]_hat; kmax and lambdaMax describe the
                                % most unstable mode of the fitted linear operator.
                                [theta_sub, kmax_sub, lambda_sub] = analyze_SH_scaled(u_sub, t_sub, xvec_full);
                            catch
                                continue
                            end

                            theta_sub = reshape(theta_sub(:),[],1);
                            if numel(theta_sub) < nThetaComp
                                theta_sub(end+1:nThetaComp,1) = NaN;
                            end

                            theta_store{it} = theta_sub(1:nThetaComp);
                            kmax_store(it)  = kmax_sub;
                            lambda_store(it) = lambda_sub;
                            validTimes(it)  = true;
                            local_maxIt     = max(local_maxIt, it);

                        end

                        % Store the result for this seed.
                        if any(validTimes)
                            results(ib).valid      = true;
                            results(ib).seed       = seed;
                            results(ib).maxIt      = local_maxIt;
                            results(ib).theta      = theta_store;
                            results(ib).kmax       = kmax_store;
                            results(ib).lambdaMax  = lambda_store;
                            results(ib).validTimes = validTimes;
                        end

                    end % parfor

                    % ---- Merge new valid results into the pool (skip duplicate seeds) ----
                    % This is an extra protection line for if the seeds are
                    % chosen differently so that duplicates exist, currently it
                    % is not relevant
                    existingSeeds = [successful.seed];

                    for ib = 1:numel(results)
                        if ~results(ib).valid
                            continue
                        end
                        if ismember(results(ib).seed, existingSeeds)
                            continue
                        end
                        successful(end+1)    = results(ib); %#ok<SAGROW>
                        existingSeeds(end+1) = results(ib).seed; %#ok<SAGROW>
                    end

                    % ---- Progress under the current horizon ----
                    collected = sum([successful.maxIt] >= activeIt);
                    fprintf('Batch %d done | collected %d / %d (activeIt=%d)\n', ...
                        batch, collected, n_realizations, activeIt);

                    % ---- Adaptive reduction of the observation horizon ----
                    % If fewer than 10% of the requested realisations reach the
                    % active horizon, shorten it by one step.
                    minAccept = ceil(0.1 * n_realizations);

                    if collected < minAccept
                        if activeIt > 1
                            activeIt = activeIt - 1;

                            beforeN    = numel(successful);
                            successful = successful([successful.maxIt] >= activeIt);
                            fprintf('Reducing active tmax %.3g -> %.3g (kept seeds: %d -> %d)\n', ...
                                tmax_max, tmax_values(activeIt), beforeN, numel(successful));

                            tmax_max = tmax_values(activeIt);
                            [tsteps_max, tidx_ref] = horizonGrid(tmax_values, tmax_max, dt_target);
                        else
                            warning('Cannot reduce tmax further. Stopping.');
                            break
                        end
                    end

                end % while

                % ---- Store the accepted fits directly in the output arrays ----
                % 'successful' is only needed while the horizon is being
                % determined. Once the final active horizon is known, keep the
                % accepted seeds, copy their fits into the saved arrays and
                % discard the temporary structure.
                if ~isempty(successful)
                    successful = successful([successful.maxIt] >= activeIt);

                    if numel(successful) > n_realizations
                        successful = successful(1:n_realizations);
                    end

                    for iRun = 1:numel(successful)
                        runResult = successful(iRun);

                        for it = 1:nT
                            if isempty(runResult.theta) || numel(runResult.theta) < it
                                continue
                            end

                            theta_it = runResult.theta{it};
                            if isempty(theta_it)
                                continue
                            end

                            for kc = 1:nThetaComp
                                if kc <= numel(theta_it) && isfinite(theta_it(kc))
                                    per_theta{kc,iMu,iS,iL,it}(end+1)       = theta_it(kc);
                                    per_theta_seeds{kc,iMu,iS,iL,it}(end+1) = runResult.seed;
                                end
                            end

                            if numel(runResult.kmax) >= it && isfinite(runResult.kmax(it))
                                per_kmax{iMu,iS,iL,it}(end+1) = runResult.kmax(it);
                            end

                            if numel(runResult.lambdaMax) >= it && isfinite(runResult.lambdaMax(it))
                                per_lambda_max{iMu,iS,iL,it}(end+1) = runResult.lambdaMax(it);
                            end
                        end
                    end
                end

            end
        end
    end

else
    fprintf('Proceeding directly to analysis and plotting using loaded fitted data.\n');
end

%% ========================================================================
%% 5. CLASSIFICATION PERFORMANCE (AUC)
%% ========================================================================
% The fitted cubic coefficient is the score for telling the two signs of the
% nonlinear term apart; s > 0 is the positive ROC class. AUC is computed
% separately for every mu, tau_max and domain length so that conditions can
% be compared. AUC = 1: perfect separation; AUC = 0.5: chance level.
%
% NOTE: every available fit is used, however few there are (>= 1 per class).
% This differs from Pitchfork_loops_results.m and the gamma-sweep script,
% which require n_realizations fits in BOTH classes before an AUC is
% computed. For horizons longer than the final active horizon of a (mu, s)
% setting only the few early realisations are available, so those AUC values
% rest on small samples.
fprintf('\nComputing AUC...\n');

auc_all_all = nan(nMu, nT, nL);

for iMu = 1:nMu
    for iL = 1:nL
        for it = 1:nT

            scores_all = [];
            labels_all = [];

            for iS = 1:nS
                vals = per_theta{kc_shat, iMu, iS, iL, it};

                if isempty(vals)
                    continue
                end

                vals = vals(:)';

                scores_all = [scores_all, vals]; %#ok<AGROW>
                labels_all = [labels_all, repmat(s_values(iS) > 0, 1, numel(vals))]; %#ok<AGROW>
            end

            % Both classes are needed.
            if numel(unique(labels_all)) < 2 || isempty(scores_all)
                continue
            end

            [~,~,~,auc_all_all(iMu,it,iL)] = perfcurve(labels_all, scores_all, 1);

        end
    end
end

%% ========================================================================
%% 6. SAVE PROCESSED RESULTS
%% ========================================================================
% Fitted coefficients and derived quantities are saved rather than the full
% stochastic trajectories. This keeps the file compact while preserving
% everything needed to reproduce the analysis and figures.
if ~loadedFromFile
    if ~exist(outDataDir,'dir'), mkdir(outDataDir); end

    % MAT-file v7 is used on purpose: the files are only a few MB, so
    % v7.3/HDF5 is unnecessary and can be problematic when writing directly
    % into a synchronised cloud folder (e.g. OneDrive).
    save(savePath, ...
        'per_theta', ...
        'per_theta_seeds', ...
        'per_kmax', ...
        'per_lambda_max', ...
        'mu_values', ...
        'tmax_values', ...
        'L_values', ...
        's_values', ...
        'kc_shat', ...
        'nThetaComp', ...
        'n_realizations', ...
        'gamma_value', ...
        'dx_target', ...
        'dt_target', ...
        'auc_all_all', ...
        'base_dir_root', ...
        '-v7');
    fprintf('Saved results to: %s\n', savePath);
else
    fprintf('Loaded data present: skipping save to avoid overwriting existing file.\n');
end

%% ========================================================================
%% 7. FIGURE STYLE
%% ========================================================================
% Everything below is presentation code: the numbers have already been
% computed above.
paperFont = 'Helvetica';

fsTick        = 13;    % Tick labels
fsLabel       = 15;    % Axis labels
fsPanelTitle  = 14;    % Panel titles and row labels
fsFigureTitle = 17;    % Figure titles
fsLegend      = 12;
fsTextBox     = 10;    % Statistics boxes inside panels

axisLineWidth = 1.1;

% Apply the common axes style; extra name/value pairs are passed to set().
styleAxes = @(ax,varargin) set(ax,'FontName',paperFont,'FontSize',fsTick, ...
    'LineWidth',axisLineWidth,'TickDir','in',varargin{:});

%% ========================================================================
%% 8. FIGURE: AUC vs mu (one line per tau_max)
%% ========================================================================
% Uses the first domain length (iL = 1); change the index below to show another.
fig = figure('Name','AUC vs mu - all tmax (fixed L)','Color','w','Units','normalized', ...
    'Position',[0.15 0.15 0.7 0.6]);
ax = axes(fig); hold(ax,'on');
styleAxes(ax,'Box','on');

colorsT = lines(nT);    % One colour per tau_max

% MATLAB log axes need positive coordinates, so |mu| is plotted and the
% ticks are labelled with the negative physical values (all mu < 0 here).
xpos   = -mu_values;
useLog = (max(xpos)/min(xpos) > 10);
if useLog
    xplot = xpos;
else
    xplot = mu_values;
end

hPlot  = gobjects(nT,1);
labels = cell(nT,1);

% One colour-coded line per tau_max, plus dots at the simulated mu values.
for it = 1:nT
    auc_vec = auc_all_all(:,it,1);
    if all(isnan(auc_vec)), continue; end

    color = colorsT(it,:);
    valid = isfinite(auc_vec) & isfinite(xplot(:));
    if ~any(valid), continue; end

    if useLog
        hPlot(it) = semilogx(ax, xplot, auc_vec, 'Color', color, 'LineWidth', 1.8);
    else
        hPlot(it) = plot(ax, xplot, auc_vec, 'Color', color, 'LineWidth', 1.8);
    end

    scatter(ax, xplot(valid), auc_vec(valid), 36, 'Marker', 'o', 'MarkerFaceColor', color, ...
        'MarkerEdgeColor', color*0.6, 'MarkerFaceAlpha', 0.9);

    labels{it} = sprintf('$\\tau_{\\mathrm{max}}=%.3g$', tmax_values(it));
end

% Keep only the lines that were actually drawn for the legend.
validIdx = arrayfun(@isgraphics, hPlot);
hPlot  = hPlot(validIdx);
labels = labels(validIdx);

% Axis: mu from -10^2 on the left to -10^-1 on the right.
if useLog
    setLogMuAxis(ax,xpos);
else
    xlim(ax, [min(mu_values) max(mu_values)]);
end

xlabel(ax, '$\mu$', 'FontSize', fsLabel, 'Interpreter','latex');
ylabel(ax, 'AUC', 'FontSize', fsLabel);
title(ax, sprintf('AUC vs $\\mu$ $(L = %.3g \\pi)$', L_values(1)/pi), ...
    'FontSize', fsFigureTitle, 'Interpreter','latex');
ylim(ax, [0 1]);
set(ax, 'Layer', 'top');

if ~isempty(hPlot)
    lg = legend(ax, hPlot, labels, 'Location','best', 'Interpreter','latex', ...
        'FontSize', fsLegend);
    lg.Box = 'on';
    lg.NumColumns = min(3, numel(hPlot));
end

grid(ax,'off');
hold(ax,'off');
drawnow;

outDirAUC = fullfile(base_dir_root,'auc_plots');
if ~exist(outDirAUC,'dir'), mkdir(outDirAUC); end

exportgraphics(fig, fullfile(outDirAUC, sprintf('AUC_vs_mu_L=%g.jpg', L_values(1))), ...
    'Resolution',300);
savefig(fig, fullfile(outDirAUC, sprintf('AUC_vs_mu_L=%g.fig', L_values(1))));

%% ========================================================================
%% 9. FIGURES: distributions (ridge panels) of the fitted coefficients
%% ========================================================================
% One figure per component, domain length and tau_max. Rows = a few mu
% values spread over the sweep, columns = the two true classes. Each panel
% shows the kernel density estimate (scaled to peak 1) of the fitted
% coefficient, its mean (dashed line) and a box with n, mean and standard
% deviation. They show how much the classes overlap, complementing the AUC.
kc_list = unique([1, kc_shat]);
minSamplesPlot = 3;                     % Minimum samples needed to draw a panel

cmapPos = [0 0.4470 0.7410];            % s > 0 (blue)
cmapNeg = [0.8500 0.3250 0.0980];       % s < 0 (orange)

maxMuPlot   = 5;                        % Show at most 5 mu values (rows)
idx_mu_plot = unique(round(linspace(1, nMu, min(nMu, maxMuPlot))));
nMu_plot    = numel(idx_mu_plot);

outDir = fullfile(base_dir_root, 'ridge_plots_all');
if ~exist(outDir,'dir'), mkdir(outDir); end

for kc = kc_list
    if kc == 1
        % theta(1) estimates mu - gamma, so mu_hat - gamma is compared with mu.
        compLabel = sprintf('$\\hat{\\mu}-%g$', gamma_value);
        plotTag   = sprintf('mu_minus_%g', gamma_value);
    else
        compLabel = '$\hat{s}$';
        plotTag   = 's_hat';
    end

    for iL_plot = 1:nL
        for it_plot = 1:nT

            figRidge = figure('Name', sprintf('Ridge %s L=%g t=%g', plotTag, ...
                L_values(iL_plot)/pi, tmax_values(it_plot)), ...
                'Color','w','Units','normalized','Position',[0.08 0.07 0.82 0.82]);

            t = tiledlayout(figRidge, nMu_plot, nS, ...
                'TileSpacing','compact','Padding','compact');

            for iMu_plot = 1:nMu_plot
                iMu = idx_mu_plot(iMu_plot);
                mu  = mu_values(iMu);

                leftHasSamples = false;   % Reset for each row (see y tick labels below)

                for iS = 1:nS
                    ax = nexttile(t, (iMu_plot-1)*nS + iS);
                    cla(ax); hold(ax,'on');
                    styleAxes(ax,'TickLabelInterpreter','latex');

                    scores = per_theta{kc, iMu, iS, iL_plot, it_plot};
                    scores = scores(isfinite(scores));

                    hasData = numel(scores) >= minSamplesPlot;   % Fewer samples: no panel

                    % Column header: the true s.
                    if iMu_plot == 1
                        title(ax, sprintf('$s=%g$',s_values(iS)), ...
                            'Interpreter','latex','FontSize',fsPanelTitle);
                    end

                    % y tick labels are kept on the left; the right panel hides
                    % them only if the left panel has data.
                    if iS == 1
                        leftHasSamples = hasData;
                    end

                    % Row label: mu, written vertically on the outer right.
                    if iS == nS
                        addRowLabel(ax, mu, fsPanelTitle);
                    end

                    % Column label under the bottom row: the fitted coefficient.
                    % The numerical x tick labels stay in all rows because the
                    % mu_hat panels use row-dependent x limits.
                    if iMu_plot == nMu_plot
                        xlabel(ax, compLabel, ...
                            'Interpreter','latex','FontSize',fsLabel);
                    end

                    if ~hasData
                        showNoSamples(ax,fsTick);
                        continue
                    end

                    scores = scores(:);

                    % Fixed x windows so panels can be compared within a row.
                    if kc == 1
                        theoreticalValue = mu - gamma_value;   % True value of theta(1) = mu - gamma
                        xLimits = theoreticalValue + [-10 10];
                    else
                        xLimits = [-30 30];
                    end
                    xgrid = linspace(xLimits(1), xLimits(2), 400);

                    % Kernel density estimate only (no histograms). If it cannot be
                    % computed, f stays empty and the panel is marked "no samples".
                    try
                        [f, xi] = ksdensity(scores, xgrid, 'Support','unbounded');
                    catch
                        f = [];
                    end

                    if isempty(f) || all(~isfinite(f)) || max(f,[],'omitnan') <= 0
                        % The density could not be evaluated on this window.
                        text(ax,0.5,0.5,'no samples', ...
                            'Units','normalized','HorizontalAlignment','center', ...
                            'VerticalAlignment','middle','FontSize',fsTick);
                        xlim(ax,xLimits); ylim(ax,[0 1]);
                        if iS ~= 1 && leftHasSamples
                            yticklabels(ax,[]);
                        end
                        box(ax,'on');
                        hold(ax,'off');
                        continue
                    end

                    f = f / max(f,[],'omitnan');   % Normalise the peak to 1

                    color = cmapPos;
                    if s_values(iS) < 0, color = cmapNeg; end

                    fill(ax, [xi fliplr(xi)], [f zeros(size(f))], color, ...
                        'FaceAlpha',0.6, 'EdgeColor','none');
                    plot(ax, xi, f, 'k-', 'LineWidth',0.6);

                    % Mean (dashed line) and summary statistics.
                    m     = mean(scores,'omitnan');
                    sd    = std(scores,'omitnan');
                    nSamp = numel(scores);

                    plot(ax, [m m], [0 1], '--', 'Color',0.15*[1 1 1]);

                    ridgeText = sprintf('n=%d\nm=%.3g\nsd=%.3g', nSamp, m, sd);
                    text(ax, 0.02, 0.96, ridgeText, ...
                        'Units','normalized','VerticalAlignment','top', ...
                        'BackgroundColor','w','EdgeColor','k', ...
                        'FontSize',fsTextBox,'Margin',1);

                    xlim(ax, xLimits);
                    ylim(ax, [0 1.05]);

                    % Hide the numerical y tick labels of the right panel (the
                    % rows share their scale); ticks and grid stay.
                    if iS ~= 1 && leftHasSamples
                        yticklabels(ax,[]);
                    end

                    grid(ax,'on'); box(ax,'on');
                    hold(ax,'off');
                end
            end

            % One shared density label for the whole figure.
            ylabel(t, 'Normalized density', 'FontSize',fsLabel);

            titleStr = sprintf('%s distributions $(L=%.3g\\pi,\\;\\tau_{max}=%.3g)$', ...
                compLabel, L_values(iL_plot)/pi, tmax_values(it_plot));
            title(t, titleStr, 'Interpreter','latex','FontSize',fsFigureTitle);

            fname = fullfile(outDir, sprintf('ridge_%s_L=%g_t=%g.jpg', ...
                plotTag, L_values(iL_plot), tmax_values(it_plot)));
            try
                exportgraphics(figRidge, fname, 'Resolution',300);
            catch
                warning('Could not save ridge figure to %s', fname);
            end
        end
    end
end

%% ========================================================================
%% 10. FIGURES: ROC curves per observation horizon
%% ========================================================================
% ROC curves show the full classification trade-off for a representative
% set of mu values (one figure per tau_max, first domain length).
iL_plot = 1;
kc_plot = kc_shat;

outDirROC = fullfile(base_dir_root, 'roc_plots_per_time');
if ~exist(outDirROC,'dir'), mkdir(outDirROC); end

% Representative mu values: the grid points closest (in log10|mu|) to ten
% log-spaced targets between -100 and -0.1.
nROCcurves    = 10;
target_mu_roc = -logspace(2,-1,nROCcurves);

idx_mu_roc = zeros(1,nROCcurves);
for j = 1:nROCcurves
    [~,idx_mu_roc(j)] = min(abs( ...
        log10(-mu_values) - log10(-target_mu_roc(j))));
end
idx_mu_roc = unique(idx_mu_roc,'stable');

for it_plot = 1:nT

    % ---- Pool scores and labels of both classes for every mu ----
    all_scores_by_mu = cell(nMu,1);
    all_labels_by_mu = cell(nMu,1);

    for iMu = 1:nMu
        scores = [];
        labels = [];

        for iS = 1:nS
            vals = per_theta{kc_plot, iMu, iS, iL_plot, it_plot};
            if isempty(vals), continue; end

            vals = vals(:);
            lbl  = repmat(s_values(iS)>0, numel(vals), 1);
            scores = [scores; vals]; %#ok<AGROW>
            labels = [labels; lbl];  %#ok<AGROW>
        end

        all_scores_by_mu{iMu} = scores;
        all_labels_by_mu{iMu} = labels;
    end

    % ---- ROC figure for this tau_max ----
    figROC = figure('Name',sprintf('ROC t=%g', tmax_values(it_plot)), ...
        'Color','w','Units','normalized','Position',[0.18 0.18 0.52 0.58]);
    axROC = axes(figROC); hold(axROC,'on');
    styleAxes(axROC,'TickLabelInterpreter','latex');

    plot(axROC, [0 1], [0 1], 'k--','LineWidth',1,'HandleVisibility','off');   % Chance line

    colors = lines(numel(idx_mu_roc));

    for j = 1:numel(idx_mu_roc)

        iMu = idx_mu_roc(j);

        scores = all_scores_by_mu{iMu};
        labels = all_labels_by_mu{iMu};

        if isempty(scores) || numel(unique(labels)) < 2
            continue
        end

        [fpr,tpr,~,auc] = perfcurve(labels,scores,1);

        plot(axROC,fpr,tpr, ...
            'Color',colors(j,:), ...
            'LineWidth',1.6, ...
            'DisplayName',sprintf( ...
            '\\mu=%.3g (AUC=%.2f)',mu_values(iMu),auc));
    end

    xlabel(axROC,'False Positive Rate','FontSize',fsLabel);
    ylabel(axROC,'True Positive Rate','FontSize',fsLabel);
    title(axROC, sprintf('ROC $(L=%.3g\\pi, \\tau_{max}=%g)$', ...
        L_values(iL_plot)/pi, tmax_values(it_plot)), 'Interpreter','latex', ...
        'FontSize',fsFigureTitle);

    axis(axROC,[0 1 0 1]);
    box(axROC,'on');
    grid(axROC,'on');
    lgROC = legend(axROC,'Location','southeast','Interpreter','latex', ...
        'FontSize',fsLegend+1);
    lgROC.Box = 'on';
    drawnow;

    fn = fullfile(outDirROC, sprintf('ROC_L%g_t%g.jpg', ...
        L_values(iL_plot), tmax_values(it_plot)));
    try
        exportgraphics(figROC, fn, 'Resolution',300);
    catch
        warning('Could not save %s', fn);
    end
end

%% ========================================================================
%% 11. FIGURES: estimated vs theoretical dispersion relation
%% ========================================================================
% One figure per domain length and tau_max. Rows = a few mu values,
% columns = the two true classes. Each panel shows
%   * the theoretical dispersion lambda(k) = mu - gamma*(1-k^2)^2 (black),
%     with its true maximum (k*, lambda*) = (1, mu) as a yellow star;
%   * the dispersion curves lambda(k) = a0 - a2*k^2 + a4*k^4 built from the
%     fitted coefficients of every realisation: their mean (blue) and the
%     5th / 95th percentiles (green / red dashed);
%   * the fitted maxima (k*, lambda*) of all realisations (blue dots) and
%     their 1-sigma covariance ellipse (orange);
%   * a box with the mean +/- std of the errors Delta k and Delta lambda.
ksize = 2;                              % Plot k in [-ksize, ksize]
kvec  = linspace(-ksize,ksize,400);

outDir = fullfile(base_dir_root,'dispersion_summary_panels');
if ~exist(outDir,'dir'), mkdir(outDir); end

% Reserve space for the shared title and the legend strip at the bottom.
legendStripHeight = 0.105;
topStripHeight    = 0.080;

for iL = 1:nL
    for it = 1:nT

        fig = figure('Color','w','Units','normalized', ...
            'Position',[0.10 0.08 0.78 0.80]);

        TL = tiledlayout(fig, nMu_plot, nS, ...
            'TileSpacing','compact','Padding','compact');

        % Shrink only the tiled plotting region so that the title and the
        % shared legend stay inside the exported figure.
        defaultPos = TL.InnerPosition;
        newBottom  = defaultPos(2) + legendStripHeight;
        newHeight  = max(defaultPos(4) - legendStripHeight - topStripHeight, 0.01);
        TL.InnerPosition = [defaultPos(1), newBottom, defaultPos(3), newHeight];

        drawnow;
        hasDispersionData = false;

        for iMu_plot = 1:nMu_plot

            iMu = idx_mu_plot(iMu_plot);
            mu  = mu_values(iMu);

            % Theoretical dispersion relation and its maximum.
            eigfun_real = @(k) mu - gamma_value*(1-k.^2).^2;
            k_max_real = 1;
            max_eigenvalue_real = mu;

            leftHasSamples = false;   % Reset for each row (see y tick labels below)

            for iS = 1:nS

                ax = nexttile(TL, (iMu_plot-1)*nS + iS);
                hold(ax,'on');
                styleAxes(ax,'TickLabelInterpreter','latex');

                % Column header: the true s.
                if iMu_plot == 1
                    title(ax, sprintf('$s=%g$',s_values(iS)), ...
                        'Interpreter','latex','FontSize',fsPanelTitle);
                end

                % Row label: mu, written vertically on the outer right.
                if iS == nS
                    addRowLabel(ax, mu, fsPanelTitle);
                end

                % The k axis is the same everywhere: label it once per column.
                if iMu_plot == nMu_plot
                    xlabel(ax,'$k$','Interpreter','latex','FontSize',fsLabel);
                end

                % ---- Data for this panel ----
                vals1 = per_theta{1,iMu,iS,iL,it};     % a0_hat
                vals2 = per_theta{2,iMu,iS,iL,it};     % a2_hat
                vals3 = per_theta{3,iMu,iS,iL,it};     % a4_hat
                k_max_data      = per_kmax{iMu,iS,iL,it};
                eigenvalue_data = per_lambda_max{iMu,iS,iL,it};

                nSim = min([numel(vals1),numel(vals2),numel(vals3), ...
                    numel(k_max_data),numel(eigenvalue_data)]);

                hasData = nSim >= minSamplesPlot;      % Fewer samples: no panel

                % y tick labels are kept on the left; the right panel hides
                % them only if the left panel has data.
                if iS == 1
                    leftHasSamples = hasData;
                end

                if ~hasData
                    showNoSamples(ax,fsTick);
                    continue
                end

                hasDispersionData = true;

                vals1 = vals1(1:nSim);
                vals2 = vals2(1:nSim);
                vals3 = vals3(1:nSim);
                k_max_data      = k_max_data(1:nSim)';
                eigenvalue_data = eigenvalue_data(1:nSim)';

                % Dispersion curve of every realisation (rows) on the k grid.
                plot_data = evalDispersion(vals1, vals2, vals3, kvec);

                meanValues   = mean(plot_data,1,'omitnan');
                percentile95 = prctile(plot_data,95,1);
                percentile5  = prctile(plot_data,5,1);

                % Fitted maxima (k*, lambda*) of all realisations.
                scatter(ax, k_max_data, eigenvalue_data, 20, 'filled', ...
                    'MarkerFaceColor',[0 0.4470 0.7410], ...
                    'MarkerEdgeColor','none', ...
                    'MarkerFaceAlpha',0.30);

                % 1-sigma covariance ellipse of those maxima.
                data = [k_max_data, eigenvalue_data];

                if size(data,1) > 2 && rank(cov(data)) == 2
                    dataMean = mean(data);
                    C = cov(data);
                    [V,D] = eig(C);
                    phi = linspace(0,2*pi,100);
                    ellipse = V*sqrt(D)*[cos(phi);sin(phi)] + dataMean';

                    plot(ax, ellipse(1,:), ellipse(2,:), ...
                        'Color',[1 .5 0],'LineWidth',2);
                end

                % Theoretical curve, mean and percentile bands of the fits.
                fplot(ax, eigfun_real, [-ksize,ksize], 'k','LineWidth',1.5);
                plot(ax, kvec, meanValues,   'b',  'LineWidth',2);
                plot(ax, kvec, percentile95, 'r--','LineWidth',1.5);
                plot(ax, kvec, percentile5,  'g--','LineWidth',1.5);

                % True theoretical maximum.
                plot(ax, k_max_real, max_eigenvalue_real, 'p', ...
                    'MarkerSize',18, ...
                    'MarkerFaceColor','yellow', ...
                    'MarkerEdgeColor','black');

                % Errors of the fitted maxima relative to the true maximum.
                delta_k      = k_max_data - k_max_real;
                delta_lambda = eigenvalue_data - max_eigenvalue_real;

                delta_k_mean = mean(delta_k,'omitnan');
                delta_k_std  = std(delta_k,'omitnan');
                delta_l_mean = mean(delta_lambda,'omitnan');
                delta_l_std  = std(delta_lambda,'omitnan');

                textboxString = sprintf(['n=%d\n', ...
                    '$\\Delta k = %.3g \\pm %.3g$\n', ...
                    '$\\Delta\\lambda = %.3g \\pm %.3g$'], ...
                    nSim, delta_k_mean, delta_k_std, delta_l_mean, delta_l_std);

                text(ax, 0.02, 0.98, textboxString, ...
                    'Units','normalized', ...
                    'Interpreter','latex', ...
                    'VerticalAlignment','top', ...
                    'BackgroundColor','w', ...
                    'EdgeColor','k', ...
                    'Margin',0.5, ...
                    'FontSize',fsTextBox);

                xlim(ax, [-ksize ksize]);

                % Keep at most three y ticks: start from MATLAB's automatic
                % (usually nicely rounded) ticks and retain three spread
                % across the range if there are more.
                currentYTicks = yticks(ax);
                if numel(currentYTicks) > 3
                    keepIdx = unique(round(linspace(1,numel(currentYTicks),3)));
                    yticks(ax, currentYTicks(keepIdx));
                end

                % Only the bottom row carries numerical x tick labels.
                if iMu_plot ~= nMu_plot
                    xticklabels(ax,[]);
                end

                box(ax,'on');
                grid(ax,'on');

                % Last: hide the y tick labels of the right panel when the
                % left panel has data. Only the labels go; ticks and grid stay.
                if iS ~= 1 && leftHasSamples
                    yticklabels(ax,[]);
                end
                hold(ax,'off');
            end
        end

        % One shared lambda(k) label for the whole figure.
        ylabel(TL,'$\lambda(k)$','Interpreter','latex','FontSize',fsLabel);

        % Shared legend at the bottom (only if at least one panel has data).
        if hasDispersionData
            ax_legend = axes(fig,'Position',[0 0 1 1], ...
                'Visible','off','HandleVisibility','off');
            hold(ax_legend,'on');

            % Dummy (NaN) graphics objects that only carry the legend entries.
            h_star = plot(ax_legend,NaN,NaN,'p', ...
                'MarkerSize',12,'MarkerFaceColor','yellow','MarkerEdgeColor','k');
            h_mean = plot(ax_legend,NaN,NaN,'b-','LineWidth',2);
            h_95   = plot(ax_legend,NaN,NaN,'r--','LineWidth',1.5);
            h_5    = plot(ax_legend,NaN,NaN,'g--','LineWidth',1.5);
            h_theo = plot(ax_legend,NaN,NaN,'k','LineWidth',1.2);
            h_dots = scatter(ax_legend,NaN,NaN,20,'filled', ...
                'MarkerFaceColor',[0 0.4470 0.7410], ...
                'MarkerFaceAlpha',0.3,'MarkerEdgeColor','none');

            leg = legend(ax_legend,[h_theo,h_mean,h_95,h_5,h_star,h_dots], ...
                {'Theoretical curve','Mean','95\%','5\%', ...
                 'True $(k_{*},\lambda_{*})$', ...
                 'Data $(k_{*},\lambda_{*})$'}, ...
                'Orientation','horizontal', ...
                'Location','southoutside', ...
                'Box','on', ...
                'FontSize',fsLegend, ...
                'Interpreter','latex');

            drawnow;
            pos = leg.Position;
            pos(2) = 0.015;
            leg.Position = pos;
            hold(ax_legend,'off');
        end

        title(TL, ...
            sprintf('Dispersion summaries $(L=%.3g\\pi,\\;\\tau_{max}=%.3g)$', ...
            L_values(iL)/pi,tmax_values(it)), ...
            'Interpreter','latex','FontSize',fsFigureTitle);

        exportgraphics(fig, ...
            fullfile(outDir,sprintf('dispersion_panel_L=%g_t=%g.jpg', ...
            L_values(iL),tmax_values(it))), ...
            'Resolution',300);

        savefig(fig, ...
            fullfile(outDir,sprintf('dispersion_panel_L=%g_t=%g.fig', ...
            L_values(iL),tmax_values(it))));
    end
end

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function [tsteps_max, tidx_ref, dt_actual] = horizonGrid(tmax_values, tmax_max, dt_target)
% HORIZONGRID  Time-grid size and truncation indices for an observation horizon.
%   A trajectory is simulated up to tmax_max with step dt_target.
%     tsteps_max  number of stored time points (including t = 0)
%     tidx_ref    for each tmax_values(it) <= tmax_max, the index of the time
%                 point t = tmax_values(it); NaN for longer horizons
%     dt_actual   realised time step (equal to dt_target if the check passes)
%   Errors if tmax_max or any tmax_values(it) is not an integer multiple of
%   dt_target, because then the truncation indices would not be exact.
    tsteps_max = round(tmax_max/dt_target) + 1;
    dt_actual  = tmax_max/(tsteps_max-1);
    if abs(dt_actual-dt_target) > 1e-12*max(1,abs(dt_target))
        error('tmax must be an integer multiple of dt_target.');
    end

    nT = numel(tmax_values);
    tidx_ref = nan(nT,1);
    for it = 1:nT
        if tmax_values(it) <= tmax_max
            idx_float = tmax_values(it)/dt_target;
            if abs(idx_float-round(idx_float)) > 1e-10
                error('Each tmax value must be an integer multiple of dt_target.');
            end
            tidx_ref(it) = round(idx_float) + 1;
        end
    end
end

function n = nextBatchSize(numToGo, numCores, batchSize)
% NEXTBATCHSIZE  Number of new seeds to launch in the next parallel batch.
%   Ask for 25% more seeds than are still needed (some runs are rejected),
%   but at least one per core, round up to a multiple of the core count so
%   that all workers are busy, and never exceed batchSize.
    overProvision = 1.25;
    tentative = ceil(overProvision * numToGo);
    tentative = max(numCores, tentative);
    tentative = min(batchSize, tentative);
    n = ceil(tentative/numCores)*numCores;
    n = min(n, batchSize);
end

function setLogMuAxis(ax, xpos)
% SETLOGMUAXIS  Log-scaled |mu| axis, shown reversed with ticks labelled -10^k.
%   xpos = -mu (positive). Decade ticks cover the range of xpos.
    set(ax,'XScale','log');

    emin = floor(log10(min(xpos)));
    emax = ceil(log10(max(xpos)));
    if emax-emin < 1
        emin = emin-1;
        emax = emax+1;
    end

    xlim(ax,[10^emin 10^emax]);

    exps = emin:emax;
    set(ax,'XTick',10.^exps);

    % The plotted coordinate is |mu|, but the physical parameter is negative.
    xtick_labels = arrayfun(@(e) sprintf('$-10^{%d}$',e), ...
        exps,'UniformOutput',false);

    set(ax,'XTickLabel',xtick_labels, ...
        'TickLabelInterpreter','latex', ...
        'XDir','reverse');
end

function addRowLabel(ax, mu, fontSize)
% ADDROWLABEL  Write "mu = value" vertically along the right edge of a panel.
    text(ax,1.025,0.5, ...
        sprintf('$\\mu=%.1f$',mu), ...
        'Units','normalized', ...
        'Interpreter','latex', ...
        'FontSize',fontSize, ...
        'Rotation',90, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'Clipping','off');
end

function showNoSamples(ax, fontSize)
% SHOWNOSAMPLES  Turn a panel into an empty frame labelled "no samples".
    text(ax,0.5,0.5,'no samples', ...
        'Units','normalized', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',fontSize);
    xlim(ax,[0 1]);
    ylim(ax,[0 1]);
    xticks(ax,[]);
    yticks(ax,[]);
    box(ax,'on');
    hold(ax,'off');
end

function D = evalDispersion(a0, a2, a4, k)
% EVALDISPERSION  Fitted dispersion curves lambda(k) = a0 - a2*k^2 + a4*k^4.
%   a0, a2, a4 are vectors with one fitted coefficient per realisation and k
%   is a row vector of wavenumbers. D has one row per realisation and one
%   column per wavenumber.
    D = a0(:) - a2(:).*k.^2 + a4(:).*k.^4;
end
