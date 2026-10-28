% SH_loops_results_gamma.m
% -------------------------------------------------------------------------
% MODEL
%   Data are generated with the stochastic Swift-Hohenberg equation
%       dq = [ mu*q - gamma*(1 + d_xx)^2 q + s*q^3 ] dtau + dW ,
%   and fitted to its deterministic part (the noise is not fitted, see
%   FITTED QUANTITIES below). The linear dispersion relation of the
%   deterministic part is lambda(k) = mu - gamma*(1 - k^2)^2.
%
% MODEL
%   dq = [ mu*q - gamma*(1 + d_xx)^2 q + s*q^3 ] dtau + dW
%   with linear dispersion relation lambda(k) = mu - gamma*(1 - k^2)^2.
%
% FITTED QUANTITIES (analyze_SH_scaled)
%   theta = [a0_hat; a2_hat; a4_hat; s_hat] for the fit
%       q_tau = a0*q + a2*q_XX + a4*q_XXXX + s*q^3,
%   with expected values a0 = mu - gamma, a2 = -2*gamma, a4 = -gamma, and s
%   the cubic coefficient. k_max and lambda_max are the wavenumber and
%   growth rate of the most unstable mode of the fitted dispersion relation.
%
% WORKFLOW
%   1. For every (mu, gamma, s) combination, simulate independent noise
%      realizations until n_realizations COMPLETE runs are collected (a run
%      is complete only if every fitted quantity is finite at every
%      requested horizon). If fewer than 10% of the first batch is usable,
%      that combination is skipped and left empty.
%   2. Compute the classification AUC of s_hat for every
%      (mu, tau_max, L, gamma) with complete data for both classes.
%   3. Save the fitted coefficients (not the raw trajectories).
%   4. Produce figures: AUC vs mu (one line per gamma), distributions of the
%      fitted coefficients, ROC curves and dispersion-relation summaries
%      (one figure per gamma).
%
% REQUIREMENTS
%   * simulate_scaled_SH_Fourier.m and analyze_SH_scaled.m on the MATLAB path.
%   * MATLAB R2020a or later (tiledlayout, exportgraphics, orderedcolors).
%   * Parallel Computing Toolbox (parfor).
%   * Statistics and Machine Learning Toolbox (perfcurve, ksdensity, prctile).
%
% OUTPUTS
%   A folder named after the key settings is created in the current
%   directory, containing
%       saved_data/                   processed MAT-file
%       plots/auc/                    AUC vs mu, all gamma
%       ridge_plots_all/              distributions of fitted coefficients
%       roc_plots_per_gamma/          ROC curves, one figure per gamma
%       dispersion_summary_panels/    estimated vs theoretical dispersion
%
% REUSING SAVED DATA
%   If a MAT-file for these settings (n_realizations, dx, dt) already
%   exists, it is loaded instead of re-simulating, and the settings stored in
%   that file (mu_values, gamma_values, ...) replace the ones set below.
%   Set rerunSimulations = true to ignore the file and simulate again.
% -------------------------------------------------------------------------

clearvars; clc;

numCores = feature('numCores');   % used to size the parallel batches

%% ------------------------
%% User settings
%% ------------------------

% --- Noise and initial condition ---
noiseType          = "White";   % "White", "Correlated", "Uniform" or "None"
correlation_length = 0;         % only used for "Correlated" noise
InitialScale       = 0;         % std. dev. of the random initial condition (0 = start at q = 0)
NoiseScale         = 1;         % noise amplitude

% --- Parameter sweep ---
mu_values    = -logspace(-1,2,20);                    % negative mu: below the instability
gamma_values = [0,0.1,0.2,0.5,1,2,5,10,20,50];       % gamma sweep
s_values     = [1, -1];                               % cubic coefficient / sign s (the two classes)

% --- Domain and discretisation ---
L_values  = 20*pi;              % periodic domain length(s)
dx_target = 0.2*pi;             % exact grid spacing: Nx = L/dx = 100
% The time step must satisfy the explicit-Euler stability limit
% dt*max(-lambda_k) < 2. The largest damping rate grows with gamma and |mu|,
% hence the smaller step than in SH_loops_results.m (largest case:
% gamma = 50, mu = -100 gives dt*max(-lambda_k) = 1.45).
dt_target = 5e-5;

% --- Observation horizons ---
% Each value must be an integer multiple of dt_target.
tmax_values = 2;

% --- Ensemble ---
n_realizations        = 100;    % complete realizations per (mu, gamma, s)
u_explosion_threshold = 1e10;   % reject trajectories that exceed this magnitude

% --- Classification score ---
% theta(4) = s_hat, the estimated cubic coefficient, is the ROC score.
kc_shat    = 4;
nThetaComp = 4;

% --- Run control ---
rerunSimulations = false;   % true: ignore any saved MAT-file and simulate again

nT = numel(tmax_values);
nS = numel(s_values);
nG = numel(gamma_values);

%% ------------------------
%% Output folders and saved-data file
%% ------------------------
base_dir_root = fullfile(pwd, sprintf( ...
    'SH_ROC_scaled_FFT_mu_s_real=%d_dx=%.3gpi_dt=%.3g_gamma', ...
    n_realizations, dx_target/pi, dt_target));
if ~exist(base_dir_root,'dir'), mkdir(base_dir_root); end

outDataDir = fullfile(base_dir_root, 'saved_data');
saveName = sprintf('SH_ROC_data_FFT_mu_s_real%d_dx%.3gpi_dt%.3g_gamma.mat', ...
    n_realizations, dx_target/pi, dt_target);
savePath = fullfile(outDataDir, saveName);

%% ------------------------
%% Load saved results if available (skips the expensive simulation stage)
%% ------------------------
loadedFromFile = false;
if exist(savePath, 'file') && ~rerunSimulations
    fprintf('Found saved data at %s. Loading and skipping simulations...\n', savePath);
    S = load(savePath);

    requiredFields = {'per_theta','per_theta_seeds','per_kmax','per_lambda_max', ...
        'mu_values','tmax_values','L_values','s_values','gamma_values', ...
        'kc_shat','nThetaComp','n_realizations'};
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
    gamma_values    = S.gamma_values;
    kc_shat         = S.kc_shat;
    nThetaComp      = S.nThetaComp;
    n_realizations  = S.n_realizations;

    nMu = numel(mu_values);
    nT  = numel(tmax_values);
    nL  = numel(L_values);
    nS  = numel(s_values);
    nG  = numel(gamma_values);
    loadedFromFile = true;
end

if ~loadedFromFile
    %% ============================================================
    %% SIMULATION AND ESTIMATION
    %% ============================================================
    nMu = numel(mu_values);
    nL  = numel(L_values);

    % Storage for the fitted quantities. Each cell holds a row vector with
    % exactly n_realizations entries, or stays empty if the case was skipped:
    %   per_theta{comp, iMu, iS, iL, it, ig}   fitted coefficient 'comp' (1..nThetaComp)
    %   per_theta_seeds{...}                   rng seeds of those realizations
    %   per_kmax{iMu, iS, iL, it, ig}          k_max of the fitted dispersion relation
    %   per_lambda_max{iMu, iS, iL, it, ig}    lambda_max of the fitted dispersion relation
    per_theta       = cell(nThetaComp, nMu, nS, nL, nT, nG);
    per_theta_seeds = cell(nThetaComp, nMu, nS, nL, nT, nG);
    per_kmax        = cell(nMu, nS, nL, nT, nG);
    per_lambda_max  = cell(nMu, nS, nL, nT, nG);

    batchSize = n_realizations;      % maximum number of seeds per batch

    for iMu = 1:nMu
        mu = mu_values(iMu);
        fprintf('mu (%d/%d) = %g\n', iMu, nMu, mu);

        for iL = 1:nL
            L = L_values(iL);

            % Periodic FFT grid: Nx unique points on [0,L), so x = L is not duplicated.
            xsteps = round(L/dx_target);
            dx_actual = L/xsteps;
            if abs(dx_actual-dx_target) > 1e-12*max(1,abs(dx_target))
                error('L must be an integer multiple of dx_target.');
            end

            % All requested horizons come from one run of the longest length.
            tmax_max = max(tmax_values);
            [tsteps_max, tidx_ref] = time_indices(tmax_values, tmax_max, dt_target);

            fprintf('L=%.6g, Nx=%d, dx=%.12g, Nt=%d, dt=%.12g\n', ...
                L, xsteps, dx_actual, tsteps_max, dt_target);

            for ig = 1:nG
                gamma_now = gamma_values(ig);
                fprintf('\ngamma = %g\n', gamma_now);

                for iS = 1:nS
                    s_true = s_values(iS);
                    fprintf('s = %g\n', s_true);

                    seedStart    = 0;
                    batch        = 0;
                    collectedIdx = 0;

                    % Complete realizations collected so far (in increasing seed order).
                    collectedTheta  = nan(n_realizations, nT, nThetaComp);
                    collectedKmax   = nan(n_realizations, nT);
                    collectedLambda = nan(n_realizations, nT);
                    collectedSeeds  = nan(n_realizations, 1);

                    %% Batch loop: simulate until n_realizations complete runs are collected
                    while collectedIdx < n_realizations
                        batch = batch + 1;

                        % Over-provision seeds to allow for rejected runs.
                        numSeedsThisBatch = batch_size(n_realizations-collectedIdx, numCores, batchSize);
                        seeds = seedStart + (1:numSeedsThisBatch);
                        seedStart = seedStart + numSeedsThisBatch;

                        theta_batch  = nan(numel(seeds), nT, nThetaComp);
                        kmax_batch   = nan(numel(seeds), nT);
                        lambda_batch = nan(numel(seeds), nT);

                        fprintf('Trying gamma=%g seeds %d -- %d\n', gamma_now, seeds(1), seeds(end));

                        %% Parallel simulation + fitting of one realization per seed
                        parfor ib = 1:numel(seeds)
                            seed = seeds(ib);
                            rng(seed,'combRecursive');   % independent, reproducible stream per seed

                            try
                                [tvec_full,xvec_full,u_full] = simulate_scaled_SH_Fourier( ...
                                    mu, s_true, gamma_now, tmax_max, tsteps_max, L, xsteps, ...
                                    noiseType, correlation_length, NoiseScale, InitialScale);
                            catch
                                % Any error skips this seed silently. This includes the
                                % explicit-Euler stability check in simulate_scaled_SH_Fourier:
                                % if dt_target is too large for the largest gamma and |mu|,
                                % EVERY seed fails and the case is skipped with the "only
                                % collected 0/..." warning below. Test one call without
                                % try/catch to see the real error message.
                                continue
                            end

                            if isempty(u_full)
                                continue
                            end

                            % Make sure time runs along the rows.
                            [n1,n2] = size(u_full);
                            if n1 == numel(tvec_full)
                                u_sim = u_full;
                            elseif n2 == numel(tvec_full)
                                u_sim = u_full';
                            else
                                continue
                            end

                            for it = 1:nT
                                idx_end = tidx_ref(it);
                                if isnan(idx_end) || idx_end < 3
                                    continue
                                end

                                t_sub = tvec_full(1:idx_end);
                                u_sub = u_sim(1:idx_end,:);

                                % Reject non-finite or diverging trajectories.
                                if any(~isfinite(u_sub(:))) || max(abs(u_sub(:))) > u_explosion_threshold
                                    continue
                                end

                                try
                                    [theta_sub,kmax_sub,lambda_sub] = ...
                                        analyze_SH_scaled(u_sub,t_sub,xvec_full);
                                catch
                                    continue
                                end

                                theta_sub = theta_sub(:);
                                if numel(theta_sub) < nThetaComp
                                    theta_sub(end+1:nThetaComp,1) = NaN;
                                end

                                theta_batch(ib,it,:) = theta_sub(1:nThetaComp);
                                kmax_batch(ib,it)    = kmax_sub;
                                lambda_batch(ib,it)  = lambda_sub;
                            end
                        end % parfor

                        % Keep only complete realizations, in increasing seed order.
                        % A realization counts only when every fitted coefficient,
                        % k_max and lambda_max is finite at every requested horizon.
                        for ib = 1:numel(seeds)
                            if collectedIdx >= n_realizations
                                break
                            end

                            theta_this  = squeeze(theta_batch(ib,:,:));
                            kmax_this   = kmax_batch(ib,:);
                            lambda_this = lambda_batch(ib,:);

                            if any(~isfinite(theta_this(:))) || ...
                                    any(~isfinite(kmax_this(:))) || ...
                                    any(~isfinite(lambda_this(:)))
                                continue
                            end

                            collectedIdx = collectedIdx + 1;
                            collectedTheta(collectedIdx,:,:) = theta_batch(ib,:,:);
                            collectedKmax(collectedIdx,:)    = kmax_batch(ib,:);
                            collectedLambda(collectedIdx,:)  = lambda_batch(ib,:);
                            collectedSeeds(collectedIdx)     = seeds(ib);
                        end

                        fprintf('gamma %g batch %d: collected %d / %d\n', ...
                            gamma_now, batch, collectedIdx, n_realizations);

                        % If fewer than 10% are collected in the first batch,
                        % skip this (mu, gamma, s) case.
                        if batch == 1 && collectedIdx < 0.1*n_realizations
                            warning(['mu=%g, L=%g, s=%g, gamma=%g: only collected %d/%d ', ...
                                'in the first batch; skipping this case.'], ...
                                mu, L, s_true, gamma_now, collectedIdx, n_realizations);
                            break
                        end
                    end % while

                    if collectedIdx < n_realizations
                        % This case was deliberately skipped; leave its cells empty.
                        continue
                    end

                    % Store exactly the same n_realizations complete runs for
                    % every fitted coefficient and dispersion quantity.
                    seedVals = collectedSeeds(1:n_realizations);

                    for it = 1:nT
                        for kc = 1:nThetaComp
                            vals = squeeze(collectedTheta(1:n_realizations,it,kc));
                            per_theta{kc,iMu,iS,iL,it,ig}       = vals(:)';
                            per_theta_seeds{kc,iMu,iS,iL,it,ig} = seedVals(:)';
                        end

                        valsK = collectedKmax(1:n_realizations,it);
                        per_kmax{iMu,iS,iL,it,ig} = valsK(:)';

                        valsL = collectedLambda(1:n_realizations,it);
                        per_lambda_max{iMu,iS,iL,it,ig} = valsL(:)';
                    end
                end % s
            end % gamma
        end % L
    end % mu

else
    fprintf('Proceeding to postprocessing and plotting using loaded data.\n');
end

%% ============================================================
%% CLASSIFICATION PERFORMANCE (AUC)
% The fitted cubic coefficient s_hat is the score that separates the two
% signs of the nonlinear term; s > 0 is the positive ROC class. AUC = 1 is
% perfect separation, AUC = 0.5 is chance. A (mu, tau_max, L, gamma) point
% is evaluated only if BOTH classes have n_realizations complete fits;
% otherwise its AUC stays NaN.
%% ============================================================
fprintf('\nComputing AUC...\n');

if kc_shat > nThetaComp
    error('kc_shat (%d) exceeds nThetaComp (%d).', kc_shat, nThetaComp);
end

auc_all_all = nan(nMu, nT, nL, nG);

for iMu = 1:nMu
    for iL = 1:nL
        for it = 1:nT
            for ig = 1:nG
                [scores, labels, complete] = collect_complete_scores( ...
                    per_theta, kc_shat, iMu, iL, it, ig, s_values, n_realizations);
                if ~complete
                    continue
                end

                [~,~,~,auc_all_all(iMu,it,iL,ig)] = perfcurve(labels, scores, 1);
            end
        end
    end
end

%% ============================================================
%% SAVE PROCESSED RESULTS
% (MAT-file v7 is used: the files are only a few MB, and v7.3/HDF5 can be
% problematic when writing inside a synchronised cloud folder.)
%% ============================================================
if ~loadedFromFile
    if ~exist(outDataDir,'dir'), mkdir(outDataDir); end

    save(savePath, ...
        'per_theta','per_theta_seeds','per_kmax','per_lambda_max', ...
        'mu_values','tmax_values','L_values','s_values','gamma_values', ...
        'kc_shat','nThetaComp','n_realizations','auc_all_all', ...
        'dx_target','dt_target','base_dir_root','-v7');

    fprintf('Saved results to: %s\n', savePath);
else
    fprintf('Loaded data present: skipping save to avoid overwriting existing file.\n');
end

%% ============================================================
%% FIGURE STYLE
% Everything below is visualisation of the quantities computed above.
%% ============================================================
paperFont     = 'Helvetica';

fsTick        = 13;
fsLabel       = 15;
fsPanelTitle  = 14;
fsFigureTitle = 17;
fsLegend      = 12;
fsTextBox     = 10;

axisLineWidth = 1.1;

maxMuPlot = 5;   % number of mu values shown in the panel figures

% mu values used in the ridge and dispersion panel figures (spread across the sweep)
idx_mu_plot = unique(round(linspace(1, nMu, min(maxMuPlot, nMu))));
nMu_plot    = numel(idx_mu_plot);

cmapPos = [0 0.4470 0.7410];      % colour used for s > 0
cmapNeg = [0.8500 0.3250 0.0980]; % colour used for s < 0

%% ============================================================
%% AUC vs mu: one coloured line per gamma
%% ============================================================
iL_plot_idx = 1;   % domain length shown
it_plot_idx = 1;   % observation horizon shown

fig = figure('Name','AUC vs mu - all gamma','Color','w','Units','normalized', ...
    'Position',[0.15 0.15 0.7 0.6]);
ax = axes(fig); hold(ax,'on');
set(ax,'FontName',paperFont,'FontSize',fsTick,'Box','on', ...
    'LineWidth',axisLineWidth,'TickDir','in');

colorsG = get_colors(nG);

% MATLAB log axes cannot show negative values, so |mu| is plotted on the
% axis and the ticks are labelled with the negative mu values (see format_mu_axis).
if mu_axis_is_log(mu_values)
    xplot = -mu_values(:);
else
    xplot = mu_values(:);
end

hPlot  = gobjects(nG,1);
labels = cell(nG,1);

for ig = 1:nG
    auc_vec = squeeze(auc_all_all(:,it_plot_idx,iL_plot_idx,ig));
    valid   = isfinite(auc_vec(:)) & isfinite(xplot);
    if ~any(valid), continue; end

    color = colorsG(ig,:);

    % Line through all mu, plus dots marking the mu values actually computed.
    hPlot(ig) = plot(ax, xplot, auc_vec(:), 'Color', color, 'LineWidth', 1.8);
    scatter(ax, xplot(valid), auc_vec(valid), 36, 'o', ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', color*0.6, ...
        'MarkerFaceAlpha', 0.9, 'HandleVisibility', 'off');

    labels{ig} = sprintf('$\\gamma=%.3g$', gamma_values(ig));
end

% Keep only the lines that were actually drawn.
validIdx = isgraphics(hPlot);
hPlot  = hPlot(validIdx);
labels = labels(validIdx);

format_mu_axis(ax, mu_values);

xlabel(ax, '$\mu$', 'FontSize', fsLabel, 'Interpreter','latex');
ylabel(ax, 'AUC', 'FontSize', fsLabel);
title(ax, sprintf('AUC vs $\\mu$ $(L=%.3g\\pi,\\;\\tau_{max}=%.3g)$', ...
    L_values(iL_plot_idx)/pi, tmax_values(it_plot_idx)), ...
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

outDirAUC = fullfile(base_dir_root,'plots','auc');
if ~exist(outDirAUC,'dir'), mkdir(outDirAUC); end
save_figure(fig, fullfile(outDirAUC,'AUC_vs_mu_all_gamma.jpg'), true);

%% ============================================================
%% RIDGE / DENSITY PANELS: one figure per (coefficient, gamma, L, tau_max)
% Distribution of a fitted coefficient for each true class s, in a grid of
% panels (rows: mu, columns: s). Each density is a kernel-density estimate
% scaled to unit height, with the sample mean as a dashed line. These panels
% show the overlap between the two classes directly, complementing the
% single-number AUC summary.
%
% The coefficients shown are a0_hat (compared with mu - gamma) and s_hat
% (the classification score). Only complete datasets (n_realizations
% samples) are drawn.
%% ============================================================
outDirRidge = fullfile(base_dir_root,'ridge_plots_all');
if ~exist(outDirRidge,'dir'), mkdir(outDirRidge); end

kc_list = unique([1, kc_shat]);

for kc = kc_list

    if kc == 1
        compLabel = '$\hat{\mu}-\gamma$';
        plotTag   = 'mu_minus_gamma';
    else
        compLabel = '$\hat{s}$';
        plotTag   = 's_hat';
    end

    for ig = 1:nG
        gamma_now = gamma_values(ig);

        for iL_plot = 1:nL
            for it_plot = 1:nT

                figRidge = figure('Name', sprintf('Ridge %s L=%g t=%g gamma=%g', ...
                    plotTag, L_values(iL_plot)/pi, tmax_values(it_plot), gamma_now), ...
                    'Color','w','Units','normalized','Position',[0.08 0.07 0.82 0.82]);

                t = tiledlayout(figRidge, nMu_plot, nS, ...
                    'TileSpacing','compact','Padding','compact');

                for iMu_plot = 1:nMu_plot
                    iMu = idx_mu_plot(iMu_plot);
                    mu  = mu_values(iMu);

                    % Fixed x-window per row, so the two classes can be compared.
                    if kc == 1
                        xLimits = (mu - gamma_now) + [-10 10];   % centred on the theoretical value
                    else
                        xLimits = [-30 30];
                    end

                    leftHasSamples = false;   % reset for every row

                    for iS = 1:nS
                        ax = nexttile(t, (iMu_plot-1)*nS + iS);
                        hold(ax,'on');
                        set(ax,'FontName',paperFont,'FontSize',fsTick, ...
                            'LineWidth',axisLineWidth,'TickDir','in', ...
                            'TickLabelInterpreter','latex');

                        scores = per_theta{kc, iMu, iS, iL_plot, it_plot, ig};
                        scores = scores(isfinite(scores));
                        scores = scores(:);
                        hasData = numel(scores) >= n_realizations;
                        if hasData
                            scores = scores(1:n_realizations);
                        end

                        % Column headers show s; row labels (outer right) show mu.
                        if iMu_plot == 1
                            title(ax, sprintf('$s=%g$', s_values(iS)), ...
                                'Interpreter','latex','FontSize',fsPanelTitle);
                        end
                        if iS == nS
                            text(ax, 1.025, 0.5, sprintf('$\\mu=%.1f$', mu), ...
                                'Units','normalized','Interpreter','latex', ...
                                'FontSize',fsPanelTitle,'Rotation',90, ...
                                'HorizontalAlignment','center', ...
                                'VerticalAlignment','middle','Clipping','off');
                        end
                        % The variable label goes below each column in the bottom row.
                        if iMu_plot == nMu_plot
                            xlabel(ax, compLabel, 'Interpreter','latex','FontSize',fsLabel);
                        end

                        % Draw the density (if the dataset is complete).
                        if hasData
                            if s_values(iS) < 0
                                color = cmapNeg;
                            else
                                color = cmapPos;
                            end
                            hasData = draw_density_panel(ax, scores, color, xLimits, fsTextBox);
                        end

                        if hasData
                            grid(ax,'on'); box(ax,'on');
                        else
                            draw_no_samples(ax, fsTick);
                        end

                        if iS == 1
                            leftHasSamples = hasData;
                        end

                        % Show y tick labels only in the left column; the right
                        % column hides them only when the left panel has data.
                        % (Ticks and grid are kept.)
                        if iS > 1 && leftHasSamples
                            yticklabels(ax, []);
                        end

                        hold(ax,'off');
                    end
                end

                % Shared labels and title.
                ylabel(t, 'Normalized density', 'FontSize', fsLabel);
                title(t, sprintf( ...
                    '%s distributions $(L=%.3g\\pi,\\;\\tau_{max}=%.3g,\\;\\gamma=%.3g)$', ...
                    compLabel, L_values(iL_plot)/pi, tmax_values(it_plot), gamma_now), ...
                    'Interpreter','latex','FontSize',fsFigureTitle);

                save_figure(figRidge, fullfile(outDirRidge, sprintf( ...
                    'ridge_%s_L=%g_t=%g_gamma=%g.jpg', plotTag, ...
                    L_values(iL_plot), tmax_values(it_plot), gamma_now)), ...
                    false);
            end
        end
    end
end

%% ============================================================
%% ROC CURVES: one figure per gamma
% Each figure shows the ROC curve of s_hat for ten representative mu values
% (from -100 to -0.1), labelled with the corresponding AUC. Only mu values
% with complete data for both classes are drawn.
%% ============================================================
iL_plot = 1;   % domain length shown
it_plot = 1;   % observation horizon shown

outDirROC = fullfile(base_dir_root,'roc_plots_per_gamma');
if ~exist(outDirROC,'dir'), mkdir(outDirROC); end

% Pick the computed mu values closest (on a log scale) to ten targets
% between -100 and -0.1.
nROCcurves    = 10;
target_mu_roc = -logspace(2,-1,nROCcurves);
idx_mu_roc    = zeros(1,nROCcurves);
for j = 1:nROCcurves
    [~,idx_mu_roc(j)] = min(abs(log10(-mu_values) - log10(-target_mu_roc(j))));
end
idx_mu_roc = unique(idx_mu_roc,'stable');

colorsROC = get_colors(numel(idx_mu_roc));

for ig = 1:nG

    figROC = figure('Name', sprintf('ROC gamma=%g', gamma_values(ig)), ...
        'Color','w','Units','normalized','Position',[0.18 0.18 0.52 0.58]);
    axROC = axes(figROC); hold(axROC,'on');
    set(axROC,'FontName',paperFont,'FontSize',fsTick, ...
        'LineWidth',axisLineWidth,'TickDir','in', ...
        'TickLabelInterpreter','latex');

    % Chance line
    plot(axROC, [0 1], [0 1], 'k--', 'LineWidth', 1, 'HandleVisibility','off');

    hasROC = false;

    for j = 1:numel(idx_mu_roc)
        iMu = idx_mu_roc(j);

        [scores, labels, complete] = collect_complete_scores( ...
            per_theta, kc_shat, iMu, iL_plot, it_plot, ig, s_values, n_realizations);
        if ~complete
            continue
        end

        [fpr, tpr, ~, auc] = perfcurve(labels, scores, 1);

        plot(axROC, fpr, tpr, 'Color', colorsROC(j,:), 'LineWidth', 1.6, ...
            'DisplayName', sprintf('\\mu=%.3g (AUC=%.2f)', mu_values(iMu), auc));
        hasROC = true;
    end

    xlabel(axROC, 'False Positive Rate', 'FontSize', fsLabel);
    ylabel(axROC, 'True Positive Rate', 'FontSize', fsLabel);
    title(axROC, sprintf('ROC $(L=%.3g\\pi,\\;\\tau_{max}=%.3g,\\;\\gamma=%.3g)$', ...
        L_values(iL_plot)/pi, tmax_values(it_plot), gamma_values(ig)), ...
        'Interpreter','latex','FontSize',fsFigureTitle);

    axis(axROC, [0 1 0 1]);
    box(axROC,'on');
    grid(axROC,'on');

    if hasROC
        lgROC = legend(axROC, 'Location','southeast', 'Interpreter','latex', ...
            'FontSize', fsLegend+1);
        lgROC.Box = 'on';

        save_figure(figROC, fullfile(outDirROC, sprintf('ROC_L%g_t%g_gamma%g.jpg', ...
            L_values(iL_plot), tmax_values(it_plot), gamma_values(ig))), ...
            true);
    else
        warning('No valid ROC curves for gamma=%g.', gamma_values(ig));
    end
end

%% ============================================================
%% DISPERSION SUMMARY PANELS: one figure per (gamma, L, tau_max)
% A grid of panels (rows: mu, columns: s) comparing the estimated
% dispersion relations with the theoretical one,
%     lambda(k) = mu - gamma*(1 - k^2)^2   (maximum lambda = mu at k = 1).
% Each panel shows
%   - black line: theoretical curve; yellow star: its maximum
%   - blue line:  mean of the fitted curves a0 - a2*k^2 + a4*k^4
%   - red/green dashed: 95th/5th percentiles of the fitted curves
%   - blue dots:  (k_max, lambda_max) of each fit, with a 1-sigma covariance
%                 ellipse (orange)
%   - text box:   mean +/- std of the errors of k_max and lambda_max
% Only complete datasets (n_realizations samples) are drawn.
%% ============================================================
ksize = 2;                          % wavenumbers plotted: [-ksize, ksize]
kvec  = linspace(-ksize, ksize, 400);

minSamplesPlot = 3;

outDirDisp = fullfile(base_dir_root,'dispersion_summary_panels');
if ~exist(outDirDisp,'dir'), mkdir(outDirDisp); end

%% Common y-axis limits per mu, taken from the gamma = 1 data
% All gamma figures of the same mu share one y range (the 5th--95th
% percentile band of the fitted curves at the gamma closest to 1, plus a 5%
% margin), so the panels can be compared across gamma. Curves at other gamma
% that leave this range are clipped by the axes.
[~, ig_yref] = min(abs(gamma_values - 1));

yLimitsPerMu = nan(nMu, nL, nT, 2);

for iL_tmp = 1:nL
    for it_tmp = 1:nT
        for iMu = 1:nMu

            ymin_ref = inf;
            ymax_ref = -inf;

            for iS_tmp = 1:nS
                vals1 = per_theta{1,iMu,iS_tmp,iL_tmp,it_tmp,ig_yref};
                vals2 = per_theta{2,iMu,iS_tmp,iL_tmp,it_tmp,ig_yref};
                vals3 = per_theta{3,iMu,iS_tmp,iL_tmp,it_tmp,ig_yref};

                nSim = min([numel(vals1), numel(vals2), numel(vals3), n_realizations]);
                if nSim < minSamplesPlot
                    continue
                end

                curves = dispersion_curves(vals1(1:nSim), vals2(1:nSim), vals3(1:nSim), kvec);

                ymax_ref = max(ymax_ref, max(prctile(curves, 95, 1)));
                ymin_ref = min(ymin_ref, min(prctile(curves, 5, 1)));
            end

            if isfinite(ymin_ref) && isfinite(ymax_ref)
                yrange = ymax_ref - ymin_ref;
                if yrange <= 0
                    yrange = max(1, abs(ymax_ref));
                end
                margin = 0.05*yrange;
                yLimitsPerMu(iMu,iL_tmp,it_tmp,:) = [ymin_ref - margin, ymax_ref + margin];
            end
        end
    end
end

% Space reserved below/above the tiles for the shared legend and title.
legendStripHeight = 0.105;
topStripHeight    = 0.080;

for ig = 1:nG

    gamma_now = gamma_values(ig);

    for iL = 1:nL
        for it = 1:nT

            figDisp = figure('Color','w','Units','normalized', ...
                'Position',[0.10 0.08 0.78 0.80]);

            TL = tiledlayout(figDisp, nMu_plot, nS, ...
                'TileSpacing','compact','Padding','compact');

            % Shrink only the tile region so the large title and the shared
            % legend stay inside the exported figure.
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
                eigfun_real         = @(k) mu - gamma_now*(1 - k.^2).^2;
                k_max_real          = 1;
                max_eigenvalue_real = mu;

                leftHasSamples = false;   % reset for every row

                for iS = 1:nS

                    ax = nexttile(TL, (iMu_plot-1)*nS + iS);
                    hold(ax,'on');
                    set(ax,'FontName',paperFont,'FontSize',fsTick, ...
                        'LineWidth',axisLineWidth,'TickDir','in', ...
                        'TickLabelInterpreter','latex');

                    % Column headers show s; row labels (outer right) show mu.
                    if iMu_plot == 1
                        title(ax, sprintf('$s=%g$', s_values(iS)), ...
                            'Interpreter','latex','FontSize',fsPanelTitle);
                    end
                    if iS == nS
                        text(ax, 1.025, 0.5, sprintf('$\\mu=%.1f$', mu), ...
                            'Units','normalized','Interpreter','latex', ...
                            'FontSize',fsPanelTitle,'Rotation',90, ...
                            'HorizontalAlignment','center', ...
                            'VerticalAlignment','middle','Clipping','off');
                    end
                    % The k axis is shared: one label below each column.
                    if iMu_plot == nMu_plot
                        xlabel(ax, '$k$', 'Interpreter','latex','FontSize',fsLabel);
                    end

                    %% Fitted coefficients and (k_max, lambda_max) for this panel
                    vals1 = per_theta{1,iMu,iS,iL,it,ig};
                    vals2 = per_theta{2,iMu,iS,iL,it,ig};
                    vals3 = per_theta{3,iMu,iS,iL,it,ig};
                    k_max_data      = per_kmax{iMu,iS,iL,it,ig};
                    eigenvalue_data = per_lambda_max{iMu,iS,iL,it,ig};

                    nSimAvailable = min([numel(vals1), numel(vals2), numel(vals3), ...
                        numel(k_max_data), numel(eigenvalue_data)]);

                    if nSimAvailable < n_realizations
                        draw_no_samples(ax, fsTick);
                        hold(ax,'off');
                        continue
                    end

                    hasDispersionData = true;
                    if iS == 1
                        leftHasSamples = true;
                    end

                    nSim = n_realizations;
                    vals1 = vals1(1:nSim);
                    vals2 = vals2(1:nSim);
                    vals3 = vals3(1:nSim);
                    k_max_data      = k_max_data(1:nSim);
                    eigenvalue_data = eigenvalue_data(1:nSim);
                    k_max_data      = k_max_data(:);
                    eigenvalue_data = eigenvalue_data(:);

                    % Fitted dispersion curves and their summary statistics.
                    plot_data    = dispersion_curves(vals1, vals2, vals3, kvec);
                    meanValues   = mean(plot_data, 1, 'omitnan');
                    percentile95 = prctile(plot_data, 95, 1);
                    percentile5  = prctile(plot_data, 5, 1);

                    %% Estimated (k_max, lambda_max) of every realization
                    scatter(ax, k_max_data, eigenvalue_data, 20, 'filled', ...
                        'MarkerFaceColor',[0 0.4470 0.7410], ...
                        'MarkerEdgeColor','none','MarkerFaceAlpha',0.30);

                    % 1-sigma covariance ellipse of the scatter.
                    data = [k_max_data, eigenvalue_data];
                    if size(data,1) > 2 && rank(cov(data)) == 2
                        dataMean = mean(data);
                        [V, D]   = eig(cov(data));
                        angles   = linspace(0, 2*pi, 100);
                        ellipse  = V*sqrt(D)*[cos(angles); sin(angles)] + dataMean';
                        plot(ax, ellipse(1,:), ellipse(2,:), 'Color',[1 .5 0], 'LineWidth',2);
                    end

                    %% Dispersion curves
                    fplot(ax, eigfun_real, [-ksize, ksize], 'k', 'LineWidth', 1.5);
                    plot(ax, kvec, meanValues,   'b',   'LineWidth', 2);
                    plot(ax, kvec, percentile95, 'r--', 'LineWidth', 1.5);
                    plot(ax, kvec, percentile5,  'g--', 'LineWidth', 1.5);

                    % True maximum of the theoretical curve.
                    plot(ax, k_max_real, max_eigenvalue_real, 'p', ...
                        'MarkerSize',18,'MarkerFaceColor','yellow','MarkerEdgeColor','black');

                    %% Error statistics of the estimated maximum
                    delta_k      = k_max_data - k_max_real;
                    delta_lambda = eigenvalue_data - max_eigenvalue_real;

                    textboxString = sprintf(['n=%d\n', ...
                        '$\\Delta k = %.3g \\pm %.3g$\n', ...
                        '$\\Delta\\lambda = %.3g \\pm %.3g$'], ...
                        nSim, mean(delta_k,'omitnan'), std(delta_k,'omitnan'), ...
                        mean(delta_lambda,'omitnan'), std(delta_lambda,'omitnan'));

                    text(ax, 0.02, 0.98, textboxString, ...
                        'Units','normalized','Interpreter','latex', ...
                        'VerticalAlignment','top','BackgroundColor','w', ...
                        'EdgeColor','k','Margin',0.5,'FontSize',fsTextBox);

                    xlim(ax, [-ksize ksize]);

                    % Common y range for this mu (see above). If no reference
                    % data exist for this mu, MATLAB's automatic limits are kept.
                    yl = squeeze(yLimitsPerMu(iMu,iL,it,:))';
                    if all(isfinite(yl))
                        ylim(ax, yl);
                    end

                    % Keep at most three y ticks, spread across MATLAB's
                    % automatically chosen (usually nicely rounded) ticks.
                    currentYTicks = yticks(ax);
                    if numel(currentYTicks) > 3
                        keepIdx = unique(round(linspace(1, numel(currentYTicks), 3)));
                        yticks(ax, currentYTicks(keepIdx));
                    end

                    % Only the bottom row carries numerical x tick labels.
                    if iMu_plot ~= nMu_plot
                        xticklabels(ax, []);
                    end

                    box(ax,'on');
                    grid(ax,'on');

                    % Show y tick labels only in the left column; the right
                    % column hides them only when the left panel has data.
                    % Done last, so ticks and grid stay.
                    if iS > 1 && leftHasSamples
                        yticklabels(ax, []);
                    end

                    hold(ax,'off');
                end
            end

            % lambda(k) is a single shared vertical label for the figure.
            ylabel(TL, '$\lambda(k)$', 'Interpreter','latex', 'FontSize', fsLabel);

            %% Shared legend at the bottom (only if at least one panel has data)
            if hasDispersionData
                ax_legend = axes(figDisp, 'Position',[0 0 1 1], ...
                    'Visible','off','HandleVisibility','off');
                hold(ax_legend,'on');

                % Invisible dummy objects that only provide legend entries.
                h_star = plot(ax_legend, NaN, NaN, 'p', ...
                    'MarkerSize',12,'MarkerFaceColor','yellow','MarkerEdgeColor','k');
                h_mean = plot(ax_legend, NaN, NaN, 'b-',  'LineWidth', 2);
                h_95   = plot(ax_legend, NaN, NaN, 'r--', 'LineWidth', 1.5);
                h_5    = plot(ax_legend, NaN, NaN, 'g--', 'LineWidth', 1.5);
                h_theo = plot(ax_legend, NaN, NaN, 'k',   'LineWidth', 1.2);
                h_dots = scatter(ax_legend, NaN, NaN, 20, 'filled', ...
                    'MarkerFaceColor',[0 0.4470 0.7410], ...
                    'MarkerFaceAlpha',0.3,'MarkerEdgeColor','none');

                leg = legend(ax_legend, [h_theo,h_mean,h_95,h_5,h_star,h_dots], ...
                    {'Theoretical curve','Mean','95\%','5\%', ...
                     'True $(k_{*},\lambda_{*})$', ...
                     'Data $(k_{*},\lambda_{*})$'}, ...
                    'Orientation','horizontal','Location','southoutside', ...
                    'Box','on','FontSize',fsLegend,'Interpreter','latex');

                drawnow;
                pos = leg.Position;
                pos(2) = 0.015;
                leg.Position = pos;
                hold(ax_legend,'off');
            end

            title(TL, sprintf( ...
                'Dispersion summaries $(L=%.3g\\pi,\\;\\tau_{max}=%.3g,\\;\\gamma=%.3g)$', ...
                L_values(iL)/pi, tmax_values(it), gamma_now), ...
                'Interpreter','latex','FontSize',fsFigureTitle);

            save_figure(figDisp, fullfile(outDirDisp, sprintf( ...
                'dispersion_panel_L=%g_t=%g_gamma=%g.jpg', ...
                L_values(iL), tmax_values(it), gamma_now)), ...
                true);
        end
    end
end

fprintf('\nDone. Results are in %s\n', base_dir_root);


%% ========================================================================
%% LOCAL FUNCTIONS
%% (In a script file, local functions must come after all other code.)
%% ========================================================================

function [tsteps, tidx] = time_indices(tmax_values, tmax_active, dt_target)
% Time-grid bookkeeping for the longest horizon.
%   tsteps  number of stored time points for a run of length tmax_active
%   tidx    row index at which each requested horizon ends (NaN for horizons
%           longer than tmax_active)
% Errors if a horizon is not an integer multiple of dt_target.

tsteps    = round(tmax_active/dt_target) + 1;
dt_actual = tmax_active/(tsteps-1);
if abs(dt_actual-dt_target) > 1e-12*max(1,abs(dt_target))
    error('tmax must be an integer multiple of dt_target.');
end

tidx = nan(numel(tmax_values),1);
for it = 1:numel(tmax_values)
    if tmax_values(it) <= tmax_active
        idx_float = tmax_values(it)/dt_target;
        if abs(idx_float-round(idx_float)) > 1e-10
            error('Each tmax value must be an integer multiple of dt_target.');
        end
        tidx(it) = round(idx_float) + 1;
    end
end
end


function n = batch_size(numToGo, numCores, maxBatch)
% Number of seeds to try in the next batch: over-provision by 25% to allow
% for rejected runs, round up to a multiple of the number of cores, and
% never exceed maxBatch.
n = ceil(1.25*numToGo);
n = max(numCores, n);
n = min(maxBatch, n);
n = ceil(n/numCores)*numCores;
n = min(maxBatch, n);
end


function [scores, labels, complete] = collect_complete_scores( ...
    per_theta, kc, iMu, iL, it, ig, s_values, n)
% Pool the fitted coefficient 'kc' from all classes into one score vector
% with logical labels (true for s > 0), using exactly n finite fits per
% class. If any class has fewer than n, complete = false and scores/labels
% are empty.
scores   = [];
labels   = [];
complete = true;

for iS = 1:numel(s_values)
    vals = per_theta{kc,iMu,iS,iL,it,ig};
    vals = vals(isfinite(vals));

    if numel(vals) < n
        scores   = [];
        labels   = [];
        complete = false;
        return
    end

    vals   = vals(1:n);
    scores = [scores; vals(:)];                         %#ok<AGROW>
    labels = [labels; repmat(s_values(iS) > 0, n, 1)];  %#ok<AGROW>
end
end


function curves = dispersion_curves(a0, a2, a4, kvec)
% Fitted dispersion relations lambda(k) = a0 - a2*k^2 + a4*k^4 on the grid
% kvec (row vector); one row per realization.
curves = a0(:) - a2(:)*kvec.^2 + a4(:)*kvec.^4;
end


function colors = get_colors(n)
% n distinct colours from MATLAB's 12-colour "gem12" palette (cycled if n > 12).
base   = orderedcolors("gem12");
colors = base(mod(0:n-1, size(base,1)) + 1, :);
end


function tf = mu_axis_is_log(mu_values)
% True if the (negative) mu sweep spans more than one decade, in which case
% |mu| is shown on a logarithmic axis.
tf = all(mu_values < 0) && (max(-mu_values)/min(-mu_values) > 10);
end


function format_mu_axis(ax, mu_values)
% Format the x axis of a mu plot. For a wide negative sweep, |mu| is plotted
% on a reversed log axis (so mu increases towards 0 on the right) and the
% ticks are labelled -10^n. Otherwise the plain mu range is used.
if mu_axis_is_log(mu_values)
    xpos = -mu_values(:);
    emin = floor(log10(min(xpos)));
    emax = ceil(log10(max(xpos)));
    if emax - emin < 1
        emin = emin - 1;
        emax = emax + 1;
    end
    exps = emin:emax;

    set(ax, 'XScale', 'log');
    xlim(ax, [10^emin, 10^emax]);
    set(ax, 'XTick', 10.^exps);
    set(ax, 'XTickLabel', arrayfun(@(e) sprintf('$-10^{%d}$', e), exps, ...
        'UniformOutput', false), 'TickLabelInterpreter', 'latex');
    set(ax, 'XDir', 'reverse');
else
    xlim(ax, [min(mu_values) max(mu_values)]);
end
end


function ok = draw_density_panel(ax, scores, color, xLimits, fsTextBox)
% Draw a kernel-density estimate of 'scores' on xLimits, scaled to unit
% height, with the sample mean (dashed line) and a small box with n, mean
% and standard deviation. Returns false (and draws nothing) if the density
% cannot be computed.
ok = false;

xgrid = linspace(xLimits(1), xLimits(2), 400);
try
    [f, xi] = ksdensity(scores, xgrid, 'Support', 'unbounded');
catch
    return
end

if isempty(f) || all(~isfinite(f)) || max(f,[],'omitnan') <= 0
    return
end

f = f / max(f,[],'omitnan');

fill(ax, [xi fliplr(xi)], [f zeros(size(f))], color, ...
    'FaceAlpha', 0.6, 'EdgeColor', 'none');
plot(ax, xi, f, 'k-', 'LineWidth', 0.6);

m  = mean(scores, 'omitnan');
sd = std(scores, 'omitnan');
plot(ax, [m m], [0 1], '--', 'Color', 0.15*[1 1 1]);

text(ax, 0.02, 0.96, sprintf('n=%d\nm=%.3g\nsd=%.3g', numel(scores), m, sd), ...
    'Units','normalized','VerticalAlignment','top', ...
    'BackgroundColor','w','EdgeColor','k','FontSize',fsTextBox,'Margin',1);

xlim(ax, xLimits);
ylim(ax, [0 1.05]);
ok = true;
end


function draw_no_samples(ax, fsTick)
% Placeholder for a panel without enough data.
text(ax, 0.5, 0.5, 'no samples', 'Units','normalized', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fsTick);
xlim(ax, [0 1]);
ylim(ax, [0 1]);
xticks(ax, []);
yticks(ax, []);
box(ax,'on');
end


function save_figure(fig, jpgPath, alsoSaveFig)
% Export a figure as a 300-dpi JPEG (and optionally as a .fig file with the
% same name). A failed export gives a warning rather than an error.
try
    exportgraphics(fig, jpgPath, 'Resolution', 300);
    if alsoSaveFig
        [folder, name] = fileparts(jpgPath);
        savefig(fig, fullfile(folder, [name '.fig']));
    end
catch ME
    warning('Could not save figure %s (%s)', jpgPath, ME.message);
end
end
