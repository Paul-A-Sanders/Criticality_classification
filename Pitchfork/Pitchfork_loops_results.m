% Pitchfork_loops_results.m
% =========================================================================
% PURPOSE
%   Test how well the fitted cubic coefficient s_hat of data from a
%   stochastic pitchfork normal form separates the two criticality classes
%   (s = +1 subcritical, s = -1 supercritical):
%
%       dq = (mu*q + s*q^3) dtau,      mu < 0 .
%
%   Trajectories are simulated for a sweep of mu values and both classes,
%   the coefficients [mu_hat; s_hat] are estimated from the first tau_max
%   time units of each trajectory (for several tau_max), and the
%   classification performance of s_hat is summarised with ROC curves, AUC
%   values and optimal ROC thresholds.
%
% WORKFLOW
%   1. If a compatible processed MAT-file already exists, load it and skip
%      the simulations.
%   2. Otherwise, for every (mu, s): simulate independent noise
%      realisations in parallel (parfor), reject diverging or failed runs,
%      and fit the model with analyze_pitchfork to find the coefficients. 
%      If too few realisations survive at the longest observation horizon, 
%      the horizon is reduced step by step until enough do.
%   3. Collect the estimates, compute AUC for every (mu, tau_max), save.
%   4. Draw figures: AUC vs mu, distributions of the estimates, ROC curves,
%      optimal ROC threshold vs mu.
%
% REQUIREMENTS
%   * simulate_pitchfork.m and analyze_pitchfork.m on the MATLAB path.
%   * Parallel Computing Toolbox (parfor) and Statistics and Machine
%     Learning Toolbox (perfcurve, ksdensity).
%   * MATLAB R2020a or newer (exportgraphics, orderedcolors, tiledlayout,
%     yline); developed with R2026b.
%
% OUTPUTS
%   A folder named after the key settings is created in the current
%   directory:
%       Pitchfork_ROC_nr=..._tmin=..._tmax=..._nt=..._dt=.../
%           saved_data/           processed MAT-file (see DATA LAYOUT)
%           auc_plots/            AUC vs mu, one line per tau_max
%           ridge_plots_all/      distributions of mu_hat and s_hat
%           roc_plots_per_time/   ROC curves, one figure per tau_max
%           threshold_plots/      optimal ROC threshold vs mu
%
% DATA LAYOUT (saved in the MAT-file)
%   per_theta{comp, iMu, iS, it}        row vector of fitted values for
%                                       component comp (1 = mu_hat,
%                                       2 = s_hat) at mu_values(iMu),
%                                       s_values(iS), tmax_values(it)
%   per_theta_seeds{comp, iMu, iS, it}  RNG seeds of those realisations
%   auc_all_all(iMu, it)                AUC of s_hat for the two classes
%   roc_opt_fpr, roc_opt_tpr, roc_opt_threshold (nMu x nT)
%                                       optimal ROC operating point and the
%                                       s_hat threshold that attains it
%
% REUSING SAVED DATA
%   If the MAT-file for these settings already exists, it is loaded and the
%   simulations are skipped. The values stored in that file (mu_values,
%   tmax_values, s_values, n_realizations, ...) then REPLACE the ones set
%   below. To re-simulate after changing a setting that is not part of the
%   folder name (for example mu_values), set rerunSimulations = true or
%   delete the saved file.
%
% Review the user settings below before running.
% =========================================================================

clearvars; clc;

numCores = feature('numCores');   % Used to size the parallel batches.

%% ========================================================================
%% 1. USER SETTINGS
%% ========================================================================

% --- Model and noise ---
noiseType    = "White";   % 'White', 'Uniform' or 'None' (see simulate_pitchfork.m)
InitialScale = 0;         % Std. dev. of the random initial condition q(0); 0 starts at q = 0.
sigma        = 1;         % Noise amplitude.

% --- Parameter sweep ---
mu_values = -logspace(-1,2,50);   % mu from -0.1 to -100 (log-spaced). Values
                                  % closer to 0 are closer to the bifurcation.
s_values  = [1, -1];              % True sign of the cubic term: the two classes to distinguish.
nS        = numel(s_values);

% --- Ensemble and observation horizons ---
n_realizations        = 100;      % Target number of realisations per (mu, s, tau_max).
u_explosion_threshold = 1e10;     % Reject trajectories whose |q| exceeds this value.

tmax_values = [1,2,5,10,20,50,100,200,500,1000,2000,5000];   % Observation horizons tau_max.
                                  % Each must be an integer multiple of dt_target; at most 12
                                  % are allowed (figure colours come from orderedcolors("gem12")).
nT = numel(tmax_values);

dt_target = 0.01;                 % Euler-Maruyama time step.

% --- Classification score ---
kc_shat    = 2;   % theta = [mu_hat; s_hat]; component 2 (s_hat) is the ROC score.
nThetaComp = 2;   % Number of fitted components that are kept.

% --- Run control ---
rerunSimulations = true;   % true: ignore any saved MAT-file and simulate again

nMu = numel(mu_values);

%% ========================================================================
%% 2. OUTPUT LOCATIONS
%% ========================================================================
% The key settings are part of the folder name, so different settings never
% overwrite each other.
tmin     = min(tmax_values);
tmax_val = max(tmax_values);
base_dir_root = fullfile(pwd, sprintf( ...
    'Pitchfork_ROC_nr=%d_tmin=%.3g_tmax=%.3g_nt=%d_dt=%.3g', ...
    n_realizations, tmin, tmax_val, nT, dt_target));
if ~exist(base_dir_root,'dir'), mkdir(base_dir_root); end

outDataDir = fullfile(base_dir_root, 'saved_data');
saveName   = sprintf('Pitchfork_ROC_data_nr%d_dt%.3g.mat', n_realizations, dt_target);
savePath   = fullfile(outDataDir, saveName);

%% ========================================================================
%% 3. LOAD PROCESSED DATA (IF AVAILABLE)
%% ========================================================================
% A complete saved file lets the script skip the simulation and
% post-processing stages. Only the processed arrays needed for the figures
% are required.
loadedFromFile = false;
if exist(savePath, 'file') && ~rerunSimulations
    fprintf('Found saved data at %s. Loading and skipping simulations...\n', savePath);
    saved = load(savePath);

    requiredFields = {'per_theta','per_theta_seeds','mu_values','tmax_values', ...
        's_values','kc_shat','nThetaComp','n_realizations','dt_target', ...
        'auc_all_all'};
    missingFields = requiredFields(~isfield(saved, requiredFields));

    if isempty(missingFields)
        per_theta       = saved.per_theta;
        per_theta_seeds = saved.per_theta_seeds;
        mu_values       = saved.mu_values;
        tmax_values     = saved.tmax_values;
        s_values        = saved.s_values;
        kc_shat         = saved.kc_shat;
        nThetaComp      = saved.nThetaComp;
        n_realizations  = saved.n_realizations;
        dt_target       = saved.dt_target;
        auc_all_all     = saved.auc_all_all;

        % base_dir_root is deliberately NOT restored from the file: figures
        % are written next to the file found by this run, so the script
        % still works if the folder has been moved.
        nMu = numel(mu_values);
        nT  = numel(tmax_values);
        nS  = numel(s_values);
        loadedFromFile = true;

        fprintf('Loaded processed data successfully. Proceeding directly to plotting.\n');
    else
        warning(['Saved file exists but is missing required fields (%s). ' ...
            'The simulations will be run again.'], strjoin(missingFields, ', '));
        clear saved
    end
end

%% ========================================================================
%% 4. SIMULATION AND PARAMETER ESTIMATION
%% ========================================================================
% Skipped when processed data were loaded above.
if ~loadedFromFile

    per_theta       = cell(nThetaComp, nMu, nS, nT);   % See DATA LAYOUT in the header.
    per_theta_seeds = cell(nThetaComp, nMu, nS, nT);

    successfulCells = cell(nMu, nS);             % Accepted realisations per (mu, s)
    usedItPerCell   = zeros(nMu, nS, 'uint8');   % Longest tau_max index reached per (mu, s)

    for iMu = 1:nMu
        mu = mu_values(iMu);
        fprintf('mu (%d/%d) = %g\n', iMu, nMu, mu);

        for iS = 1:nS
            % Both classes are run at the same mu so that their fitted
            % coefficient distributions can be compared directly.
            s_true = s_values(iS);
            fprintf('\ns = %g\n', s_true);

            % Start every (mu, s) at the longest requested horizon. A
            % reduction needed for one class therefore does not carry over
            % to the other. activeIt indexes tmax_values.
            activeIt = nT;
            tmax_max = tmax_values(activeIt);
            [tsteps_max, tidx_ref] = horizonGrid(tmax_values, tmax_max, dt_target);

            batchSize = n_realizations;   % Upper limit on the number of seeds per batch
            seedStart = 0;
            batch     = 0;

            % 'successful' pools the accepted realisations. Each entry has
            %   seed        RNG seed of the realisation
            %   maxIt       index of the longest tau_max with a valid fit
            %   theta       cell(1,nT) of fitted coefficient vectors
            %   validTimes  logical(1,nT) marking horizons with a valid fit
            successful = struct('valid',false,'seed',[],'maxIt',0,'theta',{{}},'validTimes',false(1,nT));
            successful = successful([]);   % Start with an empty struct array

            % There is no fixed cap on the number of seed attempts: new
            % batches are generated until n_realizations realisations are
            % available at the active horizon (or the horizon is reduced).
            while true
                collected = sum([successful.maxIt] >= activeIt);
                if collected >= n_realizations, break; end

                batch = batch + 1;
                fprintf('\nBatch %d\n', batch);

                numSeedsThisBatch = nextBatchSize(n_realizations - collected, numCores, batchSize);
                seeds     = seedStart + (1:numSeedsThisBatch);
                seedStart = seedStart + numSeedsThisBatch;

                results = struct('valid',false,'seed',0,'maxIt',0,'theta',{{}},'validTimes',false(1,nT));
                results = repmat(results, numel(seeds), 1);
                fprintf('Trying seeds %d -- %d\n', seeds(1), seeds(end));

                % ---- Parallel simulation and fitting, one seed per iteration ----
                parfor ib = 1:numel(seeds)
                    seed = seeds(ib);
                    rng(seed,'combRecursive');   % Same seed = same independent realisation

                    theta_store = cell(1,nT);
                    validTimes  = false(1,nT);
                    local_maxIt = 0;

                    try
                        % One stochastic trajectory for this (mu, s, seed). It is simulated
                        % once up to the active horizon and then truncated for every tau_max.
                        [t_vec, u_vec] = simulate_pitchfork(mu, s_true, tmax_max, tsteps_max, noiseType, sigma, InitialScale);
                    catch
                        % Any error skips this seed silently. If EVERY seed fails
                        % (e.g. an invalid noiseType), the horizon-reduction logic
                        % below ends in "Cannot reduce tmax further".
                        continue
                    end

                    if isempty(u_vec)
                        continue
                    end

                    % Make sure the time dimension runs along the rows.
                    if numel(t_vec) == size(u_vec,1)
                        u_sim = u_vec;
                        t_sim = t_vec;
                    elseif numel(t_vec) == size(u_vec,2)
                        u_sim = u_vec';
                        t_sim = t_vec;
                    else
                        continue
                    end

                    % Fit the first tau_max time units for every horizon.
                    for it = 1:nT
                        idx_end = tidx_ref(it);
                        if isnan(idx_end) || idx_end < 3, continue; end

                        t_sub = t_sim(1:idx_end);
                        u_sub = u_sim(1:idx_end);

                        % Reject non-finite or exploding trajectory segments.
                        if any(isnan(u_sub(:))) || any(isinf(u_sub(:))) || max(abs(u_sub(:))) > u_explosion_threshold
                            continue
                        end

                        try
                            % Fit the normal form; theta = [mu_hat; s_hat].
                            [theta_sub, ~, ~] = analyze_pitchfork(u_sub, t_sub);
                        catch
                            continue
                        end

                        theta_sub = reshape(theta_sub(:),[],1);
                        if numel(theta_sub) < nThetaComp
                            theta_sub(end+1:nThetaComp,1) = NaN;
                        end

                        theta_store{it} = theta_sub(1:nThetaComp);
                        validTimes(it)  = true;
                        local_maxIt     = max(local_maxIt, it);
                    end

                    if any(validTimes)
                        results(ib).valid      = true;
                        results(ib).seed       = seed;
                        results(ib).maxIt      = local_maxIt;
                        results(ib).theta      = theta_store;
                        results(ib).validTimes = validTimes;
                    end
                end % parfor

                % ---- Merge new valid results into the pool (skip duplicate seeds) ----
                % This is an extra protection line for if the seeds are
                % chosen differently so that duplicates exist, currently it
                % is not relevant
                existingSeeds = [successful.seed];
                for ib = 1:numel(results)
                    if ~results(ib).valid, continue; end
                    if ismember(results(ib).seed, existingSeeds), continue; end
                    successful(end+1)    = results(ib); %#ok<SAGROW>
                    existingSeeds(end+1) = results(ib).seed; %#ok<SAGROW>
                end

                % ---- Progress report ----
                collected = sum([successful.maxIt] >= activeIt);
                fprintf('Status: iMu=%d/%d, iS=%d/%d, batch=%d, activeIt=%d, collected=%d/%d\n', ...
                    iMu, nMu, iS, nS, batch, activeIt, collected, n_realizations);

                if ~isempty(successful)
                    % Number of accepted realisations by the longest horizon they reached.
                    counts = histcounts([successful.maxIt], 0.5:1:(nT+0.5));
                    fprintf('  successful maxIt counts (it=1..%d): %s\n', nT, sprintf('%d ', counts));

                    seedList = unique([successful.seed]);
                    nShow = min(10, numel(seedList));
                    fprintf('  sample kept seeds (%d): %s\n', numel(seedList), sprintf('%d ', seedList(1:nShow)));
                end

                % ---- Adaptive reduction of the observation horizon ----
                % If fewer than 10% of the requested realisations reach the
                % active horizon, shorten it by one step. This stops a few
                % diverging or invalid trajectories from blocking the sweep.
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

            % Remember the accepted realisations and the final horizon for this (mu, s).
            successfulCells{iMu, iS} = successful;
            usedItPerCell(iMu, iS)   = uint8(activeIt);
        end
    end
end

%% ========================================================================
%% 5. ASSEMBLE ESTIMATES, COMPUTE AUC, SAVE
%% ========================================================================
if ~loadedFromFile

    % ---- Build per_theta from the accepted realisations ----
    fprintf('\nBuilding per_theta from successful simulations...\n');

    for iMu = 1:nMu
        for iS = 1:nS
            S = successfulCells{iMu,iS};
            if isempty(S), continue; end

            % Keep realisations that reached the final active horizon and
            % use exactly n_realizations of them, in seed order.
            finalActiveIt = double(usedItPerCell(iMu, iS));
            S = S([S.maxIt] >= finalActiveIt);
            if numel(S) > n_realizations, S = S(1:n_realizations); end

            for sidx = 1:numel(S)
                sStruct = S(sidx);
                for it = 1:nT
                    if isempty(sStruct.theta), continue; end
                    if numel(sStruct.theta) < it, continue; end
                    theta_it = sStruct.theta{it};
                    if isempty(theta_it), continue; end
                    for kc = 1:nThetaComp
                        if kc <= numel(theta_it) && ~isnan(theta_it(kc))
                            per_theta{kc, iMu, iS, it}(end+1)       = theta_it(kc);
                            per_theta_seeds{kc, iMu, iS, it}(end+1) = sStruct.seed;
                        end
                    end
                end
            end
        end
    end

    % ---- Classification performance: AUC of s_hat ----
    % For each (mu, tau_max) the fitted s_hat values of the two true classes
    % are pooled and passed to perfcurve. AUC = 1 means perfect separation;
    % AUC = 0.5 is chance level. The positive class is s > 0, and larger
    % s_hat is taken to favour it.
    fprintf('\nComputing AUC...\n');

    auc_all_all = nan(nMu, nT);

    for iMu = 1:nMu
        for it = 1:nT
            scores_all = [];
            labels_all = [];
            completeClasses = true;

            for iS = 1:nS
                vals = per_theta{kc_shat, iMu, iS, it};
                vals = vals(isfinite(vals));

                % Use a (mu, s, tau_max) dataset only if it holds at least
                % n_realizations values; other points are unaffected.
                if numel(vals) < n_realizations
                    completeClasses = false;
                    break
                end

                % Equal class sizes: keep exactly n_realizations values.
                vals = vals(1:n_realizations);
                vals = vals(:)';

                scores_all = [scores_all, vals]; %#ok<AGROW>
                labels_all = [labels_all, ...
                    repmat(s_values(iS) > 0, 1, n_realizations)]; %#ok<AGROW>
            end

            % AUC needs both classes to be complete at this (mu, tau_max).
            if ~completeClasses
                continue
            end

            [~,~,~,auc_all_all(iMu,it)] = perfcurve(labels_all, scores_all, 1);
        end
    end

    % ---- Save processed results ----
    % Only compact fitted quantities are stored; the trajectories are not
    % needed because every figure is generated from these estimates.
    if ~exist(outDataDir,'dir'), mkdir(outDataDir); end

    save(savePath, ...
        'per_theta', ...
        'per_theta_seeds', ...
        'mu_values', ...
        'tmax_values', ...
        's_values', ...
        'kc_shat', ...
        'nThetaComp', ...
        'n_realizations', ...
        'dt_target', ...
        'auc_all_all', ...
        'base_dir_root', ...
        '-v7');

    fprintf('Saved results to: %s\n', savePath);

else
    fprintf('Using saved per_theta/AUC data; postprocessing and saving are skipped.\n');
end

%% ========================================================================
%% 6. FIGURE STYLE
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
%% 7. FIGURE: AUC vs mu (one line per tau_max)
%% ========================================================================
fig = figure('Name','AUC vs mu - all tmax','Color','w','Units','normalized', ...
    'Position',[0.15 0.15 0.7 0.6]);
ax = axes(fig); hold(ax,'on');
styleAxes(ax,'Box','on');

colorsT = orderedcolors("gem12");

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

for it = 1:nT
    auc_vec = auc_all_all(:,it);
    if all(isnan(auc_vec)), continue; end

    color = colorsT(it,:);
    valid = isfinite(auc_vec) & isfinite(xplot(:));
    if ~any(valid), continue; end

    if useLog
        hPlot(it) = semilogx(ax,xplot,auc_vec,'Color',color,'LineWidth',1.8);
    else
        hPlot(it) = plot(ax,xplot,auc_vec,'Color',color,'LineWidth',1.8);
    end

    % Dots mark the mu values that were actually simulated.
    scatter(ax,xplot(valid),auc_vec(valid),36,'o', ...
        'MarkerFaceColor',color, ...
        'MarkerEdgeColor',color*0.6, ...
        'MarkerFaceAlpha',0.9);

    labels{it} = sprintf('$\\tau_{\\mathrm{max}}=%.4g$',tmax_values(it));
end

% Keep only the lines that were actually drawn for the legend.
validIdx = arrayfun(@isgraphics,hPlot);
hPlot  = hPlot(validIdx);
labels = labels(validIdx);

if useLog
    setLogMuAxis(ax,xpos);
else
    xlim(ax,[min(mu_values) max(mu_values)]);
end

xlabel(ax,'$\mu$','FontSize',fsLabel,'Interpreter','latex');
ylabel(ax,'AUC','FontSize',fsLabel);
title(ax,'AUC vs $\mu$','FontSize',fsFigureTitle,'Interpreter','latex');
ylim(ax,[0 1]);
set(ax,'Layer','top');

if ~isempty(hPlot)
    lg = legend(ax,hPlot,labels,'Location','best', ...
        'Interpreter','latex','FontSize',fsLegend);
    lg.Box = 'on';
    lg.NumColumns = min(3,numel(hPlot));
end

grid(ax,'off');
hold(ax,'off');
drawnow;

outDirAUC = fullfile(base_dir_root,'auc_plots');
if ~exist(outDirAUC,'dir'), mkdir(outDirAUC); end

fnBase = sprintf('AUC_vs_mu_nmu=%d_tmin=%g_tmax=%g.jpg', ...
    numel(mu_values),min(tmax_values),max(tmax_values));
fnPath = fullfile(outDirAUC,fnBase);

try
    exportgraphics(fig,fnPath,'Resolution',300);
    fprintf('Saved AUC figure to: %s\n',fnPath);
catch
    % Fallback if exportgraphics fails: try saveas, otherwise just warn.
    try
        saveas(fig,fullfile(outDirAUC,[fnBase(1:end-4) '.fig']));
        saveas(fig,fullfile(outDirAUC,[fnBase(1:end-4) '.jpg']));
    catch
        fprintf(2,'WARNING: Could not save AUC figure.\n');
    end
end

%% ========================================================================
%% 8. FIGURES: distributions (ridge panels) of mu_hat and s_hat
%% ========================================================================
% One figure per component and tau_max. Rows = a few mu values spread over
% the sweep, columns = the two true classes. Each panel shows the kernel
% density estimate (scaled to peak 1) of the fitted coefficient, its mean
% (dashed line) and a box with n, mean and standard deviation. Together
% these show how much the two classes overlap, complementing the AUC.
%
% The x axis is a FIXED window, the same for both classes and every tau_max,
% so panels can be compared directly (as in the Swift-Hohenberg scripts):
%     mu_hat : true mu +/- 10  (centred on the value being estimated)
%     s_hat  : [-30, 30]       (the true values are s = +1 and s = -1)
% The y axis is fixed to [0, 1.05]. Values outside the window are not shown,
% so very wide distributions (typically short tau_max) are cut off at the
% edges; change xLimits below if you need a different window.
kc_list = unique([1,kc_shat]);
minForKDE = 3;                          % Minimum samples for a density estimate

cmapPos = [0 0.4470 0.7410];            % s > 0 (blue)
cmapNeg = [0.8500 0.3250 0.0980];       % s < 0 (orange)

maxMuPlot   = 5;                        % Show at most 5 mu values (rows)
idx_mu_plot = unique(round(linspace(1,nMu,min(nMu,maxMuPlot))));
nMu_plot    = numel(idx_mu_plot);

outDir = fullfile(base_dir_root,'ridge_plots_all');
if ~exist(outDir,'dir'), mkdir(outDir); end

for kc = kc_list
    if kc == 1
        compLabel = '$\hat{\mu}$';
        plotTag   = 'mu_hat';
    else
        compLabel = '$\hat{s}$';
        plotTag   = 's_hat';
    end

    for it_plot = 1:nT

        figRidge = figure('Name',sprintf('Ridge %s t=%g', ...
            plotTag,tmax_values(it_plot)), ...
            'Color','w','Units','normalized', ...
            'Position',[0.08 0.07 0.82 0.82]);

        t = tiledlayout(figRidge,nMu_plot,nS, ...
            'TileSpacing','compact','Padding','compact');

        for iMu_plot = 1:nMu_plot
            iMu = idx_mu_plot(iMu_plot);
            mu  = mu_values(iMu);

            % Fixed x window for this row (identical for both classes).
            if kc == 1
                xLimits = mu + [-10 10];   % Centred on the true value of mu_hat (= mu)
            else
                xLimits = [-30 30];
            end
            xgrid = linspace(xLimits(1),xLimits(2),400);

            leftHasSamples = false;   % Reset for each row (see y tick labels below)

            for iS = 1:nS
                ax = nexttile(t,(iMu_plot-1)*nS+iS);
                cla(ax); hold(ax,'on');
                styleAxes(ax,'TickLabelInterpreter','latex');

                scores = per_theta{kc,iMu,iS,it_plot};
                scores = scores(isfinite(scores));

                % A panel is drawn only if this individual (mu, s, tau_max)
                % dataset is complete; other datasets do not matter.
                hasData = numel(scores) >= max(n_realizations,minForKDE);
                if hasData
                    scores = scores(1:n_realizations);
                end

                % Column header: the true s.
                if iMu_plot == 1
                    title(ax,sprintf('$s=%g$',s_values(iS)), ...
                        'Interpreter','latex','FontSize',fsPanelTitle);
                end

                % y tick labels are kept on the left; the right panel hides
                % them only if the left panel has data.
                if iS == 1
                    leftHasSamples = hasData;
                end

                % Row label: mu, written vertically on the outer right.
                if iS == nS
                    addRowLabel(ax,mu,fsPanelTitle);
                end

                % Column label under the bottom row: the fitted coefficient.
                if iMu_plot == nMu_plot
                    xlabel(ax,compLabel, ...
                        'Interpreter','latex','FontSize',fsLabel);
                end

                if ~hasData
                    showNoSamples(ax,fsTick);
                    continue
                end

                scores = scores(:);

                % Kernel density estimate on the fixed window (no histograms).
                try
                    [f,xi] = ksdensity(scores,xgrid,'Support','unbounded');
                catch
                    f = [];
                end

                if isempty(f) || all(~isfinite(f)) || max(f,[],'omitnan') <= 0
                    % The density could not be evaluated on this window (for
                    % example because all values lie far outside it).
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

                f = f/max(f,[],'omitnan');   % Normalise the peak to 1

                color = cmapPos;
                if s_values(iS) < 0
                    color = cmapNeg;
                end

                fill(ax,[xi fliplr(xi)],[f zeros(size(f))],color, ...
                    'FaceAlpha',0.6,'EdgeColor','none');
                plot(ax,xi,f,'k-','LineWidth',0.6);

                % Mean (dashed line) and summary statistics.
                m     = mean(scores,'omitnan');
                sd    = std(scores,'omitnan');
                nSamp = numel(scores);

                plot(ax,[m m],[0 1],'--','Color',0.15*[1 1 1]);

                ridgeText = sprintf('n=%d\nm=%.3g\nsd=%.3g',nSamp,m,sd);
                text(ax,0.02,0.96,ridgeText, ...
                    'Units','normalized', ...
                    'VerticalAlignment','top', ...
                    'BackgroundColor','w', ...
                    'EdgeColor','k', ...
                    'FontSize',fsTextBox, ...
                    'Margin',1);

                xlim(ax,xLimits);
                ylim(ax,[0 1.05]);

                % Hide the numerical y tick labels of the right panel (the
                % rows share their scale); ticks and grid stay.
                if iS ~= 1 && leftHasSamples
                    yticklabels(ax,[]);
                end

                grid(ax,'on');
                box(ax,'on');
                hold(ax,'off');
            end
        end

        ylabel(t,'Normalized density','FontSize',fsLabel);

        titleStr = sprintf('%s distributions $(\\tau_{max}=%.4g)$', ...
            compLabel,tmax_values(it_plot));
        title(t,titleStr,'Interpreter','latex','FontSize',fsFigureTitle);

        fname = fullfile(outDir,sprintf('ridge_%s_t=%g.jpg', ...
            plotTag,tmax_values(it_plot)));
        exportJpg(figRidge,fname);
    end
end

%% ========================================================================
%% 9. FIGURES: ROC curves per observation horizon + optimal thresholds
%% ========================================================================
% One ROC figure per tau_max, with curves for a representative set of mu
% values. In the same loop the optimal ROC operating point (and the s_hat
% threshold that attains it) is stored for EVERY valid (mu, tau_max), not
% only for the mu values drawn here; it is shown in section 10.
kc_plot = kc_shat;

outDirROC = fullfile(base_dir_root,'roc_plots_per_time');
if ~exist(outDirROC,'dir'), mkdir(outDirROC); end

% Representative mu values: the grid points closest (in log10|mu|) to ten
% log-spaced targets between -100 and -0.1.
nROCcurves    = 10;
target_mu_roc = -logspace(2,-1,nROCcurves);

idx_mu_roc = zeros(1,nROCcurves);
for j = 1:nROCcurves
    [~,idx_mu_roc(j)] = min(abs( ...
        log10(-mu_values)-log10(-target_mu_roc(j))));
end
idx_mu_roc = unique(idx_mu_roc,'stable');

% Rows correspond to mu, columns to tau_max.
roc_opt_fpr       = nan(nMu,nT);
roc_opt_tpr       = nan(nMu,nT);
roc_opt_threshold = nan(nMu,nT);

for it_plot = 1:nT

    % ---- Pool scores and labels of both classes for every mu ----
    all_scores_by_mu = cell(nMu,1);
    all_labels_by_mu = cell(nMu,1);

    for iMu = 1:nMu
        scores = [];
        labels = [];
        completeClasses = true;

        for iS = 1:nS
            vals = per_theta{kc_plot,iMu,iS,it_plot};
            vals = vals(isfinite(vals));

            % Check only this individual (mu, s, tau_max) dataset.
            if numel(vals) < n_realizations
                completeClasses = false;
                break
            end

            vals = vals(1:n_realizations);
            vals = vals(:);

            scores = [scores;vals]; %#ok<AGROW>
            labels = [labels; ...
                repmat(s_values(iS)>0,n_realizations,1)]; %#ok<AGROW>
        end

        if completeClasses
            all_scores_by_mu{iMu} = scores;
            all_labels_by_mu{iMu} = labels;
        end
        % Otherwise leave empty: only this mu is skipped at this tau_max.
    end

    % ---- Optimal ROC point and threshold for ALL valid mu values ----
    for iMu = 1:nMu
        scores = all_scores_by_mu{iMu};
        labels = all_labels_by_mu{iMu};

        if isempty(scores) || ...
                numel(scores) ~= nS*n_realizations || ...
                numel(unique(labels)) < 2
            continue
        end

        [fpr_all,tpr_all,thresholds_all,~,optROCPoint_all] = ...
            perfcurve(labels,scores,1);

        % Find the threshold belonging to MATLAB's optimal ROC point. The
        % nearest-point formulation is robust against tiny floating-point
        % differences between OPTROCPT and the returned ROC coordinates.
        [~,idxOpt_all] = min((fpr_all-optROCPoint_all(1)).^2 + ...
                             (tpr_all-optROCPoint_all(2)).^2);

        optThreshold_all = thresholds_all(idxOpt_all);

        roc_opt_fpr(iMu,it_plot) = optROCPoint_all(1);
        roc_opt_tpr(iMu,it_plot) = optROCPoint_all(2);

        % Ignore +/-Inf thresholds, which can occur at the ROC end points.
        if isfinite(optThreshold_all)
            roc_opt_threshold(iMu,it_plot) = optThreshold_all;
        end
    end

    % ---- ROC figure for this tau_max ----
    figROC = figure('Name',sprintf('ROC t=%g',tmax_values(it_plot)), ...
        'Color','w','Units','normalized', ...
        'Position',[0.18 0.18 0.52 0.58]);

    axROC = axes(figROC); hold(axROC,'on');
    styleAxes(axROC,'TickLabelInterpreter','latex');

    plot(axROC,[0 1],[0 1],'k--','LineWidth',1, ...
        'HandleVisibility','off');            % Chance line

    colors = orderedcolors("gem12");
    hasROC = false;

    for j = 1:numel(idx_mu_roc)
        iMu = idx_mu_roc(j);

        scores = all_scores_by_mu{iMu};
        labels = all_labels_by_mu{iMu};

        if isempty(scores) || ...
                numel(scores) ~= nS*n_realizations || ...
                numel(unique(labels)) < 2
            continue
        end

        % perfcurve also returns MATLAB's optimal ROC operating point as
        % [FPR, TPR]. With the balanced classes used here and equal
        % misclassification costs this maximises TPR - FPR (Youden's J).
        [fpr,tpr,~,auc,optROCPoint] = perfcurve(labels,scores,1);

        plot(axROC,fpr,tpr, ...
            'Color',colors(j,:), ...
            'LineWidth',1.6, ...
            'DisplayName',sprintf('\\mu=%.3g (AUC=%.2f)', ...
            mu_values(iMu),auc));

        % Mark the optimal operating point on the curve.
        plot(axROC,optROCPoint(1),optROCPoint(2),'o', ...
            'MarkerSize',7, ...
            'MarkerFaceColor',colors(j,:), ...
            'MarkerEdgeColor','k', ...
            'LineWidth',0.8, ...
            'HandleVisibility','off');

        hasROC = true;
    end

    xlabel(axROC,'False Positive Rate','FontSize',fsLabel);
    ylabel(axROC,'True Positive Rate','FontSize',fsLabel);

    title(axROC,sprintf('ROC $(\\tau_{max}=%.4g)$', ...
        tmax_values(it_plot)), ...
        'Interpreter','latex','FontSize',fsFigureTitle);

    axis(axROC,[0 1 0 1]);
    grid(axROC,'on');
    box(axROC,'on');

    if hasROC
        lgROC = legend(axROC,'Location','southeast', ...
            'Interpreter','latex','FontSize',fsLegend+1);
        lgROC.Box = 'on';
    end

    drawnow;

    fn = fullfile(outDirROC,sprintf('ROC_kc%d_t%g.jpg', ...
        kc_plot,tmax_values(it_plot)));
    exportJpg(figROC,fn);
end

%% ========================================================================
%% 10. FIGURE: optimal ROC threshold vs mu (all observation horizons)
%% ========================================================================
% Value of s_hat at the optimal ROC operating point for each mu and tau_max.
% A threshold of 0 (dashed line) corresponds to classifying by the sign of
% the fitted cubic coefficient.
figThresh = figure('Name','Optimal ROC threshold vs mu - all tmax', ...
    'Color','w','Units','normalized','Position',[0.15 0.15 0.7 0.6]);
axThresh = axes(figThresh); hold(axThresh,'on');
styleAxes(axThresh,'Box','on');

colorsT = orderedcolors("gem12");

% As in the AUC figure: plot |mu| on a log axis with negative tick labels.
xpos   = -mu_values;
useLog = (max(xpos)/min(xpos) > 10);
if useLog
    xplot = xpos;
else
    xplot = mu_values;
end

hThresh      = gobjects(nT,1);
labelsThresh = cell(nT,1);

% Only thresholds inside [yMin, yMax] are shown. Values outside are set to
% NaN, which breaks the line instead of drawing artificial connections
% across large excursions, while valid points can reappear closer to mu = 0.
yMin = -2;
yMax = 2;

for it = 1:nT
    threshold_vec = roc_opt_threshold(:,it);

    valid   = isfinite(threshold_vec) & isfinite(xplot(:));
    inRange = valid & threshold_vec >= yMin & threshold_vec <= yMax;

    if ~any(inRange)
        continue
    end

    yPlot = threshold_vec;
    yPlot(~inRange) = NaN;

    color = colorsT(it,:);

    if useLog
        hThresh(it) = semilogx(axThresh, xplot, yPlot, ...
            'Color', color, 'LineWidth', 1.8);
    else
        hThresh(it) = plot(axThresh, xplot, yPlot, ...
            'Color', color, 'LineWidth', 1.8);
    end

    % Markers only at the accepted points.
    scatter(axThresh, xplot(inRange), threshold_vec(inRange), 36, 'o', ...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', color*0.6, ...
        'MarkerFaceAlpha', 0.9, ...
        'HandleVisibility', 'off');

    labelsThresh{it} = sprintf('$\\tau_{\\mathrm{max}}=%.4g$', ...
        tmax_values(it));
end

ylim(axThresh,[yMin yMax]);

% Reference line: threshold = 0.
yline(axThresh,0,'k--','LineWidth',1.1,'HandleVisibility','off');

validIdx = arrayfun(@isgraphics,hThresh);
hThresh      = hThresh(validIdx);
labelsThresh = labelsThresh(validIdx);

if useLog
    setLogMuAxis(axThresh,xpos);
else
    xlim(axThresh,[min(mu_values) max(mu_values)]);
end

xlabel(axThresh,'$\mu$','FontSize',fsLabel,'Interpreter','latex');
ylabel(axThresh,'Optimal threshold for $\hat{s}$', ...
    'FontSize',fsLabel,'Interpreter','latex');
title(axThresh,'Optimal ROC threshold vs $\mu$', ...
    'FontSize',fsFigureTitle,'Interpreter','latex');

set(axThresh,'Layer','top');
grid(axThresh,'off');

if ~isempty(hThresh)
    lgThresh = legend(axThresh,hThresh,labelsThresh,'Location','best', ...
        'Interpreter','latex','FontSize',fsLegend);
    lgThresh.Box = 'on';
    lgThresh.NumColumns = min(3,numel(hThresh));
end

hold(axThresh,'off');
drawnow;

outDirThresh = fullfile(base_dir_root,'threshold_plots');
if ~exist(outDirThresh,'dir'), mkdir(outDirThresh); end

fnThresh = fullfile(outDirThresh,sprintf( ...
    'optimal_threshold_vs_mu_nmu=%d_tmin=%g_tmax=%g.jpg', ...
    numel(mu_values),min(tmax_values),max(tmax_values)));

if exportJpg(figThresh,fnThresh)
    fprintf('Saved optimal-threshold figure to: %s\n',fnThresh);
end

% Append the operating-point results to the saved MAT-file. This works both
% after a fresh simulation and when per_theta was loaded from the file.
try
    save(savePath,'roc_opt_fpr','roc_opt_tpr','roc_opt_threshold','-append');
    fprintf('Appended ROC operating points and thresholds to: %s\n',savePath);
catch
    warning('Could not append ROC operating-point data');
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

function ok = exportJpg(fig, fname)
% EXPORTJPG  Save a figure at 300 dpi; warn (instead of erroring) on failure.
    ok = true;
    try
        exportgraphics(fig,fname,'Resolution',300);
    catch
        ok = false;
        warning('Could not save %s',fname);
    end
end
