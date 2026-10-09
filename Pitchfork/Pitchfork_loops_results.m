% test_code_pitchfork_loops.m
% -------------------------------------------------------------------------
% PURPOSE
% Simulate pitchfork trajectories over a sweep of negative mu values and
% both s classes. Estimate the model parameters at several observation
% horizons, then assess classification performance using ROC curves, AUC,
% and optimal ROC thresholds.
%
% WORKFLOW
% 1. Load previously processed results when a compatible MAT-file exists.
% 2. Otherwise, simulate trajectories, analyse their parameter estimates,
%    and adaptively reduce the active observation horizon when needed.
% 3. Assemble the estimates, calculate AUC values, save the results, and
%    generate AUC, density, ROC, and threshold figures.
%
% REQUIREMENTS
% simulate_pitchfork.m and analyze_pitchfork.m must be on the MATLAB path.
% This script also uses MATLAB parallel-computing and statistics functions,
% including parfor, perfcurve, ksdensity, and exportgraphics.
%
% OUTPUTS
% A settings-specific directory is created in the current working directory.
% It contains the processed MAT-file and subfolders for generated figures.
%
% Review the user settings below before running. The code below retains the
% original calculations and plotting commands; explanatory comments have
% been added to make the workflow easier to follow.
% -------------------------------------------------------------------------

clearvars; clc;

numCores = feature('numCores'); % Used to choose parallel batch sizes.

%% ------------------------
%% User settings
%% ------------------------
noiseType = "White";
InitialScale = 0;
sigma = 1;

mu_values = -logspace(-1,2,50);    % Negative mu sweep.
% mu is the distance-to-bifurcation parameter; values closer to 0 are
% closer to the pitchfork and therefore provide the strongest early-warning regime.
n_realizations = 30;               % Target number of runs per s class.
u_explosion_threshold = 1e10;      % Reject trajectories that grow too large.

tmax_values = [1,2,5,10,20,50,100,200,500,1000,2000,5000]; % Observation horizons.
dt_target = 0.01;                   % Requested simulation time step.
nT = numel(tmax_values);

s_values = [1, -1];                % The two classes to distinguish.
% Here s is the true sign of the cubic normal-form coefficient. The ROC analysis
% asks whether its fitted value can distinguish the two criticality classes.
nS = numel(s_values);

% Select the theta component used as the score for ROC analysis.
kc_shat = 2;
% The second fitted coefficient is used as the classification score (s-hat).
nThetaComp = 2;

%% ------------------------
%% Output paths
%% ------------------------
% Include key settings in the output folder name to distinguish runs.
tmin = min(tmax_values);
tmax_val = max(tmax_values);
base_dir_root = fullfile(pwd, sprintf( ...
    'Pitchfork_ROC_nr=%d_tmin=%.3g_tmax=%.3g_nt=%d_dt=%.3g', ...
    n_realizations, tmin, tmax_val, nT, dt_target));
if ~exist(base_dir_root,'dir'), mkdir(base_dir_root); end

outDataDir = fullfile(base_dir_root, 'saved_data');
saveName = sprintf('Pitchfork_ROC_data_nr%d_dt%.3g.mat', n_realizations, dt_target);
savePath = fullfile(outDataDir, saveName);

%% ------------------------
%% Load processed data if available
%% ------------------------
% A complete saved file lets the script skip the simulation and postprocessing
% stages. The required fields are the processed arrays needed by the plots.
loadedFromFile = false;
if exist(savePath, 'file')
    fprintf('Found saved data at %s. Loading and skipping simulations...\n', savePath);
    saved = load(savePath);

    % These are the processed quantities needed to reproduce all plots.
    % Do NOT require successfulCells: those are only needed while building
    % per_theta after a fresh simulation.
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

        % Keep the CURRENT base_dir_root. This makes the script portable if
        % the folder was moved; figures are written next to the MAT-file that
        % was found by this run rather than to an old saved absolute path.
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

%% ------------------------
%% Simulation and estimation
%% ------------------------
if ~loadedFromFile
    % A fresh simulation is only performed when no compatible processed data
    % file was found above. All expensive simulation/fitting work is therefore
    % separated from the later plotting workflow.
    %% ------------------------
    %% Preallocations
    %% ------------------------
    nMu = numel(mu_values);

    % aggregated AUCs: (mu x t)
    auc_all_all = nan(nMu, nT);

    % per_theta{comp, iMu, iS, it} (cells for later plotting)
    per_theta = cell(nThetaComp, nMu, nS, nT);
    per_theta_seeds = cell(nThetaComp, nMu, nS, nT);

    %% ------------------------
    %% Main loops (synchronous, simple)
    %% ------------------------
    for iMu = 1:nMu
        % Repeat the experiment for each distance from the bifurcation.
        mu = mu_values(iMu);
        fprintf('mu (%d/%d) = %g\n', iMu, nMu, mu);

        for iS = 1:nS
            % Run both cubic-sign classes at the same mu so that their fitted
            % coefficient distributions can be compared directly.
            s_true = s_values(iS);
            fprintf('\ns = %g\n', s_true);

            % Reset the maximum observation horizon independently for each s.
            % A reduction required for one criticality class therefore does not
            % carry over to the other class.
            activeIt = nT;
            tmax_max = tmax_values(activeIt);
            tsteps_max = round(tmax_max/dt_target) + 1;
            dt_actual = tmax_max/(tsteps_max-1);
            if abs(dt_actual-dt_target) > 1e-12*max(1,abs(dt_target))
                error('tmax must be an integer multiple of dt_target.');
            end
            tvec_ref = (0:tsteps_max-1)'*dt_target;

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

            batchSize = n_realizations;

            % No fixed maximum number of seed attempts. Continue generating
            % batches until n_realizations have been collected at the current
            % active observation horizon, unless the adaptive-horizon logic
            % reduces activeIt.
            seedStart = 0;
            batch = 0;

            successful = struct('valid',false,'seed',[],'maxIt',0,'theta',{{}},'validTimes',false(1,nT));
            successful = successful([]);  % empty

            while true
                collected = sum([successful.maxIt] >= activeIt);
                if collected >= n_realizations, break; end

                batch = batch + 1;
                fprintf('\nBatch %d\n',batch);

                numToGo = n_realizations - collected;
                overProvision = 1.25;
                tentative = ceil(overProvision * numToGo);
                tentative = max(numCores, tentative);
                tentative = min(batchSize, tentative);
                numSeedsThisBatch = ceil(tentative/numCores)*numCores;
                numSeedsThisBatch = min(numSeedsThisBatch, batchSize);

                seeds = seedStart + (1:numSeedsThisBatch);
                seedStart = seedStart + numSeedsThisBatch;

                results = struct('valid',false,'seed',0,'maxIt',0,'theta',{{}},'validTimes',false(1,nT));
                results = repmat(results, numel(seeds), 1);
                fprintf('Trying seeds %d -- %d\n', seeds(1), seeds(end));

                parfor ib = 1:numel(seeds)
                    seed = seeds(ib);
                    rng(seed,'combRecursive');

                    theta_store = cell(1,nT);
                    validTimes = false(1,nT);
                    local_maxIt = 0;

                    try
                        % Simulate one stochastic trajectory for this (mu,s,seed) combination.
                        % The same seed identifies one independent realization across
                        % all observation horizons extracted from this trajectory.
                        [t_vec, u_vec] = simulate_pitchfork(mu, s_true, tmax_max, tsteps_max, noiseType, sigma, InitialScale);
                    catch
                        continue
                    end

                    if isempty(u_vec)
                        continue
                    end

                    % ensure time dimension corresponds to rows
                    if numel(t_vec) == size(u_vec,1)
                        u_sim = u_vec;
                        t_sim = t_vec;
                    elseif numel(t_vec) == size(u_vec,2)
                        u_sim = u_vec';
                        t_sim = t_vec;
                    else
                        continue
                    end

                    for it = 1:nT
                        idx_end = tidx_ref(it);
                        if isnan(idx_end) || idx_end < 3, continue; end

                        t_sub = t_sim(1:idx_end);
                        u_sub = u_sim(1:idx_end);

                        if any(isnan(u_sub(:))) || any(isinf(u_sub(:))) || max(abs(u_sub(:))) > u_explosion_threshold
                            continue
                        end

                        try
                            % Fit the normal-form model to the trajectory segment. The resulting theta
                            % contains the estimated coefficients used for classification.
                            [theta_sub, ~, ~] = analyze_pitchfork(u_sub, t_sub);
                        catch
                            continue
                        end

                        theta_sub = reshape(theta_sub(:),[],1);
                        if numel(theta_sub) < nThetaComp
                            theta_sub(end+1:nThetaComp,1) = NaN;
                        end

                        theta_store{it} = theta_sub(1:nThetaComp);
                        validTimes(it) = true;
                        local_maxIt = max(local_maxIt, it);
                    end

                    if any(validTimes)
                        results(ib).valid = true;
                        results(ib).seed = seed;
                        results(ib).maxIt = local_maxIt;
                        results(ib).theta = theta_store;
                        results(ib).validTimes = validTimes;
                    end
                end % parfor

                existingSeeds = [successful.seed];
                for ib = 1:numel(results)
                    if ~results(ib).valid, continue; end
                    if ismember(results(ib).seed, existingSeeds), continue; end
                    successful(end+1) = results(ib); %#ok<SAGROW>
                    existingSeeds(end+1) = results(ib).seed; %#ok<SAGROW>
                end

                collected = sum([successful.maxIt] >= activeIt);
                % --- Insert after: collected = sum([successful.maxIt] >= activeIt);
                % Run-status summary for debugging / monitoring
                fprintf('Status: iMu=%d/%d, iS=%d/%d, batch=%d, activeIt=%d, collected=%d/%d\n', ...
                    iMu, nMu, iS, nS, batch, activeIt, collected, n_realizations);

                % show distribution of successful simulations by achieved maxIt
                if ~isempty(successful)
                    maxIts = [successful.maxIt];
                    if ~isempty(maxIts)
                        % counts for each it level (1..nT)
                        counts = histcounts(maxIts, 0.5:1:(nT+0.5));
                        countsStr = sprintf('%d ', counts);
                        fprintf('  successful maxIt counts (it=1..%d): %s\n', nT, countsStr);
                        % show a few kept seed IDs (up to 10)
                        seedList = [successful.seed];
                        seedList = unique(seedList);
                        nShow = min(10, numel(seedList));
                        fprintf('  sample kept seeds (%d): %s\n', numel(seedList), sprintf('%d ', seedList(1:nShow)));
                    end
                end

                % If fewer than 10% of the requested realizations reach the current
                % horizon, the horizon is shortened. This prevents a small number
                % of divergent/invalid trajectories from stopping the entire sweep.
                % adaptive reduction
                minAccept = ceil(0.1 * n_realizations);
                if collected < minAccept
                    if activeIt > 1
                        activeIt = activeIt - 1;
                        % before/after sizes
                        beforeN = numel(successful);
                        successful = successful([successful.maxIt] >= activeIt);
                        afterN = numel(successful);
                        fprintf('Reducing active tmax -> %.3g. kept seeds: %d -> %d (activeIt=%d)\n', ...
                            tmax_max, beforeN, afterN, activeIt);
                        tmax_max = tmax_values(activeIt);
                        fprintf('Reducing active tmax -> %.3g\n', tmax_max);
                        tsteps_max = round(tmax_max/dt_target) + 1;
                        dt_actual = tmax_max/(tsteps_max-1);
                        if abs(dt_actual-dt_target) > 1e-12*max(1,abs(dt_target))
                            error('tmax must be an integer multiple of dt_target.');
                        end
                        tvec_ref = (0:tsteps_max-1)'*dt_target;

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
                    else
                        warning('Cannot reduce tmax further. Stopping.');
                        break
                    end
                end
            end % while

            % store successfulSeeds for this (iMu,iS)
            if ~exist('successfulCells','var') || isempty(successfulCells)
                successfulCells = cell(nMu, nS);
            end
            successfulCells{iMu, iS} = successful;

            if ~exist('usedItPerCell','var') || isempty(usedItPerCell)
                usedItPerCell = zeros(nMu, nS, 'uint8');
            end
            usedItPerCell(iMu, iS) = uint8(activeIt);
        end
    end
end

if ~loadedFromFile
    %% ============================================================
    %% POSTPROCESSING: build per_theta from successful simulations (no L)
    %% ============================================================
    fprintf('\nBuilding per_theta from successful simulations...\n');

for kc = 1:nThetaComp
    for it = 1:nT
        for iMu = 1:nMu
            for iS = 1:nS
                per_theta{kc, iMu, iS, it} = [];
                per_theta_seeds{kc, iMu, iS, it} = [];
            end
        end
    end
end

for iMu = 1:nMu
    for iS = 1:nS
        S = {};
        if exist('successfulCells','var') && ~isempty(successfulCells{iMu,iS})
            S = successfulCells{iMu,iS};
        end
        if isempty(S), continue; end

        finalActiveIt = double(usedItPerCell(iMu, iS));
        keepMaskSeeds = [S.maxIt] >= finalActiveIt;
        S = S(keepMaskSeeds);
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
                        per_theta{kc, iMu, iS, it} = [per_theta{kc, iMu, iS, it}, theta_it(kc)]; %#ok<SAGROW>
                        per_theta_seeds{kc, iMu, iS, it} = [per_theta_seeds{kc, iMu, iS, it}, sStruct.seed]; %#ok<SAGROW>
                    end
                end
            end
        end
    end
end

%% ============================================================
%% CLASSIFICATION PERFORMANCE
% For every (mu,tau_max), the fitted s-hat values from the two true classes
% are combined and passed to perfcurve. AUC=1 means perfect separation,
% whereas AUC=0.5 corresponds to chance-level discrimination.
%% AUC computation (no L)
%% ============================================================
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

            % Use this (mu,s,tau_max) dataset only when it contains at least
            % n_realizations samples. Other times and mu values are unaffected.
            if numel(vals) < n_realizations
                completeClasses = false;
                break
            end

            % Keep the class sizes equal if more than n_realizations exist.
            vals = vals(1:n_realizations);
            vals = vals(:)';

            scores_all = [scores_all, vals]; %#ok<AGROW>
            labels_all = [labels_all, ...
                repmat(s_values(iS) > 0, 1, n_realizations)]; %#ok<AGROW>
        end

        % AUC is only computed if both s classes are complete for this
        % particular (mu,tau_max) point.
        if ~completeClasses
            continue
        end

        % The positive ROC class is s>0. The score itself is the fitted cubic
        % coefficient, so larger values should favour the positive class.
        [~,~,~,auc_all_all(iMu,it)] = ...
            perfcurve(labels_all, scores_all, 1);
    end
end

%% ============================================================
%% SAVE PROCESSED RESULTS
% Only compact fitted quantities are saved. The original trajectories do not
% need to be retained because all subsequent figures are generated from these
% parameter estimates and their associated seeds.
%% SAVE RESULTS
%% ============================================================
outDataDir = fullfile(base_dir_root, 'saved_data');
if ~exist(outDataDir,'dir')
    mkdir(outDataDir);
end

saveName = sprintf('Pitchfork_ROC_data_nr%d_dt%.3g.mat', ...
    n_realizations, dt_target);
savePath = fullfile(outDataDir, saveName);

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

%% ============================================================
%% FIGURE SETTINGS
% The calculations above produce the scientific results; everything below is
% presentation code for figures used to inspect and communicate those results.
%% Paper-ready figure style (Below is all just code and settings for plotting)
%% ============================================================
paperFont = 'Helvetica';

fsTick        = 13;
fsLabel       = 15;
fsPanelTitle  = 14;
fsFigureTitle = 17;
fsLegend      = 12;
fsTextBox     = 10;

axisLineWidth = 1.1;

%% ============================================================
%% AUC vs mu
%% ============================================================
fig = figure('Name','AUC vs mu - all tmax','Color','w','Units','normalized', ...
    'Position',[0.15 0.15 0.7 0.6]);
ax = axes(fig); hold(ax,'on');

set(ax,'FontName',paperFont,'FontSize',fsTick,'Box','on', ...
    'LineWidth',axisLineWidth,'TickDir','in');

colorsT = orderedcolors("gem12");

% MATLAB logarithmic axes require positive coordinates. Plot |mu| and use
% negative tick labels because all mu values in this sweep are negative.
xpos = -mu_values;
useLog = (max(xpos)/min(xpos) > 10);

hPlot = gobjects(nT,1);
labels = cell(nT,1);

firstBad_idx = nan(nT,1);
firstBad_mu = nan(nT,1);

for it = 1:nT
    auc_vec = squeeze(auc_all_all(:,it));
    if all(isnan(auc_vec)), continue; end

    color = colorsT(it,:);

    if useLog
        xplot = xpos;
    else
        xplot = mu_values;
    end

    valid = isfinite(auc_vec) & isfinite(xplot(:));

    lastNanIdx = find(isnan(auc_vec),1,'last');
    if ~isempty(lastNanIdx)
        firstBad_idx(it) = lastNanIdx;
        firstBad_mu(it) = mu_values(lastNanIdx);
    end

    if ~any(valid), continue; end

    if useLog
        hPlot(it) = semilogx(ax,xplot,auc_vec, ...
            'Color',color,'LineWidth',1.8);
    else
        hPlot(it) = plot(ax,xplot,auc_vec, ...
            'Color',color,'LineWidth',1.8);
    end

    scatter(ax,xplot(valid),auc_vec(valid),36,'o', ...
        'MarkerFaceColor',color, ...
        'MarkerEdgeColor',color*0.6, ...
        'MarkerFaceAlpha',0.9);

    labels{it} = sprintf('$\\tau_{\\mathrm{max}}=%.4g$', ...
        tmax_values(it));
end

validIdx = arrayfun(@isgraphics,hPlot);
hPlot = hPlot(validIdx);
labels = labels(validIdx);

if useLog
    set(ax,'XScale','log');

    xmin = min(xpos);
    xmax = max(xpos);
    emin = floor(log10(xmin));
    emax = ceil(log10(xmax));

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
    try
        saveas(fig,fullfile(outDirAUC,[fnBase(1:end-4) '.fig']));
        saveas(fig,fullfile(outDirAUC,[fnBase(1:end-4) '.jpg']));
    catch ME
        % test_code_pitchfork_loops.m
clearvars; clc;

numCores = feature('numCores');
%% ------------------------
%% User settings
%% ------------------------
noiseType = "White";
InitialScale = 0;
sigma = 1;

mu_values = -logspace(-1,2,50);    % negative mu sweep
n_realizations = 100;
u_explosion_threshold = 1e10;

tmax_values = [1,2,5,10,20,50,100,200,500,1000,2000,5000];
dt_target = 0.01;
nT = numel(tmax_values);

s_values = [1, -1];
nS = numel(s_values);

% which theta component is used as s_hat for ROC
kc_shat = 2;
nThetaComp = 2;

% output folders: include t-range, number of tmax values and timestamp
tmin = min(tmax_values);
tmax_val = max(tmax_values);
base_dir_root = fullfile(pwd, sprintf( ...
    'Pitchfork_ROC_nr=%d_tmin=%.3g_tmax=%.3g_nt=%d_dt=%.3g', ...
    n_realizations, tmin, tmax_val, nT, dt_target));
if ~exist(base_dir_root,'dir'), mkdir(base_dir_root); end

outDataDir = fullfile(base_dir_root, 'saved_data');
saveName = sprintf('Pitchfork_ROC_data_nr%d_dt%.3g.mat', n_realizations, dt_target);
savePath = fullfile(outDataDir, saveName);

loadedFromFile = false;
if exist(savePath, 'file')
    fprintf('Found saved data at %s. Loading and skipping simulations...\n', savePath);
    saved = load(savePath);

    % These are the processed quantities needed to reproduce all plots.
    % Do NOT require successfulCells: those are only needed while building
    % per_theta after a fresh simulation.
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

        % Keep the CURRENT base_dir_root. This makes the script portable if
        % the folder was moved; figures are written next to the MAT-file that
        % was found by this run rather than to an old saved absolute path.
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

if ~loadedFromFile
    %% ------------------------
    %% Preallocations
    %% ------------------------
    nMu = numel(mu_values);

    % aggregated AUCs: (mu x t)
    auc_all_all = nan(nMu, nT);

    % per_theta{comp, iMu, iS, it} (cells for later plotting)
    per_theta = cell(nThetaComp, nMu, nS, nT);
    per_theta_seeds = cell(nThetaComp, nMu, nS, nT);

    %% ------------------------
    %% Main loops (synchronous, simple)
    %% ------------------------
    for iMu = 1:nMu
        mu = mu_values(iMu);
        fprintf('mu (%d/%d) = %g\n', iMu, nMu, mu);

        for iS = 1:nS
            s_true = s_values(iS);
            fprintf('\ns = %g\n', s_true);

            % Reset the maximum observation horizon independently for each s.
            % A reduction required for one criticality class therefore does not
            % carry over to the other class.
            activeIt = nT;
            tmax_max = tmax_values(activeIt);
            tsteps_max = round(tmax_max/dt_target) + 1;
            dt_actual = tmax_max/(tsteps_max-1);
            if abs(dt_actual-dt_target) > 1e-12*max(1,abs(dt_target))
                error('tmax must be an integer multiple of dt_target.');
            end
            tvec_ref = (0:tsteps_max-1)'*dt_target;

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

            batchSize = n_realizations;

            % No fixed maximum number of seed attempts. Continue generating
            % batches until n_realizations have been collected at the current
            % active observation horizon, unless the adaptive-horizon logic
            % reduces activeIt.
            seedStart = 0;
            batch = 0;

            successful = struct('valid',false,'seed',[],'maxIt',0,'theta',{{}},'validTimes',false(1,nT));
            successful = successful([]);  % empty

            while true
                collected = sum([successful.maxIt] >= activeIt);
                if collected >= n_realizations, break; end

                batch = batch + 1;
                fprintf('\nBatch %d\n',batch);

                numToGo = n_realizations - collected;
                overProvision = 1.25;
                tentative = ceil(overProvision * numToGo);
                tentative = max(numCores, tentative);
                tentative = min(batchSize, tentative);
                numSeedsThisBatch = ceil(tentative/numCores)*numCores;
                numSeedsThisBatch = min(numSeedsThisBatch, batchSize);

                seeds = seedStart + (1:numSeedsThisBatch);
                seedStart = seedStart + numSeedsThisBatch;

                results = struct('valid',false,'seed',0,'maxIt',0,'theta',{{}},'validTimes',false(1,nT));
                results = repmat(results, numel(seeds), 1);
                fprintf('Trying seeds %d -- %d\n', seeds(1), seeds(end));

                parfor ib = 1:numel(seeds)
                    seed = seeds(ib);
                    rng(seed,'combRecursive');

                    theta_store = cell(1,nT);
                    validTimes = false(1,nT);
                    local_maxIt = 0;

                    try
                        [t_vec, u_vec] = simulate_pitchfork(mu, s_true, tmax_max, tsteps_max, noiseType, sigma, InitialScale);
                    catch
                        continue
                    end

                    if isempty(u_vec)
                        continue
                    end

                    % ensure time dimension corresponds to rows
                    if numel(t_vec) == size(u_vec,1)
                        u_sim = u_vec;
                        t_sim = t_vec;
                    elseif numel(t_vec) == size(u_vec,2)
                        u_sim = u_vec';
                        t_sim = t_vec;
                    else
                        continue
                    end

                    for it = 1:nT
                        idx_end = tidx_ref(it);
                        if isnan(idx_end) || idx_end < 3, continue; end

                        t_sub = t_sim(1:idx_end);
                        u_sub = u_sim(1:idx_end);

                        if any(isnan(u_sub(:))) || any(isinf(u_sub(:))) || max(abs(u_sub(:))) > u_explosion_threshold
                            continue
                        end

                        try
                            [theta_sub, ~, ~] = analyze_pitchfork(u_sub, t_sub);
                        catch
                            continue
                        end

                        theta_sub = reshape(theta_sub(:),[],1);
                        if numel(theta_sub) < nThetaComp
                            theta_sub(end+1:nThetaComp,1) = NaN;
                        end

                        theta_store{it} = theta_sub(1:nThetaComp);
                        validTimes(it) = true;
                        local_maxIt = max(local_maxIt, it);
                    end

                    if any(validTimes)
                        results(ib).valid = true;
                        results(ib).seed = seed;
                        results(ib).maxIt = local_maxIt;
                        results(ib).theta = theta_store;
                        results(ib).validTimes = validTimes;
                    end
                end % parfor

                existingSeeds = [successful.seed];
                for ib = 1:numel(results)
                    if ~results(ib).valid, continue; end
                    if ismember(results(ib).seed, existingSeeds), continue; end
                    successful(end+1) = results(ib); %#ok<SAGROW>
                    existingSeeds(end+1) = results(ib).seed; %#ok<SAGROW>
                end

                collected = sum([successful.maxIt] >= activeIt);
                % --- Insert after: collected = sum([successful.maxIt] >= activeIt);
                % Run-status summary for debugging / monitoring
                fprintf('Status: iMu=%d/%d, iS=%d/%d, batch=%d, activeIt=%d, collected=%d/%d\n', ...
                    iMu, nMu, iS, nS, batch, activeIt, collected, n_realizations);

                % show distribution of successful simulations by achieved maxIt
                if ~isempty(successful)
                    maxIts = [successful.maxIt];
                    if ~isempty(maxIts)
                        % counts for each it level (1..nT)
                        counts = histcounts(maxIts, 0.5:1:(nT+0.5));
                        countsStr = sprintf('%d ', counts);
                        fprintf('  successful maxIt counts (it=1..%d): %s\n', nT, countsStr);
                        % show a few kept seed IDs (up to 10)
                        seedList = [successful.seed];
                        seedList = unique(seedList);
                        nShow = min(10, numel(seedList));
                        fprintf('  sample kept seeds (%d): %s\n', numel(seedList), sprintf('%d ', seedList(1:nShow)));
                    end
                end

                % adaptive reduction
                minAccept = ceil(0.1 * n_realizations);
                if collected < minAccept
                    if activeIt > 1
                        activeIt = activeIt - 1;
                        % before/after sizes
                        beforeN = numel(successful);
                        successful = successful([successful.maxIt] >= activeIt);
                        afterN = numel(successful);
                        fprintf('Reducing active tmax -> %.3g. kept seeds: %d -> %d (activeIt=%d)\n', ...
                            tmax_max, beforeN, afterN, activeIt);
                        tmax_max = tmax_values(activeIt);
                        fprintf('Reducing active tmax -> %.3g\n', tmax_max);
                        tsteps_max = round(tmax_max/dt_target) + 1;
                        dt_actual = tmax_max/(tsteps_max-1);
                        if abs(dt_actual-dt_target) > 1e-12*max(1,abs(dt_target))
                            error('tmax must be an integer multiple of dt_target.');
                        end
                        tvec_ref = (0:tsteps_max-1)'*dt_target;

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
                    else
                        warning('Cannot reduce tmax further. Stopping.');
                        break
                    end
                end
            end % while

            % store successfulSeeds for this (iMu,iS)
            if ~exist('successfulCells','var') || isempty(successfulCells)
                successfulCells = cell(nMu, nS);
            end
            successfulCells{iMu, iS} = successful;

            if ~exist('usedItPerCell','var') || isempty(usedItPerCell)
                usedItPerCell = zeros(nMu, nS, 'uint8');
            end
            usedItPerCell(iMu, iS) = uint8(activeIt);
        end
    end
end

if ~loadedFromFile
    %% ============================================================
    %% POSTPROCESSING: build per_theta from successful simulations (no L)
    %% ============================================================
    fprintf('\nBuilding per_theta from successful simulations...\n');

for kc = 1:nThetaComp
    for it = 1:nT
        for iMu = 1:nMu
            for iS = 1:nS
                per_theta{kc, iMu, iS, it} = [];
                per_theta_seeds{kc, iMu, iS, it} = [];
            end
        end
    end
end

for iMu = 1:nMu
    for iS = 1:nS
        S = {};
        if exist('successfulCells','var') && ~isempty(successfulCells{iMu,iS})
            S = successfulCells{iMu,iS};
        end
        if isempty(S), continue; end

        finalActiveIt = double(usedItPerCell(iMu, iS));
        keepMaskSeeds = [S.maxIt] >= finalActiveIt;
        S = S(keepMaskSeeds);
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
                        per_theta{kc, iMu, iS, it} = [per_theta{kc, iMu, iS, it}, theta_it(kc)]; %#ok<SAGROW>
                        per_theta_seeds{kc, iMu, iS, it} = [per_theta_seeds{kc, iMu, iS, it}, sStruct.seed]; %#ok<SAGROW>
                    end
                end
            end
        end
    end
end

%% ============================================================
%% AUC computation (no L)
%% ============================================================
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

            % Use this (mu,s,tau_max) dataset only when it contains at least
            % n_realizations samples. Other times and mu values are unaffected.
            if numel(vals) < n_realizations
                completeClasses = false;
                break
            end

            % Keep the class sizes equal if more than n_realizations exist.
            vals = vals(1:n_realizations);
            vals = vals(:)';

            scores_all = [scores_all, vals]; %#ok<AGROW>
            labels_all = [labels_all, ...
                repmat(s_values(iS) > 0, 1, n_realizations)]; %#ok<AGROW>
        end

        % AUC is only computed if both s classes are complete for this
        % particular (mu,tau_max) point.
        if ~completeClasses
            continue
        end

        [~,~,~,auc_all_all(iMu,it)] = ...
            perfcurve(labels_all, scores_all, 1);
    end
end

%% ============================================================
%% SAVE RESULTS
%% ============================================================
outDataDir = fullfile(base_dir_root, 'saved_data');
if ~exist(outDataDir,'dir')
    mkdir(outDataDir);
end

saveName = sprintf('Pitchfork_ROC_data_nr%d_dt%.3g.mat', ...
    n_realizations, dt_target);
savePath = fullfile(outDataDir, saveName);

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

%% ============================================================
%% Paper-ready figure style
%% ============================================================
paperFont = 'Helvetica';

fsTick        = 13;
fsLabel       = 15;
fsPanelTitle  = 14;
fsFigureTitle = 17;
fsLegend      = 12;
fsTextBox     = 10;

axisLineWidth = 1.1;

%% ============================================================
%% AUC vs mu
%% ============================================================
fig = figure('Name','AUC vs mu - all tmax','Color','w','Units','normalized', ...
    'Position',[0.15 0.15 0.7 0.6]);
ax = axes(fig); hold(ax,'on');

set(ax,'FontName',paperFont,'FontSize',fsTick,'Box','on', ...
    'LineWidth',axisLineWidth,'TickDir','in');

colorsT = orderedcolors("gem12");

% MATLAB logarithmic axes require positive coordinates. Plot |mu| and use
% negative tick labels because all mu values in this sweep are negative.
xpos = -mu_values;
useLog = (max(xpos)/min(xpos) > 10);

hPlot = gobjects(nT,1);
labels = cell(nT,1);

firstBad_idx = nan(nT,1);
firstBad_mu = nan(nT,1);

for it = 1:nT
    auc_vec = squeeze(auc_all_all(:,it));
    if all(isnan(auc_vec)), continue; end

    color = colorsT(it,:);

    if useLog
        xplot = xpos;
    else
        xplot = mu_values;
    end

    valid = isfinite(auc_vec) & isfinite(xplot(:));

    lastNanIdx = find(isnan(auc_vec),1,'last');
    if ~isempty(lastNanIdx)
        firstBad_idx(it) = lastNanIdx;
        firstBad_mu(it) = mu_values(lastNanIdx);
    end

    if ~any(valid), continue; end

    if useLog
        hPlot(it) = semilogx(ax,xplot,auc_vec, ...
            'Color',color,'LineWidth',1.8);
    else
        hPlot(it) = plot(ax,xplot,auc_vec, ...
            'Color',color,'LineWidth',1.8);
    end

    scatter(ax,xplot(valid),auc_vec(valid),36,'o', ...
        'MarkerFaceColor',color, ...
        'MarkerEdgeColor',color*0.6, ...
        'MarkerFaceAlpha',0.9);

    labels{it} = sprintf('$\\tau_{\\mathrm{max}}=%.4g$', ...
        tmax_values(it));
end

validIdx = arrayfun(@isgraphics,hPlot);
hPlot = hPlot(validIdx);
labels = labels(validIdx);

if useLog
    set(ax,'XScale','log');

    xmin = min(xpos);
    xmax = max(xpos);
    emin = floor(log10(xmin));
    emax = ceil(log10(xmax));

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
    try
        saveas(fig,fullfile(outDirAUC,[fnBase(1:end-4) '.fig']));
        saveas(fig,fullfile(outDirAUC,[fnBase(1:end-4) '.jpg']));
    catch
        fprintf(2,'WARNING: Could not save AUC figure.\n');
    end
end

%% ============================================================
%% Ridge / density panels
%% ============================================================
kc_list = unique([1,kc_shat]);
minForKDE = 3;

cmapPos = [0 0.4470 0.7410];
cmapNeg = [0.8500 0.3250 0.0980];

for kc = kc_list
    for it_plot = 1:nT

        maxMuPlot = 5;
        nMu_total = numel(mu_values);
        nMu_plot = min(nMu_total,maxMuPlot);

        idx_mu_plot = unique(round(linspace(1,nMu_total,nMu_plot)));
        nMu_plot = numel(idx_mu_plot);

        if kc == 1
            compLabel = '$\hat{\mu}$';
            plotTag = 'mu_hat';
        else
            compLabel = '$\hat{s}$';
            plotTag = 's_hat';
        end

        figRidge = figure('Name',sprintf('Ridge %s t=%g', ...
            plotTag,tmax_values(it_plot)), ...
            'Color','w','Units','normalized', ...
            'Position',[0.08 0.07 0.82 0.82]);

        t = tiledlayout(figRidge,nMu_plot,nS, ...
            'TileSpacing','compact','Padding','compact');

        for iMu_plot = 1:nMu_plot
            iMu = idx_mu_plot(iMu_plot);
            mu = mu_values(iMu);

            for iS = 1:nS
                ax = nexttile(t,(iMu_plot-1)*nS+iS);
                cla(ax); hold(ax,'on');

                set(ax,'FontName',paperFont,'FontSize',fsTick, ...
                    'LineWidth',axisLineWidth,'TickDir','in', ...
                    'TickLabelInterpreter','latex');

                scores = per_theta{kc,iMu,iS,it_plot};
                scores = scores(isfinite(scores));

                % Check only this individual (mu,s,tau_max) dataset.
                % An incomplete dataset at another time does not affect it.
                if numel(scores) >= n_realizations
                    scores = scores(1:n_realizations);
                end

                % s as column header.
                if iMu_plot == 1
                    title(ax,sprintf('$s=%g$',s_values(iS)), ...
                        'Interpreter','latex','FontSize',fsPanelTitle);
                end

                % Numerical density labels only on the left.
                if iS ~= 1
                    yticklabels(ax,[]);
                end

                % mu as vertical row label on the outer right.
                if iS == nS
                    text(ax,1.025,0.5, ...
                        sprintf('$\\mu=%.1f$',mu), ...
                        'Units','normalized', ...
                        'Interpreter','latex', ...
                        'FontSize',fsPanelTitle, ...
                        'Rotation',90, ...
                        'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle', ...
                        'Clipping','off');
                end

                % Fitted coefficient label below each column.
                if iMu_plot == nMu_plot
                    xlabel(ax,compLabel, ...
                        'Interpreter','latex','FontSize',fsLabel);
                end

                if numel(scores) < n_realizations
                    text(ax,0.5,0.5,'no samples', ...
                        'Units','normalized', ...
                        'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle', ...
                        'FontSize',fsTick);
                    xlim(ax,[0 1]);
                    ylim(ax,[0 1]);
                    xticks(ax,[]);
                    yticks(ax,[]);
                    box(ax,'on');
                    hold(ax,'off');
                    continue
                end

                scores = scores(:);

                if numel(scores) >= minForKDE
                    pr = prctile(scores,[1 99]);

                    if pr(1) == pr(2)
                        pr = pr+[-0.5 0.5];
                    end

                   prRange = pr(2) - pr(1);

                   xgrid = linspace( ...
                            pr(1) - 0.1*prRange, ...
                            pr(2) + 0.1*prRange, ...
                            400);

                                            [f,xi] = ksdensity(scores,xgrid,'Support','unbounded');

                    if max(f) > 0
                        f = f/max(f);
                    end

                    color = cmapPos;
                    if s_values(iS) < 0
                        color = cmapNeg;
                    end

                    fill(ax,[xi fliplr(xi)],[f zeros(size(f))],color, ...
                        'FaceAlpha',0.6,'EdgeColor','none');
                    plot(ax,xi,f,'k-','LineWidth',0.6);
                    xlim(ax,[min(xi) max(xi)]);
                else
                    color = (s_values(iS)>0)*cmapPos + ...
                        (s_values(iS)<0)*cmapNeg;
                    histogram(ax,scores,'Normalization','pdf', ...
                        'FaceAlpha',0.6,'FaceColor',color);
                end

                m = mean(scores,'omitnan');
                sd = std(scores,'omitnan');
                nSamp = numel(scores);

                yl = ylim(ax);
                if isempty(yl) || any(~isfinite(yl))
                    yl = [0 1];
                    ylim(ax,yl);
                end

                plot(ax,[m m],[0.02*yl(2) 0.98*yl(2)],'--', ...
                    'Color',0.15*[1 1 1]);

                ridgeText = sprintf('n=%d\nm=%.3g\nsd=%.3g', ...
                    nSamp,m,sd);

                text(ax,0.02,0.95,ridgeText, ...
                    'Units','normalized', ...
                    'VerticalAlignment','top', ...
                    'BackgroundColor','w', ...
                    'EdgeColor','k', ...
                    'FontSize',fsTextBox, ...
                    'Margin',1);

                if iS ~= 1
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

        outDir = fullfile(base_dir_root,'ridge_plots_all');
        if ~exist(outDir,'dir'), mkdir(outDir); end

        fname = fullfile(outDir,sprintf('ridge_%s_t=%g.jpg', ...
            plotTag,tmax_values(it_plot)));

        try
            exportgraphics(figRidge,fname,'Resolution',300);
        catch
            fprintf(2,'WARNING: Could not save ridge figure to %s\n',fname);
        end
    end
end

%% ============================================================
%% ROC curves per observation horizon
%% ============================================================
kc_plot = kc_shat;

outDirROC = fullfile(base_dir_root,'roc_plots_per_time');
if ~exist(outDirROC,'dir'), mkdir(outDirROC); end

% Match the SH presentation: show a representative set of mu values.
nROCcurves = 10;
target_mu_roc = -logspace(2,-1,nROCcurves);

idx_mu_roc = zeros(1,nROCcurves);
for j = 1:nROCcurves
    [~,idx_mu_roc(j)] = min(abs( ...
        log10(-mu_values)-log10(-target_mu_roc(j))));
end
idx_mu_roc = unique(idx_mu_roc,'stable');

% Store optimal operating-point information for EVERY valid (mu,tau_max)
% pair, not only for the representative mu values shown in the ROC figures.
% Rows correspond to mu, columns to tau_max.
roc_opt_fpr       = nan(nMu,nT);
roc_opt_tpr       = nan(nMu,nT);
roc_opt_threshold = nan(nMu,nT);

for it_plot = 1:nT
    all_scores_by_mu = cell(nMu,1);
    all_labels_by_mu = cell(nMu,1);

    for iMu = 1:nMu
        scores = [];
        labels = [];
        completeClasses = true;

        for iS = 1:nS
            vals = per_theta{kc_plot,iMu,iS,it_plot};
            vals = vals(isfinite(vals));

            % Check only this individual (mu,s,tau_max) dataset.
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
        else
            % Skip only this mu value at this observation time.
            all_scores_by_mu{iMu} = [];
            all_labels_by_mu{iMu} = [];
        end
    end

    % ------------------------------------------------------------
    % Optimal ROC point + threshold for ALL valid mu values
    % ------------------------------------------------------------
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

        % Find the threshold associated with MATLAB's optimal ROC point.
        % The distance formulation is robust to tiny floating-point
        % differences between OPTROCPT and the returned ROC coordinates.
        [~,idxOpt_all] = min((fpr_all-optROCPoint_all(1)).^2 + ...
                             (tpr_all-optROCPoint_all(2)).^2);

        optThreshold_all = thresholds_all(idxOpt_all);

        roc_opt_fpr(iMu,it_plot) = optROCPoint_all(1);
        roc_opt_tpr(iMu,it_plot) = optROCPoint_all(2);

        % Ignore +/-Inf thresholds, which can occur at the ROC endpoints.
        if isfinite(optThreshold_all)
            roc_opt_threshold(iMu,it_plot) = optThreshold_all;
        end
    end

    figROC = figure('Name',sprintf('ROC t=%g',tmax_values(it_plot)), ...
        'Color','w','Units','normalized', ...
        'Position',[0.18 0.18 0.52 0.58]);

    axROC = axes(figROC); hold(axROC,'on');

    set(axROC,'FontName',paperFont,'FontSize',fsTick, ...
        'LineWidth',axisLineWidth,'TickDir','in', ...
        'TickLabelInterpreter','latex');

    plot(axROC,[0 1],[0 1],'k--','LineWidth',1, ...
        'HandleVisibility','off');

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

        % perfcurve returns MATLAB's optimal ROC operating point as
        % [FPR, TPR]. With the balanced class sizes used here and equal
        % misclassification costs, this corresponds to maximizing TPR-FPR
        % (Youden's J statistic).
        [fpr,tpr,thresholds,auc,optROCPoint] = perfcurve(labels,scores,1);

        plot(axROC,fpr,tpr, ...
            'Color',colors(j,:), ...
            'LineWidth',1.6, ...
            'DisplayName',sprintf('\\mu=%.3g (AUC=%.2f)', ...
            mu_values(iMu),auc));

        % Mark the optimal operating point on this ROC curve.
        plot(axROC,optROCPoint(1),optROCPoint(2),'o', ...
            'MarkerSize',7, ...
            'MarkerFaceColor',colors(j,:), ...
            'MarkerEdgeColor','k', ...
            'LineWidth',0.8, ...
            'HandleVisibility','off');

        % The corresponding optimal threshold has already been computed
        % above for every valid mu value. Re-use the stored operating point
        % here so the ROC figures and threshold-vs-mu figure are consistent.

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

    try
        exportgraphics(figROC,fn,'Resolution',300);
    catch
        warning('Could not save %s',fn);
    end
end

%% ============================================================
%% OPTIMAL CLASSIFICATION THRESHOLD
% This figure shows the fitted s-hat value at the ROC operating point that
% provides the best separation of the two classes for each mu and tau_max.
%% Optimal ROC threshold vs mu for all observation horizons
%% ============================================================
figThresh = figure('Name','Optimal ROC threshold vs mu - all tmax', ...
    'Color','w','Units','normalized','Position',[0.15 0.15 0.7 0.6]);
axThresh = axes(figThresh); hold(axThresh,'on');

set(axThresh,'FontName',paperFont,'FontSize',fsTick,'Box','on', ...
    'LineWidth',axisLineWidth,'TickDir','in');

colorsT = orderedcolors("gem12");

% As in the AUC figure, plot |mu| on a logarithmic axis and display
% negative tick labels, because all mu values in the sweep are negative.
xpos = -mu_values;
useLog = (max(xpos)/min(xpos) > 10);

hThresh = gobjects(nT,1);
labelsThresh = cell(nT,1);

% Show every optimal threshold that lies inside the chosen plotting range.
% Out-of-range values are replaced by NaN, which breaks the plotted line.
% This prevents artificial connections across large threshold excursions,
% while still allowing valid points to reappear closer to mu = 0.
yMin = -2;
yMax = 2;

for it = 1:nT
    threshold_vec = roc_opt_threshold(:,it);

    if useLog
        xplot = xpos;
    else
        xplot = mu_values;
    end

    valid = isfinite(threshold_vec) & isfinite(xplot(:));
    inRange = valid & threshold_vec >= yMin & threshold_vec <= yMax;

    if ~any(inRange)
        continue
    end

    % Keep the full x-grid, but replace rejected threshold values by NaN.
    % MATLAB then breaks the line at those locations instead of connecting
    % distant accepted points across an out-of-range excursion.
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

    % Plot only the accepted points as markers.
    scatter(axThresh, xplot(inRange), threshold_vec(inRange), 36, 'o', ...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', color*0.6, ...
        'MarkerFaceAlpha', 0.9, ...
        'HandleVisibility', 'off');

    labelsThresh{it} = sprintf('$\\tau_{\\mathrm{max}}=%.4g$', ...
        tmax_values(it));
end

ylim(axThresh,[yMin yMax])

% Natural reference for the fitted cubic coefficient: threshold = 0.
yline(axThresh,0,'k--','LineWidth',1.1,'HandleVisibility','off');

validIdx = arrayfun(@isgraphics,hThresh);
hThresh = hThresh(validIdx);
labelsThresh = labelsThresh(validIdx);

if useLog
    set(axThresh,'XScale','log');

    xmin = min(xpos);
    xmax = max(xpos);
    emin = floor(log10(xmin));
    emax = ceil(log10(xmax));

    if emax-emin < 1
        emin = emin-1;
        emax = emax+1;
    end

    xlim(axThresh,[10^emin 10^emax]);
    ylim(axThresh,[-2 2])


    exps = emin:emax;
    set(axThresh,'XTick',10.^exps);

    xtick_labels = arrayfun(@(e) sprintf('$-10^{%d}$',e), ...
        exps,'UniformOutput',false);

    set(axThresh,'XTickLabel',xtick_labels, ...
        'TickLabelInterpreter','latex', ...
        'XDir','reverse');
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

try
    exportgraphics(figThresh,fnThresh,'Resolution',300);
    fprintf('Saved optimal-threshold figure to: %s\n',fnThresh);
catch ME
    warning('Could not save optimal-threshold figure');
end

% Append the newly derived operating-point quantities to the saved MAT file.
% This works both after a fresh simulation and when the original per_theta
% data were loaded from an existing file.
try
    save(savePath,'roc_opt_fpr','roc_opt_tpr','roc_opt_threshold','-append');
    fprintf('Appended ROC operating points and thresholds to: %s\n',savePath);
catch 
    warning('Could not append ROC operating-point data');
end


    end
end

%% ============================================================
%% Ridge / density panels
%% ============================================================
kc_list = unique([1,kc_shat]);
minForKDE = 3;

cmapPos = [0 0.4470 0.7410];
cmapNeg = [0.8500 0.3250 0.0980];

for kc = kc_list
    for it_plot = 1:nT

        maxMuPlot = 5;
        nMu_total = numel(mu_values);
        nMu_plot = min(nMu_total,maxMuPlot);

        idx_mu_plot = unique(round(linspace(1,nMu_total,nMu_plot)));
        nMu_plot = numel(idx_mu_plot);

        if kc == 1
            compLabel = '$\hat{\mu}$';
            plotTag = 'mu_hat';
        else
            compLabel = '$\hat{s}$';
            plotTag = 's_hat';
        end

        figRidge = figure('Name',sprintf('Ridge %s t=%g', ...
            plotTag,tmax_values(it_plot)), ...
            'Color','w','Units','normalized', ...
            'Position',[0.08 0.07 0.82 0.82]);

        t = tiledlayout(figRidge,nMu_plot,nS, ...
            'TileSpacing','compact','Padding','compact');

        for iMu_plot = 1:nMu_plot
            iMu = idx_mu_plot(iMu_plot);
            mu = mu_values(iMu);

            for iS = 1:nS
                ax = nexttile(t,(iMu_plot-1)*nS+iS);
                cla(ax); hold(ax,'on');

                set(ax,'FontName',paperFont,'FontSize',fsTick, ...
                    'LineWidth',axisLineWidth,'TickDir','in', ...
                    'TickLabelInterpreter','latex');

                scores = per_theta{kc,iMu,iS,it_plot};
                scores = scores(isfinite(scores));

                % Check only this individual (mu,s,tau_max) dataset.
                % An incomplete dataset at another time does not affect it.
                if numel(scores) >= n_realizations
                    scores = scores(1:n_realizations);
                end

                % s as column header.
                if iMu_plot == 1
                    title(ax,sprintf('$s=%g$',s_values(iS)), ...
                        'Interpreter','latex','FontSize',fsPanelTitle);
                end

                % Numerical density labels only on the left.
                if iS ~= 1
                    yticklabels(ax,[]);
                end

                % mu as vertical row label on the outer right.
                if iS == nS
                    text(ax,1.025,0.5, ...
                        sprintf('$\\mu=%.1f$',mu), ...
                        'Units','normalized', ...
                        'Interpreter','latex', ...
                        'FontSize',fsPanelTitle, ...
                        'Rotation',90, ...
                        'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle', ...
                        'Clipping','off');
                end

                % Fitted coefficient label below each column.
                if iMu_plot == nMu_plot
                    xlabel(ax,compLabel, ...
                        'Interpreter','latex','FontSize',fsLabel);
                end

                if numel(scores) < n_realizations
                    text(ax,0.5,0.5,'no samples', ...
                        'Units','normalized', ...
                        'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle', ...
                        'FontSize',fsTick);
                    xlim(ax,[0 1]);
                    ylim(ax,[0 1]);
                    xticks(ax,[]);
                    yticks(ax,[]);
                    box(ax,'on');
                    hold(ax,'off');
                    continue
                end

                scores = scores(:);

                if numel(scores) >= minForKDE
                    pr = prctile(scores,[1 99]);

                    if pr(1) == pr(2)
                        pr = pr+[-0.5 0.5];
                    end

                   prRange = pr(2) - pr(1);

                   xgrid = linspace( ...
                            pr(1) - 0.1*prRange, ...
                            pr(2) + 0.1*prRange, ...
                            400);

                                            [f,xi] = ksdensity(scores,xgrid,'Support','unbounded');

                    if max(f) > 0
                        f = f/max(f);
                    end

                    color = cmapPos;
                    if s_values(iS) < 0
                        color = cmapNeg;
                    end

                    fill(ax,[xi fliplr(xi)],[f zeros(size(f))],color, ...
                        'FaceAlpha',0.6,'EdgeColor','none');
                    plot(ax,xi,f,'k-','LineWidth',0.6);
                    xlim(ax,[min(xi) max(xi)]);
                else
                    color = (s_values(iS)>0)*cmapPos + ...
                        (s_values(iS)<0)*cmapNeg;
                    histogram(ax,scores,'Normalization','pdf', ...
                        'FaceAlpha',0.6,'FaceColor',color);
                end

                m = mean(scores,'omitnan');
                sd = std(scores,'omitnan');
                nSamp = numel(scores);

                yl = ylim(ax);
                if isempty(yl) || any(~isfinite(yl))
                    yl = [0 1];
                    ylim(ax,yl);
                end

                plot(ax,[m m],[0.02*yl(2) 0.98*yl(2)],'--', ...
                    'Color',0.15*[1 1 1]);

                ridgeText = sprintf('n=%d\nm=%.3g\nsd=%.3g', ...
                    nSamp,m,sd);

                text(ax,0.02,0.95,ridgeText, ...
                    'Units','normalized', ...
                    'VerticalAlignment','top', ...
                    'BackgroundColor','w', ...
                    'EdgeColor','k', ...
                    'FontSize',fsTextBox, ...
                    'Margin',1);

                if iS ~= 1
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

        outDir = fullfile(base_dir_root,'ridge_plots_all');
        if ~exist(outDir,'dir'), mkdir(outDir); end

        fname = fullfile(outDir,sprintf('ridge_%s_t=%g.jpg', ...
            plotTag,tmax_values(it_plot)));

        try
            exportgraphics(figRidge,fname,'Resolution',300);
        catch
            warning('Could not save ridge figure to %s',fname);
        end
    end
end

%% ============================================================
%% ROC curves per observation horizon
%% ============================================================
kc_plot = kc_shat;

outDirROC = fullfile(base_dir_root,'roc_plots_per_time');
if ~exist(outDirROC,'dir'), mkdir(outDirROC); end

% Match the SH presentation: show a representative set of mu values.
nROCcurves = 10;
target_mu_roc = -logspace(2,-1,nROCcurves);

idx_mu_roc = zeros(1,nROCcurves);
for j = 1:nROCcurves
    [~,idx_mu_roc(j)] = min(abs( ...
        log10(-mu_values)-log10(-target_mu_roc(j))));
end
idx_mu_roc = unique(idx_mu_roc,'stable');

% Store optimal operating-point information for EVERY valid (mu,tau_max)
% pair, not only for the representative mu values shown in the ROC figures.
% Rows correspond to mu, columns to tau_max.
roc_opt_fpr       = nan(nMu,nT);
roc_opt_tpr       = nan(nMu,nT);
roc_opt_threshold = nan(nMu,nT);

for it_plot = 1:nT
    all_scores_by_mu = cell(nMu,1);
    all_labels_by_mu = cell(nMu,1);

    for iMu = 1:nMu
        scores = [];
        labels = [];
        completeClasses = true;

        for iS = 1:nS
            vals = per_theta{kc_plot,iMu,iS,it_plot};
            vals = vals(isfinite(vals));

            % Check only this individual (mu,s,tau_max) dataset.
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
        else
            % Skip only this mu value at this observation time.
            all_scores_by_mu{iMu} = [];
            all_labels_by_mu{iMu} = [];
        end
    end

    % ------------------------------------------------------------
    % Optimal ROC point + threshold for ALL valid mu values
    % ------------------------------------------------------------
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

        % Find the threshold associated with MATLAB's optimal ROC point.
        % The distance formulation is robust to tiny floating-point
        % differences between OPTROCPT and the returned ROC coordinates.
        [~,idxOpt_all] = min((fpr_all-optROCPoint_all(1)).^2 + ...
                             (tpr_all-optROCPoint_all(2)).^2);

        optThreshold_all = thresholds_all(idxOpt_all);

        roc_opt_fpr(iMu,it_plot) = optROCPoint_all(1);
        roc_opt_tpr(iMu,it_plot) = optROCPoint_all(2);

        % Ignore +/-Inf thresholds, which can occur at the ROC endpoints.
        if isfinite(optThreshold_all)
            roc_opt_threshold(iMu,it_plot) = optThreshold_all;
        end
    end

    figROC = figure('Name',sprintf('ROC t=%g',tmax_values(it_plot)), ...
        'Color','w','Units','normalized', ...
        'Position',[0.18 0.18 0.52 0.58]);

    axROC = axes(figROC); hold(axROC,'on');

    set(axROC,'FontName',paperFont,'FontSize',fsTick, ...
        'LineWidth',axisLineWidth,'TickDir','in', ...
        'TickLabelInterpreter','latex');

    plot(axROC,[0 1],[0 1],'k--','LineWidth',1, ...
        'HandleVisibility','off');

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

        % perfcurve returns MATLAB's optimal ROC operating point as
        % [FPR, TPR]. With the balanced class sizes used here and equal
        % misclassification costs, this corresponds to maximizing TPR-FPR
        % (Youden's J statistic).
        [fpr,tpr,thresholds,auc,optROCPoint] = perfcurve(labels,scores,1);

        plot(axROC,fpr,tpr, ...
            'Color',colors(j,:), ...
            'LineWidth',1.6, ...
            'DisplayName',sprintf('\\mu=%.3g (AUC=%.2f)', ...
            mu_values(iMu),auc));

        % Mark the optimal operating point on this ROC curve.
        plot(axROC,optROCPoint(1),optROCPoint(2),'o', ...
            'MarkerSize',7, ...
            'MarkerFaceColor',colors(j,:), ...
            'MarkerEdgeColor','k', ...
            'LineWidth',0.8, ...
            'HandleVisibility','off');

        % The corresponding optimal threshold has already been computed
        % above for every valid mu value. Re-use the stored operating point
        % here so the ROC figures and threshold-vs-mu figure are consistent.

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

    try
        exportgraphics(figROC,fn,'Resolution',300);
    catch
        warning('Could not save %s',fn);
    end
end

%% ============================================================
%% Optimal ROC threshold vs mu for all observation horizons
%% ============================================================
figThresh = figure('Name','Optimal ROC threshold vs mu - all tmax', ...
    'Color','w','Units','normalized','Position',[0.15 0.15 0.7 0.6]);
axThresh = axes(figThresh); hold(axThresh,'on');

set(axThresh,'FontName',paperFont,'FontSize',fsTick,'Box','on', ...
    'LineWidth',axisLineWidth,'TickDir','in');

colorsT = orderedcolors("gem12");

% As in the AUC figure, plot |mu| on a logarithmic axis and display
% negative tick labels, because all mu values in the sweep are negative.
xpos = -mu_values;
useLog = (max(xpos)/min(xpos) > 10);

hThresh = gobjects(nT,1);
labelsThresh = cell(nT,1);

% Show every optimal threshold that lies inside the chosen plotting range.
% Out-of-range values are replaced by NaN, which breaks the plotted line.
% This prevents artificial connections across large threshold excursions,
% while still allowing valid points to reappear closer to mu = 0.
yMin = -2;
yMax = 2;

for it = 1:nT
    threshold_vec = roc_opt_threshold(:,it);

    if useLog
        xplot = xpos;
    else
        xplot = mu_values;
    end

    valid = isfinite(threshold_vec) & isfinite(xplot(:));
    inRange = valid & threshold_vec >= yMin & threshold_vec <= yMax;

    if ~any(inRange)
        continue
    end

    % Keep the full x-grid, but replace rejected threshold values by NaN.
    % MATLAB then breaks the line at those locations instead of connecting
    % distant accepted points across an out-of-range excursion.
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

    % Plot only the accepted points as markers.
    scatter(axThresh, xplot(inRange), threshold_vec(inRange), 36, 'o', ...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', color*0.6, ...
        'MarkerFaceAlpha', 0.9, ...
        'HandleVisibility', 'off');

    labelsThresh{it} = sprintf('$\\tau_{\\mathrm{max}}=%.4g$', ...
        tmax_values(it));
end

ylim(axThresh,[yMin yMax])

% Natural reference for the fitted cubic coefficient: threshold = 0.
yline(axThresh,0,'k--','LineWidth',1.1,'HandleVisibility','off');

validIdx = arrayfun(@isgraphics,hThresh);
hThresh = hThresh(validIdx);
labelsThresh = labelsThresh(validIdx);

if useLog
    set(axThresh,'XScale','log');

    xmin = min(xpos);
    xmax = max(xpos);
    emin = floor(log10(xmin));
    emax = ceil(log10(xmax));

    if emax-emin < 1
        emin = emin-1;
        emax = emax+1;
    end

    xlim(axThresh,[10^emin 10^emax]);
    ylim(axThresh,[-2 2])


    exps = emin:emax;
    set(axThresh,'XTick',10.^exps);

    xtick_labels = arrayfun(@(e) sprintf('$-10^{%d}$',e), ...
        exps,'UniformOutput',false);

    set(axThresh,'XTickLabel',xtick_labels, ...
        'TickLabelInterpreter','latex', ...
        'XDir','reverse');
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

try
    exportgraphics(figThresh,fnThresh,'Resolution',300);
    fprintf('Saved optimal-threshold figure to: %s\n',fnThresh);
catch ME
    warning('Could not save optimal-threshold figure');
end

% Append the newly derived operating-point quantities to the saved MAT file.
% This works both after a fresh simulation and when the original per_theta
% data were loaded from an existing file.
try
    save(savePath,'roc_opt_fpr','roc_opt_tpr','roc_opt_threshold','-append');
    fprintf('Appended ROC operating points and thresholds to: %s\n',savePath);
catch ME
    warning('Could not append ROC operating-point data');
end

