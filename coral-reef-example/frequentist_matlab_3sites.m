% Frequentist MLE fitting for coral reef model (3 sites)
% REQUIRES: Optimization Toolbox (lsqnonlin, MultiStart)

clear; clc; close all;

csvPath  = './datasets/new_reef_data.xlsx';
saveDir  = './results/';
% if ~exist(saveDir, 'dir'), mkdir(saveDir); end

datasets = [
    ds("North Direction Reef", 1993, 2012, [1 2 3])
    ds("Chicken Reef",         2009, 2023, [1 2 3])
    ds("Turner Reef",          2002, 2017, [1 2 3])
    ds("Chinaman Reef",        1992, 2004, [1 2 3])
    ds("SNAKE (22088)",        1994, 2008, [1 2 3])
    ds("HORSESHOE",            1999, 2012, [1 2 3])
    ds("Broomfield Reef",      2008, 2022, [1 2 3])
    ds("One Tree Reef",        1992, 2006, [1 2 3])
    ds("Lady Musgrave Reef",   1992, 2007, [1 2 3])
];

nStarts   = 200;
nDraws    = 4000;  % match Bayesian (4 chains x 1000 samples)
odeSolver = @ode45;

LB = [0, 0, 0, 0];
UB = [5, 5, 5, 5];

[~,~,ext] = fileparts(csvPath);
if strcmpi(ext, '.xlsx') || strcmpi(ext, '.xls')
    T = readtable(csvPath);
else
    T = readtable(csvPath, 'FileType','text', 'Delimiter',';');
end

if any(strcmpi(T.Properties.VariableNames, 'GROUP_CODE'))
    T = T(strcmpi(string(T.GROUP_CODE), 'Hard Coral'), :);
end

aliases.reef  = lower(strsplit("Reef,ReefName,reef_name,Location,REEF_NAME"));
aliases.site  = lower(strsplit("Site,SiteNo,Site_Number,site_no,SITE_NO"));
aliases.year  = lower(strsplit("Year,SurveyYear,Time,SAMPLE_DATE,YEAR_CODE"));
aliases.cover = lower(strsplit("CoralCover,coral_cover,LiveCoral,COVER,Live.Coral"));

T = standardize_table(T, aliases);

results = table('Size',[0 16], 'VariableTypes', ...
    {'string','double','double','string','double','double','double','double', ...
     'double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'Reef','YearStart','YearEnd','SitesIncluded', ...
     'alpha','beta','gamma','mu','alpha_lo','alpha_hi','beta_lo','beta_hi', ...
     'gamma_lo','gamma_hi','mu_lo','mu_hi'});

for k = 1:numel(datasets)
    spec = datasets(k);
    [D_fit, D_sites] = subset_dataset(T, spec);
    
    if height(D_fit) < 4
        continue;
    end

    D = sortrows(D_fit, 'Year');
    t_obs = D.Year;
    t0 = t_obs(1);
    tRel = t_obs - t0;
    N_obs = D.N;

    N0 = N_obs(1);
    x0 = [N0*N0; N0*(1-N0)];

    fun = @(p) residuals(p, tRel, N_obs, x0, odeSolver);

    X0 = LB + rand(nStarts, 4).*(UB - LB);

    opts = optimoptions('lsqnonlin', 'Display','off', ...
        'MaxFunctionEvaluations',2e5, 'StepTolerance',1e-12);

    problem = createOptimProblem('lsqnonlin', 'objective',fun, ...
        'x0',X0(1,:), 'lb',LB, 'ub',UB, 'options',opts);

    ms = MultiStart('UseParallel',true, 'Display','off');
    startPts = CustomStartPointSet(X0);

    [pHat0, ~] = run(ms, problem, startPts);

    [pHat, resnorm, ~, ~, ~, ~, J] = lsqnonlin(fun, pHat0, LB, UB, opts);

    dof = numel(N_obs) - numel(pHat);
    s2 = resnorm / max(dof, 1);
    JTJ = full(J)' * full(J);
    JTJ = (JTJ + JTJ')/2 + 1e-10*eye(numel(pHat));
    covP = s2 * pinv(JTJ);
    
    se = sqrt(max(diag(covP), 0));
    CI = [pHat(:) - 1.96*se, pHat(:) + 1.96*se];
    CI(:,1) = max(0, CI(:,1));  % clip lower bounds at 0 (parameters are non-negative)

    Pdraws = mvnrnd_local(pHat(:)', covP, nDraws);
    Pdraws = Pdraws(all(Pdraws > 0, 2), :);
    if size(Pdraws,1) < nDraws
        Pdraws(end+1:nDraws,:) = repmat(pHat(:)', nDraws-size(Pdraws,1), 1);
    end

    [tFine, Nmean, Cmean, Bmean, Nall] = compute_predictions(Pdraws, x0, odeSolver, max(tRel));

    sdN = std(Nall, 0, 2);
    NmeanP = 100*Nmean; 
    CmeanP = 100*Cmean; 
    BmeanP = 100*Bmean;
    N1loP = max(0, 100*(Nmean - sdN));  N1hiP = min(100, 100*(Nmean + sdN));
    N2loP = max(0, 100*(Nmean - 2*sdN)); N2hiP = min(100, 100*(Nmean + 2*sdN));
    N3loP = max(0, 100*(Nmean - 3*sdN)); N3hiP = min(100, 100*(Nmean + 3*sdN));

    figure('Color','w', 'Position',[120 120 900 520]); hold on;

    tFine_years = tFine + t0;

    fill([tFine_years; flipud(tFine_years)], [N3loP; flipud(N3hiP)], ...
        0.92*[1 1 1], 'EdgeColor','none');
    fill([tFine_years; flipud(tFine_years)], [N2loP; flipud(N2hiP)], ...
        0.85*[1 1 1], 'EdgeColor','none');
    h1 = fill([tFine_years; flipud(tFine_years)], [N1loP; flipud(N1hiP)], ...
        0.78*[1 1 1], 'EdgeColor','none');

    hN = plot(tFine_years, NmeanP, 'Color',[0.2 0.6 0.2], 'LineWidth',1.8);
    hC = plot(tFine_years, CmeanP, '--', 'LineWidth',1.4, 'Color',[0.15 0.35 0.85]);
    hB = plot(tFine_years, BmeanP, '--', 'LineWidth',1.4, 'Color',[0.85 0.25 0.25]);

    siteIDs = unique(D_sites.Site);
    mk = {'o','s','^','d','v','>'};
    col = lines(max(3, numel(siteIDs)));
    hSites = gobjects(numel(siteIDs), 1);
    siteLabels = cell(numel(siteIDs), 1);
    for si = 1:numel(siteIDs)
        mask = D_sites.Site == siteIDs(si);
        hSites(si) = plot(D_sites.Year(mask), 100*D_sites.N(mask), ...
            mk{mod(si-1, numel(mk))+1}, 'MarkerFaceColor','w', ...
            'MarkerEdgeColor','k', 'Color',col(si,:), 'LineStyle','none', ...
            'LineWidth',1.0, 'MarkerSize',6);
        siteLabels{si} = sprintf('Site %d', siteIDs(si));
    end

    xlabel('Year'); ylabel('% Coral Cover');
    xlim([spec.y0 spec.y1]);
    ylim([0 100]);
    grid on; box on;
    set(gca, 'FontSize',11, 'LineWidth',1.2);

    text(0.01, 0.98, char(spec.name), 'Units','normalized', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'FontSize',14, 'FontWeight','bold');

    % move legend to north for last 3 reefs to avoid overlap
    if contains(lower(spec.name), {'broomfield', 'one tree', 'lady musgrave'})
        legend([hN hC hB hSites'], [{'Total N(t)', 'Healthy C(t)', 'Bleached B(t)'}, siteLabels'], ...
            'Location','north', 'FontSize',9);
    else
        legend([hN hC hB hSites'], [{'Total N(t)', 'Healthy C(t)', 'Bleached B(t)'}, siteLabels'], ...
            'Location','northeast', 'FontSize',9);
    end

    reefSafe = regexprep(char(spec.name), '\W+', '_');
    sitesStr = strjoin(string(spec.sites), '-');
    fname = sprintf('%s_%d-%d_sites_%s', reefSafe, spec.y0, spec.y1, sitesStr);
    exportgraphics(gcf, fullfile(saveDir, [fname '.pdf']), 'ContentType','vector');
    exportgraphics(gcf, fullfile(saveDir, [fname '.png']), 'Resolution',300);

    results = [results; {string(spec.name), spec.y0, spec.y1, ...
        strjoin(string(spec.sites),','), pHat(1), pHat(2), pHat(3), pHat(4), ...
        CI(1,1), CI(1,2), CI(2,1), CI(2,2), CI(3,1), CI(3,2), CI(4,1), CI(4,2)}];

    preds = table(tFine,tFine_years,NmeanP,sdN,N1loP,N1hiP,N2loP,N2hiP,N3loP,N3hiP,CmeanP,BmeanP, ...
    'VariableNames', {'t_rel','year','N_mean','N_sd', ...
     'N_lo1','N_hi1','N_lo2','N_hi2','N_lo3','N_hi3','C_mean','B_mean'});

    writetable(preds, fullfile(saveDir, reefSafe + "_predictions_freq.csv"));


end

writetable(results, fullfile(saveDir, 'results_table.csv'));


function S = ds(name, y0, y1, sites)
    S = struct('name',name, 'y0',y0, 'y1',y1, 'sites',sites);
end

function parts = strsplit(str)
    C = regexp(str, '\s*,\s*', 'split');
    parts = string(C);
end

function T = standardize_table(T, aliases)
    v = lower(string(T.Properties.VariableNames));
    findcol = @(cands) find(ismember(v, cands), 1);

    iReef = findcol(aliases.reef);
    iSite = findcol(aliases.site);
    iYear = findcol(aliases.year);
    iCover = findcol(aliases.cover);

    T = renamevars(T, T.Properties.VariableNames([iReef iSite iYear iCover]), ...
        {'Reef','Site','Year','N'});

    T.Reef = string(T.Reef);
    T.Site = double(T.Site);
    
    if isdatetime(T.Year)
        T.Year = year(T.Year);
    else
        T.Year = double(T.Year);
    end

    T.N = double(T.N);
    if median(T.N, 'omitnan') > 1
        T.N = T.N / 100;
    end
    T.N = max(min(T.N, 1), 0);

    T = rmmissing(T, 'DataVariables', {'Reef','Site','Year','N'});
end

function [D_fit, D_sites] = subset_dataset(T, spec)
    normfun = @(s) regexprep(upper(string(s)), '[^A-Z0-9 ]', '');
    reefNorm = normfun(T.Reef);
    keySpec = normfun(spec.name);
    
    maskReef = contains(reefNorm, keySpec);
    D = T(maskReef, :);

    availSites = unique(D.Site(~isnan(D.Site)));
    useSites = intersect(availSites, spec.sites);
    if isempty(useSites)
        useSites = availSites;
    end

    D = D(ismember(D.Site, useSites) & D.Year >= spec.y0 & D.Year <= spec.y1, :);

    D_sites = sortrows(D, {'Site','Year'});

    years = unique(D.Year);
    Nmean = arrayfun(@(y) mean(D.N(D.Year == y)), years);
    D_fit = table(repmat(string(spec.name), numel(years), 1), ...
        NaN(numel(years), 1), years, Nmean, ...
        'VariableNames', {'Reef','Site','Year','N'});
end

function r = residuals(p, tRel, N_obs, x0, odeSolver)
    if any(p <= 0)
        r = inf(size(N_obs));
        return;
    end

    try
        sol = odeSolver(@(t,x) coral_ode(t, x, p), [0 max(tRel)], x0);
        Npred = arrayfun(@(t) sum(deval(sol, t)), tRel);
        r = Npred(:) - N_obs(:);
    catch
        r = inf(size(N_obs));
    end
end

function dx = coral_ode(~, x, p)
    C = x(1);
    B = x(2);
    N = C + B;
    
    alpha = p(1);
    beta = p(2);
    gamma = p(3);
    mu = p(4);
    
    dC = alpha*C*(1 - N) - beta*C + gamma*B;
    dB = beta*C - gamma*B - mu*B;
    
    dx = [dC; dB];
end

function [tFine, Nmean, Cmean, Bmean, Nall] = compute_predictions(Pdraws, x0, odeSolver, tMax)
    tFine = linspace(0, tMax, 200)';
    nD = size(Pdraws, 1);
    Nall = zeros(numel(tFine), nD);
    Call = Nall;
    Ball = Nall;

    parfor j = 1:nD
        p = Pdraws(j, :);
        sol = odeSolver(@(t,x) coral_ode(t, x, p), [0 tMax], x0);
        C = zeros(numel(tFine), 1);
        B = zeros(numel(tFine), 1);
        for i = 1:numel(tFine)
            xi = deval(sol, tFine(i));
            C(i) = xi(1);
            B(i) = xi(2);
        end
        Nall(:,j) = C + B;
        Call(:,j) = C;
        Ball(:,j) = B;
    end

    Nmean = mean(Nall, 2);
    Cmean = mean(Call, 2);
    Bmean = mean(Ball, 2);
end

function X = mvnrnd_local(mu, Sigma, n)
    mu = mu(:)';
    d = numel(mu);
    Sigma = (Sigma + Sigma')/2;
    [R, p] = chol(Sigma);
    if p ~= 0
        [V, D] = eig(Sigma);
        D = max(diag(D), 0);
        R = (V * diag(sqrt(D)))';
    end
    Z = randn(n, d);
    X = bsxfun(@plus, Z*R, mu);
end
