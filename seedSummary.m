% =========================================================================
% seedSummary.m  [N7]  Further statistics from seedStudy.mat (no new runs).
%
% Prints, as mean +- sd over the seeds with the value of the reference run
% (seed 42), the quantities quoted in the thesis that seedStudy.m saves but
% does not print:
%   - volumes over the threshold grid D = 3, 2, 1, 0.5, 0.25;
%   - share of the volume without eavesdropper, ratios between adversaries
%     and between thresholds;
%   - transmission times and rates at D = 1, rates over the window;
%   - accumulated deflection at the end of the pass;
%   - epoch in which each adversary stops Alice (how often over the seeds).
% =========================================================================
clc; clear;
S = load('seedStudy.mat');  R = S.R;  Rf = S.Rref;  nS = numel(R);
fprintf('=== seedSummary: %d seeds, reference run = seed %d ===\n', nS, Rf.seed);
fprintf('  %-44s %10s %10s %10s %10s %10s %8s\n', 'quantity', 'mean', 'sd', 'min', 'max', 'seed 42', 'z(42)');

eve = {'Eve-A', 'Eve-B', 'Eve-C'};
thr = Rf.thr_grid;

% Volumes over the threshold grid.
fprintf('--- volume over the threshold grid [bit] ---\n');
for iE = 1:3
    for t = 1:numel(thr)
        v = arrayfun(@(r) r.volGrid(iE,t), R);
        row(sprintf('%s, D = %g', eve{iE}, thr(t)), v, Rf.volGrid(iE,t), '%10.1f');
    end
end

% Shares and ratios.
fprintf('--- shares and ratios ---\n');
noEve = arrayfun(@(r) r.vol3(1), R);
for iE = 1:3
    v = arrayfun(@(r) r.vol3(iE+1), R)./noEve*100;
    row(sprintf('share of no-eavesdropper volume, %s [%%]', eve{iE}), v, Rf.vol3(iE+1)/Rf.vol3(1)*100, '%10.2f');
end
for iE = 1:3
    v = arrayfun(@(r) r.vol3(iE+1)/r.vol1(iE), R);
    row(sprintf('volume D = 3 / D = 1, %s', eve{iE}), v, Rf.vol3(iE+1)/Rf.vol1(iE), '%10.2f');
end
v = arrayfun(@(r) r.vol3(2)/r.vol3(3), R);
row('volume Eve-A / Eve-B, D = 3', v, Rf.vol3(2)/Rf.vol3(3), '%10.2f');
v = arrayfun(@(r) r.vol3(3)/r.vol3(4), R);
row('volume Eve-B / Eve-C, D = 3 (hopping gain)', v, Rf.vol3(3)/Rf.vol3(4), '%10.2f');
v = arrayfun(@(r) r.vol1(2)/r.vol1(3), R);
row('volume Eve-B / Eve-C, D = 1', v, Rf.vol1(2)/Rf.vol1(3), '%10.2f');

% Times and rates.
fprintf('--- transmission times and rates ---\n');
for iE = 1:3
    row(sprintf('transmission time %s, D = 3 [s]', eve{iE}), arrayfun(@(r) r.tTx3(iE), R), Rf.tTx3(iE), '%10.1f');
    row(sprintf('transmission time %s, D = 1 [s]', eve{iE}), arrayfun(@(r) r.tTx1(iE), R), Rf.tTx1(iE), '%10.1f');
end
for iE = 1:3
    row(sprintf('covert rate %s, D = 1 [bit/s]', eve{iE}), arrayfun(@(r) r.rateTx1(iE), R), Rf.rateTx1(iE), '%10.2f');
    row(sprintf('rate over window %s, D = 3 [bit/s]', eve{iE}), arrayfun(@(r) r.rateWin3(iE), R), Rf.rateWin3(iE), '%10.3f');
    row(sprintf('rate over window %s, D = 1 [bit/s]', eve{iE}), arrayfun(@(r) r.rateWin1(iE), R), Rf.rateWin1(iE), '%10.3f');
end

% Evidence at the end of the pass.
fprintf('--- accumulated deflection at the end of the pass ---\n');
for iE = 1:3
    row(sprintf('%s', eve{iE}), arrayfun(@(r) r.eveEnd(iE), R), Rf.eveEnd(iE), '%10.2f');
end

% Epoch of the stop.
fprintf('--- epoch in which Alice stops (start time of the epoch: count over the seeds) ---\n');
for D = [3 1]
    for iE = 1:3
        if D == 3, k = arrayfun(@(r) r.is3(iE), R); kr = Rf.is3(iE);
        else,      k = arrayfun(@(r) r.is1(iE), R); kr = Rf.is1(iE); end
        [u, ~, ic] = unique(k);
        fprintf('  %s, D = %d: ', eve{iE}, D);
        for i = 1:numel(u), fprintf('%3.0f s (x%d)  ', Rf.t_ep(u(i)), sum(ic == i)); end
        fprintf('| seed 42: %3.0f s\n', Rf.t_ep(kr));
    end
end

function row(label, v, r, f)
% One line: mean, sd, min, max, reference value, z of the reference.
m = mean(v);  sd = std(v);
if sd > 1e-9*max(abs(m),1), z = (r - m)/sd; else, sd = 0; z = 0; end
fprintf(['  %-44s ' f ' ' f ' ' f ' ' f ' ' f ' %+8.2f\n'], label, m, sd, min(v), max(v), r, z);
end