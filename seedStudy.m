% =========================================================================
% seedStudy.m  [N7]  Variability of the results over the master seed.
%
% 1. Check A2: runSeed(42) must reproduce the reference run of
%    CovertSimFix.m (output-run-riferimento.txt). If it does not, the study
%    stops with an error.
% 2. Runs runSeed for the seeds 1-30. Each run re-optimises the packet
%    separation s, so the spread includes the choice of the operating point.
% 3. Prints mean, standard deviation, minimum and maximum of every quantity,
%    with the value of the reference run (seed 42) and its distance from the
%    mean in standard deviations.
% 4. Saves seedStudy.mat (R, Rref, seeds, logs). If the file is in the same
%    folder as CovertSimFix.m, block 9 adds the bands to Fig 2b.
%
% Needs runSeed.m in the same folder. The printout of each run is captured
% with evalc and saved in logs, so the command window shows only the study.
% =========================================================================
clc; clear; close all;

seeds   = 1:30;
seedRef = 42;

%% ------------------------------------------------------------------------
% A2. Equivalence check on the reference seed
% -------------------------------------------------------------------------
% Values printed by the reference run (claude/output-run-riferimento.txt),
% compared at the precision of the printout.
fprintf('=== A2: runSeed(%d) against the reference run ===\n', seedRef);
logRef = evalc('Rref = runSeed(seedRef);');
chk = {
  'optimised s',              sprintf('%.2f', Rref.s),        '6.65'
  'design volume vs Eve-B',   sprintf('%.3e', Rref.B_design), '5.115e+03'
  'volume, no eavesdropper',  sprintf('%.3e', Rref.vol3(1)),  '6.720e+04'
  'volume vs Eve-A, D = 3',   sprintf('%.3e', Rref.vol3(2)),  '1.204e+04'
  'volume vs Eve-B, D = 3',   sprintf('%.3e', Rref.vol3(3)),  '4.402e+03'
  'volume vs Eve-C, D = 3',   sprintf('%.3e', Rref.vol3(4)),  '3.199e+02'
  'volume vs Eve-A, D = 1',   sprintf('%.3e', Rref.vol1(1)),  '5.019e+03'
  'volume vs Eve-B, D = 1',   sprintf('%.3e', Rref.vol1(2)),  '1.308e+03'
  'volume vs Eve-C, D = 1',   sprintf('%.3e', Rref.vol1(3)),  '3.555e+01'
  'stop Eve-A, D = 3 [s]',    sprintf('%.0f', Rref.ts3(1)),   '410'
  'stop Eve-B, D = 3 [s]',    sprintf('%.0f', Rref.ts3(2)),   '316'
  'stop Eve-C, D = 3 [s]',    sprintf('%.0f', Rref.ts3(3)),   '42'
  'BER at closest approach',  sprintf('%.3f', Rref.BER_cross),'0.114'
  'Eve-B statistic, end',     sprintf('%.2f', Rref.eveEnd(2)),'53.40'
};
okAll = true;
for i = 1:size(chk,1)
    ok = strcmp(chk{i,2}, chk{i,3});
    okAll = okAll && ok;
    fprintf('  %-26s %12s  expected %12s  %s\n', chk{i,1}, chk{i,2}, chk{i,3}, ...
            char(ok*'OK  ' + ~ok*'FAIL'));
end
if ~okAll
    error('A2 failed: runSeed(42) does not reproduce the reference run. Study stopped.');
end
fprintf('A2 passed: runSeed reproduces the reference run.\n');

%% ------------------------------------------------------------------------
% Runs over the seeds
% -------------------------------------------------------------------------
nS   = numel(seeds);
logs = cell(1, nS);
fprintf('\n=== Runs over %d seeds ===\n', nS);
tStart = tic;
for j = 1:nS
    logs{j} = evalc('out = runSeed(seeds(j));');
    if j == 1, R = repmat(out, 1, nS); end
    R(j) = out;
    fprintf('  seed %2d | s = %.2f | vol vs Eve-B = %8.1f bit (D = 3), %7.1f bit (D = 1) | stop %3.0f s | %.0f s elapsed\n', ...
            seeds(j), out.s, out.vol3(3), out.vol1(2), out.ts3(2), toc(tStart));
end
save('seedStudy.mat', 'R', 'Rref', 'seeds', 'seedRef', 'logs', 'logRef');
fprintf('Saved seedStudy.mat\n');

%% ------------------------------------------------------------------------
% Summary: mean, standard deviation, range, reference run
% -------------------------------------------------------------------------
% Each row: label, values over the seeds, value of the reference run, format.
rows = {
  'packet separation s',                 [R.s],                 Rref.s,            '%.2f'
  'design volume vs Eve-B [bit]',        [R.B_design],          Rref.B_design,     '%.0f'
  'volume, no eavesdropper [bit]',       vec(R,'vol3',1),       Rref.vol3(1),      '%.0f'
  'volume vs Eve-A, D = 3 [bit]',        vec(R,'vol3',2),       Rref.vol3(2),      '%.0f'
  'volume vs Eve-B, D = 3 [bit]',        vec(R,'vol3',3),       Rref.vol3(3),      '%.0f'
  'volume vs Eve-C, D = 3 [bit]',        vec(R,'vol3',4),       Rref.vol3(4),      '%.1f'
  'volume vs Eve-A, D = 1 [bit]',        vec(R,'vol1',1),       Rref.vol1(1),      '%.0f'
  'volume vs Eve-B, D = 1 [bit]',        vec(R,'vol1',2),       Rref.vol1(2),      '%.0f'
  'volume vs Eve-C, D = 1 [bit]',        vec(R,'vol1',3),       Rref.vol1(3),      '%.2f'
  'Eve-B lower bound, D = 3 [bit]',      [R.B_Bc],              Rref.B_Bc,         '%.0f'
  'stop Eve-A, D = 3 [s]',               vec(R,'ts3',1),        Rref.ts3(1),       '%.0f'
  'stop Eve-B, D = 3 [s]',               vec(R,'ts3',2),        Rref.ts3(2),       '%.0f'
  'stop Eve-C, D = 3 [s]',               vec(R,'ts3',3),        Rref.ts3(3),       '%.0f'
  'stop Eve-A, D = 1 [s]',               vec(R,'ts1',1),        Rref.ts1(1),       '%.0f'
  'stop Eve-B, D = 1 [s]',               vec(R,'ts1',2),        Rref.ts1(2),       '%.0f'
  'stop Eve-C, D = 1 [s]',               vec(R,'ts1',3),        Rref.ts1(3),       '%.0f'
  'covert rate Eve-A, D = 3 [bit/s]',    vec(R,'rateTx3',1),    Rref.rateTx3(1),   '%.2f'
  'covert rate Eve-B, D = 3 [bit/s]',    vec(R,'rateTx3',2),    Rref.rateTx3(2),   '%.2f'
  'covert rate Eve-C, D = 3 [bit/s]',    vec(R,'rateTx3',3),    Rref.rateTx3(3),   '%.2f'
  'covert rate Eve-B, D = 1 [bit/s]',    vec(R,'rateTx1',2),    Rref.rateTx1(2),   '%.2f'
  'rate over window Eve-B, D = 3 [bit/s]', vec(R,'rateWin3',2), Rref.rateWin3(2),  '%.2f'
  'BER at closest approach',             [R.BER_cross],         Rref.BER_cross,    '%.4f'
  'capacity at closest approach',        [R.C_cross],           Rref.C_cross,      '%.4f'
  'Eve-B statistic at end of pass',      vec(R,'eveEnd',2),     Rref.eveEnd(2),    '%.2f'
  'mean pointing loss [dB]',             -10*log10([R.L_point_mean]), -10*log10(Rref.L_point_mean), '%.3f'
  'fading factor of the adversaries',    [R.eveGain],           Rref.eveGain,      '%.4f'
  'Campaign A, max relative deviation',  [R.campA_maxRel],      Rref.campA_maxRel, '%.4f'
  'Campaign A, max deviation [sd]',      [R.campA_maxZ],        Rref.campA_maxZ,   '%.2f'
};
fprintf('\n=== Variability over %d seeds (seeds %d-%d), reference run = seed %d ===\n', ...
        nS, seeds(1), seeds(end), seedRef);
fprintf('  %-38s %12s %12s %9s %12s %12s %12s %8s\n', 'quantity', 'mean', 'sd', 'sd/mean', ...
        'min', 'max', 'seed 42', 'z(42)');
for i = 1:size(rows,1)
    v = rows{i,2};  r = rows{i,3};  f = rows{i,4};
    m = mean(v);  sd = std(v);
    if sd > 1e-9*max(abs(m),1), z = (r - m)/sd; else, sd = 0; z = 0; end   % sd = 0: no spread
    fprintf(['  %-38s ' pad(f) ' ' pad(f) ' %8.1f%% ' pad(f) ' ' pad(f) ' ' pad(f) ' %+8.2f\n'], ...
            rows{i,1}, m, sd, 100*sd/abs(m), min(v), max(v), r, z);
end

% Distribution of the optimised separation.
[us, ~, ic] = unique([R.s]);
fprintf('\n  optimised s over the seeds: ');
for i = 1:numel(us), fprintf('%.2f (x%d)  ', us(i), sum(ic == i)); end
fprintf('\n');

%% ------------------------------------------------------------------------
% B2. Decision point: where the reference run sits
% -------------------------------------------------------------------------
vB  = vec(R,'vol3',3);
zB  = (Rref.vol3(3) - mean(vB))/std(vB);
fprintf('\n=== B2: reference run against the seeds ===\n');
fprintf('  volume vs Eve-B, D = 3: seed 42 = %.0f bit, seeds = %.0f +- %.0f bit (z = %+.2f)\n', ...
        Rref.vol3(3), mean(vB), std(vB), zB);
fprintf('  s: seed 42 = %.2f, seeds: mean %.2f, range %.2f - %.2f\n', ...
        Rref.s, mean([R.s]), min([R.s]), max([R.s]));
if abs(zB) > 1
    fprintf('  -> |z| > 1: decide with Giulio whether the thesis quotes the mean or seed 42.\n');
else
    fprintf('  -> |z| <= 1: seed 42 is representative; the thesis keeps its numbers and adds +- sd.\n');
end

%% ------------------------------------------------------------------------
function v = vec(R, field, i)
% Element i of a vector field, for every run.
    v = arrayfun(@(r) r.(field)(i), R);
end

function f = pad(f)
% Width-12 version of a numeric format.
    f = strrep(f, '%', '%12');
end