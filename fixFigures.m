% =========================================================================
% fixFigures.m
% Run AFTER NewCovertSimFix.m, in the same MATLAB session (same workspace).
% It does NOT modify the simulator: it redraws three panels from the
% variables the run leaves in the workspace and brings the four figures to
% the front. Nothing is saved: export them by hand (File > Save As, or
% Edit > Copy Figure / exportgraphics), with the file names of Chapter 6.
%
%   Fig 2(b): cumulative volume drawn at the END of each epoch, so that the
%             black curve meets each adversary line exactly at its stop
%             time; stop points marked.
%   Fig 3(b): part of the curve above OD = 1 dashed and shaded (outside the
%             range of the linear imprint).
%   Fig 4(a): marker at the alphabet allowed by the sampling bound alone.
% =========================================================================

% Time convention for Fig 2(b).
%   'printed'  : epoch k spans [t_ep(k), t_ep(k)+45 s), the same convention
%                as stopEpochTime, so the stop times match the printout and
%                the text (303 s, 407 s, 30 s). The last epoch, sampled at
%                the end of the window, is drawn up to 945 s.
%   'physical' : epoch k spans its real weight of 900/21 = 42.9 s, so the
%                axis ends at 900 s; the stop times shift by up to ~20 s
%                (e.g. Eve-C at ~289 s instead of 303 s).
timeConv = 'printed';

H2 = @(p) -p.*log2(max(p,eps)) - (1-p).*log2(max(1-p,eps));   % binary entropy

%% ------------------------------------------------------------------------
% FIG 2(b): cumulative covert volume
% -------------------------------------------------------------------------
f2 = findobj('Type', 'figure', 'Name', 'Fig 2 Orbital pass and eavesdropper accumulation');
ax = findPanel(f2, '(b)');
legend(ax, 'off');  cla(ax);  hold(ax, 'on');

switch timeConv
    case 'printed',  dtE = t_ep(2) - t_ep(1);       % 45 s
    case 'physical', dtE = T_window/N_epochs;       % 42.9 s
end
tEdge = (0:N_epochs)*dtE;                    % start of epoch 1 ... end of epoch 21
cumE  = [0, bits_cum]/1e3;                   % bits are uniform inside an epoch

plot(ax, tEdge, cumE, 'k-', 'LineWidth', 1.5, 'DisplayName', 'no eavesdropper');

vols   = [B_covert_A, B_covert_C, B_covert_B]/1e3;
kStop  = [isA, isC, isB];
tStopP = [tsA, tsC, tsB];
names  = {'vs Eve-A (envelope energy)', 'vs Eve-C (hop-blind)', 'vs Eve-B (key compromised)'};
cols   = {[0.47 0.67 0.19], [0 0.45 0.74], [0.85 0.33 0.10]};
styles = {'--', '-.', '--'};

for i = 1:3
    yline(ax, vols(i), styles{i}, 'Color', cols{i}, 'LineWidth', 1.3, ...
          'DisplayName', names{i});
    if isnan(tStopP(i)), continue; end               % threshold never reached
    k = kStop(i);
    % Fraction of the crossing epoch already transmitted at the stop. It is
    % the same fraction volumeAtThreshold uses, since both D^2 and the bits
    % grow linearly with the number of symbols inside an epoch.
    if k < N_epochs, span = t_ep(k+1) - t_ep(k); else, span = t_ep(k) - t_ep(k-1); end
    frac = (tStopP(i) - t_ep(k))/span;
    tS   = (k - 1 + frac)*dtE;
    plot(ax, [tS tS], [0 vols(i)], ':', 'Color', cols{i}, 'LineWidth', 1.1, ...
         'HandleVisibility', 'off');
    plot(ax, tS, vols(i), 'o', 'MarkerSize', 6, 'MarkerFaceColor', cols{i}, ...
         'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    text(ax, tS + 12, vols(i), sprintf('%.0f s', tS), 'Color', cols{i}, ...
         'FontSize', 8, 'VerticalAlignment', 'bottom');
end

xline(ax, t_sec(idx_cross), 'k:', 'Closest approach', ...
      'LabelOrientation', 'horizontal', 'HandleVisibility', 'off');
xlim(ax, [0 tEdge(end)]);  ylim(ax, [0 1.1*max(cumE)]);
xlabel(ax, 'Time from Pass Start [s]');
ylabel(ax, 'Cumulative Volume [kbit]');
title(ax, '(b) Covert Information Volume');
legend(ax, 'Location', 'northwest', 'FontSize', 7.5);
grid(ax, 'on');  hold(ax, 'off');

%% ------------------------------------------------------------------------
% FIG 3(b): imprint depth, valid range highlighted
% -------------------------------------------------------------------------
f3 = findobj('Type', 'figure', 'Name', 'Fig 3 Design optimisation and sensitivities');
ax = findPanel(f3, '(b)');
legend(ax, 'off');  cla(ax);  hold(ax, 'on');

OD_fine = logspace(log10(cfg.OD/20), log10(cfg.OD*8), 300);
d_fine  = (OD_fine/cfg.OD)*(g_mf/sigma_ref);         % same as the simulator
C_fine  = 1 - H2(0.5*erfc(d_fine/sqrt(2)));
in      = OD_fine <= 1;

patch(ax, [1 max(OD_fine) max(OD_fine) 1], [0 0 1 1], [0.9 0.9 0.9], ...
      'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, OD_fine(in),  C_fine(in),  '-',  'Color', [0 0.45 0.74], 'LineWidth', 1.6, ...
     'DisplayName', 'linear imprint (OD \leq 1)');
plot(ax, OD_fine(~in), C_fine(~in), '--', 'Color', [0 0.45 0.74], 'LineWidth', 1.2, ...
     'DisplayName', 'extrapolation');
d_dec = g_mf/sigma_ref;
plot(ax, cfg.OD, 1 - H2(0.5*erfc(d_dec/sqrt(2))), 'ko', 'MarkerFaceColor', 'k', ...
     'HandleVisibility', 'off');
xline(ax, cfg.OD, 'k:', 'declared', 'LineWidth', 1.2, 'HandleVisibility', 'off');
set(ax, 'XScale', 'log', 'Layer', 'top');
xlim(ax, [min(OD_fine) max(OD_fine)]);  ylim(ax, [0 1]);
xlabel(ax, 'Imprint Depth [OD]');
ylabel(ax, 'Capacity [bit/symbol]');
title(ax, '(b) Imprint Depth Sensitivity');
legend(ax, 'Location', 'northwest', 'FontSize', 7.5);
grid(ax, 'on');  hold(ax, 'off');

%% ------------------------------------------------------------------------
% FIG 4(a): reading at the sampling-limited alphabet
% -------------------------------------------------------------------------
f4 = findobj('Type', 'figure', 'Name', 'Fig 4 Hopping alphabet');
ax = findPanel(f4, '(a)');
delete(findobj(ax, 'Tag', 'fixMarker'));     % safe to run the script twice
hold(ax, 'on');
B_nyq = 10^interp1(log10(M_vec), log10(max(B_vec,eps)), log10(Mres_nyq));
plot(ax, Mres_nyq, B_nyq, 'kd', 'MarkerFaceColor', [0.85 0.33 0.10], ...
     'HandleVisibility', 'off', 'Tag', 'fixMarker');
text(ax, Mres_nyq, B_nyq, sprintf('  \\sigma_f = %.0f Hz', ...
     cfg.sigma_f_track*opt.F_trk/opt.F_nyq), 'FontSize', 8, ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Tag', 'fixMarker');
hold(ax, 'off');
fprintf('Fig 4(a): volume vs Eve-C read at M = %.1f -> %.2f kbit (fixed operating point)\n', ...
        Mres_nyq, B_nyq/1e3);

%% ------------------------------------------------------------------------
% FIG 3 SPLIT: each panel in its own figure, for its own section of Ch. 6
%   Fig 3a -> Section 6.1 (choice of the packet separation)
%   Fig 3b -> Section 6.9 (sensitivity to the optical depth)
% -------------------------------------------------------------------------
delete(findobj('Type', 'figure', '-regexp', 'Name', '^Fig 3[ab] '));   % rerun-safe
panelToFigure(findPanel(f3, '(a)'), 'Fig 3a Operating point optimisation');
panelToFigure(findPanel(f3, '(b)'), 'Fig 3b Imprint depth sensitivity');
f3.Visible = 'off';                          % the combined figure is no longer needed

%% ------------------------------------------------------------------------
% FIG 1 SPLIT
%   Fig 1ab -> Section 6.2 (signal and receiver output), panels (a) and (b)
%   Fig 1c  -> Section 6.3 (Campaign A)
% Panel (a) has two y-axes (yyaxis), which copyobj does not copy reliably,
% so (a) and (b) are redrawn from the workspace with the simulator's own
% code; panel (c) is a plain axes and is copied.
% -------------------------------------------------------------------------
f1 = findobj('Type', 'figure', 'Name', 'Fig 1 Spectral fringe and matched filtering');
delete(findobj('Type', 'figure', '-regexp', 'Name', '^Fig 1(ab|c) '));   % rerun-safe

f1ab = figure('Name', 'Fig 1ab Spectral fringe and matched filtering', ...
              'Color', 'w', 'Position', [80 100 900 380]);
% (a) spectra
subplot(1,2,1);
plot(f_MHz, sig.S_cl, 'LineWidth', 1.4); hold on;
plot(f_MHz, S_bit0, 'LineWidth', 0.9);
yyaxis right;
plot(f_MHz, S_bit0 - sig.S_cl, 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);
ylabel('Interference term [a.u.]');
yyaxis left;
xlabel('Detuning [MHz]');
ylabel('Modulation Depth [a.u.]');
title(sprintf('(a) Optical Spectrum (OD = %.2f)', cfg.OD));
legend({'S_{classical}','S_{quantum} (bit 0)','difference (right axis)'}, ...
       'Location', 'best', 'FontSize', 7.5);
grid on;  xlim([-1 1]*lineWidth_Hz/1e6);
% (b) matched filter output
subplot(1,2,2);
uax   = (1:cfg.N) - cfg.N/2 - 1;
Y0_mf = real(fftshift(fft(ifftshift( sig.tmpl.*weight))));
Y1_mf = real(fftshift(fft(ifftshift(-sig.tmpl.*weight))));
plot(uax, Y0_mf/max(abs(Y0_mf)), 'b-',  'LineWidth', 1.4); hold on;
plot(uax, Y1_mf/max(abs(Y0_mf)), 'r--', 'LineWidth', 1.4);
yline(0, 'k:');
xlabel('Baseband Bin Index');
ylabel('Normalized Output');
title('(b) Compressed Matched Filter Peak');
legend({'bit 0 (matched)','bit 1 (matched)'}, 'Location', 'best', 'FontSize', 8);
grid on;  xlim([-60 120]);       % wide enough to show the image peak at +54 bins

panelToFigure(findPanel(f1, '(c)'), 'Fig 1c Error rate validation');
f1.Visible = 'off';                          % the combined figure is no longer needed

%% ------------------------------------------------------------------------
% FIG 2 SPLIT
%   Fig 2a -> Section 6.4 (Campaign B: error rate and capacity along the pass)
%   Fig 2b -> Section 6.5 (cumulative volume and stop of each adversary)
% Panel (a) has two y-axes (yyaxis) and is redrawn with the simulator's own
% code; panel (b), already redrawn above, is a plain axes and is copied.
% -------------------------------------------------------------------------
delete(findobj('Type', 'figure', '-regexp', 'Name', '^Fig 2[ab] '));   % rerun-safe

figure('Name', 'Fig 2a Channel quality over the pass', 'Color', 'w', ...
       'Position', [110 80 620 460]);
yyaxis left;
errorbar(t_ep, BER_B, BER_B - BER_B_lo, BER_B_hi - BER_B, 'o-', ...
         'Color', [0 0.45 0.74], 'LineWidth', 1.3, 'MarkerSize', 4);
ylabel('Raw BER');  ylim([0 0.55]);
ax = gca;  ax.YColor = [0 0.45 0.74];
yyaxis right;
plot(t_ep, C_ep_sim, 's--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2, 'MarkerSize', 4);
ylabel('Capacity [bit/symbol]');  ylim([0 1.05]);
ax.YColor = [0.85 0.33 0.1];
xline(t_sec(idx_cross), 'k:', 'Closest approach', 'LabelOrientation', 'horizontal');
xlabel('Time from Pass Start [s]');
title('Channel Quality Over Fly-By Pass');
grid on;

panelToFigure(findPanel(f2, '(b)'), 'Fig 2b Covert information volume');
f2.Visible = 'off';                          % the combined figure is no longer needed

%% ------------------------------------------------------------------------
% SHOW: bring the four figures to the front (nothing is saved)
% -------------------------------------------------------------------------
figNames = {'Fig 1ab Spectral fringe and matched filtering', ...
            'Fig 1c Error rate validation', ...
            'Fig 2a Channel quality over the pass', ...
            'Fig 2b Covert information volume', ...
            'Fig 3a Operating point optimisation', ...
            'Fig 3b Imprint depth sensitivity', ...
            'Fig 4 Hopping alphabet'};
for i = 1:numel(figNames)
    f = findobj('Type', 'figure', 'Name', figNames{i});
    f.Visible = 'on';
    figure(f);                                   % raise the window
    fprintf('%s  ->  save as %s.pdf\n', figNames{i}, strrep(figNames{i}, ' ', '_'));
end

%% ------------------------------------------------------------------------
function ax = findPanel(f, tag)
% Axes of figure f whose title starts with tag, e.g. '(b)'.
    if isempty(f), error('Figure not found: run NewCovertSimFix.m first.'); end
    axs = findobj(f, 'Type', 'axes');
    for a = axs.'
        if startsWith(strjoin(string(a.Title.String)), tag)
            ax = a;  return;
        end
    end
    error('Panel %s not found.', tag);
end

function fNew = panelToFigure(ax, name)
% Copy one panel (with its legend, if any) into a new single-panel figure.
% The "(a) "/"(b) " prefix is removed from the title.
    fNew = figure('Name', name, 'Color', 'w', 'Position', [120 120 620 460]);
    lg = ax.Legend;
    if ~isempty(lg)
        h = copyobj([lg, ax], fNew);  axNew = h(2);
    else
        axNew = copyobj(ax, fNew);
    end
    axNew.Units = 'normalized';
    axNew.OuterPosition = [0 0 1 1];
    axNew.Title.String = regexprep(char(strjoin(string(axNew.Title.String))), '^\(\w\)\s*', '');
end