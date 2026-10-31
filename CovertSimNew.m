% =========================================================================
% QUANTUM-COVERT INTER-SATELLITE LINK SIMULATOR
% Feasibility Study of a Quantum Doppler Effect-Based Channel for Covert
% Inter-Satellite Communications - Master Thesis, Universita' di Padova.
%
% Spectral fringe encoding on the quantum Doppler effect, matched-filter
% reception and covert-volume assessment over a LEO-LEO pass.
%
% PHYSICAL BASIS
%   P. T. Grochowski, A. R. H. Smith, A. Dragan, K. Debski, "Quantum time
%   dilation in atomic spectra", Phys. Rev. Research 3, 023053 (2021).
%   Equation numbers quoted in the code refer to that work. Taken from it:
%   the coherent superposition of two momentum packets [Eq. A1], the
%   incoherent mixture [Eq. A2], the first-order Doppler line shape
%   [Eq. 26] and the centroid correction delta_Q [Eq. 32]. Added here: a
%   relative spatial displacement Delta_x between the branches, which turns
%   the interference lobe into a periodic fringe.
%
% PERIMETER
%   Channel simulator. The emitter is described by three declared inputs
%   (Section 4) and is not designed here:
%     Delta    velocity spread of the two momentum packets, in units of v/c
%     Delta_x  spatial separation between the branches
%     OD       rms optical depth of the imprint on the probe beam
%   The classical Doppler shift is assumed compensated from the ephemerides;
%   only the residual offset is simulated.
%
% MODEL CHAIN
%   1. Coherent superposition of two momentum packets [Eq. A1], relative
%      phase phi in {0, pi}: the bit.
%   2. The relative displacement Delta_x gives a phase linear in p; through
%      dw = Omega p/(mc) it becomes Phi(dw) = kappa dw, with
%      kappa = m c Delta_x/(hbar Omega): a fringe of spacing 1/kappa hertz.
%   3. Imprinted spectrum S(w) = S_class(w) + I_env(w) cos(Phi(w) + phi).
%   4. Receiver: subtracts S_class, weights by I_env, multiplies by
%      exp(+i kappa dw) and reads the sign of the zero-frequency FFT bin.
%   5. Fringe hopping: kappa is redrawn per symbol from a keyed sequence;
%      Eve-A (energy), Eve-B (key-aware) and Eve-C (scheme-aware,
%      key-unaware) are evaluated against it.
%
% BLOCKS
%   1 constants | 2 orbital geometry | 3 link budget | 4 signal model and
%   operating point | 5 Campaign A | 6 Campaign B | 7 sensitivity analyses |
%   8 capacity and covert volume | 9 figures
%
% SWITCHES
%   cfg.hopping = false  fixed-kappa baseline
%   cfg.fastMC  = false  full spectral Monte Carlo (slow, used as a check)
%
% REQUIRED TOOLBOXES
%   Satellite Communications (satelliteScenario, satellite, states)
%   Phased Array System (fspl)
%   Signal Processing (findpeaks)
%
% FIGURES
%   Block 9 draws the seven figures of Chapter 6 in their final form.
%
% Repository: https://github.com/giuliomora/QuantumSatelliteCovertSIM
% =========================================================================

clc; clear; close all;

SEED_MASTER = 42;
rng(SEED_MASTER, 'twister');

%% ========================================================================
% 1. PHYSICAL CONSTANTS
% =========================================================================
c_light    = 299792458;
lambda_0   = 267e-9;                % 27Al+ 1S0-3P1 intercombination line
f0_optical = c_light/lambda_0;
Omega_rad  = 2*pi*f0_optical;
h_planck   = 6.62607015e-34;
hbar       = 1.054571817e-34;
E_photon   = h_planck*f0_optical;
u_amu      = 1.66053906660e-27;
m_Al27     = 27*u_amu;
R_earth    = 6378137;
mu_earth   = 3.986004418e14;
Gamma0_Hz  = 520;                   % natural linewidth of the transition

fprintf('=== 1. CONSTANTS ===\n');
fprintf('27Al+ 1S0-3P1: lambda = %.1f nm, f0 = %.3f THz, Gamma0 = %d Hz\n', ...
        lambda_0*1e9, f0_optical/1e12, Gamma0_Hz);

%% ========================================================================
% 2. ORBITAL GEOMETRY (SGP4 PROPAGATION, J2-CONSISTENT INITIALISATION)
% =========================================================================
% Two co-altitude LEO satellites on planes 60 deg apart. The RAAN offset and
% the true anomalies are solved for numerically so that the closest approach
% is about 300 km at mid-window; SGP4 then propagates the resulting elements.
% The classical Doppler shift is computed and reported but never applied to
% the spectrum: its compensation is assumed.
fprintf('\n=== 2. ORBITAL GEOMETRY ===\n');

startTime  = datetime(2026,8,14,10,0,0);
stopTime   = startTime + minutes(15);
sampleTime = 0.1;
T_window   = seconds(stopTime - startTime);

altitude_orbit = 800e3;  sma = R_earth + altitude_orbit;
ecc = 0.0001; incl_deg = 60; argPeri_deg = 0; J2 = 1.08263e-3;
mu_n = sqrt(mu_earth/sma^3);
raan_rate = -1.5*mu_n*J2*(R_earth/(sma*(1-ecc^2)))^2*cosd(incl_deg);
planeAngle_deg = 60;     % angle between the two orbital planes
% Spherical trigonometry: cos(planeAngle) = cos^2 i + sin^2 i cos(dRAAN).
dRAAN_deg = acosd((cosd(planeAngle_deg)-cosd(incl_deg)^2)/sind(incl_deg)^2);
raan_Alice = 0; raan_Bob = dRAAN_deg;
T_orbit = 2*pi/mu_n; incl_rad = deg2rad(incl_deg);

posC = @(u0,raan0,t) sma*[ ...
  cos(raan0+raan_rate*t).*cos(u0+mu_n*t)-sin(raan0+raan_rate*t).*sin(u0+mu_n*t).*cos(incl_rad); ...
  sin(raan0+raan_rate*t).*cos(u0+mu_n*t)+cos(raan0+raan_rate*t).*sin(u0+mu_n*t).*cos(incl_rad); ...
  sin(u0+mu_n*t).*sin(incl_rad)];

% Coarse then fine search for the argument of latitude giving a 300 km pass.
% The analytic model is used only to choose the initial elements; SGP4 is
% the sole source of the geometry downstream.
tc = 0:2:T_orbit; pA = posC(0,0,tc); u0c = deg2rad(0:0.5:359.5);
bestOff = Inf; bestU = 0; bestT = 0;
for k = 1:numel(u0c)
    dk = vecnorm(posC(u0c(k),deg2rad(raan_Bob),tc)-pA,2,1);
    [dm,ik] = min(dk); off = abs(dm-300e3);
    if off < bestOff, bestOff = off; bestU = u0c(k); bestT = tc(ik); end
end
tf = max(0,bestT-10):0.02:min(T_orbit,bestT+10); pAf = posC(0,0,tf);
u0f = bestU + deg2rad(-1:0.002:1); bestOff = Inf;
for k = 1:numel(u0f)
    dk = vecnorm(posC(u0f(k),deg2rad(raan_Bob),tf)-pAf,2,1);
    [dm,ik] = min(dk); off = abs(dm-300e3);
    if off < bestOff, bestOff = off; bestU = u0f(k); bestT = tf(ik); end
end
trueAnom_Alice = mod(rad2deg(mu_n*bestT)-rad2deg(mu_n*T_window/2),360);
trueAnom_Bob   = mod(trueAnom_Alice+rad2deg(bestU),360);

sc   = satelliteScenario(startTime, stopTime, sampleTime);
satA = satellite(sc,sma,ecc,incl_deg,raan_Alice,argPeri_deg,trueAnom_Alice, ...
                 "OrbitPropagator","sgp4","Name","Alice");
satB = satellite(sc,sma,ecc,incl_deg,raan_Bob,argPeri_deg,trueAnom_Bob, ...
                 "OrbitPropagator","sgp4","Name","Bob");
[posAlice,velAlice,tS] = states(satA,"CoordinateFrame","inertial");
[posBob,  velBob,  ~ ] = states(satB,"CoordinateFrame","inertial");

t_sec = seconds(tS-tS(1)); N_t = numel(t_sec);
% Sign convention: uLOS points from Alice to Bob, so a receding Bob gives
% v_rel > 0 and a redshift, f_D < 0.
dPos  = posBob-posAlice; d_t = vecnorm(dPos,2,1);
uLOS  = dPos./d_t; v_rel_t = sum((velBob-velAlice).*uLOS,1);
f_D_t = -(v_rel_t/c_light)*f0_optical;
dfD_dt_t = gradient(f_D_t, sampleTime);

% Line-of-sight clearance above a 80 km atmospheric shell.
r_min = zeros(1,N_t);
for i = 1:N_t
    r_min(i) = norm(cross(posAlice(:,i),posBob(:,i)))/norm(posBob(:,i)-posAlice(:,i));
end
los_margin = r_min-(R_earth+80e3);

% Solar exclusion: Sun at least 30 deg from Bob's line of sight. The solar
% position uses the mean longitude only (error up to about 2 deg).
AU = 1.495978707e11; eps_ecl = deg2rad(23.439291);
nd = days(tS-datetime(2000,1,1,12,0,0,'TimeZone','UTC'));
lam_s = deg2rad(mod(280.460+0.9856474*nd,360));
posSun = AU*[cos(lam_s); sin(lam_s)*cos(eps_ecl); sin(lam_s)*sin(eps_ecl)];
uSun = (posSun-posBob)./vecnorm(posSun-posBob,2,1);
theta_sun = rad2deg(acos(min(max(sum((-uLOS).*uSun,1),-1),1)));
link_available = (los_margin>=0) & (theta_sun>=30);

[d_min, idx_cross] = min(d_t);
fprintf('Closest approach: %.2f km at t = %.1f s\n', d_min/1000, t_sec(idx_cross));
fprintf('Peak optical Doppler: %.3f GHz | Doppler rate at closest approach: %.3f MHz/s\n', ...
        max(abs(f_D_t))/1e9, abs(dfD_dt_t(idx_cross))/1e6);
fprintf('Link availability (LOS and solar exclusion): %.1f %%\n', 100*mean(link_available));

% --- Checks on the orbital block -----------------------------------------
v_orb   = sqrt(mu_earth/sma);                       % circular orbital speed
v_rel_max = max(abs(v_rel_t));
fD_check  = (v_rel_max/c_light)*f0_optical;         % from the velocity alone
sep_geo   = 2*sma*sind(planeAngle_deg/2);           % max plane separation
fprintf('  check: plane angle %.0f deg -> RAAN offset %.2f deg\n', ...
        planeAngle_deg, dRAAN_deg);
fprintf('  check: orbital speed %.3f km/s, peak LOS relative speed %.3f km/s (%.0f %% of it)\n', ...
        v_orb/1e3, v_rel_max/1e3, 100*v_rel_max/v_orb);
fprintf('  check: Doppler from velocity = %.3f GHz vs reported %.3f GHz (delta %.2e)\n', ...
        fD_check/1e9, max(abs(f_D_t))/1e9, abs(fD_check-max(abs(f_D_t)))/fD_check);
fprintf('  check: closest approach %.1f km against %.0f km target, %.0f km plane separation\n', ...
        d_min/1e3, 300, sep_geo/1e3);
if abs(d_min-300e3) > 50e3
    warning('Closest approach is far from the 300 km target: the element search did not converge.');
end

%% ========================================================================
% 3. FREE-SPACE OPTICAL LINK BUDGET
% =========================================================================
% The beam is expanded to theta_div = 8 urad (diffraction limit of the
% 0.20 m aperture: 1.7 urad) to tolerate 1.5 urad of pointing jitter; the
% transmit gain is that of the expanded beam. theta_div is the FULL 1/e^2
% divergence, so G = 8/(theta_div/2)^2. Pointing jitter is applied
% independently at both terminals.
fprintf('\n=== 3. OPTICAL LINK BUDGET ===\n');
D_tx = 0.20; D_rx = 0.30; P_tx = 100e-3;
eta_tx = 0.85; eta_rx = 0.80; eta_det = 0.65; rho_ap = 0.65;
sigma_jit_pat = 1.5e-6; theta_div = 8.0e-6;

G_tx = 8/(theta_div/2)^2;
G_rx = rho_ap*(pi*D_rx/lambda_0)^2;
FSPL_t = 10.^(fspl(d_t,lambda_0)/10);
patLoss = @(n) exp(-8*((sqrt((sigma_jit_pat*randn(1,n)).^2 + ...
                            (sigma_jit_pat*randn(1,n)).^2))./theta_div).^2);
Lp_t = patLoss(N_t).*patLoss(N_t);
P_rx_ideal_t = P_tx*eta_tx*eta_rx*(G_tx*G_rx)./FSPL_t;
P_rx_t = P_rx_ideal_t.*Lp_t;
fprintf('Received power at closest approach: %.2f dBm\n', 10*log10(P_rx_t(idx_cross)*1e3));

% --- Checks on the link budget -------------------------------------------
fprintf('  budget at closest approach [dB]: Ptx %.1f + Gtx %.1f + Grx %.1f - FSPL %.1f + eta %.1f = %.2f dBm\n', ...
        10*log10(P_tx*1e3), 10*log10(G_tx), 10*log10(G_rx), ...
        10*log10(FSPL_t(idx_cross)), 10*log10(eta_tx*eta_rx), ...
        10*log10(P_rx_ideal_t(idx_cross)*1e3));
fprintf('  pointing loss at that epoch: %.2f dB (mean over the pass %.2f dB)\n', ...
        -10*log10(Lp_t(idx_cross)), -10*log10(mean(Lp_t)));

% Cost of the beam expansion against the diffraction-limited aperture gain.
G_ap        = rho_ap*(pi*D_tx/lambda_0)^2;
theta_diff  = 2*lambda_0/(pi*D_tx/2);
fprintf('  BEAM: expanded to %.1f urad (gain %.1f dB) against a diffraction limit of %.2f urad (gain %.1f dB)\n', ...
        theta_div*1e6, 10*log10(G_tx), theta_diff*1e6, 10*log10(G_ap));
fprintf('  BEAM: cost of the expansion %.1f dB, bought against a pointing jitter of %.2f urad\n', ...
        10*log10(G_ap/G_tx), sigma_jit_pat*1e6);
if theta_div < 3*sigma_jit_pat
    warning('Beam width is within three times the pointing jitter: the loss model is not reliable.');
end

% Pointing-loss statistics over the pass.
Lp_sorted = sort(Lp_t);
fprintf('  pointing loss over the pass: mean %.2f dB, median %.2f dB, 5th percentile %.2f dB\n', ...
        -10*log10(mean(Lp_t)), -10*log10(median(Lp_t)), ...
        -10*log10(Lp_sorted(max(1,round(0.05*N_t)))));
% The same theta_div is used at both terminals (conservative simplification).

%% ========================================================================
% 4. QUANTUM STATE, SPECTRAL FRINGE, IMPRINTING AND OPERATING POINT
% =========================================================================
fprintf('\n=== 4. QUANTUM STATE AND OPERATING POINT ===\n');
cfg = struct();

% --- 4.1 Source: momentum spread -----------------------------------------
% Delta is declared. The expression is only a convenient way to obtain a
% plausible value, not a derivation. [P8] Delta is the AMPLITUDE width
% parameter, phi(p) ~ exp(-(p-p0)^2/(2 Delta^2)): it equals the ground-state
% spread of |phi|^2 at f_source, so each lobe (standard deviation Delta/sqrt(2))
% corresponds to a trap of f_source/2 = 0.5 MHz (printed in Section 4.14).
cfg.f_source = 1e6;
cfg.Delta = sqrt(hbar*(2*pi*cfg.f_source)/(2*m_Al27))/c_light;   % v/c

% --- 4.2 Superposition weight --------------------------------------------
% The interference term scales as sin(2 theta), the centroid correction
% delta_Q of Eq. (32) as sin(4 theta): theta = pi/4 maximises the signal and
% nulls delta_Q.
cfg.theta = pi/4;

% --- 4.3 Branch separation and fringe -------------------------------------
% Phi(dw) = kappa dw,  kappa = m c Delta_x/(hbar Omega)  (compression parameter)
cfg.Delta_x = 90e-9;
cfg.kappa   = m_Al27*c_light*cfg.Delta_x/(hbar*Omega_rad);

% --- 4.3b Fringe hopping --------------------------------------------------
% With hopping, Delta_x (hence kappa) is redrawn every symbol from a keyed
% sequence. The sequence is not generated: only the size of the alphabet
% matters, bounded in operatingPoint by Nyquist sampling, image clearance and
% Doppler tracking residual.
cfg.hopping       = true;   % false restores the fixed-kappa baseline
cfg.hopLossTol_dB = 1.0;    % loss budget fixing the tracking ceiling on kappa
cfg.imageMargin   = 3.0;    % image peak must clear this many compressed-peak
                            % half-widths (same guard as Sections 4.12, 7.3)

% --- 4.4 Imprint depth ----------------------------------------------------
% rms optical depth over the band; the signal lives in the optically thin
% valley between the line peaks, where the linear imprint is accurate.
% Swept in Section 7.2.
cfg.OD = 0.8;

% --- 4.5 Residual errors and receiver -------------------------------------
cfg.eps_sub = 1e-3;   % per-bin relative error on the known classical profile,
                      % frozen over a frame (a bias, not noise)
cfg.eps_kappa = 1e-3; % relative error on kappa; displaces the compressed
                      % peak by F*eps bins
cfg.eps_theta_probe = 1e-2;   % theta miscalibration used to linearise the
                              % centroid leakage (Section 4.10)
cfg.F_excess     = 3.0;    % [P9] factor on the noise AMPLITUDE: 9x in power (9.5 dB)
                           % = 3 dB heterodyne penalty + 6.5 dB margin for laser
                           % intensity noise; 1.0 is the pure shot-noise limit
% [P13] Structured (smooth) calibration errors on Bob's copy of S_cl, drawn
% once per frame as Gaussian with these standard deviations: height (gain,
% relative), centre (Hz) and lobe width (relative, area kept, centres fixed).
% With cfg.profileFit = true Bob estimates the three quantities from the
% average of the frame and removes them before deciding (Section 4.2 of the
% thesis); false reproduces a receiver without the estimation step.
cfg.eps_gain   = 1e-3;
cfg.eps_c_Hz   = 1e2;
cfg.eps_w      = 1e-3;
cfg.profileFit = true;
cfg.sigma_f_track = 1e3;   % residual Doppler offset error, 1-sigma, in Hz
cfg.wellDepth = 1e6;       % photoelectrons per bin per symbol handled linearly
                           % by the receiver chain (dynamic range check, 4.8)

cfg.N = 2^14;
SNR_thr_Eve = 3;           % detection threshold of the eavesdropper
thr_grid    = [3 2 1 0.5 0.25];   % [P1][P2] thresholds reported in Section 8
                                  % (same operating point; 0.25 ~ criterion of
                                  % the laboratory experiments, error >= 0.45)

%% ------------------------------------------------------------------------
% 4.6 Operating-point optimisation
% -------------------------------------------------------------------------
% The packet separation is swept and the value maximising the covert volume
% against the target adversary is kept:
%     B(sep) = sum_m n_m C(d_m), truncated where sqrt(sum_m n_m delta_m^2)
%     reaches the detection threshold.
fprintf('\n--- 4.6 Operating-point optimisation ---\n');

N_epochs = 21;
idx_ep   = round(linspace(1,N_t,N_epochs));
t_ep     = t_sec(idx_ep);
% Equal share of the window per epoch; availability at the epoch centre.
dt_ep    = (T_window/N_epochs)*ones(1,N_epochs);
avail_ep = double(link_available(idx_ep));

% Design power: received power times the pass-mean pointing loss, the same
% mean power delivered by the Monte Carlo campaigns.
L_point_mean = mean(Lp_t);
fprintf('Design power includes the pass-mean pointing loss: %.2f dB\n', -10*log10(L_point_mean));
% [P10] Eve's evidence is D^2 = sum delta^2 with delta ~ d^2 ~ 1/sigma^2 ~ L,
% so it averages L^2, not L: her effective power is P*sqrt(<L^2>). She can
% weight each frame by its fade because the probe is overt and bright (she
% measures its power). Bob keeps <L>. eveGain multiplies the adversaries' d^2.
L_eve_eff = sqrt(mean(Lp_t.^2));
eveGain   = L_eve_eff/L_point_mean;
mu_jit    = 16*sigma_jit_pat^2/theta_div^2;              % per terminal, analytic
eveGain_an = sqrt((1/(1+2*mu_jit))^2)/(1/(1+mu_jit))^2;
fprintf('  [P10] adversaries: effective power sqrt(<L^2>) -> d^2 x%.3f, D^2 x%.3f (analytic x%.3f, x%.3f)\n', ...
        eveGain, eveGain^2, eveGain_an, eveGain_an^2);

par = struct('theta',cfg.theta,'Omega_rad',Omega_rad, ...
             'eps_theta_probe',cfg.eps_theta_probe, ...
             'f0',f0_optical,'N',cfg.N,'Delta',cfg.Delta,'kappa',cfg.kappa, ...
             'OD',cfg.OD,'Gamma0',Gamma0_Hz, ...
             'F_excess',cfg.F_excess,'E_photon',E_photon,'eta_det',eta_det, ...
             'P_ep',P_rx_ideal_t(idx_ep)*L_point_mean,'dt_ep',dt_ep.*avail_ep, ...
             'SNR_thr',SNR_thr_Eve,'m',m_Al27,'c',c_light,'hbar',hbar, ...
             'hopping',cfg.hopping,'sigma_f',cfg.sigma_f_track, ...
             'hopLossTol_dB',cfg.hopLossTol_dB,'imageMargin',cfg.imageMargin, ...
             'eveGain',eveGain);

sep_vec = 4.0:0.05:9.5;
nSep = numel(sep_vec);
B_naive = zeros(1,nSep);   % covert volume against Eve-A
B_keyed = zeros(1,nSep);   % covert volume against Eve-B
B_hop   = zeros(1,nSep);   % covert volume against Eve-C
B_noEve = zeros(1,nSep);   % volume with no eavesdropper
for k = 1:nSep
    o = operatingPoint(sep_vec(k), par);
    if o.nyquistPhase < pi                 % reject an undersampled fringe
        B_naive(k) = o.volume_naive;
        B_keyed(k) = o.volume_keyed;
        B_hop(k)   = o.volume_hopped;
        B_noEve(k) = o.volume_noEve;
    end
end
% Target adversary: Eve-C with hopping, Eve-A without (against Eve-B the link
% is not covert for any separation).
if cfg.hopping
    [B_best, k_best] = max(B_hop);
    targetName = 'Eve-C (scheme-aware, key-unaware)';
else
    [B_best, k_best] = max(B_naive);
    targetName = 'Eve-A (energy detector)';
end
if B_best <= 0
    error('No feasible separation: the fringe is undersampled over the whole sweep.');
end
cfg.sep_over_Delta = sep_vec(k_best);
opt = operatingPoint(cfg.sep_over_Delta, par);

fprintf('Optimisation target: %s\n', targetName);
fprintf('Optimal separation = %.2f Delta -> covert volume = %.3e bit\n', ...
        cfg.sep_over_Delta, B_best);
fprintf('  same point, no eavesdropper        = %.3e bit\n', opt.volume_noEve);
fprintf('  same point, energy detector Eve-A  = %.3e bit\n', opt.volume_naive);
fprintf('  same point, scheme-aware Eve-C     = %.3e bit\n', opt.volume_hopped);
fprintf('  same point, key-aware eavesdropper = %.3e bit\n', opt.volume_keyed);
fprintf('  hopping alphabet: %d integer fringe counts (F in [%d, %.1f]) -> %.1f resolvable pitches\n', ...
        opt.M_alpha, opt.F_min, min(opt.F_nyq,opt.F_trk), opt.M_res);
% Design deflection (pass-mean fading), not the Monte Carlo value.
fprintf('  median deflection d = %.3f (design, mean fading) -> raw BER = %.3f, capacity = %.4f bit/symbol\n', ...
        median(opt.d_ep), median(opt.BER_ep), median(opt.C_ep));

% --- Checks on the optimisation ------------------------------------------
% Closed-form covert optimum for equal epochs, small d:
%     d_opt = (2 M_eff thr^2 / M_pass)^(1/4)
% with M_eff the cells of the targeted adversary. Agreement with the sweep is
% expected in order of magnitude only (unequal epochs, early truncation).
d_opt = (2*opt.M_eff_target*SNR_thr_Eve^2/opt.M_pass)^(1/4);
fprintf('  check: analytic covert optimum d_opt = %.3f (equal epochs), sweep median %.3f (ratio %.2f)\n', ...
        d_opt, median(opt.d_ep), median(opt.d_ep)/d_opt);
% [P10] The closed form assumes that Eve sees Bob's d; with the fading
% average of the adversaries she sees d*sqrt(eveGain).
fprintf('  [P10] same check on the deflection seen by the adversaries: median %.3f (ratio %.2f)\n', ...
        median(opt.d_ep)*sqrt(eveGain), median(opt.d_ep)*sqrt(eveGain)/d_opt);
% The optimum must be interior to the sweep.
if k_best == 1 || k_best == nSep
    warning('Optimal separation sits on the edge of sep_vec: widen the sweep.');
else
    fprintf('  check: optimum interior to the sweep (%.2f in [%.2f, %.2f])\n', ...
            cfg.sep_over_Delta, sep_vec(1), sep_vec(end));
end
% Binding constraint: detection or pass duration.
B_target_now = opt.volume_naive;
if cfg.hopping, B_target_now = opt.volume_hopped; end
if B_target_now >= 0.999*opt.volume_noEve
    fprintf('  check: BINDING CONSTRAINT IS PASS DURATION - the eavesdropper never reaches threshold,\n');
    fprintf('         so the covert model is inert at this operating point.\n');
else
    fprintf('  check: binding constraint is DETECTION - %.1f %% of the transmissible volume is usable\n', ...
            100*B_target_now/opt.volume_noEve);
end

%% ------------------------------------------------------------------------
% 4.7 Frozen operating point and constraints on the fringe
% -------------------------------------------------------------------------
cfg.p1     = -cfg.sep_over_Delta*cfg.Delta/2;
cfg.p2     = +cfg.sep_over_Delta*cfg.Delta/2;
cfg.B_IF   = opt.B_IF;
cfg.T_sym  = opt.T_sym;
cfg.dw     = opt.dw;
cfg.Phi    = opt.Phi;
cfg.key    = exp(+1i*cfg.kappa*cfg.dw);

% In-phase template (matched filter) and quadrature template (used to apply
% a fringe-phase error without rebuilding the spectrum).
sig = struct('S_cl', opt.S_cl, 'I_env', opt.I_env, ...
             'tmpl',   opt.I_env.*cos(opt.Phi), ...
             'tmpl_q', opt.I_env.*sin(opt.Phi));

R_b_raw   = 1/cfg.T_sym;
dw_max    = 2*pi*cfg.B_IF/2;
dw_step   = cfg.dw(2)-cfg.dw(1);
lineWidth_Hz  = f0_optical*(cfg.sep_over_Delta+6)*cfg.Delta;
sigma_env_rad = (cfg.Delta/sqrt(2))*Omega_rad;

fprintf('\n--- 4.7 Frozen operating point ---\n');
fprintf('Source: Delta = %.3e v/c -> velocity spread = %.4f m/s (declared)\n', ...
        cfg.Delta, cfg.Delta*c_light);
fprintf('Imprint depth OD = %.3f (declared; swept in Section 7.2)\n', cfg.OD);
fprintf('Packet separation = %.2f Delta -> velocity difference = %.4f m/s\n', ...
        cfg.sep_over_Delta, (cfg.p2-cfg.p1)*c_light);
fprintf('Derived B_IF = %.2f MHz -> T_sym = %.1f us -> raw bit rate = %.3f kbit/s\n', ...
        cfg.B_IF/1e6, cfg.T_sym*1e6, R_b_raw/1e3);
fringeSpacing_Hz = 1/cfg.kappa;
nFringe          = cfg.B_IF*cfg.kappa;
fprintf('Branch separation Delta_x = %.1f nm (declared)\n', cfg.Delta_x*1e9);
fprintf('kappa = %.4e s/rad -> fringe spacing = %.3f MHz -> %.1f fringes in the band\n', ...
        cfg.kappa, fringeSpacing_Hz/1e6, nFringe);
% Reference only: Delta_x that would give as many fringes as the packet
% separation suggests.
Delta_x_matched = cfg.Delta_x * (cfg.sep_over_Delta*(cfg.sep_over_Delta+6)/pi)/nFringe;
fprintf('  for reference, %.1f fringes would need Delta_x = %.1f nm\n', ...
        cfg.sep_over_Delta*(cfg.sep_over_Delta+6)/pi, Delta_x_matched*1e9);

% Ceiling 1: Nyquist sampling of the fringe on the spectral grid.
phasePerSample = cfg.kappa*dw_step;
fprintf('Ceiling 1 (Nyquist): fringe phase per sample = %.4f rad (limit pi), margin x%.0f\n', ...
        phasePerSample, pi/phasePerSample);

% Ceiling 2: fringe spacing against the natural linewidth.
fprintf('Ceiling 2 (natural linewidth): fringe spacing %.3f MHz vs Gamma0 = %d Hz, margin x%.0f\n', ...
        fringeSpacing_Hz/1e6, Gamma0_Hz, fringeSpacing_Hz/Gamma0_Hz);
fprintf('Ceiling 3 (coherence time of the superposition) is experimental and outside this model.\n');

% --- Hopping alphabet ----------------------------------------------------
% The three ceilings bounding the alphabet, reported side by side to show
% which one binds.
fprintf('\n--- 4.7b Fringe-hopping alphabet ---\n');
if cfg.hopping
    fprintf('Hopping ENABLED: kappa redrawn per symbol from the keyed sequence.\n');
else
    fprintf('Hopping DISABLED: fixed-kappa baseline, Eve-C degenerates to Eve-B.\n');
end
fprintf('Declared operating point sits at F = %.1f fringes in the band\n', opt.F_now);
fprintf('  ceiling (a) Nyquist                  : F <= %.1f (integers up to %d)\n', ...
        opt.F_nyq, floor(opt.F_nyq));
fprintf('  ceiling (b) image clearance          : F >= %d (compressed half-width %.1f bins)\n', ...
        opt.F_min, opt.halfW);
fprintf('  ceiling (c) tracking residual %.0f Hz  : F <= %.1f (integers up to %d, budget %.1f dB)\n', ...
        cfg.sigma_f_track, opt.F_trk, floor(opt.F_trk), cfg.hopLossTol_dB);
if opt.F_trk < opt.F_nyq
    fprintf('  BINDING CEILING: the Doppler tracking residual, not the fringe physics.\n');
    fprintf('    tightening sigma_f to %.0f Hz would open the full Nyquist alphabet\n', ...
            cfg.sigma_f_track*opt.F_trk/opt.F_nyq);
else
    fprintf('  BINDING CEILING: Nyquist sampling of the fringe.\n');
end
fprintf('Integer fringe counts in range: %d\n', opt.M_alpha);
fprintf('Resolvable pitches under the envelope: M = %.1f\n', opt.M_res);
fprintf('  template correlation for adjacent counts %.4f, for 10 counts %.3f, for 30 counts %.3f\n', ...
        exp(-(opt.rhoA)^2/2), exp(-(10*opt.rhoA)^2/2), exp(-(30*opt.rhoA)^2/2));
fprintf('  -> Eve-C spread over %.1f cells (Eve-A: %.0f, Eve-B: 1)\n', opt.M_eff, opt.N_env);
fprintf('  envelope occupancy N_env = %.0f samples: upper bound on the cells ANY keyed scheme can create\n', ...
        opt.N_env);
if cfg.hopping
    Mres_nyq = resolvablePitches(floor(opt.F_nyq)-opt.F_min+1, opt.rhoA);
    fprintf('  with the full Nyquist range the resolvable alphabet would be %.1f (x%.2f in M, x%.2f in volume\n', ...
            Mres_nyq, Mres_nyq/opt.M_res, sqrt(Mres_nyq/opt.M_res));
    fprintf('    by the asymptotic sqrt(M) scaling, not a simulated volume)\n');
end
if opt.F_now < opt.F_min || opt.F_now > min(opt.F_nyq,opt.F_trk)
    warning('The declared Delta_x lies outside the admissible hopping alphabet.');
end

% --- [P14] Alphabet used by Alice and meaning of M ----------------------------
% Alice draws the pitch uniformly from ALL the admissible integer counts
% (option A). Eve-C sums the energies of the matched filters of the same
% counts; with unit-energy templates of correlation matrix R, her statistic
% has mean shift d^2*sum_ij R_ij^2/M_alpha and standard deviation
% sqrt(2*sum_ij R_ij^2), i.e. delta_C = d^2/sqrt(2K) with
%     K = M_alpha^2 / sum_ij R_ij^2,
% the participation ratio of R. resolvablePitches evaluates it with the
% Gaussian-envelope correlation; here it is recomputed on the exact templates.
if cfg.hopping
    F_all = opt.F_min:floor(min(opt.F_nyq,opt.F_trk));
    Tall  = single(opt.I_env.*cos((F_all/opt.B_IF).*opt.dw));     % N x M_alpha
    Tall  = Tall./vecnorm(Tall,2,1);
    Rall  = double(Tall.'*Tall);
    K_exact = numel(F_all)^2/sum(Rall(:).^2);
    % Same quantity if Alice used only M_res equispaced counts (option B).
    F_eq  = unique(round(linspace(opt.F_min, F_all(end), max(1,round(opt.M_res)))));
    [~,iEq] = ismember(F_eq, F_all);
    Req   = Rall(iEq,iEq);
    K_eq  = numel(F_eq)^2/sum(Req(:).^2);
    clear Tall Rall
    fprintf('  [P14] Alice draws from all %d counts: Eve-C cells K = %.2f on the exact templates\n', ...
            numel(F_all), K_exact);
    fprintf('        (Gaussian-envelope approximation used in the volumes: %.2f, ratio %.3f)\n', ...
            opt.M_eff, K_exact/opt.M_eff);
    fprintf('        for reference, %d equispaced counts would give K = %.2f (delta_C x%.3f)\n', ...
            numel(F_eq), K_eq, sqrt(K_exact/K_eq));
    if abs(K_exact/opt.M_eff - 1) > 0.05
        warning('Exact participation ratio departs from the Gaussian approximation by more than 5 %%.');
    end
end

fprintf('Fringe-coefficient error: eps_kappa = %.1e -> peak displacement = %.3f bins (declared pitch), %.2f bins (finest pitch)\n', ...
        cfg.eps_kappa, cfg.eps_kappa*cfg.kappa*cfg.B_IF, ...
        cfg.eps_kappa*floor(min(opt.F_nyq,opt.F_trk)));

%% ------------------------------------------------------------------------
% 4.8 Detection noise and Doppler tracking residual
% -------------------------------------------------------------------------
% sigma = F_excess / sqrt(n_bin),  n_bin = P_rx T_sym eta_det/(E_photon N):
% the contrast does not attenuate with range, the photon count does.
fprintf('\n--- 4.8 Noise and tracking ---\n');
fprintf('Shot-derived noise (F_excess = %.1f): sigma = %.2e at closest approach, %.2e at pass edge\n', ...
        cfg.F_excess, min(opt.sigma_ep), max(opt.sigma_ep));

% A residual offset s rotates the compressed peak by kappa*s; the averaged
% loss for a zero-mean Gaussian offset is exp(-sigma_phi^2/2).
cfg.sigma_phi_track = cfg.kappa*2*pi*cfg.sigma_f_track;
fprintf('Tracking offset: sigma_f = %.0f Hz -> fringe phase error %.4f rad -> %.3f dB loss\n', ...
        cfg.sigma_f_track, cfg.sigma_phi_track, ...
        -20*log10(exp(-cfg.sigma_phi_track^2/2)));
fprintf('  offset tolerated for 1 dB of loss: %.1f kHz\n', ...
        sqrt(2*log(10^(1/20)))/(cfg.kappa*2*pi)/1e3);
% Same offset on the finest admissible pitch.
F_top = floor(min(opt.F_nyq, opt.F_trk));
fprintf('  same offset on the finest admissible pitch (F = %d): %.3f dB loss (budget %.1f dB)\n', ...
        F_top, -20*log10(exp(-(2*pi*cfg.sigma_f_track*F_top/cfg.B_IF)^2/2)), cfg.hopLossTol_dB);

% --- Residual Doppler rate: tolerance ------------------------------------
% A residual rate r sweeps the fringe during one symbol and scales it by
% sinc(x), x = pi*kappa*r*T_sym. The tolerated rate is computed on the finest
% pitch for the same loss budget: a requirement on the compensation, not a
% simulated residual.
lossTarget = 10^(-cfg.hopLossTol_dB/20);
x_tol   = fzero(@(x) sin(x)./x - lossTarget, [1e-6, pi-1e-6]);
kap_top = F_top/cfg.B_IF;
r_tol_top  = x_tol/(pi*kap_top*cfg.T_sym);
r_tol_decl = x_tol/(pi*cfg.kappa*cfg.T_sym);
cfg.rateTol_Hzs = r_tol_top;
fprintf('Tracking rate: tolerated residual for %.1f dB loss = %.2f MHz/s on the finest pitch (%.1f MHz/s on the declared one)\n', ...
        cfg.hopLossTol_dB, r_tol_top/1e6, r_tol_decl/1e6);
fprintf('  against a raw rate of %.1f MHz/s at closest approach: the compensation must remove %.2f %% of it\n', ...
        abs(dfD_dt_t(idx_cross))/1e6, 100*(1 - r_tol_top/abs(dfD_dt_t(idx_cross))));
fprintf('  uncompensated, the line would move %.0f kHz within one symbol against a finest pitch of %.1f kHz\n', ...
        abs(dfD_dt_t(idx_cross))*cfg.T_sym/1e3, cfg.B_IF/F_top/1e3);

% --- Checks on the noise model -------------------------------------------
% sigma recomputed from the photon count (with the pass-mean pointing loss),
% and its scaling with range.
n_bin_chk = (P_rx_ideal_t(idx_cross)*cfg.T_sym/E_photon)*eta_det/cfg.N;   % unfaded
sig_chk   = cfg.F_excess/sqrt(n_bin_chk*L_point_mean);
fprintf('  check: photons per bin at closest approach = %.2e unfaded, %.2e with mean fading -> sigma = %.2e (reported %.2e)\n', ...
        n_bin_chk, n_bin_chk*L_point_mean, sig_chk, min(opt.sigma_ep));
[~,i_far] = max(d_t(idx_ep));
fprintf('  check: sigma edge/closest = %.2f against range ratio %.2f (shot noise scales as range)\n', ...
        max(opt.sigma_ep)/min(opt.sigma_ep), d_t(idx_ep(i_far))/d_min);
fprintf('  photons per symbol at closest approach: %.2e unfaded, %.2e with mean fading, over %d bins\n', ...
        n_bin_chk*cfg.N, n_bin_chk*L_point_mean*cfg.N, cfg.N);

% --- Receiver dynamic range ----------------------------------------------
% Flux per bin at the brightest epoch (closest approach, unfaded) against the
% declared linear range of the receiver. The receiver is assumed not to
% saturate; this reports by how much that assumption is stretched.
n_bin_max = (max(P_rx_ideal_t)*cfg.T_sym/E_photon)*eta_det/cfg.N;
satRatio  = n_bin_max/cfg.wellDepth;
fprintf('  dynamic range: %.2e photoelectrons per bin per symbol vs declared linear range %.1e -> x%.0f\n', ...
        n_bin_max, cfg.wellDepth, satRatio);
if satRatio > 1
    fprintf('              absorbing it needs %.1f dB of attenuation, or an integration cut to %.2f %% of the symbol\n', ...
            10*log10(satRatio), 100/satRatio);
    fprintf('              if absorbed by attenuation alone the deflection would fall by x%.1f\n', ...
            sqrt(satRatio));
    warning(['Flux per bin exceeds the declared linear range by a factor %.0f: the ' ...
             'no-saturation assumption is not self-evident.'], satRatio);
else
    fprintf('              the receiver stays within the declared linear range; the assumption is inert here\n');
end

%% ------------------------------------------------------------------------
% 4.9 Adversary models
% -------------------------------------------------------------------------
% All adversaries receive with the same SNR as Bob (conservative).
%   Eve-A, key-unaware: envelope-weighted energy, delta_A = d^2/sqrt(2 N_env),
%          N_env = (sum I_env^2)^2/sum I_env^4.
%   Eve-B, key-aware (Kerckhoffs): same matched filter as Bob,
%          delta_B = d^2/sqrt(2). Key-compromise case.
%   Eve-C, scheme-aware, key-unaware: tests all M resolvable pitches
%          non-coherently, delta_C = d^2/sqrt(2 M).
% Ratios: delta_B/delta_A = sqrt(N_env), delta_B/delta_C = sqrt(M).
d_cross = interp1(t_ep, opt.d_ep, t_sec(idx_cross), 'linear', 'extrap');
d_cross_eve = d_cross*sqrt(eveGain);       % [P10] deflection seen by the adversaries
delta_A = d_cross_eve^2/sqrt(2*opt.N_env);
delta_B = d_cross_eve^2/sqrt(2);
delta_C = d_cross_eve^2/sqrt(2*opt.M_eff);
fprintf('\n--- 4.9 Adversary models (evaluated at closest approach) ---\n');
fprintf('Receiver deflection d = %.3f -> raw BER = %.4f\n', d_cross, 0.5*erfc(d_cross/sqrt(2)));
fprintf('  [P10] deflection seen by the adversaries (fading averaged on L^2): %.3f\n', d_cross_eve);
fprintf('Eve-A (energy)   delta = %.3e per symbol -> threshold after %.2e symbols\n', ...
        delta_A, (SNR_thr_Eve/delta_A)^2);
fprintf('Eve-C (hop-blind) delta = %.3e per symbol -> threshold after %.2e symbols\n', ...
        delta_C, (SNR_thr_Eve/delta_C)^2);
fprintf('Eve-B (matched)  delta = %.3e per symbol -> threshold after %.2e symbols\n', ...
        delta_B, (SNR_thr_Eve/delta_B)^2);
fprintf('Sensitivity ratio Eve-B/Eve-A = x%.1f (expected sqrt(N_env) = %.1f)\n', ...
        delta_B/delta_A, sqrt(opt.N_env));
fprintf('  hopping divides it by x%.2f (sqrt(M)): Eve-C/Eve-A = x%.2f\n', ...
        delta_B/max(delta_C,eps), delta_C/delta_A);
if abs(delta_B/delta_A - sqrt(opt.N_env)) > 1e-6*sqrt(opt.N_env)
    warning('Adversary ratio departs from sqrt(N_env): the detector models are inconsistent.');
end
if abs(delta_B/delta_C - sqrt(opt.M_eff)) > 1e-6*sqrt(opt.M_eff)
    warning('Eve-C ratio departs from sqrt(M): the hopping model is inconsistent.');
end
% Single-symbol pitch identification by Eve-C: she needs the matched-filter
% amplitude d to exceed the largest of M noise cells, ~sqrt(2 ln M).
% Evaluated on the unfaded power (the adversary's best case).
if cfg.hopping
    d_id   = sqrt(2*log(max(opt.M_eff,2)));
    d_pass = opt.d_ep;
    fprintf('  check: single-symbol pitch identification needs d >~ %.2f\n', d_id);
    d_cross_uf = d_cross/sqrt(L_point_mean);
    d_pass_uf  = d_pass/sqrt(L_point_mean);
    fprintf('         d at closest approach %.2f with mean fading, %.2f unfaded (margin x%.2f on the unfaded value)\n', ...
            d_cross, d_cross_uf, d_id/max(d_cross_uf,eps));
    fprintf('         pass median %.2f with mean fading, %.2f unfaded (margin x%.2f)\n', ...
            median(d_pass), median(d_pass_uf), d_id/max(median(d_pass_uf),eps));
    fprintf('         epochs above half that bound (unfaded): %.0f %% of the pass\n', ...
            100*mean(d_pass_uf > 0.5*d_id));
    if d_cross_uf > 0.5*d_id
        fprintf('         NOTE: the margin is thin only near closest approach. Capping the\n');
        fprintf('         per-symbol deflection at %.2f restores it and costs almost nothing,\n', 0.5*d_id);
        fprintf('         since capacity saturates in d while the detection statistic grows as\n');
        fprintf('         d^2: the strong epochs buy little and spend the covert budget fast.\n');
    end
end

% --- Checks on the adversary models --------------------------------------
% Eve-A measured by Monte Carlo on her energy statistic and compared with
% the closed form d^2/sqrt(2 N_env).
rsEve   = RandStream('mrg32k3a','Seed',SEED_MASTER); rsEve.Substream = 9;
nTrial  = 2000;
% [P10] The closed form delta_A above uses the deflection seen by the
% adversaries, d_cross_eve, so the measured statistic must use the same one.
s_unit  = sig.tmpl*(d_cross_eve/sqrt(sum(sig.tmpl.^2)));   % ||s_unit|| = d_cross_eve
E_null  = zeros(1,nTrial); E_sig = zeros(1,nTrial);
wA      = sig.I_env.^2/max(sig.I_env.^2);          % Eve-A weights (known envelope)
for it = 1:nTrial
    z = randn(rsEve,cfg.N,1);
    E_null(it) = sum(wA.*z.^2);                 % envelope-weighted energy
    E_sig(it)  = sum(wA.*(z+s_unit).^2);
end
delta_A_mc = (mean(E_sig)-mean(E_null))/std(E_null);
fprintf('  check: Eve-A measured on the energy statistic = %.3e vs closed form %.3e (ratio %.3f)\n', ...
        delta_A_mc, delta_A, delta_A_mc/delta_A);
if abs(delta_A_mc/delta_A - 1) > 0.15
    warning('Measured Eve-A deflection departs from the closed form by more than 15 %%.');
end
% If Eve-B crosses in less than one symbol, she detects the first symbol.
nB = (SNR_thr_Eve/delta_B)^2;
if nB < 1
    fprintf('  check: Eve-B crosses threshold within the FIRST symbol (%.3f), so the\n', nB);
    fprintf('         fractional figure above is a formal extrapolation, not a count.\n');
end
fprintf('  check: Eve-A tolerates %.2e symbols of %.2e in the pass (%.1f %%)\n', ...
        (SNR_thr_Eve/delta_A)^2, opt.M_pass, 100*(SNR_thr_Eve/delta_A)^2/opt.M_pass);

%% ------------------------------------------------------------------------
% 4.10 Leakage through the classical Doppler centroid
% -------------------------------------------------------------------------
% delta_Q vanishes at theta = pi/4; a miscalibration of theta makes the line
% centroid depend on the bit. Two estimates of the tolerance on theta are
% printed: a pessimistic one at the closest-approach photon count, and the
% epoch-by-epoch one computed in operatingPoint, which is the one used.
dp_pk  = cfg.p2 - cfg.p1;
dQ_fun = @(ph,th) cos(ph)*sin(4*th)*dp_pk / ...
                  (4*(cos(ph)*sin(2*th) + exp(dp_pk^2/(4*cfg.Delta^2))));
swing_exact_Hz = abs(dQ_fun(0,cfg.theta) - dQ_fun(pi,cfg.theta))*f0_optical;
lineWidthEff_Hz = f0_optical*cfg.Delta*(cfg.sep_over_Delta/2 + 1);
n_ph_cross = (P_rx_t(idx_cross)*cfg.T_sym/E_photon)*eta_det;
sigma_centroid_Hz = cfg.F_excess*lineWidthEff_Hz/sqrt(n_ph_cross);

eps_theta = 0.01;                    % probe miscalibration, rad
swing_eps_Hz = abs(dQ_fun(0,cfg.theta+eps_theta) - dQ_fun(pi,cfg.theta+eps_theta))*f0_optical;
snr_centroid_sym = (swing_eps_Hz/2)/sigma_centroid_Hz;
snr_centroid_acc = snr_centroid_sym*sqrt(opt.M_pass);
% delta_Q ~ 4 eps near pi/4: the tolerance scales linearly.
eps_theta_req = eps_theta*SNR_thr_Eve/snr_centroid_acc;

fprintf('\n--- 4.10 Centroid leakage (delta_Q, Eq. 32) ---\n');
fprintf('At theta = pi/4 exactly: centroid swing = %.2e Hz (numerically zero)\n', swing_exact_Hz);
fprintf('With %.3f rad of miscalibration: swing = %.3f Hz vs centroid precision %.2f Hz\n', ...
        eps_theta, swing_eps_Hz, sigma_centroid_Hz);
fprintf('  per symbol SNR = %.3f, accumulated over %.2e symbols = %.1f (threshold %.1f)\n', ...
        snr_centroid_sym, opt.M_pass, snr_centroid_acc, SNR_thr_Eve);
fprintf('  REQUIREMENT (pessimistic, closest-approach photon count): %.2e rad\n', eps_theta_req);
fprintf('  REQUIREMENT (accumulated epoch by epoch, the correct one): %.2e rad\n', ...
        opt.eps_theta_req);

%% ------------------------------------------------------------------------
% 4.11 Geometric low-probability-of-intercept argument
% -------------------------------------------------------------------------
% Under Kerckhoffs the remaining defence is geometric: the eavesdropper must
% be inside the beam. Half-angle theta_div/2 (full 1/e^2 divergence).
fprintf('\n--- 4.11 Geometric LPI ---\n');
a_half = theta_div/2;
fprintf('Divergence %.1f urad (full, 1/e^2) -> beam radius %.2f m, diameter %.2f m at %.0f km\n', ...
        theta_div*1e6, a_half*d_min, 2*a_half*d_min, d_min/1000);
fprintf('Illuminated fraction of solid angle = %.2e\n', a_half^2/4);

% --- [P12] Proximity inspector next to Bob (RPO) ----------------------------
% An inspector does not need to cross the beam mid-path: it can follow Bob and
% stay inside the beam near him. Alice points at Bob, so the beam is centred
% on Bob; an inspector at lateral offset rho stays inside while rho < w(t),
% w = (theta_div/2)*d(t). One placed a distance r behind Bob along the line of
% sight must follow the rotation of that line, at lateral speed r*omega(t).
dU = zeros(3,N_t);
for kk = 1:3, dU(kk,:) = gradient(uLOS(kk,:), sampleTime); end
omega_los = vecnorm(dU,2,1);                       % rad/s
w_beam    = a_half*d_t;                            % beam radius at Bob, m
fprintf('  [P12] beam radius at Bob: %.1f m at the window edge, %.2f m at closest approach\n', ...
        w_beam(1), w_beam(idx_cross));
fprintf('        line-of-sight rotation: %.2e rad/s at the edge, %.2e rad/s at closest approach,\n', ...
        omega_los(1), omega_los(idx_cross));
fprintf('        %.0f deg over the whole pass (the time Eve-C allows is printed in Section 8.3)\n', ...
        rad2deg(sum(omega_los)*sampleTime));

%% ------------------------------------------------------------------------
% 4.12 Consistency checks
% -------------------------------------------------------------------------
S_bit0 = sig.S_cl + sig.I_env.*cos(cfg.Phi);
S_bit1 = sig.S_cl - sig.I_env.*cos(cfg.Phi);
residual_avg = max(abs(0.5*(S_bit0+S_bit1) - sig.S_cl))/max(abs(sig.S_cl));
fprintf('\n--- 4.12 Consistency checks ---\n');
fprintf('Cancellation of the quantum term under bit averaging: relative residual = %.3e\n', residual_avg);

% Doppler peaks and envelope width against their analytic values.
[~, ipk] = findpeaks(sig.S_cl, 'MinPeakProminence', 0.1*max(sig.S_cl));
if numel(ipk) == 2
    peakPos_MHz = cfg.dw(ipk)/(2*pi)/1e6;
    peakPred_MHz = cfg.sep_over_Delta*cfg.Delta/2*f0_optical/1e6;
    fprintf('Doppler peaks measured at %+.3f / %+.3f MHz, predicted %+.3f MHz (error %.2f %%)\n', ...
            peakPos_MHz(1), peakPos_MHz(2), peakPred_MHz, ...
            100*abs(abs(peakPos_MHz(2))-peakPred_MHz)/peakPred_MHz);
else
    fprintf('Doppler peaks: %d found instead of 2 - the packets have merged.\n', numel(ipk));
end
wI     = sig.I_env/sum(sig.I_env);
sig_meas = sqrt(sum(wI.*(cfg.dw - sum(wI.*cfg.dw)).^2));
fprintf('Interference envelope width: measured %.3f MHz, analytic %.3f MHz\n', ...
        sig_meas/(2*pi)/1e6, sigma_env_rad/(2*pi)/1e6);
fprintf('Band occupancy: structure spans %.2f MHz of a %.2f MHz band\n', ...
        2*(cfg.sep_over_Delta/2+3)*cfg.Delta*f0_optical/1e6, cfg.B_IF/1e6);
fprintf('Emission weight of Eq.(26): departs from unity by %.1e over the band\n', ...
        3*cfg.B_IF/(2*f0_optical));

weight = sig.I_env.*cfg.key;                % envelope weighting and fringe key
k_expected = cfg.N/2 + 1;
Y_ref = fftshift(fft(ifftshift(sig.tmpl.*weight)));

% Read bin fixed analytically. tmpl.*weight = (1/2) I_env^2 (1 + exp(2i kappa dw)):
% a peak at zero frequency and an image of equal height at 2 kappa B_IF bins.
% Only the zero-frequency bin satisfies real(Y) = ||s||^2.
cfg.kpk = k_expected;
imageOffset = 2*cfg.kappa*cfg.B_IF;
peakCheck   = real(Y_ref(cfg.kpk))/sum(sig.tmpl.^2);
fprintf('Compression bin fixed at %d; image peak at %+.0f bins\n', ...
        cfg.kpk, imageOffset);
fprintf('Matched-filter identity real(Y(kpk))/||s||^2 = %.6f (must be 1)\n', peakCheck);
if abs(peakCheck-1) > 1e-3
    warning('Matched-filter identity violated: check kappa and the spectral grid.');
end

% Image clearance: the image must stay outside the compressed peak, whose
% half-width is N/(2 pi (sigma_env/sqrt(2))/dw_step) bins.
peakWidth_bins = cfg.N/(2*pi*((sigma_env_rad/sqrt(2))/dw_step));
fprintf('Peak/image separation: image at %.0f bins vs peak half-width %.0f bins, margin x%.1f\n', ...
        imageOffset, peakWidth_bins, imageOffset/peakWidth_bins);
if imageOffset < 3*peakWidth_bins
    Delta_x_clear = cfg.Delta_x*3*peakWidth_bins/imageOffset;
    warning(['Image peak within three half-widths of the read bin (margin %.1f). ' ...
             'Raising Delta_x to %.0f nm would clear it.'], ...
             imageOffset/peakWidth_bins, Delta_x_clear*1e9);
end

% Fringe spreading: signal collected by a signed-sum detector (no kappa)
% against the matched filter.
signedRatio = (sqrt(sum(sig.tmpl.^2))*sqrt(cfg.N))/abs(sum(sig.tmpl));
fprintf('Fringe spreading: signed-sum detector collects x%.1f less signal than the matched filter\n', ...
        signedRatio);

% Matched-filter gain over the unweighted DC-bin statistic.
g_mf   = sqrt(sum(sig.tmpl.^2));                          % ||s||
g_flat = sum(sig.tmpl.*real(cfg.key))/sqrt(cfg.N/2);      % unweighted statistic
fprintf('Matched-filter gain over the unweighted DC bin: x%.1f in deflection\n', ...
        g_mf/max(abs(g_flat),eps));
fprintf('Useful samples under the envelope: %.0f out of N = %d\n', ...
        sqrt(2*pi)*sigma_env_rad/dw_step, cfg.N);

% Deterministic part of the fast Monte Carlo path against the FFT bin.
A_sum  = sum(sig.tmpl.*real(weight));
A_fft  = real(Y_ref(cfg.kpk));
fprintf('Fast-path identity sum(tmpl.*Re w) / FFT bin = %.6f (must be 1)\n', A_sum/A_fft);
if abs(A_sum/A_fft - 1) > 1e-9
    warning('Sufficient-statistic shortcut disagrees with the transform.');
end

cfg.Nbit_frame = 1024;
cfg.N_MC       = 16;
cfg.block      = 32;

% --- [P13] Profile model and per-frame estimation basis ----------------------
% Bob's copy of S_cl may be wrong in height (gain), centre (Hz) and lobe width
% (relative, area kept, centres fixed). profileOf rebuilds the profile exactly
% for given errors. Bob fits the three quantities to the average of the frame
% on the basis X = [S_cl, dS/dcentre, dS/dwidth] (columns normalised).
sig.pp = struct('dw',cfg.dw,'Delta',cfg.Delta,'sep',cfg.sep_over_Delta, ...
                'theta',cfg.theta,'Om',Omega_rad,'kS',1);
raw0      = profileOf(sig.pp, 0, 0, 0);
sig.pp.kS = cfg.OD/sqrt(mean(raw0.^2));
profErr   = max(abs(profileOf(sig.pp,0,0,0) - sig.S_cl))/max(sig.S_cl);
if profErr > 1e-9
    warning('profileOf does not reproduce S_cl (relative error %.1e).', profErr);
end
hc = 1;  he = 1e-6;
Xb = [sig.S_cl, ...
      (profileOf(sig.pp,0, hc,0) - profileOf(sig.pp,0,-hc,0))/(2*hc), ...
      (profileOf(sig.pp,0,0, he) - profileOf(sig.pp,0,0,-he))/(2*he)];
sig.xs   = vecnorm(Xb,2,1);            % scale of each column (gain, per Hz, width)
sig.Xn   = Xb./sig.xs;
sig.XtXi = inv(sig.Xn.'*sig.Xn);
fprintf('\n--- 4.13 [P13] Per-frame estimation of the classical profile ---\n');
fprintf('profileOf reproduces S_cl to %.1e; estimation enabled: %d\n', profErr, cfg.profileFit);
% Bias that an UNCORRECTED error leaves on the statistic, as a fraction of the
% signal, (S_hat - S_cl).W / A, on several pitches (the statistic is shifted
% by minus this amount). After the fit only the part of the error outside
% the span of X remains (second order).
F_show = unique(round([opt.F_min, opt.F_now, 45, 60, 100, F_top]));
errSet = {'gain 1e-3', [1e-3 0 0]; 'centre 100 Hz', [0 100 0]; 'width 1e-3', [0 0 1e-3]; ...
          'gain 1e-2', [1e-2 0 0]; 'width 1e-2', [0 0 1e-2]};   % 1e-2: limit of the one-step linear fit
fprintf('  bias / signal for F = %s\n', mat2str(F_show));
for ie2 = 1:size(errSet,1)
    ev = errSet{ie2,2};
    dS = profileOf(sig.pp, ev(1), ev(2), ev(3)) - sig.S_cl;          % S_hat - S_cl
    dSf = dS - sig.Xn*(sig.XtXi*(sig.Xn.'*dS));                      % after a noiseless fit
    b0 = zeros(size(F_show)); b1 = b0;
    for jf = 1:numel(F_show)
        Wf = sig.I_env.*cos((F_show(jf)/cfg.B_IF)*cfg.dw);
        Af = sum(sig.I_env.*cos((F_show(jf)/cfg.B_IF)*cfg.dw).*Wf);
        b0(jf) = sum(dS.*Wf)/Af;  b1(jf) = sum(dSf.*Wf)/Af;
    end
    fprintf('  %-14s without fit: %s\n', errSet{ie2,1}, mat2str(b0,3));
    fprintf('  %-14s after fit  : %s\n', '', mat2str(b1,3));
end
% Precision of the per-frame estimates (noise of the frame average).
sdFit = @(sg) sg/sqrt(cfg.Nbit_frame)*sqrt(diag(sig.XtXi)).'./sig.xs;
sF = [sdFit(min(opt.sigma_ep)); sdFit(median(opt.sigma_ep)); sdFit(max(opt.sigma_ep))];
fprintf('  estimate precision (1 sd) at closest approach / median / edge:\n');
fprintf('    gain   %.1e / %.1e / %.1e\n', sF(:,1));
fprintf('    centre %.2f / %.2f / %.2f Hz\n', sF(:,2));
fprintf('    width  %.1e / %.1e / %.1e\n', sF(:,3));
cfg.smoothFixed = false;   % true: smooth errors applied with fixed size (tests of 7.1b)

%% ------------------------------------------------------------------------
% 4.14 [P7][P8][P11] Checks on the declared source
% -------------------------------------------------------------------------
% Deterministic orders of magnitude quoted in Section 3 of the thesis. No
% random numbers are drawn here, so the rest of the run is unchanged.
fprintf('\n--- 4.14 [P7][P8][P11] Checks on the declared source ---\n');
F_fine = [opt.F_min, opt.F_now, F_top];
% [P7] Photon recoil: known offset h/(2 m lambda^2), removed with the classical
% Doppler shift. Rounded constants of the code against exact ones.
m_exact   = (26.98153841 - 5.48579909e-4)*u_amu;   % 27Al+ ion: atom minus one electron
lam_exact = c_light/1122842857334736;           % 1S0-3P1 line
rec_code  = h_planck/(2*m_Al27*lambda_0^2);
rec_exact = h_planck/(2*m_exact*lam_exact^2);
ph_rec    = 2*pi*(F_fine/cfg.B_IF)*rec_exact;    % fringe rotation if not removed
fprintf('  [P7] recoil shift %.2f kHz (exact constants), %.2f kHz (code constants), difference %.0f Hz\n', ...
        rec_exact/1e3, rec_code/1e3, abs(rec_exact-rec_code));
fprintf('       uncompensated fringe rotation at F = %d / %.1f / %d: %.2f / %.2f / %.1f rad\n', ...
        F_fine(1), F_fine(2), F_fine(3), ph_rec);
fprintf('       error of the rounded constants on the finest fringe: %.3f rad\n', ...
        2*pi*(F_top/cfg.B_IF)*abs(rec_exact-rec_code));
fprintf('       recoil velocity hbar k/m = %.1f cm/s\n', 100*h_planck/(m_exact*lam_exact));
% [P8] Delta is the amplitude width; each lobe of |phi|^2 has sd Delta/sqrt(2).
v_lobe  = cfg.Delta*c_light/sqrt(2);
f_lobe  = 2*m_Al27*v_lobe^2/hbar/(2*pi);         % trap whose ground state has that spread
fprintf('  [P8] Delta = %.2f cm/s (amplitude width), lobe sd Delta/sqrt(2) = %.2f cm/s -> trap %.2f MHz\n', ...
        100*cfg.Delta*c_light, 100*v_lobe, f_lobe/1e6);
% [P11a] Number of emitters. Each emitter puts weight 1/2 in each lobe; the
% resonant cross-section 3 lambda^2/(2 pi) is spread over a Doppler lobe of
% sd f0*Delta/sqrt(2).
OD_pk    = max(sig.S_cl);
sig0     = 3*lambda_0^2/(2*pi);
sd_nu    = f0_optical*cfg.Delta/sqrt(2);
sigLobe  = 0.5*sig0*(pi*Gamma0_Hz/2)/(sqrt(2*pi)*sd_nu);
n_col    = OD_pk/sigLobe;
n_ph_sym = P_tx*cfg.T_sym/E_photon;
f_abs    = mean(1 - exp(-sig.S_cl));             % flat probe over the band
eps0     = 8.8541878128e-12;  B_pen = 7;
n_Bril   = eps0*B_pen^2/(2*m_Al27);               % Brillouin density limit
fprintf('  [P11a] peak OD %.2f -> column density %.2e m^-2; emitters for a %d / %d um probe radius: %.1e / %.1e\n', ...
        OD_pk, n_col, 10, 100, n_col*pi*(10e-6)^2, n_col*pi*(100e-6)^2);
fprintf('         probe photons per symbol %.2e, absorbed fraction %.3f -> at least %.1e emitters (<= 1 photon each)\n', ...
        n_ph_sym, f_abs, f_abs*n_ph_sym);
fprintf('         Brillouin limit at %d T: %.2e m^-3 -> length needed for that column %.0f m\n', ...
        B_pen, n_Bril, n_col/n_Bril);
% [P11b] Range of displacements of the alphabet and alignment tolerance for
% 90 %% of the fringe contrast, exp(-(dx*alpha)^2/(8 sigma_perp^2)) >= 0.9.
Dx_rng   = cfg.Delta_x*[opt.F_min, F_top]/opt.F_now;
sig_perp = sqrt(hbar/(2*m_Al27*2*pi*f_lobe));
alpha90  = sqrt(8*log(1/0.9))*sig_perp./Dx_rng;
fprintf('  [P11b] displacements %.0f nm - %.2f um; sigma_perp %.1f nm -> alignment %.2f rad / %.1f mrad\n', ...
        Dx_rng(1)*1e9, Dx_rng(2)*1e6, sig_perp*1e9, alpha90(1), 1e3*alpha90(2));
% [P11c] Time scales: symbol against the 3P1 lifetime, finest fringe period.
tau_3P1 = 1/(2*pi*Gamma0_Hz);
fprintf('  [P11c] T_sym %.0f us, lifetime %.0f us (ratio %.1f); finest fringe period %.1f kHz -> %.0f us to resolve\n', ...
        cfg.T_sym*1e6, tau_3P1*1e6, cfg.T_sym/tau_3P1, cfg.B_IF/F_top/1e3, 1e6*F_top/cfg.B_IF);
% Probe: one rectangular pulse per symbol; power at the band edge relative to
% the centre is sinc^2(B_IF/2 * tau_p).
tp = [60e-9 10e-9 5e-9];
edgeP = (sin(pi*cfg.B_IF/2*tp)./(pi*cfg.B_IF/2*tp)).^2;
fprintf('  [P6] probe pulse energy %.0f uJ, peak power %.1e W at 5 ns; band-edge power for 60/10/5 ns: %.2f / %.3f / %.3f\n', ...
        1e6*P_tx*cfg.T_sym, P_tx*cfg.T_sym/5e-9, edgeP);

% Monte Carlo path. The decision statistic is the zero-frequency bin of
% fft(Rx.*weight), i.e. sum(Rx.*weight):
%     stat = real(sum((S_cl - S_hat).*weight))          <- frame bias
%          + s*[A cos(dphi) - Aq sin(dphi)]             <- signal
%          + real(sum(noise.*weight)) ,                 <- noise
% and the noise term is a single Gaussian of variance sigma^2 sum(real(weight).^2).
% One variate per symbol replaces N, exactly. false = full spectral path.
cfg.fastMC = true;

% --- Fringe hopping inside the Monte Carlo -------------------------------
% [P14] Alice draws from all the admissible counts (option A, see 4.7b). For
% Bob the pitch matters only through the tracking loss, a smooth function of
% F, so the Monte Carlo samples that range with M_res representative pitches
% spread evenly over it. Campaign A keeps the declared pitch only.
if cfg.hopping
    cfg.F_alpha = unique(round(linspace(opt.F_min, F_top, max(1,round(opt.M_res)))));
else
    cfg.F_alpha = opt.F_now;
end
cfg_A = cfg;  cfg_A.F_alpha = opt.F_now;
% Template energy averaged over the alphabet, used for the adversaries'
% deflections in Section 8 (without Bob's tracking loss).
K_a   = cfg.F_alpha/cfg.B_IF;
g_eve = sqrt(mean(sum((sig.I_env.*cos(K_a.*cfg.dw)).^2, 1)));
fprintf('\nMonte Carlo alphabet: %d pitches, F from %d to %d (declared pitch F = %.1f used for Campaign A)\n', ...
        numel(cfg.F_alpha), min(cfg.F_alpha), max(cfg.F_alpha), opt.F_now);
trkLoss_a = exp(-(2*pi*cfg.sigma_f_track*K_a).^2/2);        % amplitude factor per pitch
fprintf('  tracking loss over the alphabet: mean %.3f dB, worst %.3f dB (amplitude); Bob''s d^2 falls by %.1f %% on average\n', ...
        -20*log10(mean(trkLoss_a)), -20*log10(min(trkLoss_a)), 100*(1-mean(trkLoss_a.^2)));
fprintf('  alphabet-averaged template energy vs declared pitch: %.3f\n', g_eve^2/sum(sig.tmpl.^2));

%% ========================================================================
% 5. CAMPAIGN A: BIT-ERROR-RATE WATERFALL AND ANALYTICAL VALIDATION
% =========================================================================
% The simulated BER must reproduce Q(||s||/sigma) at the declared pitch.
fprintf('\n=== 5. CAMPAIGN A: WATERFALL ===\n');
sigma_vec = logspace(log10(min(opt.sigma_ep)/6), log10(max(opt.sigma_ep)*6), 16);
nS = numel(sigma_vec); errA = zeros(1,nS); bitA = zeros(1,nS);
BER_theory_A = 0.5*erfc((g_mf./sigma_vec)/sqrt(2));
for is = 1:nS
    for mc = 1:cfg.N_MC
        rs = mcStream(SEED_MASTER,1,is,mc);
        [ne,nb] = runFrame(cfg_A, sig, sigma_vec(is), rs);
        errA(is) = errA(is)+ne; bitA(is) = bitA(is)+nb;
    end
    fprintf('  sigma = %.2e | d = %5.2f | BER simulated = %.3e | BER theory = %.3e\n', ...
            sigma_vec(is), g_mf/sigma_vec(is), errA(is)/bitA(is), BER_theory_A(is));
end
BER_A = errA./bitA; [BER_A_lo,BER_A_hi] = wilsonCI(errA,bitA,1.96);
% Comparison restricted to points with at least err_min errors; deviation
% also expressed in standard deviations of the binomial count.
cfg.err_min = 100;
valid = errA >= cfg.err_min;
if any(valid)
    relDev = abs(BER_A(valid)-BER_theory_A(valid))./BER_theory_A(valid);
    sdCnt  = sqrt(max(bitA(valid).*BER_theory_A(valid).*(1-BER_theory_A(valid)),eps));
    zDev   = abs(errA(valid)-bitA(valid).*BER_theory_A(valid))./sdCnt;
    fprintf('  maximum relative deviation over points with >= %d errors: %.1f %%\n', ...
            cfg.err_min, 100*max(relDev));
    fprintf('  maximum deviation in standard deviations of the count: %.1f sigma (%d points used)\n', ...
            max(zDev), sum(valid));
else
    fprintf('  no point reached %d errors: the waterfall grid is too coarse to validate.\n', cfg.err_min);
end
n_per_pt = cfg.Nbit_frame*cfg.N_MC;
fprintf('  campaign size: %d bits per point -> resolvable BER floor ~%.1e\n', ...
        n_per_pt, 1/n_per_pt);

%% ========================================================================
% 6. CAMPAIGN B: PERFORMANCE ALONG THE ORBITAL PASS
% =========================================================================
% Each realisation draws its own pointing loss, hence its own sigma.
% [P10] The adversaries need one noise level per epoch. Their evidence is
% D^2 = sum delta^2 with delta ~ d^2 ~ 1/sigma^2, so D^2 ~ sum 1/sigma^4 and
% the average is taken on 1/sigma^4:
%     sigma_B = (<sigma^-4>)^(-1/4)            (adversaries)
% The power-weighted mean (<sigma^-2>)^(-1/2) is kept as sigma_B2 for Bob's
% sensitivity runs and for comparison with the previous version.
fprintf('\n=== 6. CAMPAIGN B: ORBITAL PASS ===\n');
BER_B = zeros(1,N_epochs); errB = zeros(1,N_epochs); bitB = zeros(1,N_epochs);
sigma_B = zeros(1,N_epochs); sigma_B2 = zeros(1,N_epochs); sigma_B_arith = zeros(1,N_epochs);
for ie = 1:N_epochs
    ii = idx_ep(ie);
    acc = 0; accInv2 = 0; accInv4 = 0;
    for mc = 1:cfg.N_MC
        rs = mcStream(SEED_MASTER,2,ie,mc);
        th1 = sigma_jit_pat*randn(rs,1,2); th2 = sigma_jit_pat*randn(rs,1,2);
        Lp = exp(-8*(norm(th1)/theta_div)^2)*exp(-8*(norm(th2)/theta_div)^2);
        n_bin = (P_rx_ideal_t(ii)*Lp*cfg.T_sym/E_photon)*eta_det/cfg.N;
        sigma_eff = cfg.F_excess/sqrt(max(n_bin,eps));
        acc     = acc + sigma_eff;
        accInv2 = accInv2 + 1/sigma_eff^2;
        accInv4 = accInv4 + 1/sigma_eff^4;
        [ne,nb] = runFrame(cfg, sig, sigma_eff, rs);
        errB(ie) = errB(ie)+ne; bitB(ie) = bitB(ie)+nb;
    end
    sigma_B(ie)       = (accInv4/cfg.N_MC)^(-1/4);  % [P10] adversaries
    sigma_B2(ie)      = 1/sqrt(accInv2/cfg.N_MC);   % power-weighted (Bob, reference)
    sigma_B_arith(ie) = acc/cfg.N_MC;               % reference only
    BER_B(ie)   = errB(ie)/bitB(ie);
    fprintf('  t = %6.1f s | range = %7.1f km | sigma(Eve) = %.2e | sigma(power) = %.2e | d(Eve) = %5.2f | d(power) = %5.2f | BER = %.3e\n', ...
            t_ep(ie), d_t(ii)/1000, sigma_B(ie), sigma_B2(ie), g_eve/sigma_B(ie), g_eve/sigma_B2(ie), BER_B(ie));
end
% [P10] Effect of the averaging rule on the adversaries' accumulated D^2.
okB = avail_ep > 0;
fprintf('  [P10] noise for the adversaries: mean of sigma^-4. Against it, over the pass:\n');
fprintf('        power-weighted mean would understate D^2 by x%.3f (mean), x%.3f (min) to x%.3f (max)\n', ...
        mean((sigma_B2(okB)./sigma_B(okB)).^4), min((sigma_B2(okB)./sigma_B(okB)).^4), ...
        max((sigma_B2(okB)./sigma_B(okB)).^4));
fprintf('        arithmetic mean of sigma would understate D^2 by x%.3f (mean)\n', ...
        mean((sigma_B_arith(okB)./sigma_B(okB)).^4));
fprintf('        expected from the pointing statistics: x%.3f\n', eveGain^2);
[BER_B_lo,BER_B_hi] = wilsonCI(errB,bitB,1.96);
% Capacity is meaningful only where the BER is separated from one half by
% more than two Wilson half-widths.
sepFromHalf = (0.5 - BER_B)./max(BER_B_hi-BER_B, eps);
nWeak = sum(sepFromHalf < 2 & avail_ep > 0);
fprintf('  epochs whose capacity is within 2 sigma of zero: %d of %d\n', ...
        nWeak, sum(avail_ep > 0));
if nWeak > 0
    fprintf('  -> the capacity of those epochs is an upper bound, not a measurement.\n');
end

%% ========================================================================
% 7. SENSITIVITY TO CALIBRATION AND IMPRINT DEPTH
% =========================================================================
% Subtraction error on the classical profile, frozen over a frame (a bias).
% Reference: averaging M_cal symbols gives eps_sub ~ sigma/sqrt(M_cal).
fprintf('\n=== 7. SENSITIVITY ANALYSES ===\n');
fprintf('--- 7.1 Subtraction error on the classical profile ---\n');
eps_vec = [0 1e-4 1e-3 1e-2 3e-2 1e-1];
sigma_ref = median(sigma_B2);   % [P10] Bob's noise: power-weighted mean, as before
BER_eps = zeros(1,numel(eps_vec));
cfg_tmp = cfg;
cfg_tmp.eps_gain = 0; cfg_tmp.eps_c_Hz = 0; cfg_tmp.eps_w = 0;   % [P13] per-bin error only
for k = 1:numel(eps_vec)
    cfg_tmp.eps_sub = eps_vec(k);
    e = 0; b = 0;
    for mc = 1:cfg.N_MC
        rs = mcStream(SEED_MASTER,3,k,mc);
        [ne,nb] = runFrame(cfg_tmp, sig, sigma_ref, rs);
        e = e+ne; b = b+nb;
    end
    BER_eps(k) = e/b;
    fprintf('  eps_sub = %.1e (equivalent to averaging %.1e calibration symbols) -> BER = %.3e\n', ...
            eps_vec(k), (sigma_ref/max(eps_vec(k),eps))^2, BER_eps(k));
end
fprintf('  This sweep models an independent per-bin error, which averages over the\n');
fprintf('  useful samples under the envelope. Smooth errors are tested in 7.1b.\n');

% --- 7.1b [P13] Smooth calibration errors, before and after the estimation ---
% Fixed-size error on one quantity at a time (gain, centre, width), per-bin
% error off, at the median noise and at the closest-approach noise, on the
% reference pitch alone and on the hopping alphabet; receiver without and with
% the per-frame estimation of Section 4.13.
fprintf('--- 7.1b [P13] Smooth calibration errors (BER; fixed error size) ---\n');
sigTest   = [sigma_ref, min(sigma_B2)];
sigName   = {'median', 'closest'};
errType   = {'gain', 'centre [Hz]', 'width'};
errVals   = {[1e-4 1e-3 1e-2], [1e1 1e2 1e3], [1e-4 1e-3 1e-2]};
alphaSet  = {cfg_A, cfg};  alphaName = {'reference pitch', 'alphabet'};
maxDev    = zeros(2,2);   % [no fit, fit] x [median, closest]: max |BER - BER(no error)|
cnt = 0;
for ia = 1:2
    for isg = 1:2
        cb = alphaSet{ia};
        cb.eps_sub = 0; cb.eps_gain = 0; cb.eps_c_Hz = 0; cb.eps_w = 0; cb.smoothFixed = true;
        % reference: no error
        cnt = cnt + 1;  e = 0; b = 0;
        for mc = 1:cfg.N_MC
            rs = mcStream(SEED_MASTER,4,cnt,mc);
            [ne,nb] = runFrame(cb, sig, sigTest(isg), rs);  e = e+ne; b = b+nb;
        end
        BER0 = e/b;
        fprintf('  %s, %s noise (sigma %.2e): no error BER = %.4f\n', ...
                alphaName{ia}, sigName{isg}, sigTest(isg), BER0);
        for it = 1:3
            for iv = 1:numel(errVals{it})
                BERfit = zeros(1,2);
                for fitOn = 0:1
                    ct = cb;  ct.profileFit = logical(fitOn);
                    switch it
                        case 1, ct.eps_gain = errVals{it}(iv);
                        case 2, ct.eps_c_Hz = errVals{it}(iv);
                        case 3, ct.eps_w    = errVals{it}(iv);
                    end
                    cnt = cnt + 1;  e = 0; b = 0;
                    for mc = 1:cfg.N_MC
                        rs = mcStream(SEED_MASTER,4,cnt,mc);
                        [ne,nb] = runFrame(ct, sig, sigTest(isg), rs);  e = e+ne; b = b+nb;
                    end
                    BERfit(fitOn+1) = e/b;
                    maxDev(fitOn+1,isg) = max(maxDev(fitOn+1,isg), abs(e/b - BER0));
                end
                fprintf('    %-12s %8.1e : BER without estimation %.4f | with estimation %.4f\n', ...
                        errType{it}, errVals{it}(iv), BERfit(1), BERfit(2));
            end
        end
    end
end
sdBER = sqrt(0.25/(cfg.Nbit_frame*cfg.N_MC));
fprintf('  largest |BER - BER(no error)|: without estimation %.4f (median) / %.4f (closest);\n', ...
        maxDev(1,1), maxDev(1,2));
fprintf('                                 with estimation    %.4f (median) / %.4f (closest);\n', ...
        maxDev(2,1), maxDev(2,2));
fprintf('  counting noise of one BER point about %.4f (difference of two points about %.4f)\n', ...
        sdBER, sqrt(2)*sdBER);

% The deflection is linear in OD: no Monte Carlo needed. OD clipped at 1.
fprintf('--- 7.2 Imprint depth (deflection at the pass-median noise, faded) ---\n');
OD_vec = unique(min(cfg.OD*[0.1 0.25 0.5 1 1.25], 1));
for k = 1:numel(OD_vec)
    d_k = (OD_vec(k)/cfg.OD)*(g_mf/sigma_ref);
    flag = '';
    if abs(OD_vec(k)-cfg.OD) < 1e-12,  flag = '  <- design point';  end
    fprintf('  OD = %8.2e -> d = %6.3f -> raw BER = %.3e -> capacity = %.4f bit/symbol%s\n', ...
            OD_vec(k), d_k, 0.5*erfc(d_k/sqrt(2)), ...
            1-binaryEntropy(0.5*erfc(d_k/sqrt(2))), flag);
end

%% ------------------------------------------------------------------------
% 7.3 Admissible interval of the branch separation
% -------------------------------------------------------------------------
% Delta_x sets readability of the fringe, not its strength. Bounds:
%   from below - image clearance (Section 4.12);
%   from above - Nyquist sampling (Ceiling 1) and natural linewidth (Ceiling 2).
% Guard factors are declared: three half-widths for the image, ten on the
% linewidth.
fprintf('--- 7.3 Branch separation: admissible interval ---\n');
g_image = 3;     % required image-to-read-bin separation, in peak half-widths
g_gamma = 10;    % required ratio of fringe spacing to natural linewidth
kap_of  = @(dx) m_Al27*c_light*dx/(hbar*Omega_rad);
Dx_lo_img = g_image*peakWidth_bins/(2*cfg.B_IF) * (hbar*Omega_rad)/(m_Al27*c_light);
Dx_hi_nyq = (pi/dw_step)                        * (hbar*Omega_rad)/(m_Al27*c_light);
Dx_hi_gam = 1/(g_gamma*Gamma0_Hz)               * (hbar*Omega_rad)/(m_Al27*c_light);
Dx_hi     = min(Dx_hi_nyq, Dx_hi_gam);
fprintf('  lower bound (image clearance, x%d)      : %8.1f nm\n', g_image, Dx_lo_img*1e9);
fprintf('  upper bound (Nyquist on the grid)       : %8.1f nm\n', Dx_hi_nyq*1e9);
fprintf('  upper bound (natural linewidth, x%d)    : %8.1f nm\n', g_gamma, Dx_hi_gam*1e9);
fprintf('  admissible interval                     : [%.1f, %.1f] nm, declared value %.1f nm\n', ...
        Dx_lo_img*1e9, Dx_hi*1e9, cfg.Delta_x*1e9);
if cfg.Delta_x < Dx_lo_img || cfg.Delta_x > Dx_hi
    warning('Declared Delta_x lies outside the admissible interval computed in Section 7.3.');
end
fprintf('  margin to the binding bound             : x%.2f below, x%.2e above\n', ...
        cfg.Delta_x/Dx_lo_img, Dx_hi/cfg.Delta_x);
% Numerical sweep: template energy and deflection at multiples of Delta_x.
Dx_vec = cfg.Delta_x*[0.25 0.5 0.75 1 2 5 20];
fprintf('  %10s %10s %10s %10s %10s %10s %10s\n', ...
        'Dx [nm]','fringes','img/width','Nyq [rad]','sp/Gamma0','E/E_ref','d/d_ref');
d_ref_dx = g_mf/sigma_ref;
for k = 1:numel(Dx_vec)
    kap_k  = kap_of(Dx_vec(k));
    tmpl_k = sig.I_env.*cos(kap_k*cfg.dw);
    d_k    = sqrt(sum(tmpl_k.^2))/sigma_ref;
    fprintf('  %10.1f %10.1f %10.2f %10.4f %10.0f %10.4f %10.4f\n', ...
            Dx_vec(k)*1e9, cfg.B_IF*kap_k, (2*kap_k*cfg.B_IF)/peakWidth_bins, ...
            kap_k*dw_step, 1/(kap_k*Gamma0_Hz), (d_k/d_ref_dx)^2, d_k/d_ref_dx);
end
fprintf('  the deflection is flat at and above the declared value: Delta_x sets readability,\n');
fprintf('  not signal strength. Below it the column rises, because with few fringes cos^2 no\n');
fprintf('  longer averages to one half - but that gain is unusable, since it is exactly where\n');
fprintf('  the image peak walks into the read bin.\n');

%% ========================================================================
% 8. CAPACITY, COVERT VOLUME AND SUMMARY
% =========================================================================
% Per-symbol capacity 1 - H2(BER), accumulated over the pass until each
% adversary's accumulated deflection reaches threshold.
fprintf('\n=== 8. CAPACITY AND COVERT VOLUME ===\n');
n_sym_ep   = (dt_ep.*avail_ep)/cfg.T_sym;
C_ep_sim   = 1 - binaryEntropy(BER_B);
d_ep_sim   = g_eve./sigma_B;   % alphabet-averaged template energy,
                               % [P10] noise averaged on sigma^-4
delta_A_ep = d_ep_sim.^2/sqrt(2*opt.N_env);
delta_B_ep = d_ep_sim.^2/sqrt(2);
delta_C_ep = d_ep_sim.^2/sqrt(2*opt.M_eff);

bits_cum = cumsum(n_sym_ep.*C_ep_sim);
eveA_cum = sqrt(cumsum(n_sym_ep.*delta_A_ep.^2));
eveB_cum = sqrt(cumsum(n_sym_ep.*delta_B_ep.^2));
eveC_cum = sqrt(cumsum(n_sym_ep.*delta_C_ep.^2));

B_covert_A = volumeAtThreshold(bits_cum, eveA_cum, SNR_thr_Eve, n_sym_ep, C_ep_sim, delta_A_ep);
B_covert_B = volumeAtThreshold(bits_cum, eveB_cum, SNR_thr_Eve, n_sym_ep, C_ep_sim, delta_B_ep);
B_covert_C = volumeAtThreshold(bits_cum, eveC_cum, SNR_thr_Eve, n_sym_ep, C_ep_sim, delta_C_ep);

% Time and epoch at which each adversary reaches threshold.
tStop = @(ev) stopEpochTime(ev, SNR_thr_Eve, t_ep);
[tsA, isA] = tStop(eveA_cum); [tsC, isC] = tStop(eveC_cum); [tsB, isB] = tStop(eveB_cum);

% --- Second threshold ----------------------------------------------------
% thr gives the eavesdropper's optimal test a total error P_FA + P_MD =
% 2 Q(thr/2) at stop (0.13 for thr = 3). Volumes also at thr = 1, same
% operating point (not re-optimised).
thr_alt = 1;
B_covert_A1 = volumeAtThreshold(bits_cum, eveA_cum, thr_alt, n_sym_ep, C_ep_sim, delta_A_ep);
B_covert_C1 = volumeAtThreshold(bits_cum, eveC_cum, thr_alt, n_sym_ep, C_ep_sim, delta_C_ep);
B_covert_B1 = volumeAtThreshold(bits_cum, eveB_cum, thr_alt, n_sym_ep, C_ep_sim, delta_B_ep);
Pe_tot = @(t) erfc(t/(2*sqrt(2)));          % 2 Q(t/2)

% --- Eve-C: conservative lower bound -------------------------------------
% Epochs whose unfaded design deflection exceeds half the identification
% threshold sqrt(2 ln M) are counted with Eve-B's deflection, as if the
% pitch had been identified there. The Eve-C volume lies in [lower, nominal].
d_id_c   = sqrt(2*log(max(opt.M_eff,2)));
d_uf_ep  = opt.d_ep/sqrt(L_point_mean);                 % unfaded design deflection
idMask   = d_uf_ep > 0.5*d_id_c;
delta_Cc_ep = delta_C_ep;  delta_Cc_ep(idMask) = delta_B_ep(idMask);
eveCc_cum   = sqrt(cumsum(n_sym_ep.*delta_Cc_ep.^2));
B_covert_Cc = volumeAtThreshold(bits_cum, eveCc_cum, SNR_thr_Eve, n_sym_ep, C_ep_sim, delta_Cc_ep);
B_covert_Cc1 = volumeAtThreshold(bits_cum, eveCc_cum, 1, n_sym_ep, C_ep_sim, delta_Cc_ep);

% --- Full-band radiometer (reference only, not a threat model) ------------
delta_Af_ep = d_ep_sim.^2/sqrt(2*cfg.N);
eveAf_cum   = sqrt(cumsum(n_sym_ep.*delta_Af_ep.^2));
B_covert_Af = volumeAtThreshold(bits_cum, eveAf_cum, SNR_thr_Eve, n_sym_ep, C_ep_sim, delta_Af_ep);

% --- Key cost of the hopping sequence ------------------------------------
% [P14] Alice draws from all M_alpha admissible counts (option A), so each
% symbol consumes log2(M_alpha) bits of keystream. [P15-iii] Compared with the
% O(sqrt(n) log n) of Bash et al. (reference of the thesis) and with the
% tighter O(sqrt(n)) of Bloch.
key_per_sym = log2(max(opt.M_alpha,2));
key_total   = key_per_sym*opt.M_pass;
key_srl     = sqrt(opt.M_pass)*log2(max(opt.M_pass,2));   % Bash
key_bloch   = sqrt(opt.M_pass);                           % Bloch
% Uniform index in [1, M_alpha] from the keystream by rejection sampling on
% ceil(log2 M_alpha)-bit words (no modulo bias): expected bits per symbol.
nb_word     = ceil(log2(max(opt.M_alpha,2)));
key_drawn   = nb_word/(opt.M_alpha/2^nb_word);
[~,ie_cr] = min(abs(idx_ep-idx_cross));

fprintf('\n============ SUMMARY: QUANTUM DOPPLER COVERT CHANNEL ============\n');
fprintf('Branch separation Delta_x                : %8.3f um\n', cfg.Delta_x*1e6);
fprintf('Fringe spacing / fringes in band         : %8.3f MHz / %.1f\n', ...
        1/cfg.kappa/1e6, cfg.B_IF*cfg.kappa);
fprintf('Optimised packet separation              : %8.2f Delta\n', cfg.sep_over_Delta);
fprintf('Imprint depth OD (declared)              : %8.3f\n', cfg.OD);
fprintf('Theta calibration required of the source : %8.2e rad\n', opt.eps_theta_req);
fprintf('Hopping alphabet (integer fringe counts) : %8d\n', opt.M_alpha);
fprintf('Resolvable pitches M (cells for Eve-C)   : %8.1f\n', opt.M_eff);
fprintf('Derived intermediate-frequency bandwidth : %8.2f MHz\n', cfg.B_IF/1e6);
fprintf('Raw symbol rate                          : %8.3f kbit/s\n', R_b_raw/1e3);
fprintf('Matched-filter gain over unweighted sum  : x%7.1f in deflection\n', g_mf/max(abs(g_flat),eps));
fprintf('Raw BER at closest approach              : %8.3f\n', BER_B(ie_cr));
fprintf('Capacity at closest approach             : %8.4f bit/symbol\n', C_ep_sim(ie_cr));
fprintf('Median capacity over the pass            : %8.4f bit/symbol\n', median(C_ep_sim(avail_ep>0)));
fprintf('Peak information throughput              : %8.3f kbit/s\n', max(R_b_raw*C_ep_sim.*avail_ep)/1e3);
fprintf('----------------------------------------------------------------\n');
fprintf('VOLUME, no eavesdropper                  : %8.3e bit\n', bits_cum(end));
fprintf('COVERT VOLUME vs Eve-A (envelope energy)  : %8.3e bit  (N_env = %.0f)\n', B_covert_A, opt.N_env);
fprintf('  Eve-A statistic at end of pass         : %8.2f  (threshold %.1f)\n', eveA_cum(end), SNR_thr_Eve);
fprintf('COVERT VOLUME vs Eve-C (hop-blind)       : %8.3e bit  (M = %.1f)\n', B_covert_C, opt.M_eff);
fprintf('  Eve-C statistic at end of pass         : %8.2f  (threshold %.1f)\n', eveC_cum(end), SNR_thr_Eve);
fprintf('  Eve-C lower bound (identification assumed where unfaded d > %.2f, %d of %d epochs): %8.3e bit\n', ...
        0.5*d_id_c, sum(idMask), numel(idMask), B_covert_Cc);
fprintf('  -> defensible covert volume vs Eve-C lies in [%.3e, %.3e] bit\n', B_covert_Cc, B_covert_C);
fprintf('COVERT VOLUME vs Eve-B (key compromised) : %8.3e bit\n', B_covert_B);
fprintf('  Eve-B statistic at end of pass         : %8.2e (threshold %.1f)\n', eveB_cum(end), SNR_thr_Eve);
fprintf('  reference: full-band radiometer (N = %d, not a threat model): %8.3e bit\n', ...
        cfg.N, B_covert_Af);
fprintf('THRESHOLD CROSSING (thr = %.1f), epoch where it happens:\n', SNR_thr_Eve);
printStop('Eve-A', tsA, isA, idx_ep, d_t);
printStop('Eve-C', tsC, isC, idx_ep, d_t);
printStop('Eve-B', tsB, isB, idx_ep, d_t);
fprintf('----------------------------------------------------------------\n');
fprintf('THRESHOLD READING: thr = %.1f -> eavesdropper total error P_FA+P_MD = %.3f at stop\n', ...
        SNR_thr_Eve, Pe_tot(SNR_thr_Eve));
fprintf('  same operating point, thr = %.1f (P_FA+P_MD = %.3f):\n', thr_alt, Pe_tot(thr_alt));
fprintf('    vs Eve-A %8.3e bit | vs Eve-C %8.3e bit (lower bound %8.3e) | vs Eve-B %8.3e bit\n', ...
        B_covert_A1, B_covert_C1, B_covert_Cc1, B_covert_B1);
% [P1][P2] Full threshold grid at the same operating point: D = 3 (design) down
% to D = 0.25, the criterion of the laboratory experiments (P_FA+P_MD >= 0.9).
% This replaces the former prose estimate of "a few hundred bits" at D ~ 0.25.
fprintf('  [P1][P2] volume vs each threshold (same operating point, not re-optimised):\n');
fprintf('    %6s  %10s  %12s  %12s  %12s\n', 'D', 'P_FA+P_MD', 'Eve-A', 'Eve-C', 'Eve-B');
for tg = thr_grid
    fprintf('    %6.2f  %10.3f  %12.3e  %12.3e  %12.3e\n', tg, Pe_tot(tg), ...
            volumeAtThreshold(bits_cum, eveA_cum, tg, n_sym_ep, C_ep_sim, delta_A_ep), ...
            volumeAtThreshold(bits_cum, eveC_cum, tg, n_sym_ep, C_ep_sim, delta_C_ep), ...
            volumeAtThreshold(bits_cum, eveB_cum, tg, n_sym_ep, C_ep_sim, delta_B_ep));
end
% [P9] Eve's receiver carries the same 3 dB heterodyne penalty as Bob's. An Eve
% without it (homodyne or better) halves her noise power: her d^2, hence each
% delta, doubles. Volume against such an Eve, same operating point.
volNoHet = @(dl, thr) volumeAtThreshold(bits_cum, sqrt(cumsum(n_sym_ep.*(2*dl).^2)), ...
                                        thr, n_sym_ep, C_ep_sim, 2*dl);
B_C_nh3 = volNoHet(delta_C_ep, SNR_thr_Eve);  B_C_nh1 = volNoHet(delta_C_ep, 1);
fprintf('  [P9] Eve without the 3 dB heterodyne penalty (her delta x2):\n');
fprintf('       vs Eve-C %8.3e bit at D = %.0f (x%.2f lower), %8.3e bit at D = 1 (x%.2f lower)\n', ...
        B_C_nh3, SNR_thr_Eve, B_covert_C/max(B_C_nh3,eps), B_C_nh1, B_covert_C1/max(B_C_nh1,eps));
fprintf('       vs Eve-A %8.3e bit (x%.2f lower), vs Eve-B %8.3e bit (x%.2f lower), at D = %.0f\n', ...
        volNoHet(delta_A_ep, SNR_thr_Eve), B_covert_A/max(volNoHet(delta_A_ep, SNR_thr_Eve),eps), ...
        volNoHet(delta_B_ep, SNR_thr_Eve), B_covert_B/max(volNoHet(delta_B_ep, SNR_thr_Eve),eps), SNR_thr_Eve);
fprintf('  Eve-C / Eve-A volume ratio             : %8.3f (asymptotic sqrt(M/N_env) = %.3f)\n', ...
        B_covert_C/max(B_covert_A,eps), sqrt(opt.M_eff/opt.N_env));
fprintf('    The asymptotic ratio assumes equal epochs and truncation deep inside\n');
fprintf('    the pass for both adversaries. Here the epochs are strongly unequal and\n');
fprintf('    the two adversaries truncate at different points of the pass, so the\n');
fprintf('    measured ratio differs from sqrt(M/N_env). The square-root scaling itself\n');
fprintf('    is verified in Sections 8.2 and 9.1.\n');
fprintf('----------------------------------------------------------------\n');
fprintf('KEY COST: %.2f bit/symbol (%d admissible counts) -> %.2e bit per pass\n', ...
        key_per_sym, opt.M_alpha, key_total);
fprintf('  [P15-iii] Bash et al., sqrt(n) log2 n : %.2e bit -> ratio %.0fx\n', key_srl, key_total/max(key_srl,eps));
fprintf('            Bloch, sqrt(n)              : %.2e bit -> ratio %.0fx\n', key_bloch, key_total/max(key_bloch,eps));
fprintf('  [P14] keystream drawn with rejection sampling on %d-bit words: %.1f bit/symbol on average\n', ...
        nb_word, key_drawn);
fprintf('  the channel is NOT key efficient. A stream cipher seeded\n');
fprintf('  with a short pre-shared key closes the gap, at the price of\n');
fprintf('  computational rather than information-theoretic covertness.\n');
fprintf('----------------------------------------------------------------\n');
if cfg.hopping
    fprintf('With hopping enabled the defensible figure is the Eve-C volume, and it\n');
    fprintf('holds under the Kerckhoffs assumption because the secret is now a key\n');
    fprintf('and not a constant of the apparatus. Eve-A remains an upper bound and\n');
    fprintf('Eve-B the key-compromise case, where only the geometric argument of\n');
    fprintf('Section 4.11 survives. Hopping does not evade the square root law.\n');
else
    fprintf('The defensible figure is the Eve-A volume, and it holds only if covertness\n');
    fprintf('rests on the geometric argument of Section 4.11. Under the Kerckhoffs\n');
    fprintf('assumption (Eve-B) the channel is not covert.\n');
end
fprintf('================================================================\n');

%% ------------------------------------------------------------------------
% 8.3 [P12] Proximity inspector: time it must stay inside the beam
% -------------------------------------------------------------------------
% An inspector that follows Bob stays inside the beam while its lateral offset
% is below the beam radius w = (theta_div/2)*d. What the geometry costs it is
% not a launch "knowing where the link is" but keeping station on Bob while
% the line of sight rotates: over the time Alice actually transmits against
% Eve-C, the line turns by the angle below, i.e. the inspector sweeps that arc
% at radius equal to its stand-off distance behind Bob.
if ~isnan(tsC)
    mC = (t_sec >= 0) & (t_sec <= tsC);
    rotToStop = rad2deg(sum(omega_los(mC))*sampleTime);
    fprintf('\n--- 8.3 [P12] Proximity inspector (RPO) ---\n');
    fprintf('  Eve-C forces Alice to stop at t ~ %.0f s; until then the line of sight turns %.1f deg\n', ...
            tsC, rotToStop);
    fprintf('  (beam radius at Bob %.2f m at closest approach, %.1f m at the edge): the inspector must\n', ...
            w_beam(idx_cross), w_beam(1));
    fprintf('  hold station on Bob within the beam radius, not cross the beam at mid-path.\n');
    [~, iStopT] = min(abs(t_sec - tsC));
    fprintf('  while Alice transmits the beam radius at Bob falls from %.1f m to %.1f m (range %.0f -> %.0f km)\n', ...
            w_beam(1), w_beam(iStopT), d_t(1)/1e3, d_t(iStopT)/1e3);
end

%% ------------------------------------------------------------------------
% 8.4 [P12] Keystream revealed by partial identification of the pitch
% -------------------------------------------------------------------------
% Eve-C observes the M resolvable matched-filter outputs of each symbol. What
% she learns about the pitch is the mutual information between the pitch
% index and those outputs (sign of the bit unknown):
%     I(d) = log2 M - E[H(posterior)],  posterior_j ~ cosh(d y_j).
% It is the keystream information that leaks, summed over the symbols sent
% before the stop. With AES-CTR the leaked keystream says nothing about the
% key or about the other blocks, so the damage is limited to those symbols.
% Evaluated with the adversaries' fading-averaged deflection and, as an upper
% bound, with the unfaded design deflection (Eve's best case).
if cfg.hopping && ~isnan(tsC)
    Mi   = max(2, round(opt.M_eff));
    rsMI = RandStream('mrg32k3a','Seed',SEED_MASTER); rsMI.Substream = 11;
    d_uf_eve = opt.d_ep/sqrt(L_point_mean);           % unfaded, Eve's best case
    dGrid = linspace(0, 1.05*max([d_ep_sim(:); d_uf_eve(:)]), 31);
    Igrid = pitchInfo(dGrid, Mi, 4000, rsMI);
    I_avg = interp1(dGrid, Igrid, d_ep_sim);
    I_uf  = interp1(dGrid, Igrid, d_uf_eve);
    % symbols sent before Eve-C stops the transmission (same rule as the volume)
    n_tx = zeros(size(n_sym_ep));
    n_tx(1:isC-1) = n_sym_ep(1:isC-1);
    if isC == 1, e0 = 0; else, e0 = eveC_cum(isC-1)^2; end
    n_tx(isC) = n_sym_ep(isC)*min(max((SNR_thr_Eve^2 - e0)/max(n_sym_ep(isC)*delta_C_ep(isC)^2,eps),0),1);
    leak_avg = sum(n_tx.*I_avg);  leak_uf = sum(n_tx.*I_uf);
    fprintf('\n--- 8.4 [P12] Keystream revealed before the stop (M = %d) ---\n', Mi);
    fprintf('  symbols sent %.2e, keystream used for the pitch %.2e bit (%.2f bit/symbol)\n', ...
            sum(n_tx), sum(n_tx)*key_per_sym, key_per_sym);
    fprintf('  information on the pitch revealed to Eve-C: %.1f bit (fading-averaged), %.1f bit (unfaded bound)\n', ...
            leak_avg, leak_uf);
    fprintf('  per symbol %.1e / %.1e bit; largest single-epoch value %.3f bit/symbol (unfaded)\n', ...
            leak_avg/max(sum(n_tx),eps), leak_uf/max(sum(n_tx),eps), max(I_uf(n_tx>0)));
    fprintf('  for reference: a fully identified pitch would reveal %.2f bit (log2 M)\n', log2(Mi));
end

%% ------------------------------------------------------------------------
% 8.1 Control case: stationary geometry
% -------------------------------------------------------------------------
% Same chain re-optimised with the received power frozen at its closest-
% approach value, full availability and the same total dwell: an upper
% bound on a fixed link, used to measure what the pass geometry costs.
fprintf('\n--- 8.1 Control: stationary geometry at closest-approach range ---\n');
par_stat = par;
par_stat.P_ep  = P_rx_ideal_t(idx_cross)*L_point_mean*ones(1,N_epochs);   % same fading as the pass design
par_stat.dt_ep = (T_window/N_epochs)*ones(1,N_epochs);
B_stat = zeros(1,nSep);
for k = 1:nSep
    o_s = operatingPoint(sep_vec(k), par_stat);
    if o_s.nyquistPhase < pi
        % Same target adversary as the pass optimum.
        if cfg.hopping, B_stat(k) = o_s.volume_hopped; else, B_stat(k) = o_s.volume_naive; end
    end
end
[Bs_best, ks_best] = max(B_stat);
o_stat = operatingPoint(sep_vec(ks_best), par_stat);
fprintf('  optimal separation %.2f Delta (pass: %.2f) -> covert volume %.3e bit vs %.3e for the pass\n', ...
        sep_vec(ks_best), cfg.sep_over_Delta, Bs_best, B_best);
fprintf('  ratio stationary/pass = %.2f ; median deflection %.3f vs %.3f (design, mean fading)\n', ...
        Bs_best/B_best, median(o_stat.d_ep), median(opt.d_ep));
fprintf('  no eavesdropper: %.3e bit vs %.3e ; key-aware Eve-B: %.3e bit vs %.3e\n', ...
        o_stat.volume_noEve, opt.volume_noEve, o_stat.volume_keyed, opt.volume_keyed);
if Bs_best/B_best > 1.5
    fprintf('  -> the pass is the conservative case: a fixed link at best range would carry x%.1f more.\n', ...
            Bs_best/B_best);
else
    fprintf('  -> the variability of the pass costs little: the result is close to that of a fixed link.\n');
end

%% ------------------------------------------------------------------------
% 8.2 Scaling of the covert volume with the length of the window
% -------------------------------------------------------------------------
% Measures the square-root-law exponent. Packet separation frozen, OD
% scanned as a proxy for the deflection, per-epoch link conditions kept and
% only the dwell stretched. Volume against the target adversary.
fprintf('\n--- 8.2 Scaling of the covert volume with window length ---\n');
scale_vec = [1 4 16 64 256];
g_vec     = logspace(-4, log10(1/cfg.OD), 60);   % OD multiplier, OD clipped at 1
B_scale = zeros(1,numel(scale_vec));  d_scale = zeros(1,numel(scale_vec));
for is = 1:numel(scale_vec)
    par_s = par;  par_s.dt_ep = par.dt_ep*scale_vec(is);
    best = 0;  bestd = 0;
    for ig = 1:numel(g_vec)
        par_s.OD = cfg.OD*g_vec(ig);
        o_s = operatingPoint(cfg.sep_over_Delta, par_s);
        if cfg.hopping, v_s = o_s.volume_hopped; else, v_s = o_s.volume_naive; end
        if v_s > best
            best = v_s;  bestd = median(o_s.d_ep);
        end
    end
    B_scale(is) = best;  d_scale(is) = bestd;
    fprintf('  window x%-4d (%7.0f s) -> covert volume %.3e bit at median deflection %.3f\n', ...
            scale_vec(is), scale_vec(is)*T_window, best, bestd);
end
% Local exponents: the law is asymptotic in the small-deflection limit.
slopes = diff(log(B_scale))./diff(log(scale_vec));
fprintf('  local exponents between consecutive windows: ');  fprintf('%.3f ', slopes);  fprintf('\n');
fprintf('  at the operating point (first interval)  : %.3f\n', slopes(1));
fprintf('  in the asymptotic regime (last interval) : %.3f  (square-root law predicts 0.500)\n', ...
        slopes(end));
if abs(slopes(end)-0.5) > 0.05
    warning(['Asymptotic covert-volume exponent is %.3f rather than one half: either the ' ...
             'window sweep does not reach the small-deflection limit, or the truncation ' ...
             'model is wrong.'], slopes(end));
end
if abs(slopes(1)-slopes(end)) < 0.02
    fprintf('  note: the operating point already sits in the asymptotic regime.\n');
else
    fprintf('  the operating point departs from the asymptotic law by %+.3f in exponent.\n', ...
            slopes(1)-slopes(end));
end

%% ========================================================================
% 9. FIGURES
% =========================================================================
% Draws the seven figures of Chapter 6 directly in their final form (the
% former fixFigures.m is merged here). Nothing is saved: export each figure
% by hand with the file name printed at the end of this block.
%   Fig 1ab -> Sec. 6.2 (signal and receiver output)
%   Fig 1c  -> Sec. 6.3 (Campaign A)
%   Fig 2a  -> Sec. 6.4 (Campaign B)
%   Fig 2b  -> Sec. 6.5 (cumulative volume, stop of each adversary)
%   Fig 3a  -> Sec. 6.1 (choice of the packet separation)
%   Fig 3b  -> Sec. 6.9 (sensitivity to the optical depth)
%   Fig 4   -> Sec. 6.6 (hopping alphabet)
f_MHz = cfg.dw/(2*pi)/1e6;

% Time convention for Fig 2b.
%   'printed'  : epoch k spans [t_ep(k), t_ep(k)+45 s), the convention of
%                stopEpochTime, so the stop times match the printout and the
%                text. The last epoch is drawn up to 945 s.
%   'physical' : epoch k spans its real weight of 900/21 = 42.9 s; the stop
%                times shift by up to ~20 s.
timeConv = 'printed';

% -------------------------------------------------------------------------
% FIG 1ab: SPECTRAL FRINGE AND MATCHED FILTERING
% -------------------------------------------------------------------------
figure('Name','Fig 1ab Spectral fringe and matched filtering','Color','w','Position',[80 100 900 380]);

% (a) Spectra. The interference term is plotted on its own axis, since it is
% a small ripple on the classical line.
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
grid on;
xlim([-1 1]*lineWidth_Hz/1e6);

% (b) Matched filter response, wide enough to show the image peak.
subplot(1,2,2);
uax = (1:cfg.N) - cfg.N/2 - 1;
Y0_mf = real(fftshift(fft(ifftshift( sig.tmpl.*weight))));
Y1_mf = real(fftshift(fft(ifftshift(-sig.tmpl.*weight))));
plot(uax, Y0_mf/max(abs(Y0_mf)), 'b-', 'LineWidth', 1.4); hold on;
plot(uax, Y1_mf/max(abs(Y0_mf)), 'r--', 'LineWidth', 1.4);
yline(0, 'k:');
xlabel('Baseband Bin Index');
ylabel('Normalized Output');
title('(b) Compressed Matched Filter Peak');
legend({'bit 0 (matched)','bit 1 (matched)'}, 'Location', 'best', 'FontSize', 8);
grid on;
xlim([-60 max(120, ceil(1.2*imageOffset))]);

% -------------------------------------------------------------------------
% FIG 1c: WATERFALL, clipped at the resolvable floor 1/(Nbit_frame*N_MC)
% -------------------------------------------------------------------------
figure('Name','Fig 1c Error rate validation','Color','w','Position',[120 120 620 460]);
BER_floor = 1/(cfg.Nbit_frame*cfg.N_MC);
bp = BER_A;
bp(bp==0) = NaN;
errorbar(sigma_vec, bp, bp-BER_A_lo, BER_A_hi-bp, 'o', 'LineWidth', 1.2, 'MarkerSize', 5); hold on;
plot(sigma_vec, BER_theory_A, '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.3);
yline(BER_floor, 'k-.', 'resolvable floor', 'LineWidth', 1.1, 'LabelOrientation', 'horizontal');
set(gca, 'XScale', 'log', 'YScale', 'log');
ylim([BER_floor/3 1]);
xlim([min(sigma_vec) max(sigma_vec)]);
xline(min(opt.sigma_ep), 'g--', '\sigma closest', 'LineWidth', 1.1, 'LabelOrientation', 'horizontal');
xline(max(opt.sigma_ep), 'm--', '\sigma edge', 'LineWidth', 1.1, 'LabelOrientation', 'horizontal');
xlabel('Relative Noise Floor \sigma');
ylabel('Raw BER');
title('Error Rate Validation');
legend({'Monte Carlo', 'theoretical filter curve'}, 'Location', 'southeast', 'FontSize', 8);
grid on;

% -------------------------------------------------------------------------
% FIG 2a: BER AND CAPACITY ALONG THE PASS
% -------------------------------------------------------------------------
figure('Name','Fig 2a Channel quality over the pass','Color','w','Position',[110 80 620 460]);
yyaxis left;
errorbar(t_ep, BER_B, BER_B-BER_B_lo, BER_B_hi-BER_B, 'o-', 'Color', [0 0.45 0.74], 'LineWidth', 1.3, 'MarkerSize', 4);
ylabel('Raw BER');
ylim([0 0.55]);
ax = gca; ax.YColor = [0 0.45 0.74];
yyaxis right;
plot(t_ep, C_ep_sim, 's--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2, 'MarkerSize', 4);
ylabel('Capacity [bit/symbol]');
ylim([0 1.05]);
ax.YColor = [0.85 0.33 0.1];
xline(t_sec(idx_cross), 'k:', 'Closest approach', 'LabelOrientation', 'horizontal');
xlabel('Time from Pass Start [s]');
title('Channel Quality Over Fly-By Pass');
grid on;

% -------------------------------------------------------------------------
% FIG 2b: CUMULATIVE VOLUME AND COVERT CEILINGS
% The cumulative volume is drawn at the END of each epoch, so the black
% curve meets each adversary line exactly at its stop time.
% -------------------------------------------------------------------------
figure('Name','Fig 2b Covert information volume','Color','w','Position',[140 100 620 460]);
ax = gca; hold(ax, 'on');
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
    yline(ax, vols(i), styles{i}, 'Color', cols{i}, 'LineWidth', 1.3, 'DisplayName', names{i});
    if isnan(tStopP(i)), continue; end       % threshold never reached
    k = kStop(i);
    % Fraction of the crossing epoch already transmitted at the stop: the
    % same fraction volumeAtThreshold uses (D^2 and bits linear in symbols).
    if k < N_epochs, span = t_ep(k+1) - t_ep(k); else, span = t_ep(k) - t_ep(k-1); end
    frac = (tStopP(i) - t_ep(k))/span;
    tS   = (k - 1 + frac)*dtE;
    plot(ax, [tS tS], [0 vols(i)], ':', 'Color', cols{i}, 'LineWidth', 1.1, 'HandleVisibility', 'off');
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
title(ax, 'Covert Information Volume');
legend(ax, 'Location', 'northwest', 'FontSize', 7.5);
grid(ax, 'on');  hold(ax, 'off');

% -------------------------------------------------------------------------
% FIG 3a: OPERATING POINT SWEEP
% -------------------------------------------------------------------------
figure('Name','Fig 3a Operating point optimisation','Color','w','Position',[150 110 620 460]);
semilogy(sep_vec, max(B_noEve,1e-3), '--', 'LineWidth', 1.2); hold on;
semilogy(sep_vec, max(B_naive,1e-3), 'LineWidth', 1.2);
semilogy(sep_vec, max(B_hop,1e-3), 'LineWidth', 1.6);
semilogy(sep_vec, max(B_keyed,1e-3), 'LineWidth', 1.2);
xline(cfg.sep_over_Delta, 'k:', 'optimum', 'LineWidth', 1.2);
xlabel('Wavepacket Separation [\Delta_p]');
ylabel('Pass Volume [bit]');
title('Operating Point Optimization');
legend({'no eavesdropper', 'vs Eve-A', 'vs Eve-C (hop-blind)', 'vs Eve-B'}, ...
       'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 7);
grid on;

% -------------------------------------------------------------------------
% FIG 3b: IMPRINT DEPTH SENSITIVITY, valid range (OD <= 1) highlighted
% -------------------------------------------------------------------------
figure('Name','Fig 3b Imprint depth sensitivity','Color','w','Position',[160 120 620 460]);
ax = gca; hold(ax, 'on');
OD_fine = logspace(log10(cfg.OD/20), log10(cfg.OD*8), 300);
d_fine  = (OD_fine/cfg.OD)*(g_mf/sigma_ref);
C_fine  = 1 - binaryEntropy(0.5*erfc(d_fine/sqrt(2)));
inR     = OD_fine <= 1;
patch(ax, [1 max(OD_fine) max(OD_fine) 1], [0 0 1 1], [0.9 0.9 0.9], ...
      'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, OD_fine(inR),  C_fine(inR),  '-',  'Color', [0 0.45 0.74], 'LineWidth', 1.6, ...
     'DisplayName', 'linear imprint (OD \leq 1)');
plot(ax, OD_fine(~inR), C_fine(~inR), '--', 'Color', [0 0.45 0.74], 'LineWidth', 1.2, ...
     'DisplayName', 'extrapolation');
d_dec = g_mf/sigma_ref;
plot(ax, cfg.OD, 1 - binaryEntropy(0.5*erfc(d_dec/sqrt(2))), 'ko', 'MarkerFaceColor', 'k', ...
     'HandleVisibility', 'off');
xline(ax, cfg.OD, 'k:', 'declared', 'LineWidth', 1.2, 'HandleVisibility', 'off');
set(ax, 'XScale', 'log', 'Layer', 'top');
xlim(ax, [min(OD_fine) max(OD_fine)]);  ylim(ax, [0 1]);
xlabel(ax, 'Imprint Depth [OD]');
ylabel(ax, 'Capacity [bit/symbol]');
title(ax, 'Imprint Depth Sensitivity');
legend(ax, 'Location', 'northwest', 'FontSize', 7.5);
grid(ax, 'on');  hold(ax, 'off');

% -------------------------------------------------------------------------
% FIG 4: HOPPING ALPHABET AND ITS CEILINGS
% -------------------------------------------------------------------------
figure('Name','Fig 4 Hopping alphabet','Color','w','Position',[170 40 900 400]);

% (a) Covert volume against Eve-C versus alphabet size, from the simulated
% per-epoch deflections at the fixed operating point.
subplot(1,2,1);
M_vec = unique(round(logspace(0, log10(opt.N_env), 45)));   % beyond N_env Eve-C would use the energy detector
B_vec = zeros(size(M_vec));
for k = 1:numel(M_vec)
    dC = d_ep_sim.^2/sqrt(2*M_vec(k));
    eC = sqrt(cumsum(n_sym_ep.*dC.^2));
    B_vec(k) = volumeAtThreshold(bits_cum, eC, SNR_thr_Eve, n_sym_ep, C_ep_sim, dC);
end
loglog(M_vec, max(B_vec,eps), 'LineWidth', 1.6); hold on;
% Slope-1/2 reference (re-optimised amplitude, equal epochs).
ref = max(B_vec(1),eps)*sqrt(M_vec/M_vec(1));
loglog(M_vec, ref, 'k:', 'LineWidth', 1.1);
xline(max(opt.M_eff,1), 'b-.', 'M resolvable', 'LineWidth', 1.3, 'LabelOrientation', 'horizontal');
xline(opt.N_env, 'g--', 'Eve-A (M = N_{env})', 'LineWidth', 1.1, 'LabelOrientation', 'horizontal');
% Marker at the alphabet allowed by the sampling bound alone.
if cfg.hopping
    B_nyq = 10^interp1(log10(M_vec), log10(max(B_vec,eps)), log10(Mres_nyq));
    plot(Mres_nyq, B_nyq, 'kd', 'MarkerFaceColor', [0.85 0.33 0.10], 'HandleVisibility', 'off');
    text(Mres_nyq, B_nyq, sprintf('  \\sigma_f = %.0f Hz', cfg.sigma_f_track*opt.F_trk/opt.F_nyq), ...
         'FontSize', 8, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
end
xlabel('Resolvable Alphabet M [independent pitches]');
ylabel('Covert Volume vs Eve-C [bit]');
title('(a) Covert Volume Bought with Alphabet Size');
legend({'simulated', 'slope 1/2 reference'}, 'Location', 'southeast', 'FontSize', 8);
grid on; xlim([1 opt.N_env]);

% Fitted exponent, only where detection binds.
useFit = B_vec > 0 & B_vec < 0.95*bits_cum(end);
if sum(useFit) >= 5
    pf = polyfit(log10(M_vec(useFit)), log10(B_vec(useFit)), 1);
    fprintf('\n--- 9.1 Hopping alphabet scaling ---\n');
    fprintf('  covert volume at FIXED operating point grows as M^%.3f over %d points\n', ...
            pf(1), sum(useFit));
    fprintf('  (not expected to be exactly 0.5: the amplitude is not re-optimised for each M\n');
    fprintf('   and the statistic accumulates chronologically over very unequal epochs, so the\n');
    fprintf('   truncation point moves across epochs of different capacity. The square-root\n');
    fprintf('   law itself is tested in Section 8.2, where the amplitude is re-optimised.)\n');
    if pf(1) <= 0 || pf(1) > 1
        warning('Covert volume vs alphabet size has a non-physical slope: check the Eve-C model.');
    end
else
    fprintf('\n--- 9.1 Hopping alphabet scaling: too few unsaturated points to fit ---\n');
end
if cfg.hopping
    fprintf('  volume vs Eve-C read at the sampling-limited alphabet M = %.1f -> %.2f kbit (fixed operating point)\n', ...
            Mres_nyq, B_nyq/1e3);
end

% --- [P3] What limits the covert volume against Eve-C ---------------------
% Each factor is the gain obtained by relaxing one element, at the same
% operating point unless stated. The alphabet (set by the Doppler tracking)
% is one factor among several, not the dominant one.
fprintf('\n--- 9.2 [P3] Factors on the covert volume vs Eve-C (D = %.0f: %.3e bit) ---\n', ...
        SNR_thr_Eve, B_covert_C);
if cfg.hopping
    fprintf('  alphabet, tracking to the sampling limit (M %.1f -> %.1f) : x%.2f\n', ...
            opt.M_eff, Mres_nyq, B_nyq/max(B_covert_C,eps));
end
fprintf('  alphabet, upper limit of any keyed scheme (Eve-A)      : x%.2f\n', B_covert_A/max(B_covert_C,eps));
fprintf('  detection threshold, D = %.0f instead of D = 1           : x%.2f\n', ...
        SNR_thr_Eve, B_covert_C/max(B_covert_C1,eps));
fprintf('  pass geometry, stationary link at closest range (design): x%.2f\n', Bs_best/max(B_best,eps));
fprintf('  window four times longer (design, re-optimised)        : x%.2f\n', B_scale(2)/max(B_scale(1),eps));
fprintf('  [P9] Eve without the heterodyne penalty                : x%.2f lower\n', B_covert_C/max(B_C_nh3,eps));

% (b) Alphabet ceiling versus Doppler tracking residual.
subplot(1,2,2);
sf_vec  = logspace(1, 4, 200);
phi_tol = sqrt(2*log(10^(cfg.hopLossTol_dB/20)));
F_trk_v = (phi_tol./(2*pi*sf_vec))*cfg.B_IF;
M_v     = resolvablePitches(max(1, floor(min(F_trk_v, opt.F_nyq)) - opt.F_min + 1), opt.rhoA);
loglog(sf_vec, M_v, 'LineWidth', 1.6); hold on;
yline(resolvablePitches(floor(opt.F_nyq)-opt.F_min+1, opt.rhoA), 'k--', 'Nyquist ceiling', 'LineWidth', 1.1, 'LabelOrientation', 'horizontal');
xline(cfg.sigma_f_track, 'r:', '\sigma_f declared', 'LineWidth', 1.3, 'LabelOrientation', 'horizontal');
xlabel('Doppler Tracking Residual \sigma_f [Hz]');
ylabel('Resolvable Alphabet M');
title(sprintf('(b) Alphabet Ceiling (budget %.1f dB)', cfg.hopLossTol_dB));
grid on;

% --- File names for export ------------------------------------------------
figNames = {'Fig 1ab Spectral fringe and matched filtering', 'Fig 1c Error rate validation', ...
            'Fig 2a Channel quality over the pass', 'Fig 2b Covert information volume', ...
            'Fig 3a Operating point optimisation', 'Fig 3b Imprint depth sensitivity', ...
            'Fig 4 Hopping alphabet'};
fprintf('\nFigures to export:\n');
for i = 1:numel(figNames)
    fprintf('  %s  ->  %s.pdf\n', figNames{i}, strrep(figNames{i}, ' ', '_'));
end

fprintf('\n--- SIMULATION COMPLETE ---\n');
%% ========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function out = operatingPoint(sep, par)
% Design point for one packet separation: spectra, Bob's deflection and
% capacity per epoch, hopping alphabet, adversary deflections and volumes.
    N = par.N;  Delta = par.Delta;  kappa = par.kappa;
    % Band: peaks at -+ sep*Delta/2, six standard deviations of margin each
    % side. The Lorentzian of Eq. (26) is treated as a delta (Gamma0 is four
    % orders below the fringe spacing).
    B_IF  = 4*par.f0*(sep+6)*Delta;
    T_sym = N/B_IF;
    dw = linspace(-1,1,N).'*(2*pi*B_IF/2);
    p  = dw/par.Omega_rad;                  % velocity in units of v/c
    % Momentum amplitudes of Eq. (A1), centred at -+ sep*Delta/2.
    A1 = exp(-((p+sep*Delta/2).^2)./(2*Delta^2))./(pi^0.25*sqrt(Delta));
    A2 = exp(-((p-sep*Delta/2).^2)./(2*Delta^2))./(pi^0.25*sqrt(Delta));
    % Emission weight of Eq. (26) (a ~1e-8 correction, kept for fidelity).
    w  = 1 + 3*p;
    % Incoherent part [Eq. A2] and cross term. The bit-dependent
    % normalisation of Eq. (29) is omitted: the probe power is fixed by the
    % laser.
    S_cl  = (cos(par.theta)^2*A1.^2 + sin(par.theta)^2*A2.^2).*w;
    I_env = (2*cos(par.theta)*sin(par.theta)*A1.*A2).*w;
    % Normalised to unit rms and scaled by the declared optical depth.
    nrm   = sqrt(mean(S_cl.^2));
    OD    = par.OD;
    S_cl  = OD*S_cl/nrm;
    I_env = OD*I_env/nrm;
    Phi   = kappa*dw;

    signalEnergy = sum((I_env.*cos(Phi)).^2);          % ||s||^2 at unit sigma

    n_sym    = par.dt_ep(:).'/T_sym;
    n_bin_ep = (par.P_ep(:).'*T_sym/par.E_photon)*par.eta_det/N;
    sigma_ep = par.F_excess./sqrt(max(n_bin_ep,eps));
    d_ep     = sqrt(signalEnergy)./sigma_ep;
    BER_ep   = 0.5*erfc(d_ep/sqrt(2));
    C_ep     = 1 - binaryEntropy(BER_ep);
    % Eve-A: envelope-weighted energy; N_env is the participation ratio of
    % I_env^2.
    N_env    = sum(I_env.^2)^2/sum(I_env.^4);
    % [P10] The adversaries average the fading on L^2: d^2 x eveGain.
    d2_eve   = d_ep.^2*par.eveGain;
    delta_A  = d2_eve/sqrt(2*N_env);                   % key-unaware detector
    delta_B  = d2_eve/sqrt(2);                         % key-aware detector

    % --- Fringe-hopping alphabet ------------------------------------------
    % Three ceilings bound the range of the fringe count F = kappa*B_IF; the
    % number of distinct pitches is then computed from template correlation.
    F_now = kappa*B_IF;
    % (a) Nyquist: fringe phase per sample below pi.
    F_nyq = (N-1)/2;
    % (b) Image clearance: the image at 2F bins must clear the compressed
    %     peak (same half-width as Section 4.12).
    sig_env = (Delta/sqrt(2))*par.Omega_rad;
    dwStep  = dw(2)-dw(1);
    halfW   = N/(2*pi*((sig_env/sqrt(2))/dwStep));
    F_min   = ceil(par.imageMargin*halfW/2);
    % (c) Doppler tracking residual: averaged loss exp(-sigma_phi^2/2) equal
    %     to the declared budget.
    phi_tol = sqrt(2*log(10^(par.hopLossTol_dB/20)));
    F_trk   = (phi_tol/(2*pi*par.sigma_f))*B_IF;
    M_alpha = max(0, floor(min(F_nyq,F_trk)) - F_min + 1);   % integer counts
    % Resolvable pitches: templates F and F+k correlate as
    % rho(k) = exp(-(k sigW/B_IF)^2/2), sigW the rms width of I_env^2; the
    % number of independent hypotheses is the participation ratio of the
    % correlation matrix.
    wI    = I_env.^2/sum(I_env.^2);
    sigW  = sqrt(sum(wI.*dw.^2) - sum(wI.*dw)^2);
    rhoA  = sigW/B_IF;                      % correlation decay per unit of F
    M_res = resolvablePitches(M_alpha, rhoA);
    if par.hopping
        M_eff = min(max(M_res,1), N_env); % capped at N_env: beyond it Eve-C
    else                                  % would simply use the energy
        M_eff = 1;                        % detector instead
    end
    delta_C = d2_eve/sqrt(2*M_eff);       % scheme-aware, key-unaware

    % Centroid leakage: calibration of theta required of the source.
    dp_pk = sep*Delta;
    dQ = @(ph,th) cos(ph)*sin(4*th)*dp_pk / ...
                  (4*(cos(ph)*sin(2*th) + exp(dp_pk^2/(4*Delta^2))));
    th_e     = par.theta + par.eps_theta_probe;
    swing_Hz = abs(dQ(0,th_e) - dQ(pi,th_e))*par.f0;
    lineWidthEff_Hz = par.f0*Delta*(sep/2 + 1);
    n_ph_ep  = (par.P_ep(:).'*T_sym/par.E_photon)*par.eta_det;
    sig_cen  = par.F_excess*lineWidthEff_Hz./sqrt(max(n_ph_ep,eps));
    snr_cen  = sqrt(sum(n_sym.*((swing_Hz/2)./sig_cen).^2));
    out.eps_theta_req = par.eps_theta_probe*par.SNR_thr/max(snr_cen,eps);

    bits_cum = cumsum(n_sym.*C_ep);
    eveA     = sqrt(cumsum(n_sym.*delta_A.^2));
    eveB     = sqrt(cumsum(n_sym.*delta_B.^2));
    eveC     = sqrt(cumsum(n_sym.*delta_C.^2));

    out.B_IF = B_IF;  out.T_sym = T_sym;  out.dw = dw;
    out.S_cl = S_cl;  out.I_env = I_env;  out.Phi = Phi;
    out.sigma_ep = sigma_ep;  out.d_ep = d_ep;
    out.BER_ep = BER_ep;      out.C_ep = C_ep;
    out.M_pass = sum(n_sym);
    out.volume_noEve = bits_cum(end);
    out.volume_naive = volumeAtThreshold(bits_cum, eveA, par.SNR_thr, n_sym, C_ep, delta_A);
    out.volume_keyed = volumeAtThreshold(bits_cum, eveB, par.SNR_thr, n_sym, C_ep, delta_B);
    out.volume_hopped = volumeAtThreshold(bits_cum, eveC, par.SNR_thr, n_sym, C_ep, delta_C);
    out.M_alpha = M_alpha;   out.M_eff = M_eff;
    out.M_res   = M_res;     out.rhoA  = rhoA;
    out.N_env   = N_env;
    out.F_now   = F_now;     out.F_min = F_min;
    out.F_nyq   = F_nyq;     out.F_trk = F_trk;
    out.halfW   = halfW;
    % Cells of the targeted adversary, for the closed-form optimum (4.6).
    if par.hopping, out.M_eff_target = M_eff; else, out.M_eff_target = N_env; end
    out.kappa   = kappa;
    out.nyquistPhase = kappa*(dw(2)-dw(1));
end

function B = volumeAtThreshold(bits_cum, eve, thr, n_sym, C_ep, delta)
% Information volume delivered before the accumulated statistic crosses
% threshold; linear interpolation in accumulated deflection energy inside
% the crossing epoch.
    k = find(eve > thr, 1);
    if isempty(k)
        B = bits_cum(end);
        return;
    end
    if k == 1
        prevEnergy = 0; prevBits = 0;
    else
        prevEnergy = eve(k-1)^2; prevBits = bits_cum(k-1);
    end
    quota = n_sym(k)*delta(k)^2;
    if quota <= 0
        B = prevBits;
    else
        frac = min(max((thr^2 - prevEnergy)/quota, 0), 1);
        B = prevBits + frac*n_sym(k)*C_ep(k);
    end
end

function [tStop, k] = stopEpochTime(eve, thr, t_ep)
% Time at which the accumulated statistic crosses threshold (same rule as
% volumeAtThreshold); NaN if never reached.
    k = find(eve > thr, 1);
    if isempty(k), tStop = NaN; return; end
    if k == 1, e0 = 0; else, e0 = eve(k-1)^2; end
    frac = (thr^2 - e0)/max(eve(k)^2 - e0, eps);
    if k < numel(t_ep), dtk = t_ep(k+1)-t_ep(k); else, dtk = t_ep(k)-t_ep(k-1); end
    tStop = t_ep(k) + frac*dtk;
end

function printStop(name, ts, k, idx_ep, d_t)
    if isnan(ts)
        fprintf('  %s never reaches threshold within the pass\n', name);
    else
        fprintf('  %s: epoch %d, t ~ %.0f s, range at that epoch %.0f km\n', ...
                name, k, ts, d_t(idx_ep(k))/1e3);
    end
end

function Mr = resolvablePitches(Mint, rhoA)
% Independent pitch hypotheses in Mint consecutive fringe counts:
% participation ratio of the Toeplitz matrix rho(|i-j|) =
% exp(-((i-j) rhoA)^2/2). Accepts a vector.
    Mr = zeros(size(Mint));
    for q = 1:numel(Mint)
        M = floor(Mint(q));
        if M <= 1, Mr(q) = max(M,0); continue; end
        k  = 1:M-1;
        r2 = exp(-(k*rhoA).^2);            % rho(k)^2
        Mr(q) = M^2/(M + 2*sum((M-k).*r2));
    end
end

function [nerr,nbit] = runFrame(cfg, sig, sigma_noise, rs)
% Dispatcher: sufficient-statistic path if cfg.fastMC, full spectral path
% otherwise.
    if cfg.fastMC
        [nerr,nbit] = simulateFrameFast(cfg, sig, sigma_noise, rs);
    else
        [nerr,nbit] = simulateFrame(cfg, sig, sigma_noise, rs);
    end
end

function [nerr,nbit] = simulateFrame(cfg, sig, sigma_noise, rs)
% One frame on the full spectral path. Bit = phi in {0, pi} (antipodal).
% Each symbol draws its pitch from cfg.F_alpha; the residual Doppler offset
% rotates its fringe by 2*pi*kappa_n*s. Calibration errors are frozen over
% the frame. Decision: sign of the zero-frequency bin of
% fft(Rx .* I_env .* key).
% [P13] Bob's copy S_hat carries the per-bin error eps_sub and the smooth
% errors of frameProfile. With cfg.profileFit, Bob fits gain, centre and
% width to the average of the frame and removes them: since the statistic is
% linear, this equals subtracting the corrected profile before the FFT.
    M = cfg.Nbit_frame; N = cfg.N;
    Fa = cfg.F_alpha(:).';  Q = numel(Fa);
    bits = randi(rs,[0 1],1,M);

    S_hat = sig.S_cl.*(1 + cfg.eps_sub*randn(rs,N,1));   % per-bin error
    eps_k = cfg.eps_kappa*randn(rs,1,1);                  % error on kappa
    S_hat = S_hat + frameProfile(cfg, sig, rs) - sig.S_cl; % smooth errors

    statRaw = zeros(1,M);  XW = zeros(3,M);  Rsum = zeros(N,1);
    for i1 = 1:cfg.block:M
        i2 = min(i1+cfg.block-1,M); b = bits(i1:i2); nb = numel(b);
        s  = 1-2*b;                                   % bit 0 -> +1, bit 1 -> -1
        if Q > 1, iq = randi(rs,Q,1,nb); else, iq = ones(1,nb); end
        kap  = Fa(iq)/cfg.B_IF;                       % pitch of each symbol
        dphi = 2*pi*cfg.sigma_f_track*kap.*randn(rs,1,nb);   % Doppler rotation
        Tx   = sig.I_env.*cos(kap.*cfg.dw + dphi);    % N x nb
        W    = sig.I_env.*exp(+1i*(kap*(1+eps_k)).*cfg.dw);
        Rx = sig.S_cl + Tx.*s + sigma_noise*randn(rs,N,nb);
        Rx = Rx - S_hat;
        % ifftshift: the signal is centred at zero detuning (array centre).
        Y  = fftshift(fft(ifftshift(Rx.*W,1),[],1),1);
        statRaw(i1:i2) = real(Y(cfg.kpk,:));
        XW(:,i1:i2)    = real(sig.Xn.'*W);            % projection of the basis
        Rsum = Rsum + sum(Rx,2);
    end
    stat = statRaw;
    if cfg.profileFit
        thetaHat = sig.XtXi*(sig.Xn.'*(Rsum/M));       % fit on the frame average
        stat = statRaw - thetaHat.'*XW;
    end
    nerr = sum(double(stat<0) ~= bits);
    nbit = M;
end

function [nerr,nbit] = simulateFrameFast(cfg, sig, sigma_noise, rs)
% Same model as simulateFrame on the sufficient statistic (exact; see
% cfg.fastMC). Bias, in-phase and quadrature projections and noise scale are
% precomputed once per frame for every pitch. Uses a different number of
% random variates, so realisations differ while the statistics agree.
% [P13] The fit of the frame average is also reduced exactly: its input is
% X'(S_cl - S_hat) + (1/M) sum_n X'(signal_n) + (sigma/M) sum_n X'w_n, and
% the noise of each symbol is drawn jointly with its projection on X
% (4-dimensional Gaussian, covariance [W'W W'X; X'W X'X]), so the
% correlation between the fit and each symbol's own noise is kept.
    M = cfg.Nbit_frame;  N = cfg.N;
    Fa = cfg.F_alpha(:).';  Q = numel(Fa);
    bits = randi(rs,[0 1],1,M);

    S_hat = sig.S_cl.*(1 + cfg.eps_sub*randn(rs,N,1));
    eps_k = cfg.eps_kappa*randn(rs,1,1);
    S_hat = S_hat + frameProfile(cfg, sig, rs) - sig.S_cl;

    K    = Fa/cfg.B_IF;                                        % 1 x Q
    Wr   = real(sig.I_env.*exp(+1i*(K*(1+eps_k)).*cfg.dw));    % N x Q
    bias = sum((sig.S_cl - S_hat).*Wr, 1);
    A    = sum(sig.I_env.*cos(K.*cfg.dw).*Wr, 1);
    Aq   = sum(sig.I_env.*sin(K.*cfg.dw).*Wr, 1);
    XWq  = sig.Xn.'*Wr;                                        % 3 x Q
    Gc   = sig.Xn.'*(sig.I_env.*cos(K.*cfg.dw));               % 3 x Q
    Gs   = sig.Xn.'*(sig.I_env.*sin(K.*cfg.dw));               % 3 x Q
    XtX  = sig.Xn.'*sig.Xn;
    Lq   = zeros(4,4,Q);
    for q = 1:Q
        C = [sum(Wr(:,q).^2), XWq(:,q).'; XWq(:,q), XtX];
        [L,pf] = chol(C,'lower');
        if pf, L = chol(C + 1e-12*trace(C)*eye(4),'lower'); end
        Lq(:,:,q) = L;
    end

    s = 1-2*bits;
    if Q > 1, iq = randi(rs,Q,1,M); else, iq = ones(1,M); end
    dphi = 2*pi*cfg.sigma_f_track*K(iq).*randn(rs,1,M);
    z = randn(rs,4,M);  u = zeros(4,M);
    for q = 1:Q
        m = (iq == q);
        if any(m), u(:,m) = Lq(:,:,q)*z(:,m); end
    end
    sigTerm = s.*(A(iq).*cos(dphi) - Aq(iq).*sin(dphi));
    stat = bias(iq) + sigTerm + sigma_noise*u(1,:);
    if cfg.profileFit
        sigProj  = Gc(:,iq).*(s.*cos(dphi)) - Gs(:,iq).*(s.*sin(dphi));   % 3 x M
        meanProj = sig.Xn.'*(sig.S_cl - S_hat) + mean(sigProj,2) + sigma_noise*mean(u(2:4,:),2);
        thetaHat = sig.XtXi*meanProj;
        stat = stat - thetaHat.'*XWq(:,iq);
    end
    nerr = sum(double(stat<0) ~= bits);
    nbit = M;
end

function S = frameProfile(cfg, sig, rs)
% [P13] Bob's copy of the classical profile with the smooth errors of one
% frame: gain, centre (Hz) and lobe width. Random (Gaussian, standard
% deviations cfg.eps_*) or, with cfg.smoothFixed, of fixed size.
    g = randn(rs,1,3);
    if cfg.smoothFixed, g = [1 1 1]; end
    S = profileOf(sig.pp, cfg.eps_gain*g(1), cfg.eps_c_Hz*g(2), cfg.eps_w*g(3));
end

function S = profileOf(pp, eg, ec_Hz, ew)
% [P13] Classical profile of operatingPoint rebuilt with a gain error eg, a
% centre shift ec_Hz and a relative error ew on the lobe width (amplitude
% normalisation kept, so the area of each lobe is unchanged; centres fixed).
    p  = (pp.dw - 2*pi*ec_Hz)/pp.Om;
    D  = pp.Delta*(1+ew);
    A1 = exp(-((p+pp.sep*pp.Delta/2).^2)./(2*D^2))./(pi^0.25*sqrt(D));
    A2 = exp(-((p-pp.sep*pp.Delta/2).^2)./(2*D^2))./(pi^0.25*sqrt(D));
    S  = pp.kS*(1+eg)*(cos(pp.theta)^2*A1.^2 + sin(pp.theta)^2*A2.^2).*(1 + 3*p);
end

function rs = mcStream(seed,campaign,point,trial)
% Independent, reproducible substream per campaign, sweep point and trial.
    rs = RandStream('mrg32k3a','Seed',seed);
    rs.Substream = 1 + trial + 1000*point + 1000000*campaign;
end

function [lo,hi] = wilsonCI(k,n,z)
% Wilson score interval for a binomial proportion.
    p = k./n; d = 1+z^2./n; c = p+z^2./(2*n);
    s = z*sqrt(p.*(1-p)./n + z^2./(4*n.^2));
    lo = max((c-s)./d,0); hi = min((c+s)./d,1);
end

function I = pitchInfo(dGrid, M, nT, rs)
% [P12] Mutual information (bit) between a uniform pitch index among M
% orthogonal cells and the M real matched-filter outputs, for deflection d;
% the bit sign is unknown, so the likelihood of cell j is cosh(d y_j).
    I = zeros(size(dGrid));
    for k = 1:numel(dGrid)
        d = dGrid(k);
        y = randn(rs, nT, M);  y(:,1) = y(:,1) + d;     % true pitch in cell 1
        x = abs(d*y);
        lw = x + log1p(exp(-2*x));                      % log(2 cosh(d y))
        lw = lw - max(lw,[],2);
        w  = exp(lw);  w = w./sum(w,2);
        H  = -sum(w.*log2(max(w,realmin)),2);
        I(k) = log2(M) - mean(H);
    end
    I = max(I, 0);
end

function H = binaryEntropy(p)
    p = min(max(p,0),1); H = zeros(size(p)); m = (p>0)&(p<1);
    H(m) = -p(m).*log2(p(m)) - (1-p(m)).*log2(1-p(m));
end