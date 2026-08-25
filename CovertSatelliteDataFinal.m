% =========================================================================
% QUANTUM-COVERT INTER-SATELLITE LINK SIMULATOR
% Spectral chirp encoding on the quantum Doppler effect, fractional-domain
% de-chirping, and covert-capacity assessment over a LEO-LEO pass.
%
% PHYSICAL BASIS
%   P. T. Grochowski, A. R. H. Smith, A. Dragan, K. Debski, "Quantum time
%   dilation in atomic spectra", Phys. Rev. Research 3, 023053 (2021).
%   Equation numbers quoted in the code refer to that work.
%
%   The scheme extends the reference by introducing a differential
%   preparation delay Delta_t between the two branches of the momentum
%   superposition. Without that extension the interference term is a single
%   modulated lobe and not a chirp; see Note 2 in the closing section.
%
% MODEL CHAIN (every quantity is derived; free parameters are listed in
% Note 3 of the closing section)
%   1. A 27Al+ ion is prepared in a coherent superposition of two momentum
%      wave packets [Eq. A1], the branches being prepared Delta_t apart.
%   2. Free evolution in the momentum representation acts as
%      exp(-i p^2 t / 2 m hbar), so the branch-to-branch phase difference is
%      Phi_2(p) - Phi_1(p) = -p^2 Delta_t / (2 m hbar): quadratic in p.
%   3. The first-order Doppler map dw = Omega p/(mc) turns that into
%      Phi(dw) = -C dw^2 with C = m c^2 Delta_t / (2 hbar Omega^2). This is
%      the spectral chirp; its coefficient is not a tunable parameter.
%   4. The emitted spectrum is
%         S_quant(w) = S_class(w) + I_env(w) cos(Phi(w) + phi)
%      where S_class is the incoherent statistical mixture [Eq. A2], I_env
%      is the interference envelope, and phi in {0, pi} carries the bit.
%   5. The receiver removes the known classical profile, applies the
%      envelope-weighted de-chirp key, and reads the sign of the compressed
%      peak. This is the matched filter for the real template.
%   6. Averaged over equiprobable bits the interference term cancels,
%      cos(x) + cos(x+pi) = 0, and only S_class survives.
%
% REQUIRED TOOLBOXES
%   Satellite Communications Toolbox (satelliteScenario, satellite, states)
%   Phased Array System Toolbox (fspl)
%
% RUNTIME
%   Dominated by Gaussian variate generation in the Monte Carlo campaigns
%   (order 1e10 samples through the mrg32k3a substream generator), not by
%   the FFTs. Expect several minutes. See Note 9 for reduction strategies
%   that do not compromise the statistical validity of the capacity figure.
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
mc2_Al27   = m_Al27*c_light^2;
R_earth    = 6378137;
mu_earth   = 3.986004418e14;
Gamma0_Hz  = 520;                   % natural linewidth of the transition

fprintf('=== 1. CONSTANTS ===\n');
fprintf('27Al+ 1S0-3P1: lambda = %.1f nm, f0 = %.3f THz, Gamma0 = %d Hz\n', ...
        lambda_0*1e9, f0_optical/1e12, Gamma0_Hz);

%% ========================================================================
% 2. ORBITAL GEOMETRY (SGP4 PROPAGATION, J2-CONSISTENT INITIALISATION)
% =========================================================================
% Two co-altitude LEO satellites on planes separated so that their relative
% geometry produces a close fly-by. The RAAN offset and the true anomalies
% are solved for numerically so that the minimum separation lands near
% 300 km at mid-window; SGP4 then propagates the resulting elements.
fprintf('\n=== 2. ORBITAL GEOMETRY ===\n');

startTime  = datetime(2026,8,14,10,0,0);
stopTime   = startTime + minutes(15);
sampleTime = 0.1;
T_window   = seconds(stopTime - startTime);

altitude_orbit = 800e3;  sma = R_earth + altitude_orbit;
ecc = 0.0001; incl_deg = 60; argPeri_deg = 0; J2 = 1.08263e-3;
mu_n = sqrt(mu_earth/sma^3);
raan_rate = -1.5*mu_n*J2*(R_earth/(sma*(1-ecc^2)))^2*cosd(incl_deg);
angolo_piani_deg = 60;   % angle between the two planes
dRAAN_deg = acosd((cosd(angolo_piani_deg)-cosd(incl_deg)^2)/sind(incl_deg)^2);
raan_Alice = 0; raan_Bob = dRAAN_deg;
T_orbit = 2*pi/mu_n; incl_rad = deg2rad(incl_deg);

posC = @(u0,raan0,t) sma*[ ...
  cos(raan0+raan_rate*t).*cos(u0+mu_n*t)-sin(raan0+raan_rate*t).*sin(u0+mu_n*t).*cos(incl_rad); ...
  sin(raan0+raan_rate*t).*cos(u0+mu_n*t)+cos(raan0+raan_rate*t).*sin(u0+mu_n*t).*cos(incl_rad); ...
  sin(u0+mu_n*t).*sin(incl_rad)];

% Coarse then fine search for the argument of latitude giving a 300 km pass.
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

% Solar exclusion: reject epochs with the Sun within 30 deg of the receive
% line of sight.
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

%% ========================================================================
% 3. FREE-SPACE OPTICAL LINK BUDGET
% =========================================================================
% The probe beam that carries the imprinted spectral response is a
% conventional FSO link. Pointing, acquisition and tracking jitter is
% applied independently at both terminals.
fprintf('\n=== 3. OPTICAL LINK BUDGET ===\n');
D_tx = 0.20; D_rx = 0.30; P_tx = 100e-3;
eta_tx = 0.85; eta_rx = 0.80; eta_det = 0.65; rho_ap = 0.65;
sigma_jit_pat = 1.5e-6; theta_div = 8.0e-6;

G_tx = rho_ap*(pi*D_tx/lambda_0)^2;
G_rx = rho_ap*(pi*D_rx/lambda_0)^2;
FSPL_t = 10.^(fspl(d_t,lambda_0)/10);
patLoss = @(n) exp(-8*((sqrt((sigma_jit_pat*randn(1,n)).^2 + ...
                            (sigma_jit_pat*randn(1,n)).^2))./theta_div).^2);
Lp_t = patLoss(N_t).*patLoss(N_t);
P_rx_ideal_t = P_tx*eta_tx*eta_rx*(G_tx*G_rx)./FSPL_t;
P_rx_t = P_rx_ideal_t.*Lp_t;
fprintf('Received power at closest approach: %.2f dBm\n', 10*log10(P_rx_t(idx_cross)*1e3));

%% ========================================================================
% 4. QUANTUM STATE, SPECTRAL CHIRP, IMPRINTING AND OPERATING POINT
% =========================================================================
fprintf('\n=== 4. QUANTUM STATE AND OPERATING POINT ===\n');
cfg = struct();

% --- 4.1 Trap: the momentum spread follows from the motional ground state
cfg.f_secular = 3e6;
cfg.Delta = sqrt(m_Al27*hbar*(2*pi*cfg.f_secular)/2)/(m_Al27*c_light);   % units of p/(mc)

% --- 4.2 Superposition weight ------------------------------------------
% theta governs two quantities in opposite directions: the interference
% term read by the receiver scales as sin(2 theta) while the quantum
% correction to the classical Doppler shift, delta_Q of Eq. (32), scales as
% sin(4 theta). Here the chirp is the signal and delta_Q is a leakage
% channel, so the optimum is theta = pi/4, where the signal is maximal and
% delta_Q vanishes identically. This inverts the choice recommended in the
% reference, where delta_Q is itself the observable of interest.
cfg.theta = pi/4;

% --- 4.3 Differential preparation delay ---------------------------------
% The quadratic branch-to-branch phase requires the two branches to have
% evolved for different durations. Note that exp(-i p^2 t/2 m hbar) is the
% free-particle propagator: inside a harmonic trap a coherent state merely
% rotates in phase space and accumulates no such phase. Realising Delta_t
% therefore requires either released ions, in which case the branches
% separate by the amount reported below, or a state-dependent potential
% that holds one branch. This is an addition to the reference scheme and
% must be declared as such.
cfg.Delta_t = 20e-6;
cfg.C_chirp = mc2_Al27*cfg.Delta_t/(2*hbar*Omega_rad^2);   % s^2 rad^-2

% --- 4.4 Imprint depth --------------------------------------------------
% The ion imprints a fractional modulation on the probe beam. Its depth is
% bounded by the single-ion optical density, OD = sigma_abs/A_beam with
% sigma_abs = 3 lambda^2/(2 pi) at resonance. Treating the normalised
% classical profile as the full probe intensity would amount to assuming
% 100 % modulation by a single ion, which is not physical.
sigma_abs = 3*lambda_0^2/(2*pi);
w0_focus  = 0.20e-6;                       % beam waist at the ion (NA ~ 0.42)
OD_max    = sigma_abs/(pi*w0_focus^2);
cfg.OD    = 0.15;                          % design value, ~55 % of the bound
fprintf('Imprinting: sigma_abs = %.2e m^2, waist = %.2f um -> OD_max = %.3f | design OD = %.2f\n', ...
        sigma_abs, w0_focus*1e6, OD_max, cfg.OD);
if cfg.OD > OD_max
    warning('Design OD exceeds the single-ion bound (%.3f): reduce cfg.OD or tighten the focus.', OD_max);
end

% --- 4.5 Residual calibration errors ------------------------------------
cfg.eps_sub = 1e-3;   % per-bin relative error on the known classical profile;
                      % frozen over a frame, hence a bias rather than noise
cfg.eps_C   = 1e-3;   % relative error on the chirp coefficient, collecting
                      % both the Delta_t uncertainty and the uncompensated
                      % residual Doppler rate, which enters as a spurious
                      % quadratic phase adding to C
cfg.F_excess     = 3.0;    % 1.0 is the pure shot-noise limit
cfg.sigma_f_track = 1e3;   % residual Doppler offset error, 1-sigma, in Hz

cfg.N = 2^14;
SNR_thr_Eve = 3;           % assumed detection threshold for the eavesdropper

%% ------------------------------------------------------------------------
% 4.6 Operating-point optimisation
% -------------------------------------------------------------------------
% The packet separation is not chosen by hand. The covert volume is
%     B(sep) = sum_m n_m C(d_m),  truncated where
%     sqrt( sum_m n_m delta_m^2 ) exceeds the detection threshold,
% with d_m the matched-filter deflection at epoch m and delta_m the
% per-symbol deflection of the eavesdropper's detector. The maximum
% normally lies at the boundary where the eavesdropper reaches threshold
% exactly at the end of the pass: below it detection budget is wasted, above
% it the transmission is truncated.
fprintf('\n--- 4.6 Operating-point optimisation ---\n');

N_epochs = 21;
idx_ep   = round(linspace(1,N_t,N_epochs));
t_ep     = t_sec(idx_ep);
dt_ep    = (T_window/N_epochs)*ones(1,N_epochs);
avail_ep = double(link_available(idx_ep));

par = struct('Delta',cfg.Delta,'theta',cfg.theta,'Omega_rad',Omega_rad, ...
             'f0',f0_optical,'N',cfg.N,'C_chirp',cfg.C_chirp,'OD',cfg.OD, ...
             'F_excess',cfg.F_excess,'E_photon',E_photon,'eta_det',eta_det, ...
             'P_ep',P_rx_ideal_t(idx_ep),'dt_ep',dt_ep.*avail_ep, ...
             'SNR_thr',SNR_thr_Eve);

sep_vec = 4.0:0.05:9.5;
nSep = numel(sep_vec);
B_naive = zeros(1,nSep);    % covert volume against the energy detector
B_keyed = zeros(1,nSep);    % covert volume against the key-aware detector
B_noEve = zeros(1,nSep);    % volume with no eavesdropper present
for k = 1:nSep
    o = operatingPoint(sep_vec(k), par);
    if o.nyquistPhase < pi                  % reject an undersampled chirp
        B_naive(k) = o.volume_naive;
        B_keyed(k) = o.volume_keyed;
        B_noEve(k) = o.volume_noEve;
    end
end
[B_best, k_best] = max(B_naive);
cfg.sep_over_Delta = sep_vec(k_best);
opt = operatingPoint(cfg.sep_over_Delta, par);

fprintf('Optimal separation = %.2f Delta -> covert volume (energy detector) = %.3e bit\n', ...
        cfg.sep_over_Delta, B_best);
fprintf('  same point, no eavesdropper                = %.3e bit\n', opt.volume_noEve);
fprintf('  same point, key-aware eavesdropper         = %.3e bit\n', opt.volume_keyed);
fprintf('  median deflection d = %.3f -> raw BER = %.3f, capacity = %.4f bit/symbol\n', ...
        median(opt.d_ep), median(opt.BER_ep), median(opt.C_ep));

% Reference point: the hand-picked separation used before the optimisation
% was introduced, evaluated at the same imprint depth.
ref = operatingPoint(7.35, par);
fprintf('  reference sep = 7.35 (hand-picked): covert volume = %.3e bit, median d = %.3f\n', ...
        ref.volume_naive, median(ref.d_ep));

%% ------------------------------------------------------------------------
% 4.7 Frozen operating point and constraints on Delta_t
% -------------------------------------------------------------------------
cfg.p1     = -cfg.sep_over_Delta*cfg.Delta/2;
cfg.p2     = +cfg.sep_over_Delta*cfg.Delta/2;
cfg.B_IF   = opt.B_IF;
cfg.T_sym  = opt.T_sym;
cfg.dw     = opt.dw;
cfg.Phi    = opt.Phi;
cfg.key    = exp(+1i*cfg.C_chirp*cfg.dw.^2);

sig = struct('S_cl', opt.S_cl, 'I_env', opt.I_env, ...
             'tmpl', opt.I_env.*cos(opt.Phi));      % real matched-filter template

R_b_raw   = 1/cfg.T_sym;
dw_max    = 2*pi*cfg.B_IF/2;
dw_step   = cfg.dw(2)-cfg.dw(1);
lineWidth_Hz  = f0_optical*(cfg.sep_over_Delta+6)*cfg.Delta;
sigma_env_rad = (cfg.Delta/sqrt(2))*Omega_rad;

fprintf('\n--- 4.7 Frozen operating point ---\n');
fprintf('Trap f_sec = %.1f MHz -> Delta = %.3e p/(mc), velocity spread = %.4f m/s\n', ...
        cfg.f_secular/1e6, cfg.Delta, cfg.Delta*c_light);
fprintf('Packet separation = %.2f Delta -> velocity difference = %.4f m/s\n', ...
        cfg.sep_over_Delta, (cfg.p2-cfg.p1)*c_light);
fprintf('Derived B_IF = %.2f MHz -> T_sym = %.1f us -> raw bit rate = %.3f kbit/s\n', ...
        cfg.B_IF/1e6, cfg.T_sym*1e6, R_b_raw/1e3);
fprintf('Delta_t = %.1f us -> C = %.4e s^2/rad^2\n', cfg.Delta_t*1e6, cfg.C_chirp);
fprintf('Feasibility of Delta_t: in free flight the branches separate by %.1f um in %.1f us\n', ...
        (cfg.p2-cfg.p1)*c_light*cfg.Delta_t*1e6, cfg.Delta_t*1e6);

% Ceiling 1: Nyquist sampling of the chirp on the spectral grid.
phasePerSample = 2*cfg.C_chirp*(2*sigma_env_rad)*dw_step;
fprintf('Ceiling 1 (Nyquist): chirp phase per sample at the 2-sigma edge = %.3f rad -> Delta_t < %.1f us\n', ...
        phasePerSample, cfg.Delta_t*pi/phasePerSample*1e6);

% Ceiling 2: the chirp fringes must stay resolvable against the natural
% Lorentzian width of the transition.
fringeSpacing_Hz = 1/(2*cfg.C_chirp*(2*sigma_env_rad));
fprintf('Ceiling 2 (natural linewidth): fringe spacing at the edge = %.2f kHz vs Gamma0 = %d Hz -> Delta_t < %.0f us\n', ...
        fringeSpacing_Hz/1e3, Gamma0_Hz, cfg.Delta_t*fringeSpacing_Hz/Gamma0_Hz*1e6);
fprintf('Ceiling 3 (coherence time of the superposition) is experimental and outside this model.\n');

fprintf('Chirp time-bandwidth product BT = %.0f\n', cfg.C_chirp*(4*sigma_env_rad)^2/pi);
fprintf('Chirp-coefficient error: eps_C = %.1e -> residual edge phase = %.3f rad ', ...
        cfg.eps_C, cfg.eps_C*cfg.C_chirp*(2*sigma_env_rad)^2);
fprintf('(Delta_t must be known to %.3f %% for 1 rad)\n', 100/(cfg.C_chirp*(2*sigma_env_rad)^2));

%% ------------------------------------------------------------------------
% 4.8 Detection noise and Doppler tracking residual
% -------------------------------------------------------------------------
% The interference contrast is a fractional modulation, that is a ratio, and
% therefore does not attenuate with propagation. What changes with distance
% is the collected photon count, hence the relative shot noise:
%     sigma = F_excess / sqrt(n_bin),  n_bin = P_rx T_sym eta_det/(E_photon N)
% so that the signal-to-noise ratio scales as sqrt(P_rx) and not as P_rx.
fprintf('\n--- 4.8 Noise and tracking ---\n');
fprintf('Shot-derived noise (F_excess = %.1f): sigma = %.2e at closest approach, %.2e at pass edge\n', ...
        cfg.F_excess, min(opt.sigma_ep), max(opt.sigma_ep));

% A residual spectral offset s leaves a linear phase 2 C s dw after the
% de-chirp, displacing the compressed peak by 4 C df dw_max bins.
cfg.sigma_bin_track = 4*cfg.C_chirp*cfg.sigma_f_track*dw_max;
% With envelope weighting the compressed response is the transform of
% I_env^2, whose spectral width is sigma_env/sqrt(2); the peak is therefore
% sqrt(2) wider, and correspondingly more tolerant to tracking jitter, than
% with an unweighted sum.
sigma_k_peak = cfg.N/(2*pi*((sigma_env_rad/sqrt(2))/dw_step));
fprintf('Tracking offset: %.2f bins of jitter vs peak width sigma_k = %.1f bins -> %.2f dB loss\n', ...
        cfg.sigma_bin_track, sigma_k_peak, ...
        -20*log10(exp(-cfg.sigma_bin_track^2/(2*sigma_k_peak^2))));

%% ------------------------------------------------------------------------
% 4.9 Adversary models
% -------------------------------------------------------------------------
% Eve-A, key-unaware: knows the classical profile, subtracts it and sums
%   energies. Per-symbol deflection delta_A = d^2/sqrt(2N).
% Eve-B, key-aware (Kerckhoffs): knows C, which depends only on Delta_t, m
%   and Omega and is in any case estimable from the spectrum without
%   knowledge of the bit. Applies the same matched filter as the receiver
%   and measures |z|^2. Per-symbol deflection delta_B = d^2/sqrt(2).
% The ratio delta_B/delta_A is exactly sqrt(N). The receiver and Eve-B use
% the same filter; the only advantage the receiver holds is knowledge of
% phi, which reads the bit but is not needed to detect the transmission.
% Consequently, whenever the link decodes reliably, Eve-B detects it.
d_cross = interp1(t_ep, opt.d_ep, t_sec(idx_cross), 'linear', 'extrap');
delta_A = d_cross^2/sqrt(2*cfg.N);
delta_B = d_cross^2/sqrt(2);
fprintf('\n--- 4.9 Adversary models (evaluated at closest approach) ---\n');
fprintf('Receiver deflection d = %.3f -> raw BER = %.4f\n', d_cross, 0.5*erfc(d_cross/sqrt(2)));
fprintf('Eve-A (energy)   delta = %.3e per symbol -> threshold after %.2e symbols\n', ...
        delta_A, (SNR_thr_Eve/delta_A)^2);
fprintf('Eve-B (matched)  delta = %.3e per symbol -> threshold after %.2e symbols\n', ...
        delta_B, (SNR_thr_Eve/delta_B)^2);
fprintf('Sensitivity ratio Eve-B/Eve-A = %.1f (expected sqrt(N) = %.1f)\n', ...
        delta_B/delta_A, sqrt(cfg.N));

%% ------------------------------------------------------------------------
% 4.10 Leakage through the classical Doppler centroid
% -------------------------------------------------------------------------
% A non-zero delta_Q would shift the line centroid as a function of the bit,
% corrupting the classical velocimetry of the link and opening a third
% detection channel that exploits the entire photon count of the line. At
% theta = pi/4 it vanishes identically. Since theta cannot be calibrated
% exactly, the residual is evaluated both per symbol and accumulated over
% the whole pass, the latter being the binding requirement.
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
% delta_Q scales as sin(4 theta) ~ 4 eps near pi/4, so the tolerance scales
% linearly with the accumulated deflection.
eps_theta_req = eps_theta*SNR_thr_Eve/snr_centroid_acc;

fprintf('\n--- 4.10 Centroid leakage (delta_Q, Eq. 32) ---\n');
fprintf('At theta = pi/4 exactly: centroid swing = %.2e Hz (numerically zero)\n', swing_exact_Hz);
fprintf('With %.3f rad of miscalibration: swing = %.3f Hz vs centroid precision %.2f Hz\n', ...
        eps_theta, swing_eps_Hz, sigma_centroid_Hz);
fprintf('  per symbol SNR = %.3f, accumulated over %.2e symbols = %.1f (threshold %.1f)\n', ...
        snr_centroid_sym, opt.M_pass, snr_centroid_acc, SNR_thr_Eve);
fprintf('  REQUIREMENT: theta must be calibrated to better than %.1e rad (%.4f deg)\n', ...
        eps_theta_req, rad2deg(eps_theta_req));

%% ------------------------------------------------------------------------
% 4.11 Geometric low-probability-of-intercept argument
% -------------------------------------------------------------------------
% Because statistical covertness does not survive the Kerckhoffs assumption,
% the defensible security argument is that the eavesdropper must be
% physically inside the FSO beam. That claim is quantifiable.
fprintf('\n--- 4.11 Geometric LPI ---\n');
fprintf('Divergence %.1f urad -> beam radius %.2f m at %.0f km\n', ...
        theta_div*1e6, theta_div*d_min, d_min/1000);
fprintf('Illuminated fraction of solid angle = %.2e\n', theta_div^2/4);

%% ------------------------------------------------------------------------
% 4.12 Consistency checks
% -------------------------------------------------------------------------
S_bit0 = sig.S_cl + sig.I_env.*cos(cfg.Phi);
S_bit1 = sig.S_cl - sig.I_env.*cos(cfg.Phi);
residual_avg = max(abs(0.5*(S_bit0+S_bit1) - sig.S_cl))/max(abs(sig.S_cl));
fprintf('\n--- 4.12 Consistency checks ---\n');
fprintf('Cancellation of the quantum term under bit averaging: relative residual = %.3e\n', residual_avg);

weight = sig.I_env.*cfg.key;                % envelope weighting and de-chirp
k_expected = cfg.N/2 + 1;
Y_ref = fftshift(fft(ifftshift(sig.tmpl.*weight)));
[~, k_measured] = max(abs(real(Y_ref)));
fprintf('Compression bin: expected %d, measured %d\n', k_expected, k_measured);
if abs(k_measured-k_expected) > 2
    warning('Compression peak misplaced: check C_chirp and the spectral grid.');
end
cfg.kpk = k_measured;

% Matched-filter gain over the unweighted DC-bin statistic. The useful
% signal occupies only sqrt(2 pi) sigma_bin samples out of N, so an
% unweighted sum collects noise from all N samples for no additional signal.
g_mf   = sqrt(sum(sig.tmpl.^2));                          % ||s||
g_flat = sum(sig.tmpl.*real(cfg.key))/sqrt(cfg.N/2);      % unweighted statistic
fprintf('Matched-filter gain over the unweighted DC bin: %+.2f dB\n', ...
        20*log10(g_mf/max(abs(g_flat),eps)));
fprintf('Useful samples under the envelope: %.0f out of N = %d\n', ...
        sqrt(2*pi)*sigma_env_rad/dw_step, cfg.N);

cfg.Nbit_frame = 1024;
cfg.N_MC       = 16;
cfg.block      = 32;

%% ========================================================================
% 5. CAMPAIGN A: BIT-ERROR-RATE WATERFALL AND ANALYTICAL VALIDATION
% =========================================================================
% With a correct matched filter the simulated bit-error rate must reproduce
% Q(d) with d = ||s||/sigma, up to the tracking jitter and the calibration
% bias. A discrepancy indicates residual loss in the receiver.
fprintf('\n=== 5. CAMPAIGN A: WATERFALL ===\n');
sigma_vec = logspace(log10(min(opt.sigma_ep)/6), log10(max(opt.sigma_ep)*6), 16);
nS = numel(sigma_vec); errA = zeros(1,nS); bitA = zeros(1,nS);
BER_theory_A = 0.5*erfc((g_mf./sigma_vec)/sqrt(2));
for is = 1:nS
    for mc = 1:cfg.N_MC
        rs = mcStream(SEED_MASTER,1,is,mc);
        [ne,nb] = simulateFrame(cfg, sig, sigma_vec(is), rs);
        errA(is) = errA(is)+ne; bitA(is) = bitA(is)+nb;
    end
    fprintf('  sigma = %.2e | d = %5.2f | BER simulated = %.3e | BER theory = %.3e\n', ...
            sigma_vec(is), g_mf/sigma_vec(is), errA(is)/bitA(is), BER_theory_A(is));
end
BER_A = errA./bitA; [BER_A_lo,BER_A_hi] = wilsonCI(errA,bitA,1.96);
% The comparison is restricted to points with enough errors for the ratio to
% be statistically meaningful; below that the discrepancy is counting noise.
valid = errA >= 20;
if any(valid)
    fprintf('  maximum relative deviation over points with >= 20 errors: %.1f %%\n', ...
            100*max(abs(BER_A(valid)-BER_theory_A(valid))./BER_theory_A(valid)));
end

%% ========================================================================
% 6. CAMPAIGN B: PERFORMANCE ALONG THE ORBITAL PASS
% =========================================================================
% Each Monte Carlo realisation draws its own pointing loss, so a deep fade
% reduces the collected photon count and therefore raises sigma. This is the
% physical mechanism by which pointing degrades the error rate.
fprintf('\n=== 6. CAMPAIGN B: ORBITAL PASS ===\n');
BER_B = zeros(1,N_epochs); errB = zeros(1,N_epochs); bitB = zeros(1,N_epochs);
sigma_B = zeros(1,N_epochs);
for ie = 1:N_epochs
    ii = idx_ep(ie);
    acc = 0;
    for mc = 1:cfg.N_MC
        rs = mcStream(SEED_MASTER,2,ie,mc);
        th1 = sigma_jit_pat*randn(rs,1,2); th2 = sigma_jit_pat*randn(rs,1,2);
        Lp = exp(-8*(norm(th1)/theta_div)^2)*exp(-8*(norm(th2)/theta_div)^2);
        n_bin = (P_rx_ideal_t(ii)*Lp*cfg.T_sym/E_photon)*eta_det/cfg.N;
        sigma_eff = cfg.F_excess/sqrt(max(n_bin,eps));
        acc = acc + sigma_eff;
        [ne,nb] = simulateFrame(cfg, sig, sigma_eff, rs);
        errB(ie) = errB(ie)+ne; bitB(ie) = bitB(ie)+nb;
    end
    sigma_B(ie) = acc/cfg.N_MC;
    BER_B(ie)   = errB(ie)/bitB(ie);
    fprintf('  t = %6.1f s | range = %7.1f km | sigma = %.2e | d = %5.2f | BER = %.3e\n', ...
            t_ep(ie), d_t(ii)/1000, sigma_B(ie), g_mf/sigma_B(ie), BER_B(ie));
end
[BER_B_lo,BER_B_hi] = wilsonCI(errB,bitB,1.96);

%% ========================================================================
% 7. SENSITIVITY TO CALIBRATION AND IMPRINT DEPTH
% =========================================================================
% How accurately must the classical profile be known before subtraction?
% The error is frozen over a frame, hence it enters the decision statistic
% as a bias rather than as noise that averages away. As a reference, an
% estimate obtained by averaging M_cal symbols carries eps_sub ~
% sigma/sqrt(M_cal).
fprintf('\n=== 7. SENSITIVITY ANALYSES ===\n');
fprintf('--- 7.1 Subtraction error on the classical profile ---\n');
eps_vec = [0 1e-4 1e-3 1e-2 3e-2 1e-1];
sigma_ref = median(sigma_B);
BER_eps = zeros(1,numel(eps_vec));
cfg_tmp = cfg;
for k = 1:numel(eps_vec)
    cfg_tmp.eps_sub = eps_vec(k);
    e = 0; b = 0;
    for mc = 1:cfg.N_MC
        rs = mcStream(SEED_MASTER,3,k,mc);
        [ne,nb] = simulateFrame(cfg_tmp, sig, sigma_ref, rs);
        e = e+ne; b = b+nb;
    end
    BER_eps(k) = e/b;
    fprintf('  eps_sub = %.1e (equivalent to averaging %.1e calibration symbols) -> BER = %.3e\n', ...
            eps_vec(k), (sigma_ref/max(eps_vec(k),eps))^2, BER_eps(k));
end
fprintf('  This sweep models an independent per-bin error, which averages over the\n');
fprintf('  useful samples under the envelope. A structured error (trap drift, residual\n');
fprintf('  secular motion, local-oscillator drift) projects far more efficiently onto\n');
fprintf('  the template and requires a separate study.\n');

% The deflection is linear in OD, so this sensitivity needs no Monte Carlo.
fprintf('--- 7.2 Imprint depth ---\n');
OD_vec = sort([0.05 0.10 0.20 0.30 0.50 OD_max]);
for k = 1:numel(OD_vec)
    d_k = (OD_vec(k)/cfg.OD)*(g_mf/sigma_ref);
    flag = ''; if OD_vec(k) > OD_max, flag = '  (above the single-ion bound for this waist)'; end
    fprintf('  OD = %5.3f -> d = %6.3f -> raw BER = %.3e -> capacity = %.4f bit/symbol%s\n', ...
            OD_vec(k), d_k, 0.5*erfc(d_k/sqrt(2)), ...
            1-binaryEntropy(0.5*erfc(d_k/sqrt(2))), flag);
end

%% ========================================================================
% 8. CAPACITY, COVERT VOLUME AND SUMMARY
% =========================================================================
% The figure of merit is the per-symbol capacity and the volume integrated
% over the pass, not an uncoded error-rate threshold. A covert channel
% operates by construction at low signal-to-noise ratio, with reliability
% recovered by low-rate coding; the covert volume is what accumulates before
% the eavesdropper's detection statistic crosses threshold.
fprintf('\n=== 8. CAPACITY AND COVERT VOLUME ===\n');
n_sym_ep   = (dt_ep.*avail_ep)/cfg.T_sym;
C_ep_sim   = 1 - binaryEntropy(BER_B);
d_ep_sim   = g_mf./sigma_B;
delta_A_ep = d_ep_sim.^2/sqrt(2*cfg.N);
delta_B_ep = d_ep_sim.^2/sqrt(2);

bits_cum = cumsum(n_sym_ep.*C_ep_sim);
eveA_cum = sqrt(cumsum(n_sym_ep.*delta_A_ep.^2));
eveB_cum = sqrt(cumsum(n_sym_ep.*delta_B_ep.^2));

B_covert_A = volumeAtThreshold(bits_cum, eveA_cum, SNR_thr_Eve, n_sym_ep, C_ep_sim, delta_A_ep);
B_covert_B = volumeAtThreshold(bits_cum, eveB_cum, SNR_thr_Eve, n_sym_ep, C_ep_sim, delta_B_ep);
[~,ie_cr] = min(abs(idx_ep-idx_cross));

fprintf('\n============ SUMMARY: QUANTUM DOPPLER COVERT CHANNEL ============\n');
fprintf('Differential preparation delay Delta_t   : %8.1f us\n', cfg.Delta_t*1e6);
fprintf('Optimised packet separation              : %8.2f Delta\n', cfg.sep_over_Delta);
fprintf('Imprint depth OD                         : %8.3f  (single-ion bound %.3f)\n', cfg.OD, OD_max);
fprintf('Derived intermediate-frequency bandwidth : %8.2f MHz\n', cfg.B_IF/1e6);
fprintf('Raw symbol rate                          : %8.3f kbit/s\n', R_b_raw/1e3);
fprintf('Matched-filter gain over unweighted sum  : %+8.2f dB\n', 20*log10(g_mf/max(abs(g_flat),eps)));
fprintf('Raw BER at closest approach              : %8.3f\n', BER_B(ie_cr));
fprintf('Capacity at closest approach             : %8.4f bit/symbol\n', C_ep_sim(ie_cr));
fprintf('Median capacity over the pass            : %8.4f bit/symbol\n', median(C_ep_sim(avail_ep>0)));
fprintf('Peak information throughput              : %8.3f kbit/s\n', max(R_b_raw*C_ep_sim.*avail_ep)/1e3);
fprintf('----------------------------------------------------------------\n');
fprintf('VOLUME, no eavesdropper                  : %8.3e bit\n', bits_cum(end));
fprintf('COVERT VOLUME vs Eve-A (energy detector) : %8.3e bit\n', B_covert_A);
fprintf('  Eve-A statistic at end of pass         : %8.2f  (threshold %.1f)\n', eveA_cum(end), SNR_thr_Eve);
fprintf('COVERT VOLUME vs Eve-B (key-aware)       : %8.3e bit\n', B_covert_B);
fprintf('  Eve-B statistic at end of pass         : %8.2e (threshold %.1f)\n', eveB_cum(end), SNR_thr_Eve);
fprintf('----------------------------------------------------------------\n');
fprintf('The defensible figure is the Eve-A volume, and it holds only if covertness\n');
fprintf('rests on the geometric argument of Section 4.11. Under the Kerckhoffs\n');
fprintf('assumption (Eve-B) the channel is not covert.\n');
fprintf('================================================================\n');

%% ========================================================================
% 9. FIGURES
% =========================================================================
f_MHz = cfg.dw/(2*pi)/1e6;

% -------------------------------------------------------------------------
% FIG 1: SPECTRAL CHIRP AND MATCHED FILTERING
% -------------------------------------------------------------------------
figure('Name','Fig 1 Spectral chirp and matched filtering','Color','w','Position',[80 100 1200 380]);

% (a) Spectra
subplot(1,3,1);
plot(f_MHz, sig.S_cl, 'LineWidth', 1.4); hold on;
plot(f_MHz, S_bit0, 'LineWidth', 0.9);
xlabel('Detuning [MHz]'); 
ylabel('Modulation Depth [a.u.]');
title(sprintf('(a) Optical Spectrum (OD = %.2f)', cfg.OD));
legend({'S_{classical}','S_{quantum} (bit 0)'}, 'Location', 'best', 'FontSize', 8);
grid on; 
xlim([-1 1]*lineWidth_Hz/1e6);

% (b) Matched Filter Response
subplot(1,3,2);
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
xlim([-60 60]);

% (c) Waterfall
subplot(1,3,3);
bp = BER_A; 
bp(bp==0) = NaN;
errorbar(sigma_vec, bp, bp-BER_A_lo, BER_A_hi-bp, 'o', 'LineWidth', 1.2, 'MarkerSize', 4); hold on;
plot(sigma_vec, BER_theory_A, '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.3);
set(gca, 'XScale', 'log', 'YScale', 'log');
xline(min(opt.sigma_ep), 'g--', '\sigma closest', 'LineWidth', 1.1, 'LabelOrientation', 'horizontal');
xline(max(opt.sigma_ep), 'm--', '\sigma edge', 'LineWidth', 1.1, 'LabelOrientation', 'horizontal');
xlabel('Relative Noise Floor \sigma'); 
ylabel('Raw BER');
title('(c) Error Rate Validation');
legend({'Monte Carlo', 'theoretical filter curve'}, 'Location', 'best', 'FontSize', 8);
grid on;

% -------------------------------------------------------------------------
% FIG 2: ORBITAL PASS PERFORMANCE
% -------------------------------------------------------------------------
figure('Name','Fig 2 Orbital pass and eavesdropper accumulation','Color','w','Position',[110 80 900 380]);

% (a) BER and Capacity along pass
subplot(1,2,1);
yyaxis left;
h1 = errorbar(t_ep, BER_B, BER_B-BER_B_lo, BER_B_hi-BER_B, 'o-', 'Color', [0 0.45 0.74], 'LineWidth', 1.3, 'MarkerSize', 4);
ylabel('Raw BER'); 
ylim([0 0.55]);
ax = gca; ax.YColor = [0 0.45 0.74];

yyaxis right; 
h2 = plot(t_ep, C_ep_sim, 's--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2, 'MarkerSize', 4);
ylabel('Capacity [bit/symbol]'); 
ylim([0 1.05]);
ax.YColor = [0.85 0.33 0.1];

xline(t_sec(idx_cross), 'k:', 'Closest approach', 'LabelOrientation', 'horizontal'); 
xlabel('Time from Pass Start [s]');
title('(a) Channel Quality Over Fly-By Pass'); 
grid on;

% (b) Cumulative Volume
subplot(1,2,2);
plot(t_ep, bits_cum/1e3, 'LineWidth', 1.5); hold on;
yline(B_covert_A/1e3, 'g--', 'Covert volume vs Eve-A', 'LineWidth', 1.2);
yline(B_covert_B/1e3, 'r--', 'Covert volume vs Eve-B', 'LineWidth', 1.2);
xlabel('Time from Pass Start [s]'); 
ylabel('Cumulative Volume [kbit]');
title('(b) Covert Information Volume'); 
grid on;

% -------------------------------------------------------------------------
% FIG 3: DESIGN OPTIMIZATION AND SENSITIVITY
% -------------------------------------------------------------------------
figure('Name','Fig 3 Design optimisation and sensitivities','Color','w','Position',[140 60 900 420]);

% (a) Operating Point Sweep
subplot(1,2,1);
semilogy(sep_vec, max(B_noEve,1e-3), '--', 'LineWidth', 1.2); hold on;
semilogy(sep_vec, max(B_naive,1e-3), 'LineWidth', 1.6);
semilogy(sep_vec, max(B_keyed,1e-3), 'LineWidth', 1.2);
xline(cfg.sep_over_Delta, 'k:', 'optimum', 'LineWidth', 1.2);
xlabel('Wavepacket Separation [\Delta_p]'); 
ylabel('Pass Volume [bit]');
title('(a) Operating Point Optimization');
legend({'no eavesdropper', 'covert vs Eve-A', 'covert vs Eve-B'}, ...
       'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 7.5);
grid on;

% (b) Imprint Depth Sensitivity
subplot(1,2,2);
OD_fine = linspace(0.01, max(OD_max,0.6), 200);
d_fine  = (OD_fine/cfg.OD)*(g_mf/sigma_ref);
plot(OD_fine, 1-binaryEntropy(0.5*erfc(d_fine/sqrt(2))), 'LineWidth', 1.5); hold on;
xline(cfg.OD, 'k:', 'design OD', 'LineWidth', 1.2);
xline(OD_max, 'r:', 'single-ion limit', 'LineWidth', 1.2);
xlabel('Imprint Depth [OD]'); 
ylabel('Capacity [bit/symbol]');
title('(b) Imprint Depth Sensitivity'); 
grid on;

fprintf('\n--- SIMULATION COMPLETE ---\n');

%% ========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function out = operatingPoint(sep, par)
% Design point for a given packet separation. Returns the spectra, the
% matched-filter deflection of the receiver, the per-symbol deflections of
% both adversary models, and the resulting information volumes. Every
% quantity is derived from the physical parameters in par.
    Delta = par.Delta; N = par.N;
    B_IF  = 4*par.f0*(sep+6)*Delta;         % from the atomic spectral width
    T_sym = N/B_IF;
    dw = linspace(-1,1,N).'*(2*pi*B_IF/2);
    p  = dw/par.Omega_rad;                  % momentum in units of p/(mc)
    A1 = exp(-((p+sep*Delta/2).^2)./(2*Delta^2))./(pi^0.25*sqrt(Delta));
    A2 = exp(-((p-sep*Delta/2).^2)./(2*Delta^2))./(pi^0.25*sqrt(Delta));
    w  = 1 + 3*p;                           % emission weight of Eq. (26)
    S_cl  = (cos(par.theta)^2*A1.^2 + sin(par.theta)^2*A2.^2).*w;
    I_env = (2*cos(par.theta)*sin(par.theta)*A1.*A2).*w;
    % Normalise the classical profile to unit root-mean-square, then scale
    % both terms by the fractional imprint depth so that sigma below is
    % directly the relative shot noise.
    nrm   = sqrt(mean(S_cl.^2));
    S_cl  = par.OD*S_cl/nrm;
    I_env = par.OD*I_env/nrm;
    Phi   = -par.C_chirp*dw.^2;

    signalEnergy = sum((I_env.*cos(Phi)).^2);          % ||s||^2 at unit sigma

    n_sym    = par.dt_ep(:).'/T_sym;
    n_bin_ep = (par.P_ep(:).'*T_sym/par.E_photon)*par.eta_det/N;
    sigma_ep = par.F_excess./sqrt(max(n_bin_ep,eps));
    d_ep     = sqrt(signalEnergy)./sigma_ep;
    BER_ep   = 0.5*erfc(d_ep/sqrt(2));
    C_ep     = 1 - binaryEntropy(BER_ep);
    delta_A  = d_ep.^2/sqrt(2*N);                      % key-unaware detector
    delta_B  = d_ep.^2/sqrt(2);                        % key-aware detector

    bits_cum = cumsum(n_sym.*C_ep);
    eveA     = sqrt(cumsum(n_sym.*delta_A.^2));
    eveB     = sqrt(cumsum(n_sym.*delta_B.^2));

    out.B_IF = B_IF;  out.T_sym = T_sym;  out.dw = dw;
    out.S_cl = S_cl;  out.I_env = I_env;  out.Phi = Phi;
    out.sigma_ep = sigma_ep;  out.d_ep = d_ep;
    out.BER_ep = BER_ep;      out.C_ep = C_ep;
    out.M_pass = sum(n_sym);
    out.volume_noEve = bits_cum(end);
    out.volume_naive = volumeAtThreshold(bits_cum, eveA, par.SNR_thr, n_sym, C_ep, delta_A);
    out.volume_keyed = volumeAtThreshold(bits_cum, eveB, par.SNR_thr, n_sym, C_ep, delta_B);
    out.nyquistPhase = 2*par.C_chirp*(2*(Delta/sqrt(2))*par.Omega_rad)*(dw(2)-dw(1));
end

function B = volumeAtThreshold(bits_cum, eve, thr, n_sym, C_ep, delta)
% Information volume transmissible before the accumulated detection
% statistic of the eavesdropper crosses threshold. Within the crossing epoch
% the statistic is interpolated linearly in accumulated deflection energy.
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

function [nerr,nbit] = simulateFrame(cfg, sig, sigma_noise, rs)
% Transmission and reception of one frame. The bit is the relative phase phi
% in {0, pi} of the momentum superposition; by Eq. (32) the interference
% term varies as cos(phi), so the signalling is antipodal by construction.
%
% The transmitted contrast is a fractional modulation and is therefore not
% scaled with range; distance enters only through sigma_noise, which the
% caller derives from the collected photon count.
%
% The receiver weights by the known interference envelope before applying
% the de-chirp key. The direct-current bin of fft(Rx .* I_env .* key) is
% then exactly sum(Rx .* I_env .* cos(Phi)), the matched filter for the real
% template.
    M = cfg.Nbit_frame; N = cfg.N;
    bits = randi(rs,[0 1],1,M);

    % Calibration error on the classical profile, frozen over the frame and
    % therefore acting as a decision bias rather than as averaging noise.
    S_hat = sig.S_cl.*(1 + cfg.eps_sub*randn(rs,N,1));

    % Residual error on the chirp coefficient, collecting the Delta_t
    % uncertainty and the uncompensated Doppler rate.
    C_bob  = cfg.C_chirp*(1 + cfg.eps_C*randn(rs,1,1));
    weight = sig.I_env.*exp(+1i*C_bob*cfg.dw.^2);

    nerr = 0;
    for i1 = 1:cfg.block:M
        i2 = min(i1+cfg.block-1,M); b = bits(i1:i2); nb = numel(b);
        s  = 1-2*b;                                   % bit 0 -> +1, bit 1 -> -1
        Rx = sig.S_cl + sig.tmpl.*s + sigma_noise*randn(rs,N,nb);
        Rx = Rx - S_hat;
        % The ifftshift is required: the signal is centred at zero detuning,
        % that is at the centre of the array rather than at the first index.
        % Without it the centring is read as a phase ramp and the real part
        % of the compressed peak alternates sign from bin to bin.
        Y  = fftshift(fft(ifftshift(Rx.*weight,1),[],1),1);
        % Residual Doppler offset displaces the read bin by 4 C df dw_max.
        jit  = round(cfg.sigma_bin_track*randn(rs,1,nb));
        kuse = min(max(cfg.kpk + jit,1),N);
        stat = real(Y(sub2ind(size(Y),kuse,1:nb)));
        nerr = nerr + sum(double(stat<0) ~= b);
    end
    nbit = M;
end

function rs = mcStream(seed,campaign,point,trial)
% Independent, reproducible substreams: campaign, sweep point and trial
% index map to disjoint regions of a single mrg32k3a stream.
    rs = RandStream('mrg32k3a','Seed',seed);
    rs.Substream = 1 + trial + 1000*point + 1000000*campaign;
end

function [lo,hi] = wilsonCI(k,n,z)
% Wilson score interval for a binomial proportion.
    p = k./n; d = 1+z^2./n; c = p+z^2./(2*n);
    s = z*sqrt(p.*(1-p)./n + z^2./(4*n.^2));
    lo = max((c-s)./d,0); hi = min((c+s)./d,1);
end

function H = binaryEntropy(p)
    p = min(max(p,0),1); H = zeros(size(p)); m = (p>0)&(p<1);
    H(m) = -p(m).*log2(p(m)) - (1-p(m)).*log2(1-p(m));
end

% =========================================================================
% NOTES FOR THE ACCOMPANYING DOCUMENT
%
% The following must be stated explicitly in any write-up of this work. They
% are ordered from the claims the simulator supports, through the modelling
% choices that shape the numbers, to the limitations that bound them.
% =========================================================================
%
% NOTE 1 - WHAT IS GENUINELY QUANTUM HERE, AND WHAT IT BUYS.
%   The carrier is the cross term of |psi_sup|^2, that is the difference
%   between a coherent superposition [Eq. A1] and a statistical mixture
%   [Eq. A2] of the same two wave packets. It vanishes identically without
%   momentum coherence, and the simulator verifies at run time that
%   averaging over the two bits returns exactly the classical profile.
%   Two honest qualifications must accompany that statement. First, the
%   observable signature is a spectral fringe pattern, which is trivially
%   reproducible by classical means: an electro-optic modulator would
%   generate the same chirp with a contrast many orders of magnitude larger,
%   at higher bandwidth and with no cryogenic ion trap. From a purely
%   communications standpoint the quantum mechanism confers no advantage.
%   Second, the information is carried by phi, the relative phase of the
%   preparation laser; the quantum state is the transducer, not the message.
%   The value of the scheme is therefore as a physical-signature channel and
%   as a quantitative feasibility study, not as a competitive modulation.
%
% NOTE 2 - THE CHIRP IS AN EXTENSION OF THE REFERENCE, NOT A RESULT OF IT.
%   With the real Gaussian packets of Eq. (A1) one has Phi_1 = Phi_2 = 0 and
%   the interference term is a single lobe modulated by cos(phi), not a
%   chirp. Free dispersion does not help: exp(-i p^2 t/2 m hbar) is common
%   to both branches and cancels in the squared modulus, since for a free
%   particle |psi(p)|^2 is rigorously constant in time. A relative quadratic
%   phase requires the branches to evolve for different durations. Preparing
%   them Delta_t apart yields Phi_2 - Phi_1 = -p^2 Delta_t/(2 m hbar) and,
%   through the Doppler map, Phi(dw) = -C dw^2 with C = m c^2 Delta_t /
%   (2 hbar Omega^2). Grochowski et al. prepare the branches simultaneously;
%   this delay is an addition to their experimental scheme and must be
%   declared as such. The chirp coefficient nevertheless remains fully
%   determined by Delta_t and atomic constants.
%
% NOTE 3 - DERIVED VERSUS ASSUMED PARAMETERS.
%   Derived, and therefore not adjustable:
%     Delta, the momentum spread, from the motional ground state of the trap;
%     C, from Delta_t and atomic constants;
%     B_IF, from the width of the atomic spectral structure;
%     the antipodal character of the signalling, from cos(phi) in Eq. (32);
%     the compression bin, which the de-chirp places at the grid centre;
%     the detection noise, from the collected photon count;
%     the packet separation, from the constrained optimisation of Section 4.6;
%     the interference contrast, from the packet overlap.
%   Assumed, and to be declared:
%     F_excess, the excess-noise factor over the shot limit;
%     OD, the imprint depth, bounded above by the single-ion optical density;
%     eps_sub and eps_C, the residual calibration errors;
%     sigma_f_track, the residual Doppler offset;
%     Delta_t, chosen within the three ceilings of Note 6;
%     the optical hardware parameters (apertures, power, efficiencies);
%     the orbital scenario (altitude, inclination, fly-by geometry);
%     the eavesdropper's detection threshold, taken as 3.
%
% NOTE 4 - THE RECEIVER MUST BE A MATCHED FILTER, AND WHY.
%   After the de-chirp the useful signal is (1/2) I_env(w) exp(i phi), which
%   occupies of order sqrt(2 pi) sigma_bin samples out of N. The direct-
%   current bin of an unweighted transform sums all N samples and therefore
%   collects noise sigma sqrt(N/2) against a signal proportional to
%   sum(I_env). The ratio of the two deflections is
%       SNR_matched / SNR_unweighted = ||g|| sqrt(N/2) / sum(g),
%   which for a Gaussian envelope of sigma_bin samples is approximately
%   sqrt(N / (2 sqrt(pi) sigma_bin)). At the operating point of this study
%   the simulator measures roughly 13 dB. An earlier version of this work
%   read the unweighted bin and consequently reported a bit-error rate near
%   0.5 across the entire pass and zero goodput; that result was an artefact
%   of the receiver, not a property of the channel. Campaign A now validates
%   the receiver by reproducing Q(||s||/sigma) analytically.
%
% NOTE 5 - THE PERFORMANCE CRITERION.
%   A covert channel operates by construction at low signal-to-noise ratio.
%   Its optimal operating point lies at a raw bit-error rate of order 0.2 to
%   0.4, with reliability recovered by low-rate forward error correction.
%   Evaluating uncoded frame success, (1-BER)^L, in that regime is
%   meaningless and returns identically zero. The figure of merit adopted
%   here is the per-symbol capacity 1 - H2(BER), integrated over the symbols
%   of the pass.
%   A statistical caveat belongs alongside it. The Wilson half-width on a
%   bit-error rate near one half is approximately 1.96 sqrt(0.25/n). If the
%   operating point drifts to a raw error rate near 0.49, the corresponding
%   capacity is not distinguishable from zero at any practical n, and must
%   then be reported as an upper bound rather than as a measurement. The
%   frame length and repetition count used here place the capacity estimate
%   several standard deviations away from zero; that must be checked, not
%   assumed, whenever parameters change.
%
% NOTE 6 - CONSTRAINTS ON Delta_t, AND THE FEASIBILITY QUESTION.
%   Three ceilings apply, two of which the simulator evaluates:
%     (a) Nyquist sampling of the chirp on the spectral grid;
%     (b) the chirp fringe spacing at the envelope edge must remain above
%         the natural linewidth Gamma0, otherwise the fringes are washed out
%         by the Lorentzian of the transition;
%     (c) the coherence time of the momentum superposition, which is
%         experimental and cannot be verified from within this model.
%   The compression gain grows with Delta_t while the tolerance to residual
%   Doppler error degrades in the same proportion, since both scale with C.
%   That trade-off is a central design result and should be presented as
%   such.
%   A separate and more serious question concerns the mechanism itself. The
%   propagator exp(-i p^2 t/2 m hbar) is that of a free particle. Inside a
%   harmonic trap a coherent state rotates in phase space and accumulates no
%   quadratic phase in p, so Delta_t cannot be realised by simply waiting.
%   Either the ion is released, in which case the simulator reports the
%   resulting spatial excursion of the two branches, or a state-dependent
%   potential holds one branch while the other evolves freely. The second
%   option is the honest one to propose, and its coherence budget must be
%   argued from the cat-state literature rather than assumed.
%
% NOTE 7 - PHOTON BUDGET AND THE IMPRINTING ARCHITECTURE.
%   Direct fluorescence from a single ion, collected across hundreds of
%   kilometres by a 0.30 m aperture, is far below one photon per second and
%   cannot close the link. The model therefore assumes that the spectral
%   response is imprinted on a transmitted probe beam by absorption and
%   dispersion, consistent with the emission-to-absorption transformation
%   indicated in the reference; the photon budget crossing the link is that
%   of the laser. The depth of that imprint is bounded by the single-ion
%   optical density, sigma_abs/A_beam with sigma_abs = 3 lambda^2/(2 pi),
%   and the receiver deflection is linear in it. After Delta_t, the quality
%   of the focus on the ion is the most critical experimental parameter.
%
% NOTE 8 - THE COVERTNESS RESULT, STATED WITHOUT EMBELLISHMENT.
%   Two adversary models are evaluated. The key-unaware detector sums
%   residual energy after subtracting the known classical profile; the
%   key-aware detector applies the same matched filter as the receiver. The
%   ratio of their per-symbol deflections is exactly sqrt(N), verified at
%   run time. The receiver's only advantage over the key-aware detector is
%   knowledge of phi, which reads the bit but is not required to detect the
%   transmission. It follows that whenever the link decodes reliably, a
%   key-aware eavesdropper detects it. The covertness of this scheme cannot
%   rest on the eavesdropper's ignorance of C: that coefficient depends only
%   on Delta_t, m and Omega, and is in any case estimable from the received
%   spectrum by autocorrelation or Wigner analysis, over a one-dimensional
%   search space.
%   The defensible security argument is geometric. With a divergence of a
%   few microradians the illuminated fraction of solid angle is of order
%   1e-11: the eavesdropper must be physically inside the beam. That should
%   be presented as the primary mechanism, with the key-unaware statistic as
%   a second layer against an interceptor who is inside the beam but does
%   not know what to look for.
%   The choice theta = pi/4 closes a third channel by nulling delta_Q, the
%   quantum correction to the classical Doppler shift, which would otherwise
%   move the line centroid with the bit and expose the transmission to
%   centroid tracking. That null is exact only at exactly pi/4; the residual
%   accumulates over the pass like any other detection statistic, and
%   Section 4.10 converts the accumulated constraint into a calibration
%   tolerance on theta. That tolerance is tight and belongs in the
%   experimental requirements.
%
% NOTE 9 - RESULTS, AND HOW TO READ THEM.
%   The headline quantity is the covert volume per pass against the
%   key-unaware detector, printed in the summary block and plotted as the
%   cumulative curve of Fig. 2(b). At the operating point studied here it is
%   of order 1e5 bits, and it coincides with the volume that would be
%   transmissible with no eavesdropper present: the binding constraint is
%   the duration of the pass, not detection. Against the key-aware detector
%   the same configuration yields of order 1e3 bits.
%   Two structural observations deserve emphasis. First, repeating the
%   optimisation at different imprint depths leaves the covert volume nearly
%   unchanged, because the optimiser trades deflection against symbol count;
%   this is the square-root law of covert communication appearing as an
%   invariance rather than as an assumption. The analytical optimum,
%   d_opt = (2 N thr^2 / M_pass)^(1/4), agrees with the numerical sweep.
%   Second, the design sweep uses the unfaded received power while the Monte
%   Carlo campaign applies pointing loss, so the simulated volume falls
%   below the design prediction by roughly the mean fading factor. The
%   Monte Carlo figure is the one to quote. Aligning the two, by feeding the
%   fading-averaged power into the design function, would shift the optimum
%   slightly toward smaller separations and is a recommended refinement.
%   On runtime: the cost is dominated by Gaussian variate generation, not by
%   the transforms. It can be reduced by enlarging the block size, by
%   reducing N given that the envelope is oversampled on the present grid,
%   by using common random numbers across the calibration sweep, and by
%   sampling the sufficient statistic directly away from the waterfall knee
%   while retaining the full simulation for validation. None of these
%   changes the physics; all of them change the random-number consumption
%   order and therefore the exact realisations.
%
% NOTE 10 - KNOWN LIMITATIONS OF THE MODEL.
%   The image term is intrinsic and not an implementation defect. The
%   transported quantity is real, so cos(Phi + phi) contains both exp(i(Phi
%   + phi)) and its conjugate; the de-chirp compresses only the first, and
%   half the energy remains chirped and spread over the time-bandwidth
%   product. This costs 3 dB and is already included in the deflection.
%   The subtraction-error sweep models an independent per-bin error, which
%   averages efficiently over the samples under the envelope and therefore
%   shows a wide margin. A structured error, from trap drift, residual
%   secular motion or local-oscillator drift, projects far more efficiently
%   onto the template and is not covered. It should be studied separately
%   before the calibration margin reported here is relied upon. The sweep is
%   also evaluated at the pass-median noise level, where the deflection is
%   already small; repeating it at the closest-approach noise level would
%   make the sensitivity visible.
%   The de-chirp implemented here is the chirp-multiply and transform branch
%   of the fractional Fourier transform, not a full fractional transform.
%   For peak-polarity detection the two are equivalent, but the distinction
%   should be stated rather than glossed, and a demonstration of equivalence
%   with a full implementation would strengthen the presentation.
%   Finally, the residual Doppler model separates offset from rate: the
%   offset displaces the read bin and is modelled explicitly, while the
%   uncompensated rate is folded into eps_C as a spurious quadratic phase.
%   Given a Doppler rate of hundreds of megahertz per second at closest
%   approach, the rate specification deserves an independent budget rather
%   than absorption into a single relative coefficient.
% =========================================================================