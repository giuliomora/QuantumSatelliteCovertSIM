% =========================================================================
% QUANTUM-COVERT INTER-SATELLITE LINK SIMULATOR
% Spectral fringe encoding on the quantum Doppler effect, matched-filter
% reception, and covert-capacity assessment over a LEO-LEO pass.
%
% PHYSICAL BASIS
%   P. T. Grochowski, A. R. H. Smith, A. Dragan, K. Debski, "Quantum time
%   dilation in atomic spectra", Phys. Rev. Research 3, 023053 (2021).
%   Equation numbers quoted in the code refer to that work.
%
%   WHAT IS TAKEN FROM THE REFERENCE, AND WHAT IS ADDED HERE.
%   Taken: the coherent superposition of two momentum wave packets [Eq. A1],
%   the incoherent mixture it must be distinguished from [Eq. A2], the
%   first-order Doppler line shape [Eq. 26], and the centroid correction
%   delta_Q [Eq. 32]. Added: a relative spatial displacement between the
%   branches, which the reference does not consider and which turns the
%   smooth interference lobe into a periodic fringe. That fringe is what
%   makes the signature readable with a shared key, and it is the one
%   element of the physics that is NOT in the reference. It has been checked
%   numerically - free evolution alone produces no fringe, a relative
%   displacement produces one and it persists indefinitely, and the bit
%   inverts it exactly - but it is an extension and must be presented as one.
%
% PERIMETER
%   This is a CHANNEL simulator. It takes an emitter that writes the quantum
%   Doppler signature on a probe beam, described by three declared numbers,
%   and asks what covert channel can be built on it. It does not design the
%   emitter, and it does not choose those three numbers to flatter the
%   result: each is swept for sensitivity. See Note 3.
%
% DECLARED INPUTS DESCRIBING THE SOURCE (Section 4)
%   Delta    velocity spread of the two momentum packets, in units of v/c
%   Delta_x  spatial separation between the branches
%   OD       root-mean-square fractional modulation written on the probe
%
% MODEL CHAIN (everything below follows from the declared inputs)
%   1. The emitter is prepared in a coherent superposition of two momentum
%      wave packets [Eq. A1] with relative phase phi in {0, pi}: the bit.
%   2. The branches are separated in space by Delta_x. A displacement is a
%      phase linear in p, exp(i p Delta_x/hbar), and it is RELATIVE between
%      the branches, so it survives in |psi(p)|^2. Anything common to both
%      branches - free dispersion, or a spectral phase applied downstream to
%      the probe - cancels in the cross term and produces no fringe.
%   3. The first-order Doppler map dw = Omega p/(mc) turns that into
%      Phi(dw) = kappa dw with kappa = m c Delta_x/(hbar Omega): a fringe of
%      uniform spacing 1/kappa hertz.
%   4. The spectrum imprinted on the probe is
%         S_quant(w) = S_class(w) + I_env(w) cos(Phi(w) + phi)
%      where S_class is the incoherent statistical mixture [Eq. A2] and
%      I_env is the interference envelope.
%   5. The probe crosses the link. The classical Doppler shift of tens of
%      GHz is assumed compensated from the ephemerides; only the residual
%      offset is modelled.
%   6. The receiver subtracts the known classical profile, weights by the
%      envelope, multiplies by the fringe key exp(+i kappa dw), and reads
%      the sign of the compressed peak. The transform is an ordinary FFT:
%      with a linear fringe phase no fractional-domain step is required.
%   7. Averaged over equiprobable bits the interference term cancels,
%      cos(x) + cos(x+pi) = 0, and only S_class survives. This is checked
%      at run time in Section 4.12.
%
% REQUIRED TOOLBOXES
%   Satellite Communications Toolbox (satelliteScenario, satellite, states)
%   Phased Array System Toolbox (fspl)
%
% RUNTIME
%   The Monte Carlo campaigns run on the sufficient statistic (cfg.fastMC),
%   which replaces N Gaussian variates per symbol with one, exactly and not
%   approximately. Expect well under a minute. Setting cfg.fastMC = false
%   restores the full spectral path, which costs of order 1e10 variates and
%   several minutes; it should be run at least once whenever the signal model
%   changes. See Note 9.
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
% NOTE ON THE CLASSICAL DOPPLER SHIFT. The 26.7 GHz shift is computed and
% reported here, but it is never applied to the spectrum and then removed.
% The model ASSUMES it is compensated from the ephemerides, which is a
% declared assumption and not a result of the simulation, and models only
% what compensation leaves behind: a residual offset (Section 4.8) and a
% residual rate folded into eps_kappa (Section 4.5). The text must say this
% explicitly rather than implying that the simulation demonstrates the
% compensation.
%
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
% The analytic model below, including the J2 nodal drift, is used ONLY to
% choose the initial elements. Once those are fixed, SGP4 propagates them and
% is the sole authority for the geometry used everywhere downstream.
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

% Solar exclusion: reject epochs with the Sun within 30 deg of the receive
% line of sight. The solar position uses the mean longitude only, without the
% equation of centre, so it carries an error of up to about 2 degrees. That is
% immaterial against a 30 degree exclusion cone but would not be for a
% tighter one.
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

% --- Verification of the orbital block -----------------------------------
% Independent cross-checks. Each recomputes a printed quantity by a different
% route; a mismatch means the geometry is not what the model believes.
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
% The probe beam that carries the imprinted spectral response is a
% conventional FSO link. Pointing, acquisition and tracking jitter is
% applied independently at both terminals.
%
% BEAM DEFINITION, AND WHY THE APERTURE FORMULA IS NOT USED.
% Two descriptions of a transmitted beam are available and they describe
% different beams. The aperture gain rho (pi D/lambda)^2 belongs to a beam
% diffraction-limited by the 0.20 m aperture, whose far-field divergence is
% 1.7 urad. With 1.5 urad of pointing jitter the error would be comparable
% to the beam width and the link would not close. The beam is therefore
% deliberately EXPANDED to theta_div = 8 urad, and the transmit gain must be
% the one belonging to that beam and not to the aperture. The cost of the
% expansion is about 8.6 dB and it is paid deliberately.
% CONVENTION: theta_div is the FULL 1/e^2 divergence angle, so the half-angle
% is theta_div/2 and G = 8/(theta_div/2)^2. The text must state this: the
% alternative convention moves the budget by about 6 dB.
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

% --- Verification of the link budget -------------------------------------
% Term-by-term recomputation, so the budget can be checked against a hand
% calculation rather than trusted as a single number.
fprintf('  budget at closest approach [dB]: Ptx %.1f + Gtx %.1f + Grx %.1f - FSPL %.1f + eta %.1f = %.2f dBm\n', ...
        10*log10(P_tx*1e3), 10*log10(G_tx), 10*log10(G_rx), ...
        10*log10(FSPL_t(idx_cross)), 10*log10(eta_tx*eta_rx), ...
        10*log10(P_rx_ideal_t(idx_cross)*1e3));
fprintf('  pointing loss at that epoch: %.2f dB (mean over the pass %.2f dB)\n', ...
        -10*log10(Lp_t(idx_cross)), -10*log10(mean(Lp_t)));

% Cost of the beam expansion, reported so it is visible rather than implied.
% G_ap is the gain the aperture would give if the beam were diffraction
% limited; the difference is what is paid to make the pointing tolerable.
G_ap        = rho_ap*(pi*D_tx/lambda_0)^2;
theta_diff  = 2*lambda_0/(pi*D_tx/2);
fprintf('  BEAM: expanded to %.1f urad (gain %.1f dB) against a diffraction limit of %.2f urad (gain %.1f dB)\n', ...
        theta_div*1e6, 10*log10(G_tx), theta_diff*1e6, 10*log10(G_ap));
fprintf('  BEAM: cost of the expansion %.1f dB, bought against a pointing jitter of %.2f urad\n', ...
        10*log10(G_ap/G_tx), sigma_jit_pat*1e6);
if theta_div < 3*sigma_jit_pat
    warning('Beam width is within three times the pointing jitter: the loss model is not reliable.');
end

% Pointing statistics. The jitter is a small fraction of the beam width, but
% the loss distribution has a long tail: the fifth percentile is what sets
% the worst usable epochs, not the mean.
Lp_sorted = sort(Lp_t);
fprintf('  pointing loss over the pass: mean %.2f dB, median %.2f dB, 5th percentile %.2f dB\n', ...
        -10*log10(mean(Lp_t)), -10*log10(median(Lp_t)), ...
        -10*log10(Lp_sorted(max(1,round(0.05*N_t)))));
% Modelling note: the same theta_div is used for the loss at both terminals.
% At the receiver the relevant quantity is the field of view rather than the
% transmit beam width, so this is a simplification and, with equal jitter at
% both ends, a conservative one.

%% ========================================================================
% 4. QUANTUM STATE, SPECTRAL CHIRP, IMPRINTING AND OPERATING POINT
% =========================================================================
fprintf('\n=== 4. QUANTUM STATE AND OPERATING POINT ===\n');
cfg = struct();

% --- 4.1 Emitter: momentum spread of the source -------------------------
% The two branches are Gaussian momentum packets of spread Delta, in units
% of v/c. Delta is DECLARED. The expression below is only a convenient way to
% land on a plausible value - it is the ground-state spread of a source
% confined at f_source - and must not be read as a derivation: nothing in
% this model requires the source to be confined, and a free ensemble with the
% same Delta behaves identically. What an atomic source can deliver is a
% question of apparatus and is an input here, not a subject.
cfg.f_source = 1e6;
cfg.Delta = sqrt(hbar*(2*pi*cfg.f_source)/(2*m_Al27))/c_light;   % v/c

% --- 4.2 Superposition weight ------------------------------------------
% theta governs two quantities in opposite directions: the interference
% term read by the receiver scales as sin(2 theta) while the quantum
% correction to the classical Doppler shift, delta_Q of Eq. (32), scales as
% sin(4 theta). Here the fringe is the signal and delta_Q is a leakage
% channel, so the optimum is theta = pi/4, where the signal is maximal and
% delta_Q vanishes identically. This inverts the choice recommended in the
% reference, where delta_Q is itself the observable of interest.
cfg.theta = pi/4;

% --- 4.3 Branch separation and the spectral fringe ----------------------
% Only phases RELATIVE between the branches survive in |psi(p)|^2. Anything
% common to both - free dispersion, or a spectral phase applied downstream
% to the probe - is a multiplicative factor that cancels in the cross term
% and produces no fringe. The relative phase that does survive is a spatial
% displacement: branches separated by Delta_x carry exp(i p Delta_x/hbar),
% which through the Doppler map dw = Omega p/(mc) becomes
%     Phi(dw) = kappa dw ,   kappa = m c Delta_x/(hbar Omega) ,
% a fringe of uniform spacing 1/kappa hertz whose sign is flipped by phi.
% Delta_x is declared for the same reason as Delta.
cfg.Delta_x = 90e-9;
cfg.kappa   = m_Al27*c_light*cfg.Delta_x/(hbar*Omega_rad);

% --- 4.4 Imprint depth --------------------------------------------------
% OD is the root-mean-square fractional modulation the emitter writes on the
% probe beam. Declared, not designed: Section 7.2 sweeps it, so every result
% can be read at whatever value a given platform actually delivers.
cfg.OD = 0.8;

% --- 4.5 Residual calibration errors ------------------------------------
cfg.eps_sub = 1e-3;   % per-bin relative error on the known classical profile;
                      % frozen over a frame, hence a bias rather than noise
cfg.eps_kappa = 1e-3; % relative error on the fringe coefficient, collecting
                      % the uncertainty on the branch separation and on the
                      % trap frequency; it displaces the compression peak by
                      % kappa*eps*B_IF bins, which is sub-bin here
% Miscalibration of theta used to linearise the centroid leakage of Section
% 4.10. That section reports the calibration the emitter must achieve: a
% requirement placed on the apparatus, not a quantity optimised here.
cfg.eps_theta_probe = 1e-2;
cfg.F_excess     = 3.0;    % 1.0 is the pure shot-noise limit
cfg.sigma_f_track = 1e3;   % residual Doppler offset error, 1-sigma, in Hz
% Full-well capacity assumed for a spectral bin of the receiving array. It is
% not used to model the detector, only to make the no-saturation assumption of
% Section 4.8 checkable rather than tacit.
cfg.wellDepth = 1e6;       % electrons per bin, declared

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
% Epochs sample the pass uniformly and each is assigned an equal share of the
% window, so the total time is conserved exactly. Note that linspace places
% the first and last epochs at the window edges, where the link is worst, and
% weights them as full epochs: the discretisation is therefore mildly
% conservative. Availability is evaluated at the epoch centre only.
dt_ep    = (T_window/N_epochs)*ones(1,N_epochs);
avail_ep = double(link_available(idx_ep));

par = struct('theta',cfg.theta,'Omega_rad',Omega_rad, ...
             'eps_theta_probe',cfg.eps_theta_probe, ...
             'f0',f0_optical,'N',cfg.N,'Delta',cfg.Delta,'kappa',cfg.kappa, ...
             'OD',cfg.OD,'Gamma0',Gamma0_Hz, ...
             'F_excess',cfg.F_excess,'E_photon',E_photon,'eta_det',eta_det, ...
             'P_ep',P_rx_ideal_t(idx_ep),'dt_ep',dt_ep.*avail_ep, ...
             'SNR_thr',SNR_thr_Eve,'m',m_Al27,'c',c_light,'hbar',hbar);

% Only the packet separation is optimised. It is a property of the SIGNAL,
% not of the apparatus: it sets the spectral width of the line and hence the
% symbol rate, and the covert volume has an interior maximum in it. Every
% quantity describing the emitter itself - Delta, Delta_x, OD - is declared
% in Section 4 and swept for sensitivity in Section 7.
sep_vec = 4.0:0.05:9.5;
nSep = numel(sep_vec);
B_naive = zeros(1,nSep);   % covert volume against the energy detector
B_keyed = zeros(1,nSep);   % covert volume against the key-aware detector
B_noEve = zeros(1,nSep);   % volume with no eavesdropper present
for k = 1:nSep
    o = operatingPoint(sep_vec(k), par);
    if o.nyquistPhase < pi                 % reject an undersampled fringe
        B_naive(k) = o.volume_naive;
        B_keyed(k) = o.volume_keyed;
        B_noEve(k) = o.volume_noEve;
    end
end
[B_best, k_best] = max(B_naive);
if B_best <= 0
    error('No feasible separation: the fringe is undersampled over the whole sweep.');
end
cfg.sep_over_Delta = sep_vec(k_best);
opt = operatingPoint(cfg.sep_over_Delta, par);

fprintf('Optimal separation = %.2f Delta -> covert volume = %.3e bit\n', ...
        cfg.sep_over_Delta, B_best);
fprintf('  same point, no eavesdropper        = %.3e bit\n', opt.volume_noEve);
fprintf('  same point, key-aware eavesdropper = %.3e bit\n', opt.volume_keyed);
% This deflection is the DESIGN value: it uses the unfaded received power,
% because the optimisation must not depend on a particular pointing
% realisation. The Monte Carlo campaigns apply pointing loss and therefore
% report a lower value; Section 7.2 quotes the pass-median one. The two are
% both correct and must never be confused in the text.
fprintf('  median deflection d = %.3f (design, unfaded) -> raw BER = %.3f, capacity = %.4f bit/symbol\n', ...
        median(opt.d_ep), median(opt.BER_ep), median(opt.C_ep));

% --- Verification of the optimisation ------------------------------------
% The covert optimum has a closed form. Maximising sum_m n_m C(d) subject to
% sqrt(sum_m n_m d^4/(2N)) = thr, with C ~ d^2/(pi ln2) at small d, gives
%     d_opt = (2 N thr^2 / M_pass)^(1/4) .
% The numerical sweep must land near it; a large disagreement means the
% truncation or the deflection is wrong, not that the sweep found something
% clever.
d_opt = (2*cfg.N*SNR_thr_Eve^2/opt.M_pass)^(1/4);
fprintf('  check: analytic covert optimum d_opt = %.3f, sweep sits at %.3f (ratio %.2f)\n', ...
        d_opt, median(opt.d_ep), median(opt.d_ep)/d_opt);
% The optimum must be interior. If it sits on an end of sep_vec the sweep is
% reporting a boundary, not a maximum, and the range must be widened.
if k_best == 1 || k_best == nSep
    warning('Optimal separation sits on the edge of sep_vec: widen the sweep.');
else
    fprintf('  check: optimum interior to the sweep (%.2f in [%.2f, %.2f])\n', ...
            cfg.sep_over_Delta, sep_vec(1), sep_vec(end));
end
% Which constraint binds. If the covert volume equals the unconstrained one
% the eavesdropper never reaches threshold and the covert framing is inert:
% the optimiser is then maximising plain capacity and the Eve model has no
% influence on the design.
if opt.volume_naive >= 0.999*opt.volume_noEve
    fprintf('  check: BINDING CONSTRAINT IS PASS DURATION - the eavesdropper never reaches threshold,\n');
    fprintf('         so the covert model is inert at this operating point.\n');
else
    fprintf('  check: binding constraint is DETECTION - %.1f %% of the transmissible volume is usable\n', ...
            100*opt.volume_naive/opt.volume_noEve);
end

%% ------------------------------------------------------------------------
% 4.7 Frozen operating point and constraints on the fringe
% -------------------------------------------------------------------------
% Delta, Delta_x, kappa and OD are declared inputs (Section 4) and are not
% re-read from the optimisation; only the quantities the separation actually
% determines are frozen here.
cfg.p1     = -cfg.sep_over_Delta*cfg.Delta/2;
cfg.p2     = +cfg.sep_over_Delta*cfg.Delta/2;
cfg.B_IF   = opt.B_IF;
cfg.T_sym  = opt.T_sym;
cfg.dw     = opt.dw;
cfg.Phi    = opt.Phi;
cfg.key    = exp(+1i*cfg.kappa*cfg.dw);

% The in-phase template is the matched filter; the quadrature template is
% needed to apply a fringe-phase error without rebuilding the spectrum.
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
% Delta_x and sep are now independent declared inputs, so the fringe count is
% not fixed by sep alone. The equivalent branch separation that WOULD make the
% fringe span the same number of periods as the packet separation is reported
% for reference: it says whether the declared Delta_x is coarse or fine
% relative to the spectral structure it has to fill.
Delta_x_matched = cfg.Delta_x * (cfg.sep_over_Delta*(cfg.sep_over_Delta+6)/pi)/nFringe;
fprintf('  for reference, %.1f fringes would need Delta_x = %.1f nm\n', ...
        cfg.sep_over_Delta*(cfg.sep_over_Delta+6)/pi, Delta_x_matched*1e9);

% Ceiling 1: Nyquist sampling of the fringe on the spectral grid.
phasePerSample = cfg.kappa*dw_step;
fprintf('Ceiling 1 (Nyquist): fringe phase per sample = %.4f rad (limit pi), margin x%.0f\n', ...
        phasePerSample, pi/phasePerSample);

% Ceiling 2: the fringes must stay resolvable against the natural Lorentzian
% width of the transition, otherwise the convolution washes them out.
fprintf('Ceiling 2 (natural linewidth): fringe spacing %.3f MHz vs Gamma0 = %d Hz, margin x%.0f\n', ...
        fringeSpacing_Hz/1e6, Gamma0_Hz, fringeSpacing_Hz/Gamma0_Hz);
fprintf('Ceiling 3 (coherence time of the superposition) is experimental and outside this model.\n');

fprintf('Fringe-coefficient error: eps_kappa = %.1e -> peak displacement = %.3f bins\n', ...
        cfg.eps_kappa, cfg.eps_kappa*cfg.kappa*cfg.B_IF);

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

% A residual spectral offset s does NOT displace the compressed peak here:
% with a linear fringe phase it rotates the peak by kappa*s, so the real
% part is scaled by cos(kappa*s). This is the main robustness gain over the
% quadratic-phase version of this work, where the same offset moved the read
% bin. The averaged loss for a zero-mean Gaussian offset is exp(-sigma^2/2).
cfg.sigma_phi_track = cfg.kappa*2*pi*cfg.sigma_f_track;
fprintf('Tracking offset: sigma_f = %.0f Hz -> fringe phase error %.4f rad -> %.3f dB loss\n', ...
        cfg.sigma_f_track, cfg.sigma_phi_track, ...
        -20*log10(exp(-cfg.sigma_phi_track^2/2)));
fprintf('  offset tolerated for 1 dB of loss: %.1f kHz\n', ...
        sqrt(2*log(10^(1/20)))/(cfg.kappa*2*pi)/1e3);

% --- Verification of the noise model -------------------------------------
% sigma is recomputed here from the photon count by an independent route, and
% its scaling with range is checked: the collected power falls as 1/d^2 and
% the relative shot noise therefore rises linearly with range.
n_bin_chk = (P_rx_ideal_t(idx_cross)*cfg.T_sym/E_photon)*eta_det/cfg.N;
sig_chk   = cfg.F_excess/sqrt(n_bin_chk);
fprintf('  check: photons per bin at closest approach = %.2e -> sigma = %.2e (reported %.2e)\n', ...
        n_bin_chk, sig_chk, min(opt.sigma_ep));
[~,i_far] = max(d_t(idx_ep));
fprintf('  check: sigma edge/closest = %.2f against range ratio %.2f (shot noise scales as range)\n', ...
        max(opt.sigma_ep)/min(opt.sigma_ep), d_t(idx_ep(i_far))/d_min);
% Photons per symbol, the quantity an experimentalist recognises.
fprintf('  photons per symbol at closest approach: %.2e over %d bins\n', ...
        n_bin_chk*cfg.N, cfg.N);
% DECLARED ASSUMPTION ON THE RECEIVER. At this flux the excess-noise factor
% of 3 is conservative rather than generous: shot noise at 0.6 mW received
% corresponds to a relative intensity noise of -146 dBc/Hz, which sits above
% the RIN of an ordinary laser, so the fundamental term still dominates.
% What the flux does constrain is dynamic range. Each spectral bin collects
% of order 1e7 to 1e8 photons per symbol, tens of times the well depth of a
% conventional detector array, so the receiver is ASSUMED to handle the
% received flux without saturating - by dynamic range, by deliberate
% attenuation, or by a shorter integration. Which of these is a question of
% apparatus and is outside the perimeter of this model; the assumption
% belongs in the list of declared hypotheses in the text.

% --- Detector dynamic range: the assumption made checkable ---------------
% An assumption that is never evaluated is indistinguishable from an
% oversight, so the paragraph above is turned into a test. The binding epoch
% is the one collecting the MOST light, that is the closest approach taken
% unfaded: a pointing fade only relieves saturation, so ignoring it is the
% conservative choice. The test does not model the detector. It reports the
% factor by which the flux exceeds a conventional array, which is the factor
% that dynamic range, deliberate attenuation or a shorter integration has to
% absorb, and the price that absorbing it by attenuation alone would carry.
n_bin_max = (max(P_rx_ideal_t)*cfg.T_sym/E_photon)*eta_det/cfg.N;
satRatio  = n_bin_max/cfg.wellDepth;
fprintf('  saturation: %.2e electrons per bin per symbol vs declared well depth %.1e -> x%.0f\n', ...
        n_bin_max, cfg.wellDepth, satRatio);
if satRatio > 1
    fprintf('              absorbing it needs %.1f dB of attenuation, or an integration cut to %.2f %% of the symbol\n', ...
            10*log10(satRatio), 100/satRatio);
    fprintf('              if absorbed by attenuation alone the deflection would fall by x%.1f\n', ...
            sqrt(satRatio));
    warning(['Flux per bin exceeds the declared well depth by a factor %.0f: the ' ...
             'no-saturation assumption is not self-evident and must be stated in the text.'], satRatio);
else
    fprintf('              the array is not saturated at the declared well depth; the assumption is inert here\n');
end

%% ------------------------------------------------------------------------
% 4.9 Adversary models
% -------------------------------------------------------------------------
% Two adversaries are modelled. Both are assumed to receive the signal with
% the same signal-to-noise ratio as the intended receiver, which is the
% conservative choice: Section 4.11 argues they must in fact be inside the
% beam to receive anything at all.
%
% Eve-A, KEY-UNAWARE. She knows the classical profile, subtracts it, and sums
%   the residual energy over the band. Under the null hypothesis the residual
%   is N bins of variance sigma^2, so the summed energy has mean N sigma^2 and
%   standard deviation sigma^2 sqrt(2N). The signal adds sum(tmpl^2) =
%   d^2 sigma^2. Hence
%       delta_A = d^2 sigma^2 / (sigma^2 sqrt(2N)) = d^2/sqrt(2N).
%   She pays the full bandwidth in noise for a signal that occupies a small
%   part of it, because without kappa she cannot compress it.
%
% Eve-B, KEY-AWARE (Kerckhoffs). She knows kappa - it depends only on the
%   branch separation and atomic constants, and is in any case estimable from
%   the received spectrum by autocorrelation - and applies the same matched
%   filter as the receiver, measuring |z|^2. That collapses the N bins to one:
%       delta_B = d^2/sqrt(2).
%
% The ratio is therefore exactly sqrt(N), and the simulator checks it. The
% receiver's only advantage over Eve-B is knowledge of phi, which reads the
% bit but is not needed to detect that a transmission is taking place. It
% follows that whenever the link decodes reliably, Eve-B detects it. This is
% the central negative result of the work and must not be softened.
%
% WHAT THIS PAIR DOES NOT COVER. Both are quadratic detectors on the same
% observable. An adversary using a different observable - the line centroid,
% treated separately in Section 4.10 - or a different statistic, or one who
% accumulates across several passes, is outside the model. The two chosen
% here bracket the key-knowledge axis, not the space of all strategies.
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
if abs(delta_B/delta_A - sqrt(cfg.N)) > 1e-6*sqrt(cfg.N)
    warning('Adversary ratio departs from sqrt(N): the detector models are inconsistent.');
end

% --- Verification of the adversary models --------------------------------
% Eve-A's deflection is measured by Monte Carlo on the energy statistic
% itself, not rebuilt from the closed form, which would be circular. Noise of
% unit variance is drawn with and without a signal scaled to the deflection
% d_cross; the shift of the summed energy divided by its null spread must
% reproduce d^2/sqrt(2N).
rsEve   = RandStream('mrg32k3a','Seed',SEED_MASTER); rsEve.Substream = 9;
nTrial  = 2000;
s_unit  = sig.tmpl*(d_cross/sqrt(sum(sig.tmpl.^2)));   % ||s_unit|| = d_cross
E_null  = zeros(1,nTrial); E_sig = zeros(1,nTrial);
for it = 1:nTrial
    z = randn(rsEve,cfg.N,1);
    E_null(it) = sum(z.^2);
    E_sig(it)  = sum((z+s_unit).^2);
end
delta_A_mc = (mean(E_sig)-mean(E_null))/std(E_null);
fprintf('  check: Eve-A measured on the energy statistic = %.3e vs closed form %.3e (ratio %.3f)\n', ...
        delta_A_mc, delta_A, delta_A_mc/delta_A);
if abs(delta_A_mc/delta_A - 1) > 0.15
    warning('Measured Eve-A deflection departs from the closed form by more than 15 %%.');
end
% Symbols to threshold. If Eve-B crosses in less than one symbol the
% accumulated-statistic model is being extrapolated below the regime where
% it means anything: the honest statement is then that she detects the very
% first symbol, not that she needs a fraction of one.
nB = (SNR_thr_Eve/delta_B)^2;
if nB < 1
    fprintf('  check: Eve-B crosses threshold within the FIRST symbol (%.3f), so the\n', nB);
    fprintf('         fractional figure above is a formal extrapolation, not a count.\n');
end
% Sanity of the covert framing: how much of the pass survives Eve-A.
fprintf('  check: Eve-A tolerates %.2e symbols of %.2e in the pass (%.1f %%)\n', ...
        (SNR_thr_Eve/delta_A)^2, opt.M_pass, 100*(SNR_thr_Eve/delta_A)^2/opt.M_pass);

%% ------------------------------------------------------------------------
% 4.10 Leakage through the classical Doppler centroid
% -------------------------------------------------------------------------
% A non-zero delta_Q would shift the line centroid as a function of the bit,
% corrupting the classical velocimetry of the link and opening a third
% detection channel that exploits the entire photon count of the line. At
% theta = pi/4 it vanishes identically. Since theta cannot be calibrated
% exactly, the residual is evaluated both per symbol and accumulated over
% the whole pass, the latter being the binding requirement.
%
% WARNING ON TWO ESTIMATORS. The version below evaluates the centroid
% precision at CLOSEST APPROACH and then multiplies by sqrt(M_pass), as if
% every symbol of the pass were transmitted at maximum received power. That
% overstates the leakage, because most epochs collect far fewer photons. The
% constraint enforced inside operatingPoint accumulates epoch by epoch and
% is the correct one; the two disagreed by a factor of about two, with the
% pessimistic version reported here. Both are printed so the discrepancy is
% visible rather than buried, and the per-epoch figure is the one used for
% the design and the one to quote.
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
fprintf('  REQUIREMENT (pessimistic, closest-approach photon count): %.2e rad\n', eps_theta_req);
fprintf('  REQUIREMENT (accumulated epoch by epoch, the correct one): %.2e rad\n', ...
        opt.eps_theta_req);

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

% --- Verification of the signal construction -----------------------------
% The spectrum is rebuilt from first principles and compared with what
% operatingPoint produced, so a wrong centre, width or mapping shows here.
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
% Envelope width, from the second moment of the interference term, against
% the analytic value Delta/sqrt(2) mapped through the Doppler relation.
wI     = sig.I_env/sum(sig.I_env);
sig_meas = sqrt(sum(wI.*(cfg.dw - sum(wI.*cfg.dw)).^2));
fprintf('Interference envelope width: measured %.3f MHz, analytic %.3f MHz\n', ...
        sig_meas/(2*pi)/1e6, sigma_env_rad/(2*pi)/1e6);
% Fraction of the band actually occupied, and the emission weight of Eq.(26).
fprintf('Band occupancy: structure spans %.2f MHz of a %.2f MHz band\n', ...
        2*(cfg.sep_over_Delta/2+3)*cfg.Delta*f0_optical/1e6, cfg.B_IF/1e6);
fprintf('Emission weight of Eq.(26): departs from unity by %.1e over the band\n', ...
        3*cfg.B_IF/(2*f0_optical));

weight = sig.I_env.*cfg.key;                % envelope weighting and fringe key
k_expected = cfg.N/2 + 1;
Y_ref = fftshift(fft(ifftshift(sig.tmpl.*weight)));

% The read bin is fixed analytically, not found by search. Writing
% cos(Phi) = (exp(i Phi) + exp(-i Phi))/2, the product tmpl.*weight is
% (1/2) I_env^2 (1 + exp(2i kappa dw)): a peak at zero frequency and an
% image of EQUAL height displaced by 2 kappa B_IF bins. With a quadratic
% phase the image stayed spread and was harmless; with a linear fringe it is
% a second compressed peak, and an argmax over |real(Y)| picks one of the two
% at random. Only the zero-frequency bin satisfies
%     real(Y(k_expected)) = sum(tmpl.^2) = ||s||^2,
% which is the matched-filter statistic the deflection is derived from.
cfg.kpk = k_expected;
imageOffset = 2*cfg.kappa*cfg.B_IF;
peakCheck   = real(Y_ref(cfg.kpk))/sum(sig.tmpl.^2);
fprintf('Compression bin fixed at %d; image peak at %+.0f bins\n', ...
        cfg.kpk, imageOffset);
fprintf('Matched-filter identity real(Y(kpk))/||s||^2 = %.6f (must be 1)\n', peakCheck);
if abs(peakCheck-1) > 1e-3
    warning('Matched-filter identity violated: check kappa and the spectral grid.');
end

% The image must stay clear of the read bin. The compressed peak has the
% width of the transform of I_env^2, whose spectral width is sigma_env/sqrt(2);
% in bins that is N/(2 pi sigma_env/(sqrt(2) dw_step)). If the image, at
% 2 kappa B_IF bins, falls inside that width the two overlap and the read
% bin picks up energy that does not belong to it. The margin shrinks as
% Delta_x is reduced, so it must be checked and not assumed.
peakWidth_bins = cfg.N/(2*pi*((sigma_env_rad/sqrt(2))/dw_step));
fprintf('Peak/image separation: image at %.0f bins vs peak half-width %.0f bins, margin x%.1f\n', ...
        imageOffset, peakWidth_bins, imageOffset/peakWidth_bins);
if imageOffset < 3*peakWidth_bins
    Delta_x_clear = cfg.Delta_x*3*peakWidth_bins/imageOffset;
    warning(['Image peak within three half-widths of the read bin (margin %.1f). ' ...
             'Raising Delta_x to %.0f nm would clear it.'], ...
             imageOffset/peakWidth_bins, Delta_x_clear*1e9);
end

% Spreading check. An adversary who subtracts S_cl and sums the SIGNED
% residual, without knowing kappa, collects sum(tmpl); the matched filter
% collects ||tmpl||. The ratio measures how effectively the fringe hides the
% signature from a coherent but key-unaware detector. With a smooth
% single-lobe interference term the two would be comparable.
signedSum_dB = 20*log10(abs(sum(sig.tmpl))/(sqrt(sum(sig.tmpl.^2))*sqrt(cfg.N)));
fprintf('Fringe spreading: signed-sum detector is %.1f dB below the matched filter\n', ...
        -signedSum_dB);

% Matched-filter gain over the unweighted DC-bin statistic. The useful
% signal occupies only sqrt(2 pi) sigma_bin samples out of N, so an
% unweighted sum collects noise from all N samples for no additional signal.
g_mf   = sqrt(sum(sig.tmpl.^2));                          % ||s||
g_flat = sum(sig.tmpl.*real(cfg.key))/sqrt(cfg.N/2);      % unweighted statistic
fprintf('Matched-filter gain over the unweighted DC bin: %+.2f dB\n', ...
        20*log10(g_mf/max(abs(g_flat),eps)));
fprintf('Useful samples under the envelope: %.0f out of N = %d\n', ...
        sqrt(2*pi)*sigma_env_rad/dw_step, cfg.N);

% Equivalence of the two Monte Carlo paths. The deterministic parts of the
% statistic are compared directly: if the shortcut had a sign or convention
% error it would show here rather than as a silently wrong error rate.
A_sum  = sum(sig.tmpl.*real(weight));
A_fft  = real(Y_ref(cfg.kpk));
fprintf('Fast-path identity sum(tmpl.*Re w) / FFT bin = %.6f (must be 1)\n', A_sum/A_fft);
if abs(A_sum/A_fft - 1) > 1e-9
    warning('Sufficient-statistic shortcut disagrees with the transform.');
end

cfg.Nbit_frame = 1024;
cfg.N_MC       = 16;
cfg.block      = 32;

% RUNTIME. The full simulation builds a 16384-point noisy spectrum for every
% symbol, which costs of order 1e10 Gaussian variates and dominates the run.
% It is not necessary. The decision statistic is the zero-frequency bin of
% fft(Rx.*weight), which is exactly sum(Rx.*weight): the transform is a
% permutation and does not change the sum. Writing the statistic out,
%     stat = real(sum((S_cl - S_hat).*weight))          <- frame bias
%          + s*[A cos(dphi) - Aq sin(dphi)]             <- signal
%          + real(sum(noise.*weight)) ,                 <- noise
% with A = real(sum(tmpl.*weight)) and Aq = real(sum(tmpl_q.*weight)), the
% noise term is a weighted sum of independent Gaussians and is therefore
% itself Gaussian, with variance sigma^2 sum(real(weight).^2). One variate
% per symbol replaces N of them, exactly and not approximately.
% The full path is retained and both are compared in Section 4.12, because
% an exact shortcut that is silently wrong is worse than a slow one.
cfg.fastMC = true;

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
        [ne,nb] = runFrame(cfg, sig, sigma_vec(is), rs);
        errA(is) = errA(is)+ne; bitA(is) = bitA(is)+nb;
    end
    fprintf('  sigma = %.2e | d = %5.2f | BER simulated = %.3e | BER theory = %.3e\n', ...
            sigma_vec(is), g_mf/sigma_vec(is), errA(is)/bitA(is), BER_theory_A(is));
end
BER_A = errA./bitA; [BER_A_lo,BER_A_hi] = wilsonCI(errA,bitA,1.96);
% The comparison is restricted to points with enough errors for the ratio to
% be statistically meaningful. A twenty-error threshold was too loose: the
% Poisson relative deviation there is 22 %, so the reported maximum was
% dominated by counting noise at the single sparsest point and swung between
% 2 % and 44 % across runs with no change in the physics. The threshold is
% raised, and the deviation is also expressed in standard deviations of the
% binomial count, which is the quantity that actually indicates a mismatch.
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
% Statistical reach of this campaign, so the validation is not over-read.
% With Nbit_frame*N_MC bits per point the smallest error rate that can be
% resolved at all is of order 1/n; below that a zero count is expected and
% carries no information.
n_per_pt = cfg.Nbit_frame*cfg.N_MC;
fprintf('  campaign size: %d bits per point -> resolvable BER floor ~%.1e\n', ...
        n_per_pt, 1/n_per_pt);

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
        [ne,nb] = runFrame(cfg, sig, sigma_eff, rs);
        errB(ie) = errB(ie)+ne; bitB(ie) = bitB(ie)+nb;
    end
    sigma_B(ie) = acc/cfg.N_MC;
    BER_B(ie)   = errB(ie)/bitB(ie);
    fprintf('  t = %6.1f s | range = %7.1f km | sigma = %.2e | d = %5.2f | BER = %.3e\n', ...
            t_ep(ie), d_t(ii)/1000, sigma_B(ie), g_mf/sigma_B(ie), BER_B(ie));
end
[BER_B_lo,BER_B_hi] = wilsonCI(errB,bitB,1.96);
% Note 5 requires this check rather than assuming it. Near BER = 0.5 the
% capacity is not distinguishable from zero at any practical sample size, so
% each epoch is tested: the capacity is meaningful only where the measured
% error rate is separated from one half by more than its Wilson half-width.
sepFromHalf = (0.5 - BER_B)./max(BER_B_hi-BER_B, eps);
nWeak = sum(sepFromHalf < 2 & avail_ep > 0);
fprintf('  epochs whose capacity is within 2 sigma of zero: %d of %d\n', ...
        nWeak, sum(avail_ep > 0));
if nWeak > 0
    fprintf('  -> the capacity of those epochs is an upper bound, not a measurement (Note 5).\n');
end

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
        [ne,nb] = runFrame(cfg_tmp, sig, sigma_ref, rs);
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
fprintf('--- 7.2 Imprint depth (deflection at the pass-median noise, faded) ---\n');
% The sweep is clipped at unity: OD is a fractional modulation depth and a
% value above one has no physical meaning. With a declared OD near the top of
% the range the multipliers are chosen to bracket it from below.
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
% Delta_x is the third declared input, and the perimeter promises a
% sensitivity analysis on all three. It is not swept like Delta and OD,
% because it barely moves the deflection: once several fringes fit in the
% band the matched-filter energy sum((I_env cos(kappa dw))^2) approaches half
% of sum(I_env^2) and stops depending on kappa. What Delta_x decides is
% whether the fringe is READABLE at all, and it is squeezed between three
% constraints closing from opposite sides:
%   from below - the image peak sits at 2 kappa B_IF bins, and as Delta_x
%                shrinks it walks into the read bin (Section 4.12);
%   from above - the fringe phase per spectral sample approaches pi and the
%                grid can no longer sample the fringe (Ceiling 1);
%   from above - the fringe spacing approaches the natural linewidth and the
%                Lorentzian convolution washes the contrast out (Ceiling 2).
% The sweep locates the interval and reports where the declared value sits
% inside it. The guard factors are declared, not derived: three half-widths
% for the image, and a factor ten on the linewidth.
fprintf('--- 7.3 Branch separation: admissible interval ---\n');
g_image = 3;     % required image-to-read-bin separation, in peak half-widths
g_gamma = 10;    % required ratio of fringe spacing to natural linewidth
kap_of  = @(dx) m_Al27*c_light*dx/(hbar*Omega_rad);
% Analytic bounds, each inverted from its constraint.
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
% Numerical sweep. The deflection is recomputed at each Delta_x so the claim
% that it is insensitive is measured rather than asserted.
Dx_vec = cfg.Delta_x*[0.25 0.5 0.75 1 2 5 20];
fprintf('  %10s %10s %10s %10s %10s %10s\n', ...
        'Dx [nm]','fringes','img/width','Nyq [rad]','sp/Gamma0','d/d_ref');
d_ref_dx = g_mf/sigma_ref;
for k = 1:numel(Dx_vec)
    kap_k  = kap_of(Dx_vec(k));
    tmpl_k = sig.I_env.*cos(kap_k*cfg.dw);
    d_k    = sqrt(sum(tmpl_k.^2))/sigma_ref;
    fprintf('  %10.1f %10.1f %10.2f %10.4f %10.0f %10.4f\n', ...
            Dx_vec(k)*1e9, cfg.B_IF*kap_k, (2*kap_k*cfg.B_IF)/peakWidth_bins, ...
            kap_k*dw_step, 1/(kap_k*Gamma0_Hz), d_k/d_ref_dx);
end
fprintf('  the deflection is flat at and above the declared value: Delta_x sets readability,\n');
fprintf('  not signal strength. Below it the column rises, because with few fringes cos^2 no\n');
fprintf('  longer averages to one half - but that gain is unusable, since it is exactly where\n');
fprintf('  the image peak walks into the read bin.\n');

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
fprintf('Branch separation Delta_x                : %8.3f um\n', cfg.Delta_x*1e6);
fprintf('Fringe spacing / fringes in band         : %8.3f MHz / %.1f\n', ...
        1/cfg.kappa/1e6, cfg.B_IF*cfg.kappa);
fprintf('Optimised packet separation              : %8.2f Delta\n', cfg.sep_over_Delta);
fprintf('Imprint depth OD (declared)              : %8.3f\n', cfg.OD);
fprintf('Theta calibration required of the source : %8.2e rad\n', opt.eps_theta_req);
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

%% ------------------------------------------------------------------------
% 8.1 Control case: stationary geometry
% -------------------------------------------------------------------------
% A co-orbiting pair flying in formation would hold the link permanently at a
% fixed range. That is not a second scenario here but a CONTROL: the same
% chain is re-optimised with the received power frozen at its closest-approach
% value, full availability, and the same total dwell, so the ONLY difference
% from the simulated pass is the variability of the geometry. The comparison
% says how much of the result is owed to that variability and how much would
% have been obtained anyway. The stationary case is deliberately generous - it
% grants the best epoch of the pass for the whole window - so it is an upper
% bound on what a formation-flying pair would deliver, not a prediction.
fprintf('\n--- 8.1 Control: stationary geometry at closest-approach range ---\n');
par_stat = par;
par_stat.P_ep  = P_rx_ideal_t(idx_cross)*ones(1,N_epochs);
par_stat.dt_ep = (T_window/N_epochs)*ones(1,N_epochs);
B_stat = zeros(1,nSep);
for k = 1:nSep
    o_s = operatingPoint(sep_vec(k), par_stat);
    if o_s.nyquistPhase < pi, B_stat(k) = o_s.volume_naive; end
end
[Bs_best, ks_best] = max(B_stat);
o_stat = operatingPoint(sep_vec(ks_best), par_stat);
fprintf('  optimal separation %.2f Delta (pass: %.2f) -> covert volume %.3e bit vs %.3e for the pass\n', ...
        sep_vec(ks_best), cfg.sep_over_Delta, Bs_best, B_best);
fprintf('  ratio stationary/pass = %.2f ; median deflection %.3f vs %.3f (design, unfaded)\n', ...
        Bs_best/B_best, median(o_stat.d_ep), median(opt.d_ep));
fprintf('  no eavesdropper: %.3e bit vs %.3e ; key-aware Eve-B: %.3e bit vs %.3e\n', ...
        o_stat.volume_noEve, opt.volume_noEve, o_stat.volume_keyed, opt.volume_keyed);
% Reading the control. A ratio near unity would say that the geometry of the
% pass is incidental and the result is essentially that of a fixed link; a
% ratio well above unity says the pass is paying for its worst epochs, and
% that the cross-plane fly-by is the conservative choice claimed in the text.
if Bs_best/B_best > 1.5
    fprintf('  -> the pass is the conservative case: a fixed link at best range would carry x%.1f more.\n', ...
            Bs_best/B_best);
else
    fprintf('  -> the variability of the pass costs little: the result is close to that of a fixed link.\n');
end

%% ------------------------------------------------------------------------
% 8.2 Scaling of the covert volume with the length of the window
% -------------------------------------------------------------------------
% Covert communication over an additive-noise channel obeys a square-root
% law: the transmissible volume grows as the square root of the number of
% channel uses, not linearly, because the deflection must be dialled DOWN as
% the window lengthens in order to keep the accumulated detection statistic at
% threshold. The closed form checked in Section 4.6,
% d_opt = (2 N thr^2/M_pass)^(1/4), carries that law inside it. This block
% measures the exponent instead of quoting it, and the measurement matters: at
% the actual operating point it is NOT one half, and the text must be explicit
% about which of the two numbers it is using.
%
% HOW THE MEASUREMENT IS SET UP. The law concerns the trade between amplitude
% and duration, so the amplitude must be the free variable. The packet
% separation is frozen at its design value, because letting the optimiser move
% it as well would confound the fringe bandwidth with the covert trade, and
% the imprint depth is scanned as a proxy for the deflection, which is linear
% in OD. The per-epoch link conditions are held exactly as the pass delivers
% them and only the DWELL is stretched, so this is a hypothetical window of a
% different length under the same geometry and not a different orbit. The law
% applies against the key-unaware detector, which is the defensible figure of
% this work; against Eve-B the threshold is crossed almost at once and the
% length of the window never enters.
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
        if o_s.volume_naive > best
            best = o_s.volume_naive;  bestd = median(o_s.d_ep);
        end
    end
    B_scale(is) = best;  d_scale(is) = bestd;
    fprintf('  window x%-4d (%7.0f s) -> covert volume %.3e bit at median deflection %.3f\n', ...
            scale_vec(is), scale_vec(is)*T_window, best, bestd);
end
% Local exponents rather than one global fit. The law is asymptotic in the
% small-deflection limit, so a single fit over a range that begins outside
% that limit would report a number belonging to neither end.
slopes = diff(log(B_scale))./diff(log(scale_vec));
fprintf('  local exponents between consecutive windows: ');  fprintf('%.3f ', slopes);  fprintf('\n');
fprintf('  at the operating point (first interval)  : %.3f\n', slopes(1));
fprintf('  in the asymptotic regime (last interval) : %.3f  (square-root law predicts 0.500)\n', ...
        slopes(end));
% The asymptotic value is the one the law predicts and the one that must
% agree. The value at the operating point is larger, and legitimately so: the
% deflection there is of order unity, so the capacity has not yet entered its
% quadratic regime, and the epochs of the pass are strongly unequal in
% collected power. Both are reported because quoting the asymptotic exponent
% as if it described this link would overstate how quickly a longer window
% stops paying.
if abs(slopes(end)-0.5) > 0.05
    warning(['Asymptotic covert-volume exponent is %.3f rather than one half: either the ' ...
             'window sweep does not reach the small-deflection limit, or the truncation ' ...
             'model is wrong.'], slopes(end));
end
if slopes(1)-slopes(end) < 0.02
    fprintf('  note: the operating point already sits in the asymptotic regime.\n');
else
    fprintf('  the operating point sits ABOVE the asymptotic law by %.3f in exponent, so a longer\n', ...
            slopes(1)-slopes(end));
    fprintf('  window pays somewhat better here than the square-root law alone would suggest.\n');
end

%% ========================================================================
% 9. FIGURES
% =========================================================================
f_MHz = cfg.dw/(2*pi)/1e6;

% -------------------------------------------------------------------------
% FIG 1: SPECTRAL CHIRP AND MATCHED FILTERING
% -------------------------------------------------------------------------
figure('Name','Fig 1 Spectral fringe and matched filtering','Color','w','Position',[80 100 1200 380]);

% (a) Spectra
subplot(1,3,1);
% The classical profile and the transmitted spectrum are indistinguishable at
% this scale, which is the point: the interference term is a small ripple on a
% large line. The difference is therefore plotted on its own axis, otherwise
% the figure reads as though there were no signal at all.
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

% (c) Waterfall.
% The vertical range is clipped at the resolvable floor, 1/(Nbit_frame*N_MC).
% Left unclipped it follows the analytic curve down to 1e-300, which pushes
% every measured point into a thin band at the top and makes an agreement of
% a few per cent look like a total mismatch. Points below the floor carry no
% information: a zero count is what the campaign size predicts there.
subplot(1,3,3);
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
title('(c) Error Rate Validation');
legend({'Monte Carlo', 'theoretical filter curve'}, 'Location', 'southeast', 'FontSize', 8);
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
OD_fine = logspace(log10(cfg.OD/20), log10(cfg.OD*8), 200);
d_fine  = (OD_fine/cfg.OD)*(g_mf/sigma_ref);
semilogx(OD_fine, 1-binaryEntropy(0.5*erfc(d_fine/sqrt(2))), 'LineWidth', 1.5); hold on;
xline(cfg.OD, 'k:', 'declared', 'LineWidth', 1.2);
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
    N = par.N;  Delta = par.Delta;  kappa = par.kappa;
    % Band: the two Doppler peaks sit at -+ sep*Delta/2 and each has width
    % Delta, so 4(sep+6)Delta spans the structure with six standard
    % deviations of margin on each side. The Lorentzian of Eq. (26) is
    % treated as a delta function: its width Gamma0 is four orders below the
    % fringe spacing, so the convolution would reduce the fringe contrast by
    % less than 0.2 per cent. Ceiling 2 in Section 4.7 guards this.
    B_IF  = 4*par.f0*(sep+6)*Delta;
    T_sym = N/B_IF;
    dw = linspace(-1,1,N).'*(2*pi*B_IF/2);
    p  = dw/par.Omega_rad;                  % velocity in units of v/c
    % kappa is a declared input (Section 4.3) and is independent of sep: the
    % fringe spacing is set by the branch separation, the width of the line
    % by the packet separation.
    % Momentum amplitudes of Eq. (A1), centred at -+ sep*Delta/2.
    A1 = exp(-((p+sep*Delta/2).^2)./(2*Delta^2))./(pi^0.25*sqrt(Delta));
    A2 = exp(-((p-sep*Delta/2).^2)./(2*Delta^2))./(pi^0.25*sqrt(Delta));
    % Emission weight of Eq. (26). At the velocities of a cold source p is of
    % order 1e-8 in units of v/c, so this weight is a 1e-8 correction and is
    % numerically inert; it is retained for fidelity to the reference, not
    % because it changes anything.
    w  = 1 + 3*p;
    % Splitting |psi_sup|^2 into the incoherent part [Eq. A2] and the cross
    % term. The normalisation N of Eq. (29) is omitted: it depends on cos(phi)
    % and therefore on the bit, whereas the physical probe power is fixed by
    % the laser and not by the atomic state. The residual bit-dependent
    % imbalance this leaves is exactly what the signed-sum check of Section
    % 4.12 measures, and the fringe drives it down.
    S_cl  = (cos(par.theta)^2*A1.^2 + sin(par.theta)^2*A2.^2).*w;
    I_env = (2*cos(par.theta)*sin(par.theta)*A1.*A2).*w;
    % The profile is normalised to unit root-mean-square and scaled by the
    % declared imprint depth, so sigma below is directly the relative shot
    % noise seen by the receiver.
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
    delta_A  = d_ep.^2/sqrt(2*N);                      % key-unaware detector
    delta_B  = d_ep.^2/sqrt(2);                        % key-aware detector

    % Centroid leakage: the calibration of theta the emitter must achieve.
    % Reported, not enforced - it is a requirement on the apparatus.
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

    out.B_IF = B_IF;  out.T_sym = T_sym;  out.dw = dw;
    out.S_cl = S_cl;  out.I_env = I_env;  out.Phi = Phi;
    out.sigma_ep = sigma_ep;  out.d_ep = d_ep;
    out.BER_ep = BER_ep;      out.C_ep = C_ep;
    out.M_pass = sum(n_sym);
    out.volume_noEve = bits_cum(end);
    out.volume_naive = volumeAtThreshold(bits_cum, eveA, par.SNR_thr, n_sym, C_ep, delta_A);
    out.volume_keyed = volumeAtThreshold(bits_cum, eveB, par.SNR_thr, n_sym, C_ep, delta_B);
    out.kappa   = kappa;
    out.nyquistPhase = kappa*(dw(2)-dw(1));
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

function [nerr,nbit] = runFrame(cfg, sig, sigma_noise, rs)
% Dispatcher: exact sufficient-statistic path when cfg.fastMC is set, full
% spectral path otherwise.
    if cfg.fastMC
        [nerr,nbit] = simulateFrameFast(cfg, sig, sigma_noise, rs);
    else
        [nerr,nbit] = simulateFrame(cfg, sig, sigma_noise, rs);
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
% the fringe key. The direct-current bin of fft(Rx .* I_env .* key) is then
% exactly sum(Rx .* I_env .* cos(Phi)), the matched filter for the real
% template. The transform is an ordinary FFT: with a linear fringe phase no
% fractional-domain step is required.
    M = cfg.Nbit_frame; N = cfg.N;
    bits = randi(rs,[0 1],1,M);

    % Calibration error on the classical profile, frozen over the frame and
    % therefore acting as a decision bias rather than as averaging noise.
    S_hat = sig.S_cl.*(1 + cfg.eps_sub*randn(rs,N,1));

    % Residual error on the fringe coefficient, frozen over the frame. It
    % displaces the compression peak by a fraction of a bin, which the fixed
    % read bin then pays for; the effect is produced by the transform itself
    % and needs no separate model.
    kappa_bob = cfg.kappa*(1 + cfg.eps_kappa*randn(rs,1,1));
    weight    = sig.I_env.*exp(+1i*kappa_bob*cfg.dw);

    nerr = 0;
    for i1 = 1:cfg.block:M
        i2 = min(i1+cfg.block-1,M); b = bits(i1:i2); nb = numel(b);
        s  = 1-2*b;                                   % bit 0 -> +1, bit 1 -> -1
        % Residual Doppler offset rotates the fringe phase rather than
        % displacing the peak: cos(Phi + dphi) is expanded on the in-phase
        % and quadrature templates so that no spectrum has to be rebuilt.
        dphi = cfg.sigma_phi_track*randn(rs,1,nb);
        Tx   = sig.tmpl.*cos(dphi) - sig.tmpl_q.*sin(dphi);
        Rx = sig.S_cl + Tx.*s + sigma_noise*randn(rs,N,nb);
        Rx = Rx - S_hat;
        % The ifftshift is required: the signal is centred at zero detuning,
        % that is at the centre of the array rather than at the first index.
        % Without it the centring is read as a phase ramp and the real part
        % of the compressed peak alternates sign from bin to bin.
        Y  = fftshift(fft(ifftshift(Rx.*weight,1),[],1),1);
        stat = real(Y(cfg.kpk,:));
        nerr = nerr + sum(double(stat<0) ~= b);
    end
    nbit = M;
end

function [nerr,nbit] = simulateFrameFast(cfg, sig, sigma_noise, rs)
% Same model as simulateFrame, evaluated on the sufficient statistic. Exact:
% see the derivation next to cfg.fastMC. Consumes a different number of
% random variates, so individual realisations differ from the full path
% while the statistics do not.
    M = cfg.Nbit_frame;  N = cfg.N;
    bits = randi(rs,[0 1],1,M);

    S_hat  = sig.S_cl.*(1 + cfg.eps_sub*randn(rs,N,1));
    kappa_bob = cfg.kappa*(1 + cfg.eps_kappa*randn(rs,1,1));
    weight = sig.I_env.*exp(+1i*kappa_bob*cfg.dw);
    wR     = real(weight);

    bias   = sum((sig.S_cl - S_hat).*wR);
    A      = sum(sig.tmpl  .*wR);
    Aq     = sum(sig.tmpl_q.*wR);
    nStd   = sigma_noise*sqrt(sum(wR.^2));

    s    = 1-2*bits;
    dphi = cfg.sigma_phi_track*randn(rs,1,M);
    stat = bias + s.*(A*cos(dphi) - Aq*sin(dphi)) + nStd*randn(rs,1,M);
    nerr = sum(double(stat<0) ~= bits);
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
%   generate the same fringe with a contrast many orders of magnitude larger,
%   at higher bandwidth and with no cryogenic ion trap. From a purely
%   communications standpoint the quantum mechanism confers no advantage.
%   Second, the information is carried by phi, the relative phase of the
%   preparation laser; the quantum state is the transducer, not the message.
%   The value of the scheme is therefore as a physical-signature channel and
%   as a quantitative feasibility study, not as a competitive modulation.
%
% NOTE 2 - WHY THE FRINGE IS LINEAR AND NOT QUADRATIC.
%   With the real Gaussian packets of Eq. (A1) one has Phi_1 = Phi_2 = 0 and
%   the interference term is a single smooth lobe modulated by cos(phi),
%   with no structure for a key to hide. Some relative phase between the
%   branches must therefore be supplied. Three candidates were examined.
%   (a) Free dispersion. exp(-i p^2 t/2 m hbar) is common to both branches
%       and cancels in the squared modulus; for a free particle |psi(p)|^2
%       is rigorously constant in time. It contributes nothing.
%   (b) A differential preparation delay, used in earlier versions of this
%       work, would give Phi_2 - Phi_1 = -p^2 Delta_t/(2 m hbar) and hence a
%       chirp. But the quadratic-in-p phase accumulated in free flight
%       depends only on the elapsed time, not on the branch momentum, so two
%       branches that exist between the same two instants acquire the same
%       coefficient and it cancels. Worse, inside a harmonic trap there is no
%       free dispersion at all: a coherent state rotates in phase space. A
%       genuine differential quadratic phase needs a branch-selective shear,
%       which is beyond the state of the art. Supplying the chirp downstream
%       on the probe beam does not rescue it either: a spectral phase applied
%       to the field is a common factor and leaves |E(w)|^2 unchanged, so it
%       creates no fringe and hides nothing from an intensity measurement.
%   (c) A relative spatial displacement. This is the one that works. Two
%       branches separated by Delta_x carry exp(i p Delta_x/hbar), a phase
%       linear in p and genuinely relative between them, so it survives in
%       |psi(p)|^2 as a uniform fringe. This is ordinary matter-wave
%       interferometry. How a given platform produces Delta_x - ballistic
%       separation, a pulse sequence, the turning point of a trapped mode -
%       is a question of apparatus: here Delta_x is a declared input and
%       kappa follows from it and from atomic constants alone.
%
% NOTE 3 - THE PERIMETER: WHAT THIS SIMULATOR DOES AND DOES NOT DECIDE.
%   This is a CHANNEL simulator. It takes an emitter that writes the quantum
%   Doppler signature on a probe beam and asks what covert channel can be
%   built on it. It does not design the emitter.
%   Declared inputs describing the source (Section 4), each swept for
%   sensitivity rather than optimised:
%     Delta, the velocity spread of the momentum packets;
%     Delta_x, the spatial separation between the branches;
%     OD, the depth of the imprint written on the probe.
%   Everything else is derived from those inputs and from the scenario:
%     kappa and the fringe spacing, from Delta_x and atomic constants;
%     B_IF and the symbol rate, from the width of the spectral structure;
%     the antipodal character of the signalling, from cos(phi) in Eq. (32);
%     the compression bin, which the fringe key places at the grid centre;
%     the detection noise, from the collected photon count;
%     the packet separation, from the optimisation of Section 4.6;
%     the interference contrast, from the packet overlap.
%   Also assumed. These are the complete list; nothing else is taken for
%   granted anywhere in the file, and each must appear in the text.
%   ON THE SOURCE
%     that a coherent superposition of two momentum packets with controllable
%       relative phase can be prepared at all;
%     that its coherence survives one symbol, here about one millisecond.
%       This is the most severe assumption in the chain and is NOT quantified:
%       it is the binding ceiling of the whole scheme;
%     that Delta_x can be produced. Kinematically a momentum difference
%       separates the branches, and a pulse sequence can set Delta_x
%       independently of the momentum splitting through
%       Delta_x = sum_j dk_j t_j, but no platform is chosen here and no
%       sequence is shown to be realisable under its constraints;
%     that the emitter transmits continuously over the window, with no duty
%       cycle lost to preparation, cooling or recalibration. A duty cycle
%       below one reduces the volume in direct proportion;
%     that theta is calibrated to the tolerance Section 4.10 reports.
%   ON THE RECEIVER
%     that detection is coherent - amplitude and phase, not intensity alone -
%       which requires a heterodyne receiver with a local oscillator locked to
%       the line after the classical Doppler shift is removed;
%     that it handles the received flux without saturating;
%     that the classical profile is known well enough to be subtracted, with
%       eps_sub the residual; how it is calibrated is not described;
%     that kappa is shared between transmitter and receiver, which presumes a
%       prior secure channel not treated here;
%     that symbol timing and clock recovery are perfect.
%   ON THE LINK AND THE ADVERSARY
%     F_excess, the excess-noise factor over the shot limit;
%     eps_sub and eps_kappa, the residual calibration errors;
%     sigma_f_track, the residual Doppler offset, and that the 26.7 GHz shift
%       itself is compensated from the ephemerides;
%     the optical hardware parameters (apertures, power, efficiencies) and
%       the deliberate beam expansion of Section 3;
%     the orbital scenario (altitude, inclination, fly-by geometry);
%     the eavesdropper's detection threshold, taken as 3, that she receives
%       with the same signal-to-noise ratio as the intended receiver, and
%       that the statistic accumulates over a SINGLE pass. Repeated passes
%       resume the accumulation and the per-pass figure is then not the
%       right quantity.
%
% NOTE 4 - THE RECEIVER MUST BE A MATCHED FILTER, AND WHY.
%   After the fringe key the useful signal is (1/2) I_env(w) exp(i phi), which
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
% NOTE 6 - CONSTRAINTS ON THE FRINGE, AND THE FEASIBILITY QUESTION.
%   Three ceilings apply, two of which the simulator evaluates:
%     (a) Nyquist sampling of the fringe on the spectral grid;
%     (b) the fringe spacing must remain above the natural linewidth Gamma0,
%         otherwise the Lorentzian of the transition washes the fringes out;
%     (c) the coherence time of the momentum superposition, which is
%         experimental and cannot be verified from within this model.
%   Both evaluated ceilings clear by wide margins at the operating point,
%   Nyquist by two orders and the natural linewidth by three. The binding
%   limit is therefore (c), the coherence of the momentum superposition,
%   which is a property of the source and is outside this model.
%   Note that the fringe count in the band is the product of two independent
%   declared inputs: the branch separation Delta_x sets the fringe spacing,
%   the packet separation sets the width of the line. A source delivering a
%   smaller Delta_x gives coarser fringes over the same band, which costs
%   spreading against a key-unaware detector but nothing in deflection.
%
% NOTE 7 - PHOTON BUDGET AND THE IMPRINTING ARCHITECTURE.
%   Direct fluorescence from a single emitter, collected across hundreds of
%   kilometres by a 0.30 m aperture, is far below one photon per second and
%   cannot close the link. The spectral response is therefore imprinted on a
%   transmitted probe beam by absorption; the photon budget crossing the link
%   is that of the laser, not of the atom.
%   The depth of that imprint, OD, is the single number by which the atomic
%   apparatus enters this model, and the receiver deflection is linear in it.
%   Two remarks belong in the text rather than here. First, the peak
%   absorption cross-section 3 lambda^2/(2 pi) applies to a lifetime-limited
%   line, while the line used here is broadened by the very momentum
%   distribution that carries the signal; since the integrated absorption is
%   conserved, the usable cross-section is diluted in the same proportion,
%   and any estimate of OD must account for it. Second, OD is not bounded by
%   a single emitter: absorption is additive over an ensemble, so the depth
%   achievable is a property of the source chosen, not of this model.
%   Section 7.2 therefore sweeps OD, and every result can be read off at
%   whatever value the platform actually delivers.
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
%   rest on the eavesdropper's ignorance of kappa: that coefficient depends
%   only on the branch separation and atomic constants, and is in any case
%   estimable from the received spectrum by autocorrelation over a
%   one-dimensional search space.
%   What the fringe does buy, and the simulator now measures it in Section
%   4.12, is suppression of a coherent but key-unaware detector: an
%   adversary who subtracts the classical profile and sums the signed
%   residual collects sum(tmpl) instead of ||tmpl||, and the fringe drives
%   that sum towards zero. Against a smooth single-lobe interference term
%   that detector would perform as well as the matched filter.
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
%   of order 1e5 bits, against about 4e5 bits transmissible with no
%   eavesdropper present. The binding constraint is therefore DETECTION, not
%   the duration of the pass: the key-unaware statistic reaches roughly 20
%   against a threshold of 3 by the end of the pass, so the transmission has
%   to stop early. This reverses the conclusion reached with the chirped
%   template, where the two volumes coincided, and the reversal must be
%   carried into the text. Against the key-aware detector the same
%   configuration yields of order 1e2 bits.
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
%   On runtime: the cost was dominated by Gaussian variate generation, not by
%   the transforms, and is now removed rather than mitigated. The decision
%   statistic is the zero-frequency bin of the weighted transform, which is
%   exactly the weighted sum of the spectrum; the noise contribution is a
%   weighted sum of independent Gaussians and is therefore a single Gaussian
%   of known variance. One variate per symbol replaces sixteen thousand. The
%   reduction is exact, not approximate, and Section 4.12 verifies the
%   deterministic part of the identity at run time. The full spectral path is
%   retained behind cfg.fastMC = false and should be run at least once
%   whenever the signal model changes. Note that the two paths consume
%   different numbers of random variates, so individual realisations differ
%   while the statistics agree.
%
% NOTE 10 - KNOWN LIMITATIONS OF THE MODEL.
%   The image term behaves differently from the chirped version and the
%   distinction matters. The transported quantity is real, so cos(Phi + phi)
%   contains both exp(i(Phi + phi)) and its conjugate. With a quadratic phase
%   the conjugate stayed chirped and its energy spread over the time-
%   bandwidth product. With a linear fringe it does not spread: it forms a
%   second compressed peak of equal height, displaced by 2 kappa B_IF bins,
%   about thirty-five at the present operating point against a compressed
%   half-width of sixteen. That margin of 2.2 trips the guard in Section 4.12
%   and is the one open item in the model: raising Delta_x to about 90 nm
%   moves the image clear and simultaneously restores the fringe count. This
%   is not a loss of deflection. The
%   zero-frequency bin still satisfies real(Y) = sum(tmpl.^2) exactly, so the
%   deflection ||s||/sigma is the true matched-filter value with no 3 dB
%   penalty; the image carries no independent information and combining the
%   two peaks would not improve it. The practical consequence is that the
%   read bin must be fixed analytically rather than located by argmax, since
%   the two peaks are degenerate in magnitude.
%   The subtraction-error sweep models an independent per-bin error, which
%   averages efficiently over the samples under the envelope and therefore
%   shows a wide margin. A structured error, from trap drift, residual
%   secular motion or local-oscillator drift, projects far more efficiently
%   onto the template and is not covered. It should be studied separately
%   before the calibration margin reported here is relied upon. The sweep is
%   also evaluated at the pass-median noise level, where the deflection is
%   already small; repeating it at the closest-approach noise level would
%   make the sensitivity visible.
%   The receiver is now an ordinary FFT correlator. The fractional Fourier
%   transform, which the earlier chirped version required, is the degenerate
%   case of angle zero here and is retained only as theoretical framing.
%   Finally, the residual Doppler model separates offset from rate: the
%   offset rotates the fringe phase and is modelled explicitly, while the
%   uncompensated rate is folded into eps_kappa. Given a Doppler rate of
%   hundreds of megahertz per second at closest approach, the rate
%   specification deserves an independent budget rather than absorption into
%   a single relative coefficient. Note that the sensitivity to both is far
%   lower than in the chirped version, since a frequency error no longer
%   moves the compression peak.
% =========================================================================