% =========================================================================
% SIMULATORE ISL QUANTUM-COVERT BASATO SULL'EFFETTO DOPPLER QUANTISTICO
% VERSIONE 2 - revisione del ricevitore, del criterio di prestazione,
%              del punto operativo, del modello di imprinting e del
%              modello di avversario.
%
% FISICA DI RIFERIMENTO
%   P. T. Grochowski, A. R. H. Smith, A. Dragan, K. Debski,
%   "Quantum time dilation in atomic spectra", PRR 3, 023053 (2021)
%   + estensione: ritardo differenziale di preparazione Delta_t fra i due
%     rami della sovrapposizione (NOTA Q2). Senza questa estensione il
%     termine di interferenza NON e' un chirp.
%
% CHE COSA CAMBIA RISPETTO ALLA VERSIONE 1 (le cinque correzioni)
%   (1) RICEVITORE. Bob non legge piu' il bin DC di fft(Rx.*chiave), che e'
%       la somma NON PESATA di tutti gli N campioni: raccoglieva rumore da
%       16384 campioni per estrarre segnale da ~544. Ora pesa per
%       l'inviluppo noto I_env prima del de-chirp, cioe' esegue il vero
%       filtro adattato. Guadagno ~10 dB in SNR (vedi 5Q, validazione
%       BER simulata contro Q(d) analitica).
%   (2) CRITERIO. La soglia "BER <= 1e-3" e' stata rimossa come figura di
%       merito: un canale covert DEVE girare a BER grezza alta con codifica
%       forte. Il criterio e' ora la CAPACITA' per simbolo e il volume
%       integrato sul passaggio. BER_soglia resta solo come diagnostica.
%   (3) PUNTO OPERATIVO. sep_su_Delta non e' piu' scelto a mano: viene
%       trovato risolvendo numericamente il problema di progetto
%           max_sep  sum_m n_m C(d_m)   s.v.  sqrt(sum_m n_m delta_m^2) <= SNR_thr
%       cioe' massimizzare il volume informativo del passaggio sotto il
%       vincolo che la statistica di rivelazione ACCUMULATA di Eve resti
%       sotto soglia (legge della radice quadrata).
%   (4) IMPRINTING E CALIBRAZIONE. Introdotto cfg.OD, la profondita' di
%       modulazione frazionaria del fascio di sonda (prima implicitamente
%       assunta = 1, cioe' modulazione al 100%: irrealistico per un singolo
%       ione). Introdotto cfg.eps_sub, errore relativo di conoscenza del
%       profilo classico da sottrarre, con sweep dedicato (FASE 6Q-bis).
%   (5) MODELLO DI AVVERSARIO. Due Eve dichiarate:
%         Eve-A "ignara": rivelatore d'energia senza chiave di de-chirp.
%         Eve-B "informata" (Kerckhoffs): conosce C e applica lo stesso
%               filtro adattato di Bob, misurando |z|^2. E' piu' sensibile
%               di Eve-A di un fattore esatto sqrt(N).
%       Sotto Eve-B lo schema NON e' covert: se Bob decodifica, Eve rileva.
%       L'argomento di sicurezza difendibile e' GEOMETRICO (Eve deve stare
%       fisicamente dentro il fascio FSO): quantificato in 4Q.12.
%
% ADD-ON MATLAB: Satellite Communications Toolbox (satelliteScenario,
% satellite, states) e Phased Array System Toolbox (fspl). Nessun altro.
% =========================================================================

clc; clear; close all;

SEED_MASTER = 42;
rng(SEED_MASTER, 'twister');

%% ========================================================================
% SEZIONE 0: COSTANTI FISICHE
% =========================================================================
c_light   = 299792458;
lambda_0  = 267e-9;                 % 27Al+ 1S0-3P1 (intercombinazione)
f0_optical= c_light / lambda_0;
Omega_rad = 2*pi*f0_optical;
h_planck  = 6.62607015e-34;
hbar      = 1.054571817e-34;
E_fotone  = h_planck * f0_optical;
u_amu     = 1.66053906660e-27;
m_Al27    = 27 * u_amu;
mc2_Al27  = m_Al27 * c_light^2;
Re        = 6378137;
mu_earth  = 3.986004418e14;
Gamma0_Hz = 520;                    % larghezza naturale della transizione

disp('=== SEZIONE 0 ===');
disp(['27Al+ 1S0-3P1: lambda = ',num2str(lambda_0*1e9,'%.1f'),' nm, f0 = ', ...
      num2str(f0_optical/1e12,'%.3f'),' THz, Gamma0 = ',num2str(Gamma0_Hz),' Hz']);

%% ========================================================================
% FASI 1-2: GEOMETRIA ORBITALE SGP4/J2 (invariata)
% =========================================================================
disp('--- FASI 1-2: Geometria orbitale SGP4/J2 ---');

startTime = datetime(2026,8,14,10,0,0);
stopTime  = startTime + minutes(15);
sampleTime= 0.1;
T_window  = seconds(stopTime - startTime);

altitude_orbit = 800e3;  sma = Re + altitude_orbit;
ecc = 0.0001; incl_deg = 60; argPeri_deg = 0; J2 = 1.08263e-3;
mu_n = sqrt(mu_earth/sma^3);
raan_rate = -1.5*mu_n*J2*(Re/(sma*(1-ecc^2)))^2*cosd(incl_deg);
dRAAN_deg = acosd((cosd(60)-cosd(incl_deg)^2)/sind(incl_deg)^2);
raan_Alice = 0; raan_Bob = dRAAN_deg;
T_orbit = 2*pi/mu_n; incl_rad = deg2rad(incl_deg);

posC = @(u0,raan0,t) sma*[ ...
  cos(raan0+raan_rate*t).*cos(u0+mu_n*t)-sin(raan0+raan_rate*t).*sin(u0+mu_n*t).*cos(incl_rad); ...
  sin(raan0+raan_rate*t).*cos(u0+mu_n*t)+cos(raan0+raan_rate*t).*sin(u0+mu_n*t).*cos(incl_rad); ...
  sin(u0+mu_n*t).*sin(incl_rad)];

tc = 0:2:T_orbit; pA = posC(0,0,tc); u0c = deg2rad(0:0.5:359.5);
bo = Inf; bu = 0; bt = 0;
for k = 1:numel(u0c)
    dk = vecnorm(posC(u0c(k),deg2rad(raan_Bob),tc)-pA,2,1);
    [dm,ik] = min(dk); ok = abs(dm-300e3);
    if ok < bo, bo = ok; bu = u0c(k); bt = tc(ik); end
end
tf = max(0,bt-10):0.02:min(T_orbit,bt+10); pAf = posC(0,0,tf);
u0f = bu + deg2rad(-1:0.002:1); bof = Inf; buf = bu; btf = bt;
for k = 1:numel(u0f)
    dk = vecnorm(posC(u0f(k),deg2rad(raan_Bob),tf)-pAf,2,1);
    [dm,ik] = min(dk); ok = abs(dm-300e3);
    if ok < bof, bof = ok; buf = u0f(k); btf = tf(ik); end
end
trueAnom_Alice = mod(rad2deg(mu_n*btf)-rad2deg(mu_n*T_window/2),360);
trueAnom_Bob   = mod(trueAnom_Alice+rad2deg(buf),360);

sc = satelliteScenario(startTime, stopTime, sampleTime);
satA = satellite(sc,sma,ecc,incl_deg,raan_Alice,argPeri_deg,trueAnom_Alice,"OrbitPropagator","sgp4","Name","Alice");
satB = satellite(sc,sma,ecc,incl_deg,raan_Bob,  argPeri_deg,trueAnom_Bob,  "OrbitPropagator","sgp4","Name","Bob");
[posAlice,velAlice,tS] = states(satA,"CoordinateFrame","inertial");
[posBob,  velBob,  ~ ] = states(satB,"CoordinateFrame","inertial");

t_sec = seconds(tS-tS(1)); N_t = numel(t_sec);
dPos = posBob-posAlice; d_t = vecnorm(dPos,2,1);
uLOS = dPos./d_t; v_rel_t = sum((velBob-velAlice).*uLOS,1);
f_D_t = -(v_rel_t/c_light)*f0_optical;
dfD_dt_t = gradient(f_D_t, sampleTime);

r_min = zeros(1,N_t);
for i = 1:N_t
    r_min(i) = norm(cross(posAlice(:,i),posBob(:,i)))/norm(posBob(:,i)-posAlice(:,i));
end
los_margin = r_min-(Re+80e3);

AU = 1.495978707e11; eps_ecl = deg2rad(23.439291);
nd = days(tS-datetime(2000,1,1,12,0,0,'TimeZone','UTC'));
lam_s = deg2rad(mod(280.460+0.9856474*nd,360));
posSun = AU*[cos(lam_s); sin(lam_s)*cos(eps_ecl); sin(lam_s)*sin(eps_ecl)];
uSun = (posSun-posBob)./vecnorm(posSun-posBob,2,1);
theta_sun = rad2deg(acos(min(max(sum((-uLOS).*uSun,1),-1),1)));
link_available = (los_margin>=0) & (theta_sun>=30);

[d_min, idx_cross] = min(d_t);
disp(['Fly-by minimo: ',num2str(d_min/1000,'%.2f'),' km a t = ',num2str(t_sec(idx_cross),'%.1f'),' s']);
disp(['Doppler ottico max: ',num2str(max(abs(f_D_t))/1e9,'%.3f'),' GHz | Doppler-rate al crossing: ', ...
      num2str(abs(dfD_dt_t(idx_cross))/1e6,'%.3f'),' MHz/s']);
disp(['Disponibilita'' link (LOS+Sole): ',num2str(100*mean(link_available),'%.1f'),' %']);

%% ========================================================================
% FASE 3: LINK BUDGET OTTICO (invariata)
% =========================================================================
disp('--- FASE 3: Link budget FSO ---');
D_tx=0.20; D_rx=0.30; P_tx=100e-3;
eta_tx=0.85; eta_rx=0.80; eta_det=0.65; rho_ap=0.65;
sigma_jit_pat=1.5e-6; theta_div=8.0e-6;

G_tx = rho_ap*(pi*D_tx/lambda_0)^2;
G_rx = rho_ap*(pi*D_rx/lambda_0)^2;
FSPL_t = 10.^(fspl(d_t,lambda_0)/10);
Lp_t = exp(-8*((sqrt((sigma_jit_pat*randn(1,N_t)).^2+(sigma_jit_pat*randn(1,N_t)).^2))./theta_div).^2) .* ...
       exp(-8*((sqrt((sigma_jit_pat*randn(1,N_t)).^2+(sigma_jit_pat*randn(1,N_t)).^2))./theta_div).^2);
P_rx_ideal_t = P_tx*eta_tx*eta_rx*(G_tx*G_rx)./FSPL_t;
P_rx_t = P_rx_ideal_t.*Lp_t;
disp(['P_rx al crossing: ',num2str(10*log10(P_rx_t(idx_cross)*1e3),'%.2f'),' dBm']);

%% ========================================================================
% FASE 4Q: STATO QUANTISTICO, IMPRINTING E PUNTO OPERATIVO
% =========================================================================
disp('--- FASE 4Q: stato quantistico, chirp derivato, punto operativo ---');
cfg = struct();

% --- 4Q.1 Trappola: Delta discende dallo stato fondamentale del moto -----
cfg.f_secolare = 3e6;
w_sec  = 2*pi*cfg.f_secolare;
cfg.Delta = sqrt(m_Al27*hbar*w_sec/2)/(m_Al27*c_light);   % in p/(mc)

% --- 4Q.2 theta = pi/4 (scelta invertita rispetto al paper, NOTA Q6) -----
% Segnale di Bob ~ sin(2*theta) MASSIMO; delta_Q (Eq.32) ~ sin(4*theta)
% ESATTAMENTE NULLO: il centroide di riga non si sposta col bit e il canale
% di rivelazione via tracking del centroide si chiude.
cfg.theta = pi/4;

% --- 4Q.3 Ritardo differenziale: origine fisica del chirp ---------------
% ATTENZIONE (limite fisico, non numerico): exp(-i p^2 t/2 m hbar) e'
% l'evoluzione di PARTICELLA LIBERA. In trappola armonica uno stato
% coerente ruota nello spazio delle fasi e NON accumula questa fase
% quadratica. Delta_t richiede quindi o ione rilasciato (i due rami si
% separano di sep*Delta*c*Delta_t, calcolato piu' sotto) o un potenziale
% state-dependent che congela un ramo. Va dichiarato in tesi.
cfg.Delta_t = 20e-6;                 % ritardo differenziale (s)
cfg.C_chirp = mc2_Al27*cfg.Delta_t/(2*hbar*Omega_rad^2);   % [s^2/rad^2]

% --- 4Q.4 (NUOVO, punto 4) PROFONDITA' DI IMPRINTING --------------------
% La v1 normalizzava S_cl a RMS 1 e derivava sigma dal conteggio fotonico:
% equivaleva ad assumere che lo ione modulasse il fascio di sonda al 100%.
% La modulazione frazionaria vera e' fissata dalla densita' ottica del
% singolo ione: OD ~ sigma_abs / A_fascio, con sigma_abs = 3*lambda^2/(2*pi)
% risonante e A_fascio l'area del fuoco sullo ione.
sigma_abs   = 3*lambda_0^2/(2*pi);        % sezione d'urto risonante
w0_fuoco    = 0.20e-6;                     % waist del fuoco sullo ione [m] (NA ~ 0.42)
A_fascio    = pi*w0_fuoco^2;
OD_max      = sigma_abs/A_fascio;          % limite superiore ideale
cfg.OD      = 0.15;                        % progetto: ~55% del limite di singolo ione
disp(['IMPRINTING: sigma_abs = ',num2str(sigma_abs,'%.2e'),' m^2, waist = ', ...
      num2str(w0_fuoco*1e6,'%.2f'),' um -> OD_max = ',num2str(OD_max,'%.3f'), ...
      '  | OD di progetto = ',num2str(cfg.OD,'%.2f')]);
if cfg.OD > OD_max
    warning(['OD di progetto superiore al limite di singolo ione (',num2str(OD_max,'%.3f'), ...
             '): ridurre cfg.OD o stringere il fuoco.']);
end

% --- 4Q.5 Calibrazioni residue (punto 4 e nota sul rate Doppler) --------
cfg.eps_sub = 1e-3;   % errore relativo per bin sul profilo classico noto
                      % (congelato sul frame: e' una calibrazione, non rumore)
cfg.eps_C   = 1e-3;   % errore relativo sul coefficiente di chirp. Raccoglie
                      % l'incertezza su Delta_t E il residuo di RATE Doppler
                      % non compensato, che entra come fase quadratica
                      % spuria sommandosi a C.
cfg.F_eccesso = 3.0;  % 1.0 = limite shot puro; 3.0 = caso realistico
cfg.sigma_f_track = 1e3;   % errore residuo di OFFSET Doppler, 1-sigma (Hz)

cfg.N = 2^14;
SNR_soglia_Eve = 3;        % soglia di decisione assunta per Eve

%% ------------------------------------------------------------------------
% 4Q.6 (NUOVO, punto 3) SWEEP DI PROGETTO SU sep_su_Delta
% -------------------------------------------------------------------------
% sep NON e' piu' scelto a mano. Il volume covert e'
%     B(sep) = sum_m n_m * C(d_m)      troncato quando
%     sqrt( sum_m n_m * delta_m^2 ) supera SNR_soglia_Eve
% con d_m = ||s||/sigma_m la deflessione del filtro adattato all'epoca m e
% delta_m la deflessione per simbolo del rivelatore di Eve. Il massimo si
% ottiene tipicamente al confine "Eve rileva esattamente a fine passaggio":
% sotto quel punto si spreca budget di rivelazione, sopra si tronca.
disp('--- 4Q.6: sweep di progetto su sep_su_Delta ---');

N_epoche = 21;
idx_ep   = round(linspace(1,N_t,N_epoche));
t_ep     = t_sec(idx_ep);
dt_ep    = (T_window/N_epoche)*ones(1,N_epoche);
disp_av  = double(link_available(idx_ep));

par = struct('Delta',cfg.Delta,'theta',cfg.theta,'Omega_rad',Omega_rad, ...
             'f0',f0_optical,'N',cfg.N,'C_chirp',cfg.C_chirp,'OD',cfg.OD, ...
             'F_eccesso',cfg.F_eccesso,'E_fotone',E_fotone,'eta_det',eta_det, ...
             'P_ep',P_rx_ideal_t(idx_ep),'dt_ep',dt_ep.*disp_av, ...
             'SNR_thr',SNR_soglia_Eve);

sep_vec = 4.0:0.05:9.5;
nSep = numel(sep_vec);
B_en = zeros(1,nSep); B_mf = zeros(1,nSep); B_noEve = zeros(1,nSep);
d_med = zeros(1,nSep); nyq_v = zeros(1,nSep); Tsym_v = zeros(1,nSep);
for k = 1:nSep
    o = puntoOperativo(sep_vec(k), par);
    nyq_v(k)  = o.nyq;
    Tsym_v(k) = o.T_sym;
    d_med(k)  = median(o.d_ep);
    if o.nyq >= pi
        B_en(k) = 0; B_mf(k) = 0; B_noEve(k) = 0;   % chirp sottocampionato
    else
        B_en(k) = o.bits_covert_en;
        B_mf(k) = o.bits_covert_mf;
        B_noEve(k) = o.bits_no_eve;
    end
end
[B_best, k_best] = max(B_en);
cfg.sep_su_Delta = sep_vec(k_best);
opt = puntoOperativo(cfg.sep_su_Delta, par);

fprintf('   sep ottimo = %.2f  ->  volume covert (Eve-A energia) = %.3e bit\n', ...
        cfg.sep_su_Delta, B_best);
fprintf('   allo stesso punto: volume se Eve NON ascolta = %.3e bit\n', opt.bits_no_eve);
fprintf('   allo stesso punto: volume covert contro Eve-B (chiave nota) = %.3e bit\n', ...
        opt.bits_covert_mf);
fprintf('   d mediana sul passaggio = %.3f -> BER grezza mediana = %.3f, C = %.4f bit/simb\n', ...
        median(opt.d_ep), median(opt.BER_ep), median(opt.C_ep));

% Confronto esplicito con la scelta manuale della v1
o_v1 = puntoOperativo(7.35, par);
fprintf('   [confronto, stesso OD] sep = 7.35 (scelta manuale v1): volume covert = %.3e bit, ', ...
        o_v1.bits_covert_en);
fprintf('d mediana = %.3f\n', median(o_v1.d_ep));

% --- 4Q.7 Congelamento del punto operativo ------------------------------
cfg.p1 = -cfg.sep_su_Delta*cfg.Delta/2;
cfg.p2 = +cfg.sep_su_Delta*cfg.Delta/2;
cfg.B_IF  = opt.B_IF;
cfg.T_sym = opt.T_sym;
cfg.dw    = opt.dw;
cfg.p_grid= cfg.dw/Omega_rad;
cfg.chiave= exp(+1i*cfg.C_chirp*cfg.dw.^2);
S_classico = opt.S_cl;      % gia' scalati per OD
I_env      = opt.I_env;
Phi_chirp  = opt.Phi;
R_b_raw    = 1/cfg.T_sym;
dw_max     = 2*pi*cfg.B_IF/2;
larghezza_Hz = f0_optical*(cfg.sep_su_Delta+6)*cfg.Delta;

disp(['Trappola f_sec = ',num2str(cfg.f_secolare/1e6,'%.1f'),' MHz -> Delta = ', ...
      num2str(cfg.Delta,'%.3e'),' (p/mc), v_spread = ',num2str(cfg.Delta*c_light,'%.4f'),' m/s']);
disp(['Separazione OTTIMIZZATA = ',num2str(cfg.sep_su_Delta,'%.2f'),' Delta -> dv = ', ...
      num2str((cfg.p2-cfg.p1)*c_light,'%.4f'),' m/s']);
disp(['B_IF derivata = ',num2str(cfg.B_IF/1e6,'%.2f'),' MHz -> T_sym = ', ...
      num2str(cfg.T_sym*1e6,'%.1f'),' us -> R_b = ',num2str(R_b_raw/1e3,'%.3f'),' kbit/s']);
disp(['Delta_t = ',num2str(cfg.Delta_t*1e6,'%.1f'),' us -> C = ',num2str(cfg.C_chirp,'%.4e'),' s^2/rad^2']);

% Escursione spaziale dei due rami durante Delta_t (verifica di fattibilita')
dx_rami = (cfg.p2-cfg.p1)*c_light*cfg.Delta_t;
disp(['VERIFICA FATTIBILITA'' Delta_t: in volo libero i due rami si separano di ', ...
      num2str(dx_rami*1e6,'%.1f'),' um in ',num2str(cfg.Delta_t*1e6,'%.1f'),' us.']);

% --- 4Q.8 Vincoli su Delta_t: Nyquist, coerenza, larghezza naturale -----
sigma_env_rad = (cfg.Delta/sqrt(2))*Omega_rad;
dw_step  = cfg.dw(2)-cfg.dw(1);
fase_per_campione = 2*cfg.C_chirp*(2*sigma_env_rad)*dw_step;
disp(['Nyquist del chirp: fase per campione al bordo (2 sigma) = ', ...
      num2str(fase_per_campione,'%.3f'),' rad -> Delta_t max ~ ', ...
      num2str(cfg.Delta_t*pi/fase_per_campione*1e6,'%.1f'),' us']);
% NUOVO: secondo tetto su Delta_t, indipendente da Nyquist. Le frange del
% chirp al bordo dell'inviluppo hanno passo 2*pi/(2*C*dw_bordo); se scendono
% sotto Gamma0 vengono lavate dalla lorentziana naturale della transizione.
passo_frange_Hz = (2*pi/(2*cfg.C_chirp*(2*sigma_env_rad)))/(2*pi);
disp(['Larghezza naturale: passo delle frange al bordo = ', ...
      num2str(passo_frange_Hz/1e3,'%.2f'),' kHz contro Gamma0 = ',num2str(Gamma0_Hz),' Hz ', ...
      '-> margine x',num2str(passo_frange_Hz/Gamma0_Hz,'%.0f'),' (Delta_t max ~ ', ...
      num2str(cfg.Delta_t*passo_frange_Hz/Gamma0_Hz*1e6,'%.0f'),' us)']);

BT = cfg.C_chirp*(4*sigma_env_rad)^2/pi;
disp(['Prodotto tempo-banda del chirp: BT = ',num2str(BT,'%.0f')]);

% Sensibilita' del de-chirp all'errore su C (Delta_t + rate Doppler residuo)
fase_res_max = cfg.eps_C*cfg.C_chirp*(2*sigma_env_rad)^2;
disp(['Errore su C: eps_C = ',num2str(cfg.eps_C,'%.1e'),' -> fase residua al bordo = ', ...
      num2str(fase_res_max,'%.3f'),' rad  (requisito: Delta_t noto a meglio di ', ...
      num2str(100/(cfg.C_chirp*(2*sigma_env_rad)^2),'%.3f'),' %% per 1 rad)']);

% --- 4Q.9 Rumore derivato dal conteggio fotonico -------------------------
n_bin_cross = (P_rx_t(idx_cross)*cfg.T_sym/E_fotone)*eta_det/cfg.N;
sigma_cross = cfg.F_eccesso/sqrt(n_bin_cross);
n_bin_edge  = (P_rx_ideal_t(1)*cfg.T_sym/E_fotone)*eta_det/cfg.N;
sigma_edge  = cfg.F_eccesso/sqrt(n_bin_edge);
disp(['Rumore derivato (F_eccesso = ',num2str(cfg.F_eccesso,'%.1f'),'): sigma crossing = ', ...
      num2str(sigma_cross,'%.2e'),', sigma bordo = ',num2str(sigma_edge,'%.2e')]);

% Tracking di OFFSET Doppler -> spostamento del bin di lettura
cfg.sigma_bin_track = 4*cfg.C_chirp*cfg.sigma_f_track*dw_max;
% Con la pesatura per I_env la funzione compressa e' la FT di I_env^2:
% larghezza sigma_b/sqrt(2) nel dominio spettrale -> picco piu' LARGO di
% sqrt(2) rispetto alla v1, quindi PIU' tollerante al jitter.
sigma_k_picco = cfg.N/(2*pi*((sigma_env_rad/sqrt(2))/dw_step));
disp(['Tracking offset: ',num2str(cfg.sigma_bin_track,'%.2f'),' bin di jitter contro ', ...
      'larghezza del picco sigma_k = ',num2str(sigma_k_picco,'%.1f'),' bin -> degrado ', ...
      num2str(-20*log10(exp(-cfg.sigma_bin_track^2/(2*sigma_k_picco^2))),'%.2f'),' dB']);

%% ------------------------------------------------------------------------
% 4Q.10 (NUOVO, punto 5) DUE MODELLI DI EVE, DICHIARATI
% -------------------------------------------------------------------------
% Eve-A "ignara": non possiede la chiave di de-chirp. Conosce il profilo
%   classico, lo sottrae e somma energie: T = sum (x - S_cl)^2.
%   Deflessione per simbolo: delta_A = d^2 / sqrt(2N).
% Eve-B "informata" (ipotesi di Kerckhoffs): conosce C - che dipende solo da
%   Delta_t, m, Omega, cioe' da grandezze fisiche, non da un segreto - e
%   puo' comunque stimarlo dallo spettro (autocorrelazione / Wigner) senza
%   conoscere il bit. Applica lo stesso filtro adattato di Bob e misura
%   |z|^2. Deflessione per simbolo: delta_B = d^2 / sqrt(2).
% RAPPORTO ESATTO: delta_B / delta_A = sqrt(N). Con N = 16384, Eve-B e' 128
% volte piu' sensibile. Conseguenza strutturale: se Bob decodifica (d ~ 3),
% Eve-B rileva su UN SOLO simbolo. Sotto Kerckhoffs lo schema non e' covert.
d_cross = interp1(t_ep, opt.d_ep, t_sec(idx_cross), 'linear', 'extrap');
delta_A = d_cross^2/sqrt(2*cfg.N);
delta_B = d_cross^2/sqrt(2);
fprintf('\n--- 4Q.10: modelli di avversario (valutati al crossing) ---\n');
fprintf('   deflessione di Bob (filtro adattato)     d = %.3f -> BER grezza = %.4f\n', ...
        d_cross, 0.5*erfc(d_cross/sqrt(2)));
fprintf('   Eve-A (ignara, energia)      delta/simbolo = %.3e -> rileva dopo %.2e simboli\n', ...
        delta_A, (SNR_soglia_Eve/delta_A)^2);
fprintf('   Eve-B (chiave nota, matched) delta/simbolo = %.3e -> rileva dopo %.2e simboli\n', ...
        delta_B, (SNR_soglia_Eve/delta_B)^2);
fprintf('   rapporto di sensibilita'' Eve-B/Eve-A     = %.1f  (atteso sqrt(N) = %.1f)\n', ...
        delta_B/delta_A, sqrt(cfg.N));

% --- 4Q.11 delta_Q: il Doppler classico e' perturbato? (invariato) ------
dp_pk  = cfg.p2 - cfg.p1;
dQ_fun = @(ph,th) cos(ph)*sin(4*th)*dp_pk / ...
                  (4*( cos(ph)*sin(2*th) + exp(dp_pk^2/(4*cfg.Delta^2)) ));
swing_centroide_Hz = abs(dQ_fun(0,cfg.theta) - dQ_fun(pi,cfg.theta))*f0_optical;
W_riga_Hz  = f0_optical*cfg.Delta*(cfg.sep_su_Delta/2 + 1);
n_ph_cross = (P_rx_t(idx_cross)*cfg.T_sym/E_fotone)*eta_det;
sigma_centroide_Hz = cfg.F_eccesso*W_riga_Hz/sqrt(n_ph_cross);
snr_eve_centroide  = (swing_centroide_Hz/2)/sigma_centroide_Hz;
% Robustezza del nullo: theta non e' calibrabile esattamente a pi/4
eps_theta = 0.01;   % rad
swing_eps = abs(dQ_fun(0,cfg.theta+eps_theta) - dQ_fun(pi,cfg.theta+eps_theta))*f0_optical;
fprintf('   delta_Q a theta = pi/4 esatto: swing centroide = %.2e Hz (SNR Eve = %.2e)\n', ...
        swing_centroide_Hz, snr_eve_centroide);
fprintf('   robustezza: con errore di %.2f rad su theta -> swing = %.3f Hz contro sigma = %.2f Hz\n', ...
        eps_theta, swing_eps, sigma_centroide_Hz);

%% ------------------------------------------------------------------------
% 4Q.12 (NUOVO, punto 5) ARGOMENTO LPI GEOMETRICO
% -------------------------------------------------------------------------
% Poiche' sotto Kerckhoffs la covertness statistica non regge, l'argomento
% di sicurezza difendibile e' che Eve deve stare FISICAMENTE dentro il
% fascio FSO. Questo e' verificabile e va quantificato, non assunto.
r_spot_cross = theta_div*d_min;               % raggio dello spot a Bob
A_spot       = pi*r_spot_cross^2;
frazione_ang = theta_div^2/4;                 % A_spot / (4 pi d^2)
fprintf('\n--- 4Q.12: LPI geometrico (argomento di sicurezza principale) ---\n');
fprintf('   divergenza = %.1f urad -> raggio dello spot a %.0f km = %.2f m\n', ...
        theta_div*1e6, d_min/1000, r_spot_cross);
fprintf('   frazione di angolo solido illuminata = %.2e\n', frazione_ang);
fprintf('   -> Eve deve occupare quel volume per accedere al fascio; fuori da esso\n');
fprintf('      non riceve NULLA, indipendentemente dal fatto che conosca C.\n');

%% ------------------------------------------------------------------------
% 4Q.13 Verifiche di consistenza fisica (elisione, bin di compressione)
% -------------------------------------------------------------------------
S_bit0 = S_classico + I_env.*cos(Phi_chirp + 0);
S_bit1 = S_classico + I_env.*cos(Phi_chirp + pi);
residuo_Eve = max(abs(0.5*(S_bit0+S_bit1) - S_classico))/max(abs(S_classico));
disp(['Elisione del termine quantistico nella media di Eve: residuo relativo = ', ...
      num2str(residuo_Eve,'%.3e'),'  (deve essere ~0)']);

% Bin di compressione CON il filtro adattato (pesatura per I_env)
tmpl  = I_env.*cos(Phi_chirp);
peso  = I_env.*cfg.chiave;
k_atteso = cfg.N/2 + 1;
Y_ref = fftshift(fft(ifftshift(tmpl.*peso)));
[~, k_misurato] = max(abs(real(Y_ref)));
disp(['Bin di compressione: atteso = ',num2str(k_atteso),', misurato = ',num2str(k_misurato)]);
if abs(k_misurato-k_atteso) > 2
    warning('Il picco non cade dove previsto: verificare C_chirp e la griglia.');
end
cfg.kpk = k_misurato;

% Guadagno del filtro adattato rispetto alla somma non pesata della v1
g_mf   = sqrt(sum(tmpl.^2));                       % ||s||
g_flat = sum(tmpl.*real(cfg.chiave))/sqrt(cfg.N/2);% somma non pesata (statistica v1)
fprintf('GUADAGNO DEL FILTRO ADATTATO rispetto al bin DC non pesato (v1): %+.2f dB\n', ...
        20*log10(g_mf/max(abs(g_flat),eps)));

cfg.Nbit_frame = 1024;    % alzato: serve statistica per misurare capacita' basse
cfg.N_MC       = 16;
cfg.blocco     = 32;
L_frame_vec    = [256 1024];
BER_soglia     = 1e-3;    % SOLO diagnostica, non piu' figura di merito

%% ========================================================================
% FASE 5Q: CAMPAGNA A - WATERFALL E VALIDAZIONE CONTRO Q(d)
% =========================================================================
% Validazione: con il filtro adattato la BER simulata deve coincidere con
% Q(d), d = ||s||/sigma, a meno del jitter di tracking e dell'errore di
% calibrazione. Se non coincide, c'e' ancora perdita nel ricevitore.
disp('--- FASE 5Q: Campagna A - waterfall BER e validazione analitica ---');
sigma_vec = logspace(log10(min(opt.sigma_ep)/6), log10(max(opt.sigma_ep)*6), 16);
nS = numel(sigma_vec); errA = zeros(1,nS); bitA = zeros(1,nS);
BER_teo_A = 0.5*erfc((g_mf./sigma_vec)/sqrt(2));
for is = 1:nS
    for mc = 1:cfg.N_MC
        rs = streamMC(SEED_MASTER,1,is,mc);
        [ne,nb] = simulaFrameQ(cfg, S_classico, I_env, sigma_vec(is), rs);
        errA(is) = errA(is)+ne; bitA(is) = bitA(is)+nb;
    end
    fprintf('   sigma = %.2e -> d = %5.2f | BER sim = %.3e | BER teo = %.3e\n', ...
            sigma_vec(is), g_mf/sigma_vec(is), errA(is)/bitA(is), BER_teo_A(is));
end
BER_A = errA./bitA; [BER_A_lo,BER_A_hi] = wilsonCI(errA,bitA,1.96);
scarto = max(abs(BER_A - BER_teo_A)./max(BER_teo_A,1e-6));
fprintf('   scarto massimo simulazione/teoria: %.1f %% (atteso piccolo se il MF e'' corretto)\n', ...
        100*scarto);

%% ========================================================================
% FASE 6Q: CAMPAGNA B - LUNGO IL PASSAGGIO ORBITALE
% =========================================================================
disp('--- FASE 6Q: Campagna B - BER, capacita'' e accumulo di Eve vs epoca ---');
BER_B = zeros(1,N_epoche); errB = zeros(1,N_epoche); bitB = zeros(1,N_epoche);
sigma_B = zeros(1,N_epoche);
for ie = 1:N_epoche
    ii = idx_ep(ie);
    acc = 0;
    for mc = 1:cfg.N_MC
        rs = streamMC(SEED_MASTER,2,ie,mc);
        th1 = sigma_jit_pat*randn(rs,1,2); th2 = sigma_jit_pat*randn(rs,1,2);
        Lp = exp(-8*(norm(th1)/theta_div)^2)*exp(-8*(norm(th2)/theta_div)^2);
        n_bin     = (P_rx_ideal_t(ii)*Lp*cfg.T_sym/E_fotone)*eta_det/cfg.N;
        sigma_eff = cfg.F_eccesso/sqrt(max(n_bin,eps));
        acc = acc + sigma_eff;
        [ne,nb] = simulaFrameQ(cfg, S_classico, I_env, sigma_eff, rs);
        errB(ie) = errB(ie)+ne; bitB(ie) = bitB(ie)+nb;
    end
    sigma_B(ie) = acc/cfg.N_MC;
    BER_B(ie)   = errB(ie)/bitB(ie);
    fprintf('   t = %6.1f s | d = %7.1f km | sigma = %.2e | d_MF = %5.2f | BER = %.3e\n', ...
            t_ep(ie), d_t(ii)/1000, sigma_B(ie), g_mf/sigma_B(ie), BER_B(ie));
end
[BER_B_lo,BER_B_hi] = wilsonCI(errB,bitB,1.96);

%% ========================================================================
% FASE 6Q-bis: (NUOVO, punto 4) SWEEP SULL'ERRORE DI SOTTRAZIONE
% =========================================================================
% Quanto accuratamente Bob deve conoscere il profilo classico da sottrarre?
% L'errore e' congelato sul frame (e' una calibrazione), quindi entra come
% BIAS sulla statistica di decisione, non come rumore che media a zero.
% Riferimento: se S_cl e' stimato mediando M_cal simboli, eps_sub ~
% sigma/sqrt(M_cal).
disp('--- FASE 6Q-bis: sensibilita'' all''errore di sottrazione del profilo classico ---');
eps_vec = [0 1e-4 1e-3 1e-2 3e-2 1e-1];
sigma_rif = median(sigma_B);
BER_eps = zeros(1,numel(eps_vec));
cfg_tmp = cfg;
for k = 1:numel(eps_vec)
    cfg_tmp.eps_sub = eps_vec(k);
    e = 0; b = 0;
    for mc = 1:cfg.N_MC
        rs = streamMC(SEED_MASTER,3,k,mc);
        [ne,nb] = simulaFrameQ(cfg_tmp, S_classico, I_env, sigma_rif, rs);
        e = e+ne; b = b+nb;
    end
    BER_eps(k) = e/b;
    M_cal_equiv = (sigma_rif/max(eps_vec(k),eps))^2;
    fprintf('   eps_sub = %.1e (equivale a mediare %.1e simboli di calibrazione) -> BER = %.3e\n', ...
            eps_vec(k), M_cal_equiv, BER_eps(k));
end
fprintf('   [NOTA] questo sweep modella un errore per bin INDIPENDENTE, che media su ~%.0f\n', ...
        sqrt(2*pi)*(sigma_env_rad/dw_step));
fprintf('          campioni utili. Un errore STRUTTURATO (deriva di trappola, moto secolare)\n');
fprintf('          proietta molto piu'' efficacemente sul template e va studiato a parte.\n');

% Sensibilita' a OD (analitica, non serve Monte Carlo: d e'' lineare in OD)
disp('--- Sensibilita'' alla profondita'' di imprinting OD ---');
OD_vec = [0.05 0.10 0.30 0.50 OD_max];
for k = 1:numel(OD_vec)
    d_k = (OD_vec(k)/cfg.OD)*(g_mf/sigma_rif);
    fprintf('   OD = %5.3f -> d = %6.3f -> BER grezza = %.3e -> C = %.4f bit/simbolo\n', ...
            OD_vec(k), d_k, 0.5*erfc(d_k/sqrt(2)), 1-entropiaBinaria(0.5*erfc(d_k/sqrt(2))));
end

%% ========================================================================
% FASE 7Q: (NUOVO, punto 2) METRICHE BASATE SULLA CAPACITA'
% =========================================================================
% La figura di merito e' la capacita' per simbolo e il volume integrato,
% non "BER <= 1e-3". Il volume covert e' quello accumulato FINO AL PUNTO in
% cui la statistica di rivelazione di Eve supera la soglia.
disp('--- FASE 7Q: capacita'', throughput e volume covert ---');
n_sym_ep = (dt_ep.*disp_av)/cfg.T_sym;
C_ep_sim = 1 - entropiaBinaria(BER_B);
d_ep_sim = g_mf./sigma_B;
delta_A_ep = d_ep_sim.^2/sqrt(2*cfg.N);
delta_B_ep = d_ep_sim.^2/sqrt(2);

bits_cum   = cumsum(n_sym_ep.*C_ep_sim);
eveA_cum   = sqrt(cumsum(n_sym_ep.*delta_A_ep.^2));
eveB_cum   = sqrt(cumsum(n_sym_ep.*delta_B_ep.^2));

B_covert_A = volumeAllaSoglia(bits_cum, eveA_cum, SNR_soglia_Eve, ...
                              n_sym_ep, C_ep_sim, delta_A_ep);
B_covert_B = volumeAllaSoglia(bits_cum, eveB_cum, SNR_soglia_Eve, ...
                              n_sym_ep, C_ep_sim, delta_B_ep);

R_shan  = R_b_raw*C_ep_sim.*disp_av;
nL = numel(L_frame_vec); R_good = zeros(nL,N_epoche); Vol_good = zeros(1,nL);
for iL = 1:nL
    R_good(iL,:) = R_b_raw*((1-BER_B).^L_frame_vec(iL)).*disp_av;
    Vol_good(iL) = trapz(t_ep,R_good(iL,:));
end
[~,ie_cr] = min(abs(idx_ep-idx_cross));

fprintf('\n=========== SINTESI v2 - CANALE DOPPLER QUANTISTICO ===========\n');
fprintf('Delta_t (ritardo differenziale)         : %8.1f us\n', cfg.Delta_t*1e6);
fprintf('Separazione OTTIMIZZATA (era 7.35)      : %8.2f Delta\n', cfg.sep_su_Delta);
fprintf('Profondita'' di imprinting OD            : %8.3f  (max singolo ione %.3f)\n', cfg.OD, OD_max);
fprintf('B_IF derivata                           : %8.2f MHz\n', cfg.B_IF/1e6);
fprintf('Bit rate lordo                          : %8.3f kbit/s\n', R_b_raw/1e3);
fprintf('Guadagno del filtro adattato (vs v1)    : %+8.2f dB\n', 20*log10(g_mf/max(abs(g_flat),eps)));
fprintf('BER grezza al crossing                  : %8.3f\n', BER_B(ie_cr));
fprintf('Capacita'' al crossing                   : %8.4f bit/simbolo\n', C_ep_sim(ie_cr));
fprintf('Capacita'' mediana sul passaggio         : %8.4f bit/simbolo\n', median(C_ep_sim(disp_av>0)));
fprintf('Throughput informativo di picco (BSC)   : %8.3f kbit/s\n', max(R_shan)/1e3);
fprintf('[diagnostica] epoche con BER <= %.0e     : %8d su %d\n', BER_soglia, sum(BER_B<=BER_soglia), N_epoche);
fprintf('---------------------------------------------------------------\n');
fprintf('VOLUME se Eve non ascolta                : %8.3e bit\n', bits_cum(end));
fprintf('VOLUME COVERT contro Eve-A (energia)     : %8.3e bit\n', B_covert_A);
fprintf('   statistica di Eve-A a fine passaggio  : %8.2f  (soglia %.1f)\n', eveA_cum(end), SNR_soglia_Eve);
fprintf('VOLUME COVERT contro Eve-B (chiave nota) : %8.3e bit\n', B_covert_B);
fprintf('   statistica di Eve-B a fine passaggio  : %8.2e (soglia %.1f)\n', eveB_cum(end), SNR_soglia_Eve);
fprintf('---------------------------------------------------------------\n');
for iL = 1:nL
    fprintf('[diagnostica] goodput non codificato %4d bit: %8.3e kbit/s\n', ...
            L_frame_vec(iL), max(R_good(iL,:))/1e3);
end
fprintf('CONCLUSIONE: il volume difendibile e'' quello contro Eve-A, e vale solo\n');
fprintf('se la covertness poggia sull''argomento GEOMETRICO di 4Q.12. Sotto\n');
fprintf('ipotesi di Kerckhoffs (Eve-B) il canale non e'' covert.\n');
fprintf('===============================================================\n\n');

%% ========================================================================
% FASE 8Q: FIGURE
% =========================================================================
disp('--- FASE 8Q: Figure ---');
f_MHz = cfg.dw/(2*pi)/1e6;

figure('Name','FigQ1 - Chirp quantistico e filtro adattato','Color','w','Position',[80 60 1150 700]);

subplot(2,2,1);
plot(f_MHz,S_classico,'LineWidth',1.4); hold on;
plot(f_MHz,S_bit0,'LineWidth',0.9);
xlabel('Detuning [MHz]'); ylabel('Modulazione frazionaria del fascio di sonda');
title(sprintf('(a) Spettro classico e con perturbazione quantistica (OD = %.2f)',cfg.OD));
legend({'S_{classico}','S_{quant} (bit 0)'},'Location','best','FontSize',8);
grid on; xlim([-1 1]*larghezza_Hz/1e6);

subplot(2,2,2);
plot(f_MHz, tmpl, 'LineWidth',1.0); hold on;
plot(f_MHz, I_env,'--','LineWidth',1.3); plot(f_MHz,-I_env,'--','LineWidth',1.3);
xlabel('Detuning [MHz]'); ylabel('Termine di interferenza');
title(sprintf('(b) Perturbazione CHIRP, \\Phi = -C\\Delta\\omega^2 da \\Deltat = %.0f \\mus', cfg.Delta_t*1e6));
legend({'I_{env}cos(\Phi)','inviluppo \pm I_{env}'},'Location','best','FontSize',8);
grid on; xlim([-1 1]*2*sigma_env_rad/(2*pi)/1e6);

subplot(2,2,3);
uax = (1:cfg.N)-cfg.N/2-1;
Y0_mf = real(fftshift(fft(ifftshift((I_env.*cos(Phi_chirp+0)).*peso))));
Y1_mf = real(fftshift(fft(ifftshift((I_env.*cos(Phi_chirp+pi)).*peso))));
Y0_v1 = real(fftshift(fft(ifftshift((I_env.*cos(Phi_chirp+0)).*cfg.chiave))));
plot(uax,Y0_mf/max(abs(Y0_mf)),'LineWidth',1.4); hold on;
plot(uax,Y1_mf/max(abs(Y0_mf)),'LineWidth',1.4);
plot(uax,Y0_v1/max(abs(Y0_v1)),':','LineWidth',1.1);
yline(0,'k--'); xlabel('bin'); ylabel('Ampiezza compressa (norm.)');
title('(c) Filtro adattato: picco piu'' largo e piu'' tollerante al jitter');
legend({'bit 0 (MF)','bit 1 (MF)','bit 0 (somma non pesata, v1)'},'Location','best','FontSize',7);
grid on; xlim([-80 80]);

subplot(2,2,4);
bp = BER_A; bp(bp==0) = NaN;
errorbar(sigma_vec,bp,bp-BER_A_lo,BER_A_hi-bp,'o','LineWidth',1.3,'MarkerSize',5); hold on;
plot(sigma_vec,BER_teo_A,'-','LineWidth',1.2);
set(gca,'XScale','log','YScale','log');
xline(min(opt.sigma_ep),'g--','\sigma crossing','LineWidth',1.1,'LabelOrientation','horizontal');
xline(max(opt.sigma_ep),'m--','\sigma bordo','LineWidth',1.1,'LabelOrientation','horizontal');
xlabel('\sigma rumore relativo'); ylabel('BER grezza');
title('(d) Waterfall: simulazione contro Q(d) analitica');
legend({'Monte Carlo','Q(\Vert s\Vert/\sigma)'},'Location','best','FontSize',8);
grid on;

figure('Name','FigQ2 - Passaggio orbitale e accumulo di Eve','Color','w','Position',[110 40 1150 430]);
subplot(1,3,1);
yyaxis left; plot(t_ep,BER_B,'o-','LineWidth',1.3,'MarkerSize',4); ylabel('BER grezza'); ylim([0 0.55]);
yyaxis right; plot(t_ep,C_ep_sim,'s--','LineWidth',1.1,'MarkerSize',4);
ylabel('Capacita'' [bit/simbolo]');
xline(t_sec(idx_cross),'k:','Crossing'); xlabel('Tempo [s]');
title('(a) BER grezza e capacita'' lungo il passaggio'); grid on;

subplot(1,3,2);
plot(t_ep,bits_cum/1e3,'LineWidth',1.5); hold on;
yline(B_covert_A/1e3,'g--','Volume covert vs Eve-A','LineWidth',1.2);
yline(B_covert_B/1e3,'r--','Volume covert vs Eve-B','LineWidth',1.2);
xlabel('Tempo [s]'); ylabel('Volume cumulato [kbit]');
title('(b) Volume informativo cumulato'); grid on;

subplot(1,3,3);
semilogy(t_ep,max(eveA_cum,1e-12),'LineWidth',1.4); hold on;
semilogy(t_ep,max(eveB_cum,1e-12),'LineWidth',1.4);
yline(SNR_soglia_Eve,'k--','soglia','LineWidth',1.2);
xlabel('Tempo [s]'); ylabel('Statistica di rivelazione accumulata');
title('(c) Accumulo di Eve (legge della radice quadrata)');
legend({'Eve-A (ignara)','Eve-B (chiave nota)'},'Location','best','FontSize',8); grid on;

figure('Name','FigQ3 - Progetto: ottimo di sep e sensibilita''','Color','w','Position',[140 20 1150 420]);
subplot(1,3,1);
semilogy(sep_vec,max(B_noEve,1e-3),'--','LineWidth',1.2); hold on;
semilogy(sep_vec,max(B_en,1e-3),'LineWidth',1.6);
semilogy(sep_vec,max(B_mf,1e-3),'LineWidth',1.2);
xline(cfg.sep_su_Delta,'k:','ottimo','LineWidth',1.2);
xline(7.35,'r:','v1 (manuale)','LineWidth',1.1);
xlabel('sep / \Delta'); ylabel('Volume per passaggio [bit]');
title('(a) Ottimizzazione del punto operativo');
legend({'senza Eve','covert vs Eve-A','covert vs Eve-B'},'Location','best','FontSize',7); grid on;

subplot(1,3,2);
semilogx(max(eps_vec,1e-5),BER_eps,'o-','LineWidth',1.4);
xlabel('\epsilon_{sub} (errore relativo sul profilo classico)'); ylabel('BER grezza');
title('(b) Sensibilita'' alla sottrazione di S_{classico}'); grid on;

subplot(1,3,3);
OD_fine = linspace(0.01,max(OD_max,0.6),200);
d_fine  = (OD_fine/cfg.OD)*(g_mf/sigma_rif);
plot(OD_fine, 1-entropiaBinaria(0.5*erfc(d_fine/sqrt(2))),'LineWidth',1.5); hold on;
xline(cfg.OD,'k:','progetto','LineWidth',1.2);
xline(OD_max,'r:','limite ione','LineWidth',1.2);
xlabel('OD (profondita'' di imprinting)'); ylabel('Capacita'' [bit/simbolo]');
title('(c) Sensibilita'' a OD'); grid on;

disp('--- SIMULAZIONE COMPLETATA ---');

%% ========================================================================
% FUNZIONI LOCALI
% =========================================================================

function out = puntoOperativo(sep, par)
% Punto di progetto per una data separazione dei pacchetti. Restituisce
% spettri, deflessioni di Bob e di entrambe le Eve, e i volumi risultanti.
% Tutto e' derivato: nulla e' scelto a mano dentro questa funzione.
    Delta = par.Delta; N = par.N;
    larghezza_Hz = par.f0*(sep+6)*Delta;
    B_IF  = 4*larghezza_Hz;
    T_sym = N/B_IF;
    dw = linspace(-1,1,N).'*(2*pi*B_IF/2);
    p  = dw/par.Omega_rad;
    p1 = -sep*Delta/2;  p2 = +sep*Delta/2;
    A1 = exp(-((p-p1).^2)./(2*Delta^2))./(pi^0.25*sqrt(Delta));
    A2 = exp(-((p-p2).^2)./(2*Delta^2))./(pi^0.25*sqrt(Delta));
    w  = 1 + 3*p;                                   % peso di Eq.(26)
    S_cl  = (cos(par.theta)^2*A1.^2 + sin(par.theta)^2*A2.^2).*w;
    I_env = (2*cos(par.theta)*sin(par.theta)*A1.*A2).*w;
    nrm = sqrt(mean(S_cl.^2));
    % OD: profondita' di modulazione frazionaria del fascio di sonda.
    S_cl  = par.OD*S_cl/nrm;
    I_env = par.OD*I_env/nrm;
    Phi   = -par.C_chirp*dw.^2;

    E_s2 = sum((I_env.*cos(Phi)).^2);      % ||s||^2 con sigma = 1

    n_sym    = par.dt_ep(:).'/T_sym;
    n_bin_ep = (par.P_ep(:).'*T_sym/par.E_fotone)*par.eta_det/N;
    sigma_ep = par.F_eccesso./sqrt(max(n_bin_ep,eps));
    d_ep     = sqrt(E_s2)./sigma_ep;                % deflessione del MF
    BER_ep   = 0.5*erfc(d_ep/sqrt(2));
    C_ep     = 1 - entropiaBinaria(BER_ep);
    dl_A     = d_ep.^2/sqrt(2*N);                   % Eve ignara (energia)
    dl_B     = d_ep.^2/sqrt(2);                     % Eve con chiave (matched)

    bits_cum = cumsum(n_sym.*C_ep);
    eveA     = sqrt(cumsum(n_sym.*dl_A.^2));
    eveB     = sqrt(cumsum(n_sym.*dl_B.^2));

    out.sep=sep; out.B_IF=B_IF; out.T_sym=T_sym; out.dw=dw;
    out.S_cl=S_cl; out.I_env=I_env; out.Phi=Phi;
    out.sigma_ep=sigma_ep; out.d_ep=d_ep; out.BER_ep=BER_ep; out.C_ep=C_ep;
    out.M_pass=sum(n_sym);
    out.bits_no_eve = bits_cum(end);
    out.bits_covert_en = volumeAllaSoglia(bits_cum, eveA, par.SNR_thr, n_sym, C_ep, dl_A);
    out.bits_covert_mf = volumeAllaSoglia(bits_cum, eveB, par.SNR_thr, n_sym, C_ep, dl_B);
    out.nyq = 2*par.C_chirp*(2*(Delta/sqrt(2))*par.Omega_rad)*(dw(2)-dw(1));
end

function B = volumeAllaSoglia(bits_cum, eve, thr, n_sym, C_ep, dl)
% Volume trasmissibile prima che la statistica accumulata di Eve superi la
% soglia. Interpolazione lineare in energia dentro l'epoca di superamento.
    k = find(eve > thr, 1);
    if isempty(k)
        B = bits_cum(end);
        return;
    end
    if k == 1
        prevE2 = 0; prevB = 0;
    else
        prevE2 = eve(k-1)^2; prevB = bits_cum(k-1);
    end
    quota = n_sym(k)*dl(k)^2;
    if quota <= 0
        B = prevB;
    else
        frac = min(max((thr^2 - prevE2)/quota, 0), 1);
        B = prevB + frac*n_sym(k)*C_ep(k);
    end
end

function [nerr,nbit] = simulaFrameQ(cfg, S_cl, I_env, sigma_noise, rs)
% Trasmissione di un frame. Il bit e' la fase relativa phi in {0,pi}: per
% Eq.(32) il termine di interferenza va come cos(phi), quindi la codifica e'
% ANTIPODALE per costruzione fisica.
%
% RICEVITORE (correzione principale rispetto alla v1): Bob pesa per
% l'inviluppo noto I_env PRIMA del de-chirp. Il bin DC di
% fft(Rx.*I_env.*chiave) e' allora esattamente sum(Rx.*I_env.*cos(Phi)),
% cioe' il filtro adattato al template reale. La v1 usava fft(Rx.*chiave),
% cioe' la somma NON pesata su tutti gli N campioni: raccoglieva rumore da
% N campioni per estrarre segnale da ~sqrt(2*pi)*sigma_bin campioni utili.
%
% Il segnale NON viene scalato con la distanza: il contrasto e' una
% modulazione frazionaria (un rapporto non si attenua propagando). La
% distanza entra tramite sigma_noise, derivato dal conteggio fotonico.
    M = cfg.Nbit_frame; N = cfg.N;
    bits = randi(rs,[0 1],1,M);
    Phi  = -cfg.C_chirp*cfg.dw.^2;
    tmpl = I_env.*cos(Phi);

    % Errore di calibrazione del profilo classico: CONGELATO sul frame,
    % quindi agisce come bias sulla decisione, non come rumore che media.
    S_hat = S_cl.*(1 + cfg.eps_sub*randn(rs,N,1));

    % Errore residuo sul coefficiente di chirp: raccoglie incertezza su
    % Delta_t e residuo di RATE Doppler non compensato.
    C_bob  = cfg.C_chirp*(1 + cfg.eps_C*randn(rs,1,1));
    peso   = I_env.*exp(+1i*C_bob*cfg.dw.^2);

    nerr = 0;
    for i1 = 1:cfg.blocco:M
        i2 = min(i1+cfg.blocco-1,M); b = bits(i1:i2); nb = numel(b);
        segno = 1-2*b;                              % bit0 -> +1, bit1 -> -1
        Sig = S_cl + tmpl.*segno;
        Rx  = Sig + sigma_noise*randn(rs,N,nb);
        Rx  = Rx - S_hat;                           % sottrazione imperfetta
        % De-chirp + filtro adattato. L'ifftshift e' ESSENZIALE: il segnale
        % e' centrato a dw = 0, cioe' al CENTRO dell'array, non all'indice 1.
        Y   = fftshift(fft(ifftshift(Rx.*peso,1),[],1),1);
        % Residuo di OFFSET Doppler: sposta il bin di lettura di
        % k = 4*C*df*dw_max bin (derivato, non calibrato a mano).
        jit = round(cfg.sigma_bin_track*randn(rs,1,nb));
        kuse= min(max(cfg.kpk + jit,1),N);
        stat= real(Y(sub2ind(size(Y),kuse,1:nb)));
        nerr= nerr + sum(double(stat<0) ~= b);
    end
    nbit = M;
end

function rs = streamMC(seed,idc,idp,idt)
    rs = RandStream('mrg32k3a','Seed',seed);
    rs.Substream = 1 + idt + 1000*idp + 1000000*idc;
end

function [lo,hi] = wilsonCI(k,n,z)
    p = k./n; d = 1+z^2./n; c = p+z^2./(2*n);
    s = z*sqrt(p.*(1-p)./n + z^2./(4*n.^2));
    lo = max((c-s)./d,0); hi = min((c+s)./d,1);
end

function H = entropiaBinaria(p)
    p = min(max(p,0),1); H = zeros(size(p)); m = (p>0)&(p<1);
    H(m) = -p(m).*log2(p(m)) - (1-p(m)).*log2(1-p(m));
end

% =========================================================================
% NOTE TECNICHE - VERSIONE 2
% =========================================================================
% NOTA V2.1 - PERCHE' IL FILTRO ADATTATO CAMBIA TUTTO.
%   Dopo il de-chirp il segnale utile e' (1/2) I_env(w) e^{i phi}: vive su
%   ~sqrt(2 pi) sigma_bin campioni (poche centinaia su N = 16384). Il bin DC
%   della FFT non pesata somma TUTTI gli N campioni, quindi raccoglie rumore
%   sigma*sqrt(N/2) contro un segnale proporzionale a sum(I_env). Il
%   rapporto fra le due deflessioni vale
%       SNR_MF / SNR_flat = ||g|| * sqrt(N/2) / sum(g)
%   che per una gaussiana di sigma_b bin su N campioni da' circa
%       sqrt( N / (2 sqrt(pi) sigma_b) )
%   cioe' ~10 dB per i parametri di questo scenario. Il codice lo verifica a
%   runtime (variabile g_mf/g_flat) e la Campagna A lo valida confrontando la
%   BER simulata con Q(||s||/sigma).
%
% NOTA V2.2 - PERCHE' LA SOGLIA BER = 1e-3 ERA LA METRICA SBAGLIATA.
%   Un canale covert opera per costruzione a SNR bassa: il punto operativo
%   ottimo sta a BER grezza dell'ordine di 0.2-0.3, non 1e-3, e l'affidabilita'
%   si recupera con codifica a rate basso. Valutare (1-BER)^L su frame non
%   codificati e' quindi una metrica priva di significato in questo regime, e
%   produceva goodput identicamente nullo. In v2 la figura di merito e' la
%   capacita' 1 - H2(BER) integrata sui simboli del passaggio; BER_soglia
%   resta solo come diagnostica.
%   AVVERTENZA STATISTICA: per misurare capacita' dell'ordine di 1e-3
%   bit/simbolo servirebbero ~1e6 bit per punto. Con Nbit_frame*N_MC bit per
%   epoca l'intervallo di Wilson su BER ha semiampiezza ~1.96*sqrt(0.25/n):
%   se il punto operativo finisce a BER ~0.49, la capacita' misurata NON e'
%   distinguibile da zero e va riportata come limite superiore, non come
%   risultato. Il sweep di progetto serve anche a evitare quel regime.
%
% NOTA V2.3 - STRUTTURA DEL PROBLEMA DI PROGETTO.
%   Bob: BER = Q(d), d = ||s||/sigma, C(d) = 1 - H2(Q(d)).
%   Eve-A: deflessione per simbolo d^2/sqrt(2N); su M simboli sqrt(M) d^2/sqrt(2N).
%   Eve-B: deflessione per simbolo d^2/sqrt(2);  su M simboli sqrt(M) d^2/sqrt(2).
%   Massimizzando M*C(d) sotto sqrt(M) d^2/sqrt(2N) <= thr e usando
%   C(d) ~ d^2/(pi ln2) a piccolo d, si trova che il volume covert va come
%   1/d^2 finche' Eve resta il vincolo attivo: conviene quindi d PICCOLO e
%   M grande, e l'ottimo cade esattamente dove Eve raggiunge la soglia a fine
%   passaggio, cioe' d_opt = (2 N thr^2 / M_pass)^(1/4). Il volume covert
%   corrispondente scala come sqrt(2 N M_pass), cioe' come la radice del
%   numero totale di campioni: e' la legge della radice quadrata delle
%   comunicazioni covert, qui ottenuta come risultato e non assunta.
%
% NOTA V2.4 - OD E CALIBRAZIONE: DUE PARAMETRI CHE PRIMA ERANO NASCOSTI.
%   (a) OD. Normalizzare S_cl a RMS 1 e derivare sigma dal conteggio fotonico
%       equivale ad assumere modulazione al 100% del fascio di sonda. Per un
%       singolo ione la modulazione frazionaria e' limitata da
%       sigma_abs/A_fascio, con sigma_abs = 3 lambda^2/(2 pi). d e' LINEARE in
%       OD, quindi la capacita' e' molto sensibile alla qualita' del fuoco
%       sullo ione: e' il parametro sperimentale piu' critico dopo Delta_t.
%   (b) eps_sub. La v1 sottraeva S_cl esattamente noto. In v2 l'errore e'
%       congelato sul frame (e' una calibrazione) e produce un BIAS. Con
%       errore per bin indipendente il bias media su ~sqrt(2 pi) sigma_bin
%       campioni e il margine risulta ampio; un errore STRUTTURATO (deriva
%       della trappola, moto secolare residuo, deriva del LO) proietta molto
%       piu' efficacemente sul template e richiede uno studio separato. Va
%       dichiarato come limite del modello.
%
% NOTA V2.5 - IL RISULTATO ONESTO SULLA COVERTNESS.
%   delta_B / delta_A = sqrt(N) esattamente. Bob e Eve-B usano lo STESSO
%   filtro: l'unica cosa che Bob ha in piu' e' la conoscenza di phi, che
%   serve a leggere il bit, non a rilevare la trasmissione. Quindi ogni volta
%   che Bob decodifica in modo affidabile, Eve-B rileva. La covertness di
%   questo schema NON puo' poggiare sull'ignoranza di C da parte di Eve: C
%   dipende solo da Delta_t, m e Omega, e comunque e' stimabile dallo spettro
%   senza conoscere il bit (autocorrelazione, distribuzione di Wigner), con
%   uno spazio di ricerca monodimensionale.
%   L'argomento difendibile e' quello geometrico di 4Q.12: con divergenza di
%   8 urad la frazione di angolo solido illuminata e' theta_div^2/4 ~ 1.6e-11.
%   Eve deve stare dentro il fascio. Questo va presentato come il meccanismo
%   di sicurezza primario, e la statistica di Eve-A come difesa di secondo
%   livello contro un intercettatore che si trovi dentro il fascio ma non
%   sappia che cosa cercare.
%
% NOTA V2.6 - VINCOLI SU Delta_t, ORA TRE.
%   (a) NYQUIST sulla griglia spettrale (verificato a runtime).
%   (b) COERENZA della sovrapposizione (non verificabile dal simulatore).
%   (c) LARGHEZZA NATURALE: il passo delle frange del chirp al bordo
%       dell'inviluppo vale 2 pi/(2 C dw_bordo); se scende sotto Gamma0 le
%       frange vengono lavate dalla lorentziana della transizione. Calcolato
%       a runtime in 4Q.8.
%   In piu', il meccanismo stesso di Delta_t richiede evoluzione LIBERA:
%   in trappola armonica lo stato coerente ruota nello spazio delle fasi e
%   non accumula la fase quadratica in p. Il codice stampa l'escursione
%   spaziale dei due rami durante Delta_t come verifica di fattibilita'.
% =========================================================================