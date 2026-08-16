% =========================================================================
% SIMULATORE DI THROUGHPUT - LINK INTERSATELLITARE QUANTUM-COVERT (LEO-LEO)
% Estensione multi-bit + campagna Monte Carlo (BER, bit rate, volume dati)
% Modello: Ione 27Al+ | Doppler Dinamico & FrFT De-chirping
%
% NOVITA' RISPETTO ALLA VERSIONE A BIT SINGOLO
%  - Trasmissione di frame multi-bit (BPSK sul chirp: fase 0 / pi)
%  - Campagna A: BER vs livello di covertness (SNR target) al crossing
%  - Campagna B: BER / throughput vs epoca orbitale (link budget dinamico)
%  - Metriche: BER + intervalli di confidenza, bit rate, goodput, capacita'
%              di canale binario, volume dati per passaggio
%  - Riproducibilita' totale: seed master unico + substream indipendenti
%  - Rivelatore corretto: decisione sul bin di picco NOTO (non su argmax),
%    vedi NOTA TECNICA 1 in fondo al file
% =========================================================================

clc; clear; close all;

% -------------------------------------------------------------------------
% SEED MASTER UNICO PER TUTTA LA CAMPAGNA
% -------------------------------------------------------------------------
% Ogni realizzazione Monte Carlo usa un SUBSTREAM distinto dello STESSO
% generatore inizializzato con SEED_MASTER. Questo garantisce:
%   (1) riproducibilita' bit-per-bit dell'intera campagna fra esecuzioni;
%   (2) indipendenza statistica fra le realizzazioni (condizione necessaria
%       perche' la media Monte Carlo abbia senso).
% Riusare letteralmente lo stesso stato del generatore per ogni trial
% produrrebbe N copie identiche e varianza nulla: non sarebbe Monte Carlo.
SEED_MASTER = 42;
rng(SEED_MASTER, 'twister');   % stream di default (jitter PAT di Fase 3)

%% ========================================================================
% SEZIONE 0: COSTANTI FISICHE E PARAMETRI GLOBALI
% =========================================================================
c_light     = 299792458;          % Velocita' della luce nel vuoto (m/s)
lambda_0    = 267e-9;             % Transizione UV ione 27Al+ (m)
f0_optical  = c_light / lambda_0; % Frequenza portante ottica (Hz)
h_planck    = 6.626e-34;          % Costante di Planck (J*s)
E_fotone    = h_planck * f0_optical; % Energia del singolo fotone (J)
Re          = 6378137;            % Raggio equatoriale Terra, WGS84 (m)
mu_earth    = 3.986004418e14;     % Parametro gravitazionale terrestre (m^3/s^2)

disp('=== SEZIONE 0: Costanti caricate ===');
disp(['Frequenza portante ottica f0: ', num2str(f0_optical/1e12,'%.3f'), ' THz']);

%% ========================================================================
% FASE 1: GEOMETRIA ORBITALE CON EFFETTO J2 E TARGETING DEL CROSSING
% =========================================================================
disp('--- FASE 1: Ricerca geometria di crossing con perturbazioni J2 ---');

startTime  = datetime(2026,8,14,10,0,0);
stopTime   = startTime + minutes(15);
sampleTime = 0.1;                       % Risoluzione temporale (s)
T_window   = seconds(stopTime - startTime);
t_target_middle = T_window / 2;         % Centro finestra (450 s)

altitude_orbit = 800e3;                 % Quota orbitale (m)
sma            = Re + altitude_orbit;   % Semiasse maggiore (m)
ecc            = 0.0001;                % Quasi-circolare (compatibile con SGP4/J2)
incl_deg       = 60;                    % Inclinazione (deg)
argPeri_deg    = 0;                     % Argomento del perigeo

J2 = 1.08263e-3;                        % Armonica zonale terrestre

mu_n      = sqrt(mu_earth / sma^3);     % Moto medio (rad/s)
p_semi    = sma * (1 - ecc^2);
raan_rate = -1.5 * mu_n * J2 * (Re / p_semi)^2 * cosd(incl_deg); % rad/s
disp(['Precessione nodale J2 (dRAAN/dt): ', num2str(rad2deg(raan_rate)*86400,'%.4f'), ' deg/giorno']);

theta_cross_deg = 60;
cos_dRAAN = (cosd(theta_cross_deg) - cosd(incl_deg)^2) / sind(incl_deg)^2;
dRAAN_deg = acosd(cos_dRAAN);
raan_Alice = 0;
raan_Bob   = dRAAN_deg;

T_orbit   = 2*pi / mu_n;
incl_rad  = deg2rad(incl_deg);
raanA_rad = deg2rad(raan_Alice);
raanB_rad = deg2rad(raan_Bob);

posCircularJ2 = @(u0, raan0, t) sma * [ ...
    cos(raan0 + raan_rate*t).*cos(u0 + mu_n*t) - sin(raan0 + raan_rate*t).*sin(u0 + mu_n*t).*cos(incl_rad); ...
    sin(raan0 + raan_rate*t).*cos(u0 + mu_n*t) + cos(raan0 + raan_rate*t).*sin(u0 + mu_n*t).*cos(incl_rad); ...
    sin(u0 + mu_n*t).*sin(incl_rad) ];

target_distance = 300e3;                % 300 km al fly-by

t_full_coarse = 0:2:T_orbit;
posA_full_coarse = posCircularJ2(0, raanA_rad, t_full_coarse);
u0_candidates = deg2rad(0:0.5:359.5);

best_obj = Inf; best_u0 = 0; best_t_period = 0;
for k = 1:numel(u0_candidates)
    posB_k = posCircularJ2(u0_candidates(k), raanB_rad, t_full_coarse);
    d_k = vecnorm(posB_k - posA_full_coarse, 2, 1);
    [dmin_k, idx_k] = min(d_k);
    obj_k = abs(dmin_k - target_distance);
    if obj_k < best_obj
        best_obj = obj_k; best_u0 = u0_candidates(k); best_t_period = t_full_coarse(idx_k);
    end
end

t_fine = max(0,best_t_period-10):0.02:min(T_orbit, best_t_period+10);
posA_fine = posCircularJ2(0, raanA_rad, t_fine);
u0_fine = best_u0 + deg2rad(-1:0.002:1);

best_obj_f = Inf; best_u0_f = best_u0; best_t_f = best_t_period; best_dmin_f = Inf;
for k = 1:numel(u0_fine)
    posB_k = posCircularJ2(u0_fine(k), raanB_rad, t_fine);
    d_k = vecnorm(posB_k - posA_fine, 2, 1);
    [dmin_k, idx_k] = min(d_k);
    obj_k = abs(dmin_k - target_distance);
    if obj_k < best_obj_f
        best_obj_f = obj_k; best_u0_f = u0_fine(k); best_t_f = t_fine(idx_k); best_dmin_f = dmin_k;
    end
end

Delta_phase        = best_u0_f;
u_Alice_at_min_ref = mu_n * best_t_f;

trueAnom_Alice = mod(rad2deg(u_Alice_at_min_ref) - rad2deg(mu_n*t_target_middle), 360);
trueAnom_Bob   = mod(trueAnom_Alice + rad2deg(Delta_phase), 360);

disp(['Anomalia vera iniziale: Alice = ', num2str(trueAnom_Alice,'%.3f'), ...
      '°, Bob = ', num2str(trueAnom_Bob,'%.3f'), '°']);

%% ========================================================================
% FASE 2: PROPAGAZIONE SGP4/J2, GEOMETRIA ISL E CONGIUNZIONE SOLARE
% =========================================================================
disp('--- FASE 2: Propagazione SGP4 e Analisi di Congiunzione Solare ---');

sc = satelliteScenario(startTime, stopTime, sampleTime);

satAlice = satellite(sc, sma, ecc, incl_deg, raan_Alice, argPeri_deg, ...
    trueAnom_Alice, "OrbitPropagator", "sgp4", "Name", "Alice");
satBob   = satellite(sc, sma, ecc, incl_deg, raan_Bob, argPeri_deg, ...
    trueAnom_Bob, "OrbitPropagator", "sgp4", "Name", "Bob");

[posAlice, velAlice, tSamples] = states(satAlice, "CoordinateFrame", "inertial");
[posBob,   velBob,   ~]        = states(satBob,   "CoordinateFrame", "inertial");

t_sec = seconds(tSamples - tSamples(1));
N_t   = numel(t_sec);

deltaPos = posBob - posAlice;
d_t      = vecnorm(deltaPos, 2, 1);
tau_t    = d_t / c_light;

uLOS_AliceToBob = deltaPos ./ d_t;
uLOS_BobToAlice = -uLOS_AliceToBob;

deltaVel = velBob - velAlice;
v_rel_t  = sum(deltaVel .* uLOS_AliceToBob, 1);

f_D_t    = -(v_rel_t / c_light) * f0_optical;
dfD_dt_t = gradient(f_D_t, sampleTime);

% Occlusione terrestre
h_atmosphere_limit = 80e3;
r_min_los = zeros(1, N_t);
for i = 1:N_t
    r_A = posAlice(:,i);
    r_B = posBob(:,i);
    r_min_los(i) = norm(cross(r_A, r_B)) / norm(r_B - r_A);
end
los_clearance_margin = r_min_los - (Re + h_atmosphere_limit);
los_blocked = any(los_clearance_margin < 0);

% Posizione del Sole e angolo di esclusione
AU = 1.495978707e11;
eps_ecl = deg2rad(23.439291);
n_days = days(tSamples - datetime(2000,1,1,12,0,0, 'TimeZone', 'UTC'));
lambda_sun = deg2rad(mod(280.460 + 0.9856474 * n_days, 360));

posSun = AU * [cos(lambda_sun); ...
               sin(lambda_sun)*cos(eps_ecl); ...
               sin(lambda_sun)*sin(eps_ecl)];

sun_dir_Bob = posSun - posBob;
uSun_Bob    = sun_dir_Bob ./ vecnorm(sun_dir_Bob, 2, 1);
cos_theta_sun = sum(uLOS_BobToAlice .* uSun_Bob, 1);
theta_sun_rx_deg = rad2deg(acos(min(max(cos_theta_sun, -1), 1)));

theta_sun_crit = 30;
sun_blindness_flag = theta_sun_rx_deg < theta_sun_crit;

[d_min, idx_cross] = min(d_t);
disp(['Distanza minima di fly-by: ', num2str(d_min/1000,'%.2f'), ' km a t = ', num2str(t_sec(idx_cross),'%.1f'), ' s']);
disp(['Doppler massimo: ', num2str(max(abs(f_D_t))/1e9,'%.3f'), ' GHz, Doppler-rate: ', num2str(abs(dfD_dt_t(idx_cross))/1e6,'%.3f'), ' MHz/s']);

if los_blocked
    warning('ATTENZIONE: La linea di vista (LOS) è parzialmente occlusa dalla Terra/Atmosfera!');
else
    disp(['Line-of-Sight (LOS) libera: clearance minima sopra atmosfera = ', ...
          num2str(min(los_clearance_margin)/1000,'%.1f'), ' km']);
end

disp(['Angolo di esclusione solare (Sun-Rx angle) al crossing: ', num2str(theta_sun_rx_deg(idx_cross),'%.2f'), '°']);
if any(sun_blindness_flag)
    warning('ATTENZIONE: Il ricevitore entra nella zona di abbagliamento solare (angolo < 30°)!');
else
    disp('Condizione solare OK: ricevitore non abbagliato per tutta la finestra.');
end

% Maschera di disponibilita' del link (usata nel calcolo del throughput)
link_available = (los_clearance_margin >= 0) & (~sun_blindness_flag);

%% ========================================================================
% FASE 3: MODELLO DI CANALE OTTICO SATELLITARE STANDARDIZZATO (FSO / PAT)
% =========================================================================
disp('--- FASE 3: Standardizzazione Canale Ottico (FSO Link Budget & PAT) ---');

D_tx        = 0.20;         % Apertura Tx Alice (m)
D_rx        = 0.30;         % Apertura Rx Bob (m)
P_tx_laser  = 100e-3;       % Potenza media trasmessa (W)

eta_tx      = 0.85;
eta_rx      = 0.80;
eta_det     = 0.65;
DCR         = 50;           % Dark Count Rate (counts/s)

sigma_jitter_pat = 1.5e-6;  % Jitter PAT 1-sigma per asse (rad)
theta_div_beam   = 8.0e-6;  % Divergenza nominale del fascio (rad)

% NOTA TECNICA 2 (vedi in fondo): G = (pi*D/lambda)^2 presuppone un fascio
% al limite di diffrazione (lambda/D ~ 1.34 urad), incoerente con la
% divergenza dichiarata di 8 urad. Il flag seguente permette di usare il
% guadagno coerente con la divergenza effettiva.
USA_GUADAGNO_COERENTE = false;

if USA_GUADAGNO_COERENTE
    G_tx = 32 / theta_div_beam^2;       % Guadagno per fascio gaussiano reale
else
    G_tx = (pi * D_tx / lambda_0)^2;    % Formula originale (diffraction-limited)
end
G_rx = (pi * D_rx / lambda_0)^2;

FSPL_t = (4 * pi * d_t / lambda_0).^2;

% Jitter di puntamento (Rayleigh) - realizzazione nominale per i grafici
theta_az_tx = sigma_jitter_pat * randn(1, N_t);
theta_el_tx = sigma_jitter_pat * randn(1, N_t);
theta_err_tx_radiale = sqrt(theta_az_tx.^2 + theta_el_tx.^2);

theta_az_rx = sigma_jitter_pat * randn(1, N_t);
theta_el_rx = sigma_jitter_pat * randn(1, N_t);
theta_err_rx_radiale = sqrt(theta_az_rx.^2 + theta_el_rx.^2);

L_point_tx = exp(-8 * (theta_err_tx_radiale ./ theta_div_beam).^2);
L_point_rx = exp(-8 * (theta_err_rx_radiale ./ theta_div_beam).^2);
L_point_t  = L_point_tx .* L_point_rx;

P_rx_ideal_t = P_tx_laser * eta_tx * eta_rx * (G_tx * G_rx) ./ FSPL_t;
P_rx_t       = P_rx_ideal_t .* L_point_t;

P_rx_crossing_dBm = 10*log10(P_rx_t(idx_cross)*1e3);
L_point_avg_dB    = 10*log10(mean(L_point_t));

disp(['Potenza Rx al crossing (con jitter ed efficienze): ', num2str(P_rx_crossing_dBm,'%.2f'), ' dBm']);
disp(['Perdita media di puntamento (PAT Loss): ', num2str(L_point_avg_dB,'%.2f'), ' dB']);

%% ========================================================================
% FASE 4: CONFIGURAZIONE DELLA CAMPAGNA MULTI-BIT E DELLA FORMA D'ONDA
% =========================================================================
disp('--- FASE 4: Configurazione campagna multi-bit ---');

cfg = struct();

% --- Forma d'onda (parametri ereditati dal modello analitico originale) ---
cfg.N            = 2^14;         % Campioni per simbolo (fast-time)
cfg.c_chirp      = 300 * pi;     % Chirp rate normalizzato
cfg.offset       = 100 * pi;     % Disassamento atomico -> posizione del picco
cfg.sigma_doppler= 0.15;         % Larghezza del fondo classico
cfg.sigma_noise  = 0.02;         % Rumore di canale (fondo cosmico/solare/shot)
cfg.K_doppler    = 1.578e4;      % Mappa (v/c) -> banda base normalizzata

% --- Banda IF fisica: chiude il legame fra N campioni e durata simbolo ---
% T_sym = N / B_IF. Con N = 16384 e B_IF = 16.384 GHz si riottiene
% esattamente T_sym = 1 us, il valore di progetto della versione originale.
cfg.B_IF         = 16.384e9;     % Banda del ricevitore IF (Hz)
cfg.T_sym        = cfg.N / cfg.B_IF;
R_b_raw          = 1 / cfg.T_sym;    % Bit rate lordo (bit/s), BPSK 1 bit/simbolo

% --- Griglia e operatori precalcolati ---
cfg.omega        = linspace(-1, 1, cfg.N).';    % colonna N x 1
cfg.domega       = cfg.omega(2) - cfg.omega(1);
cfg.win_clutter  = max(3, round(0.01 * cfg.N)); % finestra movmean scalata su N
cfg.chiave       = exp(-1i * (cfg.c_chirp * (cfg.omega.^2)));  % chiave FrFT
cfg.blocco       = 64;           % bit elaborati per blocco (controllo memoria)

% --- Parametri della campagna Monte Carlo ---
cfg.Nbit_frame   = 256;          % bit per frame (per trial)
cfg.N_MC         = 12;           % realizzazioni Monte Carlo per punto
SNR_anchor_dB    = -41;          % covertness di progetto al crossing
snr_vec_dB       = -41:-2.5:-73.5;   % sweep di covertness (Campagna A)
N_epoche         = 21;           % epoche orbitali campionate (Campagna B)
L_frame_goodput  = 1024;         % lunghezza frame per il calcolo del FER

% --- Guadagno di processo e soglia teorica ---
% Dopo il de-chirp il termine utile e' un tono puro: la FFT concentra il
% segnale in 1 bin su N -> guadagno di processo = 10*log10(N).
G_proc_dB = 10*log10(cfg.N);
disp(['Campioni per simbolo N = ', num2str(cfg.N), ...
      '  ->  guadagno di processo FrFT = +', num2str(G_proc_dB,'%.1f'), ' dB']);
disp(['Banda IF = ', num2str(cfg.B_IF/1e9,'%.3f'), ' GHz  ->  T_sym = ', ...
      num2str(cfg.T_sym*1e6,'%.3f'), ' us  ->  bit rate lordo = ', ...
      num2str(R_b_raw/1e3,'%.1f'), ' kbit/s']);

% Fotoni per simbolo al crossing (verifica di validita' del modello analogico)
n_ph_cross = (P_rx_t(idx_cross) * cfg.T_sym / E_fotone) * eta_det;
n_dark     = DCR * cfg.T_sym;
disp(['Fotoni di segnale per simbolo al crossing: ', num2str(n_ph_cross,'%.3e')]);
disp(['Dark counts per simbolo: ', num2str(n_dark,'%.3e')]);
if n_ph_cross < 10
    warning(['Regime di conteggio fotonico (n < 10 fotoni/simbolo): il modello ' ...
             'gaussiano analogico non e'' piu'' valido, servirebbe una statistica di Poisson.']);
else
    disp('Regime classico (n >> 1): il modello analogico gaussiano e'' applicabile.');
end

%% ========================================================================
% FASE 5: CAMPAGNA A - BER vs LIVELLO DI COVERTNESS (AL CROSSING)
% =========================================================================
disp('--- FASE 5: Campagna A - BER vs covertness (Monte Carlo) ---');

Shift_cross = (v_rel_t(idx_cross) / c_light) * cfg.K_doppler;
op_cross    = operatoreShift(cfg, Shift_cross);
[kpk_cross, tmpl_cross] = trovaBinPicco(cfg, op_cross, Shift_cross);

k_teorico = round(cfg.N/2) + 1 + round(cfg.offset/pi * cfg.N/(cfg.N-1));
if abs(kpk_cross - k_teorico) > 2
    warning(['Bin di picco misurato (', num2str(kpk_cross), ') distante da quello ' ...
             'teorico (', num2str(k_teorico), '): verificare offset/chirp.']);
end
disp(['Bin di decisione (noto a Bob dalle effemeridi): k = ', num2str(kpk_cross)]);

N_snr   = numel(snr_vec_dB);
err_A   = zeros(1, N_snr);
bit_A   = zeros(1, N_snr);
stat_nom = []; bits_nom = [];

for is = 1:N_snr
    for mc = 1:cfg.N_MC
        % Substream deterministico e univoco: (campagna, punto, trial)
        rs = streamMC(SEED_MASTER, 1, is, mc);
        [ne, nb, st, bt] = simulaFrame(cfg, op_cross, Shift_cross, ...
                                       snr_vec_dB(is), kpk_cross, rs);
        err_A(is) = err_A(is) + ne;
        bit_A(is) = bit_A(is) + nb;
        if is == 1 && mc == 1
            stat_nom = st; bits_nom = bt;   % per l'istogramma dello spazio di decisione
        end
    end
    fprintf('   SNR covert = %6.2f dB  ->  BER = %8.2e  (%d errori su %d bit)\n', ...
            snr_vec_dB(is), err_A(is)/bit_A(is), err_A(is), bit_A(is));
end

BER_A = err_A ./ bit_A;
[BER_A_lo, BER_A_hi] = wilsonCI(err_A, bit_A, 1.96);

% Curva analitica di riferimento (validata numericamente):
%   statistica di decisione antipodale  mu = (a/2)*N,  sigma = s*sqrt(N/2)
%   con a = ampiezza del coseno covert  ->  BER = Q( sqrt(P_cov*N)/s )
P_classica  = mean(exp(-((cfg.omega - Shift_cross).^2)/(2*cfg.sigma_doppler^2)).^2);
snr_fine_dB = linspace(min(snr_vec_dB)-2, max(snr_vec_dB)+2, 400);
P_cov_fine  = P_classica * 10.^(snr_fine_dB/10);
arg_fine    = sqrt(P_cov_fine * cfg.N) / cfg.sigma_noise;
BER_teorica = 0.5 * erfc(arg_fine / sqrt(2));

%% ========================================================================
% FASE 6: CAMPAGNA B - BER E THROUGHPUT LUNGO IL PASSAGGIO ORBITALE
% =========================================================================
disp('--- FASE 6: Campagna B - BER e throughput vs epoca orbitale ---');

% Accoppiamento fra link budget fisico e covertness:
% la covertness nominale e' ancorata al crossing; nelle altre epoche la
% potenza covert ricevuta segue il link budget reale (FSPL + fading PAT),
% quindi  SNR_eff(t) = SNR_anchor + 10*log10( P_rx(t) / P_rx(t_cross) ).
idx_epoche = round(linspace(1, N_t, N_epoche));
t_epoche   = t_sec(idx_epoche);

BER_B   = zeros(1, N_epoche);
SNR_B   = zeros(1, N_epoche);
err_B   = zeros(1, N_epoche);
bit_B   = zeros(1, N_epoche);
shift_B = zeros(1, N_epoche);

P_ref_ideal = P_rx_ideal_t(idx_cross);

for ie = 1:N_epoche
    ii = idx_epoche(ie);
    shift_B(ie) = (v_rel_t(ii) / c_light) * cfg.K_doppler;

    if abs(shift_B(ie)) > 0.6
        warning(['Epoca ', num2str(ie), ': shift Doppler normalizzato = ', ...
                 num2str(shift_B(ie),'%.2f'), ' -> fuori dal dominio utile del modello.']);
    end

    op_e = operatoreShift(cfg, shift_B(ie));
    kpk_e = trovaBinPicco(cfg, op_e, shift_B(ie));

    snr_acc = 0;
    for mc = 1:cfg.N_MC
        rs = streamMC(SEED_MASTER, 2, ie, mc);

        % Fading PAT ri-estratto ad ogni realizzazione (Rayleigh a 2 assi).
        % Il frame dura ~260 us: il jitter e' congelato entro il frame.
        th_tx = sigma_jitter_pat * randn(rs, 1, 2);
        th_rx = sigma_jitter_pat * randn(rs, 1, 2);
        Lp = exp(-8*(norm(th_tx)/theta_div_beam)^2) * ...
             exp(-8*(norm(th_rx)/theta_div_beam)^2);

        snr_eff_dB = SNR_anchor_dB + 10*log10( (P_rx_ideal_t(ii) * Lp) / P_ref_ideal );
        snr_acc = snr_acc + snr_eff_dB;

        [ne, nb] = simulaFrame(cfg, op_e, shift_B(ie), snr_eff_dB, kpk_e, rs);
        err_B(ie) = err_B(ie) + ne;
        bit_B(ie) = bit_B(ie) + nb;
    end
    SNR_B(ie) = snr_acc / cfg.N_MC;
    BER_B(ie) = err_B(ie) / bit_B(ie);

    fprintf('   t = %6.1f s | d = %7.1f km | SNR_eff = %6.2f dB | BER = %8.2e\n', ...
            t_epoche(ie), d_t(ii)/1000, SNR_B(ie), BER_B(ie));
end

[BER_B_lo, BER_B_hi] = wilsonCI(err_B, bit_B, 1.96);

%% ========================================================================
% FASE 7: METRICHE DI THROUGHPUT
% =========================================================================
disp('--- FASE 7: Metriche di throughput ---');

disponibile = double(link_available(idx_epoche));   % maschera LOS + Sole

% 1. Throughput "ingenuo": bit corretti al secondo
R_naive = R_b_raw * (1 - BER_B) .* disponibile;

% 2. Throughput informativo: capacita' del canale binario simmetrico
%    C = 1 - H2(BER). E' il limite superiore raggiungibile con codifica
%    ideale: la metrica corretta per rispondere a "quanta INFORMAZIONE passa".
C_bsc   = 1 - entropiaBinaria(BER_B);
R_shan  = R_b_raw * C_bsc .* disponibile;

% 3. Goodput non codificato con frame da L bit (nessuna correzione d'errore)
FER     = 1 - (1 - BER_B).^L_frame_goodput;
R_good  = R_b_raw * (1 - FER) .* disponibile;

% 4. Volume dati per passaggio (integrazione sulla finestra di 15 min)
Vol_naive = trapz(t_epoche, R_naive);
Vol_shan  = trapz(t_epoche, R_shan);
Vol_good  = trapz(t_epoche, R_good);

% 5. Finestra utile (BER sotto la soglia operativa)
BER_soglia = 1e-3;
utile      = (BER_B <= BER_soglia) & (disponibile > 0);
if any(utile)
    T_utile = max(t_epoche(utile)) - min(t_epoche(utile));
else
    T_utile = 0;
end

fprintf('\n===================== SINTESI THROUGHPUT =====================\n');
fprintf('Bit rate lordo (1 bit/simbolo, T_sym = %.3f us) : %8.1f kbit/s\n', cfg.T_sym*1e6, R_b_raw/1e3);
fprintf('Guadagno di processo FrFT (N = %d)              : %8.1f dB\n', cfg.N, G_proc_dB);
fprintf('Covertness di ancoraggio al crossing            : %8.1f dB\n', SNR_anchor_dB);
[~, ie_cross] = min(abs(idx_epoche - idx_cross));
fprintf('BER al crossing                                 : %8.2e\n', BER_B(ie_cross));
fprintf('Disponibilità del link (LOS + Sole)             : %7.1f %% della finestra\n', ...
        100*mean(link_available));
fprintf('Finestra utile (BER <= %.0e)                    : %8.1f s su %.0f s\n', BER_soglia, T_utile, T_window);
fprintf('--------------------------------------------------------------\n');
fprintf('Throughput di picco (1-BER)                     : %8.1f kbit/s\n', max(R_naive)/1e3);
fprintf('Throughput informativo di picco (BSC)           : %8.1f kbit/s\n', max(R_shan)/1e3);
fprintf('Goodput di picco, frame %d bit non codificati : %8.1f kbit/s\n', L_frame_goodput, max(R_good)/1e3);
fprintf('--------------------------------------------------------------\n');
fprintf('Volume dati per passaggio (1-BER)               : %8.2f Mbit\n', Vol_naive/1e6);
fprintf('Volume informativo per passaggio (BSC)          : %8.2f Mbit\n', Vol_shan/1e6);
fprintf('Volume utile per passaggio (frame non codif.)   : %8.2f Mbit\n', Vol_good/1e6);
fprintf('==============================================================\n\n');

% 6. Frontiera covertness / rate: dalla formula analitica validata,
%    per BER = BER_obiettivo servono N >= (gamma*sigma)^2 / P_cov campioni
%    per simbolo, e a banda IF fissa  R_b = B_IF / N. Ne segue
%    R_b_max = B_IF * P_cov / (gamma^2 * sigma^2): il rate massimo cresce
%    LINEARMENTE con la potenza covert. Il de-chirp FrFT non crea capacita':
%    converte banda in invisibilita'.
BER_obiettivo = 1e-6;
gamma_req     = sqrt(2) * erfcinv(2 * BER_obiettivo);
P_cov_sweep   = P_classica * 10.^(snr_fine_dB/10);
R_max_sweep   = cfg.B_IF * P_cov_sweep / (gamma_req^2 * cfg.sigma_noise^2);
N_req_anchor  = (gamma_req * cfg.sigma_noise)^2 / (P_classica * 10^(SNR_anchor_dB/10));
fprintf('Frontiera covertness/rate (BER obiettivo = %.0e):\n', BER_obiettivo);
fprintf('   a %.1f dB di covertness servono N >= %.0f campioni/simbolo\n', SNR_anchor_dB, ceil(N_req_anchor));
fprintf('   -> rate massimo sostenibile = %.1f kbit/s con B_IF = %.2f GHz\n\n', ...
        cfg.B_IF/max(N_req_anchor,1)/1e3, cfg.B_IF/1e9);

%% ========================================================================
% FASE 8: INTERCETTORE EVE - SINGOLO FRAME vs MEDIA D'INSIEME MULTI-BIT
% =========================================================================
disp('--- FASE 8: Analisi intercettatore Eve ---');

rs_eve = streamMC(SEED_MASTER, 3, 1, 1);
[~, ~, ~, ~, Rx_eve] = simulaFrame(cfg, op_cross, Shift_cross, ...
                                   SNR_anchor_dB, kpk_cross, rs_eve, true);

Y_Eve_single       = fftshift(fft(Rx_eve(:,1)));
Spec_Eve_single_dB = 20*log10(abs(Y_Eve_single) + eps);

% Media d'insieme su tutti i bit del frame (0/1 equiprobabili):
% il termine covert antipodale si elide, resta solo il fondo classico.
Y_Eve_ens          = fftshift(fft(mean(Rx_eve, 2)));
Spec_Eve_ens_dB    = 20*log10(abs(Y_Eve_ens) + eps);

% Costo di rivelazione per Eve (rivelatore d'energia, ordine di grandezza):
% con SNR covert = eta, servono ~1/eta^2 campioni indipendenti per
% distinguere le due ipotesi. Ipotesi forte: Eve conosce esattamente la
% statistica del fondo classico.
eta_cov  = 10^(SNR_anchor_dB/10);
N_eve    = 1 / eta_cov^2;
T_eve    = N_eve / cfg.B_IF;
fprintf('Rivelazione d''energia di Eve a %.0f dB: ~%.2e campioni (~%.2e s a %.2f GHz)\n', ...
        SNR_anchor_dB, N_eve, T_eve, cfg.B_IF/1e9);

%% ========================================================================
% FASE 9: FIGURE
% =========================================================================
disp('--- FASE 9: Generazione figure ---');

u_axis = linspace(-1, 1, cfg.N);

% ---------------- FIGURA 1: geometria, link budget, Eve ------------------
figure('Name', 'Quantum-Covert ISL - Geometria, Link Budget, Eve', 'Color', 'w', ...
       'Position', [80, 50, 1050, 800]);

subplot(3,1,1);
yyaxis left;
h_dist = plot(t_sec, d_t/1000, 'b-', 'LineWidth', 1.3);
ylabel('Distanza ISL d(t) [km]', 'Color', 'b', 'FontWeight', 'bold');
ax = gca; ax.YColor = 'b';
ylim([0, max(d_t/1000)*1.1]);
yyaxis right;
h_vrel = plot(t_sec, v_rel_t/1000, 'r-', 'LineWidth', 1.3);
hold on; yline(0, 'k:', 'LineWidth', 0.8);
ylabel('v_{rel}(t) [km/s]', 'Color', 'r', 'FontWeight', 'bold');
ax.YColor = 'r';
xlabel('Tempo di simulazione [s]');
title(['(a) Geometria Orbitale ISL: crossing a t = ', num2str(t_sec(idx_cross),'%.1f'), ...
       ' s (d_{min} = ', num2str(d_min/1000,'%.1f'), ' km)']);
legend([h_dist, h_vrel], {'Distanza d(t) [sx]', 'Velocità relativa LOS [dx]'}, ...
       'Location', 'northeast', 'FontSize', 8);
grid on;

subplot(3,1,2);
h_p_ideal = plot(t_sec, 10*log10(P_rx_ideal_t*1e3), 'Color', [0.55 0.55 0.55], 'LineStyle', '--', 'LineWidth', 1.2);
hold on;
h_p_jit   = plot(t_sec, 10*log10(P_rx_t*1e3), 'Color', [0 0.45 0.75], 'LineWidth', 1.1);
xline(t_sec(idx_cross), 'k:', 'Crossing', 'LineWidth', 1.0);
xlabel('Tempo di simulazione [s]');
ylabel('P_{rx}(t) [dBm]');
legend([h_p_ideal, h_p_jit], {'P_{rx} ideale (FSPL + efficienze)', ...
       'P_{rx} con jitter di Rayleigh (\sigma = 1.5 \murad)'}, 'Location', 'best', 'FontSize', 8);
title('(b) Link Budget Ottico FSO Dinamico con Perdite PAT');
grid on;

subplot(3,1,3);
h_eve_single = plot(u_axis, Spec_Eve_single_dB, 'Color', [0.25 0.25 0.25], 'LineWidth', 0.9);
hold on;
h_eve_ens    = plot(u_axis, Spec_Eve_ens_dB, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);
xlabel('Frequenza spaziale normalizzata u');
ylabel('Densità spettrale [dB]');
legend([h_eve_single, h_eve_ens], ...
       {'Eve: singolo simbolo (FFT standard)', ...
        ['Eve: media su ', num2str(cfg.Nbit_frame), ' simboli (elisione del termine covert)']}, ...
       'Location', 'northeast', 'FontSize', 8);
title(['(c) Intercettore Eve: invisibilità energetica (SNR = ', num2str(SNR_anchor_dB), ' dB)']);
xlim([-0.25 0.25]); grid on;

% ---------------- FIGURA 2: ambiente orbitale ----------------------------
figure('Name', 'Analisi Geometria Orbitale e Ambiente Solare', 'Color', 'w', ...
       'Position', [100, 100, 900, 600]);

subplot(2,1,1);
plot(t_sec, theta_sun_rx_deg, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.4);
hold on;
yline(theta_sun_crit, 'r--', ['Soglia abbagliamento (', num2str(theta_sun_crit), '°)'], 'LineWidth', 1.2);
xline(t_sec(idx_cross), 'k:', 'Crossing', 'LineWidth', 1.0);
grid on; xlabel('Tempo di simulazione [s]'); ylabel('Angolo Sole-LOS [deg]');
title('Angolo di Esclusione Solare al Ricevitore Bob (\theta_{Sun})');

subplot(2,1,2);
plot(t_sec, los_clearance_margin / 1000, 'Color', [0 0.5 0], 'LineWidth', 1.4);
hold on;
yline(0, 'r--', 'Limite atmosfera densa (80 km)', 'LineWidth', 1.2);
xline(t_sec(idx_cross), 'k:', 'Crossing', 'LineWidth', 1.0);
grid on; xlabel('Tempo di simulazione [s]'); ylabel('Margine di clearance [km]');
title('Margine di Visibilità Ottica Intersatellitare (LOS Clearance)');

% ---------------- FIGURA 3: statistiche BER e throughput -----------------
figure('Name', 'Campagna Monte Carlo: BER e Throughput', 'Color', 'w', ...
       'Position', [120, 40, 1080, 900]);

% (a) BER vs covertness
subplot(2,2,1);
ber_plot = BER_A; ber_plot(ber_plot == 0) = NaN;   % zeri non rappresentabili in log
h_th = semilogy(snr_fine_dB, BER_teorica, 'k-', 'LineWidth', 1.3); hold on;
h_mc = semilogy(snr_vec_dB, ber_plot, 'o', 'Color', [0 0.45 0.74], ...
                'MarkerFaceColor', [0 0.45 0.74], 'MarkerSize', 5);
h_ul = semilogy(snr_vec_dB, BER_A_hi, 'v--', 'Color', [0.5 0.5 0.5], 'MarkerSize', 4);
yline(BER_soglia, 'r:', 'BER = 10^{-3}', 'LineWidth', 1.1);
xlabel('SNR covert rispetto al fondo classico [dB]');
ylabel('BER');
title(['(a) BER vs covertness - ', num2str(cfg.N_MC*cfg.Nbit_frame), ' bit/punto']);
legend([h_th, h_mc, h_ul], {'Analitica: Q(\surd(P_{cov}N)/\sigma)', ...
       'Monte Carlo', 'Limite sup. IC 95% (Wilson)'}, 'Location', 'southwest', 'FontSize', 8);
ylim([1e-5 1]); grid on; set(gca,'XDir','reverse');

% (b) Spazio di decisione
subplot(2,2,2);
if ~isempty(stat_nom)
    histogram(stat_nom(bits_nom == 0), 30, 'FaceColor', [0 0.45 0.74], 'FaceAlpha', 0.6); hold on;
    histogram(stat_nom(bits_nom == 1), 30, 'FaceColor', [0.85 0.33 0.1], 'FaceAlpha', 0.6);
    xline(0, 'k--', 'LineWidth', 1.2);
end
xlabel('Statistica di decisione  Re\{Y(k_{picco})\}');
ylabel('Occorrenze');
title(['(b) Spazio di decisione a ', num2str(snr_vec_dB(1),'%.1f'), ' dB (1 frame)']);
legend({'Bit 0 (\phi = 0)', 'Bit 1 (\phi = \pi)', 'Soglia'}, 'Location', 'best', 'FontSize', 8);
grid on;

% (c) BER lungo il passaggio
subplot(2,2,3);
ber_b_plot = BER_B; ber_b_plot(ber_b_plot == 0) = NaN;
yyaxis left;
semilogy(t_epoche, ber_b_plot, 'o-', 'LineWidth', 1.3, 'MarkerSize', 4); hold on;
semilogy(t_epoche, BER_B_hi, 'v:', 'Color', [0.5 0.5 0.5], 'MarkerSize', 4);
ylabel('BER'); ylim([1e-5 1]);
yyaxis right;
plot(t_epoche, SNR_B, 's--', 'LineWidth', 1.1, 'MarkerSize', 4);
ylabel('SNR covert efficace [dB]');
xline(t_sec(idx_cross), 'k:', 'Crossing', 'LineWidth', 1.0);
xlabel('Tempo di simulazione [s]');
title('(c) BER ed SNR efficace lungo il passaggio');
grid on;

% (d) Throughput e volume cumulato
subplot(2,2,4);
yyaxis left;
h_r1 = plot(t_epoche, R_naive/1e3, 'o-', 'LineWidth', 1.3, 'MarkerSize', 4); hold on;
h_r2 = plot(t_epoche, R_shan/1e3,  's--', 'LineWidth', 1.2, 'MarkerSize', 4);
h_r3 = plot(t_epoche, R_good/1e3,  'd:',  'LineWidth', 1.2, 'MarkerSize', 4);
ylabel('Throughput [kbit/s]');
yyaxis right;
Vol_cum = cumtrapz(t_epoche, R_shan);
h_v = plot(t_epoche, Vol_cum/1e6, 'LineWidth', 1.4);
ylabel('Volume informativo cumulato [Mbit]');
xlabel('Tempo di simulazione [s]');
title('(d) Throughput istantaneo e volume per passaggio');
legend([h_r1, h_r2, h_r3, h_v], {'R_b(1-BER)', 'R_b(1-H_2(BER)) [BSC]', ...
       ['Goodput frame ', num2str(L_frame_goodput), ' bit'], 'Volume cumulato'}, ...
       'Location', 'best', 'FontSize', 7);
grid on;

% ---------------- FIGURA 4: frontiera covertness / rate ------------------
figure('Name', 'Frontiera Covertness-Rate', 'Color', 'w', 'Position', [140, 80, 700, 420]);
loglog(10.^(snr_fine_dB/10), R_max_sweep/1e3, 'LineWidth', 1.6); hold on;
xline(10^(SNR_anchor_dB/10), 'k--', 'Punto di progetto', 'LineWidth', 1.1);
yline(R_b_raw/1e3, 'r:', 'Rate lordo attuale', 'LineWidth', 1.1);
grid on;
xlabel('SNR covert (lineare, rispetto al fondo classico)');
ylabel('Bit rate massimo sostenibile [kbit/s]');
title(['Frontiera covertness-rate a B_{IF} = ', num2str(cfg.B_IF/1e9,'%.2f'), ...
       ' GHz, BER obiettivo = ', num2str(BER_obiettivo,'%.0e')]);

disp('--- CAMPAGNA MONTE CARLO COMPLETATA CON SUCCESSO! ---');

%% ========================================================================
% FUNZIONI LOCALI
% =========================================================================

function rs = streamMC(seed_master, id_campagna, id_punto, id_trial)
% Substream deterministico e univoco a partire da un unico seed master.
% Sostituendo lo stesso SEED_MASTER si riottiene identicamente la campagna.
    sub = 1 + id_trial + 1000*id_punto + 1000000*id_campagna;
    rs  = RandStream('mrg32k3a', 'Seed', seed_master);
    rs.Substream = sub;
end

function op = operatoreShift(cfg, shift)
% Precalcola l'operatore di ri-centraggio Doppler (interpolazione lineare
% a passo costante). Equivale a interp1(...,'linear',0) ma vettorizzato su
% matrici di simboli, con un costo trascurabile.
    pos = (cfg.omega + shift - cfg.omega(1)) / cfg.domega + 1;
    i0  = floor(pos);
    fr  = pos - i0;
    ok  = (i0 >= 1) & (i0 + 1 <= cfg.N);
    op.i0 = i0(ok);
    op.fr = fr(ok);
    op.ok = ok;
    op.N  = cfg.N;
end

function Sc = applicaShift(op, S)
    Sc = zeros(op.N, size(S,2));
    Sc(op.ok, :) = (1 - op.fr) .* S(op.i0, :) + op.fr .* S(op.i0 + 1, :);
end

function Y = riceviBob(cfg, op, S)
% Catena di ricezione: ri-centraggio Doppler da effemeridi -> clutter
% rejection (rimozione dell'inviluppo classico lento) -> de-chirp FrFT ->
% compressione spettrale.
    Sc = applicaShift(op, S);
    Sc = Sc - movmean(Sc, cfg.win_clutter, 1);
    Y  = fftshift(fft(Sc .* cfg.chiave, [], 1), 1);
end

function [kpk, tmpl] = trovaBinPicco(cfg, op, shift)
% Bin di compressione atteso, ricavato dalla forma d'onda nota (senza
% rumore ne' fondo). E' l'informazione che Bob possiede per costruzione:
% conosce chirp, offset e shift Doppler dalle effemeridi.
    Phi  = cfg.c_chirp * ((cfg.omega - shift).^2) + cfg.offset * (cfg.omega - shift);
    y    = riceviBob(cfg, op, cos(Phi));
    tmpl = real(y);
    [~, kpk] = max(abs(tmpl));
end

function [nerr, nbit, stat, bits, Rx_out] = simulaFrame(cfg, op, shift, snr_db, kpk, rs, salva_rx)
% Trasmette un frame di cfg.Nbit_frame bit BPSK sul chirp, lo propaga nel
% canale (fondo classico + rumore) e lo demodula. Elaborazione a blocchi
% per limitare l'occupazione di memoria.
    if nargin < 7, salva_rx = false; end
    M    = cfg.Nbit_frame;
    bits = randi(rs, [0 1], 1, M);

    Phi   = cfg.c_chirp * ((cfg.omega - shift).^2) + cfg.offset * (cfg.omega - shift);
    Scl   = exp(-((cfg.omega - shift).^2) / (2 * cfg.sigma_doppler^2));

    % Normalizzazione di potenza per covertness energetica: la potenza del
    % simbolo covert e' fissata rispetto a quella del fondo classico.
    P_cl  = mean(Scl.^2);
    P_cov = P_cl * 10^(snr_db/10);
    a_cov = sqrt(2 * P_cov);            % ampiezza del coseno (P = a^2/2)

    stat = zeros(1, M);
    if salva_rx
        Rx_out = zeros(cfg.N, M);
    else
        Rx_out = [];
    end

    for i1 = 1:cfg.blocco:M
        i2  = min(i1 + cfg.blocco - 1, M);
        b   = bits(i1:i2);
        Stx = a_cov * cos(Phi + pi * b);            % N x blocco
        Rx  = Scl + Stx + cfg.sigma_noise * randn(rs, cfg.N, numel(b));
        if salva_rx, Rx_out(:, i1:i2) = Rx; end
        Y   = riceviBob(cfg, op, Rx);
        stat(i1:i2) = real(Y(kpk, :));
    end

    dec  = double(stat < 0);            % cos(0) = +1 -> bit 0 ; cos(pi) = -1 -> bit 1
    nerr = sum(dec ~= bits);
    nbit = M;
end

function [lo, hi] = wilsonCI(k, n, z)
% Intervallo di confidenza di Wilson per una proporzione binomiale.
% Preferibile all'approssimazione normale quando k e' piccolo o nullo.
    p = k ./ n;
    d = 1 + z^2 ./ n;
    c = p + z^2 ./ (2*n);
    s = z * sqrt( p.*(1-p)./n + z^2 ./ (4*n.^2) );
    lo = max((c - s) ./ d, 0);
    hi = min((c + s) ./ d, 1);
end

function H = entropiaBinaria(p)
% Entropia binaria H2(p), con H2(0) = H2(1) = 0.
    p = min(max(p, 0), 1);
    H = zeros(size(p));
    m = (p > 0) & (p < 1);
    H(m) = -p(m).*log2(p(m)) - (1-p(m)).*log2(1-p(m));
end

% =========================================================================
% NOTE TECNICHE
% =========================================================================
% NOTA 1 - RIVELATORE.
%   La versione a bit singolo decideva con  [~, idx] = max(abs(Picco))  e poi
%   guardava il segno in quell'indice. Con un solo bit ad alto SNR funziona,
%   ma in una campagna BER e' una trappola: sotto la soglia il massimo cade su
%   un picco di rumore in posizione casuale e la BER satura a 0.5 invece di
%   degradare con continuita'. Verifica numerica del solo blocco di segnale
%   (N = 16384, 400 bit): a SNR = -60 dB la BER con bin noto e' 5e-3, con
%   argmax e' 0.50. Qui il bin di picco e' calcolato una volta per epoca dalla
%   forma d'onda nota (Bob conosce chirp, offset e shift dalle effemeridi):
%   e' il rivelatore coerente corretto, ed e' anche piu' favorevole.
%
% NOTA 2 - GUADAGNO D'ANTENNA OTTICA.
%   G = (pi*D/lambda)^2 vale per un fascio al limite di diffrazione
%   (lambda/D = 1.34 urad per D = 0.2 m), mentre il modello PAT dichiara una
%   divergenza di 8 urad. Le due ipotesi sono incompatibili: usando la
%   divergenza reale il guadagno Tx cala di circa 15.5 dB. Il flag
%   USA_GUADAGNO_COERENTE in Fase 3 permette di verificarne l'impatto.
%
% NOTA 3 - REGIME FOTONICO.
%   Con i parametri di progetto il link riceve dell'ordine di 1e10 fotoni per
%   simbolo: e' un regime pienamente classico. La covertness qui e' di tipo
%   LPI/spread-spectrum (forma d'onda nascosta sotto il fondo), non una
%   proprieta' quantistica: nessuna garanzia di sicurezza in senso QKD
%   discende da questo modello.
%
% NOTA 4 - K_doppler e BANDA IF.
%   K_doppler = 1.578e4 e' una costante di scala grafica ereditata dal modello
%   originale, non una derivazione fisica. Il legame fisico e' invece imposto
%   qui da B_IF: T_sym = N / B_IF. Cambiando B_IF cambiano coerentemente
%   durata di simbolo, bit rate e conteggio fotonico.
%
% NOTA 5 - LIMITE FONDAMENTALE.
%   Dalla formula analitica (validata contro il Monte Carlo) il rate massimo a
%   BER fissata vale  R_b = B_IF * P_cov / (gamma^2 * sigma^2): e' lineare
%   nella potenza covert e indipendente da N. Il de-chirp FrFT non genera
%   capacita', scambia banda con invisibilita'. Ogni 3 dB di covertness in piu'
%   costano esattamente un dimezzamento del bit rate.
% =========================================================================