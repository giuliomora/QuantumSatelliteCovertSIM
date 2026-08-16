% =========================================================================
% SIMULAZIONE DINAMICA LINK INTERSATELLITARE QUANTUM-COVERT (LEO-LEO ISL)
% Estensione orbitale con Satellite Communications Toolbox
% Modello: Ione 27Al+ | Doppler Dinamico & FrFT De-chirping
% =========================================================================
 
clc; clear; close all;
rng(42);   % Riproducibilita' di rumore AWGN e pointing jitter
 
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

% Costante armonica zonale J2 terrestre
J2 = 1.08263e-3;

% Tasso secolare di precessione nodale dRAAN/dt dovuto a J2 (rad/s)
mu_n      = sqrt(mu_earth / sma^3);     % Moto medio (rad/s)
p_semi    = sma * (1 - ecc^2);
raan_rate = -1.5 * mu_n * J2 * (Re / p_semi)^2 * cosd(incl_deg); % rad/s
disp(['Precessione nodale J2 (dRAAN/dt): ', num2str(rad2deg(raan_rate)*86400,'%.4f'), ' deg/giorno']);

% Configurazione RAAN per incrocio di 60°
theta_cross_deg = 60;
cos_dRAAN = (cosd(theta_cross_deg) - cosd(incl_deg)^2) / sind(incl_deg)^2;
dRAAN_deg = acosd(cos_dRAAN);
raan_Alice = 0;
raan_Bob   = dRAAN_deg;

T_orbit   = 2*pi / mu_n;
incl_rad  = deg2rad(incl_deg);
raanA_rad = deg2rad(raan_Alice);
raanB_rad = deg2rad(raan_Bob);

% Propagatore analitico per la stima iniziale con correzione di precessione J2
posCircularJ2 = @(u0, raan0, t) sma * [ ...
    cos(raan0 + raan_rate*t).*cos(u0 + mu_n*t) - sin(raan0 + raan_rate*t).*sin(u0 + mu_n*t).*cos(incl_rad); ...
    sin(raan0 + raan_rate*t).*cos(u0 + mu_n*t) + cos(raan0 + raan_rate*t).*sin(u0 + mu_n*t).*cos(incl_rad); ...
    sin(u0 + mu_n*t).*sin(incl_rad) ];

target_distance = 300e3;                % 300 km al fly-by

% Ricerca deterministica della fase relativa
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

% Raffinamento locale fine
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

% Configurazione con propagatore SGP4 (integra J2, J3, J4 + drag)
satAlice = satellite(sc, sma, ecc, incl_deg, raan_Alice, argPeri_deg, ...
    trueAnom_Alice, "OrbitPropagator", "sgp4", "Name", "Alice");
satBob   = satellite(sc, sma, ecc, incl_deg, raan_Bob, argPeri_deg, ...
    trueAnom_Bob, "OrbitPropagator", "sgp4", "Name", "Bob");

[posAlice, velAlice, tSamples] = states(satAlice, "CoordinateFrame", "inertial");
[posBob,   velBob,   ~]        = states(satBob,   "CoordinateFrame", "inertial");

t_sec = seconds(tSamples - tSamples(1));
N_t   = numel(t_sec);

% 1. Cinematica ISL e Doppler
deltaPos = posBob - posAlice;              % Vettore LOS da Alice a Bob (3xN)
d_t      = vecnorm(deltaPos, 2, 1);        % Distanza ISL (m)
tau_t    = d_t / c_light;                  % Ritardo propagazione (s)

uLOS_AliceToBob = deltaPos ./ d_t;         % Versore di puntamento da Alice a Bob
uLOS_BobToAlice = -uLOS_AliceToBob;        % Versore di puntamento da Bob ad Alice (verso Tx)

deltaVel = velBob - velAlice;
v_rel_t  = sum(deltaVel .* uLOS_AliceToBob, 1); % Range-rate (m/s)

f_D_t    = -(v_rel_t / c_light) * f0_optical;   % Doppler ottico (Hz)
dfD_dt_t = gradient(f_D_t, sampleTime);         % Doppler rate (Hz/s)

% 2. Controllo Occlusione Terrestre (Earth Line-of-Sight Clearance)
% Calcola la quota minima sfiorata dal raggio ottico sopra la Terra
h_atmosphere_limit = 80e3;                      % Atmosfera densa (m)
r_min_los = zeros(1, N_t);
for i = 1:N_t
    r_A = posAlice(:,i);
    r_B = posBob(:,i);
    % Distanza minima della retta passante per A e B dall'origine ECI
    r_min_los(i) = norm(cross(r_A, r_B)) / norm(r_B - r_A);
end
los_clearance_margin = r_min_los - (Re + h_atmosphere_limit);
los_blocked = any(los_clearance_margin < 0);

% 3. Vettore Posizione del Sole e Angolo di Congiunzione (Sun Exclusion Angle)
% Approssimazione analitica accurata della posizione del Sole in ECI (J2000)
% (1 AU = 1.496e11 m, obliquità dell'eclittica eps_ecl = 23.44°)
AU = 1.495978707e11;
eps_ecl = deg2rad(23.439291);
n_days = days(tSamples - datetime(2000,1,1,12,0,0, 'TimeZone', 'UTC')); % Giorni da J2000.0
lambda_sun = deg2rad(mod(280.460 + 0.9856474 * n_days, 360)); % Longitudine eclittica apparente

% Coordinate geocentriche inerziali (ECI) del Sole (3xN)
posSun = AU * [cos(lambda_sun); ...
               sin(lambda_sun)*cos(eps_ecl); ...
               sin(lambda_sun)*sin(eps_ecl)];

% Calcolo del Sun-Exclusion Angle per il ricevitore Bob:
% Angolo tra la LOS ottica (direzione da cui arriva il laser, verso Alice)
% e la direzione verso il Sole (vettore Bob -> Sole)
sun_dir_Bob = posSun - posBob;
uSun_Bob    = sun_dir_Bob ./ vecnorm(sun_dir_Bob, 2, 1);

% theta_sun_exclusion: angolo tra ricezione ottica e centroide solare
cos_theta_sun = sum(uLOS_BobToAlice .* uSun_Bob, 1);
theta_sun_rx_deg = rad2deg(acos(min(max(cos_theta_sun, -1), 1)));

% Soglia minima di esclusione tipica dei telescopi ottici spaziali
theta_sun_crit = 30; % gradi
sun_blindness_flag = theta_sun_rx_deg < theta_sun_crit;

% Statistiche a console
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
 
%% ========================================================================
% FASE 3: MODELLO DI CANALE OTTICO SATELLITARE STANDARDIZZATO (FSO / PAT)
% =========================================================================
disp('--- FASE 3: Standardizzazione Canale Ottico (FSO Link Budget & PAT) ---');

% -------------------------------------------------------------------------
% 1. PARAMETRI OTTICI E OPTO-ELETTRONICI DI SISTEMA
% -------------------------------------------------------------------------
D_tx        = 0.20;         % Diametro apertura telescopio Tx Alice (m)
D_rx        = 0.30;         % Diametro apertura telescopio Rx Bob (m)
P_tx_laser  = 100e-3;       % Potenza media trasmessa dal laser (W) -> 20 dBm

% Efficienze ottiche e quantiche
eta_tx      = 0.85;         % Efficienza ottica Alice (specchi, lenti, accoppiamento)
eta_rx      = 0.80;         % Efficienza ottica Bob (baffle, filtri interferenziali)
eta_det     = 0.65;         % Efficienza quantica del rivelatore UV (SNSPD / APD)
T_sym       = 1e-6;         % Durata temporale del simbolo quantistico (1 microsecondo)
DCR         = 50;           % Dark Count Rate del rivelatore (counts/s)

% Parametri di Puntamento (PAT - Pointing, Acquisition, and Tracking)
sigma_jitter_pat = 1.5e-6;  % Jitter 1-sigma per singolo asse (rad, 1.5 urad)
theta_div_beam   = 8.0e-6;  % Divergenza nominale di progetto del fascio (rad)

% -------------------------------------------------------------------------
% 2. GUADAGNI E ATTENUAZIONE IN SPAZIO LIBERO (FSPL)
% -------------------------------------------------------------------------
% Guadagni d'antenna ottica di apertura (standard CCSDS / ITU-R P.2148)
G_tx = (pi * D_tx / lambda_0)^2;
G_rx = (pi * D_rx / lambda_0)^2;

% FSPL dinamico calcolato lungo la traiettoria orbitale SGP4
FSPL_t = (4 * pi * d_t / lambda_0).^2;

% -------------------------------------------------------------------------
% 3. MODELLO DI JITTER DI PUNTAMENTO (DISTRIBUZIONE DI RAYLEIGH)
% -------------------------------------------------------------------------
% L'errore angolare su Azimuth ed Elevazione e' disaccoppiato e Gaussiano
theta_az_tx = sigma_jitter_pat * randn(1, N_t);
theta_el_tx = sigma_jitter_pat * randn(1, N_t);
theta_err_tx_radiale = sqrt(theta_az_tx.^2 + theta_el_tx.^2); % Rayleigh

theta_az_rx = sigma_jitter_pat * randn(1, N_t);
theta_el_rx = sigma_jitter_pat * randn(1, N_t);
theta_err_rx_radiale = sqrt(theta_az_rx.^2 + theta_el_rx.^2); % Rayleigh

% Attenuazione PAT (modello di disallineamento esponenziale per fasci Gaussiani)
L_point_tx = exp(-8 * (theta_err_tx_radiale ./ theta_div_beam).^2);
L_point_rx = exp(-8 * (theta_err_rx_radiale ./ theta_div_beam).^2);
L_point_t  = L_point_tx .* L_point_rx;

% -------------------------------------------------------------------------
% 4. LINK BUDGET DINAMICO E CONTEGGIO FOTONICO
% -------------------------------------------------------------------------
% Potenza ottica ideale e con fading stocastico
P_rx_ideal_t = P_tx_laser * eta_tx * eta_rx * (G_tx * G_rx) ./ FSPL_t;
P_rx_t       = P_rx_ideal_t .* L_point_t;

% Flusso medio di fotoni di segnale incidenti sul rivelatore per simbolo
n_photons_signal_t = (P_rx_t * T_sym / E_fotone) * eta_det;

% Fotoni spuri di fondo / Dark counts per simbolo
n_photons_noise_t  = (DCR * T_sym) * ones(1, N_t);

% -------------------------------------------------------------------------
% STATISTICHE DEL CANALE A CONSOLE
% -------------------------------------------------------------------------
P_rx_crossing_dBm = 10*log10(P_rx_t(idx_cross)*1e3);
L_point_avg_dB    = 10*log10(mean(L_point_t));

disp(['Potenza Rx al crossing (con jitter ed efficienze): ', num2str(P_rx_crossing_dBm,'%.2f'), ' dBm']);
disp(['Perdita media di puntamento (PAT Loss): ', num2str(L_point_avg_dB,'%.2f'), ' dB']);
disp(['Numero medio di fotoni per simbolo al crossing: ', num2str(n_photons_signal_t(idx_cross),'%.2e'), ' fotoni/simbolo']);
disp(['Dark Counts attesi per simbolo: ', num2str(n_photons_noise_t(idx_cross),'%.2e'), ' fotoni/simbolo']);

%% ========================================================================
% FASE 4: SEGNALE QUANTISTICO - UN SOLO BIT SELEZIONABILE (ALICE)
% =========================================================================
disp('--- FASE 4: Generazione del simbolo quantistico trasmesso ---');
 
bit_to_send = 0;    % <<< MODIFICABILE A MANO: 0 oppure 1
 
% Istante di trasmissione: punto di crossing ISL (minima distanza), dove
% il link e' geometricamente valido e il Doppler e' definito dalla
% dinamica orbitale reale appena estratta.
idx_tx    = idx_cross;
v_rel_tx  = v_rel_t(idx_tx);
f_D_tx    = f_D_t(idx_tx);
dfD_dt_tx = dfD_dt_t(idx_tx);
 
N_campioni = 100000;                 % Risoluzione fast-time ad alta densita'
omega = linspace(-1, 1, N_campioni); % Banda base normalizzata
 
% Parametri chirp quantistico - IDENTICI al modello analitico originale
c_chirp            = 300 * pi;   % Chirp rate (rad, normalizzato)
offset_quantistico = 100 * pi;   % Termine di disassamento atomico
I_quant             = 1e-5;
 
% Fattore di normalizzazione IF: mappa lo shift Doppler frazionario reale
% (v_rel/c) sulla banda base normalizzata omega in [-1,1]. Costante di
% scala ereditata dal modello analitico originale (rappresenta la banda
% IF nominale del ricevitore); non e' una derivazione fisica esatta.
K_doppler_grafico = 1.578e4;
Shift_Orbitale = (v_rel_tx / c_light) * K_doppler_grafico;
 
disp(['Bit trasmesso (scelto a mano): ', num2str(bit_to_send)]);
disp(['Shift Doppler normalizzato (dinamico, dal crossing reale): ', ...
      num2str(Shift_Orbitale,'%.4f')]);
 
%% ========================================================================
% FASE 5: CANALE ISL - COVERTNESS (SNR = -41 dB) CON LE FORMULE ORIGINALI
% =========================================================================
disp('--- FASE 5: Iniezione nel canale (formule identiche al modello originale) ---');
 
% Profilo Doppler-broadening termico/classico (fondo), centrato sul
% medesimo shift Doppler dinamico usato per il simbolo quantistico
sigma_doppler = 0.15;
S_classico = exp(-((omega - Shift_Orbitale).^2) / (2 * sigma_doppler^2));
 
% Fase quantistica totale traslata in frequenza (stesso c_chirp usato poi
% nella chiave di de-chirp di Bob: nessun disallineamento del matched
% filter). Si generano ENTRAMBI i simboli di riferimento (necessari alla
% dimostrazione statistica di covertness di Eve in Fase 7), ma solo quello
% corrispondente a bit_to_send viene realmente "trasmesso" e ricevuto da
% Bob in Fase 6.
Phi_shiftata = c_chirp * ((omega - Shift_Orbitale).^2) + ...
               offset_quantistico * (omega - Shift_Orbitale);
S_quant_tx0 = I_quant * cos(Phi_shiftata + 0);    % Simbolo di riferimento Bit 0
S_quant_tx1 = I_quant * cos(Phi_shiftata + pi);   % Simbolo di riferimento Bit 1
 
% Normalizzazione di potenza per covertness energetica (SNR = -41 dB),
% IDENTICA al modello originale (scorporata dal link budget dinamico:
% quest'ultimo descrive la potenza ottica assoluta del link, mentre la
% covertness e' definita relativamente al fondo classico in banda base).
SNR_target_dB    = -41;
Potenza_Classica = mean(S_classico.^2);
Potenza_Target   = Potenza_Classica * 10^(SNR_target_dB / 10);
Fattore_Scala    = sqrt(Potenza_Target / mean(S_quant_tx0.^2));
S_quant_tx0 = S_quant_tx0 * Fattore_Scala;
S_quant_tx1 = S_quant_tx1 * Fattore_Scala;
 
% Rumore di canale (fondo cosmico / solare / shot noise)
Rumore_Spaziale = 0.02 * randn(1, N_campioni);
 
S_rx_bit0 = S_classico + S_quant_tx0 + Rumore_Spaziale;  % riferimento (per Eve)
S_rx_bit1 = S_classico + S_quant_tx1 + Rumore_Spaziale;  % riferimento (per Eve)
 
% Segnale FISICAMENTE trasmesso e ricevuto da Bob in questa esecuzione
if bit_to_send == 0
    S_rx_bob = S_rx_bit0;
else
    S_rx_bob = S_rx_bit1;
end
 
%% ========================================================================
% FASE 6: RICEVITORE DI BOB - TRACKING, CLUTTER REJECTION, FrFT DE-CHIRP
% =========================================================================
disp('--- FASE 6: Demodulazione coerente di Bob (formule originali) ---');
 
% 1. Doppler tracking basato su effemeridi (ri-centraggio dello shift noto
%    dalla soluzione orbitale, non stimato ciecamente da Bob)
S_centrato_bob = interp1(omega, S_rx_bob, omega + Shift_Orbitale, 'linear', 0);
 
% 2. Clutter rejection (rimozione inviluppo lento classico)
S_clean_bob = S_centrato_bob - movmean(S_centrato_bob, 1000);
 
% 3. Chiave di de-chirp FrFT - stesso c_chirp usato in trasmissione
Chiave_FrFT = exp(-1i * (c_chirp * (omega.^2)));
 
% 4. Compressione spettrale (Fourier frazionaria) + decisione sul bit
Y_Bob_bob    = fftshift(fft(S_clean_bob .* Chiave_FrFT));
Picco_Bob    = real(Y_Bob_bob);
 
[peak_amp, peak_idx] = max(abs(Picco_Bob));
peak_signed = Picco_Bob(peak_idx);
decoded_bit = double(peak_signed < 0);   % cos(0)=+1 -> bit 0 ; cos(pi)=-1 -> bit 1
 
rumore_rms = std(Picco_Bob([1:10000, end-10000:end]));
SNR_post_Bob_dB = 20*log10(abs(peak_amp) / rumore_rms);
 
disp(['Picco di decisione: ', num2str(peak_signed,'%.3f'), ' -> bit decodificato = ', num2str(decoded_bit)]);
disp(['Processing gain: SNR post-FrFT = +', num2str(SNR_post_Bob_dB,'%.2f'), ' dB']);
if decoded_bit == bit_to_send
    disp('DECODIFICA CORRETTA: il bit ricevuto coincide con il bit trasmesso.');
else
    disp('ERRORE DI DECODIFICA: il bit ricevuto NON coincide con il bit trasmesso.');
end
 
%% ========================================================================
% FASE 7: INTERCETTORE EVE - FFT SINGOLO SHOT E MEDIA STATISTICA
% =========================================================================
disp('--- FASE 7: Analisi intercettatore Eve ---');
 
% 1. Eve su singolo frame (spettro FFT standard, non de-chirpato) del
%    pacchetto realmente trasmesso
Y_Eve_single = fftshift(fft(S_rx_bob));
Spec_Eve_single_dB = 20*log10(abs(Y_Eve_single) + eps);
 
% 2. Eve, media statistica su molti simboli equiprobabili (0/1) -> mostra
%    la cancellazione di covertness indipendentemente da quale bit sia
%    stato inviato in QUESTA specifica trasmissione (argomento statistico
%    sull'ensemble dei pacchetti, non sulla singola osservazione)
S_Eve_ensemble = 0.5 * S_rx_bit0 + 0.5 * S_rx_bit1;
Y_Eve_ensemble = fftshift(fft(S_Eve_ensemble));
Spec_Eve_ensemble_dB = 20*log10(abs(Y_Eve_ensemble) + eps);
 
%% ========================================================================
% FASE 8: FIGURA RIASSUNTIVA A 4 PANNELLI
% =========================================================================
disp('--- FASE 8: Generazione figura riassuntiva ---');

u_axis = linspace(-1, 1, N_campioni);
figure('Name', 'Quantum-Covert ISL Dinamico', 'Color', 'w', 'Position', [80, 50, 1050, 920]);

% (a) Geometria ISL dinamica: d(t) e v_rel(t)
subplot(4,1,1);
yyaxis left;
h_dist = plot(t_sec, d_t/1000, 'b-', 'LineWidth', 1.3);
ylabel('Distanza ISL d(t) [km]', 'Color', 'b', 'FontWeight', 'bold');
ax = gca; ax.YColor = 'b';
ylim([0, max(d_t/1000)*1.1]);

yyaxis right;
h_vrel = plot(t_sec, v_rel_t/1000, 'r-', 'LineWidth', 1.3);
hold on;
yline(0, 'k:', 'LineWidth', 0.8); % semplice riferimento visivo dello zero per la velocità
ylabel('Velocità Relativa LOS v_{rel}(t) [km/s]', 'Color', 'r', 'FontWeight', 'bold');
ax.YColor = 'r';

xlabel('Tempo di simulazione [s]');
title(['(a) Geometria Orbitale ISL: Avvicinamento e Crossing a t = ', ...
       num2str(t_sec(idx_cross),'%.1f'), ' s (d_{min} = ', num2str(d_min/1000,'%.1f'), ' km)']);
legend([h_dist, h_vrel], ...
       {'Distanza Inter-satellite d(t) [scala sx]', ...
        'Velocità Relativa LOS v_{rel}(t) [scala dx]'}, ...
       'Location', 'northeast', 'FontSize', 8);
grid on;

% (b) Link budget ottico standardizzato con fading PAT di Rayleigh
subplot(4,1,2);
h_p_ideal = plot(t_sec, 10*log10(P_rx_ideal_t*1e3), 'Color', [0.55 0.55 0.55], 'LineStyle', '--', 'LineWidth', 1.2);
hold on;
h_p_jit   = plot(t_sec, 10*log10(P_rx_t*1e3), 'Color', [0 0.45 0.75], 'LineWidth', 1.1);
xline(t_sec(idx_cross), 'k:', 'Crossing Fly-by', 'LineWidth', 1.0);
xlabel('Tempo di simulazione [s]'); 
ylabel('Potenza Ricevuta P_{rx}(t) [dBm]');
legend([h_p_ideal, h_p_jit], ...
       {'P_{rx} Ideale (FSPL + Efficienze Ottiche)', ...
        'P_{rx} Effettiva (Jitter di Rayleigh \sigma = 1.5 \murad)'}, ...
       'Location', 'best', 'FontSize', 8);
title('(b) Link Budget Ottico FSO Dinamico con Perdite PAT di Rayleigh');
grid on;

% (c) Spettro Eve: singolo shot vs media di ensemble
subplot(4,1,3);
h_eve_single = plot(u_axis, Spec_Eve_single_dB, 'Color', [0.25 0.25 0.25], 'LineWidth', 0.9);
hold on;
h_eve_ens    = plot(u_axis, Spec_Eve_ensemble_dB, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);
xlabel('Frequenza Spaziale Normalizzata u'); 
ylabel('Densità Spettrale [dB]');
legend([h_eve_single, h_eve_ens], ...
       {'Eve: Singolo Frame FFT (Segnale Covert + Fondo Gaussiano)', ...
        'Eve: Media d''Insieme Statistica (Elisione Quantistica \langle\Psi\rangle \rightarrow 0)'}, ...
       'Location', 'northeast', 'FontSize', 8);
title('(c) Analisi Intercettatore Eve: Invisibilità Energetica (SNR = -41 dB)');
xlim([-0.25 0.25]); grid on;

% (d) Spazio di decisione di Bob dopo de-chirp FrFT
subplot(4,1,4);
if bit_to_send == 0
    curve_color = 'b';
    bit_str = 'Bit 0 [\phi = 0]';
else
    curve_color = 'm';
    bit_str = 'Bit 1 [\phi = \pi]';
end
h_peak = plot(u_axis, Picco_Bob, curve_color, 'LineWidth', 1.3);
hold on;
h_soglia = yline(0, 'k--', 'LineWidth', 1.0);
h_det = plot(u_axis(peak_idx), peak_signed, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
text(u_axis(peak_idx), peak_signed, ['  Decodificato: Bit ', num2str(decoded_bit)], ...
     'FontWeight', 'bold', 'FontSize', 9, 'VerticalAlignment', 'bottom');
xlabel('u (Dominio Frazionario di Fourier FrFT)'); 
ylabel('Ampiezza Compressa');
legend([h_peak, h_soglia, h_det], ...
       {['Risposta Compressa FrFT (Tx: ', bit_str, ')'], ...
        'Soglia di Decisione Coerente (Ampiezza = 0)', ...
        ['Campione di Picco Rilevato (Decisione: Bit ', num2str(decoded_bit), ')']}, ...
       'Location', 'northeast', 'FontSize', 8);
title(['(d) Ricevitore di Bob (De-Chirp FrFT): Rilevazione Coerente (SNR_{post-FrFT} \approx +', ...
       num2str(SNR_post_Bob_dB,'%.1f'), ' dB)']);
xlim([-0.1 0.1]); grid on;

%% ========================================================================
% FASE 8b: MONITORAGGIO AMBIENTALE, VISIBILITÀ LOS E CONGIUNZIONE SOLARE
% =========================================================================
figure('Name', 'Analisi Geometria Orbitale e Ambiente Solare', 'Color', 'w', 'Position', [100, 100, 900, 600]);

% (1) Angolo di esclusione solare
subplot(2,1,1);
plot(t_sec, theta_sun_rx_deg, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.4);
hold on;
yline(theta_sun_crit, 'r--', ['Soglia Critica Abbagliamento (', num2str(theta_sun_crit), '°)'], 'LineWidth', 1.2);
xline(t_sec(idx_cross), 'k:', 'Crossing Fly-by', 'LineWidth', 1.0);
grid on; 
xlabel('Tempo di simulazione [s]'); 
ylabel('Angolo Sole-LOS [deg]');
title('Angolo di Esclusione Solare al Ricevitore Bob (\theta_{Sun})');
legend({'Angolo Sole-Puntamento Ricevitore', 'Limite Cieco Solare', 'Istante di Crossing'}, ...
    'Location', 'best', 'FontSize', 8);

% (2) Margine di clearance atmosferica
subplot(2,1,2);
plot(t_sec, los_clearance_margin / 1000, 'Color', [0 0.5 0], 'LineWidth', 1.4);
hold on;
yline(0, 'r--', 'Limite Atmosfera Densa (80 km)', 'LineWidth', 1.2);
xline(t_sec(idx_cross), 'k:', 'Crossing Fly-by', 'LineWidth', 1.0);
grid on; 
xlabel('Tempo di simulazione [s]'); 
ylabel('Margine di Clearance [km]');
title('Margine di Visibilità Ottica Intersatellitare (LOS Clearance)');
legend({'Clearance sopra l''Atmosfera', 'Soglia di Occlusione Terrestre', 'Istante di Crossing'}, ...
    'Location', 'best', 'FontSize', 8);

disp('--- TUTTI I GRAFICI GENERATI CON SUCCESSO! ---');

disp('--- SIMULAZIONE DINAMICA COMPLETATA CON SUCCESSO! ---');