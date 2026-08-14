% =========================================================================
% SIMULAZIONE DINAMICA LINK INTERSATELLITARE QUANTUM-COVERT (LEO-LEO ISL)
% Estensione orbitale con Satellite Communications Toolbox
% Modello: Ione 27Al+ | Doppler Dinamico & FrFT De-chirping
% =========================================================================
%
% v2 - correzioni rispetto alla prima estensione dinamica:
%   1) Geometria orbitale: ricerca automatica (coarse->fine) della fase
%      iniziale di Bob per ottenere un vero fly-by ISL (v_rel cambia segno
%      entro la finestra) a una distanza di crossing target ~300 km,
%      invece di un "avvicinamento" troncato ai bordi della simulazione.
%   2) Pointing jitter: la divergenza di fascio usata per il pointing loss
%      e' un parametro di progetto realistico (allargata oltre il limite
%      di diffrazione), non piu' comparabile per ordine di grandezza a
%      sigma_jitter: evita le attenuazioni patologiche (-1000+ dB) viste
%      nella v1.
%   3) Catena covert: si torna a trasmettere UN SOLO simbolo selezionabile
%      a mano (bit_to_send = 0 oppure 1), generato e ricevuto con le
%      IDENTICHE formule del chirp del modello analitico originale
%      (stesso c_chirp lato Tx e lato Rx, nessuna perturbazione aggiuntiva
%      che disallinei il matched filter). Il fattore di scala di
%      covertness torna quello originale, scorporato dal link budget
%      dinamico (che nella v1 poteva far collassare l'ampiezza del
%      simbolo in coincidenza di un fading di puntamento). Risultato:
%      Bob ottiene un picco FrFT netto e correttamente decodificabile.
%
% Richiede: Satellite Communications Toolbox (satelliteScenario, satellite,
%           states). Testato concettualmente per R2022a+.
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
% FASE 1: GEOMETRIA ORBITALE E RICERCA AUTOMATICA DEL CROSSING ISL
% =========================================================================
disp('--- FASE 1: Ricerca deterministica della geometria di crossing ISL ---');
 
startTime  = datetime(2026,8,14,10,0,0);
stopTime   = startTime + minutes(15);
sampleTime = 0.1;                       % Risoluzione temporale finale (s)
T_window   = seconds(stopTime - startTime);
t_target_middle = T_window / 2;         % Centro della finestra (450 s): qui vogliamo il crossing
 
altitude_orbit = 800e3;                 % Quota (m)
sma            = Re + altitude_orbit;   % Semiasse maggiore (m)
ecc            = 0;                     % Orbita quasi-circolare
incl_deg       = 60;                    % Inclinazione comune (deg)
argPeri_deg    = 0;                     % Argomento del perigeo (orbita circolare)
 
% Offset di RAAN per un angolo di incrocio tra i piani orbitali di 60°,
% dato incl = 60° per entrambi: cos(theta) = cos(i)^2 + sin(i)^2*cos(dRAAN)
theta_cross_deg = 60;
cos_dRAAN = (cosd(theta_cross_deg) - cosd(incl_deg)^2) / sind(incl_deg)^2;
dRAAN_deg = acosd(cos_dRAAN);
raan_Alice = 0;
raan_Bob   = dRAAN_deg;
 
mu_n      = sqrt(mu_earth / sma^3);     % Moto medio (rad/s), orbita circolare
T_orbit   = 2*pi / mu_n;                % Periodo orbitale completo (s), ~100 min a 800 km
incl_rad  = deg2rad(incl_deg);
raanA_rad = deg2rad(raan_Alice);
raanB_rad = deg2rad(raan_Bob);
 
% Propagatore analitico di posizione ECI per orbita circolare (u = arg. di
% latitudine = trueAnomaly per ecc=0). Usato SOLO per la ricerca rapida
% della fase: l'estrazione dinamica finale userà satelliteScenario/states.
posCircular = @(u0, raan_rad, t) sma * [ ...
    cos(raan_rad).*cos(u0 + mu_n*t) - sin(raan_rad).*sin(u0 + mu_n*t).*cos(incl_rad); ...
    sin(raan_rad).*cos(u0 + mu_n*t) + cos(raan_rad).*sin(u0 + mu_n*t).*cos(incl_rad); ...
    sin(u0 + mu_n*t).*sin(incl_rad) ];
 
target_distance = 300e3;                % Distanza ISL target al crossing (m)
 
% -------------------------------------------------------------------
% PASSO A: trovare la FASE RELATIVA Delta = u0_Bob - u0_Alice tale che il
% minimo GLOBALE di d(t), cercato su un intero periodo orbitale (non sulla
% finestra troncata di 900 s!), valga ~300 km. Il minimo globale dipende
% solo da Delta, non dalla fase assoluta: si fissa trueAnom_Alice_ref = 0
% come riferimento arbitrario per questa ricerca.
% -------------------------------------------------------------------
t_full_coarse = 0:2:T_orbit;
posA_full_coarse = posCircular(0, raanA_rad, t_full_coarse);
u0_candidates = deg2rad(0:0.5:359.5);
 
best_obj = Inf; best_u0 = 0; best_t_period = 0;
for k = 1:numel(u0_candidates)
    posB_k = posCircular(u0_candidates(k), raanB_rad, t_full_coarse);
    d_k = vecnorm(posB_k - posA_full_coarse, 2, 1);
    [dmin_k, idx_k] = min(d_k);
    obj_k = abs(dmin_k - target_distance);
    if obj_k < best_obj
        best_obj = obj_k; best_u0 = u0_candidates(k); best_t_period = t_full_coarse(idx_k);
    end
end
 
% Raffinamento locale di Delta (attorno al minimo grossolano, ancora sul
% periodo intero nell'intorno temporale del minimo trovato)
t_fine = max(0,best_t_period-10):0.02:min(T_orbit, best_t_period+10);
posA_fine = posCircular(0, raanA_rad, t_fine);
u0_fine = best_u0 + deg2rad(-1:0.002:1);
 
best_obj_f = Inf; best_u0_f = best_u0; best_t_f = best_t_period; best_dmin_f = Inf;
for k = 1:numel(u0_fine)
    posB_k = posCircular(u0_fine(k), raanB_rad, t_fine);
    d_k = vecnorm(posB_k - posA_fine, 2, 1);
    [dmin_k, idx_k] = min(d_k);
    obj_k = abs(dmin_k - target_distance);
    if obj_k < best_obj_f
        best_obj_f = obj_k; best_u0_f = u0_fine(k); best_t_f = t_fine(idx_k); best_dmin_f = dmin_k;
    end
end
 
Delta_phase    = best_u0_f;             % Fase relativa Bob-Alice (rad), fissata
d_min_achieved = best_dmin_f;           % Distanza minima effettivamente raggiunta (m)
 
% -------------------------------------------------------------------
% PASSO B: ri-centrare la FASE ASSOLUTA (mantenendo Delta_phase costante)
% in modo che il minimo trovato al Passo A cada esattamente a
% t = t_target_middle (centro della finestra di simulazione di 15 min),
% cosi' che il crossing (e l'inversione di segno di v_rel) sia visibile
% comodamente all'interno della finestra osservata.
% -------------------------------------------------------------------
u_Alice_at_min_ref = mu_n * best_t_f;   % u di Alice al minimo, nel riferimento trueAnom_Alice_ref=0
 
trueAnom_Alice = mod(rad2deg(u_Alice_at_min_ref) - rad2deg(mu_n*t_target_middle), 360);
trueAnom_Bob   = mod(trueAnom_Alice + rad2deg(Delta_phase), 360);
 
disp(['Periodo orbitale: ', num2str(T_orbit/60,'%.2f'), ' min']);
disp(['Fase relativa Bob-Alice (Delta): ', num2str(rad2deg(Delta_phase),'%.4f'), ' deg']);
disp(['trueAnom_Alice = ', num2str(trueAnom_Alice,'%.4f'), ' deg, trueAnom_Bob = ', ...
      num2str(trueAnom_Bob,'%.4f'), ' deg']);
disp(['Crossing atteso a t ~ ', num2str(t_target_middle,'%.1f'), ' s, d_min stimato = ', ...
      num2str(d_min_achieved/1000,'%.1f'), ' km']);
 
%% ========================================================================
% FASE 2: SCENARIO DINAMICO CON SATELLITE COMMUNICATIONS TOOLBOX
% =========================================================================
disp('--- FASE 2: Propagazione dinamica con satelliteScenario ---');
 
sc = satelliteScenario(startTime, stopTime, sampleTime);
 
satAlice = satellite(sc, sma, ecc, incl_deg, raan_Alice, argPeri_deg, ...
    trueAnom_Alice, "Name", "Alice");
satBob   = satellite(sc, sma, ecc, incl_deg, raan_Bob, argPeri_deg, ...
    trueAnom_Bob, "Name", "Bob");
 
[posAlice, velAlice, tSamples] = states(satAlice, "CoordinateFrame", "inertial");
[posBob,   velBob,   ~]        = states(satBob,   "CoordinateFrame", "inertial");
 
t_sec = seconds(tSamples - tSamples(1));   % Tempo simulazione (s), riga 1xN
N_t   = numel(t_sec);
 
% Distanza ISL istantanea e ritardo di propagazione
deltaPos = posBob - posAlice;              % 3xN (m)
d_t      = vecnorm(deltaPos, 2, 1);        % Distanza ISL d(t) (m)
tau_t    = d_t / c_light;                  % Ritardo di propagazione tau(t) (s)
 
% Velocita' relativa lungo la linea di vista (LOS, range-rate)
uLOS     = deltaPos ./ d_t;                % Versore LOS Alice->Bob, 3xN
deltaVel = velBob - velAlice;              % 3xN (m/s)
v_rel_t  = sum(deltaVel .* uLOS, 1);       % v_rel(t) (m/s), >0 = allontanamento
 
% Doppler ottico non relativistico riferito alla transizione 27Al+
f_D_t    = -(v_rel_t / c_light) * f0_optical; % Shift Doppler f_D(t) (Hz)
dfD_dt_t = gradient(f_D_t, sampleTime);       % Doppler-rate df_D/dt (Hz/s)
 
[d_min, idx_cross] = min(d_t);
sign_change = any(diff(sign(v_rel_t)) ~= 0);
disp(['Distanza minima ISL: ', num2str(d_min/1000,'%.2f'), ' km al campione ', ...
      num2str(idx_cross), ' (t = ', num2str(t_sec(idx_cross),'%.1f'), ' s)']);
disp(['v_rel al crossing: ', num2str(v_rel_t(idx_cross)/1000,'%.3f'), ' km/s']);
if sign_change
    disp('Fly-by confermato: v_rel cambia segno entro la finestra di simulazione.');
else
    disp('ATTENZIONE: v_rel non cambia segno entro la finestra; ripetere la ricerca con finestra piu'' ampia.');
end
disp(['Doppler f_D al crossing: ', num2str(f_D_t(idx_cross)/1e6,'%.3f'), ' MHz']);
disp(['Doppler-rate al crossing: ', num2str(dfD_dt_t(idx_cross)/1e3,'%.3f'), ' kHz/s']);
 
%% ========================================================================
% FASE 3: LINK BUDGET OTTICO DINAMICO CON POINTING JITTER (PAT)
% =========================================================================
disp('--- FASE 3: Link budget ottico dinamico e jitter di puntamento ---');
 
D_tx        = 0.20;     % Apertura Tx Alice (m)
D_rx        = 0.30;     % Apertura Rx Bob (m)
P_tx_laser  = 100e-3;   % Potenza di pompaggio laser (W)
sigma_jitter = 1.5e-6;  % Deviazione standard jitter di puntamento (rad)
 
% Guadagni ottici dei telescopi (approssimazione ad apertura diffrattiva)
G_tx = (pi * D_tx / lambda_0)^2;
G_rx = (pi * D_rx / lambda_0)^2;
 
% Divergenza di fascio usata per il modello di pointing-loss: e' un
% parametro di PROGETTO del sistema (tipicamente allargato appositamente
% oltre il limite di diffrazione dell'apertura, cosi' da rilassare i
% requisiti di puntamento), non il limite di diffrazione stretto usato
% sopra per il guadagno d'antenna. Se si usasse qui il limite di
% diffrazione stretto (~1.3 microrad), sigma_jitter=1.5 microrad sarebbe
% dello stesso ordine di grandezza dell'intero fascio e produrrebbe fading
% patologici; il valore realistico tipico per ISL ottici e' O(10 microrad).
theta_div_beam = 8e-6;  % rad, divergenza di sistema assunta
 
% FSPL dinamico, funzione della distanza istantanea
FSPL_t = (4 * pi * d_t / lambda_0).^2;
 
% Errore di puntamento stocastico gaussiano indipendente su Alice e Bob
theta_err_tx = sigma_jitter * randn(1, N_t);
theta_err_rx = sigma_jitter * randn(1, N_t);
 
% Perdita di puntamento (modello Gaussiano del fascio, PAT loss)
L_point_tx = exp(-8 * (theta_err_tx ./ theta_div_beam).^2);
L_point_rx = exp(-8 * (theta_err_rx ./ theta_div_beam).^2);
L_point_t  = L_point_tx .* L_point_rx;
 
% Potenza ricevuta istantanea, con e senza degradazione da jitter
P_rx_ideal_t = P_tx_laser * (G_tx * G_rx) ./ FSPL_t;
P_rx_t       = P_rx_ideal_t .* L_point_t;
 
disp(['Potenza ricevuta media (con jitter): ', num2str(mean(P_rx_t),'%.3e'), ' W']);
disp(['Degradazione media da pointing jitter: ', ...
      num2str(10*log10(mean(L_point_t)),'%.2f'), ' dB']);
 
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

% (b) Link budget ottico dinamico con degradazione da pointing jitter
subplot(4,1,2);
h_p_ideal = plot(t_sec, 10*log10(P_rx_ideal_t*1e3), 'Color', [0.55 0.55 0.55], 'LineStyle', '--', 'LineWidth', 1.2);
hold on;
h_p_jit   = plot(t_sec, 10*log10(P_rx_t*1e3), 'Color', [0 0.45 0.75], 'LineWidth', 1.1);
xlabel('Tempo di simulazione [s]'); 
ylabel('Potenza Ricevuta P_{rx}(t) [dBm]');
legend([h_p_ideal, h_p_jit], ...
       {'P_{rx} Ideale (Attenuazione FSPL pura)', ...
        'P_{rx} Effettiva (con Jitter di Puntamento \sigma = 1.5 \murad)'}, ...
       'Location', 'best', 'FontSize', 8);
title('(b) Link Budget Ottico Dinamico e Fading da Pointing Jitter (PAT)');
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

disp('--- SIMULAZIONE DINAMICA COMPLETATA CON SUCCESSO! ---');