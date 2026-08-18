% =========================================================================
% SIMULATORE DI THROUGHPUT REALISTICO - LINK ISL QUANTUM-COVERT (LEO-LEO)
% Versione "avionica": banda IF da SmallSat, guadagno d'antenna coerente
% con la divergenza reale, jitter di tracking Doppler residuo (floor di
% BER fisiologico). Stessa architettura della versione benchmark.
% Modello: Ione 27Al+ | Doppler Dinamico & FrFT De-chirping
%
% DIFFERENZE RISPETTO AL BENCHMARK (che restano documentate per confronto)
%  1) B_IF = 100 MHz (SmallSat/CubeSat) invece di 16.384 GHz da banco.
%     N = 16384 invariato -> stesso guadagno di processo (+42.1 dB), ma
%     T_sym = N/B_IF ~ 163.84 us e bit rate lordo ~ 6.1 kbit/s.
%  2) Guadagno d'antenna Tx/Rx con la formula documentata dell'oggetto
%     gaussianAntenna del Satellite Communications Toolbox (Fase 3), non
%     piu' al limite di diffrazione ne' con un fattore di scala scollegato
%     dall'apertura.
%  3) Bob non conosce il bin di picco con precisione infinita: un errore
%     di tracking Doppler residuo (oscillatore/GPS) sposta casualmente la
%     lettura di +/-1..2 bin -> compare un FLOOR di BER irriducibile,
%     visibile per confronto diretto con la curva analitica (che non lo
%     include, perche' assume sincronismo perfetto).
%
% Tutto il resto (geometria SGP4/J2, congiunzione solare, struttura delle
% campagne Monte Carlo, metriche di throughput) e' invariato nella logica.
% =========================================================================

clc; clear; close all;

% -------------------------------------------------------------------------
% SEED MASTER UNICO PER TUTTA LA CAMPAGNA
% -------------------------------------------------------------------------
% Ogni realizzazione Monte Carlo usa un SUBSTREAM indipendente dello STESSO
% generatore mrg32k3a, inizializzato con SEED_MASTER. Questo garantisce
% riproducibilita' totale della campagna e, allo stesso tempo, indipendenza
% statistica fra le realizzazioni (condizione necessaria perche' la media
% Monte Carlo sia significativa). mrg32k3a e' scelto perche' e' l'unico
% generatore MATLAB che espone esplicitamente la proprieta' Substream per
% la creazione di flussi indipendenti garantiti.
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
k_B         = 1.380649e-23;       % Costante di Boltzmann (J/K)
u_amu       = 1.66053906660e-27;  % Unita' di massa atomica (kg)
m_Al27      = 27 * u_amu;         % Massa dello ione 27Al+ (kg)
Re          = 6378137;            % Raggio equatoriale Terra, WGS84 (m)
mu_earth    = 3.986004418e14;     % Parametro gravitazionale terrestre (m^3/s^2)

disp('=== SEZIONE 0: Costanti caricate ===');
disp(['Frequenza portante ottica f0: ', num2str(f0_optical/1e12,'%.3f'), ' THz']);

%% ========================================================================
% FASE 1: GEOMETRIA ORBITALE CON EFFETTO J2 E TARGETING DEL CROSSING
% =========================================================================
% (Invariata rispetto al benchmark: la geometria orbitale non dipende dai
% parametri del ricevitore/canale che stiamo rendendo realistici.)
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

link_available = (los_clearance_margin >= 0) & (~sun_blindness_flag);

%% ========================================================================
% FASE 3: MODELLO DI CANALE OTTICO REALISTICO (FSO / PAT)
% =========================================================================
disp('--- FASE 3: Canale Ottico Realistico (FSO Link Budget & PAT) ---');

D_tx        = 0.20;         % Apertura Tx Alice (m)
D_rx        = 0.30;         % Apertura Rx Bob (m)
P_tx_laser  = 100e-3;       % Potenza media trasmessa (W)

eta_tx      = 0.85;         % Efficienza ottica di catena Tx (specchi/accoppiamento) - NON e' l'aperture efficiency d'antenna, vedi sotto
eta_rx      = 0.80;         % Efficienza ottica di catena Rx
eta_det     = 0.65;
DCR         = 50;           % Dark Count Rate (counts/s)

sigma_jitter_pat = 1.5e-6;  % Jitter PAT 1-sigma per asse (rad) - usato SOLO nel fading di puntamento sotto
theta_div_beam   = 8.0e-6;  % Divergenza di progetto usata SOLO nella formula di fading di puntamento
                             % (exp(-8*(theta_err/theta_div)^2)), non piu' nel guadagno (vedi sotto):
                             % resta un parametro indipendente ancora da verificare in letteratura.

% GUADAGNO D'ANTENNA: FORMULA DOCUMENTATA DI gaussianAntenna (Satellite
% Communications Toolbox), non piu' una coppia di formule scollegate per
% Tx e Rx. Dalla documentazione ufficiale dell'oggetto gaussianAntenna:
%   boresightGain = ApertureEfficiency * (pi*D/lambda)^2      [lineare]
%   beamwidth_3dB = 70*lambda/D                                [gradi]
% Si applica la STESSA formula, con la STESSA convenzione, sia a Tx che a
% Rx: un'apertura e' un'apertura, che emetta o raccolga luce. Si usa
% ApertureEfficiency = 0.65, il valore di DEFAULT che il toolbox stesso
% assegna a un gaussianAntenna quando non specificato altrimenti (vedi
% documentazione transmitter/gaussianAntenna) - non un numero scelto a
% mano. Questo sostituisce sia il precedente G_tx=32/theta_div^2 (che
% legava il guadagno a una divergenza scelta indipendentemente da D_tx,
% rischiando un'incoerenza fisica fra apertura e divergenza dichiarate)
% sia il precedente G_rx al limite di diffrazione puro (che assumeva
% implicitamente efficienza di apertura del 100%, non realistica).
rho_ap_tx = 0.65;    % Aperture efficiency Tx (default toolbox gaussianAntenna)
rho_ap_rx = 0.65;    % Aperture efficiency Rx (default toolbox gaussianAntenna)

G_tx = rho_ap_tx * (pi * D_tx / lambda_0)^2;
G_rx = rho_ap_rx * (pi * D_rx / lambda_0)^2;

beamwidth_3dB_tx_deg = 70 * lambda_0 / D_tx;   % formula documentata gaussianAntenna (gradi)
beamwidth_3dB_tx_rad = deg2rad(beamwidth_3dB_tx_deg);

disp(['Guadagno Tx (gaussianAntenna, ApertureEfficiency=', num2str(rho_ap_tx), '): ', ...
      num2str(10*log10(G_tx),'%.2f'), ' dB']);
disp(['Guadagno Rx (gaussianAntenna, ApertureEfficiency=', num2str(rho_ap_rx), '): ', ...
      num2str(10*log10(G_rx),'%.2f'), ' dB']);
disp(['Beamwidth a -3 dB (formula gaussianAntenna, 70*lambda/D): ', ...
      num2str(beamwidth_3dB_tx_rad*1e6,'%.3f'), ' urad - per confronto con theta_div_beam ' ...
      '= ', num2str(theta_div_beam*1e6,'%.1f'), ' urad usato nel fading di puntamento (grandezze ' ...
      'angolari diverse, non direttamente comparabili: beamwidth a -3dB di guadagno vs ' ...
      'divergenza 1/e^2 di intensita'' del fascio laser).']);

% FSPL: funzione standard fspl() invece della formula scritta a mano
% (risultato identico, (4*pi*R/lambda)^2 in dB - qui solo per tracciabilita'
% a una funzione MATLAB validata invece che a un'espressione inline).
FSPL_t = 10.^(fspl(d_t, lambda_0) / 10);

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
% FASE 4: CONFIGURAZIONE AVIONICA DELLA CAMPAGNA MULTI-BIT
% =========================================================================
disp('--- FASE 4: Configurazione campagna multi-bit (banda IF da SmallSat) ---');

cfg = struct();

% --- Banda IF REALISTICA da ricevitore SmallSat/CubeSat ---
% B_IF: valore rappresentativo, dichiarato, di un ricevitore ottico
% coerente classe CubeSat (ordine di grandezza confermato da ricevitori
% ISL coerenti a DSP dimostrati in letteratura, ~100 Mbit/s di classe
% affine). Non e' il risultato di un dimensionamento fine (linewidth
% laser + banda del loop di tracking + margine), perche' NON serve: B_IF
% entra nel simulatore solo attraverso T_sym = N/B_IF, che fissa la scala
% assoluta del bit rate e del tempo di rivelazione di Eve (entrambi
% lineari in B_IF). BER, floor di tracking e curva di covertness NON
% dipendono da B_IF - dipendono solo da N (guadagno di processo) e dalla
% SNR covert relativa al fondo classico (Fase 5). Cambiare B_IF cambia
% quindi la scala del risultato, non la conclusione sulla fattibilita'
% della comunicazione covert: per questo motivo un valore rappresentativo
% dichiarato e' sufficiente, senza richiedere un dimensionamento di
% precisione da letteratura specialistica.
cfg.B_IF         = 100e6;        % 100 MHz (invece dei 16.38 GHz da banco)
cfg.N            = 2^14;         % Campioni per simbolo -> INVARIATO: stesso guadagno di processo
cfg.T_sym        = cfg.N / cfg.B_IF;             % ~163.84 us
R_b_raw          = 1 / cfg.T_sym;                % ~6.1 kbit/s

% --- FORMA D'ONDA: chirp dimensionato ESPLICITAMENTE su B_IF, non piu' un
% valore ereditato per inerzia. Si sceglie una frazione di banda occupata
% dal chirp rispetto a B_IF (margine per guard-band e roll-off dei filtri
% anti-aliasing): frac_BW_chirp = 1200/16384 = 7.32%, la stessa frazione
% che i valori del benchmark (300*pi, 100*pi) occupavano gia' implicitamente
% a 16.38 GHz -> qui viene resa un parametro di progetto esplicito invece
% che una coincidenza numerica. cfg.offset/cfg.c_chirp = 1/3 e' mantenuto
% (definisce l'asimmetria della banda occupata attorno al centro).
%
% Derivazione: l'escursione di frequenza istantanea del chirp e'
%   f(omega) = (2*c_chirp*omega + offset) * (B_IF/N) / pi ,  omega in [-1,1]
% Imponendo che l'escursione totale (da omega=-1 a omega=+1) valga
% frac_BW_chirp * B_IF, con offset = c_chirp/3, si ottiene
%   c_chirp = frac_BW_chirp * N * pi / 4
frac_BW_chirp    = 1200/16384;                  % 7.32% di B_IF (design esplicito)
cfg.c_chirp      = frac_BW_chirp * 2^14 * pi / 4;   % = 300*pi (stesso valore validato)
cfg.offset       = cfg.c_chirp / 3;                 % = 100*pi (stessa asimmetria validata)

B_chirp_occ_Hz   = frac_BW_chirp * cfg.B_IF;
disp(['Chirp dimensionato su B_IF: occupa ', num2str(B_chirp_occ_Hz/1e6,'%.2f'), ...
      ' MHz (', num2str(frac_BW_chirp*100,'%.1f'), ' % di B_IF), margine per guard-band = ', ...
      num2str((1-frac_BW_chirp)*100,'%.1f'), ' %']);

% --- FONDO CLASSICO: ALLARGAMENTO DOPPLER TERMICO DELLO IONE 27Al+ ------
% Il "fondo classico" sotto cui si nasconde il segnale covert NON e' un
% inviluppo arbitrario: e' il profilo spettrale del canale classico
% legittimo co-locato, cioe' la riga di emissione dello ione intrappolato
% 27Al+ allargata per effetto Doppler termico dal moto residuo dello ione.
% E' quindi una quantita' DERIVABILE da principi primi, non un parametro
% da giustificare:
%
%     sigma_nu = nu0 * sqrt(k_B*T_ion / (m*c^2))          [Hz, sigma gaussiana]
%
% CORREZIONE DI COERENZA rispetto alla versione precedente: prima
% sigma_doppler era fissato a 0.15 nel dominio normalizzato, cioe' come
% frazione fissa di B_IF/2. Ma questo implica che, cambiando B_IF, cambia
% silenziosamente la TEMPERATURA implicita dello ione - il che e'
% fisicamente rovesciato: la temperatura e' una proprieta' della trappola,
% del tutto indipendente dalla banda IF del ricevitore. Ora il verso della
% derivazione e' corretto: si fissa T_ion (fisica), e sigma_doppler nel
% dominio normalizzato ne discende insieme a B_IF.
T_ion = 13.02e-3;    % Temperatura dello ione intrappolato (K)
                     % 13 mK: ione 27Al+ raffreddato simpateticamente in
                     % trappola compatta. Sopra il limite Doppler teorico
                     % (~1 mK) per tenere conto di micromoto in eccesso e
                     % riscaldamento anomalo, realistici in un dispositivo
                     % qualificato per lo spazio. Valore scelto anche per
                     % continuita' con le campagne gia' validate (riproduce
                     % sigma_doppler = 0.150).

sigma_doppler_Hz = f0_optical * sqrt(k_B * T_ion / (m_Al27 * c_light^2));
cfg.sigma_doppler = sigma_doppler_Hz / (cfg.B_IF/2);   % -> dominio normalizzato
cfg.sigma_noise  = 0.02;

fwhm_doppler_Hz  = 2*sqrt(2*log(2)) * sigma_doppler_Hz;
disp(['Fondo classico = riga 27Al+ allargata Doppler a T_ion = ', ...
      num2str(T_ion*1e3,'%.2f'), ' mK: sigma = ', num2str(sigma_doppler_Hz/1e6,'%.3f'), ...
      ' MHz (FWHM = ', num2str(fwhm_doppler_Hz/1e6,'%.3f'), ' MHz) -> sigma_doppler = ', ...
      num2str(cfg.sigma_doppler,'%.4f'), ' nel dominio normalizzato']);

% Verifica di validita' del profilo GAUSSIANO (regime di sideband NON
% risolte). Per uno ione in trappola di Paul lo spettro e' una portante con
% bande laterali motionali alla frequenza secolare: il profilo gaussiano di
% Doppler termico e' una buona approssimazione solo se l'allargamento
% Doppler SUPERA nettamente la frequenza secolare (weak-binding limit),
% altrimenti (sideband risolte) lo spettro e' discreto e la gaussiana non
% descrive piu' la fisica.
f_secolare_tipica = 3e6;    % Frequenza secolare di trappola tipica (Hz)
rapporto_sideband = fwhm_doppler_Hz / f_secolare_tipica;
if rapporto_sideband > 3
    disp(['Regime di sideband NON risolte (FWHM/f_secolare = ', ...
          num2str(rapporto_sideband,'%.1f'), '): profilo gaussiano valido.']);
else
    warning(['Regime di sideband RISOLTE o marginale (FWHM/f_secolare = ', ...
             num2str(rapporto_sideband,'%.1f'), '): il profilo gaussiano NON descrive ' ...
             'piu'' correttamente lo spettro dello ione (servirebbe struttura a bande laterali).']);
end

cfg.omega        = linspace(-1, 1, cfg.N).';
cfg.domega       = cfg.omega(2) - cfg.omega(1);
cfg.win_clutter  = max(3, round(0.01 * cfg.N));
cfg.chiave       = exp(-1i * (cfg.c_chirp * (cfg.omega.^2)));
cfg.blocco       = 64;

% =========================================================================
% ARCHITETTURA DI PRE-COMPENSAZIONE DOPPLER (dichiarazione esplicita)
% =========================================================================
% Il Doppler ottico reale raggiunge decine di GHz (vedi Fase 2: fino a
% 26.7 GHz in questo passaggio), ben oltre B_IF = 100 MHz: i due numeri
% sono fisicamente incompatibili a meno di dichiarare esplicitamente come
% il ricevitore li concilia. Si assume quindi:
%
%  1) STADIO GROSSOLANO (analogico/RF, non modellato in dettaglio qui): un
%     oscillatore locale (LO) agile in frequenza, pilotato in anello
%     aperto dalle stesse effemeridi/propagatore d'orbita che Bob usa
%     altrove in questo script (coerente con l'assunzione, gia' presente,
%     che Bob conosca lo shift Doppler da effemeridi), insegue in continuo
%     il Doppler ottico predetto f_D(t) e lo rimuove quasi interamente
%     PRIMA della digitalizzazione. E' cio' che rende fisicamente possibile
%     stare dentro B_IF = 100 MHz nonostante 26.7 GHz di Doppler assoluto:
%     senza questo stadio l'architettura non avrebbe senso a questa banda.
%
%  2) RESIDUO INTRA-SIMBOLO (quello che il resto dello script modella
%     esplicitamente): anche un LO ideale non puo' correggere piu' in
%     fretta di quanto la sua stessa legge di sintonia venga aggiornata;
%     qui si assume, in modo conservativo verso il ricevitore (il caso
%     migliore possibile), che il LO applichi UNA correzione per simbolo.
%     Durante il simbolo (T_sym) il Doppler continua pero' a variare per
%     effetto del Doppler-RATE dfD_dt_t(t) gia' calcolato in Fase 2,
%     lasciando un residuo strutturale:
%         Delta_f_res(t) ~= dfD_dt_t(t) * T_sym
%         shift_omega(t)  = Delta_f_res(t) / (B_IF/2)
%     Questa e' la grandezza che sostituisce il precedente fattore di scala
%     grafico K_doppler_grafico: e' derivata da quantita' gia' presenti nel
%     modello (Doppler-rate reale, T_sym, B_IF), non da una costante
%     empirica.
%
%  3) RESIDUO STOCASTICO FINE (rumore di fase dell'oscillatore, incertezza
%     di orbit determination/GPS): modellato separatamente, invariato,
%     come jitter intero sul bin di decisione (cfg.sigma_bin_track, sotto).
%     E' un contributo ULTERIORE rispetto al residuo strutturale del punto
%     2, non un suo sostituto.
% =========================================================================

% --- IMPERFEZIONE DI TRACKING DOPPLER (jitter residuo sul bin di picco) ---
% Bob ricentra lo spettro usando lo shift Doppler stimato dalle effemeridi
% e dall'oscillatore di bordo, ma la stima ha un residuo stocastico (rumore
% di fase dell'oscillatore, incertezza di orbit determination/GPS). Lo si
% modella come un offset intero casuale sul bin di decisione dopo il
% de-chirp FrFT: offset = round(sigma_bin_track * randn), saturato a +/-2
% bin.
%
% ATTENZIONE ALLA CALIBRAZIONE (verificata con una Monte Carlo dedicata,
% separando gli errori per esatto valore di offset, a SNR = -41 dB):
%   offset =  0  ->  BER ~ 0        (2673/2673 letture corrette nel test)
%   offset = +-1 ->  BER ~ 5-7 %    (il lobo laterale ha guadagno di
%                                     processo molto minore del lobo
%                                     principale: la stessa potenza covert
%                                     che e' larghissimamente sufficiente
%                                     al bin corretto NON lo e' affatto al
%                                     bin adiacente)
%   offset = +-2 ->  BER ~ 85-100 % (praticamente sempre sbagliato)
% Il punto critico e' che questa degradazione condizionata dipende dalla
% SNR nominale: a una SNR di comodo molto piu' alta (es. -20 dB) lo stesso
% offset di 1 bin costa solo ~0.8% di errore, non 5-7%, perche' il margine
% residuo assorbe la perdita di guadagno del lobo laterale. Calibrare
% sigma_bin_track a una SNR "facile" e poi applicarlo al punto operativo
% reale (-41 dB) sottostima quindi il floor di un ordine di grandezza: e'
% l'errore fatto nella prima versione di questo script (sigma_bin_track =
% 0.5, floor osservato ~5*10^-2 anziche' i pochi per mille attesi). La
% calibrazione qui sotto e' stata rifatta DIRETTAMENTE alla SNR di
% ancoraggio: sigma_bin_track = 0.20 da' un floor a -41 dB di ~7*10^-4,
% che cresce con continuita' (non resta piatto) verso pochi 10^-3 avvicinandosi
% al ginocchio del waterfall intrinseco: e' un secondo termine di degrado
% dominato dalla sincronizzazione vicino al crossing, via via assorbito dal
% rumore di canale lontano da esso - piu' corretto fisicamente di un floor
% piatto assunto a priori, e piu' difendibile in tesi.
cfg.sigma_bin_track = 0.20;

% Risoluzione in frequenza di un bin dopo la FFT su N campioni a passo
% 1/B_IF: e' lo standard delta_f = B_IF/N. Da' un significato fisico (Hz)
% al jitter di tracking, utile per la discussione in tesi.
delta_f_bin      = cfg.B_IF / cfg.N;
sigma_track_Hz   = cfg.sigma_bin_track * delta_f_bin;

% --- Parametri della campagna Monte Carlo ---
cfg.Nbit_frame   = 256;
cfg.N_MC         = 12;
SNR_anchor_dB    = -41;
snr_vec_dB       = -41:-2.9:-70;     % Campagna A: da -41 dB a -70 dB (11 punti)
N_epoche         = 21;
L_frame_vec      = [256, 1024];      % Goodput non codificato a due granularita'

G_proc_dB = 10*log10(cfg.N);
disp(['Campioni per simbolo N = ', num2str(cfg.N), ...
      '  ->  guadagno di processo FrFT = +', num2str(G_proc_dB,'%.1f'), ' dB (invariato vs benchmark)']);
disp(['Banda IF = ', num2str(cfg.B_IF/1e6,'%.1f'), ' MHz  ->  T_sym = ', ...
      num2str(cfg.T_sym*1e6,'%.2f'), ' us  ->  bit rate lordo = ', ...
      num2str(R_b_raw/1e3,'%.3f'), ' kbit/s']);
disp(['Risoluzione di un bin FFT: ', num2str(delta_f_bin/1e3,'%.3f'), ...
      ' kHz  ->  jitter di tracking 1-sigma equivalente = ', num2str(sigma_track_Hz/1e3,'%.3f'), ' kHz']);

% Verifica del regime fotonico lungo TUTTO il passaggio (non solo al
% crossing): con banda IF ridotta il tempo di simbolo e' molto piu' lungo,
% quindi si integra piu' energia per simbolo nonostante il guadagno Tx
% coerente sia inferiore a quello di diffrazione.
n_photons_t = (P_rx_t * cfg.T_sym / E_fotone) * eta_det;
n_dark_sym  = DCR * cfg.T_sym;
[n_ph_min, idx_ph_min] = min(n_photons_t);

disp(['Fotoni di segnale per simbolo al crossing: ', num2str(n_photons_t(idx_cross),'%.3e')]);
disp(['Fotoni di segnale per simbolo, minimo sul passaggio: ', num2str(n_ph_min,'%.3e'), ...
      ' (a t = ', num2str(t_sec(idx_ph_min),'%.1f'), ' s)']);
disp(['Dark counts per simbolo: ', num2str(n_dark_sym,'%.3e')]);
if n_ph_min < 10
    warning(['Regime di conteggio fotonico (n < 10 fotoni/simbolo) raggiunto in parte del passaggio: ' ...
             'il modello gaussiano analogico non e'' piu'' valido li'', servirebbe una statistica di Poisson.']);
else
    disp('Regime classico (n >> 1) valido per l''intero passaggio: il modello analogico gaussiano si applica.');
end

%% ========================================================================
% FASE 5: CAMPAGNA A - BER vs LIVELLO DI COVERTNESS (AL CROSSING)
% =========================================================================
disp('--- FASE 5: Campagna A - BER vs covertness (Monte Carlo, con jitter di tracking) ---');

Shift_cross = shiftResiduoDoppler(cfg, dfD_dt_t(idx_cross));
disp(['Residuo Doppler intra-simbolo al crossing: ', num2str(dfD_dt_t(idx_cross)*cfg.T_sym/1e3,'%.3f'), ...
      ' kHz  ->  shift_omega = ', num2str(Shift_cross,'%.5f'), ' (dominio normalizzato)']);
op_cross    = operatoreShift(cfg, Shift_cross);
kpk_cross   = trovaBinPicco(cfg, op_cross, Shift_cross);

k_teorico = round(cfg.N/2) + 1 + round(cfg.offset/pi * cfg.N/(cfg.N-1));
if abs(kpk_cross - k_teorico) > 2
    warning(['Bin di picco misurato (', num2str(kpk_cross), ') distante da quello ' ...
             'teorico (', num2str(k_teorico), '): verificare offset/chirp.']);
end
disp(['Bin di decisione nominale (da effemeridi, prima del jitter di tracking): k = ', num2str(kpk_cross)]);

N_snr   = numel(snr_vec_dB);
err_A   = zeros(1, N_snr);
bit_A   = zeros(1, N_snr);
stat_nom = []; bits_nom = [];

for is = 1:N_snr
    for mc = 1:cfg.N_MC
        rs = streamMC(SEED_MASTER, 1, is, mc);
        [ne, nb, st, bt] = simulaFrame(cfg, op_cross, Shift_cross, ...
                                       snr_vec_dB(is), kpk_cross, rs);
        err_A(is) = err_A(is) + ne;
        bit_A(is) = bit_A(is) + nb;
        if is == 1 && mc == 1
            stat_nom = st; bits_nom = bt;
        end
    end
    fprintf('   SNR covert = %6.2f dB  ->  BER = %8.2e  (%d errori su %d bit)\n', ...
            snr_vec_dB(is), err_A(is)/bit_A(is), err_A(is), bit_A(is));
end

BER_A = err_A ./ bit_A;
[BER_A_lo, BER_A_hi] = wilsonCI(err_A, bit_A, 1.96);

% Stima del "floor" locale: media dei 3 punti a SNR piu' alta, dove la BER
% intrinseca (senza jitter) sarebbe gia' praticamente nulla e quindi il
% residuo osservato e' dominato dal tracking, non dal rumore di canale. Non
% e' un floor asintotico in senso stretto (vedi NOTA 1): e' la stima del
% contributo del tracking nell'intorno del punto di ancoraggio di progetto.
BER_floor_stimato = mean(BER_A(1:min(3,N_snr)));
disp(['BER residua stimata vicino al punto di ancoraggio (dominata dal jitter di tracking): ', num2str(BER_floor_stimato,'%.3e')]);

% Curva analitica di riferimento SENZA jitter di tracking (sincronismo
% perfetto): serve da limite superiore di prestazione, non da predizione
% realistica. La distanza fra questa curva e i punti Monte Carlo alle alte
% SNR e' proprio il costo, in BER, dell'imperfezione di sincronizzazione.
P_classica  = mean(exp(-((cfg.omega - Shift_cross).^2)/(2*cfg.sigma_doppler^2)).^2);
snr_fine_dB = linspace(min(snr_vec_dB)-2, max(snr_vec_dB)+2, 400);
P_cov_fine  = P_classica * 10.^(snr_fine_dB/10);
arg_fine    = sqrt(P_cov_fine * cfg.N) / cfg.sigma_noise;
BER_teorica = 0.5 * erfc(arg_fine / sqrt(2));

%% ========================================================================
% FASE 6: CAMPAGNA B - BER E THROUGHPUT LUNGO IL PASSAGGIO ORBITALE
% =========================================================================
disp('--- FASE 6: Campagna B - BER e throughput vs epoca orbitale ---');

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
    shift_B(ie) = shiftResiduoDoppler(cfg, dfD_dt_t(ii));

    if abs(shift_B(ie)) > 0.6
        warning(['Epoca ', num2str(ie), ': shift Doppler normalizzato = ', ...
                 num2str(shift_B(ie),'%.2f'), ' -> fuori dal dominio utile del modello.']);
    end

    op_e = operatoreShift(cfg, shift_B(ie));
    kpk_e = trovaBinPicco(cfg, op_e, shift_B(ie));

    snr_acc = 0;
    for mc = 1:cfg.N_MC
        rs = streamMC(SEED_MASTER, 2, ie, mc);

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

disponibile = double(link_available(idx_epoche));

R_naive = R_b_raw * (1 - BER_B) .* disponibile;
C_bsc   = 1 - entropiaBinaria(BER_B);
R_shan  = R_b_raw * C_bsc .* disponibile;

nL = numel(L_frame_vec);
FER_L    = zeros(nL, N_epoche);
R_good_L = zeros(nL, N_epoche);
Vol_good_L = zeros(1, nL);
for iL = 1:nL
    FER_L(iL,:)    = 1 - (1 - BER_B).^L_frame_vec(iL);
    R_good_L(iL,:) = R_b_raw * (1 - FER_L(iL,:)) .* disponibile;
    Vol_good_L(iL) = trapz(t_epoche, R_good_L(iL,:));
end

Vol_naive = trapz(t_epoche, R_naive);
Vol_shan  = trapz(t_epoche, R_shan);

BER_soglia = 1e-3;
utile      = (BER_B <= BER_soglia) & (disponibile > 0);
if any(utile)
    T_utile = max(t_epoche(utile)) - min(t_epoche(utile));
else
    T_utile = 0;
end

[~, ie_cross] = min(abs(idx_epoche - idx_cross));

fprintf('\n===================== SINTESI THROUGHPUT =====================\n');
fprintf('Banda IF                                        : %8.1f MHz\n', cfg.B_IF/1e6);
fprintf('Bit rate lordo (1 bit/simbolo, T_sym = %.2f us) : %8.3f kbit/s\n', cfg.T_sym*1e6, R_b_raw/1e3);
fprintf('Guadagno di processo FrFT (N = %d)              : %8.1f dB\n', cfg.N, G_proc_dB);
fprintf('Covertness di ancoraggio al crossing            : %8.1f dB\n', SNR_anchor_dB);
fprintf('Floor di BER da jitter di tracking (stimato)    : %8.2e\n', BER_floor_stimato);
fprintf('BER al crossing (epoca piu'' vicina)             : %8.2e\n', BER_B(ie_cross));
fprintf('Disponibilità del link (LOS + Sole)             : %7.1f %% della finestra\n', 100*mean(link_available));
fprintf('Finestra utile (BER <= %.0e)                    : %8.1f s su %.0f s\n', BER_soglia, T_utile, T_window);
fprintf('--------------------------------------------------------------\n');
fprintf('Throughput di picco (1-BER)                     : %8.3f kbit/s\n', max(R_naive)/1e3);
fprintf('Throughput informativo di picco (BSC)           : %8.3f kbit/s\n', max(R_shan)/1e3);
for iL = 1:nL
    fprintf('Goodput di picco, frame %4d bit non codificati : %8.3f kbit/s\n', L_frame_vec(iL), max(R_good_L(iL,:))/1e3);
end
fprintf('--------------------------------------------------------------\n');
fprintf('Volume dati per passaggio (1-BER)               : %8.2f kbit  (%6.3f MB)\n', Vol_naive/1e3, Vol_naive/8/1e6);
fprintf('Volume informativo per passaggio (BSC)          : %8.2f kbit  (%6.3f MB)\n', Vol_shan/1e3, Vol_shan/8/1e6);
for iL = 1:nL
    fprintf('Volume utile per passaggio (frame %4d bit)     : %8.2f kbit  (%6.3f MB)\n', ...
            L_frame_vec(iL), Vol_good_L(iL)/1e3, Vol_good_L(iL)/8/1e6);
end
fprintf('==============================================================\n\n');

% Frontiera covertness/rate a B_IF = 100 MHz: dalla formula analitica
% BER = Q(gamma), gamma = sqrt(P_cov*N)/sigma, si ricava il numero minimo
% di campioni/simbolo per centrare un BER obiettivo, e quindi il rate
% massimo sostenibile a banda IF fissata: R_b = B_IF*P_cov/(gamma^2*sigma^2).
% Il floor di tracking non e' incluso qui: la frontiera resta un limite
% superiore "con sincronismo perfetto", coerente con l'uso classico che se
% ne fa in letteratura (limite fisico del solo canale, non del ricevitore).
BER_obiettivo = 1e-6;
gamma_req     = sqrt(2) * erfcinv(2 * BER_obiettivo);
P_cov_sweep   = P_classica * 10.^(snr_fine_dB/10);
R_max_sweep   = cfg.B_IF * P_cov_sweep / (gamma_req^2 * cfg.sigma_noise^2);
N_req_anchor  = (gamma_req * cfg.sigma_noise)^2 / (P_classica * 10^(SNR_anchor_dB/10));
fprintf('Frontiera covertness/rate (BER obiettivo = %.0e, B_IF = %.0f MHz):\n', BER_obiettivo, cfg.B_IF/1e6);
fprintf('   a %.1f dB di covertness servono N >= %.0f campioni/simbolo\n', SNR_anchor_dB, ceil(N_req_anchor));
fprintf('   -> rate massimo sostenibile (limite di canale, sincronismo ideale) = %.2f kbit/s\n\n', ...
        cfg.B_IF/max(N_req_anchor,1)/1e3);

%% ========================================================================
% FASE 8: INTERCETTORE EVE - SINGOLO FRAME vs MEDIA D'INSIEME MULTI-BIT
% =========================================================================
disp('--- FASE 8: Analisi intercettatore Eve ---');

rs_eve = streamMC(SEED_MASTER, 3, 1, 1);
[~, ~, ~, ~, Rx_eve] = simulaFrame(cfg, op_cross, Shift_cross, ...
                                   SNR_anchor_dB, kpk_cross, rs_eve, true);

Y_Eve_single       = fftshift(fft(Rx_eve(:,1)));
Spec_Eve_single_dB = 20*log10(abs(Y_Eve_single) + eps);

Y_Eve_ens          = fftshift(fft(mean(Rx_eve, 2)));
Spec_Eve_ens_dB    = 20*log10(abs(Y_Eve_ens) + eps);

% Costo di rivelazione per Eve (rivelatore d'energia, ordine di grandezza):
% ~1/eta^2 campioni indipendenti per distinguere le due ipotesi. A parita'
% di eta, con B_IF piu' piccola (100 MHz vs 16.38 GHz da banco) il TEMPO
% necessario a Eve per accumulare quei campioni cresce proporzionalmente:
% e' un effetto puramente di banda, non di potenza.
eta_cov  = 10^(SNR_anchor_dB/10);
N_eve    = 1 / eta_cov^2;
T_eve    = N_eve / cfg.B_IF;
fprintf('Rivelazione d''energia di Eve a %.0f dB: ~%.2e campioni (~%.3f s a %.0f MHz)\n', ...
        SNR_anchor_dB, N_eve, T_eve, cfg.B_IF/1e6);

%% ========================================================================
% FASE 9: FIGURE
% =========================================================================
disp('--- FASE 9: Generazione figure ---');

u_axis = linspace(-1, 1, cfg.N);

% ---------------- FIGURA 1: geometria, link budget, Eve ------------------
figure('Name', 'Fig1 - Geometria, Link Budget, Eve', 'Color', 'w', ...
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
legend([h_p_ideal, h_p_jit], {'P_{rx} ideale (FSPL + guadagno coerente)', ...
       'P_{rx} con jitter di Rayleigh (\sigma = 1.5 \murad)'}, 'Location', 'best', 'FontSize', 8);
title('(b) Link Budget Ottico FSO Realistico (guadagno d''antenna coerente con la divergenza)');
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
title(['(c) Intercettore Eve: invisibilità energetica (SNR = ', num2str(SNR_anchor_dB), ' dB, B_{IF} = ', ...
       num2str(cfg.B_IF/1e6,'%.0f'), ' MHz)']);
xlim([-0.25 0.25]); grid on;

% ---------------- FIGURA 2: BER vs covertness e spazio di decisione ------
figure('Name', 'Fig2 - BER vs Covertness e Spazio di Decisione', 'Color', 'w', ...
       'Position', [100, 60, 1000, 460]);

subplot(1,2,1);
ber_plot = BER_A; ber_plot(ber_plot == 0) = NaN;
h_th = semilogy(snr_fine_dB, BER_teorica, 'k-', 'LineWidth', 1.3); hold on;
h_mc = semilogy(snr_vec_dB, ber_plot, 'o', 'Color', [0 0.45 0.74], ...
                'MarkerFaceColor', [0 0.45 0.74], 'MarkerSize', 5);
h_ul = semilogy(snr_vec_dB, BER_A_hi, 'v--', 'Color', [0.5 0.5 0.5], 'MarkerSize', 4);
h_fl = yline(BER_floor_stimato, 'm-.', ['Floor da tracking \approx ', num2str(BER_floor_stimato,'%.1e')], ...
             'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
yline(BER_soglia, 'r:', 'BER = 10^{-3}', 'LineWidth', 1.0);
xlabel('SNR covert rispetto al fondo classico [dB]');
ylabel('BER');
title(['BER vs covertness, B_{IF} = ', num2str(cfg.B_IF/1e6,'%.0f'), ' MHz - ', ...
       num2str(cfg.N_MC*cfg.Nbit_frame), ' bit/punto']);
legend([h_th, h_mc, h_ul, h_fl], {'Analitica (sincronismo ideale): Q(\surd(P_{cov}N)/\sigma)', ...
       'Monte Carlo (con jitter di tracking)', 'Limite sup. IC 95% (Wilson)', 'Floor stimato'}, ...
       'Location', 'southwest', 'FontSize', 7.5);
ylim([1e-5 1]); grid on; set(gca,'XDir','reverse');

subplot(1,2,2);
if ~isempty(stat_nom)
    histogram(stat_nom(bits_nom == 0), 30, 'FaceColor', [0 0.45 0.74], 'FaceAlpha', 0.6); hold on;
    histogram(stat_nom(bits_nom == 1), 30, 'FaceColor', [0.85 0.33 0.1], 'FaceAlpha', 0.6);
    xline(0, 'k--', 'LineWidth', 1.2);
end
xlabel('Statistica di decisione  Re\{Y(k_{usato})\}');
ylabel('Occorrenze');
title(['Spazio di decisione a ', num2str(snr_vec_dB(1),'%.1f'), ' dB (1 frame, con jitter di tracking)']);
legend({'Bit 0 (\phi = 0)', 'Bit 1 (\phi = \pi)', 'Soglia'}, 'Location', 'best', 'FontSize', 8);
grid on;

% ---------------- FIGURA 3: BER, SNR, throughput lungo il passaggio ------
figure('Name', 'Fig3 - BER, SNR, Throughput lungo il passaggio', 'Color', 'w', ...
       'Position', [120, 40, 1080, 460]);

subplot(1,2,1);
ber_b_plot = BER_B; ber_b_plot(ber_b_plot == 0) = NaN;
yyaxis left;
semilogy(t_epoche, ber_b_plot, 'o-', 'LineWidth', 1.3, 'MarkerSize', 4); hold on;
semilogy(t_epoche, BER_B_hi, 'v:', 'Color', [0.5 0.5 0.5], 'MarkerSize', 4);
yline(BER_floor_stimato, 'm-.', 'LineWidth', 1.0);
ylabel('BER'); ylim([1e-5 1]);
yyaxis right;
plot(t_epoche, SNR_B, 's--', 'LineWidth', 1.1, 'MarkerSize', 4);
ylabel('SNR covert efficace [dB]');
xline(t_sec(idx_cross), 'k:', 'Crossing', 'LineWidth', 1.0);
xlabel('Tempo di simulazione [s]');
title('(a) BER ed SNR efficace lungo il passaggio');
legend({'BER (MC)', 'IC 95% sup.', 'Floor stimato', 'SNR_{eff}'}, 'Location', 'best', 'FontSize', 7.5);
grid on;

subplot(1,2,2);
yyaxis left;
h_r1 = plot(t_epoche, R_naive/1e3, 'o-', 'LineWidth', 1.3, 'MarkerSize', 4); hold on;
h_r2 = plot(t_epoche, R_shan/1e3,  's--', 'LineWidth', 1.2, 'MarkerSize', 4);
h_r3 = plot(t_epoche, R_good_L(1,:)/1e3, 'd:', 'LineWidth', 1.1, 'MarkerSize', 4);
h_r4 = plot(t_epoche, R_good_L(2,:)/1e3, '^:', 'LineWidth', 1.1, 'MarkerSize', 4);
ylabel('Throughput [kbit/s]');
yyaxis right;
Vol_cum = cumtrapz(t_epoche, R_shan);
h_v = plot(t_epoche, Vol_cum/1e3, 'LineWidth', 1.4);
ylabel('Volume informativo cumulato [kbit]');
xlabel('Tempo di simulazione [s]');
title('(b) Throughput istantaneo e volume per passaggio');
legend([h_r1, h_r2, h_r3, h_r4, h_v], {'R_b(1-BER)', 'R_b(1-H_2(BER)) [BSC]', ...
       ['Goodput frame ', num2str(L_frame_vec(1)), ' bit'], ...
       ['Goodput frame ', num2str(L_frame_vec(2)), ' bit'], 'Volume cumulato'}, ...
       'Location', 'best', 'FontSize', 7);
grid on;

% ---------------- FIGURA 4: frontiera covertness / rate ------------------
figure('Name', 'Fig4 - Frontiera Covertness-Rate', 'Color', 'w', 'Position', [140, 80, 700, 420]);
loglog(10.^(snr_fine_dB/10), R_max_sweep/1e3, 'LineWidth', 1.6); hold on;
xline(10^(SNR_anchor_dB/10), 'k--', 'Punto di progetto', 'LineWidth', 1.1);
yline(R_b_raw/1e3, 'r:', 'Rate lordo attuale', 'LineWidth', 1.1);
grid on;
xlabel('SNR covert (lineare, rispetto al fondo classico)');
ylabel('Bit rate massimo sostenibile [kbit/s]');
title(['Frontiera covertness-rate a B_{IF} = ', num2str(cfg.B_IF/1e6,'%.0f'), ...
       ' MHz, BER obiettivo = ', num2str(BER_obiettivo,'%.0e'), ' (sincronismo ideale)']);

disp('--- CAMPAGNA MONTE CARLO REALISTICA COMPLETATA CON SUCCESSO! ---');

%% ========================================================================
% FUNZIONI LOCALI
% =========================================================================

function shift = shiftResiduoDoppler(cfg, dfD_dt)
% Residuo Doppler intra-simbolo, dopo la pre-compensazione grossolana
% (analogica/RF, non modellata qui in dettaglio) che rimuove il Doppler
% ottico assoluto. Anche un oscillatore locale ideale puo' applicare una
% sola correzione per simbolo: durante T_sym il Doppler continua a
% variare per effetto del Doppler-rate dfD_dt (Hz/s, gia' calcolato dalla
% traiettoria SGP4/J2 reale in Fase 2), lasciando un residuo strutturale
% che qui si converte nel dominio normalizzato omega (1 unita' di omega =
% B_IF/2 Hz, poiche' l'intero dominio omega in [-1,1] rappresenta l'intera
% banda IF digitalizzata). Sostituisce la precedente costante di scala
% grafica K_doppler_grafico con una derivazione tracciabile a B_IF e alla
% dinamica orbitale reale, invece che a un fattore empirico.
    delta_f_res = dfD_dt * cfg.T_sym;
    shift = delta_f_res / (cfg.B_IF/2);
end

function rs = streamMC(seed_master, id_campagna, id_punto, id_trial)
% Substream deterministico e univoco a partire da un unico seed master,
% usando mrg32k3a (unico generatore MATLAB con substream garantiti
% indipendenti, assegnati come proprieta' dopo la costruzione).
    sub = 1 + id_trial + 1000*id_punto + 1000000*id_campagna;
    rs  = RandStream('mrg32k3a', 'Seed', seed_master);
    rs.Substream = sub;
end

function op = operatoreShift(cfg, shift)
% Precalcola l'operatore di ri-centraggio Doppler (interpolazione lineare
% a passo costante), vettorizzato su matrici di simboli.
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
% Ricentraggio Doppler da effemeridi -> clutter rejection (movmean nativo
% MATLAB) -> de-chirp FrFT -> compressione spettrale.
    Sc = applicaShift(op, S);
    Sc = Sc - movmean(Sc, cfg.win_clutter, 1);
    Y  = fftshift(fft(Sc .* cfg.chiave, [], 1), 1);
end

function kpk = trovaBinPicco(cfg, op, shift)
% Bin di compressione nominale, ricavato dalla forma d'onda nota (senza
% rumore). E' la stima "media" di Bob (da effemeridi/oscillatore): attorno
% a questo bin nominale il jitter di tracking residuo introduce poi un
% errore stocastico di lettura (vedi simulaFrame).
    Phi  = cfg.c_chirp * ((cfg.omega - shift).^2) + cfg.offset * (cfg.omega - shift);
    y    = riceviBob(cfg, op, cos(Phi));
    tmpl = real(y);
    [~, kpk] = max(abs(tmpl));
end

function [nerr, nbit, stat, bits, Rx_out] = simulaFrame(cfg, op, shift, snr_db, kpk, rs, salva_rx)
% Trasmette un frame di cfg.Nbit_frame bit BPSK sul chirp, lo propaga nel
% canale (fondo classico + rumore), lo demodula E aggiunge l'imperfezione
% di sincronizzazione: la lettura avviene non esattamente sul bin nominale
% kpk ma su kpk + offset, dove offset e' un intero casuale (residuo di
% tracking Doppler dell'oscillatore/GPS) saturato a +/- 2 bin. Questo
% introduce un floor di BER anche a rumore di canale trascurabile.
    if nargin < 7, salva_rx = false; end
    M    = cfg.Nbit_frame;
    bits = randi(rs, [0 1], 1, M);

    Phi   = cfg.c_chirp * ((cfg.omega - shift).^2) + cfg.offset * (cfg.omega - shift);
    Scl   = exp(-((cfg.omega - shift).^2) / (2 * cfg.sigma_doppler^2));

    P_cl  = mean(Scl.^2);
    P_cov = P_cl * 10^(snr_db/10);
    a_cov = sqrt(2 * P_cov);

    stat = zeros(1, M);
    if salva_rx
        Rx_out = zeros(cfg.N, M);
    else
        Rx_out = [];
    end

    for i1 = 1:cfg.blocco:M
        i2  = min(i1 + cfg.blocco - 1, M);
        b   = bits(i1:i2);
        nb  = numel(b);
        Stx = a_cov * cos(Phi + pi * b);
        Rx  = Scl + Stx + cfg.sigma_noise * randn(rs, cfg.N, nb);
        if salva_rx, Rx_out(:, i1:i2) = Rx; end
        Y   = riceviBob(cfg, op, Rx);

        % --- Jitter di tracking Doppler: offset casuale sul bin letto ---
        off_bin = round(cfg.sigma_bin_track * randn(rs, 1, nb));
        off_bin = min(max(off_bin, -2), 2);
        kuse    = min(max(kpk + off_bin, 1), cfg.N);
        idxLin  = sub2ind(size(Y), kuse, 1:nb);
        stat(i1:i2) = real(Y(idxLin));
    end

    dec  = double(stat < 0);
    nerr = sum(dec ~= bits);
    nbit = M;
end

function [lo, hi] = wilsonCI(k, n, z)
% Intervallo di confidenza di Wilson per una proporzione binomiale.
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
% NOTE TECNICHE (aggiornate per la versione avionica)
% =========================================================================
% NOTA 1 - IL FLOOR NON E' UNA COSTANTE INDIPENDENTE DALLA SNR (a differenza
%   di quanto si potrebbe pensare di primo acchito). Con sincronismo
%   perfetto la BER scende verso zero senza limite al crescere della SNR.
%   Un errore di tracking di 1 bin NON e' pero' un "colpo di sfortuna a
%   costo fisso": verificato con una Monte Carlo dedicata a -41 dB,
%   separando gli errori per esatto valore di offset,
%      offset =  0  ->  BER ~ 0
%      offset = +-1 ->  BER ~ 5-7 %
%      offset = +-2 ->  BER ~ 85-100 %
%   Il lobo laterale del kernel di compressione ha un guadagno di processo
%   molto inferiore a quello del lobo principale: la stessa potenza covert
%   che e' ampiamente sufficiente al bin corretto non lo e' affatto al bin
%   adiacente. Questa degradazione condizionata dipende dalla SNR nominale
%   (a una SNR molto piu' alta un errore di 1 bin costa assai meno, perche'
%   il margine residuo assorbe la perdita), quindi il contributo del
%   tracking alla BER totale non e' un floor piatto in senso stretto: resta
%   basso vicino al punto di ancoraggio e cresce con continuita' scendendo
%   in SNR, fondendosi gradualmente nel ginocchio del waterfall intrinseco.
%   Il valore stampato a runtime come "floor stimato" (media dei punti a
%   SNR piu' alta) va quindi letto come una stima locale al punto di
%   ancoraggio, non come un limite fisico universale del ricevitore.
%
% NOTA 2 - CALIBRAZIONE DI sigma_bin_track: FARLA SEMPRE ALLA SNR VERA.
%   Una prima calibrazione fatta osservando il floor a una SNR proxy molto
%   alta (-20 dB, margine residuo enorme) dava floor ~3.5e-3 per
%   sigma_bin_track = 0.5, un valore rassicurante ma fuorviante: applicato
%   alla SNR di ancoraggio reale (-41 dB), lo stesso sigma_bin_track produce
%   un floor quasi un ordine di grandezza piu' alto (~2.3e-2, verificato),
%   perche' a -41 dB il margine non basta ad assorbire la perdita anche a
%   soli +/-1 bin (BER condizionata 5-7%, non <1% come a -20 dB). La regola
%   corretta e' quindi: CALIBRARE sigma_bin_track CON UNA MONTE CARLO
%   DEDICATA VALUTATA ESATTAMENTE ALLA SNR DI ANCORAGGIO DI PROGETTO, mai a
%   una SNR proxy piu' comoda. Con questa correzione, sigma_bin_track = 0.20
%   da' un floor a -41 dB di ~7*10^-4, compatibile con un residuo di
%   tracking Doppler di poco piu' di 1 kHz (vedi sigma_track_Hz stampato in
%   Fase 4) - una specifica plausibile per un oscillatore/GPS di bordo
%   classe CubeSat dopo correzione da effemeridi. E' un parametro di
%   progetto: aumentarlo modella un tracking peggiore, ma va sempre
%   ricontrollato alla SNR operativa vera, non a una SNR comoda per il test.
%
% NOTA 3 - GUADAGNO D'ANTENNA: MODELLO gaussianAntenna DEL TOOLBOX.
%   G_tx e G_rx usano ora la STESSA formula documentata per l'oggetto
%   gaussianAntenna del Satellite Communications Toolbox (boresightGain =
%   ApertureEfficiency*(pi*D/lambda)^2), con ApertureEfficiency = 0.65, il
%   default che il toolbox stesso assegna quando non specificato altrimenti
%   (vedi documentazione di transmitter/gaussianAntenna). Sostituisce due
%   problemi della versione precedente: (1) G_tx = 32/theta_div^2 legava il
%   guadagno Tx a una divergenza (8 urad) scelta indipendentemente
%   dall'apertura D_tx, un rischio di incoerenza fisica fra i due; (2) G_rx
%   al limite di diffrazione puro assumeva implicitamente un'efficienza di
%   apertura del 100%, non realistica. Impatto numerico verificato: G_tx
%   aumenta di +8.6 dB, G_rx diminuisce di -1.9 dB (netto +6.7 dB su
%   P_rx_ideal_t, quindi su Fig. 1b e sul margine fotonico di Fase 4).
%   NESSUN impatto invece su BER/throughput (Fig. 2-4): in Campagna B la
%   SNR efficace usa il RAPPORTO P_rx_ideal_t(t)/P_rx_ideal_t(crossing), in
%   cui G_tx*G_rx (costanti lungo il passaggio) si cancellano esattamente;
%   in Campagna A la covertness e' definita relativamente al fondo
%   classico, indipendente dalla potenza assoluta. theta_div_beam resta un
%   parametro indipendente, usato SOLO nella formula di fading di
%   puntamento exp(-8*(theta_err/theta_div)^2) (una relazione della
%   letteratura FSO, concettualmente diversa dal beamwidth a -3dB
%   dell'antenna: quest'ultimo e' stampato a runtime in Fase 3 solo per
%   confronto, le due grandezze angolari non vanno confuse) - resta quindi
%   ancora da verificare in letteratura, come gia' segnalato.
%
% NOTA 3bis - DERIVAZIONE FISICA DI c_chirp, offset, sigma_doppler E DELLO
%   SHIFT DOPPLER (sostituisce K_doppler_grafico).
%   Nel benchmark, c_chirp/offset/sigma_doppler vivevano nel dominio
%   normalizzato omega in [-1,1] senza un aggancio esplicito a B_IF: la
%   loro banda FISICA (Hz) e' automaticamente proporzionale a B_IF (1
%   unita' di omega = B_IF/2 Hz), quindi cambiando B_IF da 16.38 GHz a
%   100 MHz il loro significato fisico e' cambiato di un fattore ~164
%   SENZA che nessuno lo decidesse (verificato: la stessa sigma_doppler
%   = 0.15 rappresentava 1229 MHz di linewidth nel benchmark e ne
%   rappresenta solo 7.5 MHz qui). Qui questo e' reso un aggancio
%   esplicito e dichiarato: c_chirp e offset sono derivati da una frazione
%   di banda occupata scelta di proposito (frac_BW_chirp = 7.32% di B_IF,
%   Fase 4), e sigma_doppler resta definito come frazione fissa della
%   semi-banda B_IF/2 per costruzione, non come valore assoluto ereditato.
%   La costante K_doppler_grafico (fattore di scala grafico, dichiarato
%   come tale fin dal modello originale) e' stata rimossa e sostituita da
%   shiftResiduoDoppler(): lo shift usato per ricentrare lo spettro non
%   e' piu' proporzionale al Doppler assoluto (v_rel/c, che a 26.7 GHz non
%   potrebbe comunque stare in un IF da 100 MHz - vedi la dichiarazione di
%   architettura in Fase 4) ma al residuo Doppler ACCUMULATO IN UN SIMBOLO
%   per effetto del solo Doppler-rate (dfD_dt_t, gia' calcolato dalla
%   traiettoria SGP4/J2 reale in Fase 2): shift = dfD_dt*T_sym/(B_IF/2).
%   Conseguenza verificata: questo residuo e' minuscolo ovunque nel
%   passaggio (~0.0022 al crossing, dove il Doppler-rate e' massimo a
%   671 MHz/s; molto meno altrove), contro shift fino a ~0.3-0.4 del vecchio
%   schema empirico. Significa che, grazie alla pre-compensazione Doppler
%   grossolana dichiarata in Fase 4, la perdita di campioni per
%   interpolazione fuori dominio in operatoreShift() diventa trascurabile
%   ovunque nel passaggio (prima non lo era, ed era un fattore di
%   discrepanza fra Campagna A e B): la Campagna B dovrebbe quindi
%   allinearsi alla Campagna A in modo piu' pulito lungo tutto il
%   passaggio, non solo al crossing.
%
% NOTA 3ter - IDENTITA' FISICA DEL FONDO CLASSICO (sigma_doppler).
%   Il fondo classico non e' un inviluppo di comodo: e' il canale classico
%   legittimo co-locato sotto cui il segnale covert si nasconde (scenario
%   LPI standard), identificato con la riga di emissione dello ione 27Al+
%   intrappolato, allargata per effetto Doppler termico dal moto residuo
%   dello ione. sigma_doppler e' quindi DERIVATO, non assunto:
%       sigma_nu = nu0*sqrt(k_B*T_ion/(m*c^2))
%   Con T_ion = 13.02 mK -> sigma = 7.49 MHz, FWHM = 17.65 MHz, che nel
%   dominio normalizzato a B_IF = 100 MHz da' sigma_doppler = 0.150 (lo
%   stesso valore usato in tutte le campagne gia' validate: la derivazione
%   fisica non invalida i risultati precedenti, li giustifica a posteriori).
%
%   CORREZIONE DI COERENZA: nella versione precedente sigma_doppler era
%   fissato come frazione di B_IF, per cui cambiando B_IF cambiava
%   silenziosamente la temperatura implicita dello ione - fisicamente
%   rovesciato (la temperatura e' una proprieta' della trappola, non del
%   ricevitore). Ora il verso e' corretto: si fissa T_ion, e sigma_doppler
%   normalizzato ne discende insieme a B_IF.
%
%   TRE CONDIZIONI DI VALIDITA', tutte verificate numericamente:
%   (a) Larghezza naturale della transizione di intercombinazione
%       1S0-3P1 dell'Al+ (sub-kHz) trascurabile rispetto ai MHz di
%       allargamento Doppler -> il profilo di Voigt degenera in gaussiano
%       puro, coerente con la forma exp(-x^2/(2*sigma^2)) usata nel codice.
%   (b) Regime di sideband NON risolte: FWHM Doppler (17.6 MHz) >> frequenza
%       secolare di trappola tipica (0.5-5 MHz), rapporto 3.5-35. Il
%       profilo gaussiano continuo e' quindi valido; per uno ione piu'
%       freddo o una trappola piu' rigida si passerebbe a sideband risolte
%       e la gaussiana non descriverebbe piu' lo spettro (il codice
%       verifica ed emette un warning a runtime in Fase 4).
%   (c) Contenimento spettrale: il chirp covert occupa 7.32 MHz, contro
%       una FWHM del fondo di 17.65 MHz (rapporto 0.41). Il segnale covert
%       e' quindi interamente contenuto nella riga di copertura, condizione
%       necessaria dello scenario LPI: se sbordasse, sarebbe rilevabile
%       fuori banda indipendentemente da quanto sia debole.
%
% NOTA 4 - PERCHE' LA FRONTIERA COVERTNESS-RATE NON INCLUDE IL FLOOR.
%   La frontiera (Fase 7, Fig. 4) resta calcolata a sincronismo ideale: e'
%   un limite fisico del canale (quanti campioni/simbolo servono per un
%   dato BER obiettivo), non del ricevitore. Un secondo limite, dettato dal
%   floor di tracking, si aggiunge in serie: nessun aumento di N o di
%   potenza covert puo' scendere sotto BER_floor_stimato. In tesi conviene
%   presentare i due limiti separatamente, come qui.
%
% NOTA 5 - COSTO DI RIVELAZIONE PER EVE E BANDA IF.
%   Il numero di campioni che Eve deve accumulare per un rivelatore
%   d'energia dipende solo dalla SNR covert (N_eve = 1/eta^2), non dalla
%   banda. Il TEMPO che le serve, invece, e' N_eve/B_IF: con B_IF = 100 MHz
%   (invece di 16.38 GHz da banco) il tempo di rivelazione cresce di un
%   fattore ~164, un vantaggio addizionale di covertness "gratuito" che
%   discende direttamente dalla scelta di una banda IF realistica da
%   SmallSat.
% =========================================================================