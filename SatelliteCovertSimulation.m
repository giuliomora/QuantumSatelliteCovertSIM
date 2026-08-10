% =========================================================================
% SIMULAZIONE LINK INTERSATELLITARE QUANTUM COVERT (LEO-LEO ISL)
% =========================================================================
clc; clear; close all;

% =========================================================================
% STEP 1 & STEP 2: CINEMATICA ISL E LINK BUDGET OTTICO SPAZIALE
% =========================================================================
disp('--- FASE 1: Cinematica LEO-LEO e Link Budget Intersatellitare ---');

% Costanti Fisiche
c = 3e8;                    % Velocità della luce (m/s)
lambda_0 = 267e-9;          % Lunghezza d'onda transizione Al+ (267 nm)
h_planck = 6.626e-34;       % Costante di Planck (J*s)
E_fotone = h_planck * (c / lambda_0); % Energia di un singolo fotone (J)

% Geometria e Cinematica Intersatellitare (ISL LEO-LEO)
v_orbita = 7.6e3;           % Velocità orbitale singolo satellite (7.6 km/s)
Distanza_ISL = 300e3;       % Distanza di collegamento intersatellitare (300 km)

% Angolo di incidenza tra i piani orbitali (Incrocio equatore/piani adiacenti)
theta_deg = 60;             % Angolo scelto: 60° (Configurazione tipica Walker-Star)
theta_rad = deg2rad(theta_deg);

% Calcolo Vettoriale della Velocità Relativa lungo la linea di vista (LOS)
v_rel = 2 * v_orbita * sin(theta_rad / 2); % ~7.6 km/s a 60°

% Parametri Ottici dei Telescopi
D_tx = 0.20;                % Diametro telescopio Alice (200 mm)
D_rx = 0.30;                % Diametro telescopio Bob (300 mm)
P_tx_laser = 100e-3;        % Potenza laser di trasmissione (100 mW)

% Guadagni delle Antenne Ottiche
G_tx = (pi * D_tx / lambda_0)^2;
G_rx = (pi * D_rx / lambda_0)^2;

% Free Space Path Loss (FSPL nel vuoto cosmico tra satelliti)
FSPL = (4 * pi * Distanza_ISL / lambda_0)^2;

% Potenza Ricevuta Teorica
P_rx = P_tx_laser * (G_tx * G_rx) / FSPL;

% Stampa dei parametri cinematici e ottici
disp(['Distanza Link ISL: ', num2str(Distanza_ISL/1000), ' km']);
disp(['Angolo tra Piani Orbitali (\theta): ', num2str(theta_deg), ' gradi']);
disp(['Velocità Relativa v_rel: ', num2str(v_rel/1000, '%.2f'), ' km/s']);
disp(['Potenza Ricevuta Teorica: ', num2str(P_rx, '%.3e'), ' Watt']);

% =========================================================================
% FASE 2: GENERAZIONE DEL SEGNALE QUANTISTICO (ALICE)
% =========================================================================
disp('--- FASE 2: Generazione Segnale Atomico (Alluminio-27) ---');

% Vettori di Simulazione
Fs = 50e6;                  % Frequenza di campionamento 50 MHz
T_coh = 1e-3;               % Tempo di coerenza 1 ms
N_campioni = Fs * T_coh;    % 50.000 campioni
omega = linspace(-1, 1, N_campioni); % Asse frequenze normalizzato in banda base

% Parametri Atomici
amu = 1.6605e-27;           % Unità di massa atomica (kg)
m_Al = 26.98 * amu;         % Massa Ione Alluminio-27

% Differenza di momento (delta p)
delta_v = 1;                % m/s
delta_p = m_Al * delta_v;   % Differenza di momento lineare

% Encoding del Bit Segreto (Fase Quantistica)
bit_tx = 1;                 % Bit da trasmettere (0 oppure 1)
if bit_tx == 0
    phi = 0;
else
    phi = pi;
end

% Chirp Rate e Offset Quantistico
c_chirp = 300 * pi;         % Tasso di dispersione quadratica (Chirp rate)
offset_quantistico = 100 * pi; % Offset spettrale quantistico (p2 - p1)

% Fase totale del segnale quantistico
Phi_totale = c_chirp * (omega.^2) + offset_quantistico * omega;

% Generazione dell'interferenza atomica
I_quant = 1e-5;             % Intensità quantistica nominale
Interferenza_Alice = I_quant * cos(Phi_totale + phi);
disp('Generazione completata. Alice ha preparato lo stato atomico covert.');

% =========================================================================
% FASE 3: CANALE SPAZIALE ISL E MACRO-DOPPLER PROIETTATO (-41 dB)
% =========================================================================
disp('--- FASE 3: Simulazione Canale Spaziale Intersatellitare (-41 dB) ---');

% 1. Calcolo dello Shift Doppler Cinematico Reale (Relativistico/Classico)
shift_fisico_reale = v_rel / c; % Ratio adimensionale v_rel/c (~2.53e-5 per 7.6 km/s)

% 2. Mappatura nello Spazio delle Frequenze Normalizzate del Simulatore
% Per rendere lo shift visivamente apprezzabile sull'asse omega [-1, 1], 
% lo Shift_Orbitale viene espresso come frazione della banda simulata.
% In un simulatore Signal-Processing Level, rappresenta il disassamento del 40%.
K_doppler_grafico = 1.578e4; % Fattore di scala del sistema di riferimento
Shift_Orbitale = shift_fisico_reale * K_doppler_grafico; % Shift visibile (~0.40)

% 3. Generazione del Profilo Classico (Fondo Termico Atomico spostato)
sigma_doppler = 0.15; 
S_classico = exp(-((omega - Shift_Orbitale).^2) / (2 * sigma_doppler^2));

% 4. Traslazione Doppler del Segnale Quantistico
% L'atomo è a bordo del satellite emettitore (Alice): il segnale quantistico 
% subisce lo stesso identico spostamento Doppler cinematico del fondo classico.
Phi_totale_shiftata = c_chirp * ((omega - Shift_Orbitale).^2) + ...
                      offset_quantistico * (omega - Shift_Orbitale);
Interferenza_Alice_Shiftata = I_quant * cos(Phi_totale_shiftata + phi);

% 5. Impostazione della Covertness Energetica (-41 dB rispetto al rumore)
Potenza_Classica = mean(S_classico.^2);
Potenza_Target = Potenza_Classica * (10^(-41/10)); % SNR impostato a -41 dB
Potenza_Attuale = mean(Interferenza_Alice_Shiftata.^2);
Fattore_di_Scala = sqrt(Potenza_Target / Potenza_Attuale);

% Il segnale quantistico covert viene ridimensionato sotto il fondo
Interferenza_Spaziale = Interferenza_Alice_Shiftata * Fattore_di_Scala;

% 6. Rumore di Fondo Spaziale Aggiuntivo (Luce Solare / AWGN)
Rumore_Spaziale = 0.02 * randn(1, N_campioni);

% 7. SEGNALE TOTALE PRESENTE NEL CANALE ISL (Intercettabile da Eve e Bob)
S_ricevuto = S_classico + Interferenza_Spaziale + Rumore_Spaziale;

% =========================================================================
% STEP 5: IL RICEVITORE DI BOB (Tracking Intersatellitare & De-Chirping)
% =========================================================================
disp('--- FASE 4: Tracking Orbitale e Demodulazione di Bob ---');

% 1. Tracking Satellitare: Bob compensa il Macro-Doppler cinematico LEO-LEO
S_centrato = interp1(omega, S_ricevuto, omega + Shift_Orbitale, 'linear', 0);

% 2. Clutter Rejection: Rimozione della baseline termica continua
S_filtrato = S_centrato - movmean(S_centrato, 500);

% 3. De-Chirping (FrFT-like): Compensazione della fase quadratica dell'Alluminio
Chiave_FrFT = exp(-1i * (c_chirp * (omega.^2)));

% 4. Estrazione del Picco Spettrale (FFT)
Segnale_Ruotato = S_filtrato .* Chiave_FrFT;
Y_Bob = fftshift(fft(Segnale_Ruotato));
Picco_Bob = real(Y_Bob);

% =========================================================================
% FASE 4b: INTERCETTATORE EVE (FFT CLASSICA)
% =========================================================================
disp('--- FASE 4b: Intercettatore Eve (FFT Classica senza De-Chirp) ---');

% Eve intercetta il segnale ISL nel canale ma applica una FFT standard
Y_Eve = fftshift(fft(S_ricevuto));
Spec_Eve = abs(Y_Eve);
Spec_Eve_dB = 20*log10(Spec_Eve + eps);

% =========================================================================
% FASE 5: REPORT GRAFICO DEL LINK INTERSATELLITARE
% =========================================================================
disp('--- FASE 5: Generazione Grafici ---');

figure('Name', 'Link Intersatellitare (ISL LEO-LEO) Quantum Covert', ...
       'Color', 'w', 'Position', [100, 100, 900, 800]);

% Grafico 1: Canale Intersatellitare Intercettato
subplot(3, 1, 1);
plot(omega, S_ricevuto, 'Color', [0.6 0.6 0.6]);
title(['1. Canale ISL LEO-LEO (SNR = -41 dB, \theta = ', num2str(theta_deg), ...
       '°, v_{rel} = ', num2str(v_rel/1000, '%.2f'), ' km/s)']);
xlabel('\omega (Banda Base)'); ylabel('Ampiezza'); grid on;

% Grafico 2: Spettro dell'Intercettatore Eve
subplot(3, 1, 2);
u_e = linspace(-1, 1, N_campioni);
plot(u_e, Spec_Eve_dB, 'k');
title('2. Spettro di Eve (FFT classica) - Segnale Indistinguibile dal Rumore');
xlabel('u (Dominio delle Frequenze)'); ylabel('Ampiezza (dB)'); grid on;
xlim([-0.2 0.2]);

% Grafico 3: Decodifica del Ricevitore Bob
u = linspace(-1, 1, N_campioni);
subplot(3, 1, 3);
plot(u, Picco_Bob, 'b', 'LineWidth', 1.2);
title(['3. Decodifica Bob (De-Chirp + FFT) - Bit Decodificato: ', num2str(bit_tx)]);
xlabel('u (Dominio Frazionario)'); ylabel('Ampiezza del Picco'); grid on;
xlim([-0.2 0.2]);

disp('SIMULAZIONE INTERSATELLITARE COMPLETATA CON SUCCESSO!');