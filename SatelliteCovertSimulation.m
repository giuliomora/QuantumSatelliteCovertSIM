% =========================================================================
% SIMULAZIONE AVANZATA LINK INTERSATELLITARE QUANTUM COVERT (LEO-LEO ISL)
% Modello: Ione 27Al+ | Quantum Time Dilation & Chirp De-dispersion (FrFT)
% =========================================================================
clc; clear; close all;
%% ========================================================================
% FASE 1: CINEMATICA ISL E LINK BUDGET OTTICO SPAZIALE
% =========================================================================
disp('--- FASE 1: Cinematica LEO-LEO e Link Budget Intersatellitare ---');
% Costanti Fisiche Fondamentali
c = 3e8;                    % Velocità della luce (m/s)
lambda_0 = 267e-9;          % Transizione UV ione 27Al+ (267 nm)
h_planck = 6.626e-34;       % Costante di Planck (J*s)
E_fotone = h_planck * (c / lambda_0); % Energia singolo fotone (~7.44e-19 J)
% Parametri Orbitali Intersatellitari (Configurazione Walker-Star)
v_orbita = 7.6e3;           % Velocità orbitale singolo satellite (7.6 km/s)
Distanza_ISL = 300e3;       % Link distance ISL (300 km)
theta_deg = 60;             % Angolo di incrocio tra piani orbitali
theta_rad = deg2rad(theta_deg);
% Velocità Relativa LOS
v_rel = 2 * v_orbita * sin(theta_rad / 2); % ~7.60 km/s a 60°
% Parametri Ottici dei Terminali ISL
D_tx = 0.20;                % Telescopio Alice (200 mm)
D_rx = 0.30;                % Telescopio Bob (300 mm)
P_tx_laser = 100e-3;        % Potenza di pompaggio laser (100 mW)
% Guadagni Ottici e Attenuazione Geometrica nello Spazio Libero
G_tx = (pi * D_tx / lambda_0)^2;
G_rx = (pi * D_rx / lambda_0)^2;
FSPL = (4 * pi * Distanza_ISL / lambda_0)^2;
P_rx = P_tx_laser * (G_tx * G_rx) / FSPL;
disp(['Link ISL Distance: ', num2str(Distanza_ISL/1000), ' km']);
disp(['Velocità Relativa (v_rel): ', num2str(v_rel/1000, '%.2f'), ' km/s']);
disp(['Potenza Ottica Ricevuta Teorica: ', num2str(P_rx, '%.3e'), ' W']);
%% ========================================================================
% FASE 2: FISICA ATOMICA E GENERAZIONE SEGNALE QUANTISTICO (ALICE)
% =========================================================================
disp('--- FASE 2: Generazione Segnale Atomico e Modulazione di Fase ---');
% Griglia di Simulazione in Banda Base Normalizzata
N_campioni = 100000;         % Risoluzione ad alta densità
omega = linspace(-1, 1, N_campioni);
% Parametri Chirp Quantistico (Grochowski et al. 2021)
c_chirp = 300 * pi;          % Chirp rate (rad/s^2 normalizzato)
offset_quantistico = 100 * pi; % Termine di disassamento atomico
% Fase Quantistica Totale (Dispersione Quadratica)
Phi_totale = c_chirp * (omega.^2) + offset_quantistico * omega;
% Generazione di entrambi i simboli (Bit 0 -> phi=0, Bit 1 -> phi=pi)
I_quant = 1e-5;
S_atomico_bit0 = I_quant * cos(Phi_totale + 0);    % Simbolo Bit 0
S_atomico_bit1 = I_quant * cos(Phi_totale + pi);   % Simbolo Bit 1
%% ========================================================================
% FASE 3: TRASMISSIONE NEL CANALE ISL, DOPPLER E RUMORE (COVERAGE -41 dB)
% =========================================================================
disp('--- FASE 3: Iniezione nel Canale Spaziale ed Effetto Doppler ---');
% Shift Doppler Cinematico Relativo
shift_fisico_reale = v_rel / c; 
K_doppler_grafico = 1.578e4; 
Shift_Orbitale = shift_fisico_reale * K_doppler_grafico; % ~0.40
% Profilo Termico/Gaussiano di Fondo (Doppler Broadening classico)
sigma_doppler = 0.15; 
S_classico = exp(-((omega - Shift_Orbitale).^2) / (2 * sigma_doppler^2));
% Traslazione del segnale quantistico (trasportato dalla sorgente mobile)
Phi_shiftata = c_chirp * ((omega - Shift_Orbitale).^2) + ...
               offset_quantistico * (omega - Shift_Orbitale);
S_quant_tx0 = I_quant * cos(Phi_shiftata + 0);
S_quant_tx1 = I_quant * cos(Phi_shiftata + pi);
% Normalizzazione Potenza per Covertness Energetica (SNR = -41 dB)
SNR_target_dB = -41;
Potenza_Classica = mean(S_classico.^2);
Potenza_Target = Potenza_Classica * 10^(SNR_target_dB / 10);
Fattore_Scala = sqrt(Potenza_Target / mean(S_quant_tx0.^2));
S_quant_tx0 = S_quant_tx0 * Fattore_Scala;
S_quant_tx1 = S_quant_tx1 * Fattore_Scala;
% Rumore di canale (Fondo cosmico / Solare / Shot Noise)
rng(42); % Riproducibilità del rumore AWGN
Rumore_Spaziale = 0.02 * randn(1, N_campioni);
% Segnali totali ricevuti nel canale per Bit 0 e Bit 1
S_rx_bit0 = S_classico + S_quant_tx0 + Rumore_Spaziale;
S_rx_bit1 = S_classico + S_quant_tx1 + Rumore_Spaziale;
%% ========================================================================
% FASE 4: RICEVITORE DI BOB (Tracking, Clutter Rejection, FrFT De-Chirp)
% =========================================================================
disp('--- FASE 4: Demodulazione Coerente di Bob ---');
% 1. Doppler Tracking (Spostamento del centro banda tramite effemeridi)
S_centrato_bit0 = interp1(omega, S_rx_bit0, omega + Shift_Orbitale, 'linear', 0);
S_centrato_bit1 = interp1(omega, S_rx_bit1, omega + Shift_Orbitale, 'linear', 0);
% 2. Clutter Rejection (Rimozione inviluppo lento classico)
S_clean_bit0 = S_centrato_bit0 - movmean(S_centrato_bit0, 1000);
S_clean_bit1 = S_centrato_bit1 - movmean(S_centrato_bit1, 1000);
% 3. De-Chirping Operatore FrFT: exp(-i * c_chirp * omega^2)
Chiave_FrFT = exp(-1i * (c_chirp * (omega.^2)));
% 4. Compressione Spettrale (Fourier Frazionaria)
Y_Bob_bit0 = fftshift(fft(S_clean_bit0 .* Chiave_FrFT));
Y_Bob_bit1 = fftshift(fft(S_clean_bit1 .* Chiave_FrFT));
Picco_Bob_bit0 = real(Y_Bob_bit0);
Picco_Bob_bit1 = real(Y_Bob_bit1);
% Metriche di decisione (Processing Gain)
[val_max_0, idx_max_0] = max(Picco_Bob_bit0);
[val_min_1, idx_min_1] = min(Picco_Bob_bit1);
rumore_rms = std(Picco_Bob_bit0([1:10000, end-10000:end]));
SNR_post_Bob_dB = 20*log10(abs(val_max_0) / rumore_rms);
disp(['Processing Gain: Segnale recuperato a SNR = +', num2str(SNR_post_Bob_dB, '%.2f'), ' dB']);
%% ========================================================================
% FASE 5: INTERCETTORE EVE (FFT Classica e Analisi Statistica di Ensemble)
% =========================================================================
disp('--- FASE 5: Analisi Intercettatore Eve ---');
% 1. Eve su singolo frame (spettro FFT standard non de-chirpato)
Y_Eve_single = fftshift(fft(S_rx_bit0));
Spec_Eve_single_dB = 20*log10(abs(Y_Eve_single) + eps);
% 2. Eve media statistica su molti simboli (Covertness Quantistica)
% Con bit equiprobabili (p=0.5), i termini quantistici opposti si elidono esattamente:
S_Eve_ensemble = 0.5 * S_rx_bit0 + 0.5 * S_rx_bit1;
Y_Eve_ensemble = fftshift(fft(S_Eve_ensemble));
Spec_Eve_ensemble_dB = 20*log10(abs(Y_Eve_ensemble) + eps);
%% ========================================================================
% FASE 6: REPORT GRAFICO AD ALTA DEFINIZIONE PER LA TESI
% =========================================================================
disp('--- FASE 6: Generazione Grafici ---');
figure('Name', 'Analisi Completa Quantum Covert ISL', 'Color', 'w', 'Position', [80, 50, 1000, 850]);
u_axis = linspace(-1, 1, N_campioni);
% Sottografico 1: Canale Spaziale Grezzo
subplot(4, 1, 1);
plot(omega, S_rx_bit0, 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8);
title(['(a) Segnale nel Canale ISL (SNR_{in} = -41 dB, v_{rel} = ', num2str(v_rel/1000, '%.2f'), ' km/s, \theta = 60°)']);
xlabel('\omega (Banda Base Normalizzata)'); ylabel('Ampiezza');
legend('Segnale Ricevuto S_{rx}(\omega)', 'Location', 'northeast');
grid on;
% Sottografico 2: Spettro Intercettatore Eve (Singolo Shot)
subplot(4, 1, 2);
plot(u_axis, Spec_Eve_single_dB, 'Color', [0.2 0.2 0.2], 'LineWidth', 0.9);
title('(b) Intercettatore Eve: FFT Standard su Singolo Pacchetto (Nessuna Traccia di Chirp Quantistico)');
xlabel('Frequenza Spaziale u'); ylabel('Densità Spettrale (dB)');
xlim([-0.25 0.25]); grid on;
% Sottografico 3: Spettro Intercettatore Eve (Media di Ensemble)
subplot(4, 1, 3);
plot(u_axis, Spec_Eve_ensemble_dB, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.0);
title('(c) Covertness Quantistica: Spettro di Eve mediato su Bit Equiprobabili \rightarrow \langle\Psi\rangle = 0');
xlabel('Frequenza Spaziale u'); ylabel('Densità Spettrale (dB)');
xlim([-0.25 0.25]); grid on;
% Sottografico 4: Decodifica di Bob (Confronto Bit 0 vs Bit 1)
subplot(4, 1, 4);
plot(u_axis, Picco_Bob_bit0, 'b', 'LineWidth', 1.3);
hold on;
plot(u_axis, Picco_Bob_bit1, 'm', 'LineWidth', 1.3);
yline(0, 'k:');
title(['(d) Ricevitore Bob (De-Chirp FrFT + FFT): Picco Coerente (SNR_{out} \approx +', ...
       num2str(SNR_post_Bob_dB, '%.1f'), ' dB)']);
xlabel('u (Dominio Frazionario FrFT)'); ylabel('Ampiezza Compressa');
legend('Bit Decodificato "0" (cos(0) = +1)', 'Bit Decodificato "1" (cos(\pi) = -1)', 'Location', 'northeast');
xlim([-0.1 0.1]); grid on;
disp('--- SIMULAZIONE AVANZATA COMPLETATA CON SUCCESSO! ---');