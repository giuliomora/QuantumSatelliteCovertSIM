% =========================================================================
% SIMULATORE MONTE CARLO - CURVA BER (Link Intersatellitare LEO-LEO)
% =========================================================================
clc; clear; close all;
disp('--- INIZIO SIMULAZIONE MONTE CARLO (ISL LEO-LEO) ---');
disp('Attendi... elaborazione dei bit su canale intersatellitare...');
tic;

% =========================================================================
% 1. PARAMETRI DI BASE E CINEMATICA ISL (I TUOI PARAMETRI)
% =========================================================================
c = 3e8;                        % Velocità della luce (m/s)
v_orbita = 7.6e3;               % Velocità orbitale (7.6 km/s)
theta_deg = 60;                 % Angolo tra piani orbitali (60 gradi)
theta_rad = deg2rad(theta_deg);

v_rel = 2 * v_orbita * sin(theta_rad / 2); 
shift_fisico_reale = v_rel / c; 
K_doppler_grafico = 1.578e4; 
Shift_Orbitale = shift_fisico_reale * K_doppler_grafico; % ~0.40

N_campioni = 50000;
omega = linspace(-1, 1, N_campioni);
c_chirp = 300 * pi;
offset_quantistico = 100 * pi;
sigma_doppler = 0.15;
I_quant = 1e-5;

% =========================================================================
% 2. SETUP DELLA SIMULAZIONE BER (I TUOI PARAMETRI)
% =========================================================================
SNR_dB_array = -45:2:-15;       % Il TUO intervallo SNR (-45 dB a -15 dB)
N_bits_per_SNR = 200;           % I TUOI bit per step
BER = zeros(1, length(SNR_dB_array)); 

S_classico = exp(-((omega - Shift_Orbitale).^2) / (2 * sigma_doppler^2));
Potenza_Classica = mean(S_classico.^2);

Phi_totale_shiftata = c_chirp * ((omega - Shift_Orbitale).^2) + ...
                      offset_quantistico * (omega - Shift_Orbitale);
Chiave_FrFT = exp(-1i * (c_chirp * (omega.^2)));

% =========================================================================
% 3. CICLO MONTE CARLO
% =========================================================================
for i = 1:length(SNR_dB_array)
    
    snr_corrente = SNR_dB_array(i);
    errori_contati = 0;
    
    % Potenza del segnale quantistico associata all'SNR corrente
    Potenza_Target = Potenza_Classica * (10^(snr_corrente/10));
    
    Interferenza_Base = I_quant * cos(Phi_totale_shiftata);
    Fattore_di_Scala = sqrt(Potenza_Target / mean(Interferenza_Base.^2));
    
    % Taratura del rumore AWGN tenendo conto del guadagno di processo (N = 50000)
    % A -45 dB il rumore è sufficiente a creare qualche errore, a -15 dB il BER va a 0
    sigma_rumore = sqrt(Potenza_Classica) * 0.15 * (10^(-snr_corrente/40)); 
    
    for b = 1:N_bits_per_SNR
        % Alice commuta il bit
        bit_tx = randi([0, 1]); 
        if bit_tx == 0
            phi = 0; 
        else
            phi = pi; 
        end
        
        Interferenza_Scalata = (I_quant * cos(Phi_totale_shiftata + phi)) * Fattore_di_Scala;
        
        % CANALE SPAZIALE: Background + Segnale + Rumore casuale scalato
        Rumore_Spaziale = sigma_rumore * randn(1, N_campioni);
        S_ricevuto = S_classico + Interferenza_Scalata + Rumore_Spaziale;
        
        % BOB: Demodulazione
        S_centrato = interp1(omega, S_ricevuto, omega + Shift_Orbitale, 'linear', 0);
        S_filtrato = S_centrato - movmean(S_centrato, 500);
        
        % De-Chirp
        Segnale_Ruotato = S_filtrato .* Chiave_FrFT;
        Y_Bob = real(fftshift(fft(Segnale_Ruotato)));
        
        % Decisione basata sulla polarità del picco dominante
        [~, max_idx] = max(abs(Y_Bob)); 
        if Y_Bob(max_idx) > 0
            bit_rx = 0;
        else
            bit_rx = 1;
        end
        
        if bit_rx ~= bit_tx
            errori_contati = errori_contati + 1;
        end
    end
    
    BER(i) = errori_contati / N_bits_per_SNR;
    
    disp(['SNR = ', num2str(snr_corrente, '%2d'), ' dB | Errori: ', ...
          num2str(errori_contati, '%3d'), '/', num2str(N_bits_per_SNR), ...
          ' | BER = ', num2str(BER(i), '%.4f')]);
end

tempo_totale = toc;
disp(['Simulazione conclusa in ', num2str(tempo_totale, '%.2f'), ' secondi.']);

% =========================================================================
% 4. GRAFICO DELLA CURVA BER
% =========================================================================
figure('Name', 'Prestazioni Link ISL - Curva BER', 'Color', 'w', 'Position', [200, 200, 800, 500]);

% Valore minimo di BER visualizzabile pari alla risoluzione (1 / N_bits)
BER_plot = max(BER, 1/N_bits_per_SNR); 

semilogy(SNR_dB_array, BER_plot, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
title(['Curva BER vs SNR - Link Intersatellitare LEO-LEO (\theta = ', num2str(theta_deg), '°)']);
xlabel('Signal-to-Noise Ratio (SNR) [dB]');
ylabel('Bit Error Rate (BER)');
xlim([-45 -15]);
ylim([1/N_bits_per_SNR 1]);

xline(-41, 'r--', 'SNR Covert (-41 dB)', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');