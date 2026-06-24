% =========================================================================
% SIMULATORE MONTE CARLO - CURVA BER (Bit Error Rate)
% =========================================================================
clc; clear; close all;

disp('--- INIZIO SIMULAZIONE MONTE CARLO ---');
disp('Attendi... elaborazione di migliaia di fotoni spaziali in corso.');
tic; % Avvia il cronometro

% 1. Parametri di Base (Gli stessi già validati)
N_campioni = 50000;
omega = linspace(-1, 1, N_campioni);
c_chirp = 300 * pi;
offset_quantistico = 100 * pi;
Shift_Orbitale = 0.4;
sigma_doppler = 0.15;
I_quant = 1e-7;

% 2. Setup della Simulazione BER
SNR_dB_array = -45:2:-15;       % Variazione del rumore: da -45 dB a -15 dB
N_bits_per_SNR = 200;           % Quanti bit testare per ogni livello di rumore
BER = zeros(1, length(SNR_dB_array)); % Array vuoto per salvare i risultati

% 3. Pre-calcoli per ottimizzare la velocità del Loop
S_classico = exp(-((omega - Shift_Orbitale).^2) / (2 * sigma_doppler^2));
Potenza_Classica = mean(S_classico.^2);
Phi_totale_shiftata = c_chirp * ((omega - Shift_Orbitale).^2) + offset_quantistico * (omega - Shift_Orbitale);
Chiave_FrFT = exp(-1i * (c_chirp * (omega.^2)));

% 4. Il Cuore: Il Ciclo Monte Carlo
for i = 1:length(SNR_dB_array)
    
    snr_corrente = SNR_dB_array(i);
    errori_contati = 0;
    
    % Calcolo della potenza target per questo step di rumore
    Potenza_Target = Potenza_Classica * (10^(snr_corrente/10));
    
    for b = 1:N_bits_per_SNR
        % --- ALICE: Generazione Bit Casuale ---
        bit_tx = randi([0, 1]); % Genera 0 o 1 a caso
        if bit_tx == 0; phi = 0; else; phi = pi; end
        
        Interferenza = I_quant * cos(Phi_totale_shiftata + phi);
        Fattore_di_Scala = sqrt(Potenza_Target / mean(Interferenza.^2));
        Interferenza_Scalata = Interferenza * Fattore_di_Scala;
        
        % --- CANALE SPAZIALE: Aggiunta Rumore ---
        Rumore_Spaziale = 0.5 * randn(1, N_campioni);
        S_ricevuto = S_classico + Interferenza_Scalata + Rumore_Spaziale;
        
        % --- BOB: Ricezione e Decodifica ---
        % Tracking orbitale
        S_centrato = interp1(omega, S_ricevuto, omega + Shift_Orbitale, 'linear', 0);
        
        % Clutter Rejection
        S_filtrato = S_centrato - movmean(S_centrato, 500);
        
        % FrFT ed Estrazione
        Y_Bob = real(fftshift(fft(S_filtrato .* Chiave_FrFT)));
        
        % Logica di Decisione Automatica: Bob cerca il picco massimo assoluto.
        % Se il picco è positivo, deduce 0. Se negativo, deduce 1.
        [~, max_idx] = max(abs(Y_Bob)); 
        if Y_Bob(max_idx) > 0
            bit_rx = 0;
        else
            bit_rx = 1;
        end
        
        % Contatore Errori
        if bit_rx ~= bit_tx
            errori_contati = errori_contati + 1;
        end
    end
    
    % Calcola la percentuale di errore per questo livello di SNR
    BER(i) = errori_contati / N_bits_per_SNR;
    
    % Stampa a schermo l'avanzamento
    disp(['SNR = ', num2str(snr_corrente), ' dB | Errori: ', num2str(errori_contati), '/', num2str(N_bits_per_SNR), ' | BER = ', num2str(BER(i))]);
end

tempo_totale = toc; % Ferma il cronometro
disp(['Simulazione conclusa in ', num2str(tempo_totale), ' secondi.']);

% =========================================================================
% FASE 5: GRAFICO DELLA CURVA A CASCATA (Waterfall)
% =========================================================================

figure('Name', 'Prestazioni Ricevitore', 'Color', 'w');
% In telecomunicazioni, il BER si disegna SEMPRE con asse Y logaritmico (semilogy)
semilogy(SNR_dB_array, BER, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
title('Analisi delle Prestazioni: Bit Error Rate vs SNR Satellitare');
xlabel('Signal-to-Noise Ratio (SNR) [dB]');
ylabel('Bit Error Rate (BER)');

% Limiti estetici per rendere il grafico perfetto per la tesi
ylim([10^-3 1]); 
xlim([-45 -15]);