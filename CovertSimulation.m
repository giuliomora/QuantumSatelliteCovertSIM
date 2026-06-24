% =========================================================================
% SIMULAZIONE QUANTUM COVERT COMMUNICATION - FASE 1 e 2
% =========================================================================

% -------------------------------------------------------------------------
% OPERAZIONI DI PULIZIA (Si mettono sempre all'inizio di ogni script)
% -------------------------------------------------------------------------
clc;        % Pulisce la Command Window (lo schermo centrale)
clear;      % Cancella tutte le variabili vecchie dal Workspace
close all;  % Chiude tutti i grafici o finestre rimasti aperti da prima

% =========================================================================
% FASE 1: INIZIALIZZAZIONE (Le Costanti Fisiche e l'Asse di Simulazione)
% =========================================================================

% 1. Costanti Fisiche Generali
c = 3e8;                    % Velocità della luce in m/s
amu = 1.6605e-27;           % Unità di massa atomica in kg
m_Al = 26.98 * amu;         % Massa di uno ione di Alluminio-27 (kg)
lambda_0 = 267e-9;          % Lunghezza d'onda transizione Al+ (267 nm)

% 2. Momenti dei pacchetti d'onda dell'atomo
p1 = 2e-8 * m_Al * c;       % Momento p1
p2 = 4e-8 * m_Al * c;       % Momento p2

% 3. Parametri del Vettore di Simulazione
% Invece di usare un asse in "secondi" o "Hertz" grezzi, in teoria dei segnali 
% per applicare la FrFT si usa spesso un asse "normalizzato" (senza unità di 
% misura fisse) che va da -1 a 1. Rappresenta la nostra finestra di scansione.
Fs = 50e6;                  % Frequenza di campionamento (50 MHz)
T_obs = 1e-3;               % Tempo di osservazione (1 ms)
N_campioni = Fs * T_obs;    % MATLAB calcolerà da solo che sono 50000 punti

% Creiamo l'asse omega. 
% 'linspace' genera N_campioni distanziati equamente tra -1 e 1.
omega = linspace(-1, 1, N_campioni); 

% =========================================================================
% FASE 2: COSTRUZIONE DEL SEGNALE QUANTISTICO (Alice)
% =========================================================================

% 1. Impostiamo il Bit Segreto
bit_da_trasmettere = 0;     % Scegli 0 oppure 1. Puoi cambiarlo per fare test!

if bit_da_trasmettere == 0
    phi = 0;
else
    phi = pi;               % Pi greco in MATLAB si scrive semplicemente 'pi'
end

% 2. La Dinamica di Dispersione (Il Chirp)
% La variazione di fase dipende dalla differenza dei momenti al quadrato.
% Per la simulazione numerica, definiamo un "rate" (velocità) del chirp 
% proporzionale per rendere l'effetto visibile all'algoritmo.
c_chirp = 150 * pi;         % Fattore di variazione quadratica (scelto per test)

% Calcoliamo Phi_chirp. 
% ATTENZIONE: Il punto '.' prima del '^2' è vitale in MATLAB! 
% Significa "eleva al quadrato ogni singolo numero dell'array omega", 
% e non "fai la moltiplicazione tra matrici".
c_chirp = 150 * pi;         % Fattore di variazione quadratica
u_0 = 3000 * pi;            % Offset lineare (Differenza di momento)

% Aggiungiamo l'offset: spingerà il picco fuori dal rumore centrale!
Phi_chirp = c_chirp * (omega.^2) + u_0 * omega;

% 3. L'Inviluppo (Quanto è forte il segnale quantistico?)
% Sappiamo che è debolissimo (ordine v/c). Lo settiamo a un valore minuscolo.
I_quant = 1e-4;             % Intensità bassissima (0.0001)

% 4. Creazione dell'Interferenza (L'equazione fenomenologica)
% S_quant = I_quant * cos(Phi_chirp + phi)
Interferenza_Quantistica = I_quant * cos(Phi_chirp + phi);

% Messaggio di conferma a schermo
disp('Fase 1 e 2 completate! Segnale quantistico generato e pronto in memoria.');


% =========================================================================
% FASE 3: IL CANALE COVERT E IL RUMORE (Eve)
% =========================================================================

% 1. Generazione dello Spettro Classico (Effetto Doppler termico)
% Creiamo una grande curva a campana (Gaussiana) che simula il Doppler.
sigma_doppler = 0.2; % Definisce quanto è larga la campana termica
S_classico = exp(-(omega.^2) / (2 * sigma_doppler^2));

% 2. Scaliamo il segnale quantistico per ottenere un SNR di -41 dB
% Calcoliamo la potenza del rumore classico
Potenza_Classica = mean(S_classico.^2);

% La formula del dB in potenza è: SNR = 10 * log10(P_segnale / P_rumore)
% Quindi calcoliamo la potenza target che deve avere il segnale quantistico:
Potenza_Target_Quant = Potenza_Classica * (10^(-41/10));

% Regoliamo (scaliamo) la nostra interferenza della Fase 2
Potenza_Attuale_Quant = mean(Interferenza_Quantistica.^2);
Fattore_di_Scala = sqrt(Potenza_Target_Quant / Potenza_Attuale_Quant);

% Ecco la nostra increspatura microscopica
Interferenza_Covert = Interferenza_Quantistica * Fattore_di_Scala;

% 3. Aggiungiamo il Rumore di Fondo Spaziale (AWGN)
% Simula i fotoni spuri del sole, del deep space e delle lenti del telescopio.
% 'randn' genera numeri casuali (rumore bianco).
Rumore_Spaziale = 0.05 * randn(1, N_campioni);

% 4. IL SEGNALE RICEVUTO (Quello che transita nello spazio)
% Sommiamo il rumore enorme, l'increspatura invisibile e il rumore spaziale.
S_ricevuto = S_classico + Interferenza_Covert + Rumore_Spaziale;


% =========================================================================
% FASE 4: LA DECODIFICA CON FrFT (Bob)
% =========================================================================

% 0. IL FILTRO PASSA-ALTO (Rimozione del "Clutter" Doppler)
% Calcoliamo il profilo "liscio" della campana di rumore e lo sottraiamo.
% Questo cancella la montagna classica e lascia solo le micro-oscillazioni.
S_filtrato = S_ricevuto - movmean(S_ricevuto, 500);

% 1. La Chiave di Bob (L'angolo ottimale / De-chirp)
% Bob conosce la massa e i momenti, quindi sa generare il chirp opposto
% per compensare esattamente la dispersione spaziale (l'equivalente di alpha_opt).
% Nota l'uso di '1i', che in MATLAB è l'unità immaginaria complessa.
% Bob cancella SOLO la dispersione spaziale (il termine quadratico)
Chiave_FrFT = exp(-1i * (c_chirp * (omega.^2)));

% 2. Applicazione dell'Operatore (Compressione del Chirp)
% Moltiplichiamo il segnale ricevuto per la chiave. Il termine quantistico 
% viene fermato, mentre il rumore viene spalmato.
Segnale_Ruotato = S_filtrato .* Chiave_FrFT;

% 3. Estrazione nel dominio frazionario 'u'
% Applichiamo la Fast Fourier Transform (fft) e la centriamo (fftshift).
Y_Bob = fftshift(fft(Segnale_Ruotato));

% Calcoliamo la potenza estratta per la visualizzazione (modulo quadro)
% Conserviamo il segno reale per capire se il bit è 0 (positivo) o 1 (negativo)
Picco_Bob = real(Y_Bob);

disp('Fase 3 e 4 completate! Canale simulato e segnale decodificato da Bob.');

% =========================================================================
% FASE 5: VISUALIZZAZIONE (I Risultati Finali)
% =========================================================================

% Creiamo l'asse 'u' per il dominio frazionario di Bob
u = linspace(-1, 1, N_campioni);

% Apriamo una nuova finestra "Figure" per i grafici
figure('Name', 'Simulazione Quantum Covert', 'NumberTitle', 'off', 'Color', 'w');

% -------------------------------------------------------------------------
% Grafico 1: Quello che vede Eve (Il Canale Intercettato)
% -------------------------------------------------------------------------
subplot(2, 1, 1); % Divide la finestra in 2 righe, 1 colonna, e usa la 1^ posizione
plot(omega, S_ricevuto, 'Color', [0.6 0.6 0.6]); % Disegna in grigio
title('1. Cosa vede Eve: Spettro Intercettato (SNR = -41 dB)');
xlabel('\omega (Frequenza Normalizzata)');
ylabel('Ampiezza Spettrale');
grid on;

% -------------------------------------------------------------------------
% Grafico 2: Quello che vede Bob (La Decodifica FrFT)
% -------------------------------------------------------------------------
subplot(2, 1, 2); % Usa la 2^ posizione nella finestra
plot(u, Picco_Bob, 'b', 'LineWidth', 1.5); % Disegna in blu con linea più spessa
title(['2. Cosa vede Bob: Decodifica (Bit Trasmesso: ', num2str(bit_da_trasmettere), ')']);
xlabel('u (Dominio Frazionario ruotato)');
ylabel('Ampiezza del Picco Estrazione');
grid on;

% Imposta un limite sull'asse X per fare lo zoom al centro e vedere bene il picco
xlim([-0.2 0.2]);