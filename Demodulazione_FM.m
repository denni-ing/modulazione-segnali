fs = 100e3;        % Frequenza di campionamento (Hz)
T = 2;             % Durata del segnale (s)
t = 0:1/fs:T-1/fs; % Vettore dei tempi
fm = 1e3;          % Frequenza della modulante
Am = 1;            % Ampiezza della modulante sinusoidale
fc = 10e3;         % Frequenza della portante (Hz)
Ac = 100;          % Ampiezza della portante
kf = 4e3;          % Costante di deviazione di frequenza
SNR_db = 10:5:40;  % SNR in ingresso (dB)
flagAWGN = 1;      % Flag

% Segnale modulante
m_t = Am*cos(2*pi*fm*t);

% Segnale modulato
theta = 2*pi*(fc*t + kf * cumsum(m_t)/fs); % normalizzato
s_fm = Ac*cos(theta);
signal_power = mean(s_fm.^2);

beta_f = kf*Am/fm; %indice di modulazione

if flagAWGN == 1
    for i = 1:length(SNR_db)
        [noisy_signal, noise_power] = awgn(s_fm, SNR_db(i), 'measured');
        
        % Demodulazione
        % Filtro passa banda
        BW = 2 * (kf*Am + fm); % Larghezza di banda
        fpass = [(fc-BW/2), (fc+BW/2)]/(fs/2); 
        [b_band, a_band] = butter(6, fpass, 'bandpass');
        s_fm_filtered = filtfilt(b_band, a_band, noisy_signal);
        
        % Derivazione numerica
        s_fm_diff = diff(s_fm_filtered);

        % Rivelazione di inviluppo
        s_fm_env = abs(s_fm_diff);
        [b_low, a_low] = butter(3, 1.2e3/(fs/2), 'low'); % Filtro passa-basso
        s_fm_env_filtered = filtfilt(b_low, a_low, s_fm_env);
        
        s_fm_env_filtered = s_fm_env_filtered - mean(s_fm_env_filtered);
        signal_power_out = mean(s_fm_env_filtered.^2)*max(s_fm_env_filtered)^2;
        noise_power_out = noise_power/4;
    
        % SNR_out sperimentale
        SNR_out_exp(i) = 10 * log10(signal_power_out / noise_power_out);
        
        t_len = t(1:length(s_fm_env_filtered));
        figure;
        plot(t_len, s_fm_env_filtered, 'r', 'LineWidth', 1.0);
        xlabel('Tempo (s)');
        ylabel('Ampiezza');
        title(['SNR_{in} = ' num2str(SNR_db(i)) ' dB']);
        legend('Segnale demodulato');
        grid on;
        xlim([0 10/fm]);
        
    end

else
    
end

%confronto con SNR_out teorico
SNR_out_t = 10 * log10(3 * beta_f^2 * 10.^(SNR_db / 10));

figure;
plot(SNR_db, SNR_out_exp, 'ro-', 'LineWidth', 2);
hold on;
plot(SNR_db, SNR_out_t, 'b--o', 'LineWidth', 2);
plot(SNR_db, SNR_db, 'go-', 'LineWidth', 1.5);
xlabel('SNR_{in} (dB)');
ylabel('SNR_{out} (dB)');
title('Confronto SNR_{out} teorico vs sperimentale');
legend('SNR_{out} sperimentale', 'SNR_{out} teorico', 'SNR_{in}');
grid on;