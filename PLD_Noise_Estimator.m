% Power Level Different Noise Estimator

clear

% Select Test Data
speech = "01";
noise_type = "babble";
snr_level = ["0" "5" "10" "15"];

%Add Audio Path
addpath('/Users/napatr.tan/Documents/MATLAB/SeniorProject/Input')

% Simulation Parameter
frameLength = 256;                              % frame length
npoint = 256;                                   % n-point of fft
sp = 0.5;                                       % percent of overlap
Length_Frame = sp*frameLength;                  % overlap frame
Length_Freq_Index  = sp*npoint + 1;             % no overlap frame
segment_window = rectwin(frameLength)';         % Rectangle Window
mic_len = 0.1;                                  % Distance between 2 Mic
sound_speed = 340;                              % Speed of Sound

% Smoothing Factor
a2 = 0.9;                   % Forgetting Factor Phi_Noise
a3 = 0.8;                   % Forgetting Factor Phi_Noise
q = 0.75;                   % HR Coeff
  
% Import Clean Speech
im_file = "sp"+speech+".wav";
[signal_c, Fs_c] = audioread(im_file);

% Create Discrete Frequency Vector
fsignal = Fs_c*(0:(npoint/2))/npoint;

% Every SNR Level
for j = 1:size(snr_level,2)

    % Import Noisy Signal
    noisy_file = "sp"+speech+"_"+noise_type+"_sn"+snr_level(j)+".wav";
    [signal, Fs] = audioread(noisy_file);
    
    % Create Dual Channel Signal
    [speech1,speech2,noise1,noise2] = dual_mic(signal_c,signal);
    signal1 = speech1 + noise1;
    signal2 = speech2 + noise2;

    % Select Gramma Value for Wiener Filter
    snr_12 = 20*log10(norm(signal1)/norm(signal2));
    if snr_12 < 3
        gramma = 3.5;
    elseif snr_12 < 6
        gramma = 3.8;
    elseif snr_12 < 8
        gramma = 1.5;
    else
        gramma = 0.2;
    end

    % Calculate Signal Parameter
    nframe = floor(length(signal1)/(frameLength*(1-sp)))-1;

    % Create Output Data
    en_signal = zeros(1,(nframe+1)*Length_Frame);
    output_HR = zeros(1,(nframe+1)*Length_Frame);

    % Create Processing Data
    phi_noise = zeros(nframe,Length_Freq_Index);
    G = zeros(nframe,Length_Freq_Index);

    %% Signal Processing
    % Initial Frame Window
    currentFrame = 1;
    nextFrame = currentFrame + frameLength-1;
    
    % Initial Err Value
    Log_Err_C = 0;

    % For Every Frame
    for k = 1:nframe

        % Import Current Frame Signal
        signal1_d = signal1(currentFrame:nextFrame)';
        signal2_d = signal2(currentFrame:nextFrame)';
        clean_signal = speech1(currentFrame:nextFrame)';

        % Fourier Transform
        signal1_ft_double = fft(signal1_d,npoint);
        signal1_ft = signal1_ft_double(1:npoint/2+1);
        signal1_ft(2:end-1) = 2*signal1_ft(2:end-1);
        signal1_ft_m = abs(signal1_ft);             % Magnitude Data
        signal1_ft_ph = angle(signal1_ft);          % Phase Data

        signal2_ft_double = fft(signal2_d,npoint);
        signal2_ft = signal2_ft_double(1:npoint/2+1);
        signal2_ft(2:end-1) = 2*signal2_ft(2:end-1);
        signal2_ft_m = abs(signal2_ft);             % Magnitude Data
        signal2_ft_ph = angle(signal2_ft);          % Phase Data

        % PSD Calculation
        psdsignal1 = periodogram(signal1_d,segment_window,npoint,Fs,'power','reassigned')';
        psdsignal2 = periodogram(signal2_d,segment_window,npoint,Fs,'power','reassigned')';

        % Cross PSD Calculation
        cross_psd = abs(cpsd(signal1_d,signal2_d,segment_window,0.9*frameLength,npoint,Fs))';

        % PLDNE Calculation Normalize
        phi_PLDNE = abs((psdsignal1 - psdsignal2)./(psdsignal1 + psdsignal2));
        % Delta PLD Calculation
        phi_PLD = max((psdsignal1 - psdsignal2),0);

        if k == 1
            phi_noise(k,:) = (1-a2)*psdsignal1;
        else
            % For Every Freq. Bin
            for l = 4:length(fsignal)-10
                % Calculate Noise PSD Estimator
                if phi_PLDNE(l) < 0.2
                    phi_noise(k,l) = a2*phi_noise(k-1,l)+(1-a2)*psdsignal1(l);
                elseif phi_PLDNE(l) > 0.8
                    phi_noise(k,l) = phi_noise(k-1,l);
                else
                    phi_noise(k,l) = a3*phi_noise(k-1,l)+(1-a3)*psdsignal2(l);
                end

                % Calculate Transfer Function
                Tn1n2 = sinc(2*pi*fsignal(l)*mic_len/sound_speed);
                H12 = (cross_psd(l) - Tn1n2*phi_noise(k,l))/(psdsignal1(l) - phi_noise(k,l));
                H12 = min(H12,1);                   % Limit Value

                % Create Weiner Filter
                G_c = phi_PLD(l)/(phi_PLD(l)+gramma*(1-abs(H12)^2)*phi_noise(k,l));
                G_c = min(max(G_c,0),1);            % Limit Value
                G(k,l) = G_c;                       % Append Data
            end

        end
        
        % Specctral Subpression 
        en_signal_f_m = G(k,:).*signal1_ft_m;
        en_signal_f = en_signal_f_m.*exp(complex(0,1).*signal1_ft_ph);

        % Inverse Fourier Transform
        en_signal_d = ifft(en_signal_f,npoint);
        en_signal_c = en_signal_d(1:frameLength);
        en_signal_c = real(en_signal_c);
        
        % Export PLD Enhance Signal
        en_signal(currentFrame:nextFrame) = en_signal_c;
        
        %% Harmonic Regeneration
        harmo = max(en_signal_c,0);

        % Harmonic Regeneration Fourier Transform
        HR_fft = fft(harmo.*segment_window,npoint);
        HR_F = HR_fft(1:npoint/2+1);                % Single-Side Fourier
        HR_F(2:end-1) = 2*HR_F(2:end-1);
        HR_fft_abs = abs(HR_F);                     % Magnitude
        HR_fft_ph = angle(HR_F);                    % Phase

        % Smoothing Recursive
        HR_fft_M = q.*en_signal_f_m +2.*(1-q).*(HR_fft_abs);
        HR_fft_M(1) = 0;

        % Export Harmonic Regeneration Signal
        output_HR_fft = HR_fft_M.*exp(complex(0,1).*signal1_ft_ph);
        output_HR_ifft = ifft(output_HR_fft,npoint);
        output_HR_c = real(output_HR_ifft(1:frameLength));
        output_HR(currentFrame:nextFrame) = output_HR_c;

        

        %% Measurement Section
        Noise_PSD = periodogram(signal1_d - clean_signal,segment_window,npoint,Fs,'power','reassigned');
        EN_Noise_PSD = periodogram(signal1_d - output_HR_c,segment_window,npoint,Fs,'power','reassigned');
        Log_Err = 10*log10(Noise_PSD./EN_Noise_PSD);
        Log_Err(isnan(Log_Err)) = 0;
        Log_Err(isinf(Log_Err)) = 0;
        Log_Err_C = Log_Err_C+sum(abs(Noise_PSD - EN_Noise_PSD));

        % Update to Next Frame
        currentFrame = nextFrame-sp*frameLength+1;
        nextFrame = currentFrame + frameLength-1;

    end
    
    % Output SNR Value
    measured_signal = speech1(1:length(output_HR))';
    out_SNR(j) = 20*log10(norm(measured_signal)/norm(measured_signal - en_signal));
%         disp("Output SNR : "+num2str(out_SNR))

    logErr(j) = Log_Err_C;
%         disp("LOG Err : "+num2str(logErr))

    %% Plot Result
    % Plot Spectrogram of Signal
        spectrogram_graph = figure('Position', [10 10 1200 600]);
        sgtitle("Spectrogram of Speech "+speech+", "+noise_type+" noise SNR "+snr_level(j)+" dB")
        subplot(2,2,1)
        spectrogram_plot_sp(signal_c,npoint,Fs_c,"(a) Clean Signal");
        subplot(2,2,2)
        spectrogram_plot_sp(signal1,npoint,Fs_c,"(b) Signal from Mic 1");
        subplot(2,2,3)
        spectrogram_plot_sp(en_signal,npoint,Fs_c,"(c) Enhance Signal");
        subplot(2,2,4)
        spectrogram_plot_sp(output_HR,npoint,Fs_c,"(d) Harmonic Regeneration");
        saveas(spectrogram_graph,"/Users/napatr.tan/Documents/MATLAB/SeniorProject/Result/PLD_Result/ResultSpectrogramsp"+speech+noise_type+snr_level(j)+".jpg")

    %% Export Enhance Speech
        audiowrite(speech+"_"+noise_type+"_"+snr_level(j)+"dB_enhance.wav",output_HR,Fs)

end