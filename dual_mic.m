function [speech1,speech2,noise1,noise2]=dual_mic(signal_c,signal_n)

    % System Parameter
    SNRNoise = 0;                   % SNR between 2 noise
    SNRSpeech = 10;                 % SNR bwtween 2 speech

    %Add Audio Path
    addpath('/Users/napatr.tan/Documents/MATLAB/SeniorProject/Data')
    
    % Assign Signal
    speech1 = signal_c;             % Clean Speech
    noise2 = signal_n-signal_c;     % Noise
    
    % Random Impulse
    h = load('imp256.mat');
    h1 = h.h1_256;
    h2 = h.h2_256;
    
    % Filter and Scale Speech Signal
    speech2 = filter(h1,1,speech1);
    coeffSpeech = sqrt(var(speech1)/(var(speech2)*10^(SNRSpeech/10)));
    speech2 = speech2*coeffSpeech;

    % Filter and Scale Noise Signal
    noise1 = filter(h2,1,noise2);
    coeffNoise = sqrt(var(noise2)/(var(noise1)*10^(SNRNoise/10)));
    noise1 = noise1*coeffNoise;
    noise1 = 0.1*noise1 + 0.9*noise2;
    
end