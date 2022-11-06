clear 
close all

%% Load filter from SOE
load SOE_output.mat
%% read audio
fs_new = 8000;
mic_length = 10;
L = mic_length*fs_new;
[speech1,fs1] = audioread('audio files/speech1.wav');
speech1      = resample(speech1,fs_new,fs1);
speech1      = speech1(1:L);

%% Analysis
noverlap = 2;
nfft = 2048;
win_analysis  = sqrt(hann(nfft,'periodic'));
win_synthesis = sqrt(hann(nfft,'periodic'));
x = speech1;
x = repmat(x,1,5);

[X,f] = WOLA_analysis_skeleton(x,fs_new,win_analysis,nfft,noverlap,g);

X_pow = abs(X).^2; % *2 to consider 2-sided power
X_pow_2sided = X_pow;
X_pow_2sided(2:end,:,:) = X_pow(2:end,:,:)*2;
X_pow_2sided_dB = 10*log10(X_pow_2sided);
X_pow_2sided_dB(X_pow_2sided_dB <= 10^(-10)) = 0;

[freq_bins, time_bins] = meshgrid(1:size(X,2),f);
figure,mesh(freq_bins,time_bins,X_pow_2sided_dB)
grid("off")
colormap jet
xlabel('Time window')
ylabel('Frequency (Hz)')
title('Output of WOLA analysis')
cb = colorbar;
cb.Title.String = "Power (dB)";

%% Synthesis
x_synth = WOLA_synthesis_skeleton(X,win_synthesis,nfft,noverlap);

figure,
plot(x_synth,'b')
hold on
plot(x,'r')
legend('x synth','x');

synth_error=norm(x_synth-x(1:numel(x_synth)));
disp('compare original speech signal and reconstructed one by WOLA_synthesis:')
disp(synth_error)

soundsc(x,fs_new)
soundsc(x_synth,fs_new)

%% Perfect reconstruction condition: 
% y(n) = sum_m[x(n).*win_anal(n-mN/2).*win_synth(n-mN/2)]
% p(n)+p(n+N/2)=1  n=0,...,N/2 -1

for n=1:nfft/2
    S(n)= win_analysis(n)*win_synthesis(n) + win_analysis(n+nfft/2)*win_synthesis(n+nfft/2); 
end

if sum(S)==numel(S)
    disp('Perfect reconstruction condition is satisfied')
else
    disp('Perfect reconstruction condition is not satisfied')
end

