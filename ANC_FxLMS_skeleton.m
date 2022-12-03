% Session 4 - Exercise 1 - Single-channel FxNLMS ANC
% 
% Main points:
% (1) Generate the noisy microphone signal at the left ear
% (2) Implement the FxNLMS ANC.
% (3) Compute input/output SNRs and SNR improvement.
% (4) Implement the filter for left and right ears. Hear the synthesised signal

clear all
close all
warning off 

% Load RIRs
load Computed_RIRs.mat

% Set length
sigLenSec = 10;

% Select the speaker number
speaker = 1;


%%
% Plot the RIRs of the noise
figure(1); clf;
subplot(2,1,1)
plot(RIR_noise(:,1))
title('RIR of noise to left mic')
subplot(2,1,2)
plot(RIR_noise(:,2))
title('RIR of noise to right mic')

%%
% Read in the noise source (resample if necessary)
[noise,fs_noise] = audioread('audio files/White_noise1.wav');
noise      = resample(noise,fs_RIR,fs_noise);
noise      = noise(1:sigLenSec*fs_RIR);

noise_mic_L = fftfilt(RIR_noise(:,1),noise);
noise_mic_R = fftfilt(RIR_noise(:,2),noise);

filt_noise = noise_mic_L; % just consoder the left mic

% Plot the noisy signal
figure(2); clf;
subplot(2,1,1)
plot(noise_mic_L)
title('noisy signal in left mic')
subplot(2,1,2)
plot(noise_mic_R)
title('noisy signal in right mic')

%%
% speech signal
[speech,fs_speech] = audioread('audio files/speech1.wav');
speech      = resample(speech,fs_RIR,fs_speech);
speech      = speech(1:sigLenSec*fs_RIR);

speech_mic_L = fftfilt(RIR_sources(:,1,speaker),speech);
speech_mic_R = fftfilt(RIR_sources(:,2,speaker),speech);


mic = speech_mic_L + filt_noise; 
% soundsc(mic,fs_RIR)

%% FxLMS

M = 100;  % Length of secondary path (RIR)
L = 600;  % Adaptive filter length
sigLenSample = sigLenSec*fs_RIR;

mu = 0.2;   % Step size
delta = 5*10^(-5);

W = zeros(L,1); % Initialize adaptive filter
W_just_noise = zeros(L,1);
e = zeros(sigLenSample,1);
e_just_noise = zeros(sigLenSample,1);
x = [zeros(L+M-1,1);noise];
speech_zeropad = [zeros(L+M-1,1);speech];

h = RIR_sources(1:M,1,speaker);% consider first mic and speaker

tic
for n = 1:sigLenSample

    % STEP 1 : Arrange the previous L + M − 1 samples of x(n) up until the
    % current sample, n, in an [M x L] Hankel matrix X_Hmat (this will be
    % use for the filtering operation)

    xseg = flip(x(n:n+L+M-1)); %noise segment
    X_Hmat = hankel(xseg(1:M),xseg(M+1:end));

    % STEP 2: Filter the signal x(n) with the adaptive filter w(n) to
    % generate loudspeaker signals y(n). Store the samples of y(n) into the
    % vector y

    y = X_Hmat*W + speech_zeropad(n:n+M-1); %Mx1
    y_just_noise = X_Hmat*W_just_noise ; %just noise, without desired speech

    % STEP 3: Compute the error signal e(n)
    e(n) = filt_noise(n) + h'*y;
    e_just_noise(n) = filt_noise(n) + h'*y_just_noise;

    % STEP 4: Pre-filter x(n) with the (estimated) RIRs to obtain the
    % filtered x(n), \bar{x}(n). Store the samples of \bar{x}(n) into
    % the vector xf
    xf = X_Hmat'*h; 
    
    % STEP 5: Update the filter w(n). Store the samples of w(n) into
    % the vector w
    W = W - (mu/(norm(xf)^2+delta))*xf*e(n);
    W_just_noise = W_just_noise - (mu/(norm(xf)^2+delta))*xf*e_just_noise(n);
end
toc


%%
% Calculate the noise suppression
NSP = 10*log10(mean(e.^2)/mean(filt_noise.^2));


%%
% In the existing plot of the noisy signal, superimpose the error signal to
% appreciate the noise cancellation

figure; hold on;
plot(filt_noise)
plot(e)
legend('d(n)','e(n)')
title('Compare d and e (with desired speech)')

figure; hold on;
plot(filt_noise)
plot(e_just_noise)
legend('d(n)','e(n)')
title('Compare d and e (without desired speech)')

% soundsc(e,fs_RIR)

%% Compare SNR 
SNR_before_channel = 10*log10(var(speech)/var(noise));
SNR_mic_no_filter  = 10*log10(var(speech_mic_L)/var(filt_noise));
SNR_mic_filtered   = 10*log10(var(speech_mic_L)/var(e_just_noise));

