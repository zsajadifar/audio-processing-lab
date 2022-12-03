% Session 4 - Exercise 2 - Multi-channel FxNLMS ANC
%           - Exercise 3 - Multi-channel FxNLMS ANC in 3D audio scenario
% 
%
% Main points:
% (1) Generate the noisy microphone signals
% (2) Implement the ANC.
% (3) Compute input/output SNRs and SNR improvement
% (4) Implement the filter for left and right ears. Hear the synthesised signal

clear all;
close all
warning off 

% Load RIRs
load Computed_RIRs.mat

% Set length
sigLenSec = 5;

%%
% Read in the noise source (resample if necessary)
[noise,fs_noise] = audioread('audio files/White_noise1.wav');
noise      = resample(noise,fs_RIR,fs_noise);
noise      = noise(1:sigLenSec*fs_RIR);

noise_mic_L = fftfilt(RIR_noise(:,1),noise);
noise_mic_R = fftfilt(RIR_noise(:,2),noise);

filt_noise = [noise_mic_L,noise_mic_R]; 

% Plot the noisy signal
figure,
hold on
plot(noise_mic_L)
plot(noise_mic_R)
title('noise signal')
legend('left mic','right mic')

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
% Boolean flag 3D audio scenario
AUDIO_3D = 1;
if AUDIO_3D
    load HRTF.mat

    % Define the lengths of RIRs and g filters
    J = size(RIR_sources,3);
    Lh = 400; 
    Lg = ceil((2*Lh-2)/(J-2));  % Length of filters g
    
    % Truncate impulse response to reduce computational complexity
    RIR_sources = RIR_sources(1:Lh,:,:);
    
    % Calculate delay for SOE
    Delta=ceil(sqrt(room_dim(1)^2 + room_dim(2)^2)*fs_RIR/340);
    
    % Define the Toeplitz matrices for left and right ear (left side of SOE)
    HL=[];
    HR=[];
    for j = 1:J
        HL_ = toeplitz([RIR_sources(:,1,j);zeros(Lg-1,1)],zeros(Lg,1));
        HL = [HL HL_];
        HR_ = toeplitz([RIR_sources(:,2,j);zeros(Lg-1,1)],zeros(Lg,1));
        HR = [HR HR_];
    end
    
    % Define the HRTFs for left and right ear (right side of SOE) from the
    % loaded HRTF, make sure that length of HRTF is smaller than Lh
    Lx = Lh/2;
    xL = HRTF(:,1); % Left ear
    xR = HRTF(:,2); % Right ear
    
    % Construct H (from HL and HR) and x (from xL and xR) and remove all-zero 
    % rows in H, and the corresponding elements in x. Truncate the length of xL
    % and xR to the Lx and add zeros at the end to match their size to the
    % Lh+Lg-1
    xL_delayed = [zeros(Delta,1);xL(1:Lx);zeros(((Lh+Lg-1)-(Lx+Delta)),1)];
    xR_delayed = [zeros(Delta,1);xR(1:Lx);zeros(((Lh+Lg-1)-(Lx+Delta)),1)];
    
    H = [HL;HR];
    x = [xL_delayed;xR_delayed];
    
    % remove zero rows
    idx_nonzeros = sum(abs(H),2)> 0;
    H_rmv0 = H(idx_nonzeros,:);
    x_rmv0 = x(idx_nonzeros,:);
    
    % solve the SOE
    g = H_rmv0\x_rmv0;
    
    L = sigLenSec*fs_RIR;
    [speech1,fs1] = audioread('audio files/speech1.wav');
    speech1      = resample(speech1,fs_RIR,fs1);
    speech1      = speech1(1:L);
    
    %synthesized with H and g
    HRTF_estimated       = H*g;
    HRTF_estimated_left  = HRTF_estimated(1:floor(numel(HRTF_estimated)/2));
    HRTF_estimated_right = HRTF_estimated(floor(numel(HRTF_estimated)/2)+1:end);
    
    speech_L  = fftfilt(HRTF_estimated_left,speech1);
    speech_R  = fftfilt(HRTF_estimated_right,speech1);
    binaural_sig    = [speech_L,speech_R];

    % 0 dB SNR between left binaural_sig and left filt_noise
    SNR_babbel = 0;
    scale = var(speech_L) / ((10.^(SNR_babbel/10)).*var(filt_noise(:,1)));
    binaural_sig_scaled = binaural_sig/sqrt(scale);
end





%% MFxLMS

M = 100;  % Length of secondary path (RIR)
L = 600;  % Adaptive filter length
J = size(RIR_sources,3); %number of speakers
sigLenSample = sigLenSec*fs_RIR;

mu = 0.2;   % Step size
delta = 5*10^(-5);

W = zeros(L,J); % Initialize adaptive filter

x = [zeros(L+M-1,1);noise];

h_L = squeeze(RIR_sources(1:M,1,:));
h_R = squeeze(RIR_sources(1:M,2,:));

e_L = zeros(sigLenSample,1);
e_R = zeros(sigLenSample,1);

tic
for n=1:sigLenSample
    
    % STEP 1 : Arrange the previous L + M − 1 samples of x(n) up until the
    % current sample, n, in an [M x L] Hankel matrix X_Hmat (this will be
    % use for the filtering operation). Store the samples in the
    % appropriate form.
    
    xseg = flip(x(n:n+L+M-1)); %noise segment
    X_Hmat = hankel(xseg(1:M),xseg(M+1:end));

    % STEP 2: Filter the signal x(n) with the adaptive filter w(n) to
    % generate loudspeaker signals y(n). Store the samples in
    % the appropriate form.
    
    y = zeros(M,J);
    for j = 1:J
        y(:,j) = X_Hmat*W(:,j);
    end

    % STEP 3: Compute the error signals e_L(n) and e_R(n). Store them in
    % the appropriate form.

    e_L(n)=filt_noise(n,1);
    e_R(n)=filt_noise(n,2);

    for j = 1:J
        e_L(n)=e_L(n)+ h_L(:,j)'*y(:,j);
        e_R(n)=e_R(n)+ h_R(:,j)'*y(:,j);
    end

    if AUDIO_3D
        e_L(n) = e_L(n)+ binaural_sig_scaled(n,1);
        e_R(n) = e_R(n)+ binaural_sig_scaled(n,2);
    end
    
    % STEP 4: Pre-filter x(n) with the (estimated) RIRs to obtain the
    % filtered x(n), \bar{x}(n). Store the samples in the appropriate form.
    for j = 1:J
        xf_L(:,j)=X_Hmat'*h_L(:,j);
        xf_R(:,j)=X_Hmat'*h_R(:,j);
    end

    % STEP 5: Update the filter w(n). Store the samples in the appropriate
    % form.
    E = [e_L(n),e_R(n)];
    for j=1:J
        X = [xf_L(:,j),xf_R(:,j)];
        W(:,j) = W(:,j) - (mu/(norm(X)^2+delta))*X*E';
    end 
end
toc

%%
% Calculate the noise suppression
% Calculate the noise suppression
NSP_L = 10*log10(mean(e_L.^2)/mean(filt_noise(:,1).^2));
NSP_R = 10*log10(mean(e_R.^2)/mean(filt_noise(:,2).^2));

%%
% In the existing plot of the noisy signals, superimpose the corresponding
% error signals to appreciate the noise cancellation

figure; hold on;
plot(filt_noise(:,1))
plot(e_L)
title('Compare d and e (left ear)')
if AUDIO_3D
    plot(binaural_sig_scaled(:,1))
    legend('d(n)','e(n)','binaural left')

else 
    legend('d(n)','e(n)')
end



figure; hold on;
plot(filt_noise(:,2))
plot(e_R)
title('Compare d and e (right ear)')
if AUDIO_3D
    plot(binaural_sig_scaled(:,2))
    legend('d(n)','e(n)','binaural right')

else 
    legend('d(n)','e(n)')
end


