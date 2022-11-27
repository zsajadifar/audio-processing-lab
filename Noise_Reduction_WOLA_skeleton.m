% Lab 3 for Digital Audio Signal Processing Lab Sessions
% Session 3: Noise reduction in the STFT domain
% R.Ali, G. Bernardi, J.Schott, A. Bertrand
% 2020
%
% The following is the skeleton code for doing the noise reduction in the
% STFT domain. Complete the missing sections accordingly


%% Exercise 3.1: ## Obtain the noisy microphone signals as outlined in the session document.
%
% SOME CONVENTIONS:
%
% Your clean speech signal should be called "speech" and is a binaural signal such that:
% speech = [speech_L speech_R]  % speech_L and speech_R are the clean binaurally synthesised speech signals from session 1
% 
% Your noises should be stored in one binaural variable called "noise"


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Section of code to complete %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all

%% Exercise 3.1

% Load ATFs
load Computed_RIRs

% Load measured HRTFs
load HRTF 

% Define the signal length
mic_length = 10;

% Load the speech signal and resample
source_filename{1} = 'speech1.wav';
noise_filename{1} = 'Babble_noise1.wav';

% Number of loudspeakers
J = size(RIR_sources,3);

% Define the lengths of RIRs and g filters
Lh = 400; % Length of impulse responses of listeniLg room
Lg = ceil((2*Lh-2)/(J-2));      % Length of filters g_j

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
Delta=ceil(sqrt(room_dim(1)^2+room_dim(2)^2)*fs_RIR/340);
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

fs_new = 8000;
mic_length = 10;
L = mic_length*fs_new;
[speech1,fs1] = audioread('audio files/speech1.wav');
speech1      = resample(speech1,fs_new,fs1);
speech1      = speech1(1:L);

%synthesized with H and g
HRTF_estimated       = H*g;
HRTF_estimated_left  = HRTF_estimated(1:floor(numel(HRTF_estimated)/2));
HRTF_estimated_right = HRTF_estimated(floor(numel(HRTF_estimated)/2)+1:end);

speech_L  = fftfilt(HRTF_estimated_left,speech1);
speech_R = fftfilt(HRTF_estimated_right,speech1);
speech  = [speech_L,speech_R];

%Babble noise
[babble_noise,fs_noise] = audioread('audio files/Babble_noise1.wav');
babble_noise      = resample(babble_noise,fs_new,fs_noise);
babble_noise      = babble_noise(1:L);

babble_noise_mic_L = fftfilt(RIR_noise(:,1),babble_noise);
babble_noise_mic_R = fftfilt(RIR_noise(:,2),babble_noise);

SNR_babbel = 10*log10(var(speech_L)/var(2.3*babble_noise_mic_L));
binaural_babble_noise = 2.3*[babble_noise_mic_L,babble_noise_mic_R];

%uncorrelated noise
uncorr_noise = randn(L,1);
uncorr_noise_mic_L = fftfilt(RIR_noise(:,1),uncorr_noise);
uncorr_noise_mic_R = fftfilt(RIR_noise(:,2),uncorr_noise);

SNR_uncorr = 10*log10(var(speech_L)/var(0.0045*uncorr_noise_mic_L));
binaural_uncorr_noise = 0.0045*[uncorr_noise_mic_L,uncorr_noise_mic_R];

% combine speech signal with 2 noises
noise = binaural_uncorr_noise + binaural_babble_noise;


% Create noisy mic signals in the time domain:
y_TD = speech + noise;  % stacked microphone signals
% soundsc(y_TD,fs_new)


%% Apply WOLA analysis to observe signals in the STFT domain, Apply the SPP.

fs = fs_RIR;    % sampling freq
nfft = 512;    % number of DFT points
window = sqrt(hann(nfft,'periodic')); % analysis window
noverlap = 2;   % factor for overlap. noverlap = 2 corresponds to 50% overlap
% time = 0:1/fs:((length(x)-1)/fs);


% ## Apply the WOLA analysis to the noisy mic. signals, the speech, and the noise.
g=[];
[y_STFT,f] = WOLA_analysis_skeleton(y_TD,fs,window,nfft,noverlap,g);
[n_STFT,~] = WOLA_analysis_skeleton(noise,fs,window,nfft,noverlap,g);
[x_STFT,~] = WOLA_analysis_skeleton(speech,fs,window,nfft,noverlap,g);


% Observe the STFT
clow = -60; chigh = 10; % lower and upper limits for signal power in spectrogram (can change for different resolutions)
[N_freqs, N_frames] = size(y_STFT(:,:,1));

figure; 
imagesc(1:N_frames,f/1000, mag2db(abs(x_STFT(:,:,1))), [clow, chigh]); 
colorbar; 
axis xy; 
set(gca,'fontsize', 14);
set(gcf,'color','w'); 
xlabel('Time Frame'); 
ylabel('Frequency (kHz)')
title('STFT, only speech, left mic')

figure; 
imagesc(1:N_frames,f/1000, mag2db(abs(n_STFT(:,:,1))), [clow, chigh]); 
colorbar; 
axis xy; 
set(gca,'fontsize', 14);
set(gcf,'color','w'); 
xlabel('Time Frame'); 
ylabel('Frequency (kHz)')
title('STFT, only noise(Babble+uncorr), left mic')

figure; 
imagesc(1:N_frames,f/1000, mag2db(abs(y_STFT(:,:,1))), [clow, chigh]); 
colorbar; 
axis xy; 
set(gca,'fontsize', 14);
set(gcf,'color','w'); 
xlabel('Time Frame'); 
ylabel('Frequency (kHz)')
title('STFT, noisy speech signal, left mic')


% ## Compute the Speech Presence Probability on the reference microphone
% (you can try the speech-only signal or the noisy-speech in one of the microphones)
% Use the attached spp_calc.m function
y_TD_babble = speech + binaural_babble_noise;
y_TD_uncorr = speech + binaural_uncorr_noise;

[noisePowMat, SPP_speech]      = spp_calc(speech(:,1),nfft,nfft/noverlap);
[noisePowMat, SPP_noise]       = spp_calc(noise(:,1),nfft,nfft/noverlap);
[noisePowMat, SPP_y_TD]        = spp_calc(y_TD(:,1),nfft,nfft/noverlap);
[noisePowMat, SPP_y_TD_babble] = spp_calc(y_TD_babble(:,1),nfft,nfft/noverlap);
[noisePowMat, SPP_y_TD_uncorr] = spp_calc(y_TD_uncorr(:,1),nfft,nfft/noverlap);


% Observe the SPP
figure; imagesc(1:N_frames, f,SPP_speech); colorbar; axis xy; set(gcf,'color','w');  
set(gca,'Fontsize',14), xlabel('Time Frames'), ylabel('Frequency (Hz)'), title('SPP for ref mic,speech');

figure; imagesc(1:N_frames, f,SPP_noise); colorbar; axis xy; set(gcf,'color','w');  
set(gca,'Fontsize',14), xlabel('Time Frames'), ylabel('Frequency (Hz)'), title('SPP for ref mic,noise');

figure; imagesc(1:N_frames, f,SPP_y_TD); colorbar; axis xy; set(gcf,'color','w');  
set(gca,'Fontsize',14), xlabel('Time Frames'), ylabel('Frequency (Hz)'), title('SPP for ref mic,speech+noise');

figure; imagesc(1:N_frames, f,SPP_y_TD_babble); colorbar; axis xy; set(gcf,'color','w');  
set(gca,'Fontsize',14), xlabel('Time Frames'), ylabel('Frequency (Hz)'), title('SPP for ref mic,speech+babble noise');

figure; imagesc(1:N_frames, f,SPP_y_TD_uncorr); colorbar; axis xy; set(gcf,'color','w');  
set(gca,'Fontsize',14), xlabel('Time Frames'), ylabel('Frequency (Hz)'), title('SPP for ref mic,speech+uncorr noise');


%%  Exercise 3.2: ## Implementation of the MWF
num_mics = 2;

Rnn = cell(N_freqs,1);  Rnn(:) = {1e-6*ones(num_mics,num_mics)};      % Noise Only (NO) corr. matrix. Initialize to small random values
Ryy = cell(N_freqs,1);  Ryy(:) = {1e-6*ones(num_mics,num_mics)};      % Speech + Noise (SPN) corr. matrix. Initialize to small random values
lambda = 0.995;                                                       % Forgetting factors for correlation matrices - can change
SPP_thr = 0.95;                                                       % Threshold for SPP - can change

% For MWF filter for left ear
S_mvdr_mwfL_stft = zeros(N_freqs,N_frames);         
X_mvdr_mwfL_stft = zeros(N_freqs,N_frames);         
N_mvdr_mwfL_stft = zeros(N_freqs,N_frames);  
W_mvdr_mwfL = (1/num_mics)*ones(num_mics,N_freqs); 


% STFT Processing
% Looping through each time frame and each frequency bin
tic
for l=2:N_frames % Time index
    
    for k = 1:N_freqs % Freq index
        
        
        % Create a vector of mic. signals
        Y_kl = squeeze(y_STFT(k,l,1:num_mics));  % M x 1 noisy mic sigs for this freq and time frame
        X_kl = squeeze(x_STFT(k,l,1:num_mics));
        N_kl = squeeze(n_STFT(k,l,1:num_mics));

        % ## Update the correlation matrices using the forgetting factor.
        % Threshold the SPP in order to distinguish between periods of speech and non-speech activity
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Section of code to complete (3 lines) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % chapter 4, slide 32
        if SPP_y_TD > SPP_thr
            Ryy{k} = ((lambda.^2) * Ryy{k}) + (1-lambda.^2)*(Y_kl*Y_kl');
        else
            Rnn{k} = ((lambda.^2) * Rnn{k}) + (1-lambda.^2)*(N_kl*N_kl');
        end
        
        % Computing the MWF filter using a generalised eigenvalue
        % decomposition of the correlation matrices.
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Section of code to complete ~ 10 lines %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % chapter 4, slide 36,37
        [V,D] = eig(Ryy{k},Rnn{k});
        [B,I] = sort(diag(D),'descend');
        Q_inv_H = V(:,I);% V = Q^(-H)
        %Q_H = inv(Q_inv_H);
       
        sigma_SPN_NO = B(1);
        B(1) = 1-(1/sigma_SPN_NO);
        B(2:end)=0;

        F = (Q_inv_H * diag(B))/Q_inv_H; %Q_inv_H * diag([(1-(1/sigma_SPN_NO)),0])*Q_H


       
        W_mvdr_mwfL(:,k) = F(:,1); % Final expression for filter
        
        % Filtering the noisy speech, the speech-only, and the noise-only.
        S_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* Y_kl(1:num_mics);
        X_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* X_kl(1:num_mics);
        N_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* N_kl(1:num_mics);
        

        
    end % end freqs
end % end time frames
toc


%% Observe processed STFTs

figure; imagesc(1:N_frames,f/1000,mag2db(abs(y_STFT(:,:,1))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('microphne signal, 1st mic');
figure; imagesc(1:N_frames,f/1000,mag2db(abs(S_mvdr_mwfL_stft(:,:,1))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'),title('Enhanced Signal L - MWF');


%% Apply the synthesis stage of the WOLA framework to obtain the time domain equivalents:

s_mwfL = WOLA_synthesis_skeleton(S_mvdr_mwfL_stft,window,nfft,noverlap); %To complete (time-domain version of S_mvdr_mwfL_stft)
x_mwfL = WOLA_synthesis_skeleton(X_mvdr_mwfL_stft,window,nfft,noverlap); % To complete (time-domain version of X_mvdr_mwfL_stft) 
n_mwfL = WOLA_synthesis_skeleton(N_mvdr_mwfL_stft,window,nfft,noverlap); % To complete (time-domain version of N_mvdr_mwfL_stft)

% PLOT SIGNALS
% LISTEN TO SIGNALS!
figure, hold on
plot(s_mwfL,'r')
plot(x_mwfL,'g')
plot(n_mwfL,'b')
legend('s','x','n')



%% EVALUATION

% SNR_in = % Compute input SNR
% SNR_out = % Compute output SNR


