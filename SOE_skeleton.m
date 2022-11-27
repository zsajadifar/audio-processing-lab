% Lab 1 for Digital Audio Signal Processing Lab Sessions
% Exercise 1-4: 3D audio
% 
% In this lab, we derive a set of filters g that can be used, along with
% the measured RIRs H, to produce the proper psychocoustic rendition of 3D
% audio
%
%

clear;
close all

% Load ATFs
load Computed_RIRs

% Load measured HRTFs
load HRTF 

% Define the signal length
siglength = 10;

% Load the speech signal and resample
source_filename{1} = 'speech1.wav';

% Noise flag for the noise perturbing SOE
noiseFlag = 1;
% Noise flag for sweetspot tests
sweetspotFlag = 0;

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

% Solve the SOE
if ~noiseFlag
    % Without noise
    g = H_rmv0\x_rmv0;
else
    % With noise
    noise_std = std(H_rmv0(:,1))*0.05;
    noise = randn(size(H_rmv0))*noise_std;
    g = (H_rmv0+noise)\x_rmv0;
end

% Plot estimated and real HRTFs
figure,
plot(1:size(x,1),x,'r');
hold on
plot(1:size(x,1),H*g,'b');
legend('x','H*g');

% Calculate synthesis error
synth_error=norm(H*g-x);
disp('synthesis error:')
disp(synth_error)

save('SOE_output.mat','g','H')

%% 1.4.6
% Synthethize the binaural speech using H and g and compare it
% (psychoacoustically) to the binaural speech synthetized with x

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

speech1_left  = fftfilt(HRTF_estimated_left,speech1);
speech1_right = fftfilt(HRTF_estimated_right,speech1);
binaural_sig  = [speech1_left,speech1_right];
soundsc(binaural_sig,fs_new)

%synthesized with x
HRTF_actual       = x;
HRTF_actual_left  = HRTF_actual(1:floor(numel(HRTF_actual)/2));
HRTF_actual_right = HRTF_actual(floor(numel(HRTF_actual)/2)+1:end);

speech1_left  = fftfilt(HRTF_actual_left,speech1);
speech1_right = fftfilt(HRTF_actual_right,speech1);
binaural_sig_synth_x  = [speech1_left,speech1_right];
soundsc(binaural_sig_synth_x,fs_new)

synth_error=norm(binaural_sig-binaural_sig_synth_x);
disp('compare binaural synthesized by actual and estimated HRTF:')
disp(synth_error)

%% 1.4.8
fs_new = 8000;
mic_length = 10;
L = mic_length*fs_new;
[speech1,fs1] = audioread('audio files/speech1.wav');
[speech2,fs2] = audioread('audio files/speech2.wav');
speech1      = resample(speech1,fs_new,fs1);
speech2      = resample(speech2,fs_new,fs2);

x1 = speech1(1:L);
x2 = speech2(1:L);

binaural_sig1_1 = [x1,x1];
binaural_sig2_1 = [x1,0.5*x1];
binaural_sig3_1 = [x1,delayseq(x1,3)];
binaural_sig4_1(:,1) =fftfilt(HRTF_estimated_left,binaural_sig1_1(:,1));
binaural_sig4_1(:,2) =fftfilt(HRTF_estimated_right,binaural_sig1_1(:,2));

binaural_sig1_2 = [x2,x2];
binaural_sig2_2 = [0.5*x2,x2];
binaural_sig3_2 = [delayseq(x2,3),x2];
binaural_sig4_2(:,1) =fftfilt(HRTF_estimated_right,binaural_sig1_2(:,1));
binaural_sig4_2(:,2) =fftfilt(HRTF_estimated_left,binaural_sig1_2(:,2));

binaural_sig1 = binaural_sig1_1 + binaural_sig1_2;
binaural_sig2 = binaural_sig2_1 + binaural_sig2_2;
binaural_sig3 = binaural_sig3_1 + binaural_sig3_2;
binaural_sig4 = binaural_sig4_1 + binaural_sig4_2;

soundsc(binaural_sig1,fs_new)
soundsc(binaural_sig2,fs_new)
soundsc(binaural_sig3,fs_new)
soundsc(binaural_sig4,fs_new)

%% 1.4.10
load Computed_RIRs.mat

fs_new = 8000;
mic_length = 10;
L = mic_length*fs_new;
[speech1,fs1] = audioread('audio files/speech1.wav');
speech1      = resample(speech1,fs_new,fs1);
speech1      = speech1(1:L);

sweetspotFlag=1;
xL = HRTF(:,1); % Left ear
xR = HRTF(:,2); % Right ear
Delta=ceil(sqrt(room_dim(1)^2+room_dim(2)^2)*fs_RIR/340);
xL_delayed = [zeros(Delta,1);xL(1:Lx);zeros(((Lh+Lg-1)-(Lx+Delta)),1)];
xR_delayed = [zeros(Delta,1);xR(1:Lx);zeros(((Lh+Lg-1)-(Lx+Delta)),1)];

x = [xL_delayed;xR_delayed];

error = zeros(10,1);

if sweetspotFlag==1
    for i=0:10
        m_pos(:,2)=m_pos(:,2)+([0.01;0.01]*i);
%         m_pos(:,1)=m_pos(:,1)+([0.01;0.01]*i);
        [RIR_sources_new,RIR_noise_new]=create_rirs(m_pos,s_pos,v_pos,room_dim,rev_time,fs_new,fs_new/2);
        RIR_sources_new = RIR_sources_new(1:Lh,:,:);
        HL=[];
        HR=[];
        for j = 1:J
            HL_ = toeplitz([RIR_sources_new(:,1,j);zeros(Lg-1,1)],zeros(Lg,1));
            HL = [HL HL_];
            HR_ = toeplitz([RIR_sources_new(:,2,j);zeros(Lg-1,1)],zeros(Lg,1));
            HR = [HR HR_];
        end
        
        H = [HL;HR];

        synth_error=norm(H*g-x);
        error(i+1)=synth_error;
       
        HRTF_estimated       = H*g;
        HRTF_estimated_left  = HRTF_estimated(1:floor(numel(HRTF_estimated)/2));
        HRTF_estimated_right = HRTF_estimated(floor(numel(HRTF_estimated)/2)+1:end);
        
        speech1_left  = fftfilt(HRTF_estimated_left,speech1);
        speech1_right = fftfilt(HRTF_estimated_right,speech1);
        binaural_sig_sweetspot(:,:,i+1)  = [speech1_left,speech1_right];
    end        
end

figure,
plot(0:10,error)
xlabel('move mics vertically [cm]')
ylabel('synthesis error')

% figure,
% plot(0:10,error)
% xlabel('move mics horizontally [cm]')
% ylabel('synthesis error')

soundsc(binaural_sig_sweetspot(:,:,10),fs_new)




