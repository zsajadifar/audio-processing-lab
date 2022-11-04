dirac_IR = zeros(100,1);
dirac_IR(1)= 1;
nfft=128;
[speech1,fs1] = audioread('audio files/speech1.wav');
speech1 = speech1(1:10:fs1);
y = OLA_skeleton(speech1,dirac_IR,nfft);

error=norm(y-speech1);
disp('error:')
disp(error)

figure,
plot(speech1,'r');
hold on
plot(y,'b');
legend('actual signal','OLA signal');

%% compare OLA with conv
load SOE_output.mat

% read audio
fs_new = 8000;
mic_length = 10;
L = mic_length*fs_new;
[speech1,fs1] = audioread('audio files/speech1.wav');
speech1      = resample(speech1,fs_new,fs1);
speech1      = speech1(1:L);

%conv
HRTF_estimated       = H*g;
HRTF_estimated_left  = HRTF_estimated(1:floor(numel(HRTF_estimated)/2));
HRTF_estimated_right = HRTF_estimated(floor(numel(HRTF_estimated)/2)+1:end);
disp('conv: ')
tic
conv_L = conv(speech1,HRTF_estimated_left);
conv_R = conv(speech1,HRTF_estimated_right);
toc
binaural_sig_conv = [conv_L conv_R];

%OLA
nfft = 1024;% nfft should be larger than Lh
disp('OLA: ')
tic
OLA_L = OLA_skeleton(speech1,HRTF_estimated_left ,nfft);
OLA_R = OLA_skeleton(speech1,HRTF_estimated_right,nfft);
toc
binaural_sig_OLA = [OLA_L,OLA_R];

soundsc(binaural_sig_conv,fs_new)
soundsc(binaural_sig_OLA,fs_new)

