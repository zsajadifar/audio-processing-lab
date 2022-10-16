close all

load('HRTF.mat')
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
binaural_sig4_1(:,1) =fftfilt(HRTF(:,1),binaural_sig1_1(:,1));
binaural_sig4_1(:,2) =fftfilt(HRTF(:,2),binaural_sig1_1(:,2));

binaural_sig1_2 = [x2,x2];
binaural_sig2_2 = [0.5*x2,x2];
binaural_sig3_2 = [delayseq(x2,3),x2];
binaural_sig4_2(:,1) =fftfilt(HRTF(:,2),binaural_sig1_2(:,1));
binaural_sig4_2(:,2) =fftfilt(HRTF(:,1),binaural_sig1_2(:,2));

binaural_sig1 = binaural_sig1_1 + binaural_sig1_2;
binaural_sig2 = binaural_sig2_1 + binaural_sig2_2;
binaural_sig3 = binaural_sig3_1 + binaural_sig3_2;
binaural_sig4 = binaural_sig4_1 + binaural_sig4_2;

soundsc(binaural_sig1,fs_new)
soundsc(binaural_sig2,fs_new)
soundsc(binaural_sig3,fs_new)
soundsc(binaural_sig4,fs_new)