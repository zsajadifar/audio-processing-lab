close all
clear 

fs_new = 8000;
mic_length = 10;
L = mic_length*fs_new;
[speech,fs] = audioread('audio files/speech1.wav');
speech      = resample(speech,fs_new,fs);
x = speech(1:L);
binaural_sig1_1 = [x,x];
binaural_sig2_1 = [x,0.5*x];
binaural_sig3_1 = []