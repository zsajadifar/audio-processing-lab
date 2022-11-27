close all
clear

load('Computed_RIRs.mat');
speechfilename{1} = 'audio files/speech1.wav';
%speechfilename{2} = 'audio files/speech2.wav';
noisefilename{1}  = 'audio files/Babble_noise1.wav';

mic_length = 10; % desired length of microphone signals in Sec
mic_num = size(RIR_sources,2);
[mic] = create_micsigs(speechfilename,noisefilename,mic_length,mic_num,fs_RIR,RIR_sources,RIR_noise);
figure,
plot(mic(:,1),'r')
hold on
plot(mic(:,2),'b')

soundsc(mic(:,1),fs_RIR)


