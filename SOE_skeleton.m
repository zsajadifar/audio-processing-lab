% Lab 1 for Digital Audio Signal Processing Lab Sessions
% Exercise 1-4: 3D audio
% 
% In this lab, we derive a set of filters g that can be used, along with
% the measured RIRs H, to produce the proper psychocoustic rendition of 3D
% audio
%
%

clear;
% close all

% Load ATFs
load Computed_RIRs

% Load measured HRTFs
load HRTF 

% Define the signal length
siglength = 10;

% Load the speech signal and resample
source_filename{1} = 'speech1.wav';

% Noise flag for the noise perturbing SOE
noiseFlag = 0;
% Noise flag for sweetspot tests
sweetspotFlag = 0;

% Number of loudspeakers
J = size(RIR_sources,3);

% Define the lengths of RIRs and g filters
Lh = 400; % Length of impulse responses of listeniLg room
% Lg = ;      % Length of filters g_j

% Truncate impulse response to reduce computational complexity

% Calculate delay for SOE

% Define the Toeplitz matrices for left and right ear (left side of SOE)
HL=[];
HR=[];

% Define the HRTFs for left and right ear (right side of SOE) from the
% loaded HRTF
% xL = ; % Left ear
% xR = ; % Right ear

% Construct H (from HL and HR) and x (from xL and xR) and remove all-zero rows in H, and the corresponding elements in x


% Solve the SOE
if ~noiseFlag
    % Without noise
else
    % With noise
end

% Plot estimated and real HRTFs
figure(1); clf;

% Calculate synthesis error

% Synthethize the binaural speech using H and g and compare it
% (psychoacoustically) to the binaural speech synthetized with x


