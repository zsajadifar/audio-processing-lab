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
Lg = 266;      % Length of filters g_j

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
% loaded HRTF
xL = HRTF(:,1); % Left ear
xR = HRTF(:,2); % Right ear

% Construct H (from HL and HR) and x (from xL and xR) and remove all-zero rows in H, and the corresponding elements in x
Delta=ceil(sqrt(room_dim(1)^2+room_dim(2)^2)*fs_RIR/340);
xL_delayed = [zeros(Delta,1);xL(Delta+1:(Lg+Lh-1),1)];
xR_delayed = [zeros(Delta,1);xR(Delta+1:(Lg+Lh-1),1)];

H = [HL;HR];
x = [xL_delayed;xR_delayed];

% remove zero rows
idx_nonzeros = sum(abs(H),2)> 0;
H = H(idx_nonzeros,:);
x = x(idx_nonzeros,:);

% Solve the SOE
if ~noiseFlag
    % Without noise
    g = H\x;
else
    % With noise
    g = H\x;
end

% Plot estimated and real HRTFs
figure,
hold on
plot(1:size(x,1),x,'r');
plot(1:size(x,1),H*g,'b');
legend('x','H*g');

% Calculate synthesis error
synth_error=norm(H*g-x);
disp(synth_error);

% Synthethize the binaural speech using H and g and compare it
% (psychoacoustically) to the binaural speech synthetized with x


