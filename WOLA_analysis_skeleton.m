% Lab 2 for Digital Audio Signal Processing Lab Sessions
% Session 2: Binaural synthesis and 3D audio: OLA and WOLA frameworks
% R.Ali, G. Bernardi, J.Schott, A. Bertrand
% 2021
%
% The following is the skeleton code for the analysis stage of the WOLA method, which you need to
% complete


function [X,f] = WOLA_analysis_skeleton(x,fs,window,nfft,noverlap,g)
%WOLA_analysis  short-time fourier transform
% INPUT:
%   x           : input time signal(s) (samples x channels)
%   fs          : sampling rate
%   window      : window function
%   nfft        : FFT size
%   noverlap    : frame overlap; default: 2 (50%)
%   g           : filter(a vector with size (Lg*M)*1)
%
% OUTPUT:
%   X           : STFT matrix (bins x frames x channels)
%   f           : frequency vector for bins


% use only half FFT spectrum
N_half = nfft / 2 + 1;

% get frequency vector
f = 0:(fs / 2) / (N_half - 1):fs / 2;

% init
L = floor((length(x) - nfft + (nfft / noverlap)) / (nfft / noverlap));
M = size(x,2); % number of speakers
X = zeros(N_half, L, M);
Lg = numel(g)/M;

for m = 0:M-1
    
    G = fft(g(m*Lg+1:(m+1)*Lg),nfft);

    for l = 0:L-1 % Frame index

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Section of code to complete (3 - 5 lines) %
        x_segment = x((l*nfft/noverlap)+1:(l*nfft/noverlap)+nfft,m+1).*window;
        X_fft = fft(x_segment,nfft);
        % Add filtering with G
        X_filtered = X_fft.*G;

        X(:,l+1,m+1)=X_filtered(1:N_half);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end

end
