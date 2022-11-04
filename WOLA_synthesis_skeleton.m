% Lab 2 for Digital Audio Signal Processing Lab Sessions
% Session 2: Binaural synthesis and 3D audio: OLA and WOLA frameworks
% R.Ali, G. Bernardi, J.Schott, A. Bertrand
% 2021
%
% The following is the skeleton code for the synthesis stage of the WOLA method, which you need to
% complete


function x = WOLA_synthesis_skeleton(X,window,nfft,noverlap)
%WOLA_synthesis inverse short-time fourier transform.
%
% INPUT:
%   X           : input matrix (bins x frames x channels)
%   window      : window function
%   nfft        : FFT size
%   noverlap    : frame overlap; default: 2 (50%)
%
% OUTPUT:
%   x           : output time signal(s)


L = size(X,2); % number of frames
M = size(X,3);
x_synth = zeros(nfft,L,M);
% ## Perform IFFT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section of code to complete (1 - 3 lines) %
for m=1:M
    x_synth(:,:,m) = ifft([zeros(1,size(X,2));X(2:end-1,:,m);zeros(1,size(X,2));conj(flipud(X(2:end-1,:,m)))],nfft);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% ## Apply synthesis window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section of code to complete (1 - 3 lines) %
for m=1:M
    x_synth(:,:,m) = x_synth(:,:,m).*repmat(window,1,size(x_synth,2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% ## Obtain re-synthesised signals
x_length = (nfft/noverlap)*(L-1)+nfft;
x = zeros(x_length,M);

for m = 1:M

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Section of code to complete (1 - 3 lines) %
    for l=0:L-1
        x((l*nfft/noverlap)+1:(l*nfft/noverlap)+nfft,m)=...
           x((l*nfft/noverlap)+1:(l*nfft/noverlap)+nfft,m)+x_synth(:,l+1,m);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

end
