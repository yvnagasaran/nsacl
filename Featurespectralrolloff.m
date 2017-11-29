% =======================================================================
%> @brief computes the spectral rolloff from the magnitude spectrum
%> called by ::ComputeFeature
%>
%> @param X: spectrogram (dimension FFTLength X Observations)
%> @param f_s: sample rate of audio data 
%>
%> @retval v spectral rolloff (in Hz)
% ======================================================================
function [vsr] = Featurespectralrolloff (X, f_s, kappa)
 
    % initialize parameters
    if (nargin < 3)
        kappa   = 0.85;
    end
 
    % allocate memory
    vsr     = zeros(1,size(X,2));
 
    %compute rolloff
    afSum   = sum(X,1);
    for (n = 1:length(vsr))
        vsr(n)  = find(cumsum(X(:,n)) >= kappa*afSum(n), 1); 
    end
 
    % convert from index to Hz
    vsr1     = vsr / size(X,1) * f_s/2;
    vrms =  Featuretimerootmeansquare(X);
    vsr2=vsr1.*vrms;
    vsr=sum(vsr2)./sum(vrms);
end
% ======================================================================
function [vrms] = Featuretimerootmeansquare(x)
  
        % calculate the rms
        vrms     = sqrt(mean(x.^2));
    
end