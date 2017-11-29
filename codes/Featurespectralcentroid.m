% ======================================================================
%> @brief computes the spectral centroid from the (squared) magnitude spectrum
%> called by ::ComputeFeature
%>
%> @param X: spectrogram (dimension FFTLength X Observations)
%> @param f_s: sample rate of audio data 
%>
%> @retval v spectral centroid (in Hz)
% ======================================================================
function [vsc] = Featurespectralcentroid (X, f_s)
 
    X       = X.^2;
    vsc     = ([0:size(X,1)-1]*X)./sum(X,1);
 
    % avoid NaN for silence frames
    vsc (sum(X,1) == 0) = 0;
 
    % convert from index to Hz
    vsc1     = vsc / size(X,1) * f_s/2;
    vrms =  Featuretimerootmeansquare(X);
    vsc2=vsc1.*vrms;
   
   vsc=sum(vsc2)./sum(vrms);
end
% ======================================================================
function [vrms] = Featuretimerootmeansquare(x)
  
        % calculate the rms
        vrms     = sqrt(mean(x.^2));
    
end