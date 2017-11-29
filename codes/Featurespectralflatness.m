% ======================================================================
%> @brief computes the spectral flatness from the magnitude spectrum
%> called by ::ComputeFeature
%>
%> @param X: spectrogram (dimension FFTLength X Observations)
%> @param f_s: sample rate of audio data (unused)
%>
%> @retval vtf spectral flatness
% ======================================================================
function [vtf] = Featurespectralflatness (X, f_s)
 
    XLog    = log(X+1e-20);
    vtf1     = exp(mean(XLog,1)) ./ mean(X,1);
    vrms =  Featuretimerootmeansquare(X);
    vtf2=vtf1.*vrms;
   
   vtf=sum(vtf2)./sum(vrms);
end
% ======================================================================
function [vrms] = Featuretimerootmeansquare(x)
  
        % calculate the rms
        vrms     = sqrt(mean(x.^2));
    
end