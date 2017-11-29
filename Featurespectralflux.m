% ======================================================================
%> @brief computes the spectral flux from the magnitude spectrum
%> called by ::ComputeFeature
%>
%> @param X: spectrogram (dimension FFTLength X Observations)
%> @param f_s: sample rate of audio data (unused)
%>
%> @retval v spectral flux
% ======================================================================
function [vsf] = Featurespectralflux (X, f_s)
 
    % difference spectrum (set first diff to zero)
    afDeltaX    = diff([X(:,1), X],1,2);
 
    % flux
    vsf1         = sqrt(sum(afDeltaX.^2))/size(X,1);
    vrms =  Featuretimerootmeansquare(X);
    vsf2=vsf1.*vrms;
   
   vsf=sum(vsf2)./sum(vrms);
end
% ======================================================================
function [vrms] = Featuretimerootmeansquare(x)
  
        % calculate the rms
        vrms     = sqrt(mean(x.^2));
    
end