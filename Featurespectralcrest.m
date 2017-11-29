% ======================================================================
%> @brief computes the spectral crest from the magnitude spectrum
%> called by ::ComputeFeature
%>
%> @param X: spectrogram (dimension FFTLength X Observations)
%> @param f_s: sample rate of audio data (unused)
%>
%> @retval vtsc spectral crest factor
% ======================================================================
function [vtsc] = Featurespectralcrest (X, f_s)
 
   vtsc1 = max(X,[],1) ./ sum(X,1);
    vrms =  Featuretimerootmeansquare(X);
    vtsc2=vtsc1.*vrms;
   
   vtsc=sum(vtsc2)./sum(vrms);
end
% ======================================================================
function [vrms] = Featuretimerootmeansquare(x)
  
        % calculate the rms
        vrms     = sqrt(mean(x.^2));
    
end