% ======================================================================
%> @brief computes the spectral decrease from the magnitude spectrum
%> called by ::ComputeFeature
%>
%> @param X: spectrogram (dimension FFTLength X Observations)
%> @param f_s: sample rate of audio data (unused)
%>
%> @retval vsk spectral decrease
% ======================================================================
function [vsd] = Featurespectraldecrease (X, f_s)

    % compute index vector
    k       = [0:size(X,1)-1];
    k(1)    = 1;
    kinv    = 1./k;
    
    % compute slope
    vsd1     = (kinv*(X-repmat(X(1,:),size(X,1),1)))./sum(X(2:end,:),1);
    vrms =  Featuretimerootmeansquare(X);
    vsd2=vsd1.*vrms;
   
   vsd=sum(vsd2)./sum(vrms);
end
% ======================================================================
function [vrms] = Featuretimerootmeansquare(x)
  
        % calculate the rms
        vrms     = sqrt(mean(x.^2));
    
end
