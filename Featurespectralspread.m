% ======================================================================
%> @brief computes the spectral spread from the magnitude spectrum
%> called by ::ComputeFeature
%>
%> @param X: spectrogram (dimension FFTLength X Observations)
%> @param f_s: sample rate of audio data 
%>
%> @retval v spectral spread (in Hz)
% ======================================================================
function [vss] = Featurespectralspread (X, f_s)
 
    % get spectral centroid as index
    vsc     = centroid (X, f_s)*2/f_s * size(X,1);
 
    % allocate memory
    vss     = zeros(size(vsc));
 
    % compute spread
    X       = X.^2;
    for (n = 1:size(X,2))
        vss(n)  = (([0:size(X,1)-1]-vsc(n)).^2*X(:,n))./sum(X(:,n));
    end
    vss     = sqrt(vss);
 
    % convert from index to Hz
    vss1     = vss / size(X,1) * f_s/2;
    vrms =  Featuretimerootmeansquare(X);
    vss2=vss1.*vrms;
   
   vss=sum(vss2)./sum(vrms);
end
% ======================================================================
function [vrms] = Featuretimerootmeansquare(x)
  
        % calculate the rms
        vrms     = sqrt(mean(x.^2));
    
end
function [vsc] = centroid (X, f_s)
 
    X       = X.^2;
    vsc     = ([0:size(X,1)-1]*X)./sum(X,1);
 
    % avoid NaN for silence frames
    vsc (sum(X,1) == 0) = 0;
 
    % convert from index to Hz
    vsc     = vsc / size(X,1) * f_s/2;
end