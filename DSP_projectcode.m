
function [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11] = spectralfeaturecomputation (afAudioData, f_s, afWindow, iBlockLength, iHopLength)
 clc
 
        % set default parameters if necessary
    if (nargin < 5)
        iHopLength      = 160;
    end
    if (nargin < 4)
        iBlockLength    = 320;
    end
     % pre-processing: down-mixing
    if (size(afAudioData,2)> 1)
        afAudioData = mean(afAudioData,2);
    end
    % pre-processing: normalization (not necessary for many features)
    if (length(afAudioData)> 1)
        afAudioData = afAudioData/max(abs(afAudioData));
          if (nargin < 3 || isempty(afWindow))
            afWindow    = hann(iBlockLeth,'periodic');
        end
 
        % compute FFT window function
        if (length(afWindow) ~= iBlockLength)
            error('window length mismatch');
        end        
         % in the real world, we would do this block by block...
        [x,f,t]     = spectrogram(  afAudioData,...
                                    afWindow,...
                                    iBlockLength-iHopLength,...
                                    iBlockLength,...
                                    f_s);
 
        % magnitude spectrum
        x           = abs(x)*2/iBlockLength;
 
        % compute feature
v1           = Featurespectraltonalpowerratio(x, f_s);
v2           = Featurespectralspread(x, f_s);
v3           = Featurespectralslope(x, f_s);
v4           = Featurespectralskewness(x, f_s);
v5           = Featurespectralrolloff(x, f_s);
v6           = Featurespectralkurtosis(x, f_s);
v7           = Featurespectralflux(x, f_s);
v8           = Featurespectralflatness(x, f_s);        
v9           = Featurespectraldecrease(x, f_s);     
v10          = Featurespectralcrest(x, f_s);    
v11          = Featurespectralcentroid(x, f_s);

    end
% ======================================================================
%> computes the spectral centroid from the (squared) magnitude spectrum
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
    
end% ======================================================================
%> computes the spectral skewness from the magnitude spectrum
% ======================================================================
function [vssk] = Featurespectralskewness (X, f_s)
 
    UseBookDefinition = true;
 
    if (UseBookDefinition)
        % compute mean and standard deviation
        mu_x    = mean(abs(X), 1);
        std_x   = std(abs(X), 1);
 
        % compute skewness
        X       = X - repmat(mu_x, size(X,1), 1);
        vssk1    = sum ((X.^3)./(repmat(std_x, size(X,1), 1).^3*size(X,1)));
    else
        % interpret the spectrum as pdf, not as signal
        f       = linspace(0, f_s/2, size(X,1));
        % compute mean and standard deviation
        mu_X    = (f * X) ./ sum(X,1);
        tmp     = repmat(f, size(X,2),1) - repmat(mu_X, size(X,1),1)';
        var_X   = diag (tmp.^2 * X) ./ (sum(X,1)'*size(X,1));
 
        vssk1    = diag (tmp.^3 * X) ./ (var_X.^(3/2) .* sum(X,1)'*size(X,1));
    end
    vrms =  Featuretimerootmeansquare(X);
    vssk2=vssk1.*vrms;
    vssk=sum(vssk2)./sum(vrms);
end
% ======================================================================
function [vrms] = Featuretimerootmeansquare(x)
  
        % calculate the rms
        vrms     = sqrt(mean(x.^2));
    
end
% ======================================================================
%> computes the tonal power ratio from the magnitude spectrum
% ======================================================================
function [vtpr] = Featurespectraltonalpowerratio(X, f_s, G_T)
 
    % initiliaze
    if (nargin < 3)
        G_T = 5e-4;
    end
 
    % allocate memory
    vtpr1    = zeros(1,size(X,2));
 
    X       = X.^2;
    fSum    = sum(X,1);
 
    for (n = 1:size(X,2))
        if (fSum == 0)
            % do nothing for 0-blocks
            continue;
        end
        % find local maxima
        [afPeaks]   = findpeaks(X(:,n));
 
        % find peaks above the threshold
        k_peak      = find(afPeaks > G_T);
 
        % calculate the ratio
        vtpr1(n)     = sum(afPeaks(k_peak))/fSum(n);
    end
    vrms =  Featuretimerootmeansquare(X);
    vtpr2=vtpr1.*vrms;
   
   vtpr=sum(vtpr2)./sum(vrms);
end
% ======================================================================
function [vrms] = Featuretimerootmeansquare(x)
  
        % calculate the rms
        vrms     = sqrt(mean(x.^2));
    
end
% ======================================================================
%>  computes the spectral spread from the magnitude spectrum
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
% ======================================================================
%> computes the spectral slope from the magnitude spectrum
% ======================================================================
function [vssl] = Featurespectralslope (X, f_s)
 
    % compute mean
    mu_x    = mean(abs(X), 1);
 
    % compute index vector
    kmu     = [0:size(X,1)-1] - size(X,1)/2;
 
    % compute slope
    X       = X - repmat(mu_x, size(X,1), 1);
    vssl1    = (kmu*X)/(kmu*kmu');
    vrms =  Featuretimerootmeansquare(X);
    vssl2=vssl1.*vrms;
    vssl=sum(vssl2)./sum(vrms);
end
% ======================================================================
function [vrms] = Featuretimerootmeansquare(x)
  
        % calculate the rms
        vrms     = sqrt(mean(x.^2));
    
end
% =======================================================================
%>  computes the spectral rolloff from the magnitude spectrum
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
% =====================================================================
% computes the spectral kurtosis from the magnitude spectrum
% ======================================================================
function [vsk] = Featurespectralkurtosis (X, f_s)
 
 
    UseBookDefinition = true;
 
    if (UseBookDefinition)
        % compute mean and standard deviation
        mu_x    = mean(abs(X), 1);
        std_x   = std(abs(X), 1);
 
        % remove mean
        X       = X - repmat(mu_x, size(X,1), 1);
 
        % compute kurtosis
        vsk1     = sum ((X.^4)./(repmat(std_x, size(X,1), 1).^4*size(X,1)));
    else
        % interpret the spectrum as pdf, not as signal
        f       = linspace(0, f_s/2, size(X,1));
        % compute mean and standard deviation
        mu_X    = (f * X) ./ sum(X,1);
        tmp     = repmat(f, size(X,2),1) - repmat(mu_X, size(X,1),1)';
        var_X   = diag (tmp.^2 * X) ./ (sum(X,1)'*size(X,1));
 
        vsk1    = diag (tmp.^4 * X) ./ (var_X.^2 .* sum(X,1)'*size(X,1));
    end
    vsk1     = vsk1-3;
    vrms =  Featuretimerootmeansquare(X);
    vsk2=vsk1.*vrms;
    vsk=sum(vsk2)./sum(vrms);
end
% ======================================================================
function [vrms] = Featuretimerootmeansquare(x)
  
        % calculate the rms
        vrms     = sqrt(mean(x.^2));
    
end
% ======================================================================
%computes the spectral flux from the magnitude spectrum

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
% ======================================================================
 %computes the spectral flatness from the magnitude spectrum

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
% ======================================================================
 %computes the spectral decrease from the magnitude spectrum
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
% ======================================================================
 %computes the spectralcrest from the magnitude spectrum
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
% ======================================================================
%K-Means Clustering
% ======================================================================
[data,fs] = wavread('audiofile.wav');

data = data/abs(max(data));

%framing for silence removal
f_d = 0.025; f_size = floor(f_d*fs); n = length(data); n_f = floor(n/f_size);

temp = 0;

for i = 1 : n_f
 frames(i,:)= data(temp+1: temp + f_size);
 temp = temp + f_size;
end


%silence removal based on max amplitude
m_amp = abs(max(frames,[],2));
id = find(m_amp > 0.5); % finding ID of frames with max amp > 0.05
fr_ws = frames(id,:); % frames without silence

% reconstruct signal
y = reshape(fr_ws',1,[]);

data = data/abs(max(data));
b_l = 1000;
h_l = 500;
win    = hann(b_l,'periodic');
[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11] = spectralfeaturecomputation (data,fs,win,b_l, h_l);
s = [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11]
%labeling of dataset
H= ['://dataset of 80*11 matrix/'];
[p,c] = kmeans(H,k);
end

