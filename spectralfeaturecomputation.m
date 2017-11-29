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