%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% freq_phase_scramble
% by Til Ole Bergmann 2013 based on code from Bernhard Staresina
% last modified 2016/11/22 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% phase scramble time-frequency data 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = sa_freq_phase_scramble(cfg, data)

segmentlength = numel(data.time{1});
window        = ones(segmentlength,1); % hann(segmentlength);

for itrl = 1:numel(data.trial)    
    
    for ichan = 1:numel(data.label)      
        
        trldat  = data.trial{itrl}(ichan,:);
        dummy   = window.*trldat';
        fourier = fft(dummy);
        fft_ang = angle(fourier);
        fft_abs = abs(fourier);
        
        % add noise to angles
        randval         = rand(size(fourier)).*(2*pi);
        ang_new         = fft_ang+randval;
        helper          = fourier(1);
        fouriernew      = fft_abs.*(cos(ang_new)+1i*sin(ang_new));
        fouriernew(1)   = helper;
        
        % recreat symmetrie of fourier transformed data
        for sample = 2:segmentlength/2
            fouriernew(sample) = real(fouriernew(segmentlength-sample+2)) ...
                - imag(fouriernew(segmentlength-sample+2))*1i;
        end
        
        % transform back
        trldatnew = ifft(fouriernew);
        
        % replace original trial values
        data.trial{itrl}(ichan,:) = real(trldatnew);
        
    end
end

end % of function