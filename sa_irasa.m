%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sa_irasa
% by Til Ole Bergmann 2018
% last modified 2018/01/05 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Separating oscillatory and fractal components of power spectra using the 
% method of Irregular-Resampling Auto-Spectral Analysis (IRASA)¹, based on 
% code available at https://purr.purdue.edu/publications/1987/1. This wrapper
% works on data structures in FieldTrip format and additionally calculates 
% the percent change of oscillatory power relative to the fractal component 
% to effectively normalize across frequency bins per channel and subject.
%
% Reference:
% ¹ Wen, H., and Liu, Z. (2016). Separating Fractal and Oscillatory Components 
% in the Power Spectrum of Neurophysiological Signal. Brain Topogr 29, 13-26.
%
% Use as:
%
% FUNCTION
% freq = sa_irasa(cfg,data)
%
% INPUTS
% data = data structure in FieldTrip format
% cfg = configuration structure defing parameters
% cfg.channel = Nx1 cell array with selection of channels (default = 'all'), see ft_channelselection for detaiils
% cfg.keeptrials = whether or not to keep single trial data, 'yes' = keep single trials, 'no' = average across trials (default = 'no')
% cfg.foi = output frequency range as [fmin fmax] (default [<lowest possible freq> : <smallest possible step size> : samplingrate/5])
% cfg.detrend = detrending data before FFT, 1 = yes, 0 = no (default 1) 
% cfg.filter = filtering before downsampling to avoid aliasing, 1 = yes, 0 = no (default 1)
% cfg.hset = array containing scaling factors >1, (default 1.1:0.05:1.9)
% cfg.mainoutput = IRASA output field, the content of freq.powspctrm is replaced with: 'mixd','frac','osci','perc' (default is 'perc')
% cfg.resamplefreqs = whether ('yes') or not ('no') to resample odd frequency bins resulting from nextpow2 to the requested foi (default is 'yes')
% 
% OUTPUTS
% freq = data structure for frequency domain data in FieldTrip format 
% freq.mixd = uncorrected powerspectrum 
% freq.frac = fractal part of the power spectrum (1/f part)
% freq.osci = 1/f corrected oscillatory part
% freq.perc = percent change of the oscillatory relative to the fractal part 
% freq.Beta = power-law exponent
% freq.Cons = power intersect of power-law line in log-log scale
% freq.Plaw = power-law spectrum
% freq.powspctrm = one of the following fields: 'mixd','frac','osci','perc'  
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function freq = sa_irasa(gcfg,data)
tic;

%% check and adjust config
if ~isfield(gcfg,'channel'), gcfg.channel = {'all'}; end
if ~isfield(gcfg,'keeptrials'), gcfg.keeptrials = 'no'; end
if ~isfield(gcfg,'foi'), gcfg.foi = [1/(size(data.trial{1},2)/data.fsample):1/(size(data.trial{1},2)/data.fsample):floor(data.fsample/5)]; end
if ~isfield(gcfg,'detrend'), gcfg.detrend = 1; end
if ~isfield(gcfg,'filter'), gcfg.filter = 1; end
if ~isfield(gcfg,'hset'), gcfg.hset = [1.1:0.05:1.9]; end
if ~isfield(gcfg,'mainoutput'), gcfg.mainoutput = 'perc'; end
if ~isfield(gcfg,'resamplefreqs'), gcfg.resamplefreqs = 'yes'; end

%% preparation
display('Runing IRASA...');

%% select channels of interest
if ~strcmp(gcfg.channel,'all')
    cfg = [];
    cfg.channel = gcfg.channel;
    data = ft_selectdata(cfg,data);
end

%% set up freq structure
freq = struct; 
switch gcfg.keeptrials
    case 'yes'
        freq.dimord = 'chan_freq_rep';        
    case 'no'
        freq.dimord = 'chan_freq';        
end
freq.label = data.label;
if isfield(data,'cfg')
    gcfg.previous = data.cfg;
end
freq.cfg = gcfg;


%% call IRASA

% rearrange data to chan_rep_time dimorder to allow IRASA computation for all trials of one channel at once
dat = permute(reshape([data.trial{:}],[numel(data.label) size(data.trial{1},2) size(data.trial,2)]),[1,3,2]);
frange = [0 data.fsample/4];

switch gcfg.keeptrials
    case 'yes' % keeptrials
        display('Single trials are kept, therefore percent change and resampling are based on single-trial data.');
    case 'no' % keeptrials
        display('Single trials are NOT kept, therefore percent change and resampling are based on averaged data.');
end

for iChan = 1:size(dat,1)    
    
    display(['Computing channel ' data.label{iChan} ', number ' num2str(iChan) ' of ' num2str(length(data.label)) ' channels.']);              
    spec = amri_sig_fractal(squeeze(dat(iChan,:,:))',data.fsample,'frange',frange,'detrend',gcfg.detrend,'filter',gcfg.filter,'hset',gcfg.hset);          
    
    switch gcfg.keeptrials        
        case 'yes' % keeptrials            
            spec = amri_sig_plawfit(spec,frange);  
            spec.perc = spec.osci .* (ones(size(spec.osci)) ./ spec.frac); % percent change of osci from frac            
            switch gcfg.resamplefreqs
                case 'no' % resamplefreqs
                    freq.freq = spec.freq';
                    for iTrial = 1:size(dat,2)
                        freq.mixd(iChan,:,:) = spec.mixd;
                        freq.frac(iChan,:,:) = spec.frac;
                        freq.osci(iChan,:,:) = spec.osci;
                        freq.perc(iChan,:,:) = spec.perc;
                        freq.Plaw(iChan,:,:) = spec.Plaw;
                    end
                case 'yes' % resamplefreqs
                    freq.freq = gcfg.foi;
                    for iTrial = 1:size(dat,2)
                        freq.mixd(iChan,:,iTrial) = interpn(spec.freq,spec.mixd(:,iTrial),gcfg.foi,'pchip')';
                        freq.frac(iChan,:,iTrial) = interpn(spec.freq,spec.frac(:,iTrial),gcfg.foi,'pchip')';
                        freq.osci(iChan,:,iTrial) = interpn(spec.freq,spec.osci(:,iTrial),gcfg.foi,'pchip')';
                        freq.perc(iChan,:,iTrial) = interpn(spec.freq,spec.perc(:,iTrial),gcfg.foi,'pchip')';
                        freq.Plaw(iChan,:,iTrial) = interpn(spec.freq,spec.Plaw(:,iTrial),gcfg.foi,'pchip')';
                    end                    
            end
            clear spec;
            
        case 'no' % keeptrials
            switch gcfg.resamplefreqs                    
                case 'no' % resamplefreqs
                    freq.freq = spec.freq';
                    freq.mixd(iChan,:) = squeeze(mean(spec.mixd,2));
                    freq.frac(iChan,:) = squeeze(mean(spec.frac,2));
                    freq.osci(iChan,:) = squeeze(mean(spec.osci,2));
                    freq.perc(iChan,:) = squeeze(mean(spec.osci,2)) .* (ones(size(squeeze(mean(spec.osci,2)))) ./ squeeze(mean(spec.frac,2))); % percent change of osci from frac
                    temp = struct;
                    temp.freq = spec.freq;
                    temp.mixd = squeeze(mean(spec.mixd,2));
                    temp.frac = squeeze(mean(spec.frac,2));
                    temp.osci = squeeze(mean(spec.osci,2));
                    temp.perc = squeeze(mean(spec.osci,2)) .* (ones(size(squeeze(mean(spec.osci,2)))) ./ squeeze(mean(spec.frac,2))); % percent change of osci from frac
                    temp = amri_sig_plawfit(temp,frange);
                    freq.Beta(iChan,:) = temp.Beta;
                    freq.Cons(iChan,:) = temp.Cons;
                    freq.Plaw(iChan,:) = temp.Plaw;
                    clear temp;                    
                    
                case 'yes' % resamplefreqs
                    freq.freq = gcfg.foi;                    
                    freq.mixd(iChan,:) = interpn(spec.freq,squeeze(mean(spec.mixd,2)),gcfg.foi,'pchip')';
                    freq.frac(iChan,:) = interpn(spec.freq,squeeze(mean(spec.frac,2)),gcfg.foi,'pchip')';
                    freq.osci(iChan,:) = interpn(spec.freq,squeeze(mean(spec.osci,2)),gcfg.foi,'pchip')';                                        
                    freq.perc(iChan,:) = interpn(spec.freq,(squeeze(mean(spec.osci,2)) .* (ones(size(squeeze(mean(spec.osci,2)))) ./ squeeze(mean(spec.frac,2)))),gcfg.foi,'pchip')'; % percent change of osci from frac                    
                    temp = struct;
                    temp.freq = spec.freq;
                    temp.mixd = squeeze(mean(spec.mixd,2));
                    temp.frac = squeeze(mean(spec.frac,2));
                    temp.osci = squeeze(mean(spec.osci,2));
                    temp.perc = squeeze(mean(spec.osci,2)) .* (ones(size(squeeze(mean(spec.osci,2)))) ./ squeeze(mean(spec.frac,2))); % percent change of osci from frac
                    temp = amri_sig_plawfit(temp,frange);
                    freq.Beta(iChan,:) = temp.Beta;
                    freq.Cons(iChan,:) = temp.Cons;                    
                    freq.Plaw(iChan,:) = interpn(spec.freq,temp.Plaw,gcfg.foi,'pchip')'; 
                    clear temp;                                 
            end
            clear spec;
    end
    
    
    
end % of loop over channels
clear dat;

%% add powspctrm field to freq structure to maintain FieldTrip format
switch gcfg.mainoutput
    case 'mixd'
        freq.powspctrm = freq.mixd;
    case 'frac'
        freq.powspctrm = freq.frac;
    case 'osci'
        freq.powspctrm = freq.osci;
    case 'perc'
        freq.powspctrm = freq.perc;
end


toc
end % of main function sa_irasa()



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunction amri_sig_fractal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Separate the spectra of fractal component and oscillatory component from mixed time series
%   amri_sig_fractal()
%
% Usage
%   spec = amri_sig_fractal(sig,srate,...)
%
% Inputs
%   sig   - a time-series vector. If sig is a matrix, then separate spectra for each column
%   srate - sampling rate
%
% Outputs
%   spec  - spectrum
%           .freq = a vector of frequency points
%           .srate= sample rate;
%           .mixd = spectrum of mixed time series
%           .frac = spectrum of fractal component
%           .osci = spectrum of oscillatory component
%
% Keywords
%   frange  - [fmin fmax](default [0 srate/4]), the output frequency range.
%   detrend - 1 or 0 (default 1): 1 means detrending data before fft, otherwise 0
%   filter  - 1 or 0 (default 1): 1 means filtering before downsampling to avoid aliasing.
%   hset    - (default 1.1:0.05:1.9) an array containing scaling factors (> 1).
%
% See also
%   amri_sig_genfrac
%   amri_sig_plawfit
%
% Version
%   0.10
%
% Reference
%   -Wen H. and Liu Z. (2015) Separating Fractal and Oscillatory Components in
%    the Power Spectrum of Neurophysiological Signals

%% History
% 0.01 - HGWEN - 12/15/2013 - Use resample function instead of interp1
%                           - upsample signal before extracting fractals
% 0.02 - HGWEN - 01/15/2014 - add a new method 'nivcgsa'
%                           - fit power-law line after resampling data in euqal space
%                           - set nfft=2^nextpow2(2*Ndata);
% 0.03 - HGWEN - 02/27/2014 - Change the name 'nivcgsa' to "IRASA".
% 0.04 - HGWEN - 03/01/2014 - In IRASA, use median instead of min operator in the final step
% 0.05 - HGWEN - 05/14/2014 - If sig is a matrix, separate spectra for each column
% 0.06 - HGWEN - 08/20/2014 - Use multiple h values in CGSA
%                           - remove the power-law fitting section, and add a new function "amri_sig_plawfit"
%                           - Only return freq, srate, mixd, and frac.
% 0.07 - HGWEN - 10/11/2014 - Add a keyword "upsample"
% 0.08 - HGWEN - 10/19/2014 - Reorganized the code
% 0.09 - HCWEN - 04/11/2015 - Added keywords 'hset' and 'filter', and removed the keyword 'upsample'
% 0.10 - HGWEN - 09/26/2015 - Reorganized the structure.

%%

    function spec = amri_sig_fractal(sig,srate,varargin)

        if nargin<2
            eval('help amri_sig_fractal');
            return
        end
        
        %% defaults
        flag_detrend = 1;
        fmin = 0;
        fmax = srate/4;
        flag_filter = 1;
        hset=1.1:0.05:1.9;
        
        %% Keywords
        for i = 1:2:size(varargin,2)
            Keyword = varargin{i};
            Value   = varargin{i+1};
            if strcmpi(Keyword,'frange')
                fmin = max(Value(1),fmin);
                fmax = min(Value(2),fmax);
            elseif strcmpi(Keyword,'detrend')
                flag_detrend = Value;
            elseif strcmpi(Keyword,'filter')
                flag_filter = Value;
            elseif strcmpi(Keyword,'hset')
                hset = Value;
            else
                warning(['amri_sig_fractal(): unknown keyword ' Keyword]);
            end
        end
        
        %% preprocessing
        sig = double(sig);
        if isvector(sig)
            sig = sig(:);
        end
        
        % detrend signal
        if flag_detrend >= 1
            sig = detrend(sig,'linear');
        end
        
        %% apply IRASA method to separate fractal and oscillatory components
        [Smixd, Sfrac, freq] = irasa(sig,srate,hset,flag_filter);
        
        %% only keep the given frequency range
        ff = (freq>=fmin & freq<=fmax & freq>0);
        freq = freq(ff);
        Smixd = Smixd(ff,:);
        Sfrac = Sfrac(ff,:);
        
        %% outputs
        spec.freq  = freq;
        spec.srate = srate;
        spec.mixd  = Smixd;
        spec.frac  = Sfrac;
        spec.osci  = Smixd - Sfrac;
        
    end

%% IRASA Irregular-Resampling Auto-Spectral Analysis

    function [Smixd, Sfrac, freq] = irasa(sig,srate,hset,flag_filter)
        % Given a discrete time series (sig) of length (Ntotal)
        Ntotal = size(sig,1);
        dim = size(sig,2);
        
        % Ndata is the power of 2 that does not exceed 90% of Ntotal.
        Ndata = 2^floor(log2(Ntotal*0.9));
        
        % Nsubset is fixed to 15
        Nsubset = 15;
        
        % compute the auto-power spectrum of the originally sampled time series
        L = floor((Ntotal-Ndata)/(Nsubset-1));
        
        % set nfft greater than ceil(hset(end))*Ndata, asure that do fft without truncating
        nfft = 2^nextpow2(ceil(hset(end))*Ndata);
        
        % set output data length Nfrac
        Nfrac = nfft/2 + 1;
        freq = srate/2*linspace(0,1,Nfrac); freq = freq(:);
        
        % compute the spectrum of mixed data
        Smixd = zeros(Nfrac,dim);
        taper = gettaper([Ndata dim]);
        for k = 0:Nsubset-1
            i0 = L*k+1;
            x1 = sig(i0:1:i0+Ndata-1,:);
            p1 = fft(x1.*taper,nfft)/min(nfft,size(x1,1));
            p1(2:end,:) = p1(2:end,:)*2;
            Smixd = Smixd+abs(p1(1:Nfrac,:)).^2;
        end
        Smixd = Smixd/Nsubset;
        
        % filter the input signal to avoid alising when downsampling
        if flag_filter == 1
            sig_filtered = sig;
            for i = 1 : size(sig,2)
                sig_filtered(:,i) = amri_sig_filtfft(sig(:,i),srate,0,srate/(2*ceil(hset(end))));
            end
        end
        
        % compute fractal component.
        Sfrac = zeros(Nfrac,dim,length(hset));
        for ih = 1:length(hset)
            % compute the auto-power spectrum of xh
            h = hset(ih);
            [n, d] = rat(h); % n > d
            Sh = zeros(Nfrac,dim);
            for k = 0 : Nsubset-1
                i0 = L*k + 1;
                x1 = sig(i0:i0+Ndata-1,:);
                xh = myresample(x1, n, d);
                taperh = gettaper(size(xh));
                ph = fft(xh.*taperh,nfft)/min(nfft,size(xh,1));
                ph(2:end,:) = ph(2:end,:)*2;
                tmp = (abs(ph)).^2;
                Sh = Sh + tmp(1:Nfrac,:);
            end
            Sh = Sh / Nsubset;
            
            % compute the auto-power spectrum of X1h
            S1h = zeros(Nfrac, dim);
            for k = 0 : Nsubset - 1
                i0 = L*k + 1;
                if (flag_filter==1)
                    x1 = sig_filtered(i0:1:i0+Ndata-1,:);
                else
                    x1 = sig(i0:1:i0+Ndata-1,:);
                end
                x1h = myresample(x1,d,n);
                taper1h = gettaper(size(x1h));
                p1h = fft(x1h.*taper1h,nfft)/min(nfft,size(x1h,1));
                p1h(2:end,:) = p1h(2:end,:)*2;
                tmp = (abs(p1h)).^2;
                S1h = S1h + tmp(1:Nfrac,:);
            end
            S1h = S1h / Nsubset;
            Sfrac(:,:,ih)= sqrt(Sh.*S1h);
        end
        
        % taking median
        Sfrac = median(Sfrac,3);
    end


%% subfunctions
    function taper = gettaper(S)
        % get a tapering function for power spectrum density calculation
        if license('test','signal_toolbox')
            taper = hanning(S(1),'periodic');
        else
            taper = 0.5*(1-cos(2*pi*(1:S(1))/(S(1)-1)));
        end
        taper = taper(:);
        taper = repmat(taper,1,S(2));
    end

    function data_out = myresample(data,L,D)
        % resample signal with upsample L and downsample D
        if license('test','signal_toolbox')
            data_out = resample(data,L,D);
        else
            N = size(data,1);
            x0 = linspace(0,1,N);
            x1 = linspace(0,1,round(N*L/D));
            data_out = interp1(x0,data,x1);
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunction amri_sig_plawfit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% fitting power-law function to scale-free power-spectrum
%   amri_sig_plawfit()
%
% Usage
%   spec = amri_sig_plawfit(spec, frange)
%
% Inputs
%   spec   - spec.freq: frequency points
%          - spec.frac: the scale-free power spectrum.
%   frange - given frequency range for power-law fitting
%
% Outputs
%   spec  - spectrum
%           .Freq = a vector of frequency points with given frequency range
%           .Plaw = power-law spectrum
%           .Beta = the power-law exponent within the given frequency range
%           .Cons = the power intersect of power-law line in log-log scale.
%
% Version
%   0.01
%
% Reference
%   -Wen H. and Liu Z. (2015) Separating Fractal and Oscillatory Components in
%    the Power Spectrum of Neurophysiological Signals
%
%% History
% 0.01 - HGWEN - 12/20/2013 - original file
%

%%
    function spec = amri_sig_plawfit(spec, frange)

        % define frequency range
        ff = spec.freq >= frange(1) & spec.freq <= frange(2);
        freq = spec.freq(ff);
        frac = spec.frac(ff,:,:);
        
        % convert to log-log scale
        logfreq = log10(freq);
        y1 = log10(frac);
        
        % resample frac in equal space
        x2 = linspace(min(logfreq),max(logfreq),length(logfreq)); x2 = x2(:);
        y2 = interp1(logfreq,y1,x2);
        
        % fitting power-law function
        Nt = size(y2,2);
        Nc = size(y2,3);
        beta = zeros(Nt,Nc);
        cons = zeros(Nt,Nc);
        plaw = zeros(size(frac));
        
        for j = 1 : Nc
            for i = 1 : Nt
                % ordinary least square
                p = polyfit(x2,y2(:,i,j),1);
                
                beta(i,j) = -p(1);
                cons(i,j) = p(2);
                powlaw = 10.^(polyval(p,logfreq));
                plaw(:,i,j) = powlaw(:);
            end
        end
        
        % outputs
        spec.Beta = beta;
        spec.Cons = cons;
        spec.Plaw = plaw;
        spec.Freq = freq;
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

