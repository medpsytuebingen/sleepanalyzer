%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sa_inspectsleep
% by Til Ole Bergmann 2016
% last modified 2018/02/11 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% inspect sleep architecture in EEG data
%
% si = sleep information structure
% si.freq = results from freqanalyses separately for sleep stages
% si.architect. = sleep architecture
% si.figurehandle.hypno = handle to hypnogramme figure
% si.figurehandle.topo = handle to topoplot figure
%
% cfg = [];
% cfg.channel = Nx1 cell array with selection of channels (default = 'all'), see ft_channelselection for detaiils
% cfg.sleepscoringStruct = 'scoring', or 'trialinfo' indicating field containing sleep scoring (defaule = 'scoring')    
% cfg.sleepscoringInfoRow = row number of sleepscoringStruct containing sleep scorings (default = 1)
% cfg.loadartifactfilter = 1 (load artifacts from external mat file), 0 (use artifact field from data file (default)), 2 (neglect all artifacts)
% cfg.artifactfiltername = string indicating the .mat file containing artifact information without channel prefix
% cfg.artifactfreechannels = {'Cz'}; % channels for which artifact files need to be loaded (default = {'Cz'})
% cfg.artifactPaddingPre = time period to ignore before an artifact (in seconds), only when no artifactfilter loaded (default = 1)
% cfg.artifactPaddingPost = time period to ignore after an artifact (in seconds), only when no artifactfilter loaded (default = 1)
% cfg.layout = string indicating layout file, e.g., 'easycapM3.mat';
% cfg.absfreqwin = frequency window to plot for FFTs with absolute scaling, specified as [low high] in Hz (default = [0 45]);
% cfg.relfreqwin = frequency window to plot for FFTs with relative scaling, specified as [low high] in Hz (default = [4 20])
% cfg.maxadjfreqwin = frequency window to determine maximum for scaling of y-axis, specified as [low high] in Hz (default = [8 20]);
% cfg.topofreq = frequency ranges, specified as {[low high], [low high],... } in Hz (default = {[]}, i.e. no topolpots)              
% cfg.fftplotdim = string indicating dimension along which to spearate lots, 'channel' or 'sleepstage' (default = 'channel')
% cfg.normalizepower = string indicating normalization method, 'none', '1overf', 'irasa'* (default = 'none')
% cfg.continuousTFR = 'yes or 'no' calculate continuous TFR (default = 'no')
% cfg.fftwinlength = 4; % in seconds for FFT (default = 4)
% cfg.fftwinoverlap = 0.5; % relative window overlap for FFT, values between 0 and <1 (default = 0.5)
% cfg.foi = vector of frequency bins (default = [0.25:0.25:45]) 
% [si] = sa_inspectsleep(cfg,data);
%
%  * the method 'irasa' outputs, in additon to 'powspctrm' (the 1/f corrected oscillatory part), several additional fields, namely
% 'mixd' (the uncorrected powerspectrum), 'frac' (the fractal, 1/f part), 'rel' the percent change of the oscillatory from the faractal part, 
% 'Beta' (the power-law exponent), 'Cons' (the power intersect of power-law line in log-log scale), 'Plaw (power-law spectrum)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [si] = sa_inspectsleep(gcfg,data)
tic;

%% check and adjust config
if ~isfield(gcfg,'channel'), gcfg.channel = {'all'}; end
if ~isfield(gcfg,'sleepscoringStruct'), gcfg.sleepscoringStruct = 'scoring'; end
if ~isfield(gcfg,'sleepscoringInfoRow'), gcfg.sleepscoringInfoRow = 1; end
if ~isfield(gcfg,'loadartifactfilter'), gcfg.loadartifactfilter = 0; end
if ~isfield(gcfg,'artifactfiltername'), gcfg.artifactfiltername = ''; end
if ~isfield(gcfg,'artifactfreechannels'), gcfg.artifactfreechannels = {'Cz'}; end
if ~isfield(gcfg,'artifactPaddingPre'), gcfg.artifactPaddingPre = 1; end
if ~isfield(gcfg,'artifactPaddingPost'), gcfg.artifactPaddingPost = 1; end
if ~isfield(gcfg,'absfreqwin'), gcfg.absfreqwin = [0 45]; end
if ~isfield(gcfg,'relfreqwin'), gcfg.relfreqwin = [4 20]; end
if ~isfield(gcfg,'maxadjfreqwin'), gcfg.maxadjfreqwin = [8 20]; end
if ~isfield(gcfg,'topofreq'), gcfg.topofreq = {}; end
if ~isfield(gcfg,'fftplotdim'), gcfg.fftplotdim = 'channel'; end
if ~isfield(gcfg,'normalizepower'), gcfg.normalizepower = 'none'; end
if ~isfield(gcfg,'continuousTFR'), gcfg.continuousTFR = 'no'; end
if ~isfield(gcfg,'fftwinlength'), gcfg.fftwinlength = 4; end
if ~isfield(gcfg,'fftwinoverlap'), gcfg.fftwinoverlap = 0.5; end
if ~isfield(gcfg,'foi'), gcfg.foi = [0.25:0.25:45]; end

%% preparation
display('Inspecting sleep architecture...');
si = []; % sleep
si.cfg = gcfg;

%% select channels of interest
cfg = [];
cfg.channel = gcfg.channel;
data = ft_selectdata(cfg,data);

% %% generate scoring array
% switch gcfg.sleepscoringStruct % make sure scoring information is availabvle in data.scoring
%     case 'trialinfo', data.scoring = int8(data.trialinfo(gcfg.sleepscoringInfoRow,:));
%         data.trialinfo(gcfg.sleepscoringInfoRow) = [];
%     case 'scoring'
%         data.scoring = int8(data.scoring);
% end

% generate gcfg.scoring array
gcfg.scoring = int8([]); 

% column 1: sleep stages 
% 0 = wake, 1 = stage N1 (S1), 2 = stage N2 (S2), 3 = stage N3 (S3), 4 = stage S4, 5 = REM, 6 = MT (movement time), 7 = unknown
if isfield(gcfg,'sleepscoringInfoRow'), gcfg.scoring(1,:) = int8(data.scoring(gcfg.sleepscoringInfoRow,:));
else gcfg.scoring(1,:) = ones(1,length(data.scoring))*7; end

% column 2: movement arousals
% 0 = clean, 1 = MA (movement arousal)
if isfield(gcfg,'movementarousalinfo'), gcfg.scoring(2,:) = int8(data.scoring(gcfg.movementarousalinfo,:));
else gcfg.scoring(2,:) = zeros(1,length(data.scoring)); end

% column 3: stages to gauge search criteria
if isfield(gcfg,'gaugestages'), gcfg.scoring(3,:) = int8(ismember(gcfg.scoring(1,:),gcfg.gaugestages) & gcfg.scoring(2,:)==0);
else gcfg.scoring(3,:) = ones(1,length(data.scoring)); end

% column 4: stages to search events
if isfield(gcfg,'searchstages'), gcfg.scoring(4,:) = int8(ismember(gcfg.scoring(1,:),gcfg.searchstages) & gcfg.scoring(2,:)==0);
else gcfg.scoring(4,:) = ones(1,length(data.scoring)); end


%% prepare artifactFilter
if gcfg.loadartifactfilter == 1
    for i = 1:length(gcfg.artifactfreechannels) % over all channels that need to be artifact free
        temp = load(fullfile(gcfg.datapath,['artifactFilter_' [gcfg.artifactfreechannels{i}] gcfg.artifactfiltername]));
        filterArray(i,:) = temp.artifactFilter;
        clear temp;
    end 
    combinedFilter = all(filterArray,1); % combined filter is only true when filters are true for all channels 
    for i = 1:size(data.trial,2) % loop over trials
        data.artifactFilter{i} = [logical(ones(size(data.trial{i})))]; 
        for j = 1:size(data.trial{i},1) % loop over ALL channels
            data.artifactFilter{i}(j,:) = combinedFilter; 
        end
    end
elseif gcfg.loadartifactfilter == 0
    for i = 1:size(data.trial,2) % loop over trials
        data.artifactFilter{i} = [logical(ones(size(data.trial{i})))]; 
        for j = 1:size(data.trial{i},1) % loop over ALL channels
            for k = 1:size(data.artifact{j},1) % loop over artifacts 
                startArtifact = max(1, data.artifact{j}(k,1) - gcfg.artifactPaddingPre * data.fsample); % artifact time point minus pre padding (but minimally first data point)
                endArtifact = min(size(data.trial{i},2), data.artifact{j}(k,2) + gcfg.artifactPaddingPost * data.fsample); % artifact time point minus pre padding (but maximally last data point)
                data.artifactFilter{i}(j,startArtifact:endArtifact) = 0; 
            end
        end
    end
elseif gcfg.loadartifactfilter == 2
     for i = 1:size(data.trial,2) % loop over trials
        data.artifactFilter{i} = [logical(ones(size(data.trial{i})))]; 
     end
end


%% epoch data for FFT
% select epochs
wake_ok = 0; NREM_ok = 0; REM_ok = 0;
for i = 1:size(data.trial,2) % loop over trials
    channums = find(ismember(data.label, data.label))'; % gauge channels
    stages = [0];
    wakeGaugeFilter{i} = all([data.artifactFilter{i}(channums,:); ismember(gcfg.scoring(1,:),stages)],1); % only artefact free time points in gauge channels and gauge sleep stages 
    if any(wakeGaugeFilter{i}), wake_ok = 1; else display('No artifact-free wake data available for analysis!'); end
    stages = [2 3 4];
    NREMGaugeFilter{i} = all([data.artifactFilter{i}(channums,:); ismember(gcfg.scoring(1,:),stages)],1); % only artefact free time points in gauge channels and gauge sleep stages     
    if any(NREMGaugeFilter{i}), NREM_ok = 1; else display('No artifact-free NREM data available for analysis!'); end
    stages = [5];
    REMGaugeFilter{i} = all([data.artifactFilter{i}(channums,:); ismember(gcfg.scoring(1,:),stages)],1); % only artefact free time points in gauge channels and gauge sleep stages     
    if any(REMGaugeFilter{i}), REM_ok = 1; else display('No artifact-free REM data available for analysis!'); end
end

% extract trials for wake FFT
if wake_ok
    onsets = uint32(find(diff([0 wakeGaugeFilter{i}])==1));
    offsets = uint32(find(diff(wakeGaugeFilter{i})==-1));
    numepochs = uint32(min(numel(onsets),numel(offsets)));
    cfg = [];
    cfg.trlx = uint32([onsets(1:numepochs)', offsets(1:numepochs')', ones(numepochs,1)]); % segments of differnt (maximal) length
    cfg.trl = uint32([]);
    shift = uint32(gcfg.fftwinlength * (1-gcfg.fftwinoverlap) * data.fsample);
    tl = uint32(gcfg.fftwinlength * data.fsample);
    for k = 1:size(cfg.trlx,1)
        newonset = cfg.trlx(k,1);
        newoffset = uint32(newonset + tl - 1);
        while cfg.trlx(k,2) - newonset - 1 >= tl
            cfg.trl = [cfg.trl; newonset newoffset 1]; % add segment
            newonset = newonset + shift;% shift data
            newoffset = uint32(newonset + tl - 1);
        end
    end
    wakeFFTdata = ft_redefinetrial(cfg,data);
end

% extract trials for NREM FFT
if NREM_ok
    onsets = uint32(find(diff(NREMGaugeFilter{i})==1));
    offsets = uint32(find(diff(NREMGaugeFilter{i})==-1));
    numepochs = uint32(min(numel(onsets),numel(offsets)));
    cfg = [];
    cfg.trlx = uint32([onsets(1:numepochs)', offsets(1:numepochs')', ones(numepochs,1)]); % segments of differnt (maximal) length
    cfg.trl = uint32([]);
    shift = uint32(gcfg.fftwinlength * (1-gcfg.fftwinoverlap) * data.fsample);
    tl = uint32(gcfg.fftwinlength * data.fsample);
    for k = 1:size(cfg.trlx,1)
        newonset = cfg.trlx(k,1);
        newoffset = uint32(newonset + tl - 1);
        while cfg.trlx(k,2) - newonset - 1 >= tl
            cfg.trl = [cfg.trl; newonset newoffset 1]; % add segment
            newonset = newonset + shift;% shift data
            newoffset = uint32(newonset + tl - 1);
        end
    end
    NREMFFTdata = ft_redefinetrial(cfg,data);
end

% extract trials for REM FFT
if REM_ok
    onsets = uint32(find(diff(REMGaugeFilter{i})==1));
    offsets = uint32(find(diff(REMGaugeFilter{i})==-1));
    numepochs = uint32(min(numel(onsets),numel(offsets)));
    cfg = [];
    cfg.trlx = uint32([onsets(1:numepochs)', offsets(1:numepochs')', ones(numepochs,1)]); % segments of differnt (maximal) length
    cfg.trl = uint32([]);
    shift = uint32(gcfg.fftwinlength * (1-gcfg.fftwinoverlap) * data.fsample);
    tl = uint32(gcfg.fftwinlength * data.fsample);
    for k = 1:size(cfg.trlx,1)
        newonset = cfg.trlx(k,1);
        newoffset = uint32(newonset + tl - 1);
        while cfg.trlx(k,2) - newonset - 1 >= tl
            cfg.trl = [cfg.trl; newonset newoffset 1]; % add segment
            newonset = newonset + shift;% shift data
            newoffset = uint32(newonset + tl - 1);
        end
    end
    REMFFTdata = ft_redefinetrial(cfg,data);
end

%% Irregular resampling auto-spectral analysis (IRASA)
if strcmp(gcfg.normalizepower,'irasa')

    % parameters
    spec = struct;  temp = struct; mspec = struct;
    FrangePlaxfit = [gcfg.foi(1) gcfg.foi(end)]; % define frequency range for power-law fitting
    hsetVec = [1.1:0.05:1.9];
    
    % wake
    if wake_ok
        for iChan = 1:length(data.label)
            for iTrial = 1:length(wakeFFTdata.trial)                
                sig.wake(iChan,:,iTrial) = wakeFFTdata.trial{iTrial}(iChan,:);
            end            
            display(['Computing channel ' data.label{iChan} ', number ' num2str(iChan) ' of ' num2str(length(data.label)) ' channels for wake.']);
            spec.wake(iChan,:,:) = amri_sig_fractal(squeeze(sig.wake(iChan,:,:)),data.fsample,'detrend',1,'frange',[gcfg.foi(1) gcfg.foi(end)+1],'hset',hsetVec);
        end
        for iChan = 1:length(data.label)
            for iTrial = 1:length(wakeFFTdata.trial)
                temp.wake(iChan).mixd(:,iTrial) = interpn(spec.wake(iChan).freq,spec.wake(iChan).mixd(:,iTrial),gcfg.foi,'cubic')';
                temp.wake(iChan).frac(:,iTrial) = interpn(spec.wake(iChan).freq,spec.wake(iChan).frac(:,iTrial),gcfg.foi,'cubic')';
                temp.wake(iChan).osci(:,iTrial) = interpn(spec.wake(iChan).freq,spec.wake(iChan).osci(:,iTrial),gcfg.foi,'cubic')';
            end            
            temp.wake(iChan).freq = gcfg.foi;            
            temp.wake(iChan).srate = spec.wake.srate;
            T = amri_sig_plawfit(temp.wake(iChan),FrangePlaxfit);               
            mspec.wake.freq = gcfg.foi;
            mspec.wake.srate = spec.wake.srate;
            mspec.wake.mixd(iChan,:) = squeeze(squeeze(mean(T.mixd,2)));
            mspec.wake.frac(iChan,:) = squeeze(squeeze(mean(T.frac,2)));
            mspec.wake.osci(iChan,:) = squeeze(squeeze(mean(T.osci,2)));
            mspec.wake.Beta(iChan,:) = squeeze(squeeze(mean(T.Beta,1)));
            mspec.wake.Cons(iChan,:) = squeeze(squeeze(mean(T.Cons,1)));
            mspec.wake.Plaw(iChan,:) = squeeze(squeeze(mean(T.Plaw,2)));
            clear T;
        end
    end
    
    % NREM
    if NREM_ok
        for iChan = 1:length(data.label)
            for iTrial = 1:length(NREMFFTdata.trial)                
                sig.NREM(iChan,:,iTrial) = NREMFFTdata.trial{iTrial}(iChan,:);
            end            
            display(['Computing channel ' data.label{iChan} ', number ' num2str(iChan) ' of ' num2str(length(data.label)) ' channels for NREM.']);
            spec.NREM(iChan,:,:) = amri_sig_fractal(squeeze(sig.NREM(iChan,:,:)),data.fsample,'detrend',1,'frange',[gcfg.foi(1) gcfg.foi(end)+1],'hset',hsetVec);
        end
        for iChan = 1:length(data.label)
            for iTrial = 1:length(NREMFFTdata.trial)
                temp.NREM(iChan).mixd(:,iTrial) = interpn(spec.NREM(iChan).freq,spec.NREM(iChan).mixd(:,iTrial),gcfg.foi,'cubic')';
                temp.NREM(iChan).frac(:,iTrial) = interpn(spec.NREM(iChan).freq,spec.NREM(iChan).frac(:,iTrial),gcfg.foi,'cubic')';
                temp.NREM(iChan).osci(:,iTrial) = interpn(spec.NREM(iChan).freq,spec.NREM(iChan).osci(:,iTrial),gcfg.foi,'cubic')';
            end            
            temp.NREM(iChan).freq = gcfg.foi;            
            temp.NREM(iChan).srate = spec.NREM.srate;
            T = amri_sig_plawfit(temp.NREM(iChan),FrangePlaxfit);               
            mspec.NREM.freq = gcfg.foi;
            mspec.NREM.srate = spec.NREM.srate;
            mspec.NREM.mixd(iChan,:) = squeeze(squeeze(mean(T.mixd,2)));
            mspec.NREM.frac(iChan,:) = squeeze(squeeze(mean(T.frac,2)));
            mspec.NREM.osci(iChan,:) = squeeze(squeeze(mean(T.osci,2)));
            mspec.NREM.Beta(iChan,:) = squeeze(squeeze(mean(T.Beta,1)));
            mspec.NREM.Cons(iChan,:) = squeeze(squeeze(mean(T.Cons,1)));
            mspec.NREM.Plaw(iChan,:) = squeeze(squeeze(mean(T.Plaw,2)));
            clear T;
        end
    end    
    
     % REM
    if REM_ok
        for iChan = 1:length(data.label)
            for iTrial = 1:length(REMFFTdata.trial)                
                sig.REM(iChan,:,iTrial) = REMFFTdata.trial{iTrial}(iChan,:);
            end            
            display(['Computing channel ' data.label{iChan} ', number ' num2str(iChan) ' of ' num2str(length(data.label)) ' channels for REM.']);
            spec.REM(iChan,:,:) = amri_sig_fractal(squeeze(sig.REM(iChan,:,:)),data.fsample,'detrend',1,'frange',[gcfg.foi(1) gcfg.foi(end)+1],'hset',hsetVec);
        end
        for iChan = 1:length(data.label)
            for iTrial = 1:length(REMFFTdata.trial)
                temp.REM(iChan).mixd(:,iTrial) = interpn(spec.REM(iChan).freq,spec.REM(iChan).mixd(:,iTrial),gcfg.foi,'cubic')';
                temp.REM(iChan).frac(:,iTrial) = interpn(spec.REM(iChan).freq,spec.REM(iChan).frac(:,iTrial),gcfg.foi,'cubic')';
                temp.REM(iChan).osci(:,iTrial) = interpn(spec.REM(iChan).freq,spec.REM(iChan).osci(:,iTrial),gcfg.foi,'cubic')';
            end            
            temp.REM(iChan).freq = gcfg.foi;            
            temp.REM(iChan).srate = spec.REM.srate;
            T = amri_sig_plawfit(temp.REM(iChan),FrangePlaxfit);               
            mspec.REM.freq = gcfg.foi;
            mspec.REM.srate = spec.REM.srate;
            mspec.REM.mixd(iChan,:) = squeeze(squeeze(mean(T.mixd,2)));
            mspec.REM.frac(iChan,:) = squeeze(squeeze(mean(T.frac,2)));
            mspec.REM.osci(iChan,:) = squeeze(squeeze(mean(T.osci,2)));
            mspec.REM.Beta(iChan,:) = squeeze(squeeze(mean(T.Beta,1)));
            mspec.REM.Cons(iChan,:) = squeeze(squeeze(mean(T.Cons,1)));
            mspec.REM.Plaw(iChan,:) = squeeze(squeeze(mean(T.Plaw,2)));
            clear T;
        end
    end        
    
clear spec temp; 

end % of IRASA

%% calculate FFTs
if ~strcmp(gcfg.normalizepower,'irasa')
    
    % wake FFT
    if wake_ok
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.output = 'pow';
        cfg.channel = data.label;
        cfg.keeptrials = 'no';
        cfg.foi = gcfg.foi;
        cfg.taper = 'hanning';
        si.freq.wake = ft_freqanalysis(cfg,wakeFFTdata);
    end
    
    % NREM FFT
    if NREM_ok
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.output = 'pow';
        cfg.channel = data.label;
        cfg.keeptrials = 'no';
        cfg.foi = gcfg.foi;
        cfg.taper = 'hanning';
        si.freq.NREM = ft_freqanalysis(cfg,NREMFFTdata);
    end
    
    % REM FFT
    if REM_ok
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.output = 'pow';
        cfg.channel = data.label;
        cfg.keeptrials = 'no';
        cfg.foi = gcfg.foi;
        cfg.taper = 'hanning';
        si.freq.REM = ft_freqanalysis(cfg,REMFFTdata);
    end
    
elseif strcmp(gcfg.normalizepower,'irasa')
    if wake_ok
        si.freq.wake.label = data.label;
        si.freq.wake.dimord = 'chan_freq';
        si.freq.wake.freq = mspec.wake.freq;
        si.freq.wake.powspctrm = mspec.wake.mixd';
    end
    if NREM_ok
        si.freq.NREM.label = data.label;
        si.freq.NREM.dimord = 'chan_freq';
        si.freq.NREM.freq = mspec.NREM.freq;
        si.freq.NREM.powspctrm = mspec.NREM.mixd';
    end
    if REM_ok
        si.freq.REM.label = data.label;
        si.freq.REM.dimord = 'chan_freq';
        si.freq.REM.freq = mspec.REM.freq;
        si.freq.REM.powspctrm = mspec.REM.mixd';
    end
end

%% normalize power for 1/f distribution 
switch gcfg.normalizepower  
    case 'none'  
        if wake_ok
            si.freq.wake_norm = si.freq.wake;
        end
        if NREM_ok
            si.freq.NREM_norm = si.freq.NREM;
        end
        if REM_ok
            si.freq.REM_norm= si.freq.REM;
        end
    case '1overf'    
        if wake_ok
            si.freq.wake_norm = si.freq.wake;
            si.freq.wake_norm.powspctrm = si.freq.wake.powspctrm .*  repmat(si.freq.wake.freq,length(si.freq.wake.label),1).^2;
        end
        if NREM_ok
            si.freq.NREM_norm = si.freq.NREM;
            si.freq.NREM_norm.powspctrm = si.freq.NREM.powspctrm .*  repmat(si.freq.NREM.freq,length(si.freq.NREM.label),1).^2;
        end
        if REM_ok
            si.freq.REM_norm = si.freq.REM;
            si.freq.REM_norm.powspctrm = si.freq.REM.powspctrm .*  repmat(si.freq.REM.freq,length(si.freq.REM.label),1).^2;
        end
    case 'irasa'        
        if wake_ok
            si.freq.wake_norm = si.freq.wake;
            si.freq.wake_norm.powspctrm = mspec.wake.osci;
            si.freq.wake_norm.mixd = mspec.wake.mixd;
            si.freq.wake_norm.frac = mspec.wake.frac;
            si.freq.wake_norm.rel = mspec.wake.osci .* (ones(size(mspec.wake.osci)) ./ mspec.wake.frac);
            si.freq.wake_norm.Beta = mspec.wake.Beta;
            si.freq.wake_norm.Cons = mspec.wake.Cons;
            si.freq.wake_norm.Plaw = mspec.wake.Plaw;
        end
        if NREM_ok
            si.freq.NREM_norm = si.freq.NREM;
            si.freq.NREM_norm.powspctrm = mspec.NREM.osci;
            si.freq.NREM_norm.mixd = mspec.NREM.mixd;
            si.freq.NREM_norm.frac = mspec.NREM.frac;
            si.freq.NREM_norm.rel = mspec.NREM.osci .* (ones(size(mspec.NREM.osci)) ./ mspec.NREM.frac);
            si.freq.NREM_norm.Beta = mspec.NREM.Beta;
            si.freq.NREM_norm.Cons = mspec.NREM.Cons;
            si.freq.NREM_norm.Plaw = mspec.NREM.Plaw;
        end
        if REM_ok
            si.freq.REM_norm = si.freq.REM;
            si.freq.REM_norm.powspctrm = mspec.REM.osci;
            si.freq.REM_norm.mixd = mspec.REM.mixd;
            si.freq.REM_norm.frac = mspec.REM.frac;
            si.freq.REM_norm.rel = mspec.REM.osci .* (ones(size(mspec.REM.osci)) ./ mspec.REM.frac);
            si.freq.REM_norm.Beta = mspec.REM.Beta;
            si.freq.REM_norm.Cons = mspec.REM.Cons;
            si.freq.REM_norm.Plaw = mspec.REM.Plaw;
        end
end

%% add dummy strucures if one sleep stage is missing
if ~wake_ok
    si.freq.wake.label = data.label;
    si.freq.wake.dimord = 'chan_freq';
    si.freq.wake.freq = NaN;
    si.freq.wake.powspctrm = NaN;
    si.freq.wake_norm.label = data.label;
    si.freq.wake_norm.dimord = 'chan_freq';
    si.freq.wake_norm.freq = NaN;
    si.freq.wake_norm.powspctrm = NaN;    
end
if ~NREM_ok
    si.freq.NREM.label = data.label;
    si.freq.NREM.dimord = 'chan_freq';
    si.freq.NREM.freq = NaN;
    si.freq.NREM.powspctrm = NaN;
    si.freq.NREM_norm.label = data.label;
    si.freq.NREM_norm.dimord = 'chan_freq';
    si.freq.NREM_norm.freq = NaN;
    si.freq.NREM_norm.powspctrm = NaN;    
end
if ~REM_ok
    si.freq.REM.label = data.label;
    si.freq.REM.dimord = 'chan_freq';
    si.freq.REM.freq = NaN;
    si.freq.REM.powspctrm = NaN;
    si.freq.REM_norm.label = data.label;
    si.freq.REM_norm.dimord = 'chan_freq';
    si.freq.REM_norm.freq = NaN;
    si.freq.REM_norm.powspctrm = NaN;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nchan = length(data.label);
si.figurehandle.hypno = figure;
switch gcfg.fftplotdim
    case 'channel'
        plotcols = nchan;
    case 'sleepstage'
        plotcols = 3; % at the moment only wake, NREM and REM are considered
end

%% hypnogramme
subplot(4,plotcols,1:plotcols);
cla;

plot_scoring = single(gcfg.scoring(1,:));
plot_scoring = plot_scoring(1:10*data.fsample:end); % one value per 10 s epoch
plot_scoring(plot_scoring == 5) = 0.5;
plot_time = (1:length(plot_scoring))/(6*60); % in hours

plot_wake = plot_scoring;
plot_wake(~ismember(plot_wake,[0])) = NaN;
plot_timewake = plot_time; 
plot_timewake(~ismember(plot_wake,[0])) = NaN;

plot_REM = plot_scoring;
plot_REM(~ismember(plot_REM,[0.5])) = NaN;
plot_timeREM = plot_time; 
plot_timeREM(~ismember(plot_REM,[0.5])) = NaN;

plot_NREM = plot_scoring;
plot_NREM(~ismember(plot_NREM,[2 3 4])) = NaN;
plot_timeNREM = plot_time; 
plot_timeNREM(~ismember(plot_NREM,[2 3 4])) = NaN;

hold on
plot(plot_time,plot_scoring','k','LineWidth',2);
plot(plot_timewake,plot_wake','g','LineWidth',2);
plot(plot_timeNREM,plot_NREM','b','LineWidth',2);
plot(plot_timeREM,plot_REM','r','LineWidth',2);
hold off;
ylim([-0.5 4.5]);
set(gca,'YDir','reverse');
set(gca,'ytick',[0 0.5 1 2 3 4],'yticklabel',{'wake','REM','S1','S2','S3','S4'});
xlim([0 ceil(plot_time(end))]);
set(gca,'xtick',0:1:ceil(plot_time(end)));
ylabel('sleep stage');
xlabel('time (hours)');


%% characterize sleep architecture

sleepOnset = find(ismember(gcfg.scoring(1,:),[1 2 3 4 5]),1,'first');
awakening = find(ismember(gcfg.scoring(1,:),[1 2 3 4 5]),1,'last');

si.architect.SleepOnset = find(ismember(gcfg.scoring(1,:),[1 2 3 4 5]),1,'first')/data.fsample/60;
si.architect.NREMOnset = find(ismember(gcfg.scoring(1,:),[2 3 4 5]),1,'first')/data.fsample/60;
si.architect.REMOnset = find(ismember(gcfg.scoring(1,:),[5]),1,'first')/data.fsample/60;
si.architect.Awakening = find(ismember(gcfg.scoring(1,:),[1 2 3 4 5]),1,'last')/data.fsample/60;

si.architect.minTotal = (length(gcfg.scoring(1,:))/data.fsample/60);
si.architect.minWake = sum(ismember(gcfg.scoring(1,:),[0]))/data.fsample/60;
si.architect.minSleep = sum(ismember(gcfg.scoring(1,:),[1 2 3 4 5]))/data.fsample/60;
si.architect.minNREM = sum(ismember(gcfg.scoring(1,:),[1 2 3 4]))/data.fsample/60;
si.architect.minS1 = sum(ismember(gcfg.scoring(1,:),[1]))/data.fsample/60;
si.architect.minS2 = sum(ismember(gcfg.scoring(1,:),[2]))/data.fsample/60;
si.architect.minS3 = sum(ismember(gcfg.scoring(1,:),[3]))/data.fsample/60;
si.architect.minS4 = sum(ismember(gcfg.scoring(1,:),[4]))/data.fsample/60;
si.architect.minSWS = sum(ismember(gcfg.scoring(1,:),[3 4]))/data.fsample/60; 
si.architect.minREM = sum(ismember(gcfg.scoring(1,:),[5]))/data.fsample/60;

si.architect.minSleepOnset2wakeOnset = si.architect.Awakening - si.architect.SleepOnset;
si.architect.minWakeAfterSleepOnset = sum(ismember(gcfg.scoring(1,(si.architect.SleepOnset*data.fsample*60):(si.architect.Awakening*data.fsample*60)),[5]))/data.fsample/60;
si.architect.minWakeDuringSleep = si.architect.minSleepOnset2wakeOnset - si.architect.minSleep;
si.architect.percWake = si.architect.minWake/si.architect.minTotal*100;
si.architect.percSleep = si.architect.minSleep/si.architect.minTotal*100;
si.architect.percNREM = si.architect.minNREM/si.architect.minTotal*100;
si.architect.percS1 = si.architect.minS1/si.architect.minTotal*100;
si.architect.percS2 = si.architect.minS2/si.architect.minTotal*100;
si.architect.percS3 = si.architect.minS3/si.architect.minTotal*100;
si.architect.percS4 = si.architect.minS4/si.architect.minTotal*100;
si.architect.percSWS = si.architect.minSWS/si.architect.minTotal*100;
si.architect.percREM = si.architect.minREM/si.architect.minTotal*100;

si.analyzedTimeInSec.wake = sum(wakeGaugeFilter{1},2)/data.fsample;
si.analyzedTimeInSec.NREM = sum(NREMGaugeFilter{1},2)/data.fsample;
si.analyzedTimeInSec.REM = sum(REMGaugeFilter{1},2)/data.fsample;
 
% write table
caseNames = {'total','wake','sleep','NREM','S1','S2','S3','S4','SWS','REM','sleep_onset','NREM_onset','REM_onset', 'awakening','sleep_onset_to_wake_onset','wake_after_sleep_onset','wake_during_sleep'}';
varNames = {'min','percent'};

T = cell2table(...
    {si.architect.minTotal, '100%';...
    si.architect.minWake, si.architect.percWake;...
    si.architect.minSleep, si.architect.percSleep;...
    si.architect.minNREM, si.architect.percNREM;...
    si.architect.minS1, si.architect.percS1;...
    si.architect.minS2, si.architect.percS2;...
    si.architect.minS3, si.architect.percS3;...
    si.architect.minS4, si.architect.percS4;...
    si.architect.minSWS, si.architect.percSWS;...
    si.architect.minREM, si.architect.percREM;...
    si.architect.SleepOnset, 'N.A.';...
    si.architect.NREMOnset, 'N.A.';...
    si.architect.REMOnset, 'N.A.';...
    si.architect.Awakening, 'N.A.';...
    si.architect.minSleepOnset2wakeOnset, 'N.A.';...
    si.architect.minWakeAfterSleepOnset, 'N.A.';...
    si.architect.minWakeDuringSleep, 'N.A.'},...
    'RowNames', caseNames, 'VariableNames',varNames');
% 
% 
% T = table(...
%     {si.architect.minTotal, '100%'},...
%     {si.architect.minWake, si.architect.percWake},...
%     {si.architect.minSleep, si.architect.percSleep},...
%     {si.architect.minNREM, si.architect.percNREM},...
%     {si.architect.minS1, si.architect.percS1},...
%     {si.architect.minS2, si.architect.percS2},...
%     {si.architect.minS3, si.architect.percS3},...
%     {si.architect.minS4, si.architect.percS4},...
%     {si.architect.minSWS, si.architect.percSWS},...
%     {si.architect.minREM, si.architect.percREM},...
%     {si.architect.SleepOnset, 'N.A.'},...
%     {si.architect.NREMOnset, 'N.A.'},...
%     {si.architect.REMOnset, 'N.A.'},...
%     {si.architect.Awakening, 'N.A.'},...
%     {si.architect.minSleepOnset2wakeOnset, 'N.A.'},...
%     {si.architect.minWakeAfterSleepOnset, 'N.A.'},...
%     {si.architect.minWakeDuringSleep, 'N.A.'},...
%     'RowNames', caseNames, 'VariableNames',varNames');

% plot table
subplot(4,plotcols,plotcols+1:2*plotcols);
% posT = get(subplot(3,plotcols,plotcols+1:plotcols+3),'Position');
posT = get(subplot(4,plotcols,plotcols+1:2*plotcols),'Position');
% posT = posT .* [1 1.5 1 0.5];
set(subplot(4,plotcols,plotcols+1:2*plotcols),'YTick',[],'XTick',[]);
% uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',posT);
uitable_handle = uitable('Data',table2cell(T)','RowName',T.Properties.VariableNames,'ColumnName',T.Properties.RowNames,'Units','Normalized','Position',posT);
% set(uitable_handle, 'ColumnWidth', {'auto','auto'});


%% plot FFTs 
switch gcfg.normalizepower  
    case 'none'
        plotpow.wake = si.freq.wake.powspctrm; 
        plotpow.NREM = si.freq.NREM.powspctrm;
        plotpow.REM = si.freq.REM.powspctrm;  
    case '1overf'
        plotpow.wake = si.freq.wake_norm.powspctrm; 
        plotpow.NREM = si.freq.NREM_norm.powspctrm;
        plotpow.REM = si.freq.REM_norm.powspctrm;        
    case 'irasa'
        plotpow.wake = si.freq.wake_norm.rel; 
        plotpow.NREM = si.freq.NREM_norm.rel;
        plotpow.REM = si.freq.REM_norm.rel;
end
switch gcfg.fftplotdim
    case 'channel'
               
        % plot absolute FFTs
        firstplot = 1;
        for j = find(ismember(data.label, data.label))'
            sph = subplot(4,nchan,2*nchan+j);
            cla;
            hold on;
            plot(si.freq.wake_norm.freq(si.freq.wake_norm.freq >= gcfg.absfreqwin(1) & si.freq.wake_norm.freq <= gcfg.absfreqwin(2)), plotpow.wake(j,si.freq.wake_norm.freq >= gcfg.absfreqwin(1) & si.freq.wake_norm.freq <= gcfg.absfreqwin(2)),'g');
            plot(si.freq.NREM_norm.freq(si.freq.NREM_norm.freq >= gcfg.absfreqwin(1) & si.freq.NREM_norm.freq <= gcfg.absfreqwin(2)), plotpow.NREM(j,si.freq.NREM_norm.freq >= gcfg.absfreqwin(1) & si.freq.NREM_norm.freq <= gcfg.absfreqwin(2)),'b');
            plot(si.freq.REM_norm.freq(si.freq.REM_norm.freq >= gcfg.absfreqwin(1) & si.freq.REM_norm.freq <= gcfg.absfreqwin(2)), plotpow.REM(j,si.freq.REM_norm.freq >= gcfg.absfreqwin(1) & si.freq.REM_norm.freq <= gcfg.absfreqwin(2)),'r');
            hold off;
            xlabel({'Hz'});
            if firstplot
                switch gcfg.normalizepower
                    case 'none'
                        ylabel({'absolute power','(same scaling per channel)'});
                    case '1overf'
                        ylabel({'1/f corrected power','(same scaling per channel)'});
                    case 'irasa'
                        ylabel({'% change from 1/f','(same scaling per channel)'});
                end                
            else
                set(gca,'yticklabels','');
            end
            %     maxPow =  max(max([plotpow.wake(:,si.freq.wake_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.wake_norm.freq <= gcfg.maxadjfreqwin(2)) plotpow.NREM(:,si.freq.NREM_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.NREM_norm.freq <= gcfg.maxadjfreqwin(2)) plotpow.REM(:,si.freq.REM_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.REM_norm.freq <= gcfg.maxadjfreqwin(2))],[],2))*1.1;
            maxPow =  max(max([plotpow.NREM(:,si.freq.NREM_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.NREM_norm.freq <= gcfg.maxadjfreqwin(2)) plotpow.REM(:,si.freq.REM_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.REM_norm.freq <= gcfg.maxadjfreqwin(2))],[],2))*1.1; % without wake as it is usually full of artifacts
            if ~any(maxPow), maxPow = 1; end
            ylim([0 maxPow]);
            xlim( gcfg.absfreqwin);
            title(data.label{j});
%             pos = get(sph,'position');
%             pos(1) = pos(1)-pos(3)*0.2;
%             pos(3) = pos(3)+pos(3)*0.2;
%             set(sph, 'position', pos);
            firstplot = 0;
        end
        
        % plot relative FFTs
        firstplot = 1;
        for j = find(ismember(data.label, data.label))'
            sph = subplot(4,nchan,3*nchan+j);
            cla;
            hold on;
            plot(si.freq.wake_norm.freq(si.freq.wake_norm.freq >= gcfg.relfreqwin(1) & si.freq.wake_norm.freq <= gcfg.relfreqwin(2)), plotpow.wake(j,si.freq.wake_norm.freq >= gcfg.relfreqwin(1) & si.freq.wake_norm.freq <= gcfg.relfreqwin(2)),'g');
            plot(si.freq.NREM_norm.freq(si.freq.NREM_norm.freq >= gcfg.relfreqwin(1) & si.freq.NREM_norm.freq <= gcfg.relfreqwin(2)), plotpow.NREM(j,si.freq.NREM_norm.freq >= gcfg.relfreqwin(1) & si.freq.NREM_norm.freq <= gcfg.relfreqwin(2)),'b');
            plot(si.freq.REM_norm.freq(si.freq.REM_norm.freq >= gcfg.relfreqwin(1) & si.freq.REM_norm.freq <= gcfg.relfreqwin(2)), plotpow.REM(j,si.freq.REM_norm.freq >= gcfg.relfreqwin(1) & si.freq.REM_norm.freq <= gcfg.relfreqwin(2)),'r');
            hold off;
            xlabel({'Hz'});
            if firstplot                
                switch gcfg.normalizepower
                    case 'none'
                        ylabel({'absolute power','(different scaling per channel)'});
                    case '1overf'
                        ylabel({'1/f corrected power','(different scaling per channel)'});
                    case 'irasa'
                        ylabel({'% change from 1/f','(different scaling per channel)'});
                end                     
            else
                set(gca,'yticklabels','');
            end
            %     maxPow =  max([plotpow.wake(:,si.freq.wake_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.wake_norm.freq <= gcfg.maxadjfreqwin(2)) plotpow.NREM(:,si.freq.NREM_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.NREM_norm.freq <= gcfg.maxadjfreqwin(2)) plotpow.REM(:,si.freq.REM_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.REM_norm.freq <= gcfg.maxadjfreqwin(2))],[],2)*1.1;
            maxPow =  max([plotpow.NREM(:,si.freq.NREM_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.NREM_norm.freq <= gcfg.maxadjfreqwin(2)) plotpow.REM(:,si.freq.REM_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.REM_norm.freq <= gcfg.maxadjfreqwin(2))],[],2)*1.1;
            if ~any(maxPow), maxPow = 1; end
            ylim([0 maxPow(j)]);
            xlim(gcfg.relfreqwin);
            title(data.label{j});
%             pos = get(sph,'position');
%             pos(1) = pos(1)-pos(3)*0.2;
%             pos(3) = pos(3)+pos(3)*0.2;
%             set(sph, 'position', pos);
            firstplot = 0;
        end
        
    case 'sleepstage'       
        % plot absolute FFTs        
        
        maxPow =  max(max([plotpow.wake(:,si.freq.wake_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.wake_norm.freq <= gcfg.maxadjfreqwin(2)) plotpow.NREM(:,si.freq.NREM_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.NREM_norm.freq <= gcfg.maxadjfreqwin(2)) plotpow.REM(:,si.freq.REM_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.REM_norm.freq <= gcfg.maxadjfreqwin(2))],[],2))*1.1;
        if ~any(maxPow), maxPow = 1; end
        cc = distinguishable_colors(length(data.label));
        
        % wake
        sph = subplot(4,plotcols,2*plotcols + 1);
        cla;
        hold on;
        for i = 1:length(data.label)
            plot(si.freq.wake_norm.freq(si.freq.wake_norm.freq >= gcfg.absfreqwin(1) & si.freq.wake_norm.freq <= gcfg.absfreqwin(2)), plotpow.wake(i,si.freq.wake_norm.freq >= gcfg.absfreqwin(1) & si.freq.wake_norm.freq <= gcfg.absfreqwin(2)),'color',cc(i,:));
        end
        hold off;
        xlabel({'Hz'});
        switch gcfg.normalizepower
            case 'none'
                ylabel({'absolute power','(same scaling per sleep stage)'});
            case '1overf'
                ylabel({'1/f corrected power','(same scaling per sleep stage)'});
            case 'irasa'
                ylabel({'% change from 1/f','(same scaling per sleep stage)'});
        end
        ylim([0 maxPow]);
        xlim(gcfg.absfreqwin);
        title('wake');
        legend(data.label,'Location','northeast','Orientation','vertical');
        
        % NREM
        sph = subplot(4,plotcols,2*plotcols + 2);
        cla;
        hold on;
        for i = 1:length(data.label)
            plot(si.freq.NREM_norm.freq(si.freq.NREM_norm.freq >= gcfg.absfreqwin(1) & si.freq.NREM_norm.freq <= gcfg.absfreqwin(2)), plotpow.NREM(i,si.freq.NREM_norm.freq >= gcfg.absfreqwin(1) & si.freq.NREM_norm.freq <= gcfg.absfreqwin(2)),'color',cc(i,:));
        end
        hold off;
        xlabel({'Hz'});
        switch gcfg.normalizepower
            case 'none'
                ylabel({'absolute power','(same scaling per sleep stage)'});
            case '1overf'
                ylabel({'1/f corrected power','(same scaling per sleep stage)'});
            case 'irasa'
                ylabel({'% change from 1/f','(same scaling per sleep stage)'});
        end
        ylim([0 maxPow]);
        xlim(gcfg.absfreqwin);
        title('NREM');
        legend(data.label,'Location','northeast','Orientation','vertical');
        
        % REM
        sph = subplot(4,plotcols,2*plotcols + 3);
        cla;
        hold on;
        for i = 1:length(data.label)
            plot(si.freq.REM_norm.freq(si.freq.REM_norm.freq >= gcfg.absfreqwin(1) & si.freq.REM_norm.freq <= gcfg.absfreqwin(2)), plotpow.REM(i,si.freq.REM_norm.freq >= gcfg.absfreqwin(1) & si.freq.REM_norm.freq <= gcfg.absfreqwin(2)),'color',cc(i,:));
        end
        hold off;
        xlabel({'Hz'});
        switch gcfg.normalizepower
            case 'none'
                ylabel({'absolute power','(same scaling per sleep stage)'});
            case '1overf'
                ylabel({'1/f corrected power','(same scaling per sleep stage)'});
            case 'irasa'
                ylabel({'% change from 1/f','(same scaling per sleep stage)'});
        end
        ylim([0 maxPow]);
        xlim( gcfg.absfreqwin);
        title('REM');
        legend(data.label,'Location','northeast','Orientation','vertical');
        
        % plot relative FFTs
        % wake
        sph = subplot(4,plotcols,3*plotcols + 1);
        cla;
        hold on;
        for i = 1:length(data.label)
            plot(si.freq.wake_norm.freq(si.freq.wake_norm.freq >= gcfg.relfreqwin(1) & si.freq.wake_norm.freq <= gcfg.relfreqwin(2)), plotpow.wake(i,si.freq.wake_norm.freq >= gcfg.relfreqwin(1) & si.freq.wake_norm.freq <= gcfg.relfreqwin(2)),'color',cc(i,:));
        end
        hold off;
        xlabel({'Hz'});
        switch gcfg.normalizepower
            case 'none'
                ylabel({'absolute power','(different scaling per sleep stage)'});
            case '1overf'
                ylabel({'1/f corrected power','(different scaling per sleep stage)'});
            case 'irasa'
                ylabel({'% change from 1/f','(different scaling per sleep stage)'});
        end
        maxPow = max(max([plotpow.wake(:,si.freq.wake_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.wake_norm.freq <= gcfg.maxadjfreqwin(2))],[],2))*1.1;
        if ~any(maxPow), maxPow = 1; end
        ylim([0 maxPow]);
        xlim(gcfg.relfreqwin);
        title('wake');
        legend(data.label,'Location','northeast','Orientation','vertical');
        
        % NREM
        sph = subplot(4,plotcols,3*plotcols + 2);
        cla;
        hold on;
        for i = 1:length(data.label)
            plot(si.freq.NREM_norm.freq(si.freq.NREM_norm.freq >= gcfg.relfreqwin(1) & si.freq.NREM_norm.freq <= gcfg.relfreqwin(2)), plotpow.NREM(i,si.freq.NREM_norm.freq >= gcfg.relfreqwin(1) & si.freq.NREM_norm.freq <= gcfg.relfreqwin(2)),'color',cc(i,:));
        end
        hold off;
        xlabel({'Hz'});
        switch gcfg.normalizepower
            case 'none'
                ylabel({'absolute power','(different scaling per sleep stage)'});
            case '1overf'
                ylabel({'1/f corrected power','(different scaling per sleep stage)'});
            case 'irasa'
                ylabel({'% change from 1/f','(different scaling per sleep stage)'});
        end
        maxPow = max(max([plotpow.NREM(:,si.freq.NREM_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.NREM_norm.freq <= gcfg.maxadjfreqwin(2))],[],2))*1.1;
        if ~any(maxPow), maxPow = 1; end
        ylim([0 maxPow]);
        xlim(gcfg.relfreqwin);
        title('NREM');
        legend(data.label,'Location','northeast','Orientation','vertical');
        
        % REM
        sph = subplot(4,plotcols,3*plotcols + 3);
        cla;
        hold on;
        for i = 1:length(data.label)
            plot(si.freq.REM_norm.freq(si.freq.REM_norm.freq >= gcfg.relfreqwin(1) & si.freq.REM_norm.freq <= gcfg.relfreqwin(2)), plotpow.REM(i,si.freq.REM_norm.freq >= gcfg.relfreqwin(1) & si.freq.REM_norm.freq <= gcfg.relfreqwin(2)),'color',cc(i,:));
        end
        hold off;
        xlabel({'Hz'});
        switch gcfg.normalizepower
            case 'none'
                ylabel({'absolute power','(different scaling per sleep stage)'});
            case '1overf'
                ylabel({'1/f corrected power','(different scaling per sleep stage)'});
            case 'irasa'
                ylabel({'% change from 1/f','(different scaling per sleep stage)'});
        end
        maxPow = max(max([plotpow.REM(:,si.freq.REM_norm.freq >= gcfg.maxadjfreqwin(1) & si.freq.REM_norm.freq <= gcfg.maxadjfreqwin(2))],[],2))*1.1;
        if ~any(maxPow), maxPow = 1; end
        ylim([0 maxPow]);
        xlim( gcfg.relfreqwin);
        title('REM');
        legend(data.label,'Location','northeast','Orientation','vertical');
end


%% plot topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(gcfg,'topofreq') && ~isempty(gcfg.topofreq) 
    
    si.figurehandle.topo = figure;
    for i = 1:length(gcfg.topofreq)
        
        % Wake
        subplot(3,length(gcfg.topofreq),i);
        if wake_ok
            cfg = [];
            cfg.frequency = gcfg.topofreq{i};
            plotsi.freq.wake= ft_selectdata(cfg,si.freq.wake_norm);
            title(['wake (' num2str(cfg.frequency(1)) ' to ' num2str(cfg.frequency(2)) ' Hz)']);
            cfg = [];
            cfg.layout = gcfg.layout;
            ft_topoplotTFR(cfg,plotsi.freq.wake);
        end
        
        % NREM
        subplot(3,length(gcfg.topofreq),length(gcfg.topofreq) + i);
        if NREM_ok
            cfg = [];
            cfg.frequency = gcfg.topofreq{i};
            plotsi.freq.NREM= ft_selectdata(cfg,si.freq.NREM_norm);
            title(['NREM (' num2str(cfg.frequency(1)) ' to ' num2str(cfg.frequency(2)) ' Hz)']);
            cfg = [];
            cfg.layout = gcfg.layout;
            ft_topoplotTFR(cfg,plotsi.freq.NREM);
        end
        
        % REM
        subplot(3,length(gcfg.topofreq),2 * length(gcfg.topofreq) + i);
        if REM_ok
            cfg = [];
            cfg.frequency = gcfg.topofreq{i};
            plotsi.freq.REM= ft_selectdata(cfg,si.freq.REM_norm);
            title(['REM (' num2str(cfg.frequency(1)) ' to ' num2str(cfg.frequency(2)) ' Hz)']);
            cfg = [];
            cfg.layout = gcfg.layout;
            ft_topoplotTFR(cfg,plotsi.freq.REM);
        end
        
    end % of loop over topoplots
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finishing
si.cfg = gcfg;
ttoc = toc;
display('Sleep architecture inspected.');
display(['Inspecting sleep took ' num2str(ttoc) ' seconds.']);

end % of function
