%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sa_freqanalyzeoscillations
% by Til Ole Bergmann 2013
% last modified 2017/01/16 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% analyzes time-frequency-representation of oscillation events based on candidate epochs
%
%
% [gcfg, data] = sa_freqanalyzeoscillations(cfg,data)
%
%  cfg.eventname = string indicating the name of the event
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [oscillation_eventdata, oscillation_freq] = sa_freqanalyzeoscillations(gcfg, oscillation_eventdata)
tic;

%% check config
if ~isfield(gcfg,'eventname')
    gcfg.eventname = 'event';
end

%% preparation
display(['Analyzing time-frequency structure of ' gcfg.eventname ' events...']);

%%
% select raw channels only
procchannels = {'_bp','_rms'};
goodchan =[];
for i = 1:length(oscillation_eventdata.label)
    for j = 1:length(procchannels) % loop over channels postfixes to exclude
        k = strfind(oscillation_eventdata.label{i},procchannels{j});  
        goodchan(i,j) = isempty(k);
    end
end
goodchan = logical(prod(goodchan,2));
rawchannels = oscillation_eventdata.label(logical(goodchan));

if ismember('all', gcfg.nnmfchannel)    
    gcfg.nnmfchannel = rawchannels;
end


%% TFR
% windowlength = (size(oscillation_eventdata.trial{1},2)-1)/oscillation_eventdata.fsample; % length of time window
windowlength = size(oscillation_eventdata.trial{1},2)/oscillation_eventdata.fsample; % length of time window

% standard
cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.channel = rawchannels; 
cfg.keeptrials = 'yes';
% cfg.toi = [-windowlength/2:0.005:windowlength/2];
stepsize = 0.005;
cfg.toi = [(-(windowlength/2)+stepsize):stepsize:windowlength/2];
minfoi = max(1,1/(windowlength/5));
maxfoi = min(200,oscillation_eventdata.fsample/5);
% cfg.foi = [5:1:maxfoi]; 
cfg.foi = [minfoi:1:maxfoi]; 
for i = 1:maxfoi, cycles(i) = floor(100/(1000/i)); end % keep windows about 100ms with integer number of cycles (NN2015)
cycles(cycles < 5) = 5; % but at least 5 cycles! (NN2015)
% for i = 1:maxfoi, cycles(i) = floor(75/(1000/i)); end % keep windows about 75ms with integer number of cycles
% cycles(cycles < 3) = 3; % but at least 3 cycles!
cfg.t_ftimwin = cycles(cfg.foi(1):cfg.foi(end))./cfg.foi;
% cfg.t_ftimwin = 5./cfg.foi;
cfg.output = 'pow';	
cfg.polyremoval = 1; % 0 = mean (default), 1 = linear
oscillation_freq = ft_freqanalysis(cfg, oscillation_eventdata);


% cfg = [];
% cfg.method =  'wavelet';
% cfg.taper = 'hanning';
% cfg.channel = rawchannels; 
% cfg.keeptrials = 'yes';
% cfg.toi = [-windowlength/2:0.01:windowlength/2];
% cfg.foi = [4:1:200]; 
% cfg.width = 9;
% cfg.gwidth = 3;
% cfg.output = 'pow';	
% cfg.polyremoval = 0; % 0 = mean (default), 1 = linear
% oscillation_freq = ft_freqanalysis(cfg, oscillation_eventdata);



% % low frequencies (with hanning tapers)
% cfg = [];
% cfg.method = 'mtmconvol';
% cfg.taper = 'hanning';
% cfg.channel = rawchannels; 
% cfg.keeptrials = 'yes';
% cfg.foi = [4:1:29]; 
% cfg.toi = [-windowlength/2:0.01:windowlength/2];
% cfg.t_ftimwin = 5./cfg.foi;
% cfg.output = 'pow';	
% oscillation_freq = ft_freqanalysis(cfg, oscillation_eventdata);

% % high frequencies (with multitapers)
% cfg = [];
% cfg.method = 'mtmconvol';
% cfg.taper = 'dpss';
% cfg.channel = rawchannels; 
% cfg.keeptrials = 'yes';
% cfg.foi = [30:5:200]; 
% cfg.toi = [-windowlength/2:0.01:windowlength/2];
% cfg.t_ftimwin = ones(1,length(cfg.foi))*0.2;
% cfg.tapsmofrq = ones(1,length(cfg.foi))*10;
% cfg.output = 'pow';	
% oscillation_freq_high = ft_freqanalysis(cfg, oscillation_eventdata);
% 
% % concatenate low and high frequency TFRs
% oscillation_freq.powspctrm = cat(3,oscillation_freq.powspctrm,oscillation_freq_high.powspctrm); % concatenate along dimension 3 (freq)
% oscillation_freq.cumtapcnt = cat(2,oscillation_freq.cumtapcnt,oscillation_freq_high.cumtapcnt); % concatenate along dimension 2 (freq)
% oscillation_freq.freq = [oscillation_freq.freq oscillation_freq_high.freq]; % concatenate vector (freq)
% oscillation_freq.cfg2 = oscillation_freq_high.cfg; % add second cfg field


%% NNMF
if gcfg.nnmf
% [A,Y,numIter,tElapsed,finalResidual]=nmfrule(X,k);
% nnmf.mat works with matlab2009 but not matlab2012!

cfg = [];
cfg.freqsize = size(oscillation_freq.powspctrm); % [trials, channels, freqbins, timebins]
cfg.channel = [1];
cfg.freqwindow = [21:101];
cfg.timewindow = [150:250];
cfg.k = 2; % number of components

input = reshape(oscillation_freq.powspctrm(:,cfg.channel,cfg.freqwindow,cfg.timewindow),cfg.freqsize(1),[]); % reshape

[w,h] = nnmf(input,cfg.k,'algorithm','als');

% output = rehape(input,cfg.freqsize(1),length(cfg.freqwindow),length(cfg.timewindow)); % re-reshape


end

%% finishing
oscillation_freq.cfg = gcfg;
ttoc = toc;
display(['Time-frequency analysis of ' gcfg.eventname ' events took ' num2str(ttoc) ' seconds.']);

end % of function


