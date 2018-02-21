%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% classifyspindles
% by Til Ole Bergmann 2013
% last modified 2016/11/22 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% analyzes time-frequency-representation of spindle events based on candidate epochs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spindle_eventdata, spindle_freq] = sa_freqanalyzespindles(gcfg, spindle_eventdata)
display('Classifying spindles...');

% select raw channels only
procchannels = {'_bp','_rms'};
goodchan =[];
for i = 1:length(spindle_eventdata.label)
    for j = 1:length(procchannels) % loop over channels postfixes to exclude
        k = strfind(spindle_eventdata.label{i},procchannels{j});  
        goodchan(i,j) = isempty(k);
    end
end
goodchan = logical(prod(goodchan,2));
rawchannels = spindle_eventdata.label(logical(goodchan));

if ismember('all', gcfg.nnmfchannel)    
    gcfg.nnmfchannel = rawchannels;
end


%% TFR
windowlength = (size(spindle_eventdata.trial{1},2)-1)/spindle_eventdata.fsample; % length of time window

% standard
cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.channel = rawchannels; 
cfg.keeptrials = 'yes';
cfg.toi = [-windowlength/2:0.005:windowlength/2];
maxfoi = min(200,spindle_eventdata.fsample/5);
cfg.foi = [5:1:maxfoi]; 
for i = 1:maxfoi, cycles(i) = floor(100/(1000/i)); end % keep windows about 100ms with integer number of cycles
cycles(cycles < 5) = 5; % but at least 5 cycles!
cfg.t_ftimwin = cycles(cfg.foi(1):cfg.foi(end))./cfg.foi;
% cfg.t_ftimwin = 5./cfg.foi;
cfg.output = 'pow';	
cfg.polyremoval = 1; % 0 = mean (default), 1 = linear
spindle_freq = ft_freqanalysis(cfg, spindle_eventdata);


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
% spindle_freq = ft_freqanalysis(cfg, spindle_eventdata);



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
% spindle_freq = ft_freqanalysis(cfg, spindle_eventdata);

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
% spindle_freq_high = ft_freqanalysis(cfg, spindle_eventdata);
% 
% % concatenate low and high frequency TFRs
% spindle_freq.powspctrm = cat(3,spindle_freq.powspctrm,spindle_freq_high.powspctrm); % concatenate along dimension 3 (freq)
% spindle_freq.cumtapcnt = cat(2,spindle_freq.cumtapcnt,spindle_freq_high.cumtapcnt); % concatenate along dimension 2 (freq)
% spindle_freq.freq = [spindle_freq.freq spindle_freq_high.freq]; % concatenate vector (freq)
% spindle_freq.cfg2 = spindle_freq_high.cfg; % add second cfg field


%% NNMF
if gcfg.nnmf
% [A,Y,numIter,tElapsed,finalResidual]=nmfrule(X,k);
% nnmf.mat works with matlab2009 but not matlab2012!

cfg = [];
cfg.freqsize = size(spindle_freq.powspctrm); % [trials, channels, freqbins, timebins]
cfg.channel = [1];
cfg.freqwindow = [21:101];
cfg.timewindow = [150:250];
cfg.k = 2; % number of components

input = reshape(spindle_freq.powspctrm(:,cfg.channel,cfg.freqwindow,cfg.timewindow),cfg.freqsize(1),[]); % reshape

[w,h] = nnmf(input,cfg.k,'algorithm','als');

% output = rehape(input,cfg.freqsize(1),length(cfg.freqwindow),length(cfg.timewindow)); % re-reshape


end

spindles_freq.cfg = gcfg;
end % of function


