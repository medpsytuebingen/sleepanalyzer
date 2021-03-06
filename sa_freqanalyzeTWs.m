%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% classifyTWs
% by Til Ole Bergmann 2013
% last modified 2016/11/22 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% analyzes time-frequency-representation of theta wave events based on candidate epochs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TW_eventdata, TW_freq] = sa_freqanalyzeTWs(gcfg, TW_eventdata)
display('Classifying theta waves...');

% select raw channels only
procchannels = {'_bp','_rms'};
goodchan =[];
for i = 1:length(TW_eventdata.label)
    for j = 1:length(procchannels) % loop over channels postfixes to exclude
        k = strfind(TW_eventdata.label{i},procchannels{j});  
        goodchan(i,j) = isempty(k);
    end
end
goodchan = logical(prod(goodchan,2));
rawchannels = TW_eventdata.label(logical(goodchan));

if ismember('all', gcfg.nnmfchannel)    
    gcfg.nnmfchannel = rawchannels;
end

%% TFR
windowlength = (size(TW_eventdata.trial{1},2)-1)/TW_eventdata.fsample; % length of time window
cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.channel = rawchannels; 
cfg.keeptrials = 'yes';
cfg.toi = [-windowlength/2:0.01:windowlength/2];
maxfoi = min(200,TW_eventdata.fsample/5);
cfg.foi = [5:1:maxfoi]; 
for i = 1:maxfoi, cycles(i) = floor(100/(1000/i)); end % keep windows about 100ms with integer number of cycles
cycles(cycles < 5) = 5; % but at least 5 cycles!
cfg.t_ftimwin = cycles(cfg.foi(1):cfg.foi(end))./cfg.foi;
cfg.polyremoval = 0; % 0 = mean (default), 1 = linear
TW_freq = ft_freqanalysis(cfg, TW_eventdata);


%% NNMF
if gcfg.nnmf
% [A,Y,numIter,tElapsed,finalResidual]=nmfrule(X,k);
% nnmf.mat works with matlab2009 but not matlab2012!

cfg = [];
cfg.freqsize = size(TW_freq.powspctrm); % [trials, channels, freqbins, timebins]
cfg.channel = [1];
cfg.freqwindow = [21:101];
cfg.timewindow = [150:250];
cfg.k = 2; % number of components

input = reshape(TW_freq.powspctrm(:,cfg.channel,cfg.freqwindow,cfg.timewindow),cfg.freqsize(1),[]); % reshape

[w,h] = nnmf(input,cfg.k,'algorithm','als');

% output = rehape(input,cfg.freqsize(1),length(cfg.freqwindow),length(cfg.timewindow)); % re-reshape


end

TW_freq.cfg = gcfg;
end % of function


