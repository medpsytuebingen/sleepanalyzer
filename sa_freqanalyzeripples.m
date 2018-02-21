%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% classifyripples
% by Til Ole Bergmann 2013
% last modified 2016/11/22 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% analyzes time-frequency-representation of ripple events based on candidate epochs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ripple_gcfg, ripple_eventdata, ripple_freq] = sa_freqanalyzeripples(gcfg, ripple_gcfg, ripple_eventdata)
display('Classifying ripples...');

% select raw channels only
procchannels = {'_bp','_rms'};
goodchan =[];
for i = 1:length(ripple_eventdata.label)
    for j = 1:length(procchannels) % loop over channels postfixes to exclude
        k = strfind(ripple_eventdata.label{i},procchannels{j});  
        goodchan(i,j) = isempty(k);
    end
end
goodchan = logical(prod(goodchan,2));
rawchannels = ripple_eventdata.label(logical(goodchan));

if ismember('all', gcfg.nnmfchannel)    
    gcfg.nnmfchannel = rawchannels;
end


%% TFR
windowlength = (size(ripple_eventdata.trial{1},2)-1)/ripple_eventdata.fsample; % length of time window
cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.channel = rawchannels; 
cfg.keeptrials = 'yes';
cfg.toi = [-windowlength/2:0.005:windowlength/2];
maxfoi = min(200,ripple_eventdata.fsample/5);
cfg.foi = [5:1:maxfoi]; 
for i = 1:maxfoi, cycles(i) = floor(100/(1000/i)); end % keep windows about 100ms with integer number of cycles
cycles(cycles < 5) = 5; % but at least 5 cycles!
cfg.t_ftimwin = cycles(cfg.foi(1):cfg.foi(end))./cfg.foi;
cfg.output = 'pow';	
cfg.polyremoval = 1; % 0 = mean (default), 1 = linear
ripple_freq = ft_freqanalysis(cfg, ripple_eventdata);   
    

%% rule-based classification
if gcfg.rulebased 

% '%change_baseline'    
blcfreq = ripple_freq;
[Nrpt Nchan Nfreq Ntime] = size(blcfreq.powspctrm);
bl = [-1.5 -1.0];
blwin = blcfreq.time >= bl(1) & blcfreq.time <= bl(2);
baseline = repmat(trimmean(nanmean(blcfreq.powspctrm(:,:,:,blwin),4), 10,'round',1),[Nrpt,1,1,Ntime]); %  M = trimmean(X,PERCENT,FLAG,DIM)
blcfreq.powspctrm = (blcfreq.powspctrm - baseline) ./ baseline .* 100;

% z-transform over freq bins per event
eventw = [-0.01 0.01];
eventwin = blcfreq.time >= eventw(1) & blcfreq.time <= eventw(2);
ripple_event_freq = squeeze(mean(blcfreq.powspctrm(:,:,:,eventwin),4));    
zfreq = (ripple_event_freq - repmat(nanmean(ripple_event_freq,2),[1,size(ripple_event_freq,2)])) ./ repmat(nanstd(ripple_event_freq,1,2),[1,size(ripple_event_freq,2)]);

% detect maximum in search window
rw = [80:140]; % ripple window
clear ecv;
for i = 1:size(zfreq,1)
    [lmval, indd] = lmax(zfreq(i,:),3);
    lmvals{i} = lmval;
    indds{i} = indd;   
    ecv(i) = any(ismember(indds{i},rw));
%     locmax(i) = ismember(indds{i},rw);
    kurt(i) = kurtosis(zfreq(i,rw-5));
end
display(['Percentage of events selected: ' num2str(mean(ecv))]);
ripple_cfg.ecv = ecv; % save in ripple_gcfg structure

ecv = ecv .* kurt > median(kurt);

%% plotting
figure;

% explore selected events only
subplot(2,2,1);
plot(blcfreq.freq,zfreq(ecv,:));
hold on;
plot(blcfreq.freq,nanmean(zfreq(ecv,:),1),'k','LineWidth',3);
hold off;
title(['selected events (' num2str(sum(ecv)) ')']);
xlabel('Hz'); ylabel('z-value');  

% explore timelock of selected events
subplot(2,2,2);
cfg = [];
cfg.trials = find(ecv);
[selected_ripple_timelock] = ft_timelockanalysis(cfg, ripple_eventdata);
cfg = [];
ft_singleplotER(cfg,selected_ripple_timelock);

% explore unselected events only
subplot(2,2,3);
plot(blcfreq.freq,zfreq(~ecv,:));
hold on;
plot(blcfreq.freq,nanmean(zfreq(~ecv,:),1),'k','LineWidth',3);
hold off;
title(['unselected events (' num2str(sum(~ecv)) ')']);
xlabel('Hz'); ylabel('z-value');    

% explore timelock of selected events
subplot(2,2,4);
cfg = [];
cfg.trials = find(~ecv);
[unselected_ripple_timelock] = ft_timelockanalysis(cfg, ripple_eventdata);
cfg = [];
ft_singleplotER(cfg,unselected_ripple_timelock);


end


%% NNMF
if gcfg.nnmf
% nnmf.mat works with matlab2009 but not matlab2012!
keyboard;

    for chan = find(ismember(ripple_eventdata.label, gcfg.nnmfchannel))' % loop over nnmfchannels
        cfg = [];
        cfg.freqsize = size(ripple_freq.powspctrm); % [trials, channels, freqbins, timebins]
        cfg.channel = chan;
        cfg.freqwindow = [21:101];
        cfg.timewindow = [150:250];
        cfg.k = 2; % number of components

        input = reshape(ripple_freq.powspctrm(:,cfg.channel,cfg.freqwindow,cfg.timewindow),cfg.freqsize(1),[]); % reshape

        [w{chan},h{chan}] = nnmf(input,cfg.k);

    % output = reshape(input,cfg.freqsize(1),length(cfg.freqwindow),length(cfg.timewindow)); % re-reshape

    end

end

ripple_freq.cfg = gcfg;
end % of function


