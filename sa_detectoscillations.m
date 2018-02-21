%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detectoscillations
% by Til Ole Bergmann 2016
% last modified 2017/07/18 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% detects oscillations in EEG data
%
% [gcfg, data, bpdata, rmsdata, supthreshdata] = sa_detectoscillations(cfg,data)
%
%  cfg.eventname = string with the name of the oscillation
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gcfg, data, bpdata, rmsdata, supthreshdata] = sa_detectoscillations(gcfg,data)
tic;

%% check config
if ~isfield(gcfg,'eventname')
    gcfg.eventname = 'event';
end

%% preparation
display(['Detecting ' gcfg.eventname ' events...']);


%% get some values
gcfg.labels = data.label';
if ismember('all', gcfg.gaugechannel)    
    gcfg.gaugechannel = data.label;
end
if ismember('all', gcfg.searchchannel)    
    gcfg.searchchannel = data.label;
end  

%% generate scoring array
switch gcfg.sleepscoringStruct % make sure scoring information is availabvle in data.scoring
    case 'trialinfo', data.scoring = int8(data.trialinfo(gcfg.sleepscoringInfoRow,:));
        data.trialinfo(gcfg.sleepscoringInfoRow) = [];
    case 'scoring'
        data.scoring = int8(data.scoring);
end

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
    combinedFilter = all(filterArray,1); % combined filter is only true wi_osc_gaugefreqrangehen filters are true for all channels 
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
% % select epochs
% for i = 1:size(data.trial,2) % loop over trials
%     channums = find(ismember(data.label, gcfg.gaugechannel))'; % gauge channels
%     stages = [0];
%     wakeGaugeFilter{i} = all([data.artifactFilter{i}(channums,:); ismember(gcfg.scoring(1,:),stages)],1); % only artefact free time points in gauge channels and gauge sleep stages 
%     stages = [2 3 4];
%     NREMGaugeFilter{i} = all([data.artifactFilter{i}(channums,:); ismember(gcfg.scoring(1,:),stages)],1); % only artefact free time points in gauge channels and gauge sleep stages     
% end
% 
% % wake FFT
% onsets = find(diff([0 wakeGaugeFilter{i}])==1);
% offsets = find(diff(wakeGaugeFilter{i})==-1);
% numepochs = min(numel(onsets),numel(offsets));
% cfg = [];
% cfg.minlength = 3;% in seconds
% cfg.trl = [onsets(1:numepochs)', offsets(1:numepochs')', ones(numepochs,1)];
% wakefftdata = ft_redefinetrial(cfg,data);
% 
% % NREM FFT
% onsets = find(diff(NREMGaugeFilter{i})==1);
% offsets = find(diff(NREMGaugeFilter{i})==-1);
% numepochs = min(numel(onsets),numel(offsets));
% cfg = [];
% cfg.minlength = 3;% in seconds
% cfg.trl = [onsets(1:numepochs)', offsets(1:numepochs')', ones(numepochs,1)];
% NREMfftdata = ft_redefinetrial(cfg,data);


% if numel(data.trial) > 1
%     fftdata = data;
%     gcfg.epochScoring = gcfg.scoring;
% elseif numel(data.trial) == 1 
%     cfg = [];
% 	cfg.length = 20;
%     cfg.overlap = 0; 
%     epochStart = [1 : cfg.length*data.fsample : size(data.trial{:},2)];
%     gcfg.epochScoring = gcfg.scoring(:,epochStart); 
%     fftdata = ft_redefinetrial(cfg,data);
% end

%% determine oscillation peak frequency
% % FFT
% cfg = [];
% cfg.method = 'mtmfft';
% cfg.output = 'pow';
% cfg.trials = find(gcfg.epochScoring(3,:));
% cfg.channel = gcfg.gaugechannel;
% cfg.foilim = [1:1/3:35];
% cfg.taper = 'hanning';
% freq = ft_freqanalysis(cfg,fftdata);
% clear fftdata;

% % wake FFT
% cfg = [];
% cfg.method = 'mtmfft';
% cfg.output = 'pow';
% cfg.channel = gcfg.gaugechannel;
% % cfg.foilim = [1 35];
% cfg.foi = [1:1/3:35];
% cfg.taper = 'hanning';
% wakefreq = ft_freqanalysis(cfg,wakefftdata);
% 
% % NREM FFT
% cfg = [];
% cfg.method = 'mtmfft';
% cfg.output = 'pow';
% cfg.channel = gcfg.gaugechannel;
% % cfg.foilim = [1 35];
% cfg.foi = [1:1/3:35];
% cfg.taper = 'hanning';
% NREMfreq = ft_freqanalysis(cfg,NREMfftdata);
% 
% difffreq = NREMfreq;
% difffreq.powspctrm = (NREMfreq.powspctrm ./ wakefreq.powspctrm);
% 
% % smooth freq
% kernel = 5; % 1 = no smoothing
% for j = find(ismember(data.label, gcfg.gaugechannel))'
%     smoothwakefreq(j,:) = conv(wakefreq.powspctrm(j,:),ones(1,kernel),'same');
%     smoothNREMfreq(j,:) = conv(NREMfreq.powspctrm(j,:),ones(1,kernel),'same');
%     smoothdifffreq(j,:) = conv(difffreq.powspctrm(j,:),ones(1,kernel),'same');
% end
% 
% % select sleep oscillation peak automatically
% searchIndex = difffreq.freq > gcfg.gaugefreqrange(1) & difffreq.freq < gcfg.gaugefreqrange(2);
% searchInterval = difffreq.freq(searchIndex);
% % [gcfg.oscillationPeakValue,index] = max(log(freq.powspctrm(:,searchIndex)+1).*repmat(freq.freq(searchIndex),size(freq.powspctrm,1),1).^2,[],2); % correct for 1/f distribution of EEG power spectrum to unambigiously detect oscillation peak
% % [gcfg.oscillationPeakValue,index] = max(log(smoothfreq(:,searchIndex)+1).*repmat(freq.freq(searchIndex),size(smoothfreq,1),1).^3,[],2); % correct for 1/f distribution of EEG power (^3!!) spectrum to unambigiously detect oscillation peak
% % [gcfg.oscillationPeakValue,index] = max(log(smoothfreq(:,searchIndex)+1).*repmat(freq.freq(searchIndex),size(smoothfreq,1),1),[],2); % correct for 1/f distribution of EEG power spectrum to unambigiously detect oscillation peak
% [gcfg.oscillationPeakValue,index] = max(log(smoothdifffreq(:,searchIndex)+1));
% gcfg.oscillationPeakFreq = searchInterval(index);
% % save('-v7.3', fullfile(gcfg.resultspath,[gcfg.subjectName '_oscillation_FFT_stage' char(gcfg.gaugestages+'0') '_' gcfg.timestamp]),'wakefreq', 'NREMfreq', 'difffreq'); 

% 
% % plot power spectrum
% h = figure;
% for j = find(ismember(data.label, gcfg.gaugechannel))'
%     subplot(1,length(gcfg.gaugechannel),j);
%     hold on;
%     plot(NREMfreq.freq,log(smoothwakefreq(j,:)+1),'g');    
%     plot(NREMfreq.freq,log(smoothNREMfreq(j,:)+1),'b');
%     plot(NREMfreq.freq,log(smoothdifffreq(j,:)+1),'r');    
%     hold off;
%     xlabel({'Hz'});
%     ylabel('log(power+1)');
%     title([data.label{j} ', oscillation freq: ' num2str(gcfg.oscillationPeakFreq)]);
% end
% saveas(h, fullfile(gcfg.resultspath,[gcfg.subjectName '_' gcfg.eventname '_FFT_stage' char(gcfg.gaugestages+'0') '_' gcfg.timestamp '.fig'])); 
% close(h);


%% bandpass filter in oscillation range
bpdata = data;
if gcfg.searchfreqindividual == 1 % use individual oscillation fequency filter
    gcfg.bpfreq = [gcfg.oscillationPeakFreq' - repmat(gcfg.searchfreqmargin,length(gcfg.oscillationPeakFreq),1) gcfg.oscillationPeakFreq' +  repmat(gcfg.searchfreqmargin,length(gcfg.oscillationPeakFreq),1)];
    for i = 1:size(data.trial,2) % loop over trials
        for j = find(ismember(data.label, gcfg.searchchannel))' % loop over search channels
            bpdata.trial{i}(j,:) = ft_preproc_bandpassfilter(data.trial{i}(j,:), double(data.fsample), gcfg.bpfreq(j,:), 3*fix(data.fsample/gcfg.bpfreq(j,1)), 'fir', 'twopass'); % Mathilde's formula
        end
    end
elseif gcfg.searchfreqindividual == 0 % use common oscillation fequency filter
    chanfilter = find(ismember(data.label, gcfg.searchchannel))';
    gcfg.bpfreq = repmat(gcfg.gaugefreqrange, length(chanfilter),1);
    for i = 1:size(data.trial,2) % loop over trials
        bpdata.trial{i}(chanfilter,:) = ft_preproc_bandpassfilter(data.trial{i}(chanfilter,:), double(data.fsample), gcfg.bpfreq(1,:), 3*fix(data.fsample/gcfg.bpfreq(1))+1, 'fir', 'twopass'); % Mathilde's formula
    end
end


%% calculate root mean square (RMS) signal 
rmsdata = bpdata;
rmsdata.trial = {};
for i = 1:size(bpdata.trial,2) % loop over trials
    for j = find(ismember(rmsdata.label, gcfg.searchchannel))' % loop over search channels
        tempSquared = (bpdata.trial{i}(j,:)) .^2;
        tempConvolved = conv(tempSquared,ones(1,data.fsample * gcfg.rmswindowlength)/(data.fsample * gcfg.rmswindowlength),'same');
        rmsdata.trial{i}(j,:) = sqrt(tempConvolved);
    end
end

%% determine amplitude threshold
% prepare filter also selecting gauge sleep stages
for i = 1:size(rmsdata.trial,2) % loop over trials
    channums = find(ismember(rmsdata.label, gcfg.gaugechannel))'; % gauge channels
    finalGaugeFilter{i} = all([data.artifactFilter{i}(channums,:); gcfg.scoring(3,:)],1); % only artefact free time points in gauge channels and gauge sleep stages 
end

% gather data for variance threshold calculation
for i = 1:size(rmsdata.trial,2) % loop over trials
    for j = find(ismember(rmsdata.label, gcfg.gaugechannel))' % loop over gauge channels
        vardata(j,:) = rmsdata.trial{i}(j,finalGaugeFilter{i});
	end
end

switch gcfg.centerType
    case  'median'
        for j = find(ismember(rmsdata.label, gcfg.gaugechannel))' % loop over gauge channels
            gcfg.center(j) = median(vardata(j,:),2);
        end
    case  'mean'
        for j = find(ismember(rmsdata.label, gcfg.gaugechannel))' % loop over gauge channels
            gcfg.center(j) = mean(vardata(j,:),2);
        end        
end
switch gcfg.varianceType
    case 'mad'
        for j = find(ismember(rmsdata.label, gcfg.gaugechannel))' % loop over gauge channels
            gcfg.variance(j) = mad(vardata(j,:),1,2);
            gcfg.rmsampthresh(j) = gcfg.center(j) + gcfg.variance(j) .* gcfg.ampcriterion;
        end
    case 'sd'
        for j = find(ismember(rmsdata.label, gcfg.gaugechannel))' % loop over gauge channels
            gcfg.variance(j) = std(vardata(j,:),1,2);
            gcfg.rmsampthresh(j) = gcfg.center(j) + gcfg.variance(j) .* gcfg.ampcriterion;
        end
    case 'percent'
        for j = find(ismember(rmsdata.label, gcfg.gaugechannel))' % loop over gauge channels
            numelvardata(j) = numel(vardata(j,:));       
            sortvardata{j} = sort(vardata(j,:),2);
            index(j) = round(numelvardata(j)/100*(50+gcfg.ampcriterion));
            percentilevardata(j) = sortvardata{j}(index(j));
            gcfg.rmsampthresh(j) = percentilevardata(j);     

            if gcfg.findeventfreeepochs_flag
                % event-free (ef) data
                index_ef(j) = round(numelvardata(j)/100*(50-gcfg.ampcriterion));
                percentilevardata_ef(j) = sortvardata{j}(index_ef(j));
                gcfg.rmsampthresh_ef(j) = percentilevardata_ef(j);
            end
        end
    case 'times'
        for j = find(ismember(rmsdata.label, gcfg.gaugechannel))' % loop over gauge channels            
            gcfg.rmsampthresh(j) = gcfg.center(j) * gcfg.ampcriterion;            
        end
        
end

if isfield(gcfg,'loadampcritvalues') && gcfg.loadampcritvalues == 1
    
    % replace amplitude criteria by average values from previous runs
    display('Loading and averaging amplitude criteria from previously saved files.');
    rmsamthresh_loaded = []; center_loaded = [];
    for iFile = 1:length(gcfg.ampcritvalfiles)
        temp{iFile} = load(gcfg.ampcritvalfiles{iFile});
        rmsamthresh_loaded(iFile,:) = temp{iFile}.osc_gcfg.rmsampthresh;
        center_loaded(iFile,:) = temp{iFile}.osc_gcfg.center;
        
    end
    gcfg.rmsampthresh = mean(rmsamthresh_loaded,1);
    gcfg.center = mean(center_loaded,1);       
    
end % of if loadampcritvalues


switch gcfg.ampcriterionrule
    case 'channel-wise'
        % keep it as it is
    case 'uniform'
         % it is a bit problematic to average over SD-based thresholds in
         % case for more than one channel in gcfg.gaugechannel. One should
         % than rather average RMS signal over channels 
        uniformthresh = mean(gcfg.rmsampthresh(find(ismember(rmsdata.label, gcfg.gaugechannel))));
        gcfg.rmsampthresh = ones(length(gcfg.rmsampthresh),1)*uniformthresh;
end

%% find threshold crossings
% prepare filter based on search channels and search sleep stages (!! will now be identical for all channels, this can be changed in the future!!)
for i = 1:size(rmsdata.trial,2) % loop over trials
	for j = find(ismember(rmsdata.label, gcfg.searchchannel))' % loop over search channels
        channums = find(ismember(rmsdata.label, gcfg.searchchannel))'; % search channels
        finalSearchFilter{i}(j,:) = all([data.artifactFilter{i}(channums,:); gcfg.scoring(4,:)],1); % only artefact free time points in search channels and search sleep stages
    end
end

% detect crossings
for i = 1:size(rmsdata.trial,2) % loop over trials
    for j = find(ismember(rmsdata.label, gcfg.searchchannel))' % loop over search channels
        supthreshdata{i}(j,:) = all([rmsdata.trial{i}(j,:) >= gcfg.rmsampthresh(j); finalSearchFilter{i}(j,:)]); 
        if gcfg.findeventfreeepochs_flag
            supthreshdata_ef{i}(j,:) = all([rmsdata.trial{i}(j,:) <= gcfg.rmsampthresh_ef(j); finalSearchFilter{i}(j,:)]);
        end
    end
end

%% find sufficiently long suprathresold epochs and calculate metrics
gcfg.eventInfo = struct;
% for i = find(gcfg.scoring(:,4))' % loop over search epochs
for i = 1:size(rmsdata.trial,2) % loop over trials
    for j = find(ismember(rmsdata.label, gcfg.searchchannel))' % loop over search channels
        dsig = diff([0 supthreshdata{i}(j,:) 0]);
        startIndex = find(dsig > 0);
        endIndex = find(dsig < 0)-1;
        duration = endIndex-startIndex+1;
        ec = 0; % event counter
        for k = 1:length(startIndex) % loop over potential oscillations
            if duration(k) >= gcfg.durcriterion(1)*data.fsample && duration(k) <= gcfg.durcriterion(2)*data.fsample % duration criterion fullfilled
                ec = ec + 1;
                gcfg.eventInfo(i,j).stage(ec) = single(gcfg.scoring(1,startIndex(k))); % convert stage to single to have same data format in all fields (becomes imnportat later)
                gcfg.eventInfo(i,j).startTime(ec) = single(startIndex(k)); % event start (in datapoints)
                gcfg.eventInfo(i,j).midTime(ec) = single(round(mean([startIndex(k) endIndex(k)]))); % event mid (in datapoints)
                gcfg.eventInfo(i,j).endTime(ec) = single(endIndex(k)); % event end (in datapoints)
                gcfg.eventInfo(i,j).duration(ec) = single((endIndex(k)-startIndex(k))/data.fsample);  % spindel duration (in seconds)
                searchWindow = bpdata.trial{i}(j,startIndex(k):endIndex(k));
                [minAmp,minIndex] = min(searchWindow); 
                [maxAmp,maxIndex] = max(searchWindow); 
                gcfg.eventInfo(i,j).maxTime(ec) = single(startIndex(k) + maxIndex); % time of oscillation peak (in datapoints)
                gcfg.eventInfo(i,j).minTime(ec) = single(startIndex(k) + minIndex); % time of oscillation trough (in datapoints)               
                gcfg.eventInfo(i,j).minAmp(ec) = minAmp; 
                gcfg.eventInfo(i,j).maxAmp(ec) = maxAmp; 
                gcfg.eventInfo(i,j).p2pAmp(ec) = maxAmp-minAmp; 
                gcfg.eventInfo(i,j).p2pTime(ec) = single(abs(maxIndex-minIndex)/bpdata.fsample); % in seconds                
                RMSsearchWindow = rmsdata.trial{i}(j,startIndex(k):endIndex(k));               
                [RMSmaxAmp,RMSmaxIndex] = max(RMSsearchWindow);                 
                gcfg.eventInfo(i,j).RMSmaxAmp(ec) = RMSmaxAmp; % RMS max
                gcfg.eventInfo(i,j).RMSmaxTime(ec) = single(startIndex(k) + RMSmaxIndex); % time of RMS max (in datapoints)
            else % duration criterion NOT fullfilled
                supthreshdata{i}(j,startIndex(k):endIndex(k)) = 0; % remove false oscillations from suprathreshdata
            end                 
        end
    end
end


%% find sufficiently long event-free epochs and calculate metrics
if gcfg.findeventfreeepochs_flag
    gcfg.eventfreeInfo = struct;
    % for i = find(gcfg.scoring(:,4))' % loop over search epochs
    for i = 1:size(rmsdata.trial,2) % loop over trials
        for j = find(ismember(rmsdata.label, gcfg.searchchannel))' % loop over search channels
            dsig = diff([0 supthreshdata_ef{i}(j,:) 0]);
            startIndex = find(dsig > 0);
            endIndex = find(dsig < 0)-1;
            duration = endIndex-startIndex+1;
            efc = 0; % event-free counter
            for k = 1:length(startIndex) % loop over potential oscillations
                if duration(k) >= gcfg.durcriterion(1)*data.fsample && duration(k) <= gcfg.durcriterion(2)*data.fsample % duration criterion fullfilled
                    efc = efc + 1;
                    gcfg.eventfreeInfo(i,j).stage(efc) = single(gcfg.scoring(1,startIndex(k))); % convert stage to single to have same data format in all fields (becomes imnportat later)
                    gcfg.eventfreeInfo(i,j).startTime(efc) = single(startIndex(k)); % event start (in datapoints)
                    gcfg.eventfreeInfo(i,j).midTime(efc) = single(mean([startIndex(k) endIndex(k)])); % event mid (in datapoints)
                    gcfg.eventfreeInfo(i,j).endTime(efc) = single(endIndex(k)); % event end (in datapoints)
                    gcfg.eventfreeInfo(i,j).duration(efc) = single((endIndex(k)-startIndex(k))/data.fsample);  % spindel duration (in seconds)
                    searchWindow = bpdata.trial{i}(j,startIndex(k):endIndex(k));
                    [minAmp,minIndex] = min(searchWindow); 
                    [maxAmp,maxIndex] = max(searchWindow); 
                    gcfg.eventfreeInfo(i,j).maxTime(efc) = single(startIndex(k) + maxIndex); % time of oscillation peak (in datapoints)
                    gcfg.eventfreeInfo(i,j).minTime(efc) = single(startIndex(k) + minIndex); % time of oscillation trough (in datapoints)               
                    gcfg.eventfreeInfo(i,j).minAmp(efc) = minAmp; 
                    gcfg.eventfreeInfo(i,j).maxAmp(efc) = maxAmp; 
                    gcfg.eventfreeInfo(i,j).p2pAmp(efc) = maxAmp-minAmp; 
                    gcfg.eventfreeInfo(i,j).p2pTime(efc) = single(abs(maxIndex-minIndex)/bpdata.fsample); % in seconds    
                    RMSsearchWindow = rmsdata.trial{i}(j,startIndex(k):endIndex(k));               
                    [RMSmaxAmp,RMSmaxIndex] = max(RMSsearchWindow);
                    gcfg.eventfreeInfo(i,j).RMSmaxAmp(ec) = RMSmaxAmp; % in seconds
                    gcfg.eventfreeInfo(i,j).RMSmaxTime(ec) = single(startIndex(k) + RMSmaxIndex); % time of oscillation peak (in datapoints)
                else % duration criterion NOT fullfilled
                    supthreshdata_ef{i}(j,startIndex(k):endIndex(k)) = 0; % remove false oscillations from suprathreshdata
                end                 
            end
        end
    end

end


%% add summary statistics to gcfg.summary
for j = find(ismember(rmsdata.label, gcfg.gaugechannel))' % loop over gauge channels        
    gcfg.summary.median(j) = median(vardata(j,:),2);
    gcfg.summary.mean(j) = mean(vardata(j,:),2);    
    gcfg.summary.mad(j) = mad(vardata(j,:),1,2);
    gcfg.summary.std(j) = std(vardata(j,:),1,2);
    gcfg.summary.var(j) = var(vardata(j,:),1,2);    
    gcfg.summary.min(j) = min(vardata(j,:),[],2);
    gcfg.summary.max(j) = max(vardata(j,:),[],2);
    
    gcfg.summary.rmsampthresh(j) = gcfg.rmsampthresh(j);  
    
    numelvardata(j) = numel(vardata(j,:));
    sortvardata{j} = sort(vardata(j,:),2);
    
    if  ~isempty(gcfg.eventInfo(1,j).stage)
        gcfg.summary.MEAN_RMSmaxAmp(j) = mean(gcfg.eventInfo(1,j).RMSmaxAmp);
        gcfg.summary.MEDIAN_RMSmaxAmp(j) = median(gcfg.eventInfo(1,j).RMSmaxAmp);
        gcfg.summary.MIN_RMSmaxAmp(j) = min(gcfg.eventInfo(1,j).RMSmaxAmp);
        gcfg.summary.MAX_RMSmaxAmp(j) = max(gcfg.eventInfo(1,j).RMSmaxAmp);
        gcfg.summary.Perc_MEAN_RMSmaxAmp(j) = find(sortvardata{j} >= gcfg.summary.MEAN_RMSmaxAmp(j),1,'first')/numelvardata(j);
        gcfg.summary.Perc_MEDIAN_RMSmaxAmp(j) = find(sortvardata{j} >= gcfg.summary.MEDIAN_RMSmaxAmp(j),1,'first')/numelvardata(j);
        gcfg.summary.Perc_MIN_RMSmaxAmp(j) = find(sortvardata{j} >= gcfg.summary.MIN_RMSmaxAmp(j),1,'first')/numelvardata(j);
        gcfg.summary.Perc_MAX_RMSmaxAmp(j) = find(sortvardata{j} >= gcfg.summary.MAX_RMSmaxAmp(j),1,'first')/numelvardata(j);
    else
        gcfg.summary.MEAN_RMSmaxAmp(j) = NaN;
        gcfg.summary.MEDIAN_RMSmaxAmp(j) = NaN;
        gcfg.summary.MIN_RMSmaxAmp(j) = NaN;
        gcfg.summary.MAX_RMSmaxAmp(j) = NaN;
        gcfg.summary.Perc_MEAN_RMSmaxAmp(j) = NaN;
        gcfg.summary.Perc_MEDIAN_RMSmaxAmp(j) = NaN;
        gcfg.summary.Perc_MIN_RMSmaxAmp(j) =  NaN;
        gcfg.summary.Perc_MAX_RMSmaxAmp(j) =  NaN;
    end
    gcfg.summary.Perc_rmsampthresh(j) = find(sortvardata{j} >= gcfg.summary.rmsampthresh(j),1,'first')/numelvardata(j);
  
end


%% channel renaming
% rename channels
for j = 1:numel(bpdata.label) 
    bpdata.label{j} = [bpdata.label{j} '_' gcfg.eventname '_bp']; % change channel names  
end
for j = 1:numel(rmsdata.label) 
    rmsdata.label{j} = [rmsdata.label{j} '_ ' gcfg.eventname '_rms']; % change channel names  
end

%% finishing
display([gcfg.eventname ' oscillations detected.']);
ttoc = toc;
display(['Detecting ' gcfg.eventname ' events took ' num2str(ttoc) ' seconds.']);


end % of function
