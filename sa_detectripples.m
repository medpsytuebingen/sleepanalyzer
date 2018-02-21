%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detectripples
% by Til Ole Bergmann 2013
% last modified 2016/11/22 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% detects ripples in iEEG data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gcfg, data, bpdata, rmsdata, supthreshdata] = sa_detectripples(gcfg,data)
display('Detecting ripples...');

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
if gcfg.loadartifactfilter
%     temp = load(fullfile(gcfg.datapath,[gcfg.subjectName '_artifactFilter.mat']));
%     temp = load(fullfile(gcfg.datapath,'artifactFilter.mat'));
%     temp = load(fullfile(gcfg.datapath,['artifactFilter_' [gcfg.artifactfreechannels{:}] '_12sec_equated_wake_nREM_REM_ar.mat']));
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
else
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
end


%% bandpass filter in ripple range
bpdata = data;
% gcfg.bpfreq = [gcfg.ripplePeakFreq' - repmat(gcfg.searchfreqmargin,length(gcfg.ripplePeakFreq),1) gcfg.ripplePeakFreq' +  repmat(gcfg.searchfreqmargin,length(gcfg.ripplePeakFreq),1)];
gcfg.bpfreq = repmat(gcfg.gaugefreqrange,length(gcfg.searchchannel),1);

for i = 1:size(data.trial,2) % loop over trials
    for j = find(ismember(data.label, gcfg.searchchannel))' % loop over search channels
%         bpdata.trial{i}(j,:) = ft_preproc_bandpassfilter(data.trial{i}(j,:), data.fsample, gcfg.bpfreq(j,:), 4, 'but', 'twopass', 'reduce');
        bpdata.trial{i}(j,:) = ft_preproc_bandpassfilter(data.trial{i}(j,:), double(data.fsample), gcfg.bpfreq(j,:), 3*fix(data.fsample/gcfg.bpfreq(j,1))+1, 'fir', 'twopass'); % Mathilde's formula
    end
end

%% calculate root mean square (RMS) signal 
rmsdata = bpdata;
rmsdata.trial = {};
for i = 1:size(bpdata.trial,2) % loop over trials
    for j = find(ismember(rmsdata.label, gcfg.searchchannel))' % loop over search channels
        tempSquared = (bpdata.trial{i}(j,:)) .^2;
        tempConvolved = conv(tempSquared,ones(1,data.fsample * gcfg.rmswindowlength),'same');
        rmsdata.trial{i}(j,:) = sqrt(tempConvolved);
    end
end

%% determine amplitude threshold
% prepare filter
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
            gcfg.center = median(vardata,2);
        end
    case  'mean'
        for j = find(ismember(rmsdata.label, gcfg.gaugechannel))' % loop over gauge channels
            gcfg.center = mean(vardata,2);
        end
end
switch gcfg.varianceType
    case 'mad'
        for j = find(ismember(rmsdata.label, gcfg.gaugechannel))' % loop over gauge channels
            gcfg.variance = mad(vardata,1,2);
            gcfg.rmsampthresh = gcfg.center + gcfg.variance .* gcfg.ampcriterion;
            if isfield(gcfg,'gcfg.amplimit')
                gcfg.limrmsampthresh(j) = gcfg.center + gcfg.variance .* gcfg.amplimit;
            else
                gcfg.limrmsampthresh(j) = max(vardata(j,:),[],2);
            end
        end
    case 'sd'
        for j = find(ismember(rmsdata.label, gcfg.gaugechannel))' % loop over gauge channels
            gcfg.variance = std(vardata,1,2);
            gcfg.rmsampthresh = gcfg.center + gcfg.variance .* gcfg.ampcriterion;
            if isfield(gcfg,'gcfg.amplimit')
                gcfg.limrmsampthresh(j) = gcfg.center + gcfg.variance .* gcfg.amplimit;
            else
                gcfg.limrmsampthresh(j) = max(vardata(j,:),[],2);
            end
        end
    case 'percent'
        for j = find(ismember(rmsdata.label, gcfg.gaugechannel))' % loop over gauge channels
            for j = find(ismember(rmsdata.label, gcfg.gaugechannel))' % loop over gauge channels
                numelvardata(j) = numel(vardata(j,:));
                sortvardata{j} = sort(vardata(j,:),2);
                index(j) = round(numelvardata(j)/100*(50+gcfg.ampcriterion));
                gcfg.rmsampthresh(j) = sortvardata{j}(index(j));
                if isfield(gcfg,'gcfg.amplimit')
                    limindex(j) = round(numelvardata(j)/100*(50+gcfg.amplimit));
                    gcfg.limrmsampthresh(j) = sortvardata{j}(limindex(j));
                else
                    gcfg.limrmsampthresh(j) = max(vardata(j,:),[],2);
                end
            end
        end
end
    
switch gcfg.ampcriterionrule
    case 'channel-wise'
        % keep it as it is
    case 'uniform'
         % it is a bit problematic to average over SD-based thresholds in
         % case fo more than one channel in gcfg.gaugechannel. One should
         % than rather average RMS signal over channels 
        uniformthresh = mean(gcfg.rmsampthresh(find(ismember(rmsdata.label, gcfg.gaugechannel))));
        gcfg.rmsampthresh = ones(length(gcfg.rmsampthresh),1)*uniformthresh;
end

%% find threshold crossings
% prepare filter based on search channels and search sleep stages (!! will now be identical for all channels, this can be change in the future!!)
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
        for k = 1:length(startIndex) % loop over potential ripples
            if duration(k) >= gcfg.durcriterion(1)*data.fsample && duration(k) <= gcfg.durcriterion(2)*data.fsample % duration criterion fullfilled
                if max(rmsdata.trial{i}(j,startIndex(k):endIndex(k))) <= gcfg.limrmsampthresh(j) % limit RMS amplitude criterion not exceeded
                    
                    ec = ec + 1;
                    gcfg.eventInfo(i,j).stage(ec) = single(gcfg.scoring(1,startIndex(k))); % convert stage to single to have same data format in all fields (becomes imnportat later)
                    gcfg.eventInfo(i,j).startTime(ec) = single(startIndex(k)); % event start (in datapoints)
                    gcfg.eventInfo(i,j).midTime(ec) = single(mean([startIndex(k) endIndex(k)])); % event mid (in datapoints)
                    gcfg.eventInfo(i,j).endTime(ec) = single(endIndex(k)); % event end (in datapoints)
                    gcfg.eventInfo(i,j).duration(ec) = single((endIndex(k)-startIndex(k))/data.fsample);  % ripple duration (in seconds)
                    searchWindow = bpdata.trial{i}(j,startIndex(k):endIndex(k));
                    [minAmp,minIndex] = min(searchWindow);
                    [maxAmp,maxIndex] = max(searchWindow);
                    gcfg.eventInfo(i,j).maxTime(ec) = single(startIndex(k) + maxIndex); % time of ripple peak (in datapoints)
                    gcfg.eventInfo(i,j).minTime(ec) = single(startIndex(k) + minIndex); % time of ripple trough (in datapoints)
                    gcfg.eventInfo(i,j).minAmp(ec) = minAmp;
                    gcfg.eventInfo(i,j).maxAmp(ec) = maxAmp;
                    gcfg.eventInfo(i,j).p2pAmp(ec) = maxAmp-minAmp;
                    gcfg.eventInfo(i,j).p2pTime(ec) = single(abs(maxIndex-minIndex)/data.fsample); % in seconds
                    RMSsearchWindow = rmsdata.trial{i}(j,startIndex(k):endIndex(k));
                    [RMSmaxAmp,RMSmaxIndex] = max(RMSsearchWindow);
                    gcfg.eventInfo(i,j).RMSmaxAmp(ec) = RMSmaxAmp; % RMS max
                    gcfg.eventInfo(i,j).RMSmaxTime(ec) = single(startIndex(k) + RMSmaxIndex); % time of RMS max (in datapoints)

                else % limit RMS criterion exceeded
                    supthreshdata{i}(j,startIndex(k):endIndex(k)) = 0; % remove false ripples from suprathreshdata                    
                end
                
            else % duration criterion NOT fullfilled
                supthreshdata{i}(j,startIndex(k):endIndex(k)) = 0; % remove false ripples from suprathreshdata
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
    
    gcfg.summary.MEAN_RMSmaxAmp = mean(gcfg.eventInfo(1,j).RMSmaxAmp);
    gcfg.summary.MEDIAN_RMSmaxAmp = median(gcfg.eventInfo(1,j).RMSmaxAmp);
    gcfg.summary.MIN_RMSmaxAmp = min(gcfg.eventInfo(1,j).RMSmaxAmp);
    gcfg.summary.MAX_RMSmaxAmp = max(gcfg.eventInfo(1,j).RMSmaxAmp);
    
    gcfg.summary.rmsampthresh(j) = gcfg.rmsampthresh(j);
    
    numelvardata(j) = numel(vardata(j,:));
    sortvardata{j} = sort(vardata(j,:),2);
    gcfg.summary.Perc_MEAN_RMSmaxAmp(j) = find(sortvardata{j} >= gcfg.summary.MEAN_RMSmaxAmp,1,'first')/numelvardata(j);
    gcfg.summary.Perc_MEDIAN_RMSmaxAmp(j) = find(sortvardata{j} >= gcfg.summary.MEDIAN_RMSmaxAmp,1,'first')/numelvardata(j);
    gcfg.summary.Perc_MIN_RMSmaxAmp(j) = find(sortvardata{j} >= gcfg.summary.MIN_RMSmaxAmp,1,'first')/numelvardata(j);
    gcfg.summary.Perc_MAX_RMSmaxAmp(j) = find(sortvardata{j} >= gcfg.summary.MAX_RMSmaxAmp,1,'first')/numelvardata(j);
    gcfg.summary.Perc_rmsampthresh(j) = find(sortvardata{j} >= gcfg.summary.rmsampthresh(j),1,'first')/numelvardata(j);
 
end


%% channel renaming
% rename channels
for j = 1:numel(bpdata.label) 
    bpdata.label{j} = [bpdata.label{j} '_ripple_bp']; % change channel names  
end
for j = 1:numel(rmsdata.label) 
    rmsdata.label{j} = [rmsdata.label{j} '_ripple_rms']; % change channel names  
end


end % of function
