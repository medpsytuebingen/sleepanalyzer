%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detectTWs
% by Til Ole Bergmann 2013
% last modified 2016/11/25 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% detects theta waves in EEG data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gcfg, data, bpdata] = sa_detectTWs(gcfg,data)
display('Detecting theta waves...');

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


%% bandpass filter in TW range
bpdata = data;
gcfg.bpfreq = repmat(gcfg.gaugefreqrange,length(gcfg.searchchannel),1);
for i = 1:size(data.trial,2) % loop over trials
    for j = find(ismember(data.label, gcfg.searchchannel))' % loop over search channels
%         bpdata.trial{i}(j,:) = ft_preproc_bandpassfilter(data.trial{i}(j,:), double(data.fsample), gcfg.bpfreq(j,:), 4, 'but', 'twopass');
        bpdata.trial{i}(j,:) = ft_preproc_bandpassfilter(data.trial{i}(j,:), double(data.fsample), gcfg.bpfreq(j,:), 3*fix(data.fsample/gcfg.bpfreq(j,1))+1, 'fir', 'twopass'); % Mathilde's formula
    end
end


%% prepare filter based on gauge channels and gauge sleep stages (!! will now be identical for all channels, this can be changed in the future!!)
for i = 1:size(bpdata.trial,2) % loop over trials
	for j = find(ismember(bpdata.label, gcfg.gaugechannel))' % loop over gauge channels 
        channums = find(ismember(bpdata.label, gcfg.gaugechannel))'; % gauge channels
        finalGaugeFilter{i}(j,:) = all([data.artifactFilter{i}(channums,:); gcfg.scoring(3,:)],1); % only artefact free time points in gauge channels and gauge sleep stages 
    end
end


%% calculate root mean square (RMS) signal (for event-free epoch detection)
rmsdata = bpdata;
rmsdata.trial = {};
for i = 1:size(bpdata.trial,2) % loop over trials
    for j = find(ismember(rmsdata.label, gcfg.searchchannel))' % loop over search channels
        tempSquared = (bpdata.trial{i}(j,:)) .^2;
        tempConvolved = conv(tempSquared,ones(1,data.fsample * gcfg.rmswindowlength),'same');
        rmsdata.trial{i}(j,:) = sqrt(tempConvolved);
    end
end

% gather data for variance threshold calculation
for i = 1:size(rmsdata.trial,2) % loop over trials
    for j = find(ismember(rmsdata.label, gcfg.gaugechannel))' % loop over gauge channels
        rmsvardata(j,:) = rmsdata.trial{i}(j,finalGaugeFilter{i}(j,:));
	end
end

%% find zero crossings
pol = [];
for i = 1:size(bpdata.trial,2) % loop over trials
    for j = find(ismember(bpdata.label, gcfg.searchchannel))' % loop over search channels
%         if ~isempty(findstr('TL',bpdata.label{j})) || ~isempty(findstr('TR',bpdata.label{j})) % if HC channel
%         if ~isempty(findstr('T',bpdata.label{j})) % if HC channel
        if isempty(ft_channelselection({'eeg'},bpdata.label(j))) % if NOT EEG and thus iEEG channel
            switch gcfg.iEEG_polarity
                case 'original'
                    pol(j) = 1; % DO NOT switch polarity     
                case 'flip'
                    pol(j) = -1; % DO switch polarity
            end            
        elseif ~isempty(ft_channelselection({'eeg'},bpdata.label(j))) % if EEG channel
            switch gcfg.EEG_polarity
                case 'original'
                    pol(j) = 1; % DO NOT switch polarity     
                case 'flip'
                    pol(j) = -1; % DO switch polarity
            end    
        end
        signVec = sign(bpdata.trial{i}(j,:) .* pol(j));
        diffVec = diff(signVec);
        pos2neg0cross{i,j} = find(diffVec < 0);
        neg2pos0cross{i,j} = find(diffVec > 0);
        if pos2neg0cross{i,j}(1) > neg2pos0cross{i,j}(1) % make sure pos2neg0cross k always preceeds neg2pos0cross k
             neg2pos0cross{i,j}(1) = [];
        end
    end
end
% plot2check
% interval = 10010000:10040000;
% figure;hold on;plot(bpdata.trial{i}(1,interval),'k');plot(signVec(interval),'b');plot(diffVec(interval),'r');hold off;


%% determine peak-to-peak amplitudes between zero crossings
for i = 1:size(bpdata.trial,2) % loop over trials
    for j = find(ismember(bpdata.label, gcfg.searchchannel))' % loop over search channels
        for k = 1:numel(pos2neg0cross{i,j})-1 
            [TWtrough{i,j}(k) TWtroughTime{i,j}(k)]= min(bpdata.trial{i}(j,pos2neg0cross{i,j}(k):neg2pos0cross{i,j}(k)) .* pol(j));
            [TWpeak{i,j}(k) TWpeakTime{i,j}(k)] = max(bpdata.trial{i}(j,neg2pos0cross{i,j}(k):pos2neg0cross{i,j}(k+1)) .* pol(j));
            TWp2p{i,j}(k)= TWpeak{i,j}(k) - TWtrough{i,j}(k);
            TWdur{i,j}(k) = (pos2neg0cross{i,j}(k+1)-pos2neg0cross{i,j}(k))/bpdata.fsample;
            TWp2pTime{i,j}(k) = ((neg2pos0cross{i,j}(k)-pos2neg0cross{i,j}(k))-TWtroughTime{i,j}(k))+(TWpeakTime{i,j}(k))/bpdata.fsample;
        end
    end
end


%% determine amplitude threshold (only for candidate TWs of correct duration!!!)

% select parameter for amplitude criterion
switch gcfg.ampcritparameter
    case 'trough'
        TWampcrit = TWtrough;
	case 'peak'
        TWampcrit = TWpeak;
    case 'p2p'        
        TWampcrit = TWp2p;
    case 'both'
        TWampcrit = TWtrough;
        TWampcrit2 = TWpeak; % will be used from here onwards as additional criterion whenever gcfg.ampcritparameter =='both' ist tested
end

% gather data for variance threshold calculation (only for candidate TWs of correct duration!!!)
for i = 1:size(bpdata.trial,2) % loop over trials
    for j = find(ismember(bpdata.label, gcfg.gaugechannel))' % loop over gauge channels      
        ec = 0;
        for k = 1:numel(TWampcrit{i,j}) % loop over TW candidates
            if all(finalGaugeFilter{i}(j,pos2neg0cross{i,j}(k):pos2neg0cross{i,j}(k+1))) % if TW candidate is located in artifact free data
                if TWdur{i,j}(k) >= gcfg.durcriterion(1) && TWdur{i,j}(k) <= gcfg.durcriterion(2) % duration criterion fullfilled
                    ec = ec+1;
                    vardata{j}(ec) = abs(TWampcrit{i,j}(k)); % use absolute (!) values of amplitude criterion variable
                    if strcmp(gcfg.ampcritparameter, 'both')
                        vardata2{j}(ec) = abs(TWampcrit2{i,j}(k)); % use absolute (!) values of 2nd amplitude criterion variable
                    end
                end
            end
        end      
        switch gcfg.centerType
            case  'median'                
                medianTWampcrit(j) = median(vardata{j},2);
                gcfg.center(j) = medianTWampcrit(j);
                if strcmp(gcfg.ampcritparameter, 'both')
                    medianTWampcrit2(j) = median(vardata2{j},2);
                    gcfg.center2(j) = medianTWampcrit2(j);
                end
            case  'mean'
                meanTWampcrit(j) = mean(vardata{j},2);
                gcfg.center(j) = meanTWampcrit(j);
                if strcmp(gcfg.ampcritparameter, 'both')
                    meanTWampcrit2(j) = mean(vardata2{j},2);
                    gcfg.center2(j) = meanTWampcrit2(j);
                end
            case  'none'
                minTWampcrit(j) = min(vardata{j},[],2);
                gcfg.center(j) = minTWampcrit(j);
                if strcmp(gcfg.ampcritparameter, 'both')
                    minTWampcrit2(j) = min(vardata2{j},2);
                    gcfg.center2(j) = minTWampcrit2(j);
                end
        end % of switch
        switch gcfg.varianceType
            case 'mad'
                madTWampcrit(j) = mad(vardata{j},1,2);
                gcfg.variance(j) = madTWampcrit(j);
                gcfg.TWampthresh(j) = gcfg.center(j) + gcfg.variance(j) .* gcfg.ampcriterion;
                if strcmp(gcfg.ampcritparameter, 'both')
                    madTWampcrit2(j) = mad(vardata2{j},1,2);
                    gcfg.variance2(j) = madTWampcrit2(j);
                    gcfg.TWampthresh2(j) = gcfg.center2(j) + gcfg.variance2(j) .* gcfg.ampcriterion;
                end
            case 'sd'
                sdTWampcrit(j) = std(vardata{j},1,2);
                gcfg.variance(j) = sdTWampcrit(j);
                gcfg.TWampthresh(j) = gcfg.center(j) + gcfg.variance(j) .* gcfg.ampcriterion;
                if strcmp(gcfg.ampcritparameter, 'both')
                    sdTWampcrit2(j) = std(vardata2{j},1,2);
                    gcfg.variance2(j) = sdTWampcrit2(j);
                    gcfg.TWampthresh2(j) = gcfg.center2(j) + gcfg.variance2(j) .* gcfg.ampcriterion;
                end
            case 'percent'
                numelvardata(j) = numel(vardata{j});
                sortvardata{j} = sort(vardata{j},2);
                index(j) = round(numelvardata(j)/100*(50+gcfg.ampcriterion));
                percentilevardata(j) = sortvardata{j}(index(j));
                gcfg.TWampthresh(j) = percentilevardata(j);
                if strcmp(gcfg.ampcritparameter, 'both')
                    numelvardata2(j) = numel(vardata2{j});
                    sortvardata2{j} = sort(vardata2{j},2);
                    index2(j) = round(numelvardata2(j)/100*(50+gcfg.ampcriterion));
                    percentilevardata2(j) = sortvardata2{j}(index2(j));
                    gcfg.TWampthresh2(j) = percentilevardata2(j);
                end
                %
                %                 % event-free (ef) data
                %                 index_ef(j) = round(numelvardata(j)/100*(50-gcfg.ampcriterion));
                %                 percentilevardata_ef(j) = sortvardata{j}(index_ef(j));
                %                 gcfg.TWampthresh_ef(j) = percentilevardata_ef(j);
        end % of switch
        
    end % of loop over gauge channels
end % of loop over trials


%% find events fullfilling duration and amplitude cirterion and calculate metrics
% prepare filter based on search channels and search sleep stages (!! will now be identical for all channels, this can be change in the future!!)
for i = 1:size(bpdata.trial,2) % loop over trials
    for j = find(ismember(bpdata.label, gcfg.searchchannel))' % loop over search channels
        channums = find(ismember(bpdata.label, gcfg.searchchannel))'; % search channels
        finalSearchFilter{i}(j,:) = all([data.artifactFilter{i}(channums,:); gcfg.scoring(4,:)],1); % only artefact free time points in search channels and search sleep stages 
    end
end

ec = 0; % reset event counter
gcfg.eventInfo = struct;
for i = 1:size(bpdata.trial,2) % loop over trials
    for j = find(ismember(bpdata.label, gcfg.searchchannel))' % loop over search channels
        for k = 1:numel(TWampcrit{i,j}) % loop over potential TWs
            if all(finalSearchFilter{i}(j, pos2neg0cross{i,j}(k):pos2neg0cross{i,j}(k+1))) % correct sleep stage and no artifact accordig to finalSearchFilter                
                if abs(TWampcrit{i,j}(k)) >= gcfg.TWampthresh(j); % amplitude criterion fullfilled
                    if strcmp(gcfg.ampcritparameter, 'both') % 2nd amplitude criterion required 
                        ok = abs(TWampcrit2{i,j}(k)) >= gcfg.TWampthresh2(j); % 2nd amplitude criterion fullfilled 
                    else
                        ok = 1; % if 2nd amplitude criterion NOT required
                    end
                    if ok % 2nd apmlitude criterion either fullfilled or not required
                        if TWdur{i,j}(k) >= gcfg.durcriterion(1) && TWdur{i,j}(k) <= gcfg.durcriterion(2) % duration criterion fullfilledif TWdur{i,j}(k) >= gcfg.durcriterion(1) && TWdur{i,j}(k) <= gcfg.durcriterion(2) % duration criterion fullfilled
                            ec = ec + 1;
                            gcfg.eventInfo(i,j).stage(ec) = single(gcfg.scoring(1,pos2neg0cross{i,j}(k))); % convert stage to single to have same data format in all fields (becomes important later)               
                            gcfg.eventInfo(i,j).startTime(ec) = single(pos2neg0cross{i,j}(k)); % event start (in datapoints)              
                            gcfg.eventInfo(i,j).midTime(ec) = single(neg2pos0cross{i,j}(k)); % event mid (in datapoints)
                            gcfg.eventInfo(i,j).endTime(ec) = single(pos2neg0cross{i,j}(k+1)); % event end (in datapoints)
                            gcfg.eventInfo(i,j).duration(ec) = single(TWdur{i,j}(k));  % event duration (in seconds)
                            gcfg.eventInfo(i,j).maxTime(ec) = single(TWpeakTime{i,j}(k) + neg2pos0cross{i,j}(k)); % time of event peak (in datapoints)
                            gcfg.eventInfo(i,j).minTime(ec) = single(TWtroughTime{i,j}(k) + pos2neg0cross{i,j}(k)); % time of event trough (in datapoints)               
                            gcfg.eventInfo(i,j).minAmp(ec) = TWtrough{i,j}(k); 
                            gcfg.eventInfo(i,j).maxAmp(ec) = TWpeak{i,j}(k); 
                            gcfg.eventInfo(i,j).p2pAmp(ec) = TWp2p{i,j}(k); 
                            gcfg.eventInfo(i,j).p2pTime(ec) = single(TWp2pTime{i,j}(k)); % in seconds
                            RMSsearchWindow = rmsdata.trial{i}(j,single(pos2neg0cross{i,j}(k)):single(pos2neg0cross{i,j}(k+1)));
                            [RMSmaxAmp,RMSmaxIndex] = max(RMSsearchWindow);
                            gcfg.eventInfo(i,j).RMSmaxAmp(ec) = RMSmaxAmp;  % RMS max
                            gcfg.eventInfo(i,j).RMSmaxTime(ec) = single(pos2neg0cross{i,j}(k) + RMSmaxIndex); % time of RMS max (in datapoints)
                        end
                    end
                end  
            end
        end
    end
end

%% find event-free epochs
if gcfg.findeventfreeepochs_flag
    winlen = gcfg.eventfreewinlen; % duration of event-free epochs in s 
    targetnumepochs = ec; % number of epochs to be generated (shall equal number of event epochs)
    thresh = min(rmsvardata,[],2); % minimal rms value to initialize adaptive threshold 
    x = {};
    gcfg.eventfreeInfo = struct;
    runthis = 1;
    pc = 0.5; % percent change variable for staircasing
    for i = 1:size(bpdata.trial,2) % loop over trials
        for j = find(ismember(rmsdata.label, gcfg.searchchannel))' % loop over search channels
            trackchanges = []; prevnumepochs = 0;
            while runthis
                x{i}(j,:) = all([rmsdata.trial{i}(j,:) <= thresh(j); finalGaugeFilter{i}],1);
                xdsig = diff([0 x{i}(j,:) 0]);
                xstartIndex = find(xdsig > 0);
                xendIndex = find(xdsig < 0); % -1 ???!?
                xduration = xendIndex-xstartIndex; % +1 ???!
                if any(floor(xduration/(winlen*rmsdata.fsample)) > 1)
                    newstarts = []; newends = []; delIndex = []; 
                    for k = 1:size(xduration,2) % loop over epochs
                        if floor(xduration(k)/(winlen*rmsdata.fsample)) > 1 % if epoch can be split in two or more 
                            times = floor(xduration(k)/(winlen*rmsdata.fsample));
                            newstarts = [newstarts ones(1,times)*xstartIndex(k) + (ones(1,times)*xduration(k)).*[1/times:1/times:1]-xduration(k)/times];
                            newends = [newends ones(1,times)*xstartIndex(k) + (ones(1,times)*xduration(k)).*[1/times:1/times:1]];
                            delIndex = [delIndex k];
                        end
                    end
                    xstartIndex(delIndex) = [];
                    xendIndex(delIndex) = [];
                    xstartIndex = unique(sort([xstartIndex newstarts]));
                    xendIndex = unique(sort([xendIndex newends]));
                    xduration = xendIndex-xstartIndex;
                end
                currnumepochs = sum(xduration >= winlen*rmsdata.fsample);
                display(['Threshold ' num2str(thresh) ' �V: ' num2str(currnumepochs) ' of ' num2str(targetnumepochs) ' event-free epochs found.']);
                trackchanges = [trackchanges currnumepochs-prevnumepochs];            
                if currnumepochs < targetnumepochs
                   thresh = thresh + thresh .* pc; 
                elseif currnumepochs > targetnumepochs
                   thresh = thresh - thresh .* pc; 
                elseif currnumepochs == targetnumepochs
                    runthis = 0;
                end
                if any(trackchanges) && numel(trackchanges) > 100 && ~any(trackchanges(end-100:end))
                    display(['Only ' num2str(currnumepochs) ' event-free epochs found to match ' num2str(targetnumepochs) ' event epochs!']);
                    break
                end
                prevnumepochs = currnumepochs;
                pc = pc-pc/10;
            end % of while loop
            whichepochs = xduration >= winlen*rmsdata.fsample;
            epochs = [xstartIndex(whichepochs)' xendIndex(whichepochs)']; % get epoch start and end
            epochcenter = round(mean(epochs,2))'; % derive center of epoch as row vector
            gcfg.eventfreeInfo(i,j).stage = single(gcfg.scoring(1,epochcenter)); % convert stage to single to have same data format in all fields (becomes important later)               
            gcfg.eventfreeInfo(i,j).midTime = single(epochcenter); % mid of event-free epoch (in datapoints)
            gcfg.eventfreeInfo(i,j).maxTime = single(epochcenter); % dummy
            gcfg.eventfreeInfo(i,j).minTime = single(epochcenter); % dummy
        end
    end
end


%% add summary statistics to gcfg.summary
for j = find(ismember(bpdata.label, gcfg.gaugechannel))' % loop over gauge channels        
    gcfg.summary.median(j) = median(vardata{j},2);
    gcfg.summary.mean(j) = mean(vardata{j},2);    
    gcfg.summary.mad(j) = mad(vardata{j},1,2);
    gcfg.summary.std(j) = std(vardata{j},1,2);
    gcfg.summary.var(j) = var(vardata{j},1,2);    
    gcfg.summary.min(j) = min(vardata{j},[],2);
    gcfg.summary.max(j) = max(vardata{j},[],2);
    
    gcfg.summary.MEAN_TWp2p = mean(gcfg.eventInfo(1,j).p2pAmp);
    gcfg.summary.MEDIAN_TWp2p = median(gcfg.eventInfo(1,j).p2pAmp);
    gcfg.summary.MIN_TWp2p = min(gcfg.eventInfo(1,j).p2pAmp);
    gcfg.summary.MAX_TWp2p = max(gcfg.eventInfo(1,j).p2pAmp);
    
    gcfg.summary.TWampthresh(j) = gcfg.TWampthresh(j);
    
    numelvardata(j) = numel(vardata{j});
    sortvardata{j} = sort(vardata{j},2);
    gcfg.summary.Perc_MEAN_TWp2p(j) = find(sortvardata{j} >= gcfg.summary.MEAN_TWp2p,1,'first')/numelvardata(j);
    gcfg.summary.Perc_MEDIAN_TWp2p(j) = find(sortvardata{j} >= gcfg.summary.MEDIAN_TWp2p,1,'first')/numelvardata(j);
    gcfg.summary.Perc_MIN_TWp2p(j) = find(sortvardata{j} >= gcfg.summary.MIN_TWp2p,1,'first')/numelvardata(j);
    gcfg.summary.Perc_MAX_TWp2p(j) = find(sortvardata{j} >= gcfg.summary.MAX_TWp2p,1,'first')/numelvardata(j);
    gcfg.summary.Perc_TWampthresh(j) = find(sortvardata{j} >= gcfg.summary.TWampthresh(j),1,'first')/numelvardata(j);
 
end


%% channel renaming
% rename channels
for j = 1:numel(bpdata.label) 
    bpdata.label{j} = [bpdata.label{j} '_TW_bp']; % change channel names  
end
data.cfg = gcfg;
end % of function
