%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extractTWs
% by Til Ole Bergmann 2013
% last modified 2016/11/28 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% extracts data around theta waves from EEG data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gcfg, TWdata] = sa_extractTWs(gcfg,data)
display('Extracting theta waves...');

% define channels to be extracted
if strcmp(gcfg.extractchannel,'all')
    gcfg.extractchannel = data.label;
end

% define channels to be filtered
procchannels = {'_bp','_rms'};
goodchan =[];
for i = 1:length(data.label)
    for j = 1:length(procchannels) % loop over channels postfixes to exclude
        k = strfind(data.label{i},procchannels{j});  
        goodchan(i,j) = isempty(k);
    end
end
goodchan = logical(prod(goodchan,2));
gcfg.filterchannel = data.label(logical(goodchan));

% high pass filter data before extraction
for i = size(data.trial,2) % loop over trials
    chanfilter = find(ismember(data.label, gcfg.filterchannel))';
    data.trial{i}(chanfilter,:) = ft_preproc_highpassfilter(data.trial{i}(chanfilter,:), double(data.fsample), gcfg.hpfilter, 3*fix(data.fsample/0.16)+1, 'fir', 'twopass'); % Mathilde's formula
end

TWdata.fsample = data.fsample;
TWdata.label = data.label(find(ismember(data.label, gcfg.extractchannel)));
tlchan = find(ismember(data.label, gcfg.timelockchannel))'; % timelock channel number
ec = 0; % event counter
for i = 1:size(gcfg.eventInfo,1) % loop over trials 
    for k = 1:length(gcfg.eventInfo(i,tlchan).maxTime) % loop over TWs 
        if ismember(gcfg.eventInfo(i,tlchan).stage(k),gcfg.extractstages); % if event falls in correct sleep stage
            switch gcfg.timelockevent
                case 'peak'
                    TWwindow = [int32(gcfg.eventInfo(i,tlchan).maxTime(k)) + int32(gcfg.extractwindow(1)*data.fsample), int32(gcfg.eventInfo(i,tlchan).maxTime(k)) + int32(gcfg.extractwindow(2)*data.fsample)];
                    artifactwindow = [int32(gcfg.eventInfo(i,tlchan).maxTime(k)) + int32(gcfg.artifactfreewindow(1)*data.fsample), int32(gcfg.eventInfo(i,tlchan).maxTime(k)) + int32(gcfg.artifactfreewindow(2)*data.fsample)];
                case 'trough'
                    TWwindow = [int32(gcfg.eventInfo(i,tlchan).minTime(k)) + int32(gcfg.extractwindow(1)*data.fsample), int32(gcfg.eventInfo(i,tlchan).minTime(k)) + int32(gcfg.extractwindow(2)*data.fsample)];
                    artifactwindow = [int32(gcfg.eventInfo(i,tlchan).minTime(k)) + int32(gcfg.artifactfreewindow(1)*data.fsample), int32(gcfg.eventInfo(i,tlchan).minTime(k)) + int32(gcfg.artifactfreewindow(2)*data.fsample)];
                case 'center'
                    TWwindow = [int32(gcfg.eventInfo(i,tlchan).midTime(k)) + int32(gcfg.extractwindow(1)*data.fsample), int32(gcfg.eventInfo(i,tlchan).midTime(k)) + int32(gcfg.extractwindow(2)*data.fsample)];                    
                    artifactwindow = [int32(gcfg.eventInfo(i,tlchan).midTime(k)) + int32(gcfg.artifactfreewindow(1)*data.fsample), int32(gcfg.eventInfo(i,tlchan).midTime(k)) + int32(gcfg.artifactfreewindow(2)*data.fsample)];                    
            end
            
            % check whether to reject windows with artifacts
            if gcfg.rejectartifacts
                artifactfree = all(all(data.artifactFilter{i}(ismember(data.label, gcfg.extractchannel),max([artifactwindow(1) 1]):min(artifactwindow(2),length(data.artifactFilter{i})))),2);
            else 
                artifactfree = 1;
            end
            
            % if the extraction window does not include a trial border or artifacts in any to be extracted channel
            if TWwindow(1) >= 0 && TWwindow(2) <= size(data.trial{i},2) && artifactfree
            
                ec = ec+1;
                TWdata.trial{ec} = data.trial{i}(:,[TWwindow(1):TWwindow(2)-1]);
                TWdata.time{ec} = [gcfg.extractwindow(1):1/TWdata.fsample:gcfg.extractwindow(2)-(1/TWdata.fsample)];

                % 'stage' 'startTime' 'midTime' 'endTime' 'duration' 'maxTime' 'minTime' 'minAmp' 'maxAmp' 'p2pAmp' 'p2pTime' 'RMSmaxAmp' 'RMSmaxTime'
                temp = cell2mat(struct2cell(gcfg.eventInfo(i,tlchan)))';
                TWdata.trialinfo(ec,:) = temp(k,:); 
            end
        end
    end
end
eventdata.cfg = gcfg;
end % of function