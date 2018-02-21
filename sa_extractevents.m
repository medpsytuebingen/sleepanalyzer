%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extractoscilaltions
% by Til Ole Bergmann 2016
% last modified 2016/12/09 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% extracts data around detected events from EEG data
%
% [gcfg, data] = sa_extractevents(cfg,data)
%
%  cfg.eventname = string indicating the name of the event
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gcfg, eventdata] = sa_extractevents(gcfg,data)
tic;

%% check config
if ~isfield(gcfg,'eventname')
    gcfg.eventname = 'event';
end

%% preparation
display(['Extracting ' gcfg.eventname ' events...']);

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
%     data.trial{i}(chanfilter,:) = ft_preproc_highpassfilter(data.trial{i}(chanfilter,:), double(data.fsample), gcfg.hpfilter, 3*fix(data.fsample/0.16), 'fir', 'twopass'); 
    data.trial{i}(chanfilter,:) = ft_preproc_highpassfilter(data.trial{i}(chanfilter,:), double(data.fsample), gcfg.hpfilter, 6, 'but', 'twopass','reduce'); 
end

eventdata.fsample = data.fsample;
eventdata.label = data.label(find(ismember(data.label, gcfg.extractchannel)));
tlchan = find(ismember(data.label, gcfg.timelockchannel))'; % timelock channel number
ec = 0; % event counter

for i = 1:size(gcfg.eventInfo,1) % loop over trials 
    for k = 1:length(gcfg.eventInfo(i,tlchan).maxTime) % loop over events 
        if ismember(gcfg.eventInfo(i,tlchan).stage(k),gcfg.extractstages); % if event falls in correct sleep stage
            switch gcfg.timelockevent
                case 'peak'
                    eventwindow = [int32(gcfg.eventInfo(i,tlchan).maxTime(k)) + int32(gcfg.extractwindow(1)*data.fsample), int32(gcfg.eventInfo(i,tlchan).maxTime(k)) + int32(gcfg.extractwindow(2)*data.fsample)];
                    artifactwindow = [int32(gcfg.eventInfo(i,tlchan).maxTime(k)) + int32(gcfg.artifactfreewindow(1)*data.fsample), int32(gcfg.eventInfo(i,tlchan).maxTime(k)) + int32(gcfg.artifactfreewindow(2)*data.fsample)];
                case 'trough'
                    eventwindow = [int32(gcfg.eventInfo(i,tlchan).minTime(k)) + int32(gcfg.extractwindow(1)*data.fsample), int32(gcfg.eventInfo(i,tlchan).minTime(k)) + int32(gcfg.extractwindow(2)*data.fsample)];
                    artifactwindow = [int32(gcfg.eventInfo(i,tlchan).minTime(k)) + int32(gcfg.artifactfreewindow(1)*data.fsample), int32(gcfg.eventInfo(i,tlchan).minTime(k)) + int32(gcfg.artifactfreewindow(2)*data.fsample)];
            end
            
            % check whether to reject windows with artifacts
            if gcfg.rejectartifacts
                artifactfree = all(all(data.artifactFilter{i}(ismember(data.label, gcfg.extractchannel),max([artifactwindow(1) 1]):min(artifactwindow(2),length(data.artifactFilter{i})))),2);
            else 
                artifactfree = 1;
            end
            
            % if the extraction window does not include a trial border or artifacts in any to be extracted channel
            if eventwindow(1) >= 0 && eventwindow(2) <= size(data.trial{i},2) && artifactfree
                
                ec = ec+1;
                eventdata.trial{ec} = data.trial{i}(:,[eventwindow(1):eventwindow(2)-1]);
                eventdata.time{ec} = [gcfg.extractwindow(1):1/eventdata.fsample:gcfg.extractwindow(2)-(1/eventdata.fsample)];

                % 'stage' 'startTime' 'midTime' 'endTime' 'duration' 'maxTime' 'minTime' 'minAmp' 'maxAmp' 'p2pAmp' 'p2pTime' 'RMSmaxAmp' 'RMSmaxTime'
                temp = cell2mat(struct2cell(gcfg.eventInfo(i,tlchan)))';
                eventdata.trialinfo(ec,:) = temp(k,:); 
            end
        end
    end
end

%% finishing
eventdata.cfg = gcfg;
ttoc = toc;
display(['Extracting ' gcfg.eventname ' events took ' num2str(ttoc) ' seconds.']);

end % of function