%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualeventvalidation
% by Til Ole Bergmann 2015 (based on a script by Bernhard Staresina)
% last modified 2016/11/22 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% allows to visually validate events that were detected automatically
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sa_visualeventvalidation(gcfg)

switch gcfg.targetchan
    case 'Cz'
        gcfg.targetchannel = gcfg.Cztargetchannel; 
        gcfg.spindle_timelockevent = 'peak';
    case 'HC'
        gcfg.targetchannel = gcfg.HCtargetchannel; 
        gcfg.spindle_timelockevent = 'trough';
end


%% load individual data


gcfg.subjectname = sprintf('%s%0.2d','S',gcfg.subjectno);

% load raw data
load(fullfile(gcfg.paths.data,gcfg.subjectname,['dataone_' gcfg.targetchannel{gcfg.subjectno} '.mat']));

% load artifact filter
load(fullfile(gcfg.paths.data,gcfg.subjectname,['artifactFilter_' gcfg.targetchannel{gcfg.subjectno} '_3sec_unequated_wake_nREM_REM_ar.mat']));
artbeg = find(diff([1 artifactFilter]) == -1); 
artend = find(diff([artifactFilter 1]) == 1); 

% load SO events
load(fullfile(gcfg.paths.results,gcfg.subjectname,[gcfg.subjectname '_SO_[0.50-1.25Hz]_trough_' gcfg.targetchannel{gcfg.subjectno} '_eventdata.mat']));
events_SO = [SO_eventdata.trialinfo(:,2) SO_eventdata.trialinfo(:,4)];

% load spindle events
load(fullfile(gcfg.paths.results,gcfg.subjectname,[gcfg.subjectname '_spindle_' gcfg.spindle_timelockevent '_' gcfg.targetchannel{gcfg.subjectno} '_eventdata.mat']));
events_spindle = [spindle_eventdata.trialinfo(:,2) spindle_eventdata.trialinfo(:,4)];

% select stages
% badstagefilt = ~ismember(data.scoring,[2 3 4]);
% badstagebeg = find(diff([0 badstagefilt]) == 1); 
% badstageend = find(diff([badstagefilt 0]) == -1); 
% badstage = [badstagebeg' badstageend'];

SO_goodstagefilt = ismember(data.scoring,[2 3 4]);
SO_goodstagebeg = find(diff([0 SO_goodstagefilt]) == 1); 
SO_goodstageend = find(diff([SO_goodstagefilt 0]) == -1);
SO_goodstage = [SO_goodstagebeg' SO_goodstageend'];
SO_comb_goodstage_artifact_filt = SO_goodstagefilt .* artifactFilter;
SO_badcombbeg = find(diff([1 SO_comb_goodstage_artifact_filt]) == -1); 
SO_badcombend = find(diff([SO_comb_goodstage_artifact_filt 1]) == 1); 
SO_badcomb = [SO_badcombbeg' SO_badcombend'];
SO_badcomb(2:end,1) = SO_badcomb(2:end,1) - 2500; % extra artifact padding because of extract artifact padding 
SO_badcomb(1:end-1,2) = SO_badcomb(1:end-1,2) + 2500; % extra artifact padding because of extract artifact padding 

spindle_goodstagefilt = ismember(data.scoring,[2 3 4]);
spindle_goodstagebeg = find(diff([0 spindle_goodstagefilt]) == 1); 
spindle_goodstageend = find(diff([spindle_goodstagefilt 0]) == -1);
spindle_goodstage = [spindle_goodstagebeg' spindle_goodstageend'];
spindle_comb_goodstage_artifact_filt = spindle_goodstagefilt .* artifactFilter;
spindle_badcombbeg = find(diff([1 spindle_comb_goodstage_artifact_filt]) == -1); 
spindle_badcombend = find(diff([spindle_comb_goodstage_artifact_filt 1]) == 1); 
spindle_badcomb = [spindle_badcombbeg' spindle_badcombend'];
spindle_badcomb(2:end,1) = spindle_badcomb(2:end,1) - 2500; % extra artifact padding because of extract artifact padding 
spindle_badcomb(1:end-1,2) = spindle_badcomb(1:end-1,2) + 2500; % extra artifact padding because of extract artifact padding 


%% ft_databrowser
% 
% 
% % SO filtered data
% cfg = [];
% cfg.bpfilter    = 'yes';
% cfg.bpfreq  = [0.16 1.25];
% cfg.bpfiltord       = 3*fix(data.fsample/cfg.bpfreq(1))+1;
% cfg.bpfilttype      = 'fir';
% SOfilt_data = ft_preprocessing(cfg,data);
% SOfilt_data.label = {[gcfg.targetchan '_SOfilt']};
% 
% % spindle filtered data
% cfg = [];
% cfg.bpfilter    = 'yes';
% cfg.bpfreq  = [12 16];
% cfg.bpfiltord       = 3*fix(data.fsample/cfg.bpfreq(1))+1;
% cfg.bpfilttype      = 'fir';
% spindlefilt_data = ft_preprocessing(cfg,data);
% spindlefilt_data.label = {[gcfg.targetchan '_spindlefilt']};
% 
% %a append channels
% comb_data = ft_appenddata([],data, SOfilt_data, spindlefilt_data);
% 
% cfg = [];
% cfg.bpfilter    = 'yes';
% cfg.bpfreq  = [0.3 35];
% cfg.bpfiltord       = 3*fix(data.fsample/cfg.bpfreq(1))+1;
% cfg.bpfilttype      = 'fir';
% filt_data = ft_preprocessing(cfg,data);

% run databrowser
% cfg = [];
% cfg.channel                 = {gcfg.targetchan, [gcfg.targetchan '_SOfilt'], [gcfg.targetchan '_spindlefilt']};
% cfg.channel                 = {gcfg.targetchan, [gcfg.targetchan '_SOfilt'], [gcfg.targetchan '_spindlefilt']};
% cfg.preproc.demean          = 'no';
% cfg.preproc.hpfilter        = 'no';
% cfg.preproc.bpfilter        = 'no';
% cfg.preproc.hpfreq          = 250;
% cfg.preproc.bpfreq          = [0.3 35]; % AASM recomendations (Iber et al., 2007)
% cfg.bpfiltord               = 3*fix(data.fsample/cfg.preproc.bpfreq(1))+1;
% cfg.bpfilttype              = 'fir';

% cfg.ylim                    = [-150 150]; % Cz 75�V SO ampcrit
% % cfg.ylim                    = [-600 600]; % HC 300�V SO ampcrit
% cfg.blocksize               = 10; % 30
% cfg.continuous              = 'yes';
% cfg.viewmode                =  'vertical' ;  % 'butterfly' 'vertical'
% % cfg.chanscale               = [1 2 1.5];

% cfg.artfctdef.auto.artifact     = [artbeg' artend'];
% cfg.artfctdef.SO.artifact       = events_SO;
% cfg.artfctdef.spindle.artifact  = events_spindle;



%% filter data
cfg = [];
cfg.bpfilter    = 'yes';
cfg.bpfreq  = [0.5 35];
cfg.bpfiltord       = 3; % 4th order butterworth filter is instable!
cfg.bpfilttype      = 'but';
filt_data = ft_preprocessing(cfg,data);

keyboard;


% check out hypnogramme in 10 s steps
hypnogramme = downsample(data.scoring,10*1000);
plot(hypnogramme);


%% SOs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SO validation

cfg = [];
switch gcfg.targetchan
    case 'Cz'
        cfg.ylim                    = [-150 150]; % Cz 75�V SO ampcrit
    case 'HC'
        cfg.ylim                    = [-600 600]; % HC 300�V SO ampcrit
end
cfg.blocksize               = 10;
cfg.continuous              = 'yes';
cfg.artfctdef.auto.artifact     = SO_badcomb; 
cfg.artfctdef.SO.artifact       = [];
SOx = ft_databrowser(cfg,filt_data);

keyboard;

% SO RE-validation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cfg.artfctdef.auto.artifact     = SO_badcomb; 
cfg.artfctdef.SO.artifact       = SOx.artfctdef.SO.artifact; % !!!!!!!!
cfg.artfctdef.hit_SO.artifact  = SOhit_vec;
cfg.artfctdef.a_SO.artifact = events_SO; 
cfg.artfctdef.miss_SO.artifact  = SOmiss_vec;
cfg.artfctdef.false_SO.artifact  = SOfalse_vec;
cfg.artfctdef.aselect_SO.artifact = events_SO_select; 
SOx = ft_databrowser(cfg,filt_data);
keyboard;

SOmtimewin = [3690 7780] *1000; % provide analyzed data epoch here


%% SO evaluate match

% select automatic SOs from manually scored epoch only
events_SO_select = events_SO(events_SO(:,1)>SOmtimewin(1) & events_SO(:,2)<SOmtimewin(2),:);

% remove automatic SOs from artifact and bad stage epochs
SO_badcomb_vec = [];
for i = 1:length(SO_badcomb)
    SO_badcomb_vec = [SO_badcomb_vec SO_badcomb(i,1):SO_badcomb(i,2)];
end
remove_idx = [];
for i = 1:size(events_SO_select,1)
    if any(ismember([events_SO_select(i,1):events_SO_select(i,2)], SO_badcomb_vec))
        remove_idx = [remove_idx i];
    end
end
events_SO_select(remove_idx,:) = []; 

% construct continous event vectors
mSO = [];
for i = 1:length(SOx.artfctdef.SO.artifact)
    mSO= [mSO SOx.artfctdef.SO.artifact(i,1):SOx.artfctdef.SO.artifact(i,2)];
end
aSO = [];
for i = 1:length(events_SO_select)
    aSO = [aSO events_SO_select(i,1):events_SO_select(i,2)];
end

% check matching
mSOtab = []; SOhit = 0; SOfalse = 0; SOmiss = 0; SOhit_vec = []; SOmiss_vec = [];
for i = 1:length(SOx.artfctdef.SO.artifact)   
    if any(ismember([SOx.artfctdef.SO.artifact(i,1):SOx.artfctdef.SO.artifact(i,2)], aSO))
        mSOtab(i) = 1;
        SOhit = SOhit +1;
        % find automatic event overlapping with the hit manual event
        hitmevent = events_SO_select(SOx.artfctdef.SO.artifact(i,1) <= events_SO_select(:,2) & SOx.artfctdef.SO.artifact(i,2) >= events_SO_select(:,1),:);
        [mindist mindistloc] = min(mean(hitmevent,2) - mean(SOx.artfctdef.SO.artifact(i,:))); % the center of which of multiple overlapping events is closest
        SOhit_vec(SOhit,:) = hitmevent(mindistloc,:); % use longest event in case of multiple hit events              
%         [maxdur, maxdurloc] = max(hitmevent(:,2)-hitmevent(:,1)); % which event is of longest  duration for multiple hits
%         SOhit_vec(SOhit,:) = hitmevent(maxdurloc,:); % use longest event in case of multiple hit events
%         SOhit_vec(SOhit,:) = hitmevent(1,:); % use first event in case of multiple hit events
    else
        mSOtab(i) = 0;
        SOmiss = SOmiss +1;
        SOmiss_vec(SOmiss,:) = SOx.artfctdef.SO.artifact(i,:);
    end
end
aSOtab = []; SOfalse_vec = [];
for i = 1:length(events_SO_select)  
    if any(ismember([events_SO_select(i,1):events_SO_select(i,2)], mSO))
        aSOtab(i) = 1;
    else
        aSOtab(i) = 0;
        SOfalse = SOfalse +1;
        SOfalse_vec(SOfalse,:) = events_SO_select(i,:);
    end
end


%% SO display results
aSOtotal = length(events_SO_select);
mSOtotal = length(SOx.artfctdef.SO.artifact);
display(' ');
display(['SO VALIDATION for ' gcfg.subjectname ' in ' gcfg.targetchan]);
display(['SO epoch: [' num2str(SOmtimewin(1)) ' ' num2str(SOmtimewin(2)) ']']);
display(['SO number automatic/manual: ' num2str(aSOtotal) '/' num2str(mSOtotal)]);
display(['SO hits: ' num2str(SOhit) ' (' num2str(round(SOhit/mSOtotal*100)) '%)']);
display(['SO misses: ' num2str(SOmiss) ' (' num2str(round(SOmiss/mSOtotal*100)) '%)']);
display(['SO false alarms: ' num2str(SOfalse)]);
save(fullfile(gcfg.paths.results,gcfg.subjectname,[gcfg.subjectname '_SO_validation_' gcfg.targetchannel{gcfg.subjectno} '.mat']), 'mSOtotal', 'aSOtotal', 'SOhit','SOmiss','SOfalse','SOfalse_vec','SOhit_vec','SOmiss_vec','SOx','SOmtimewin');





%% SPINDLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% spindle validation

cfg = [];
switch gcfg.targetchan
    case 'Cz'
        cfg.ylim                    = [-150 150]; % Cz 75�V SO ampcrit
    case 'HC'
        cfg.ylim                    = [-600 600]; % HC 300�V SO ampcrit
end
cfg.blocksize               = 10;
cfg.continuous              = 'yes';
cfg.artfctdef.auto.artifact     = spindle_badcomb; 
cfg.artfctdef.spindle.artifact  = [];
spindlex = ft_databrowser(cfg,filt_data);

keyboard;

% spindle RE-validation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cfg.artfctdef.auto.artifact     = spindle_badcomb; 
cfg.artfctdef.spindle.artifact  = spindlex.artfctdef.spindle.artifact; % !!
cfg.artfctdef.hit_spindle.artifact  = spindlehit_vec;
cfg.artfctdef.a_spindle.artifact  = events_spindle; 
cfg.artfctdef.miss_spindle.artifact  = spindlemiss_vec;
cfg.artfctdef.false_spindle.artifact  = spindlefalse_vec;
cfg.artfctdef.aselect_spindle.artifact  = events_spindle_select; 
spindlex = ft_databrowser(cfg,filt_data);
keyboard;

spindlemtimewin = [3850 6180] *1000; % provide analyzed data epoch here


%% spindle evaluate match

% select automatic spindles from manually scored epoch only
events_spindle_select = events_spindle(events_spindle(:,1)>spindlemtimewin(1) & events_spindle(:,2)<spindlemtimewin(2),:);

% remove automatic spindles from artifact and bad stage epochs
spindle_badcomb_vec = [];
for i = 1:length(spindle_badcomb)
    spindle_badcomb_vec = [spindle_badcomb_vec spindle_badcomb(i,1):spindle_badcomb(i,2)];
end
remove_idx = [];
for i = 1:size(events_spindle_select,1)
    if any(ismember([events_spindle_select(i,1):events_spindle_select(i,2)], spindle_badcomb_vec))
        remove_idx = [remove_idx i];
    end
end
events_spindle_select(remove_idx,:) = []; 

% construct continous event vectors
mspindle = [];
for i = 1:length(spindlex.artfctdef.spindle.artifact)
    mspindle= [mspindle spindlex.artfctdef.spindle.artifact(i,1):spindlex.artfctdef.spindle.artifact(i,2)];
end
aspindle = [];
for i = 1:length(events_spindle_select)
    aspindle = [aspindle events_spindle_select(i,1):events_spindle_select(i,2)];
end

% check matching
mspindletab = []; spindlehit = 0; spindlefalse = 0; spindlemiss = 0; spindlehit_vec = []; spindlemiss_vec = [];
for i = 1:length(spindlex.artfctdef.spindle.artifact)   
    if any(ismember([spindlex.artfctdef.spindle.artifact(i,1):spindlex.artfctdef.spindle.artifact(i,2)], aspindle))
        mspindletab(i) = 1;
        spindlehit = spindlehit +1;        
        % find automatic event overlapping with the hit manual event
        hitmevent = events_spindle_select(spindlex.artfctdef.spindle.artifact(i,1) <= events_spindle_select(:,2) & spindlex.artfctdef.spindle.artifact(i,2) >= events_spindle_select(:,1),:);
        [mindist mindistloc] = min(mean(hitmevent,2) - mean(spindlex.artfctdef.spindle.artifact(i,:))); % the center of which of multiple overlapping events is closest
        spindlehit_vec(spindlehit,:) = hitmevent(mindistloc,:); % use longest event in case of multiple hit events   
%         [maxdur, maxdurloc] = max(hitmevent(:,2)-hitmevent(:,1)); % which event is of longest  duration for multiple hits
%         spindlehit_vec(spindlehit,:) = hitmevent(maxdurloc,:); % use longest event in case of multiple hit events
%         spindlehit_vec(spindlehit,:) = hitmevent(1,:); % use first event in case of multiple hit events
    else
        mspindletab(i) = 0;
        spindlemiss = spindlemiss +1;
        spindlemiss_vec(spindlemiss,:) = spindlex.artfctdef.spindle.artifact(i,:);        
    end
end
aspindletab = [];  spindlefalse_vec = []; 
for i = 1:length(events_spindle_select)  
    if any(ismember([events_spindle_select(i,1):events_spindle_select(i,2)], mspindle))
        aspindletab(i) = 1;
    else
        aspindletab(i) = 0;
        spindlefalse = spindlefalse +1;
        spindlefalse_vec(spindlefalse,:) = events_spindle_select(i,:);        
    end
end


%% spindle display results 
aspindletotal= length(events_spindle_select);
mspindletotal = length(spindlex.artfctdef.spindle.artifact);
display(' ');
display(['SPINDLE VALIDATION for ' gcfg.subjectname ' in ' gcfg.targetchan]);
display(['spindle epoch: [' num2str(spindlemtimewin(1)) ' ' num2str(spindlemtimewin(2)) ']']);
display(['spindle number automatic: ' num2str(aspindletotal) '/' num2str(mspindletotal)]);
display(['spindle hits: ' num2str(spindlehit) ' (' num2str(round(spindlehit/mspindletotal*100)) '%)']);
display(['spindle misses: ' num2str(spindlemiss) ' (' num2str(round(spindlemiss/mspindletotal*100)) '%)']);
display(['spindle false alarms: ' num2str(spindlefalse)]);
save(fullfile(gcfg.paths.results,gcfg.subjectname,[gcfg.subjectname '_spindle_validation_' gcfg.targetchannel{gcfg.subjectno} '.mat']), 'mspindletotal', 'aspindletotal', 'spindlehit','spindlemiss','spindlefalse','spindlefalse_vec','spindlehit_vec','spindlemiss_vec','spindlex','spindlemtimewin');

keyboard;

end % of main function

%%