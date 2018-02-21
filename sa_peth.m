%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peth (peri-event time histograms)
% by Til Ole Bergmann 2013
% last modified 2016/11/22 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calculates peri-event time histograms
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sa_peth(gcfg)
display('Generating peri-event time histograms...');
SO_peth_already_plotted = 0; % reset
spindle_peth_already_plotted = 0; % reset


for tc = 1:numel(gcfg.tc) % loop over trial classes
    
    %% define events for current trial class
    switch gcfg.tc{tc}
        % 'stage' 'startTime' 'midTime' 'endTime' 'duration' 'maxTime' 'minTime' 'minAmp' 'maxAmp' 'p2pAmp' 'p2pTime'
               
        
        case {'SO_all','SO_w_spindle','SO_wo_spindle'}
            
            %%% spindle center relative to SW down-state
            % event 1 (reference event)
            gcfg.event1.targetchannel = gcfg.event1channelname; % 'HC' or 'Cz'
            gcfg.event1.chan = gcfg.event1channel;
            gcfg.event1.spec = ['SO_' gcfg.SOspec '_' gcfg.SO_timelockevent];
            gcfg.event1.eventdata = 'SO_eventdata';
            gcfg.event1.tic = 7; % minTime (trialinfo column containing reference time point for event 1)
            
            % event 2 (target event)
            gcfg.event2.targetchannel = gcfg.event2channelname; % 'HC'or 'Cz'
            gcfg.event2.chan = gcfg.event2channel;
            gcfg.event2.spec = ['spindle_' gcfg.spindle_timelockevent];
            gcfg.event2.eventdata = 'spindle_eventdata';
            gcfg.event2.tic = 3; % midTime (trialinfo column containing reference time point for event 2)
            gcfg.eow = [0.2 1.0]; % event occurence window in seconds
            
            
        case {'spindle_all', 'spindle_w_SO', 'spindle_wo_SO'}
            
            %%% SW down-state relative to spindle center
            % event 1 (reference event)
            gcfg.event1.targetchannel = gcfg.event1channelname; % 'HC' or 'Cz'
            gcfg.event1.chan = gcfg.event1channel;
            gcfg.event1.spec = ['spindle_' gcfg.spindle_timelockevent];
            gcfg.event1.eventdata = 'spindle_eventdata';
            gcfg.event1.tic = 3; % midTime (trialinfo column containing reference time point for event 1)
            
            % event 2 (counted event)
            gcfg.event2.targetchannel = gcfg.event2channelname; % 'HC' or 'Cz'
            gcfg.event2.chan = gcfg.event2channel;
            gcfg.event2.spec = ['SO_' gcfg.SOspec '_' gcfg.SO_timelockevent];
            gcfg.event2.eventdata = 'SO_eventdata';
            gcfg.event2.tic = 7; % minTime (trialinfo column containing reference time point for event 2)
            gcfg.eow = [-1 -0.2]; % event occurence window in seconds
            
            
        case {'SO_w_SO','SO_wo_SO'}            
            %%% SW down-state relative to SW down-state
            % event 1 (reference event)
            gcfg.event1.targetchannel = gcfg.event1channelname; % 'HC' or 'Cz'
            gcfg.event1.chan = gcfg.event1channel;
            gcfg.event1.spec = ['SO_' gcfg.SOspec '_' gcfg.SO_timelockevent];
            gcfg.event1.eventdata = 'SO_eventdata';
            gcfg.event1.tic = 7; % minTime (trialinfo column containing reference time point for event 1)
            
            % event 2 (counted event)
            gcfg.event2.targetchannel = gcfg.event2channelname; % 'HC' or 'Cz'
            gcfg.event2.chan = gcfg.event2channel;
            gcfg.event2.spec = ['SO_' gcfg.SOspec '_' gcfg.SO_timelockevent];
            gcfg.event2.eventdata = 'SO_eventdata';
            gcfg.event2.tic = 7; % minTime (trialinfo column containing reference time point for event 2)
            gcfg.eow = [-2.5 2.5]; % event occurence window in seconds
            
    end % of switch

    
    %% load individual data
    for i = 1:numel(gcfg.subjects) % loop over subjects
        temp1 = load(fullfile(gcfg.resultsPath,gcfg.subjectNames{i},[gcfg.subjectNames{i} '_' gcfg.event1.spec '_' gcfg.event1.chan{gcfg.subjects(i)} '_eventdata.mat']), gcfg.event1.eventdata);
        temp2 = load(fullfile(gcfg.resultsPath,gcfg.subjectNames{i},[gcfg.subjectNames{i} '_' gcfg.event2.spec '_' gcfg.event2.chan{gcfg.subjects(i)} '_eventdata.mat']), gcfg.event2.eventdata);
        
        % 'stage' 'startTime' 'midTime' 'endTime' 'duration' 'maxTime' 'minTime' 'minAmp' 'maxAmp' 'p2pAmp' 'p2pTime'
        fsample = temp1.(gcfg.event1.eventdata).fsample;
        event1{i} = temp1.(gcfg.event1.eventdata).trialinfo;
        event2{i} = temp2.(gcfg.event2.eventdata).trialinfo;
        
        
        % apply esv to event1
        if strcmp(gcfg.event1.eventdata, 'SO_eventdata') && isfield(gcfg,'SO_esv')
            event1{i} = event1{i}(logical(gcfg.SO_esv{i}),:);
        elseif strcmp(gcfg.event1.eventdata, 'spindle_eventdata') && isfield(gcfg,'spindle_esv')
            event1{i} = event1{i}(logical(gcfg.spindle_esv{i}),:);            
        end
        
        % apply esv to event2
        if strcmp(gcfg.event2.eventdata, 'SO_eventdata') && isfield(gcfg,'SO_esv')
            event2{i} = event2{i}(logical(gcfg.SO_esv{i}),:);            
        elseif strcmp(gcfg.event2.eventdata, 'spindle_eventdata') && isfield(gcfg,'spindle_esv')
            event2{i} = event2{i}(logical(gcfg.spindle_esv{i}),:);            
        end
        
    end
    
    
    %% setup search window
    edges = gcfg.edges * fsample; % gcfg.edges in datapoints
    eow = gcfg.eow * fsample; % gcfg.eow in datapoints
    eowidx = find(ismember(edges, eow));    
    eowidx(2) = eowidx(2)-1; % to ensure the right limit of eow does not cause inclusion of another bin starting at eow(2)!!!
    zerobin = find(ismember(edges, 0));    

    %% calculate PETH, tsv and PEO
    for i = 1:numel(event1) % loop over subjects
        
        % count target events in peri-event neighbourhood around reference event
        tsv{i} = ones(1,size(event1{i},1)); % preset tsv
        AllTrial_peth{i} = zeros(size(edges,2),size(event1{i},1)); % preset PETH

        for j = 1:size(event1{i},1) % loop over events of event type 1
            
            [AllTrial_peth{i}(:,j),bin] = histc(event2{i}(:,round(gcfg.event2.tic)),[event1{i}(j,round(gcfg.event1.tic)) + edges]);
            
            switch gcfg.tc{tc} % set tsv value for current trial
                case 'SO_all'
                    tsv{i}(j) = tsv{i}(j);
                case 'SO_w_spindle'
                    tsv{i}(j) = any(AllTrial_peth{i}(eowidx(1):eowidx(2),j)); % set tsv value for current trial
                case 'SO_wo_spindle'
                    tsv{i}(j) = ~any(AllTrial_peth{i}(eowidx(1):eowidx(2),j)); % set tsv value for current trial
                case 'SO_w_SO'
                    tsv{i}(j) = any(AllTrial_peth{i}([eowidx(1):zerobin-1, zerobin+1:eowidx(2)],j)); % set tsv value for current trial (auto-correlation removed)
                case 'SO_wo_SO'
                    tsv{i}(j) = ~any(AllTrial_peth{i}([eowidx(1):zerobin-1, zerobin+1:eowidx(2)],j)); % set tsv value for current trial (auto-correlation removed)                    
                case 'spindle_all'
                    tsv{i}(j) = tsv{i}(j);
                case 'spindle_w_SO'
                    tsv{i}(j) = any(AllTrial_peth{i}(eowidx(1):eowidx(2),j)); % set tsv value for current trial
                case 'spindle_wo_SO'
                    tsv{i}(j) = ~any(AllTrial_peth{i}(eowidx(1):eowidx(2),j)); % set tsv value for current trial
            end
            
        end % of loop over events of event type 1
        
        
        switch gcfg.tc{tc} % calculate PEO 
            case 'SO_all'
               All_PEO{i} = 1;
            case 'SO_w_spindle'
               All_PEO{i} = mean(any(AllTrial_peth{i}(eowidx(1):eowidx(2),:),1));
            case 'SO_wo_spindle'
               All_PEO{i} = 1-mean(any(AllTrial_peth{i}(eowidx(1):eowidx(2),:),1));
            case 'SO_w_SO'
               All_PEO{i} = mean(any(AllTrial_peth{i}([eowidx(1):zerobin-1, zerobin+1:eowidx(2)],:),1)); % probability of target event to occur in specified event occurence window (eow) around reference event
            case 'SO_wo_SO'
               All_PEO{i} = 1-mean(any(AllTrial_peth{i}([eowidx(1):zerobin-1, zerobin+1:eowidx(2)],:),1)); % probability of target event to occur in specified event occurence window (eow) around reference event
            case 'spindle_all'
               All_PEO{i} = 1;
            case 'spindle_w_SO'
               All_PEO{i} = mean(any(AllTrial_peth{i}(eowidx(1):eowidx(2),:),1));
            case 'spindle_wo_SO'
               All_PEO{i} = 1 - mean(any(AllTrial_peth{i}(eowidx(1):eowidx(2),:),1));
        end   
        All_peth{i} = mean(AllTrial_peth{i},2);
        

        %% save tsv (trial selection vector)
        peo = All_PEO{i};
        peth = All_peth{i};
        switch gcfg.tc{tc}
            case 'SO_all'
                tsv_SO_all = tsv{i};
                save('-v7.3', fullfile(gcfg.resultsPath,gcfg.subjectNames{i},[gcfg.subjectNames{i} '_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_targetchannel{gcfg.subjects(i)} '_tsv_SO_all.mat']), 'tsv_SO_all', 'peo', 'peth');
            case 'SO_w_spindle'
                tsv_SO_w_spindle = tsv{i};
                save('-v7.3', fullfile(gcfg.resultsPath,gcfg.subjectNames{i},[gcfg.subjectNames{i} '_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_targetchannel{gcfg.subjects(i)} '_tsv_SO_w_spindle.mat']), 'tsv_SO_w_spindle', 'peo', 'peth');
            case 'SO_wo_spindle'
                tsv_SO_wo_spindle = tsv{i};
                save('-v7.3', fullfile(gcfg.resultsPath,gcfg.subjectNames{i},[gcfg.subjectNames{i} '_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_targetchannel{gcfg.subjects(i)} '_tsv_SO_wo_spindle.mat']), 'tsv_SO_wo_spindle', 'peo', 'peth');
            case 'SO_w_SO'
                tsv_SO_w_SO = tsv{i};
                save('-v7.3', fullfile(gcfg.resultsPath,gcfg.subjectNames{i},[gcfg.subjectNames{i} '_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_targetchannel{gcfg.subjects(i)} '_tsv_SO_w_SO.mat']), 'tsv_SO_w_SO', 'peo', 'peth');                
            case 'SO_wo_SO'
                tsv_SO_wo_SO = tsv{i};
                save('-v7.3', fullfile(gcfg.resultsPath,gcfg.subjectNames{i},[gcfg.subjectNames{i} '_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_targetchannel{gcfg.subjects(i)} '_tsv_SO_wo_SO.mat']), 'tsv_SO_wo_SO', 'peo', 'peth');                                
            case 'spindle_all'
                tsv_spindle_all = tsv{i};
                save('-v7.3', fullfile(gcfg.resultsPath,gcfg.subjectNames{i},[gcfg.subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_targetchannel{gcfg.subjects(i)} '_tsv_spindle_all.mat']), 'tsv_spindle_all', 'peo', 'peth');
            case 'spindle_w_SO'
                tsv_spindle_w_SO = tsv{i};
                save('-v7.3', fullfile(gcfg.resultsPath,gcfg.subjectNames{i},[gcfg.subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_targetchannel{gcfg.subjects(i)} '_tsv_spindle_w_SO.mat']), 'tsv_spindle_w_SO', 'peo', 'peth');
            case 'spindle_wo_SO'
                tsv_spindle_wo_SO = tsv{i};
                save('-v7.3', fullfile(gcfg.resultsPath,gcfg.subjectNames{i},[gcfg.subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_targetchannel{gcfg.subjects(i)} '_tsv_spindle_wo_SO.mat']), 'tsv_spindle_wo_SO', 'peo', 'peth');
        end
        
    end % of loop over subjects

    
    %% calculate group results
    % PETH (peri-event time histograms)
    GA_peth = mean([All_peth{:}],2); % average over trials AND subjects...
    
    % PEO (probability of event occurence
    GA_PEO = mean([All_PEO{:}]); % average over trials AND subjects...
    GSD_PEO = std([All_PEO{:}]); % sd over trials AND subjects...
    
    
    %% save group results
    switch gcfg.tc{tc}
        case 'SO_all'
            save('-v7.3', fullfile(gcfg.resultsPath,'group',['GA_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_targetchannel{gcfg.subjects(i)} '_tsv_SO_all.mat']), 'tsv_SO_all', 'GA_PEO', 'GSD_PEO', 'GA_peth');
        case 'SO_w_spindle'
            save('-v7.3', fullfile(gcfg.resultsPath,'group',['GA_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_targetchannel{gcfg.subjects(i)} '_tsv_SO_w_spindle.mat']), 'tsv_SO_w_spindle', 'GA_PEO', 'GSD_PEO', 'GA_peth');
        case 'SO_wo_spindle'
            save('-v7.3', fullfile(gcfg.resultsPath,'group',['GA_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_targetchannel{gcfg.subjects(i)} '_tsv_SO_wo_spindle.mat']), 'tsv_SO_wo_spindle', 'GA_PEO', 'GSD_PEO', 'GA_peth');
        case 'SO_w_SO'
            save('-v7.3', fullfile(gcfg.resultsPath,'group',['GA_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_targetchannel{gcfg.subjects(i)} '_tsv_SO_w_SO.mat']), 'tsv_SO_w_SO', 'GA_PEO', 'GSD_PEO', 'GA_peth');
        case 'SO_wo_SO'
            save('-v7.3', fullfile(gcfg.resultsPath,'group',['GA_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_targetchannel{gcfg.subjects(i)} '_tsv_SO_wo_SO.mat']), 'tsv_SO_wo_SO', 'GA_PEO', 'GSD_PEO', 'GA_peth');            
        case 'spindle_all'
            save('-v7.3', fullfile(gcfg.resultsPath,'group',['GA_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_targetchannel{gcfg.subjects(i)} '_tsv_spindle_all.mat']), 'tsv_spindle_all', 'GA_PEO', 'GSD_PEO', 'GA_peth');
        case 'spindle_w_SO'
            save('-v7.3', fullfile(gcfg.resultsPath,'group',['GA_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_targetchannel{gcfg.subjects(i)} '_tsv_spindle_w_SO.mat']), 'tsv_spindle_w_SO', 'GA_PEO', 'GSD_PEO', 'GA_peth');
        case 'spindle_wo_SO'
            save('-v7.3', fullfile(gcfg.resultsPath,'group',['GA_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_targetchannel{gcfg.subjects(i)} '_tsv_spindle_wo_SO.mat']), 'tsv_spindle_wo_SO', 'GA_PEO', 'GSD_PEO', 'GA_peth');
    end
    
    
    %% display results
    
    filename = fullfile(gcfg.resultsPath, 'group', ['GA_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.event1.targetchannel '_PETH_' gcfg.timestamp '.txt']);
    diary(filename);
    display([gcfg.tc{tc} ':']);
    display(['Probability of counting event2 occuring (w) / not occuring (wo) within [' num2str(gcfg.eow(1)) ' ' num2str(gcfg.eow(2)) '] ms releative to reference event1:']);
    display('Subject-wise:');
    display(num2str([All_PEO{:}]));
    display(' ');
    display('Group:');
    display(['Mean: ' num2str(GA_PEO) ' and SD: ' num2str(GSD_PEO)]);
    display(' ');
    display(' ');
    diary off;       
    
    
    %% plot PETH
    if gcfg.plot
        
        % apply tsv to PETH
        for i = 1:numel(event1) % loop over subjects
            AllTrial_peth_plot{i} = AllTrial_peth{i}(:,logical(tsv{i})); % REcalculate Trial results for condition-wise plotting
            All_peth_plot{i} = mean(AllTrial_peth_plot{i},2); % REcalculate subject results for condition-wise plotting
        end
        GA_peth_plot = mean([All_peth_plot{:}],2); % REcalculate group results for condition-wise plotting
        
        
        if any(ismember(gcfg.tc{tc},{'SO_all','SO_w_spindle','SO_wo_spindle','SO_w_SO','SO_wo_SO'})) % && ~SO_peth_already_plotted
            
            % remove counts for event itself 
            if any(ismember(gcfg.tc{tc},{'SO_w_SO','SO_wo_SO'}))
                GA_peth_plot(ceil(size(GA_peth_plot,1)/2)) = 0;
            end
            
            % plot group peth
            x = edges/fsample;            
            figure;
            bar(x,GA_peth_plot); % relative frequency
            ylabel(['probability of ' gcfg.event2.spec ' in ' gcfg.event2.targetchannel], 'interpreter', 'none');
            xlabel(['time (s) relative to ' gcfg.event1.spec ' in ' gcfg.event1.targetchannel], 'interpreter', 'none');
            title(['Peri-Event Time Histogram (PETH) for group of ' gcfg.tc{tc} ], 'interpreter', 'none');
                        
            % save group figure
            set(gcf, 'Name',['GA_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.event1.targetchannel '_PETH_' gcfg.tc{tc} '_' gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
            filename = fullfile(gcfg.resultsPath, 'group', ['GA_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.event1.targetchannel '_PETH_' gcfg.tc{tc} '_' gcfg.timestamp]);
            saveas(gcf,[filename '.fig']);
            set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [25.4 19.05], 'PaperPosition', [-3 -1 30 19.5]);
            print(gcf, '-dpdf', '-r300', [filename '.pdf']);
            
            % plot individual PETH
            for i = 1:numel(All_peth) % loop over subjects
                
                % remove counts for event itself 
                if any(ismember(gcfg.tc{tc},{'SO_w_SO','SO_wo_SO'}))
                    All_peth_plot{i}(ceil(size(All_peth_plot{i},1)/2)) = 0;
                end
                
                x = edges/fsample;
                figure;
                bar(x,All_peth_plot{i}); % relative frequency
                ylabel(['probability of ' gcfg.event2.spec ' in ' gcfg.event2.targetchannel], 'interpreter', 'none');
                xlabel(['time (s) relative to ' gcfg.event1.spec ' in ' gcfg.event1.targetchannel], 'interpreter', 'none');
                title(['Peri-Event Time Histogram (PETH) for ' gcfg.subjectNames{i} ' of ' gcfg.tc{tc} ], 'interpreter', 'none');
                
                % save individual figure
                set(gcf, 'Name',[gcfg.subjectNames{i} '_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.event1.chan{gcfg.subjects(i)} '_PETH_' gcfg.tc{tc} '_' gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
                filename = fullfile(gcfg.resultsPath, gcfg.subjectNames{i}, [gcfg.subjectNames{i} '_SO_' gcfg.SOspec '_'  gcfg.spindle_timelockevent '_' gcfg.event1.chan{gcfg.subjects(i)} '_PETH_' gcfg.tc{tc} '_' gcfg.timestamp]);
                saveas(gcf,[filename '.fig']);
                set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [25.4 19.05], 'PaperPosition', [-3 -1 30 19.5]);
                print(gcf, '-dpdf', '-r300', [filename '.pdf']);
                close(gcf); % close current individual figures (would be way too many open figures!)
            end
            
%             SO_peth_already_plotted = 1;
        end
         
        
        if any(ismember(gcfg.tc{tc},{'spindle_all', 'spindle_w_SO', 'spindle_wo_SO'})) % && ~spindle_peth_already_plotted
            x = edges/fsample;
            figure;
            bar(x,GA_peth_plot); % relative frequency
            ylabel(['probability of ' gcfg.event2.spec ' in ' gcfg.event2.targetchannel], 'interpreter', 'none');
            xlabel(['time (s) relative to ' gcfg.event1.spec ' in ' gcfg.event1.targetchannel], 'interpreter', 'none');
            title(['Peri-Event Time Histogram (PETH) for group of ' gcfg.tc{tc} ], 'interpreter', 'none');
            
            % save group figure
            set(gcf, 'Name',['GA_spindle_' gcfg.spindle_timelockevent '_' gcfg.event1.targetchannel '_PETH_' gcfg.tc{tc} '_' gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
            filename = fullfile(gcfg.resultsPath, 'group', ['GA_spindle_' gcfg.spindle_timelockevent '_' gcfg.event1.targetchannel '_PETH_' gcfg.tc{tc} '_' gcfg.timestamp]);
            saveas(gcf,[filename '.fig']);
            set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [25.4 19.05], 'PaperPosition', [-3 -1 30 19.5]);
            print(gcf, '-dpdf', '-r300', [filename '.pdf']);            
            
            % plot individual PETH
            for i = 1:numel(All_peth_plot) % loop over subjects
                x = edges/fsample;
                figure;
                bar(x,All_peth_plot{i}); % relative frequency
                ylabel(['probability of ' gcfg.event2.spec ' in ' gcfg.event2.targetchannel], 'interpreter', 'none');
                xlabel(['time (s) relative to ' gcfg.event1.spec ' in ' gcfg.event1.targetchannel], 'interpreter', 'none');
                title(['Peri-Event Time Histogram (PETH) for ' gcfg.subjectNames{i} ' of ' gcfg.tc{tc} ], 'interpreter', 'none');
                
                % save individual figure
                set(gcf, 'Name',[gcfg.subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.event1.chan{gcfg.subjects(i)} '_PETH_' gcfg.tc{tc} '_' gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
                filename = fullfile(gcfg.resultsPath, gcfg.subjectNames{i}, [gcfg.subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.event1.chan{gcfg.subjects(i)} '_PETH_' gcfg.tc{tc} '_' gcfg.timestamp]);
                saveas(gcf,[filename '.fig']);
                set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [25.4 19.05], 'PaperPosition', [-3 -1 30 19.5]);
                print(gcf, '-dpdf', '-r300', [filename '.pdf']);
                close(gcf); % close current individual figures (would be way too many open figures!)
                
            end
%             spindle_peth_already_plotted = 1;
        end
    end
    

end % of loop over trial classes

end % of function




