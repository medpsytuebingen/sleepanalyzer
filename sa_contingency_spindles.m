%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spindle contingencies 
% by Til Ole Bergmann 2016
% last modified 2016/11/22 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% performs contingency analyses across channels for spindles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spindle_contingencies] = sa_contingency_spindles(gcfg, All_spindle_info)


if  gcfg.loadPreviousContingencyResults ~= 1 % calculate anew or load existing contingencies?
    
    h1 = figure;

    for i = 1:numel(gcfg.subjects) % loop over subjects to generate logidx
        
        contMat{i} = zeros(length(All_spindle_info{i}.eventInfo),length(All_spindle_info{i}.eventInfo)); % setup contingency matrix
        logidx{i} = logical(zeros(length(All_spindle_info{i}.eventInfo),length(All_spindle_info{i}.scoring)));
        
        for ch = 1:length(All_spindle_info{i}.eventInfo) % first loop over channels
            for e = 1:length(All_spindle_info{i}.eventInfo(ch).startTime) % loop over events
                logidx{i}(ch, All_spindle_info{i}.eventInfo(ch).startTime(e):All_spindle_info{i}.eventInfo(ch).endTime(e)) = 1;
            end
        end % of first loop over channels
        
        
        for ch = 1:length(All_spindle_info{i}.eventInfo) % second loop over channels to test for condingencies (i.e. overlapping events)
            tempContMat{i}{ch} = zeros(length(All_spindle_info{i}.eventInfo(ch).startTime),length(All_spindle_info{i}.eventInfo));
            for e = 1:length(All_spindle_info{i}.eventInfo(ch).startTime) % loop over events
                tempContMat{i}{ch}(e,:) = sum(logidx{i}(:,All_spindle_info{i}.eventInfo(ch).startTime(e):All_spindle_info{i}.eventInfo(ch).endTime(e)),2) > 0';
            end % of loop over events
            contMat{i}(ch,:) = mean(tempContMat{i}{ch},1);
            
        end % of second loop over channels
        
        
        % plot single-subject figures
        h2 = figure;
        imagesc(contMat{i});
        colormap hot;
        set(gca,'clim',[0 1]);
        colorbar;
        set(gca,'XTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'XTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',10);
        set(gca,'YTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'YTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',10);
        set(gca,'yDir','normal');
        title({[gcfg.subjectNames{i} ': Spindle contingency matrix']});
        ylabel('Given a spindle in channel');
        xlabel('there is a spindle in channel');
        
        set(h2, 'Name',[gcfg.subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_contingencies_' gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
        filename = fullfile(gcfg.paths.results, gcfg.subjectNames{i}, [gcfg.subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_contingencies_' gcfg.timestamp]);
        saveas(h2,[filename '.fig']);
        set(h2, 'PaperUnits', 'centimeters', 'PaperSize', [33.867 19.05], 'PaperPosition', [0 -1 33.867 19.05]); % [left bottom width height]
        print(h2, '-dpdf', '-r300', [filename '.pdf']);
        
        % plot multi-subject figure
        figure(h1);
        subplot(3,4,i);
        imagesc(contMat{i});
        colormap hot;
        set(gca,'clim',[0 1]);
        colorbar;
        set(gca,'XTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'XTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',4);
        set(gca,'YTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'YTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',4);
        set(gca,'yDir','normal');
        title({[gcfg.subjectNames{i} ': Spindle contingency matrix']});
        ylabel('Given a spindle in channel');
        xlabel('there is a spindle in channel');
        
    end % of loop over subjects
    
    figure(h1);
    set(gcf, 'Name',['All_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_contingencies_' gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
    filename = fullfile(gcfg.paths.results, 'group', ['All_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_contingencies_' gcfg.timestamp]);
    saveas(gcf,[filename '.fig']);
    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [33.867 19.05], 'PaperPosition', [0 -1 33.867 19.05]); % [left bottom width height]
    print(gcf, '-dpdf', '-r300', [filename '.pdf']);
    
    
    
    % test against surrogate spindle contingencies %%%%%%%%%%%%%%%%%%%%%%%%    
    
    for i = 1:numel(gcfg.subjects) % loop over subjects
    tic; % for subjectToc
        % first set of loops to create shuffle cell array
        for ch = 1:length(All_spindle_info{i}.eventInfo) % loop over channels
            
            % define good epochs (sleepstages) only
            goodEpoch.idx = All_spindle_info{i}.scoring(4,:);
            goodEpoch.start = find(diff(goodEpoch.idx)>0);
            goodEpoch.end = find(diff(goodEpoch.idx)<0);
            st = All_spindle_info{i}.eventInfo(ch).startTime'; % start times
            et = All_spindle_info{i}.eventInfo(ch).endTime'; % end times
            
%                   % plot actual logidx
%                     dat = [logidx{i}; goodEpoch.idx]; % column 4 is 1 for 'searchstages'
%                     figure;
%                     imagesc(dat);
%                     set(gca,'YTick',1:numel(All_spindle_info{i}.extractchannel));
%                     set(gca,'YTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
            
            for ge = 1:length(goodEpoch.start) % loop over good epochs (GE, ge)
                
                % get distribution of event start and end times for good epoch (sleepstages) only
                GE{ge}.st = st(st >= goodEpoch.start(ge) & st < goodEpoch.end(ge));
                GE{ge}.et = et(et > goodEpoch.start(ge) & et <= goodEpoch.end(ge));
                eventDist{ge} =  [GE{ge}.et - GE{ge}.st];
                noEventDist{ge} = [GE{ge}.st; goodEpoch.end(ge)] - [goodEpoch.start(ge); GE{ge}.et];
                
                % build shuffle cell (SC) array to shuffle noEvents (0) and events(1) separately
                SC0{ch}{ge} = cell(1,length(noEventDist{ge}));
                SC1{ch}{ge} = cell(1,length(eventDist{ge}));
                for k = 1:length(SC0{ch}{ge}) % loop over noEvents
                    SC0{ch}{ge}{k} = logical(zeros(1,noEventDist{ge}(k)));
                end
                for k = 1:length(SC1{ch}{ge}) % loop over vents
                    SC1{ch}{ge}{k} = logical(ones(1,eventDist{ge}(k)));
                end
            end % of loop over good epochs
        end % of loop over channels
        

        % second set of loops for actual permutation test
        numIter = 10000; % number of iterations for permutation test
        randLogidx{i} = logical(zeros(size(logidx{i}))); % reset randLogidx
        randContMat{i} = zeros(numIter,length(All_spindle_info{i}.eventInfo),length(All_spindle_info{i}.eventInfo)); % setup surrogate contingency matrix
        for j = 1:numIter % surrogate permutation loop
            for ch = 1:length(All_spindle_info{i}.eventInfo) % loop over channels
                display([gcfg.subjectNames{i} ', channel ' All_spindle_info{i}.extractchannel{ch} ', surrogate iteration' num2str(j)]);
                for ge = 1:length(goodEpoch.start) % loop over good epochs (GE, ge)
                    % shuffle indices for SC0 and SC1
                    randIdx_SC0{ch}{ge} = randperm(length(SC0{ch}{ge}));
                    randIdx_SC1{ch}{ge} = randperm(length(SC1{ch}{ge}));
                    
                    % rebuild surrogate logindex
                    SCN{ch}{ge} = cell(1,length(SC0{ch}{ge})+length(SC1{ch}{ge}));
                    SCN{ch}{ge}(1:2:end) = SC0{ch}{ge}(randIdx_SC0{ch}{ge});
                    SCN{ch}{ge}(2:2:end) = SC1{ch}{ge}(randIdx_SC1{ch}{ge});
                    
                    % make logical array and add missing sample point...(?)
                    zerofill = zeros(1,length(randLogidx{i}(ch,goodEpoch.start(ge):goodEpoch.end(ge)))-length([SCN{ch}{ge}{:}]));
                    randLogidx{i}(ch,goodEpoch.start(ge):goodEpoch.end(ge)) = [[SCN{ch}{ge}{:}] zerofill];
                end % of loop over good epochs
            end % of loop over channels
            
            % build surrogate contingency table for
            %         randContMat{i}{j} = zeros(length(All_spindle_info{i}.eventInfo),length(All_spindle_info{i}.eventInfo)); % setup contingency matrix
            
            for ch = 1:length(All_spindle_info{i}.eventInfo) % second loop over channels to test for condingencies (i.e. overlapping events)
                tempContMat{i}{ch} = zeros(length(All_spindle_info{i}.eventInfo(ch).startTime),length(All_spindle_info{i}.eventInfo));
                for e = 1:length(All_spindle_info{i}.eventInfo(ch).startTime) % loop over events
                    tempContMat{i}{ch}(e,:) = sum(randLogidx{i}(:,All_spindle_info{i}.eventInfo(ch).startTime(e):All_spindle_info{i}.eventInfo(ch).endTime(e)),2) > 0';
                end % of loop over events
                randContMat{i}(j,ch,:) = mean(tempContMat{i}{ch},1);
            end % of second loop over channels
        end % of surrogate permutation loop
        
        
        % permutation test
        sortRandContMat{i} = sort(randContMat{i},1);
        for m = 1:size(sortRandContMat{i},2) % loop over rows
            for n = 1:size(sortRandContMat{i},3) % loop over columns
                pContMat{i}(n,m) = mean(squeeze(sortRandContMat{i}(:,m,n)) > contMat{i}(n,m)); % p-value
            end
        end
        
        % save results
        spindle_contingencies.contMat = contMat{i};
        spindle_contingencies.pContMat = pContMat{i};
        spindle_contingencies.randContMat = randContMat{i};
        save('-v7.3', fullfile(gcfg.resultsPath,gcfg.subjectNames{i},[gcfg.subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_targetchannel{gcfg.subjects(i)} '_contingencies.mat']), 'spindle_contingencies');
        
        
        % plot relevant contingency parameters
        h3 = figure;
        colormap hot;
        
        subplot(2,2,1);
        imagesc(squeeze(sortRandContMat{i}(1,:,:)));
        title({[gcfg.subjectNames{i} ': MIN SURROGATE Spindle contingency matrix']});
        set(gca,'clim',[0 1]);
        colorbar;
        set(gca,'XTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'XTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'YTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'YTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'yDir','normal');
        ylabel('Given a spindle in channel');
        xlabel('there is a spindle in channel');
        
        subplot(2,2,2);
        imagesc(squeeze(sortRandContMat{i}(end,:,:)));
        title({[gcfg.subjectNames{i} ': MAX SURROGATE Spindle contingency matrix']});
        set(gca,'clim',[0 1]);
        colorbar;
        set(gca,'XTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'XTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'YTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'YTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'yDir','normal');
        ylabel('Given a spindle in channel');
        xlabel('there is a spindle in channel');
        
        subplot(2,2,3);
        imagesc(pContMat{i});
        title({[gcfg.subjectNames{i} ': P-values for Spindle contingency matrix']});
        set(gca,'clim',[0 1]);
        colorbar;
        set(gca,'XTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'XTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'YTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'YTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'yDir','normal');
        ylabel('Given a spindle in channel');
        xlabel('there is a spindle in channel');
        
        subplot(2,2,4);
        imagesc(contMat{i} .* (pContMat{i}<0.0001));
        title({[gcfg.subjectNames{i} ': Significant (p < 0.0001) values in Spindle contingency matrix']});
        set(gca,'clim',[0 1]);
        colorbar;
        set(gca,'XTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'XTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'YTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'YTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'yDir','normal');
        ylabel('Given a spindle in channel');
        xlabel('there is a spindle in channel');
        
        suptitle([gcfg.subjectNames{i} ': Contingency statistics with' num2str(numIter) ' iterations']);
        
        set(h3, 'Name',[gcfg.subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_contingencySTATS_' gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
        filename = fullfile(gcfg.paths.results, gcfg.subjectNames{i}, [gcfg.subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_contingencySTATS_' gcfg.timestamp]);
        saveas(h3,[filename '.fig']);
        set(h3, 'PaperUnits', 'centimeters', 'PaperSize', [33.867 19.05], 'PaperPosition', [0 -1 33.867 19.05]); % [left bottom width height]
        print(h3, '-dpdf', '-r300', [filename '.pdf']);
        
        
        %     % plot single-subject sorted surrogate figures
        %     figure;
        %     for sp = 1:min(12,numIter)
        %         subplot(3,4,sp)
        %         imagesc(squeeze(sortRandContMat{i}(sp,:,:)));
        %         colormap hot;
        %         set(gca,'clim',[0 1]);
        %         colorbar;
        %         set(gca,'XTick',1:numel(All_spindle_info{i}.extractchannel));
        %         set(gca,'XTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',4);
        %         set(gca,'YTick',1:numel(All_spindle_info{i}.extractchannel));
        %         set(gca,'YTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',4);
        %         set(gca,'yDir','normal');
        %         title({[gcfg.subjectNames{i} ': Spindle contingency matrix']});
        %         ylabel('Given a spindle in channel');
        %         xlabel('there is a spindle in channel');
        %     end
        
        subjectToc{i} = toc;
        display(['Contingency permutation test for subject ' gcfg.subjectNames{i} ' took ' num2str(subjectToc{i}) ' seconds.']);
    end % of loop over subjects

elseif gcfg.loadPreviousContingencyResults == 1
    for i = 1:numel(gcfg.subjects)
        load(fullfile(gcfg.resultsPath,gcfg.subjectNames{i},[gcfg.subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_targetchannel{gcfg.subjects(i)} '_contingencies.mat']), 'spindle_contingencies');
        contMat{i} = spindle_contingencies.contMat;        
        pContMat{i} = spindle_contingencies.pContMat;
        randContMat{i} = spindle_contingencies.randContMat;
        clear spindle_contingencies;        
    end
    keyboard;
    
    for i = 1:numel(gcfg.subjects)
        numHCTargetChan(i) = find(strcmp(All_spindle_info{i}.timelockchannel, All_spindle_info{i}.extractchannel));
        numCzTargetChan(i) = find(strcmp('Cz', All_spindle_info{i}.extractchannel));
        
        for ch = 1:length(All_spindle_info{i}.eventInfo)        
            Cz2HCP(i) = contMat{i}(numCzTargetChan(i),numHCTargetChan(i));
            Hc2CzP(i) = contMat{i}(numHCTargetChan(i),numCzTargetChan(i));
            
        end
    end
        
    % paired t-test
    [h,p,ci,stats] = ttest(Cz2HCP,Hc2CzP)
    
    
    %% TEMP re-plotting for conference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % plot relevant contingency parameters
        h3r = figure;
        colormap jet;
        
        subplot(2,2,1);
        imagesc(squeeze(sortRandContMat{i}(1,:,:)));
        title({[gcfg.subjectNames{i} ': MIN SURROGATE Spindle contingency matrix']});
        set(gca,'clim',[0 1]);
        colorbar;
        set(gca,'XTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'XTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'YTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'YTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'yDir','normal');
        ylabel('Given a spindle in channel');
        xlabel('there is a spindle in channel');
        
        subplot(2,2,2);
        imagesc(squeeze(sortRandContMat{i}(end,:,:)));
        title({[gcfg.subjectNames{i} ': MAX SURROGATE Spindle contingency matrix']});
        set(gca,'clim',[0 1]);
        colorbar;
        set(gca,'XTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'XTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'YTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'YTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'yDir','normal');
        ylabel('Given a spindle in channel');
        xlabel('there is a spindle in channel');
        
        subplot(2,2,3);
        imagesc(pContMat{i});
        title({[gcfg.subjectNames{i} ': P-values for Spindle contingency matrix']});
        set(gca,'clim',[0 1]);
        colorbar;
        set(gca,'XTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'XTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'YTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'YTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'yDir','normal');
        ylabel('Given a spindle in channel');
        xlabel('there is a spindle in channel');
        
        subplot(2,2,4);
        imagesc(contMat{i} .* (pContMat{i}<0.0001));
        title({[gcfg.subjectNames{i} ': Significant (p < 0.0001) values in Spindle contingency matrix']});
        set(gca,'clim',[0 1]);
        colorbar;
        set(gca,'XTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'XTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'YTick',1:numel(All_spindle_info{i}.extractchannel));
        set(gca,'YTickLabel',All_spindle_info{i}.extractchannel, 'FontSize',8);
        set(gca,'yDir','normal');
        ylabel('Given a spindle in channel');
        xlabel('there is a spindle in channel');
        
        suptitle([gcfg.subjectNames{i} ': Contingency statistics with' num2str(numIter) ' iterations']);
        
        set(h3r, 'Name',[gcfg.subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_contingencySTATSjet_' gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
        filename = fullfile(gcfg.paths.results, gcfg.subjectNames{i}, [gcfg.subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_contingencySTATSjet_' gcfg.timestamp]);
        saveas(h3r,[filename '.fig']);
        set(h3r, 'PaperUnits', 'centimeters', 'PaperSize', [33.867 19.05], 'PaperPosition', [0 -1 33.867 19.05]); % [left bottom width height]
        print(h3r, '-dpdf', '-r300', [filename '.pdf']);
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    

end % of if statement deciding whether or not to load vs. calcualte contingencies

end % of function