%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SO CFC statistics
% by Til Ole Bergmann 2013
% last modified 2016/11/22 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% performs statistics for slow oscillation CFC 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results] = sa_CFC_stats_SOs(gcfg, GA_SO_freq, GA_SO_timelock, All_SO_freq, All_SO_timelock, SO_tsv)  

results = struct; % initialize results structure

%% check fields
if ~isfield(gcfg,'thisCondition')
    gcfg.thisCondition = 1;
end    
if ~isfield(gcfg,'conditionLabels')
    if ~isfiel(gcfg,'conditions')
        gcfg.conditionLabels = {''};
    elseif isfiel(gcfg,'conditions')
        gcfg.conditionLabels = repmat({''},size(gcfg.conditions));
    end
end    

%% options
opt.clusterstats = 0;
opt.plot = 0;
opt.saveplot = 0;

opt.tfrstats = 0;
opt.pacpt = 1; % phase amplitude coupling per trial

%% TFR non-parametric t-statistics against zero
if opt.tfrstats  

    All_zero_freq = All_SO_freq;
    for i = 1:length(All_SO_freq) % create fake zero data
        All_zero_freq{i}.powspctrm(~isnan(All_zero_freq{i}.powspctrm)) = 0;
    end    

    cfg = [];
    cfg.channel = 'all';
    cfg.latency = gcfg.plot_SO_xlim; % [begin end] in seconds; 'all'
    cfg.trials = 'all';
    cfg.frequency = gcfg.plot_SO_ylim; % [begin end]; 'all'
    cfg.parameter = 'powspctrm';   
    cfg.method = 'montecarlo';
    cfg.statistic = 'depsamplesT'; % dependent samples T-statistic	
    cfg.design(1,:) = [repmat(1:numel(All_SO_freq),1,2)]; % subject vector        
    cfg.design(2,:) = [ones(1,length(All_SO_freq)) ones(1,length(All_SO_freq))*2]; % condition vector        
    cfg.uvar = 1; % number or list with indices, unit variable(s)
    cfg.ivar = 2; % number or list with indices, independent variable(s)   
    cfg.numrandomization = 'all'; % 'all'
    cfg.correctm = 'no'; % 'no' , 'max', cluster', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')
    cfg.alpha = 0.05; 
    cfg.tail = 0; %, -1, 1 or 0 (default = 0)   
    cfg.correcttail = 'prob'; % correct p-values or alpha-values when doing a two-sided test, 'alpha','prob' or 'no' (default = 'no')
    cfg.feedback = 'text'; % 'gui', 'text', 'textbar' or 'no' (default = 'textbar')    
    cfg.randomseed = 'yes'; % string, 'yes', 'no' or a number (default = 'yes')       
    [stat] = ft_freqstatistics(cfg, All_SO_freq{:}, All_zero_freq{:}); 
    %     [stat] = ft_freqstatistics(cfg, All_SO_freq{:}, All_SO_surrogate_freq{:}); 
    %     [stat] = ft_freqstatistics(cfg, All_SO_freq{:}, All_SO_free_freq{:}); 

    if opt.clusterstats
        cfg.correctm = 'cluster'; % 'no' , 'max', cluster', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')
        cfg.clusterstatistic = 'maxsum'; % how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
        cfg.clusterthreshold = 'parametric'; % method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
        cfg.clusteralpha = 0.05; % for either parametric or nonparametric thresholding per tail (default = 0.05)
        cfg.clustertail = 0; % -1, 1 or 0 (default = 0)        
        [stat_cor] = ft_freqstatistics(cfg, All_SO_freq{:}, All_zero_freq{:});     
        %     [stat_cor] = ft_freqstatistics(cfg, All_SO_freq{:}, All_SO_surrogate_freq{:}); 
        %     [stat_cor] = ft_freqstatistics(cfg, All_SO_freq{:}, All_SO_free_freq{:}); 
    end

    % cluster extent threshold correction
    clustersig = 1;
    k = 200;
    if clustersig    
        x = squeeze(stat.mask);
        L = bwlabel(x,4);
        clust = unique(L);
        nclust = numel(unique(L));    
        M = zeros(size(L));
        for iclust = 2:nclust
            if numel(find(L==clust(iclust))) >= k
                M(L==clust(iclust)) = 1;
            end
        end
        stat.mask_extthresh = permute(M, [3,1,2]);
    end


    %% plot stats
    if opt.plot
        figure;
        subplot(4,1,1);
        plotdata = squeeze(GA_SO_freq.powspctrm(1,GA_SO_freq.freq >= cfg.frequency(1) & GA_SO_freq.freq <= cfg.frequency(2), GA_SO_freq.time >= cfg.latency(1) & GA_SO_freq.time <= cfg.latency(2)));
        imagesc(stat.time, stat.freq, plotdata, 'CDataMapping', 'scaled');
        set(gca,'YDir','normal');
        colormap(jet); colorbar('location', 'EastOutside');
        set(gca,'clim', [-max(max(abs(plotdata))) max(max(abs(plotdata)))]);
        title('% change from baseline'); 

        subplot(4,1,2);
        plotdata = squeeze(stat.stat);
        imagesc(stat.time, stat.freq, plotdata , 'CDataMapping', 'scaled');
        set(gca,'YDir','normal'); 
        colormap(jet); colorbar('location', 'EastOutside');
        set(gca,'clim', [-max(max(abs(plotdata))) max(max(abs(plotdata)))]);
        title('t-values'); 

        subplot(4,1,3);
        plotdata = squeeze(stat.stat.*stat.mask);
        imagesc(stat.time, stat.freq, plotdata, 'CDataMapping', 'scaled');
        set(gca,'YDir','normal');
        colormap(jet); colorbar('location', 'EastOutside');    
        crange = [-max(max(abs(plotdata))) max(max(abs(plotdata)))];
        set(gca,'clim', crange);
        title('thresholded t-values (p<0.05 uncorrected)'); 

        subplot(4,1,4);
        plotdata = squeeze(stat.stat.*stat.mask_extthresh);
        imagesc(stat.time, stat.freq, plotdata, 'CDataMapping', 'scaled');
        set(gca,'YDir','normal');
        colormap(jet); colorbar('location', 'EastOutside');    
        crange = [-max(max(abs(plotdata))) max(max(abs(plotdata)))];
        set(gca,'clim', crange);
        title(['thresholded t-values (p<0.05 cluster extent threshold = ' num2str(k) ' corrected)']);     


    %     subplot(4,1,4);
    %     plotdata = squeeze(stat_cor.stat.*stat_cor.mask);
    %     imagesc(stat_cor.time, stat_cor.freq, plotdata, 'CDataMapping', 'scaled');
    %     set(gca,'YDir','normal');
    %     colormap(jet); colorbar('location', 'EastOutside');    
    %     set(gca,'clim', crange);
    %     title('thresholded t-values (p<0.05 cluster-corrected)'); 

        % save plot
        if opt.saveplot
            set(gcf, 'Name',[gcfg.conditionLabels{gcfg.thisCondition} '_GA_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_tsv_flag '_' SO_timelockchannel '_stats_' gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
            filename = fullfile(gcfg.paths.results, 'group', [gcfg.conditionLabels{gcfg.thisCondition} '_GA_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_tsv_flag '_' SO_timelockchannel '_stats_' gcfg.timestamp]);
            saveas(gcf,[filename '.fig']);        
            set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [25.4 19.05], 'PaperPosition', [-3 -1 30 19.5]);
            print(gcf, '-dpdf', '-r300', [filename '.pdf']);   
        end
    CFC_results_SO.stat = stat;    
    end

end 

%% phase amplitude coupling per trial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opt.pacpt % phase amplitude coupling per trial

    subjectNames = gcfg.subjectNames;

    for i = 1:numel(gcfg.subjects) % loop over sujects
       
        if isfield(gcfg,'freqwinIndividual')
            gcfg.freqwin = gcfg.freqwinIndividual(gcfg.subjects(i),:);
        end
        
        % load trialwise data of subject i
        load(fullfile(gcfg.paths.results,subjectNames{i},[subjectNames{i} '_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_targetchannel{gcfg.subjects(i)} '_eventdata.mat']), 'SO_eventdata'); % load freq data
        load(fullfile(gcfg.paths.results,subjectNames{i},[subjectNames{i} '_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_targetchannel{gcfg.subjects(i)} '_freq.mat']), 'SO_freq'); % load freq data

        % select trials        
        cfg = [];
        cfg.trials = find(SO_tsv{i});
        cfg.tolerance   = 1e-3;
        SO_eventdata = ft_selectdata(cfg,SO_eventdata); 
        SO_freq = ft_selectdata(cfg,SO_freq);         
        
        % define data 
        freqindex = SO_freq.freq >= gcfg.freqwin(1) & SO_freq.freq <= gcfg.freqwin(2);     
        timeindex = unique(int16(find(SO_eventdata.time{1} >= gcfg.timewin(1) & SO_eventdata.time{1} <= gcfg.timewin(2))-1));
    
        trialN{i} = size(SO_eventdata.trial,2); % # trials        
        for j = 1:trialN{i} % loop over trials
            
            % define data for trial j
            LF = SO_eventdata.trial{j};
            HF = squeeze(mean(SO_freq.powspctrm(j,1,freqindex,:),3))';
            
            % preprocess LF data
            LFhp = ft_preproc_bandpassfilter(LF, double(SO_eventdata.fsample), [gcfg.hpfilter gcfg.lpfilter ], 3*fix(SO_eventdata.fsample/gcfg.hpfilter)+1, 'fir', 'twopass'); % Mathilde's formula
            
            % preprocess HF data
%             usfactor = double(SO_eventdata.fsample / round(1/mean(diff(SO_freq.time)))); %  upsampling factor to interpolate HF (TFR) signal to LF signal
            HF(isnan(HF)) = 0; % replace NaN by 0
%             HFus = interp(HF,usfactor); % upsample to fit LF time resolution
            HFus = transpose(resample(transpose(HF),SO_eventdata.fsample,round(1/mean(diff(SO_freq.time))))); % upsample to fit LF time resolution
            HFhp = ft_preproc_highpassfilter(HFus, double(SO_eventdata.fsample), gcfg.hpfilter, 3*fix(SO_eventdata.fsample/gcfg.hpfilter)+1, 'fir', 'twopass'); % Mathilde's formula        
            
            % LF event-locked HF power change (difference)
            bl = [-2.5 -1.5]; % for spindles it would be bl = [-2.0 -1.0]
            blindex = unique(int16(find(SO_eventdata.time{1} >= bl(1) & SO_eventdata.time{1} <= bl(2))));
            ERPow{i}(j) = mean(HFus(:,timeindex),2) - mean(HFus(:,blindex),2);
            blPow{i}(j) = mean(HFus(:,blindex),2);
            
            % get phase angles (hilbert transform)
            LFph = angle(hilbert(LFhp(1,:))); 
            HFph = angle(hilbert(HFhp(1,:))); 
        
            % synchronization index (cf. Cohen et al., J Neurosci Methods 2008)
            phi1{i}(j,:) = LFph(:,timeindex); % phase values lf modulating signal
            phi2{i}(j,:) = HFph(:,timeindex); % phase values hf power modulation 
            SI{i}(j) = mean(exp(1i*[phi1{i}(j,:) - phi2{i}(j,:)])); % Synchronization Index (SI)
            SIm{i}(j) = abs(SI{i}(j)); % SI magnitude (0 = perfect desychronization to 1 = perfect synchronization)
            SIp{i}(j) = (atan2(imag(SI{i}(j)),real(SI{i}(j)))); % preferred phase of the synchronization            
                        
        end % of loop over trials
        
        
        % calculate mean LF event-locked HF power change (percent change!)
        AERPow{i} = mean(ERPow{i})/mean(blPow{i}) * 100;        
                
        % two-sided one-sample t-test of percent LF event-timelocked HF power change against zero
        [AERPowh{i},AERPowp{i},AERPowci{i},AERPowstats{i}] = ttest([ERPow{i}(:)], 0, 0.05, 'both');
        AERPowt{i} = AERPowstats{i}.tstat;
        
        
        % calculate mean vector length from preferred phase angles  (measure of consistency of SI preferred phase over LF events)
        r{i} = circ_r(SIp{i}');
            
        % calculate mean SI magnitude (measure of average SI magnitude over LF events, independent of their respective phases)
        ASIm{i} = mean(SIm{i});       

        % calculate mean preferred phase (measure of avreage phase angle of preferred phase over LF events)
        ASIp{i} = circ_mean(SIp{i}');
        ASIp_ang{i} = circ_rad2ang(ASIp{i});
        
        % generate surrogates (from uniform distribution)
        for k = 1:10000 % iterations for surrogate distribution            
            Sur_SIp{i}(k,:) = rand(1,trialN{i})*2*pi; % generate as many random 2pi values from uniform distribution as there are trials
            Sur_r(i,k) = circ_r(Sur_SIp{i}(k,:)'); % calculate mean vector length of surrogates       
            Sur_ASIp(i,k) = circ_mean(Sur_SIp{i}(k,:)'); % calculate mean preferred phase for surrogates
        end     
        
        % single-subject permutation test of mean vector length against surrogates from uniform distribution with same trial number        
        ixr = find(sort(Sur_r(i,:)) > r{i}, 1);        
        if isempty(ixr), ixr = length(Sur_r(i,:)); end           
        Sur_r_p{i} = 1-(ixr/length(Sur_r(i,:)));        
        z{i} = (r{i} - mean(Sur_r(i,:))) / std(Sur_r(i,:)); % z-transform actual r-value with respect to surrogate distribution
        z_p{i} = normcdf(-abs(z{i}),0,1); % transform z-value to p-value                       
        
        % Rayleigh test on non-uniform distribution of preferred phase angles within subject
        [Rp{i},Rz{i}] = circ_rtest([SIp{i}]');  
 
        % single-subject chiµ test on non-uniform distribution
%         [Chi_h{i},Chi_p{i},Chi_st{i}] = chi2gof(SIp{i},'edges',linspace(-pi,pi,18),'expected',length(SIp{1})*diff(linspace(-pi,pi,18)));
        
        % single-subject KS test on non-uniform distribution
%         [KSh{i},KSp{i},KSk{i},KSc{i}] = kstest([SIp{i}]',[[SIp{i}]', unifcdf([SIp{i}]',-pi,pi)]);
 

        % print SUBJECT results
        display(' '); display(['SUBJECT ' subjectNames{i} ' RESULTS']);        
        display('LF Event-Locked HF Power Change');       
        display([subjectNames{i} ' AERPow: ' num2str(AERPow{i})]); 
        display(sprintf('%s %s%0.10f', subjectNames{i}, ['T-test t(' num2str(AERPowstats{i}.df) ') (AERPowt): '], AERPowt{i}));
        display(sprintf('%s %s%0.10f', subjectNames{i}, 'T-test p (AERPowp): ', AERPowp{i}));
        display('PAC Strength');
        display([subjectNames{i} ' mean vecotr length (r): ' num2str(r{i})]); 
        display([subjectNames{i} ' magnitude of synchronization index (ASIm): ' num2str(ASIm{i})]); 
        display(sprintf('%s %s%0.10f', subjectNames{i} , 'Rayleigh''s z (Rz): ', Rz{i}));
        display(sprintf('%s %s%0.10f', subjectNames{i} , 'Rayleigh p (Rp): ', Rp{i}));  
        display(sprintf('%s %s%0.10f', subjectNames{i} , 'Mean vector length permutation p (Sur_r_p{i}): ', Sur_r_p{i}));
        display(sprintf('%s %s%0.10f', subjectNames{i} , 'Mean vector length permutation z (z{i}): ', z{i}));
        display(sprintf('%s %s%0.10f', subjectNames{i} , 'Mean vector length permutation z-to-p (z_p{i}): ', z_p{i}));                
        display('PAC Preferred Phase');
        display([subjectNames{i} ' preferred phase (ASIp): ' num2str(ASIp{i})]); 
        
        % plot circular stats & surrogates
        figure; set(gcf,'name',subjectNames{i});        
        subplot(2,4,1); circ_plot([SIp{i}]','pretty','o',true,'linewidth',2,'color','r');        
        title({['distribution of preferred phase angles of ' num2str(trialN{i}) ' trials'], ['mean angle: ' num2str(ASIp{i}) ' rad, ' num2str(circ_rad2ang(ASIp{i})) 'degree']});
        subplot(2,4,2); circ_plot([SIp{i}]','hist',[],18,true,true,'linewidth',2,'color','r');
        title({['distribution of preferred phase angles of ' num2str(trialN{i}) ' trials'], ['mean angle: ' num2str(ASIp{i}) ' rad, ' num2str(circ_rad2ang(ASIp{i})) 'degree']});
        subplot(2,4,3); circ_plot([Sur_ASIp(i,:)]','pretty','o',true,'linewidth',2,'color','r');
        title({'surrogate distribution of preferred phase angles'});
        subplot(2,4,4); circ_plot([Sur_ASIp(i,:)]','hist',[],18,true,true,'linewidth',2,'color','r');
        title({'surrogate distribution of preferred phase angles'});
        subplot(2,4,5:8); hist([Sur_r(i,:)]',round(k/100));
        hold on; ylim = get(gca,'ylim'); line([r{i} r{i}], ylim,'color','r'); hold off;
        title(['surrogate distribution of mean vector lengths with actual r-value ' num2str(r{i}) ' (red)']);            
        
        % save figure        
        set(gcf, 'Name',[subjectNames{i} '_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_tsv_flag '_' gcfg.SO_timelockchannel '_CFCtrialwise_' sprintf('[%2.2f %2.2fs]',gcfg.timewin(1),gcfg.timewin(2)) '_' sprintf('[%2.2f %2.2fHz]',gcfg.freqwin(1),gcfg.freqwin(2)) '_'  gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
        filename = fullfile(gcfg.paths.results, subjectNames{i}, [subjectNames{i} '_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_tsv_flag '_' gcfg.SO_timelockchannel '_CFCtrialwise_' sprintf('[%2.2f %2.2fs]',gcfg.timewin(1),gcfg.timewin(2)) '_' sprintf('[%2.2f %2.2fHz]',gcfg.freqwin(1),gcfg.freqwin(2)) '_' gcfg.timestamp]);
        saveas(gcf,[filename '.fig']);   
        set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [25.4 19.05], 'PaperPosition', [-3 -1 30 19.5]);
        print(gcf, '-dpdf', '-r300', [filename '.pdf']);   
        
    end % of loop over subjects           

    
    % calculate grand average LF event-locked HF percent power change 
    GAERPow = mean([AERPow{:}]);
    SEMAERPow = std([AERPow{:}])/sqrt(length(AERPow));
    
    % two-sided one-sample t-test of percent LF event-timelocked HF power change against zero
    [GAERPowh,GAERPowp,GAERPowci,GAERPowstats] = ttest([AERPow{:}], 0, 0.05, 'both'); 
    GAERPowt = GAERPowstats.tstat;
    
    % calculate grand mean direction of preferred phase 
    GASIp = circ_mean([ASIp{:}]');
    GASIp_ang = circ_rad2ang(GASIp);
    [SASIp sASIp] = circ_var([ASIp{:}]'); 
    SASIp_ang = circ_rad2ang(SASIp); % circular variance 1-r
    sASIp_ang = circ_rad2ang(sASIp); % angular variance 2*(1-r)     
    
    % Rayleigh test on non-uniform distribution of preferred phase angles over subjects
    [GRp,GRz] = circ_rtest([ASIp{:}]');    

	% V test on non-uniform distribution of preferred phase angles over subjects    
    if ~isempty(ft_channelselection({'eeg'},{gcfg.SO_timelockchannel})) % if EEG channel
        Vtestvalue = 0;
    elseif strcmp(gcfg.SO_timelockchannel,'HC')
        Vtestvalue = 180;
    end
    
    [GVp, GV] = circ_vtest([ASIp{:}]', circ_ang2rad(Vtestvalue)); % test explicitly against a predefined direction 
    
    % generate surrogates (from uniform distribution)             
	GSur_r = mean(Sur_r,1);
    
    % permutation test against surrogates from uniform distribution with same trial number
    Gr = mean([r{:}]);
    SEMGr = std([r{:}])/sqrt(length(r));
    ixGr = find(sort(GSur_r) > Gr, 1);
    if isempty(ixGr), ixGr = length(GSur_r); end           
    GSur_r_p = 1-(ixGr/length(GSur_r));        
    Gz = (Gr - mean(GSur_r)) / std(GSur_r); % z-transform actual r-value with respect to surrogate distribution
    Gz_p = normcdf(-abs(Gz),0,1); % transform z-value to p-value       

    % save command window output from here
    filename = fullfile(gcfg.paths.results, 'group', [gcfg.conditionLabels{gcfg.thisCondition} '_GA_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_tsv_flag '_' gcfg.SO_timelockchannel '_CFCtrialwise_' sprintf('[%2.2f %2.2fs]',gcfg.timewin(1),gcfg.timewin(2)) '_' sprintf('[%2.2f %2.2fHz]',gcfg.freqwin(1),gcfg.freqwin(2)) '_' gcfg.timestamp '.txt']);
    diary(filename); % file to save command window output to filename
    
	% print SUBJECT SUMMARY results
    subjectssummaryresults = [gcfg.subjects; [trialN{:}]; 
                        [AERPow{:}]; [AERPowt{:}]; [AERPowp{:}];
                        [[r{:}]; [ASIm{:}]; [Rz{:}]; Rp{:}]; [z{:}]; [z_p{:}];
                        [ASIp{:}]; circ_rad2ang([ASIp{:}])];
    display(' '); display('SUBJECT SUMMARY RESULTS');
    display(num2str(subjectssummaryresults));
    
    % print GROUP results
    display(' '); display('GROUP RESULTS');
    display('LF Event-Locked HF Power Change');
	display(sprintf('%s%0.10f', 'Mean % power change (GAERPow): ', GAERPow));  
    display(sprintf('%s%0.10f', ['T-test t(' num2str(GAERPowstats.df) ') (GAERPowt): '], GAERPowt));    
    display(sprintf('%s%0.10f', 'T-test p (GAERPowp): ', GAERPowp));           
    display('PAC Strength');
    display(sprintf('%s%0.4f', 'Mean vector length (Gr): ', Gr));
    display(sprintf('%s%0.4f', 'Mean vector length permutation z (Gz): ', Gz));
    display(sprintf('%s%0.10f', 'Mean vector length permutation p (GSur_r_p): ', GSur_r_p));       
    display(sprintf('%s%0.10f', 'Mean vector length permutation z-to-p (Gz_p): ', Gz_p));           
    display('PAC Preferred Phase');
    display(sprintf('%s%0.10f', 'Grand mean of preferred phase (GASIp_ang): ', GASIp_ang));   
    display(sprintf('%s%0.10f', 'Angular variance of preferred phase (sASIp_ang): ', sASIp_ang));             
    display(sprintf('%s%0.10f', 'Rayleigh''s p (GRp): ', GRp));
    display(sprintf('%s%0.10f', 'Rayleigh''s z (GRz): ', GRz));
    display(sprintf('%s%u%s%0.10f', 'V Test (GVp) against ', Vtestvalue, 'µ :', GVp));
    display(sprintf('%s%u%s%0.10f', 'V Test (GV) against ', Vtestvalue, 'µ :', GV));
    
    % plot circular GROUP stats
    figure; 
    subplot(1,2,1); circ_plot([ASIp{:}]','pretty','o',true,'linewidth',2,'color','r');
    title({'distribution of individual preferred phase angles', ['mean angle: ' num2str(GASIp) ' rad, ' num2str(circ_rad2ang(GASIp)) 'degree']});
    subplot(1,2,2); hist(GSur_r,round(k/100));
    hold on; ylim = get(gca,'ylim'); line([Gr Gr], ylim,'color','r'); hold off;
    title({'surrogate distribution of group mean vector lengths',['with actual group r-value ' num2str(Gr) ' (red)']});   
    
    % save figure
    set(gcf, 'Name',[gcfg.conditionLabels{gcfg.thisCondition} '_GA_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_tsv_flag '_' gcfg.SO_timelockchannel '_CFCtrialwise_' sprintf('[%2.2f %2.2fs]',gcfg.timewin(1),gcfg.timewin(2)) '_' sprintf('[%2.2f %2.2fHz]',gcfg.freqwin(1),gcfg.freqwin(2)) '_' gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
    filename = fullfile(gcfg.paths.results, 'group', [gcfg.conditionLabels{gcfg.thisCondition} '_GA_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_tsv_flag '_' gcfg.SO_timelockchannel '_CFCtrialwise_' sprintf('[%2.2f %2.2fs]',gcfg.timewin(1),gcfg.timewin(2)) '_' sprintf('[%2.2f %2.2fHz]',gcfg.freqwin(1),gcfg.freqwin(2)) '_' gcfg.timestamp]);
    saveas(gcf,[filename '.fig']);   
    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [25.4 19.05], 'PaperPosition', [-3 -1 30 19.5]);
    print(gcf, '-dpdf', '-r300', [filename '.pdf']); 

        
    % summarize results for function output   
    subjects = gcfg.subjects;
    varList = {'subjects', 'trialN', 'AERPow', 'AERPowt', 'AERPowp', 'r','SIp', 'SIm', 'ASIm','Rz','Rp','z','z_p','ASIp','ASIp_ang',...
                'GAERPow', 'GAERPowt', 'GAERPowp', 'Gr', 'Gz', 'GSur_r_p', 'Gz_p', 'GASIp_ang', 'sASIp_ang', 'GRp', 'GRz', 'GVp', 'GV'};   
	for index = 1:numel(varList)
        results.(varList{index}) = eval(varList{index}); 
    end

    % save results !!! --> MOVED TO *_sl !!!
%     filename = fullfile(gcfg.paths.results, 'group', [gcfg.conditionLabels{gcfg.thisCondition} '_GA_SO_' gcfg.SOspec '_' gcfg.SO_timelockevent '_' gcfg.SO_tsv_flag '_' gcfg.SO_timelockchannel '_CFCtrialwise_' sprintf('[%2.2f %2.2fs]',gcfg.timewin(1),gcfg.timewin(2)) '_' sprintf('[%2.2f %2.2fHz]',gcfg.freqwin(1),gcfg.freqwin(2)) '_' gcfg.timestamp '.mat']);
%     save('-v7.3', filename, varList{:});  
    
    % stop saving command window output
    diary off;
    
    
end % of phase amplitude coupling per trial


end % of function