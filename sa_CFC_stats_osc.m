%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sa_CFC_stats_osc
% by Til Ole Bergmann 2013
% last modified 2017/01/19 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% performs cross-frequency coupling (CFC) , more specifically phase-amplitude coupling (PAC) 
%
% output arguments:
% results = CFC results structure
% results.figurehandle =  handle to figure
%
% input arguments:
% GA_osc_freq = Grand Average (with keeptrials = no) freq structure from fieldtrip time-frequency analysis (TFR) 
% GA_osc_timelock = Grand Average (with keeptrials = no) timelock structure  from fieldtrip timelock analysis 
% All_osc_freq  = 1xN cell array containing for N subject the normalized (e.g. % change from baseline) single-trial (keeptrials = yes) freq structures from fieldtrip time-frequencyanalysis (TFR) 
% All_osc_raw_freq = 1xN cell array containing for N subject the raw (unnormalized) single-trial (keeptrials = yes) freq structures from fieldtrip time-frequencyanalysis (TFR) 
% All_osc_timelock = 1xN cell array containing for N subject the single-trial (keeptrials = yes) timelock structures from fieldtrip timelocke analysis 

% cfg = [];
% cfg.subjectNames = 1xN cell array with the names of the subjects 
% cfg.subjects = 1xN cell array with the subject numbers 
%%%%%%% cfg.channel = Nx1 cell array with selection of channels (default = 'all'), see ft_channelselection for details
%%%%%%% cfg.layout = string indicating layout file, e.g., 'easycapM3.mat'
% cfg.hpfilter = frequency for high-pass filter for phase-providing frequency specified in Hz (no default)
% cfg.lpfilter = frequency for low-pass filter for phase-providing frequency specified in Hz (no default)
% cfg.timewin = time window around zero to be considered for phase-amplitude coupling, specified as [low high] in s (no default); 
% cfg.freqwin = frequency window to be considere for phase-amplitude coupling, specified as [low high] im Hz (no default);
% cfg.Vtestvalue = phase angle for which maximum power is expected in power-providing frequency, specified in angular degree (e.g., 0 = peak , 90 = falling flank, 180 = trough, 270 = rising flank) (no default)
% cfg.osc_timelockchannel = string indicating the channel for which CFC shall be calcualted, ERPs time-locked and TFR plotted (no default)
% [CFC_results_osc] = sa_CFC_stats_osc(cfg, GA_osc_freq, All_osc_freq, All_osc_timelock);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results] = sa_CFC_stats_osc(gcfg, GA_osc_freq, All_osc_freq, All_osc_timelock)
keyboard;
results = struct; % initialize results structure

%% check fields
if ~isfield(gcfg,'thisCondition')
    gcfg.thisCondition = 1;
end    
if ~isfield(gcfg,'conditionLabels')
    if ~isfield(gcfg,'conditions')
        gcfg.conditionLabels = {''};
    elseif isfield(gcfg,'conditions')
        gcfg.conditionLabels = repmat({''},size(gcfg.conditions));
    end
end    

%% options
opt.tfrstats = 0; 
opt.tfrplot = 0;
opt.clusterstats = 0;
opt.clustextthresh = 0;
opt.plot = 0;
opt.saveplot = 0;

surrogate_type = 'fft_shuffle_phase'; % 'shuffle_phase', 'cyclic_shift_phase', 'random_shift_phase', 'fft_shuffle_phase' 

opt.ptd = 0; % peak trough differences
opt.pacpt = 1; % phase amplitude coupling per trial

%% TFR non-parametric t-statistics against zero
if opt.tfrstats  
   
    % create zero data for comparison
	All_zero_freq = All_osc_freq;
    for i = 1:length(All_osc_freq) % create fake zero data
        All_zero_freq{i}.powspctrm(~isnan(All_zero_freq{i}.powspctrm)) = 0;
    end
    
    % non-parametric t-statistcs
    cfg = [];
    cfg.channel = 'all';
    cfg.latency = gcfg.plot_osc_xlim; % [begin end] in seconds; 'all'
    cfg.trials = 'all';
    cfg.frequency = gcfg.plot_osc_ylim; % [begin end]; 'all'
    cfg.parameter = 'powspctrm';   
    cfg.method = 'montecarlo';
    cfg.statistic = 'depsamplesT'; % dependent samples T-statistic	
    cfg.design(1,:) = [repmat(1:numel(All_osc_freq),1,2)]; % subject vector        
    cfg.design(2,:) = [ones(1,length(All_osc_freq)) ones(1,length(All_osc_freq))*2]; % condition vector        
    cfg.uvar = 1; % number or list with indices, unit variable(s)
    cfg.ivar = 2; % number or list with indices, independent variable(s)   
    cfg.numrandomization = 'all'; % 'all'
    cfg.correctm = 'no'; % 'no' , 'max', cluster', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')
    cfg.alpha = 0.05; 
    cfg.tail = 0; %, -1, 1 or 0 (default = 0)   
    cfg.correcttail = 'prob'; % correct p-values or alpha-values when doing a two-sided test, 'alpha','prob' or 'no' (default = 'no')
    cfg.feedback = 'text'; % 'gui', 'text', 'textbar' or 'no' (default = 'textbar')    
    cfg.randomseed = 'yes'; % string, 'yes', 'no' or a number (default = 'yes')       
    [stat] = ft_freqstatistics(cfg, All_osc_freq{:}, All_zero_freq{:}); 
%     [stat] = ft_freqstatistics(cfg, All_osc_freq{:}, All_osc_surrogate_freq{:}); 
%     [stat] = ft_freqstatistics(cfg, All_osc_freq{:}, All_osc_free_freq{:}); 

    results.stat = stat;

end


%% cluster statistics
if opt.clusterstats
    cfg.correctm = 'cluster'; % 'no' , 'max', cluster', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')
    cfg.clusterstatistic = 'maxsum'; % how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
    cfg.clusterthreshold = 'parametric'; % method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
    cfg.clusteralpha = 0.05; % for either parametric or nonparametric thresholding per tail (default = 0.05)
    cfg.clustertail = 0; % -1, 1 or 0 (default = 0)        
    [cluststat] = ft_freqstatistics(cfg, All_osc_freq{:}, All_zero_freq{:});     
%     [cluststat] = ft_freqstatistics(cfg, All_osc_freq{:}, All_osc_surrogate_freq{:}); 
%     [cluststat] = ft_freqstatistics(cfg, All_osc_freq{:}, All_osc_free_freq{:}); 

    results.cluststat = cluststat;

end


%% cluster extent threshold (cet) correction
  if opt.clustextthresh    
    k = 200; % cluster size
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
    results.stat.cetmask = permute(M, [3,1,2]);
end
       

%% plot stats
if opt.tfrplot    
    results.figurehandles.tfrplot = figure;
    subplot(4,1,1);
    plotdata = squeeze(GA_osc_freq.powspctrm(1,GA_osc_freq.freq >= cfg.frequency(1) & GA_osc_freq.freq <= cfg.frequency(2), GA_osc_freq.time >= cfg.latency(1) & GA_osc_freq.time <= cfg.latency(2)));
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
%     set(gca,'clim', [-max(max(abs(plotdata))) max(max(abs(plotdata)))]);
    set(gca,'clim', [min([-max(max(abs(plotdata))), -1e-10]) max([max(max(abs(plotdata))), 1e-10])]) % to make sure it is not [0 0]
    title('thresholded t-values (p<0.05 uncorrected)'); 

	subplot(4,1,4);
    plotdata = squeeze(stat.stat.*stat.cetmask);
    imagesc(stat.time, stat.freq, plotdata, 'CDataMapping', 'scaled');
    set(gca,'YDir','normal');
    colormap(jet); colorbar('location', 'EastOutside');    
    crange = [-max(max(abs(plotdata))) max(max(abs(plotdata)))];
    set(gca,'clim', crange);
    title(['thresholded t-values (p<0.05 cluster extent threshold = ' num2str(k) ' corrected)']);        
    
end



%% analysis of peak-trough differences (ptd)
if opt.ptd 
    
    % settings    
    peaks = {}; troughs = {}; newpeaks = {}; newtroughs = {};
    fsample = round(1/mean(diff(All_osc_timelock{1}.time))); % timelock time resolution
    tfrfsample = round(1/mean(diff(All_osc_freq{1}.time))); % tfr time resolution
	dsfactor = fsample/tfrfsample; % downsample by factor
    searchwin = find(All_osc_timelock{1}.time >= gcfg.timewin(1) & All_osc_timelock{1}.time <= gcfg.timewin(2));
	newsearchwin = unique(int16((searchwin-1)/dsfactor));                
    
	for i = 1:numel(All_osc_timelock) % loop over sujects
        
        % high-pass filter spindel timelocked individual averages
        All_osc_timelock_hpf{i}.avg = ft_preproc_highpassfilter(All_osc_timelock{i}.avg, double(fsample), gcfg.hpfilter , 3*fix(fsample/gcfg.hpfilter )+1, 'fir', 'twopass'); % Mathilde's formula
    
        % detect individual peaks and zero crossings automatically
        signal{i} = All_osc_timelock_hpf{i}.avg;
        signs{i} = sign(signal{i});
        diffs{i} = diff(signs{i});
        p2nzc{i} = find(diffs{i} < 0);
        n2pzc{i} = find(diffs{i} > 0);
        if p2nzc{i}(1) > n2pzc{i}(1) 
            pa = 1; na = 0;
        elseif p2nzc{i}(1) < n2pzc{i}(1) 
            pa = 0; na = 1;
        end
        for j = 1:size(n2pzc{i},2)-1 % loop over peaks
            [maxval, maxpos] = max(signal{i}(:,n2pzc{i}(j):p2nzc{i}(j+na)),[],2);
            peaks{i}(j) = n2pzc{i}(j) + maxpos - 1;
        end
        for j = 1:size(p2nzc{i},2)-1 % loop over troughs
            [minval, minpos] = min(signal{i}(:,p2nzc{i}(j):n2pzc{i}(j+pa)),[],2);
            troughs{i}(j) = p2nzc{i}(j) + minpos - 1;
        end
        peaks{i} = peaks{i}(:,peaks{i} > searchwin(1) & peaks{i} < searchwin(end));        
        troughs{i} = troughs{i}(:,troughs{i} > searchwin(1) & troughs{i} < searchwin(end));
        p2nzc{i} = p2nzc{i}(:,p2nzc{i} > searchwin(1) & p2nzc{i} < searchwin(end));
        n2pzc{i} = n2pzc{i}(:,n2pzc{i} > searchwin(1) & n2pzc{i} < searchwin(end));        
        
        % set interval around trials to extract
        troughwidth = [-0.01 0.01]*fsample; % in s 
        peakwidth = [-0.01 0.01]*fsample; % in s 
        newtroughs{i} = [];
        newpeaks{i} = [];
        for k = 1:size(troughs{i},2)
            newtroughs{i} = [newtroughs{i} troughs{i}(k)+troughwidth(1):troughs{i}(k)+troughwidth(2)];
        end
        for k = 1:size(peaks{i},2)
            newpeaks{i} = [newpeaks{i} peaks{i}(k)+peakwidth(1):peaks{i}(k)+peakwidth(2)];
        end
        newp2nzc = p2nzc;        
        newtroughs = troughs;
        newn2pzc = n2pzc;
        newpeaks = peaks;
        
        % extract power values from individual TFRs
        freqlim = All_osc_freq{i}.freq >= gcfg.freqwin(1) & All_osc_freq{i}.freq <= gcfg.freqwin(2); 
        [N p2nzcbin{i}] = histc(All_osc_timelock{i}.time(newp2nzc{i}), All_osc_freq{i}.time);        
        [N troughbin{i}] = histc(All_osc_timelock{i}.time(newtroughs{i}), All_osc_freq{i}.time); 
        [N n2pzcbin{i}] = histc(All_osc_timelock{i}.time(newpeaks{i}), All_osc_freq{i}.time);         
        [N peakbin{i}] = histc(All_osc_timelock{i}.time(newn2pzc{i}), All_osc_freq{i}.time); 
        avgp2nzcpow{i} = mean(mean(All_osc_freq{i}.powspctrm(:,freqlim,unique(p2nzcbin{i})),3),2);                      
        avgtroughpow{i} = mean(mean(All_osc_freq{i}.powspctrm(:,freqlim,unique(troughbin{i})),3),2);       
        avgn2pzcpow{i} = mean(mean(All_osc_freq{i}.powspctrm(:,freqlim,unique(n2pzcbin{i})),3),2);
        avgpeakpow{i} = mean(mean(All_osc_freq{i}.powspctrm(:,freqlim,unique(peakbin{i})),3),2);
    end
    
%     tn = 0; % test number
%     
%     tn = tn+1; % two-sided test of peak and trough values against each other
%     [h{tn},p{tn},ci{tn},stats{tn}] = ttest([avgpeakpow{:}],[avgtroughpow{:}],0.05,'both'); % t-test
%     p{tn}
    
%     % plot
%     dat = [avgp2nzcpow{:}; avgtroughpow{:}; avgn2pzcpow{:}; avgpeakpow{:}];
%     h = figure;
%     set(h,'DefaultAxesColorOrder',varycolor(13)); hold on;
%     plot(dat, 'linewidth', 2);
%     plot(mean(dat,2),'k','linewidth',3);
%     legend(cellstr(num2str(gcfg.subjects(:)))','Location','EastOutside');
%     
end


%% phase amplitude coupling per trial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opt.pacpt % phase amplitude coupling per trial  
    
    for i = 1:numel(gcfg.subjects) % loop over sujects
       
        % load trialwise data of subject i
        load(gcfg.files.eventdata{i}, 'osc_eventdata'); % load freq data;
        load(gcfg.files.freq{i}, 'osc_freq'); % load freq data

        % select trials        
        cfg = [];
        cfg.trials = gcfg.trials{i};
        cfg.tolerance   = 1e-3;
        osc_eventdata = ft_selectdata(cfg,osc_eventdata); 
        osc_freq = ft_selectdata(cfg,osc_freq);
        
        % define data 
        freqindex = osc_freq.freq >= gcfg.freqwin(1) & osc_freq.freq <= gcfg.freqwin(2);     
        timeindex = unique(int16(find(osc_eventdata.time{1} >= gcfg.timewin(1) & osc_eventdata.time{1} <= gcfg.timewin(2))));
    
        trialN{i} = size(osc_eventdata.trial,2); % # trials        
        for j = 1:trialN{i} % loop over trials
            
            % define data for trial j
            LF = osc_eventdata.trial{j};
            HF = squeeze(mean(osc_freq.powspctrm(j,1,freqindex,:),3))';
            
            % preprocess LF data
            LFhp = ft_preproc_bandpassfilter(LF, double(osc_eventdata.fsample), [gcfg.hpfilter  gcfg.lpfilter ], 3*fix(osc_eventdata.fsample/gcfg.hpfilter )+1, 'fir', 'twopass'); % Mathilde's formula
            
            % preprocess HF data
%             usfactor = double(osc_eventdata.fsample / round(1/mean(diff(osc_freq.time)))); %  upsampling factor to interpolate HF (TFR) signal to LF signal
            HF(isnan(HF)) = 0; % replace NaN by 0
%             HFus = interp(HF,usfactor); % upsample to fit LF time resolution
%             HFus = transpose(resample(transpose(HF),osc_eventdata.fsample,round(1/mean(diff(osc_freq.time))))); % upsample to fit LF time resolution
            HFus = transpose(resample(transpose(double(HF)), double(osc_freq.time), double(osc_eventdata.fsample), 'pchip')); % upsample to fit LF time resolution
            HFhp = ft_preproc_highpassfilter(HFus, double(osc_eventdata.fsample), gcfg.hpfilter , 3*fix(osc_eventdata.fsample/gcfg.hpfilter )+1, 'fir', 'twopass'); % Mathilde's formula        
            
            % LF event-locked HF power change (difference)
            bl = [-2.0 -1.0]; % for SO it would be bl = [-2.5 -1.5]
            blindex = unique(int16(find(osc_eventdata.time{1} >= bl(1) & osc_eventdata.time{1} <= bl(2))-1));
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
 
        % single-subject chi° test on non-uniform distribution
%         [Chi_h{i},Chi_p{i},Chi_st{i}] = chi2gof(SIp{i},'edges',linspace(-pi,pi,18),'expected',length(SIp{1})*diff(linspace(-pi,pi,18)));
        
        % single-subject KS test on non-uniform distribution
%         [KSh{i},KSp{i},KSk{i},KSc{i}] = kstest([SIp{i}]',[[SIp{i}]', unifcdf([SIp{i}]',-pi,pi)]);
 
        % print SUBJECT results
        display(' '); display(['SUBJECT ' gcfg.subjectNames{i} ' RESULTS']);        
        display('LF Event-Locked HF Power Change');       
        display([gcfg.subjectNames{i} ' AERPow: ' num2str(AERPow{i})]); 
        display(sprintf('%s %s%0.10f', gcfg.subjectNames{i}, ['T-test t(' num2str(AERPowstats{i}.df) ') (AERPowt): '], AERPowt{i}));
        display(sprintf('%s %s%0.10f', gcfg.subjectNames{i}, 'T-test p (AERPowp): ', AERPowp{i}));
        display('PAC Strength');
        display([gcfg.subjectNames{i} ' mean vecotr length (r): ' num2str(r{i})]); 
        display([gcfg.subjectNames{i} ' magnitude of synchronization index (ASIm): ' num2str(ASIm{i})]); 
        display(sprintf('%s %s%0.10f', gcfg.subjectNames{i} , 'Rayleigh''s z (Rz): ', Rz{i}));
        display(sprintf('%s %s%0.10f', gcfg.subjectNames{i} , 'Rayleigh p (Rp): ', Rp{i}));  
        display(sprintf('%s %s%0.10f', gcfg.subjectNames{i} , 'Mean vector length permutation p (Sur_r_p{i}): ', Sur_r_p{i}));
        display(sprintf('%s %s%0.10f', gcfg.subjectNames{i} , 'Mean vector length permutation z (z{i}): ', z{i}));
        display(sprintf('%s %s%0.10f', gcfg.subjectNames{i} , 'Mean vector length permutation z-to-p (z_p{i}): ', z_p{i}));                
        display('PAC Preferred Phase');
        display([gcfg.subjectNames{i} ' preferred phase (ASIp): ' num2str(ASIp{i})]); 
        
        % plot circular stats & surrogates
        results.figurehandles.CFCtrialwise{i} = figure; set(gcf,'name',gcfg.subjectNames{i});        
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
    if ~isempty(ft_channelselection({'eeg'},{gcfg.osc_timelockchannel})) % if EEG channel
        Vtestvalue = gcfg.Vtestvalue; % CHANGE!!!!! was 0
    elseif strcmp(gcfg.osc_timelockchannel,'HC')
        Vtestvalue = gcfg.Vtestvalue; % CHANGE!!!!! was 180
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
    display(sprintf('%s%u%s%0.10f', 'V Test (GVp) against ', Vtestvalue, '° :', GVp));
    display(sprintf('%s%u%s%0.10f', 'V Test (GVp) against ', Vtestvalue, '° :', GV));
    
    % plot circular GROUP stats
    results.figurehandles.CFCtrialwiseGA = figure; 
    subplot(2,3,1); circ_plot([ASIp{:}]','pretty','o',true,'linewidth',2,'color','r');
    title({'distribution of individual preferred phase angles', ['mean angle: ' num2str(GASIp) ' rad, ' num2str(circ_rad2ang(GASIp)) 'degree']});
    subplot(2,3,2); hist(GSur_r,round(k/100));
    hold on; ylim = get(gca,'ylim'); line([Gr Gr], ylim,'color','r'); hold off;
    title({'surrogate distribution of group mean vector lengths',['with actual group r-value ' num2str(Gr) ' (red)']});   
          
    % summarize results for function output   
    subjects = gcfg.subjects;
    varList = {'subjects', 'trialN', 'AERPow', 'AERPowt', 'AERPowp', 'r','SIp','SIm', 'ASIm','Rz','Rp','z','z_p','ASIp','ASIp_ang',...
                'GAERPow', 'GAERPowt', 'GAERPowp', 'Gr', 'Gz', 'GSur_r_p', 'Gz_p', 'GASIp_ang', 'sASIp_ang', 'GRp', 'GRz', 'GVp', 'GV'};   
	for index = 1:numel(varList)
        results.(varList{index}) = eval(varList{index}); 
    end


% write subject-wise table
caseNames = gcfg.subjectNames';
varNames = {'N_events', 'ER_Pow_change', 'ER_Pow_t', 'ER_Pow_p', 'mean_vec_leng','magnitude_SI','Rayleigh_z','Rayleigh_p','mean_vec_leng_perm_z','mean_vec_leng_perm_z2p','preferred_phase_radians','preferred_phase_degree'};
Tsubjects = table(trialN', AERPow', AERPowt', AERPowp', r', ASIm', Rz', Rp', z', z_p', ASIp', ASIp_ang',...
        'RowNames', caseNames, 'VariableNames',varNames);

% write group table
GcaseNames  = {'GA_ER_Pow_change', 'GA_ER_Pow_t', 'GA_ER_Pow_p', 'GA_mean_vec_leng', 'GA_mean_vec_leng_perm_z', 'GA_mean_vec_leng_perm_p', 'GA_mean_vec_leng_perm_z2p', 'GA_preferred_phase_degree', 'ang_var_preferred_phase', 'GA_Rayleigh_p', 'GA_Rayleigh_z', ['V_Test_against_' num2str(gcfg.Vtestvalue) '_t'], ['V_Test_against_' num2str(gcfg.Vtestvalue) '_p']}';
GvarNames = {'group_results'};
Tgroup = table({GAERPow, GAERPowt, GAERPowp, Gr, Gz, GSur_r_p, Gz_p, GASIp_ang, sASIp_ang, GRp, GRz, GV, GVp}',...
        'RowNames',GcaseNames, 'VariableNames',GvarNames);

% plot subject-wise table
posTsubjects = get(subplot(2,3,4:6),'Position');
set(subplot(2,3,4:6),'YTick',[],'XTick',[]);
uitable('Data',Tsubjects{:,:},'ColumnName',Tsubjects.Properties.VariableNames,'RowName',Tsubjects.Properties.RowNames,'Units','Normalized','Position', posTsubjects);

% plot group table
posTgroup = get(subplot(2,3,3),'Position');
set(subplot(2,3,3),'YTick',[],'XTick',[]);
uitable('Data',Tgroup{:,:},'ColumnName',Tgroup.Properties.VariableNames,'RowName',Tgroup.Properties.RowNames, 'Units','Normalized','Position', posTgroup);
    
end % of phase amplitude coupling per trial


%%%%%%%%%%%%%%%%%%%%%%%%%
end % of function