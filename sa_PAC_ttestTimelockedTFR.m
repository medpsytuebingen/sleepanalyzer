%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sa_PAC_ttestTimelockedTFR
% by Til Ole Bergmann 2017
% last modified 2017/01/19 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% performs phase-amplitude coupling (PAC) analysis comparing timelocked TFRs with a t-test 
%
% output arguments:
% results = CFC results structure
% results.figurehandle =  handle to figure
%
% input arguments:
% GA_freq = Grand Average (with keeptrials = no) freq structure from fieldtrip time-frequency analysis (TFR) 
% GA_timelock = Grand Average (with keeptrials = no) timelock structure  from fieldtrip timelock analysis 
% All_freq  = 1xN cell array containing for N subject the normalized (e.g. % change from baseline) single-trial (keeptrials = yes) freq structures from fieldtrip time-frequencyanalysis (TFR) 
% All_raw_freq = 1xN cell array containing for N subject the raw (unnormalized) single-trial (keeptrials = yes) freq structures from fieldtrip time-frequencyanalysis (TFR) 
% All_timelock = 1xN cell array containing for N subject the single-trial (keeptrials = yes) timelock structures from fieldtrip timelocke analysis 

% cfg = [];
% cfg.subjectNames = 1xN cell array with the names of the subjects 
% cfg.subjects = 1xN cell array with the subject numbers 
%%%%%%% cfg.channel = Nx1 cell array with selection of channels (default = 'all'), see ft_channelselection for details
%%%%%%% cfg.layout = string indicating layout file, e.g., 'easycapM3.mat'
% cfg.individuallowfreq = 'yes' or 'no' whether individual frequency ranges are provided for the phase-providing frequencies
% cfg.individualhighfreq = 'yes' or 'no' whether individual frequency ranges are provided for the amplitude-providing frequencies
% cfg.hpfilter = frequency for high-pass filter for phase-providing frequency specified in Hz (no default)
% cfg.lpfilter = frequency for low-pass filter for phase-providing frequency specified in Hz (no default)
% cfg.timewin = time window around zero to be considered for phase-amplitude coupling, specified as [low high] in s (no default); 
% cfg.freqwin = frequency window to be considere for phase-amplitude coupling, specified as [low high] im Hz (no default);
% [results] = sa_PAC_ttestTimelockedTFR(cfg, GA_freq, All_freq, All_timelock);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results] = sa_PAC_ttestTimelockedTFR(gcfg, All_freq, All_timelock)

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
opt.tfrstats = 1; 
opt.tfrplot = 1;
opt.clusterstats = 1;
opt.clustextthresh = 1;
opt.ptd = 0; % peak trough differences

%% TFR non-parametric t-statistics against zero
if opt.tfrstats  
   
    % create zero data for comparison
	All_zero_freq = All_freq;
    for i = 1:length(All_freq) % create fake zero data
        All_zero_freq{i}.powspctrm(~isnan(All_zero_freq{i}.powspctrm)) = 0;
    end
    
    % non-parametric t-statistcs
    cfg = [];
    cfg.channel = 'all';
    cfg.latency = gcfg.plot_xlim; % [begin end] in seconds; 'all'
    cfg.trials = 'all';
    cfg.frequency = gcfg.plot_ylim; % [begin end]; 'all'
    cfg.parameter = 'powspctrm';   
    cfg.method = 'montecarlo';
    cfg.statistic = 'depsamplesT'; % dependent samples T-statistic	
    cfg.design(1,:) = [repmat(1:numel(All_freq),1,2)]; % subject vector        
    cfg.design(2,:) = [ones(1,length(All_freq)) ones(1,length(All_freq))*2]; % condition vector        
    cfg.uvar = 1; % number or list with indices, unit variable(s)
    cfg.ivar = 2; % number or list with indices, independent variable(s)   
    cfg.numrandomization = 'all'; % 'all'
    cfg.correctm = 'no'; % 'no' , 'max', cluster', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')
    cfg.alpha = 0.05; 
    cfg.tail = 0; %, -1, 1 or 0 (default = 0)   
    cfg.correcttail = 'prob'; % correct p-values or alpha-values when doing a two-sided test, 'alpha','prob' or 'no' (default = 'no')
    cfg.feedback = 'text'; % 'gui', 'text', 'textbar' or 'no' (default = 'textbar')    
    cfg.randomseed = 'yes'; % string, 'yes', 'no' or a number (default = 'yes')       
    [stat] = ft_freqstatistics(cfg, All_freq{:}, All_zero_freq{:}); 
%     [stat] = ft_freqstatistics(cfg, All_freq{:}, All_surrogate_freq{:}); 
%     [stat] = ft_freqstatistics(cfg, All_freq{:}, All_free_freq{:}); 

    results.stat = stat;

end


%% cluster statistics
if opt.clusterstats
    cfg.correctm = 'cluster'; % 'no' , 'max', cluster', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')
    cfg.clusterstatistic = 'maxsum'; % how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
    cfg.clusterthreshold = 'nonparametric_individual'; % method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
    cfg.clusteralpha = 0.05; % for either parametric or nonparametric thresholding per tail (default = 0.05)
    cfg.clustertail = 0; % -1, 1 or 0 (default = 0)        
    [cluststat] = ft_freqstatistics(cfg, All_freq{:}, All_zero_freq{:});     
%     [cluststat] = ft_freqstatistics(cfg, All_freq{:}, All_surrogate_freq{:}); 
%     [cluststat] = ft_freqstatistics(cfg, All_freq{:}, All_free_freq{:}); 

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
       
%% calculate GA
cfg = [];
GA_freq = ft_freqgrandaverage(cfg, All_freq{:});
GA_timelock = ft_timelockgrandaverage(cfg, All_timelock{:});
% GA_freq.time = -2.495:0.005:2.500; % !!!! RRELIMINARY FIX !!!
% results.stat.time = -0.995:0.005:1; % !!!! RRELIMINARY FIX !!!

%% plot stats
if opt.tfrplot    
       
    cfg = [];
    if ischar(gcfg.plot_xlim) && strcmp(gcfg.plot_xlim,'all')
            cfg.latency = [GA_freq.time(1) GA_freq.time(end)]; % [begin end] in seconds; 'all'   
    elseif isnumeric(gcfg.plot_xlim)
            cfg.latency = gcfg.plot_xlim;
    end
    
    if ischar(gcfg.plot_ylim) && strcmp(gcfg.plot_ylim,'all')
        cfg.frequency = [GA_freq.freq(1) GA_freq.freq(end)]; % [begin end] in seconds; 'all'
    elseif isnumeric(gcfg.plot_ylim)
        cfg.frequency = gcfg.plot_ylim; % [begin end]; 'all'
    end
    
    timelocktimeidx = GA_timelock.time >= cfg.latency(1) & GA_timelock.time <= cfg.latency(2);
    GAfreqtimeidx = GA_freq.time > cfg.latency(1) & GA_freq.time <= cfg.latency(2);
    statfreqtimeidx = results.stat.time > cfg.latency(1) & results.stat.time <= cfg.latency(2);
    freqidx = GA_freq.freq >= cfg.frequency(1) & GA_freq.freq <= cfg.frequency(2);   
    timelockRange = max(GA_timelock.avg)-min(GA_timelock.avg);
    freqRange = max(results.stat.freq)-min(results.stat.freq);
    
    plottimelockdata  = GA_timelock.avg(timelocktimeidx) * ((freqRange)/timelockRange*0.8);
    plottimelockdata  = plottimelockdata + abs(min(plottimelockdata)) + timelockRange*0.1;        
    plottimelocktime = GA_timelock.time(timelocktimeidx);
    
    results.figurehandles.tfrplot = figure;
    colormap fireice(64);
    spMargin = 0.04;
    
    subplot_tight(1,5,1,spMargin);
    plotdata = squeeze(GA_freq.powspctrm(1,freqidx,GAfreqtimeidx));
    imagesc(GA_freq.time(GAfreqtimeidx), GA_freq.freq(freqidx), plotdata, 'CDataMapping', 'scaled');
    set(gca,'YDir','normal');
    colorbar('location', 'EastOutside');
    set(gca,'clim', [-max(max(abs(plotdata))) max(max(abs(plotdata)))]);
    title('% change from baseline'); 
    hold on;    
    plot(plottimelocktime , plottimelockdata, 'LineWidth', 1, 'Color',[1 1 1]);
    hold off;    
    
    subplot_tight(1,5,2,spMargin);
    plotdata = squeeze(results.stat.stat(1,:,statfreqtimeidx));
    imagesc(results.stat.time(statfreqtimeidx), results.stat.freq, plotdata, 'CDataMapping', 'scaled');
    set(gca,'YDir','normal'); 
    colorbar('location', 'EastOutside');
    set(gca,'clim', [-max(max(abs(plotdata))) max(max(abs(plotdata)))]);
    title('t-values'); 
    hold on;    
    plot(plottimelocktime , plottimelockdata, 'LineWidth', 1, 'Color',[1 1 1]);
    hold off;    
    
    subplot_tight(1,5,3,spMargin);
    plotdata = squeeze(results.stat.stat(1,:,statfreqtimeidx) .* results.stat.mask(1,:,statfreqtimeidx)); 
    imagesc(results.stat.time(statfreqtimeidx), results.stat.freq, plotdata, 'CDataMapping', 'scaled');
    set(gca,'YDir','normal');
    colorbar('location', 'EastOutside');    
    set(gca,'clim', [min([-max(max(abs(plotdata))), -1e-10]) max([max(max(abs(plotdata))), 1e-10])]) % to make sure it is not [0 0]
    title({'thresholded t-values', '(p<0.05 uncorrected)'}); 
    hold on;    
    plot(plottimelocktime , plottimelockdata, 'LineWidth', 1, 'Color',[1 1 1]);
    hold off;    

	subplot_tight(1,5,4,spMargin);
    plotdata = squeeze(results.stat.stat(1,:,statfreqtimeidx) .* results.stat.cetmask(1,:,statfreqtimeidx));
    imagesc(results.stat.time(statfreqtimeidx), results.stat.freq, plotdata, 'CDataMapping', 'scaled');
    set(gca,'YDir','normal');
    colorbar('location', 'EastOutside');    
    set(gca,'clim', [min([-max(max(abs(plotdata))), -1e-10]) max([max(max(abs(plotdata))), 1e-10])]) % to make sure it is not [0 0]
    title({'thresholded t-values', ['(p<0.05 cluster extent threshold = ' num2str(k) ' corrected)']});
    hold on;    
    plot(plottimelocktime , plottimelockdata, 'LineWidth', 1, 'Color',[1 1 1]);
    hold off;
    
	subplot_tight(1,5,5,spMargin);
    plotdata = squeeze(results.stat.stat(1,:,statfreqtimeidx) .* results.cluststat.mask(1,:,statfreqtimeidx));
    imagesc(results.cluststat.time(statfreqtimeidx), results.cluststat.freq, plotdata, 'CDataMapping', 'scaled');
    set(gca,'YDir','normal');
    colorbar('location', 'EastOutside');    
    set(gca,'clim', [min([-max(max(abs(plotdata))), -1e-10]) max([max(max(abs(plotdata))), 1e-10])]) % to make sure it is not [0 0]
    title({'thresholded t-values', '(p<0.05 cluster-corrected)'});    
    hold on;    
    plot(plottimelocktime , plottimelockdata, 'LineWidth', 1, 'Color',[1 1 1]);
    hold off;    
    
end


%% analysis of peak-trough differences (ptd)
if opt.ptd 
    
    % settings    
    peaks = {}; troughs = {}; newpeaks = {}; newtroughs = {};
    fsample = round(1/mean(diff(All_timelock{1}.time))); % timelock time resolution
    tfrfsample = round(1/mean(diff(All_freq{1}.time))); % tfr time resolution
	dsfactor = fsample/tfrfsample; % downsample by factor
    
	for i = 1:numel(All_timelock) % loop over sujects
        
        switch gcfg.individuallowfreq
            case 'no'
                searchwin = find(All_timelock{1}.time >= gcfg.timewin(1) & All_timelock{1}.time <= gcfg.timewin(2));
            case 'yes'
                searchwin = find(All_timelock{1}.time >= gcfg.timewin(i,1) & All_timelock{1}.time <= gcfg.timewin(i,2));
        end
        newsearchwin = unique(int16((searchwin-1)/dsfactor));
        
        % high-pass filter spindel timelocked individual averages
        switch gcfg.individuallowfreq
            case 'no'
                All_timelock_hpf{i}.avg = ft_preproc_highpassfilter(All_timelock{i}.avg, double(fsample), gcfg.hpfilter , 3*fix(fsample/gcfg.hpfilter ), 'fir', 'twopass');
            case 'no'
                All_timelock_hpf{i}.avg = ft_preproc_highpassfilter(All_timelock{i}.avg, double(fsample), gcfg.hpfilter(i) , 3*fix(fsample/gcfg.hpfilter(i) ), 'fir', 'twopass');
        end
        % detect individual peaks and zero crossings automatically
        signal{i} = All_timelock_hpf{i}.avg;
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
        freqlim = All_freq{i}.freq >= gcfg.freqwin(1) & All_freq{i}.freq <= gcfg.freqwin(2); 
        [N p2nzcbin{i}] = histc(All_timelock{i}.time(newp2nzc{i}), All_freq{i}.time);        
        [N troughbin{i}] = histc(All_timelock{i}.time(newtroughs{i}), All_freq{i}.time); 
        [N n2pzcbin{i}] = histc(All_timelock{i}.time(newpeaks{i}), All_freq{i}.time);         
        [N peakbin{i}] = histc(All_timelock{i}.time(newn2pzc{i}), All_freq{i}.time); 
        avgp2nzcpow{i} = mean(mean(All_freq{i}.powspctrm(:,freqlim,unique(p2nzcbin{i})),3),2);                      
        avgtroughpow{i} = mean(mean(All_freq{i}.powspctrm(:,freqlim,unique(troughbin{i})),3),2);       
        avgn2pzcpow{i} = mean(mean(All_freq{i}.powspctrm(:,freqlim,unique(n2pzcbin{i})),3),2);
        avgpeakpow{i} = mean(mean(All_freq{i}.powspctrm(:,freqlim,unique(peakbin{i})),3),2);
    end
              
    tn = 0; % test number    
    tn = tn+1; % two-sided test of peak and trough values against each other
    [h{tn},p{tn},ci{tn},stats{tn}] = ttest([avgpeakpow{:}],[avgtroughpow{:}],0.05,'both'); % t-test
    p{tn}
    
    % plot
    dat = [avgp2nzcpow{:}; avgtroughpow{:}; avgn2pzcpow{:}; avgpeakpow{:}];
    results.figurehandles.PeakTroughDiff = figure;
    set(results.figurehandles.PeakTroughDiff,'DefaultAxesColorOrder',varycolor(13)); hold on;
    plot(dat, 'linewidth', 2);
    plot(mean(dat,2),'k','linewidth',3);
    legend(cellstr(num2str(gcfg.subjects(:)))','Location','EastOutside');
%     
end

%%%%%%%%%%%%%%%%%%%%%%%%%

results.cfg = gcfg; % add cfg input info to output

%%%%%%%%%%%%%%%%%%%%%%%%%
end % of function