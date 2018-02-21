%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spindle CFC statistics
% by Til Ole Bergmann 2013
% last modified 2016/11/22 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% performs statistics for spindle CFC 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results] = sa_CFC_stats_spindles(gcfg, GA_spindle_freq, GA_spindle_timelock, All_spindle_freq, All_spindle_raw_freq, All_spindle_timelock, spindle_tsv)

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
opt.tfrstats = 0; 
opt.tfrplot = 0;
opt.clusterstats = 0;
opt.clustextthresh = 0;
opt.plot = 0;
opt.saveplot = 0;

surrogate_type = 'fft_shuffle_phase'; % 'shuffle_phase', 'cyclic_shift_phase', 'random_shift_phase', 'fft_shuffle_phase' 

opt.ptd = 0; % peak trough differences
opt.paco = 0; % phase amplitude coupling OLD
opt.pacpt = 1; % phase amplitude coupling per trial
opt.pacps = 0; % phase amplitude coupling per subject

%% TFR non-parametric t-statistics against zero
if opt.tfrstats  
   
    % create zero data for comparison
	All_zero_freq = All_spindle_freq;
    for i = 1:length(All_spindle_freq) % create fake zero data
        All_zero_freq{i}.powspctrm(~isnan(All_zero_freq{i}.powspctrm)) = 0;
    end
    
    % non-parametric t-statistcs
    cfg = [];
    cfg.channel = 'all';
    cfg.latency = gcfg.plot_spindle_xlim; % [begin end] in seconds; 'all'
    cfg.trials = 'all';
    cfg.frequency = gcfg.plot_spindle_ylim; % [begin end]; 'all'
    cfg.parameter = 'powspctrm';   
    cfg.method = 'montecarlo';
    cfg.statistic = 'depsamplesT'; % dependent samples T-statistic	
    cfg.design(1,:) = [repmat(1:numel(All_spindle_freq),1,2)]; % subject vector        
    cfg.design(2,:) = [ones(1,length(All_spindle_freq)) ones(1,length(All_spindle_freq))*2]; % condition vector        
    cfg.uvar = 1; % number or list with indices, unit variable(s)
    cfg.ivar = 2; % number or list with indices, independent variable(s)   
    cfg.numrandomization = 'all'; % 'all'
    cfg.correctm = 'no'; % 'no' , 'max', cluster', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')
    cfg.alpha = 0.05; 
    cfg.tail = 0; %, -1, 1 or 0 (default = 0)   
    cfg.correcttail = 'prob'; % correct p-values or alpha-values when doing a two-sided test, 'alpha','prob' or 'no' (default = 'no')
    cfg.feedback = 'text'; % 'gui', 'text', 'textbar' or 'no' (default = 'textbar')    
    cfg.randomseed = 'yes'; % string, 'yes', 'no' or a number (default = 'yes')       
    [stat] = ft_freqstatistics(cfg, All_spindle_freq{:}, All_zero_freq{:}); 
%     [stat] = ft_freqstatistics(cfg, All_spindle_freq{:}, All_spindle_surrogate_freq{:}); 
%     [stat] = ft_freqstatistics(cfg, All_spindle_freq{:}, All_spindle_free_freq{:}); 

    results.stat = stat;

end


%% cluster statistics
if opt.clusterstats
    cfg.correctm = 'cluster'; % 'no' , 'max', cluster', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')
    cfg.clusterstatistic = 'maxsum'; % how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
    cfg.clusterthreshold = 'parametric'; % method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
    cfg.clusteralpha = 0.05; % for either parametric or nonparametric thresholding per tail (default = 0.05)
    cfg.clustertail = 0; % -1, 1 or 0 (default = 0)        
    [cluststat] = ft_freqstatistics(cfg, All_spindle_freq{:}, All_zero_freq{:});     
%     [cluststat] = ft_freqstatistics(cfg, All_spindle_freq{:}, All_spindle_surrogate_freq{:}); 
%     [cluststat] = ft_freqstatistics(cfg, All_spindle_freq{:}, All_spindle_free_freq{:}); 

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
    figure;
    subplot(4,1,1);
    plotdata = squeeze(GA_spindle_freq.powspctrm(1,GA_spindle_freq.freq >= cfg.frequency(1) & GA_spindle_freq.freq <= cfg.frequency(2), GA_spindle_freq.time >= cfg.latency(1) & GA_spindle_freq.time <= cfg.latency(2)));
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
    
%     subplot(4,1,4);
%     plotdata = squeeze(stat_cor.stat.*stat_cor.mask);
%     imagesc(stat_cor.time, stat_cor.freq, plotdata, 'CDataMapping', 'scaled');
%     set(gca,'YDir','normal');
%     colormap(jet); colorbar('location', 'EastOutside');    
% %     set(gca,'clim', [-max(max(abs(plotdata))) max(max(abs(plotdata)))]);
%     set(gca,'clim', [min([-max(max(abs(plotdata))), -1e-10]) max([max(max(abs(plotdata))), 1e-10])]) % to make sure it is not [0 0]
%     title('thresholded t-values (p<0.05 cluster-corrected)'); 

	% save plot
	if opt.saveplot
        set(gcf, 'Name',[gcfg.conditionLabels{gcfg.thisCondition} '_GA_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' spindle_timelockchannel '_stats_' gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
        filename = fullfile(gcfg.paths.results, 'group', [gcfg.conditionLabels{gcfg.thisCondition} '_GA_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' spindle_timelockchannel '_stats_' gcfg.timestamp]);
        saveas(gcf,[filename '.fig']);        
        set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [25.4 19.05], 'PaperPosition', [-3 -1 30 19.5]);
        print(gcf, '-dpdf', '-r300', [filename '.pdf']);      
    end
end



%% analysis of peak-trough differences (ptd)
if opt.ptd 
    
    % settings    
    peaks = {}; troughs = {}; newpeaks = {}; newtroughs = {};
    fsample = round(1/mean(diff(All_spindle_timelock{1}.time))); % timelock time resolution
    tfrfsample = round(1/mean(diff(All_spindle_freq{1}.time))); % tfr time resolution
	dsfactor = fsample/tfrfsample; % downsample by factor
    searchwin = find(All_spindle_timelock{1}.time >= gcfg.timewin(1) & All_spindle_timelock{1}.time <= gcfg.timewin(2));
	newsearchwin = unique(int16((searchwin-1)/dsfactor));                
    
	for i = 1:numel(All_spindle_timelock) % loop over sujects
        
        % high-pass filter spindel timelocked individual averages
        All_spindle_timelock_hpf{i}.avg = ft_preproc_highpassfilter(All_spindle_timelock{i}.avg, double(fsample), gcfg.hpfilter , 3*fix(fsample/gcfg.hpfilter )+1, 'fir', 'twopass'); % Mathilde's formula
    
        % detect individual peaks and zero crossings automatically
        signal{i} = All_spindle_timelock_hpf{i}.avg;
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
        freqlim = All_spindle_freq{i}.freq >= gcfg.freqwin(1) & All_spindle_freq{i}.freq <= gcfg.freqwin(2); 
        [N p2nzcbin{i}] = histc(All_spindle_timelock{i}.time(newp2nzc{i}), All_spindle_freq{i}.time);        
        [N troughbin{i}] = histc(All_spindle_timelock{i}.time(newtroughs{i}), All_spindle_freq{i}.time); 
        [N n2pzcbin{i}] = histc(All_spindle_timelock{i}.time(newpeaks{i}), All_spindle_freq{i}.time);         
        [N peakbin{i}] = histc(All_spindle_timelock{i}.time(newn2pzc{i}), All_spindle_freq{i}.time); 
        avgp2nzcpow{i} = mean(mean(All_spindle_freq{i}.powspctrm(:,freqlim,unique(p2nzcbin{i})),3),2);                      
        avgtroughpow{i} = mean(mean(All_spindle_freq{i}.powspctrm(:,freqlim,unique(troughbin{i})),3),2);       
        avgn2pzcpow{i} = mean(mean(All_spindle_freq{i}.powspctrm(:,freqlim,unique(n2pzcbin{i})),3),2);
        avgpeakpow{i} = mean(mean(All_spindle_freq{i}.powspctrm(:,freqlim,unique(peakbin{i})),3),2);
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


%% circular-circular correlation of lfo timelocked signal and hfo power modulation and

%% % phase amplitude coupling old %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opt.paco  % phase amplitude coupling old

    % settings
    fsample = round(1/mean(diff(All_spindle_timelock{1}.time))); % timelock time resolution
    tfrfsample = round(1/mean(diff(All_spindle_freq{1}.time))); % tfr time resolution
	dsfactor = 1; % fsample/tfrfsample; % downsample modulating lf signal by factor  !!! OR !!!
    usfactor = fsample/tfrfsample; % upsample modulated hf power by factor  !!! OR !!!
    searchwin = find(All_spindle_timelock{1}.time >= gcfg.timewin(1) & All_spindle_timelock{1}.time <= gcfg.timewin(2));
	newsearchwin = unique(int16((searchwin-1)/dsfactor));   
    ima = sqrt(-1); % imaginary number
	eul = exp(1); % Euler number       
        
    for i = 1:numel(All_spindle_timelock) % loop over sujects
    
        freqlim = All_spindle_freq{i}.freq >= gcfg.freqwin(1) & All_spindle_freq{i}.freq <= gcfg.freqwin(2);     
        
        % high-pass filter spindel timelocked individual averages
        All_spindle_timelock_hpf{i}.avg = ft_preproc_highpassfilter(All_spindle_timelock{i}.avg, double(fsample), gcfg.hpfilter , 3*fix(fsample/gcfg.hpfilter )+1, 'fir', 'twopass'); % Mathilde's formula
    
        % CFC        
        filt{i} = All_spindle_timelock_hpf{i}.avg; 
        filt_ds{i} = downsample(filt{i},dsfactor);
        hilbert_phase{i} = angle(hilbert(filt{i}));
        hilbert_phase_ds{i} = downsample(hilbert_phase{i},dsfactor); 
        
        hfopow{i} = squeeze(mean(All_spindle_freq{i}.powspctrm(:,freqlim,:),2))';
        hfopow{i}(isnan(hfopow{i})) = 0;
        hfopow{i} = interp(hfopow{i},usfactor);
        hfopow_filt{i} = ft_preproc_highpassfilter(hfopow{i}, double(fsample/dsfactor), gcfg.hpfilter , 3*fix(fsample/dsfactor/gcfg.hpfilter )+1, 'fir', 'twopass'); % Mathilde's formula        
        hilbert_phase_hfopow_filt{i} = angle(hilbert(hfopow_filt{i}));    
        
        hforawpow{i} = squeeze(mean(All_spindle_raw_freq{i}.powspctrm(:,freqlim,:),2))';
        hforawpow{i}(isnan(hforawpow{i})) = 0;    
        hforawpow{i} = interp(hforawpow{i},usfactor);
        hforawpow_filt{i} = ft_preproc_highpassfilter(hforawpow{i}, double(fsample/dsfactor), gcfg.hpfilter , 3*fix(fsample/dsfactor/gcfg.hpfilter )+1, 'fir', 'twopass'); % Mathilde's formula
                
        % circular correlations
        [cfc_corrcc_hfopow_r{i} cfc_corrcc_hfopow_p{i}] = circ_corrcc(hilbert_phase_ds{i}(:,newsearchwin)',hilbert_phase_hfopow_filt{i}(:,newsearchwin)');           
        cfc_corrccabs_hfopow_r{i} = abs(cfc_corrcc_hfopow_r{i});
        [cfc_corrcl_hfopow_r{i} cfc_corrcl_hfopow_p{i}] = circ_corrcl(hilbert_phase_ds{i}(:,newsearchwin)',hfopow_filt{i}(:,newsearchwin)');   
        [cfc_corrcl_hforawpow_r{i} cfc_corrcl_hforawpow_p{i}] = circ_corrcl(hilbert_phase_ds{i}(:,newsearchwin)',hforawpow_filt{i}(:,newsearchwin)');                      

        % synchronization index (cf. Cohen et al., J Neurosci Methods 2008)
        phi1{i} = hilbert_phase_ds{i}(:,newsearchwin)'; % phase values lf modulating signal
        phi2{i} = hilbert_phase_hfopow_filt{i}(:,newsearchwin)'; % phase values hf power modulation 
        SI{i} = mean(exp(1i*[phi1{i} - phi2{i}]')); % Synchronization Index (cf. Cohen et al., J Neurosci Methods 2008)
        SIm{i} = abs(SI{i}); % SI magnitude (0 = perfect desychronization to 1 = perfect synchronization)
        SIp{i} = (atan2(imag(SI{i}),real(SI{i}))); % preferred phase of the synchronization

        % prepare data for 'fft_shuffle_phase' option of surrogate generation
        if strcmp(surrogate_type,'fft_shuffle_phase')
            % preprocess data
            data{i} = hfopow_filt{i}(:,newsearchwin);
            data{i} = detrend(data{i},'linear'); % demean and detrend

            % FFT
            fourier{i} = fft(data{i});
            fabs{i} = abs(fourier{i});
            fphase{i} = angle(fourier{i});
            
            L = length(data{i});
        end       
        
%         figure;      
%         subplot(2,1,1);hold on;
%             plot(filt{i},'k');
%             plot(filt_ds{i},'g'); 
%             plot(hfopow{i},'r');
%             plot(hfopow_filt{i},'b');
%         subplot(2,1,2);hold on;
%             plot(hilbert_phase_ds{i},'g'); 
%             plot(hilbert_phase_hfopow_filt{i},'b');     
    end % of loop over subjects


    % build surrogate distributions of r-values for subjects and t-values for group statistics
    clear cfc_corrcc_hfopow_surr_r cfc_corrccabs_hfopow_surr_r cfc_corrcl_hfopow_surr_r SI_surr SIm_surr SIp_surr;
        for k = 1:10000 % numel(newsearchwin) in case of cyclic shuffle                 
            % surogate distributions of r-values for subjects 
            switch surrogate_type
                case 'shuffle_phase' % random SHUFFLE of phase values           
                    permindex = newsearchwin(randperm(length(newsearchwin))); 
                case 'cyclic_shift_phase' % cyclically SHIFT phase values           
                    permindex = circshift(newsearchwin,[0 k]); 
                case 'random_shift_phase' % randomly SHIFT phase values           
                    permindex = newsearchwin + round(rand(1)*2*numel(newsearchwin)-numel(newsearchwin)); 
                case 'fft_shuffle_phase' % according to J�rgen Fell's code snippet
                    randphase = rand(1,L).*(2*pi);                        
            end 

            for i = 1:numel(gcfg.subjects) % loop over sujects
                % for circular-circular correlation
                hilbert_phase_hfopow_filt_surr{i} = hilbert_phase_hfopow_filt{i};
                switch surrogate_type 
                    case {'shuffle_phase', 'cyclic_shift_phase', 'random_shift_phase'}
                        hilbert_phase_hfopow_filt_surr{i}(:,newsearchwin) = hilbert_phase_hfopow_filt{i}(:,permindex);
                    case 'fft_shuffle_phase'  
                        fphasenew{i} = fphase{i} + randphase;
                        fouriernew{i} = fabs{i}.*(cos(fphasenew{i})+1i*sin(fphasenew{i}));
                        fouriernew{i}([1 end]) = fourier{i}([1 end]);
                        for sample = 2:L/2 % recreate symmetry of fourier transformed data
                            fouriernew{i}(sample) = real(fouriernew{i}(L)-sample+2) - imag(fouriernew{i}(L)-sample+2)*1i;
                        end
                        data_surr{i} = real(ifft(fouriernew{i})); % iFFT
                        data_surr_phase{i} = angle(hilbert(data_surr{i})); % phase
                        hilbert_phase_hfopow_filt_surr{i}(:,newsearchwin) = data_surr_phase{i};
%                         if i == 1 % check data for subject 2
%                             hold on; plot(hilbert_phase_hfopow_filt_surr{i}(:,newsearchwin));
%                         end
                end
                [cfc_corrcc_hfopow_surr_r{i}(k) cfc_corrcc_hfopow_surr_p{i}(k)] = circ_corrcc(hilbert_phase_ds{i}(:,newsearchwin)',hilbert_phase_hfopow_filt_surr{i}(:,newsearchwin)');   
                cfc_corrccabs_hfopow_surr_r{i}(k) = abs(cfc_corrcc_hfopow_surr_r{i}(k));

                % for circular-linear correlation
                hfopow_filt_surr{i} = hfopow_filt{i};
                switch surrogate_type 
                    case {'shuffle_phase', 'cyclic_shift_phase', 'random_shift_phase'}                
                        hfopow_filt_surr{i}(:,newsearchwin) = hfopow_filt{i}(:,permindex);
                    case 'fft_shuffle_phase'  
                        hfopow_filt_surr{i}(:,newsearchwin) = data_surr{i}; 
                end
                [cfc_corrcl_hfopow_surr_r{i}(k) cfc_corrcl_hfopow_surr_p{i}(k)] = circ_corrcl(hilbert_phase_ds{i}(:,newsearchwin)',hfopow_filt_surr{i}(:,newsearchwin)');               

                % for synchronization index (no switch necessary as hilbert_phase_hfopow_filt_surr{i}(:,newsearchwin) has been adjusted above)
                phi2_surr{i} = hilbert_phase_hfopow_filt_surr{i}(:,newsearchwin); % phase values hf power modulation 
                SI_surr{i}(k) = mean(exp(1i*[phi1{i}' - phi2_surr{i}]'))'; % Synchronization Index (cf. Cohen et al., J Neurosci Methods 2008)
                SIm_surr{i}(k) = abs(SI_surr{i}(k))'; % SI magnitude (0 = perfect desychronization to 1 = perfect synchronization)

            end % of loop over subjects

            
            % surrogate distributions of t-values for group statistics
            groupdata_cc = shiftdim(reshape([cfc_corrcc_hfopow_surr_r{:}]',numel(gcfg.subjects),k),1);
            [h_cc,p_cc,ci_cc,stats_cc] = ttest(groupdata_cc(k,:),zeros(1,length(groupdata_cc(k,:))),0.05,'both'); % t-test
            permdistr_cc_t(k) = stats_cc.tstat;        

            groupdata_ccabs = shiftdim(reshape([cfc_corrccabs_hfopow_surr_r{:}]',numel(gcfg.subjects),k),1);
            [h_ccabs,p_ccabs,ci_ccabs,stats_ccabs] = ttest(groupdata_cc(k,:),zeros(1,length(groupdata_cc(k,:))),0.05,'both'); % t-test
            permdistr_ccabs_t(k) = stats_ccabs.tstat;  

            groupdata_cl = shiftdim(reshape([cfc_corrcl_hfopow_surr_r{:}]',numel(gcfg.subjects),k),1);
            [h_cl,p_cl,ci_cl,stats_cl] = ttest(groupdata_cl(k,:),zeros(1,length(groupdata_cl(k,:))),0.05,'both'); % t-test
            permdistr_cl_t(k) = stats_cl.tstat;

            groupdata_SIm = shiftdim(reshape([SIm_surr{:}]',numel(gcfg.subjects),k),1);
            [h_cl,p_cl,ci_cl,stats_SIm] = ttest(groupdata_SIm(k,:),zeros(1,length(groupdata_SIm(k,:))),0.05,'both'); % t-test
            permdistr_SIm_t(k) = stats_SIm.tstat;        

    %         figure;
    %         subplot(1,4,1); hist(groupdata_cc(:,1),100);
    %         subplot(1,4,2); hist(groupdata_ccabs(:,1),100);
    %         subplot(1,4,3); hist(groupdata_cl(:,1),100);
    %         subplot(1,4,4); hist(groupdata_SIm(:,1),100)
        
        fprintf('Iteration %d \r ', k);
        end % of loop over k iterations
            
            
%         case 'fft_shuffle_phase' 
% %         for i = 1:numel(gcfg.subjects)
% %             subjectNames{i} = sprintf('%s%0.2d','S',gcfg.subjects(i));
% %         end
% %         subjectNames = sprintf('%s%0.2d','S',gcfg.subjects(i));
% %         load(fullfile(gcfg.paths.results,subjectNames{i},[subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_targetchannel{gcfg.subjects(i)} '_surrogate_freq.mat']), 'spindle_surrogate_freq'); % load surrogate data
%         fouriercheck = fft(data_surr{i});
%        
%         % plot stuff
%         figure;
%         set(gcf,'WindowStyle','docked')
% 
%         subplot(3,1,1);
%         title('data (blue = original, red = after phase-scrambling)');
%         cla;hold on;
%         plot(data,'b'); % original data
%         plot(datanew+1,'r'); % backtransformed data
%         hold off;
% 
%         subplot(3,1,2);
%         title('fourier complex values  (blue = original, red = after phases-crambling)');
%         cla;hold on;
%         plot(fourier,'b'); % original data
%         plot(fouriercheck+1,'r'); % backtransformed data
%         hold off;
% 
%         subplot(3,1,3);
%         title('power spectrum  (blue = original, red = after phase-scrambling)');
%         cla;hold on;
%         gcfg.freqwin = [0 50]; % frequency window to plot
%         freq = Fs/2 * linspace(0,1,L/2+1);
%         powspctrm_orig = 2*abs(fourier(1:L/2+1)); % original data
%         plot(freq(freq>=gcfg.freqwin(1) & freq<=gcfg.freqwin(2)), powspctrm_orig(freq>=gcfg.freqwin(1) & freq<=gcfg.freqwin(2)),'b'); % original data
%         powspctrm_new = 2*abs(fouriercheck(1:L/2+1)); % backtransformed data
%         plot(freq(freq>=gcfg.freqwin(1) & freq<=gcfg.freqwin(2)), powspctrm_new(freq>=gcfg.freqwin(1) & freq<=gcfg.freqwin(2)),'r'); % original data
%         hold off;
        
    
    % subject-wise permutation test against surrogates 
    for i = 1:numel(All_spindle_timelock) % loop over subjects        
        ixcc = find(sort(cfc_corrcc_hfopow_surr_r{i}) > cfc_corrcc_hfopow_r{i}, 1);
        ixccabs = find(sort(cfc_corrccabs_hfopow_surr_r{i}) > cfc_corrccabs_hfopow_r{i}, 1);
        ixcl = find(sort(cfc_corrcl_hfopow_surr_r{i}) > cfc_corrcl_hfopow_r{i}, 1);        
        ixSIm = find(sort(SIm_surr{i}) > SIm{i}, 1);
        if isempty(ixcc)
            ixcc = length(cfc_corrcc_hfopow_surr_r{i});
        end   
        if isempty(ixccabs)
            ixccabs = length(cfc_corrccabs_hfopow_surr_r{i});
        end          
        if isempty(ixcl)
            ixcl = length(cfc_corrcl_hfopow_surr_r{i});
        end
        if isempty(ixSIm)
            ixSIm = length(SIm_surr{i});
        end        
        cfc_corrcc_hfopow_p{i} = 1-(ixcc/length(cfc_corrcc_hfopow_surr_r{i}));      
        cfc_corrccabs_hfopow_p{i} = 1-(ixccabs/length(cfc_corrccabs_hfopow_surr_r{i}));      
        cfc_corrcl_hfopow_p{i} = 1-(ixcl/length(cfc_corrcl_hfopow_surr_r{i}));              
        SIm_p{i} = 1-(ixSIm/length(SIm_surr{i}));         
    end
% 	  [gcfg.subjects; cfc_corrcc_hfopow_r{:}; cfc_corrcc_hfopow_p{:}]
%     [gcfg.subjects; cfc_corrccabs_hfopow_r{:}; cfc_corrccabs_hfopow_p{:}]
%     [gcfg.subjects; cfc_corrcl_hfopow_r{:}; cfc_corrcl_hfopow_p{:}]
    [gcfg.subjects; SIm{:}; SIm_p{:}; SIp{:}]
    
    
	% group permutation test against surrogates 
    [h_cc,p_cc,ci_cc,stats_cc] = ttest([cfc_corrcc_hfopow_r{:}],zeros(1,length(cfc_corrcc_hfopow_r)),0.05,'both'); % t-test
    [h_ccabs,p_ccabs,ci_ccabs,stats_ccabs] = ttest(abs([cfc_corrcc_hfopow_r{:}]),zeros(1,length(cfc_corrcc_hfopow_r)),0.05,'right'); % t-test
    [h_cl,p_cl,ci_cl,stats_cl] = ttest([cfc_corrcl_hfopow_r{:}],zeros(1,length(cfc_corrcl_hfopow_r)),0.05,'right'); % t-test
    [h_cl,p_cl,ci_cl,stats_SIm] = ttest([SIm{:}],zeros(1,length(SIm)),0.05,'right'); % t-test
    ixcc = find(sort(permdistr_cc_t) > stats_cc.tstat, 1); 
	ixccabs = find(sort(permdistr_cc_t) > stats_ccabs.tstat, 1);
    ixcl = find(sort(permdistr_cl_t) > stats_cl.tstat, 1);    
    ixSIm = find(sort(permdistr_SIm_t) > stats_SIm.tstat, 1);    
    if isempty(ixcc)
        ixcc = length(permdistr_cc_t);
    end       
    if isempty(ixccabs)
        ixccabs = length(permdistr_cc_t);
    end
    if isempty(ixcl)
        ixcl = length(permdistr_cl_t);
    end
    if isempty(ixSIm)
        ixSIm = length(permdistr_SIm_t);
    end    
    group_corrcc_p = 1-(ixcc/length(permdistr_cc_t)); 
    group_corrccabs_p = 1-(ixccabs/length(permdistr_cc_t)); 
    group_corrcl_p = 1-(ixcl/length(permdistr_cl_t));
    group_SIm_p = 1-(ixSIm/length(permdistr_SIm_t));
    group_SIm = mean([SIm{:}]);
    group_SIp = mean([SIp{:}]);
    
%     display(['group_corrcc_p: ' num2str(group_corrcc_p)]); 
%     display(['group_corrccabs_p: ' num2str(group_corrccabs_p)]);           
%     display(['group_corrcl_p: ' num2str(group_corrcl_p)]);        

    display(['group_SIm: ' num2str(group_SIm)]); 
    display(['group_SIm_p: ' num2str(group_SIm_p)]); 
    display(['group_SIp: ' num2str(group_SIp)]); 

    
    % plot circular stats
    circ_plot([SIp{:}]','pretty','o',true,'linewidth',2,'color','r');
    
    % Rayleigh test on non-uniform distribution of preferred phase angles over subjects
    [p,z] = circ_rtest([SIp{:}]');
    sprintf('%s%0.10f', 'Rayleigh test p-value: ', p)
    sprintf('%s%0.10f', 'Rayleigh''s z-value: ', z)

end


%% phase amplitude coupling per trial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opt.pacpt % phase amplitude coupling per trial
   
    subjectNames = gcfg.subjectNames;

    for i = 1:numel(gcfg.subjects) % loop over sujects
       
        % load trialwise data of subject i
        load(fullfile(gcfg.paths.results,subjectNames{i},[subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_targetchannel{gcfg.subjects(i)} '_eventdata.mat']), 'spindle_eventdata'); % load freq data
        load(fullfile(gcfg.paths.results,subjectNames{i},[subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_targetchannel{gcfg.subjects(i)} '_freq.mat']), 'spindle_freq'); % load freq data
        
        % select trials        
        cfg = [];
        cfg.trials = find(spindle_tsv{i});
        cfg.tolerance   = 1e-3;
        spindle_eventdata = ft_selectdata(cfg,spindle_eventdata); 
        spindle_freq = ft_selectdata(cfg,spindle_freq);
        
        % define data 
        freqindex = spindle_freq.freq >= gcfg.freqwin(1) & spindle_freq.freq <= gcfg.freqwin(2);     
        timeindex = unique(int16(find(spindle_eventdata.time{1} >= gcfg.timewin(1) & spindle_eventdata.time{1} <= gcfg.timewin(2))));
    
        trialN{i} = size(spindle_eventdata.trial,2); % # trials        
        for j = 1:trialN{i} % loop over trials
            
            % define data for trial j
            LF = spindle_eventdata.trial{j};
            HF = squeeze(mean(spindle_freq.powspctrm(j,1,freqindex,:),3))';
            
            % preprocess LF data
            LFhp = ft_preproc_bandpassfilter(LF, double(spindle_eventdata.fsample), [gcfg.hpfilter  gcfg.lpfilter ], 3*fix(spindle_eventdata.fsample/gcfg.hpfilter )+1, 'fir', 'twopass'); % Mathilde's formula
            
            % preprocess HF data
%             usfactor = double(spindle_eventdata.fsample / round(1/mean(diff(spindle_freq.time)))); %  upsampling factor to interpolate HF (TFR) signal to LF signal
            HF(isnan(HF)) = 0; % replace NaN by 0
%             HFus = interp(HF,usfactor); % upsample to fit LF time resolution
            HFus = transpose(resample(transpose(HF),spindle_eventdata.fsample,round(1/mean(diff(SO_freq.time))))); % upsample to fit LF time resolution
            HFhp = ft_preproc_highpassfilter(HFus, double(spindle_eventdata.fsample), gcfg.hpfilter , 3*fix(spindle_eventdata.fsample/gcfg.hpfilter )+1, 'fir', 'twopass'); % Mathilde's formula        
            
            % LF event-locked HF power change (difference)
            bl = [-2.0 -1.0]; % for SO it would be bl = [-2.5 -1.5]
            blindex = unique(int16(find(spindle_eventdata.time{1} >= bl(1) & spindle_eventdata.time{1} <= bl(2))-1));
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
 
        % single-subject chi� test on non-uniform distribution
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
        set(gcf, 'Name',[subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_CFCtrialwise_' gcfg.sprintf('[%2.2f %2.2fs]', gcfg.timewin(1),gcfg.timewin(2)) '_' sprintf('[%2.2f %2.2fHz]',gcfg.freqwin(1),gcfg.freqwin(2)) '_' gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
        filename = fullfile(gcfg.paths.results, subjectNames{i}, [subjectNames{i} '_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_CFCtrialwise_' gcfg.timewin(1),gcfg.timewin(2)) '_' sprintf('[%2.2f %2.2fHz]',gcfg.freqwin(1),gcfg.freqwin(2)) '_' gcfg.timestamp]);
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
    if ~isempty(ft_channelselection({'eeg'},{gcfg.spindle_timelockchannel})) % if EEG channel
        Vtestvalue = 0;
    elseif strcmp(gcfg.spindle_timelockchannel,'HC')
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
    display(sprintf('%s%u%s%0.10f', 'V Test (GVp) against ', Vtestvalue, '� :', GVp));
    display(sprintf('%s%u%s%0.10f', 'V Test (GVp) against ', Vtestvalue, '� :', GV));
    
    % plot circular GROUP stats
    figure; 
    subplot(1,2,1); circ_plot([ASIp{:}]','pretty','o',true,'linewidth',2,'color','r');
    title({'distribution of individual preferred phase angles', ['mean angle: ' num2str(GASIp) ' rad, ' num2str(circ_rad2ang(GASIp)) 'degree']});
    subplot(1,2,2); hist(GSur_r,round(k/100));
    hold on; ylim = get(gca,'ylim'); line([Gr Gr], ylim,'color','r'); hold off;
    title({'surrogate distribution of group mean vector lengths',['with actual group r-value ' num2str(Gr) ' (red)']});   
    
    % save figure
    set(gcf, 'Name',[gcfg.conditionLabels{gcfg.thisCondition} '_GA_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_CFCtrialwise_' gcfg.timewin(1),gcfg.timewin(2)) '_' sprintf('[%2.2f %2.2fHz]',gcfg.freqwin(1),gcfg.freqwin(2)) '_' gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
    filename = fullfile(gcfg.paths.results, 'group', [gcfg.conditionLabels{gcfg.thisCondition} '_GA_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_CFCtrialwise_' gcfg.timewin(1),gcfg.timewin(2)) '_' sprintf('[%2.2f %2.2fHz]',gcfg.freqwin(1),gcfg.freqwin(2)) '_' gcfg.timestamp]);
    saveas(gcf,[filename '.fig']);   
    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [25.4 19.05], 'PaperPosition', [-3 -1 30 19.5]);
    print(gcf, '-dpdf', '-r300', [filename '.pdf']); 

        
    % summarize results for function output   
    subjects = gcfg.subjects;
    varList = {'subjects', 'trialN', 'AERPow', 'AERPowt', 'AERPowp', 'r','SIp','SIm', 'ASIm','Rz','Rp','z','z_p','ASIp','ASIp_ang',...
                'GAERPow', 'GAERPowt', 'GAERPowp', 'Gr', 'Gz', 'GSur_r_p', 'Gz_p', 'GASIp_ang', 'sASIp_ang', 'GRp', 'GRz', 'GVp', 'GV'};   
	for index = 1:numel(varList)
        results.(varList{index}) = eval(varList{index}); 
    end

    % save results !!! --> MOVED TO *_sl !!!
%     filename = fullfile(gcfg.paths.results, 'group', [gcfg.conditionLabels{gcfg.thisCondition} '_GA_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_CFCtrialwise_' gcfg.timewin(1),gcfg.timewin(2)) '_' sprintf('[%2.2f %2.2fHz]',gcfg.freqwin(1),gcfg.freqwin(2)) '_' gcfg.timestamp '.mat']);
%     save('-v7.3', filename, varList{:});  
    
    
end % of phase amplitude coupling per trial





%% phase amplitude coupling per subject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opt.pacps 

    % settings
    fsample = round(1/mean(diff(All_spindle_timelock{1}.time))); % timelock time resolution
    
    % peparation   
    for i = 1:numel(gcfg.subjects)
        subjectNames{i} = sprintf('%s%0.2d','S',gcfg.subjects(i));
    end

    for i = 1:numel(gcfg.subjects) % loop over sujects
       
        % define data 
        freqindex = All_spindle_freq{i}.freq >= gcfg.freqwin(1) & All_spindle_freq{i}.freq <= gcfg.freqwin(2);     
        timeindex = unique(int16(find(All_spindle_timelock{1}.time >= gcfg.timewin(1) & All_spindle_timelock{1}.time <= gcfg.timewin(2))-1));   
            
        % define data for subject i
        LF = All_spindle_timelock{i}.avg;
        HF = squeeze(mean(All_spindle_freq{i}.powspctrm(:,freqindex,:),2))';

        % preprocess LF data
        LFhp = ft_preproc_highpassfilter(LF, double(fsample), gcfg.hpfilter , 3*fix(fsample/gcfg.hpfilter )+1, 'fir', 'twopass'); % Mathilde's formula

        % preprocess HF data
        usfactor = fsample / round(1/mean(diff(All_spindle_freq{i}.time))); %  upsampling factor to interpolate HF (TFR) signal to LF signal
        HF(isnan(HF)) = 0; % replace NaN by 0
        HFus = interp(HF,usfactor); % upsample to fit LF time resolution
        HFhp = ft_preproc_highpassfilter(HFus, double(fsample), gcfg.hpfilter , 3*fix(fsample/gcfg.hpfilter )+1, 'fir', 'twopass'); % Mathilde's formula        

        % get phase angles (hilbert transform)
        LFph = angle(hilbert(LFhp)); 
        HFph = angle(hilbert(HFhp)); 

        % synchronization index (cf. Cohen et al., J Neurosci Methods 2008)
        Aphi1{i} = LFph(:,timeindex); % phase values lf modulating signal
        Aphi2{i} = HFph(:,timeindex); % phase values hf power modulation 
        ASI{i} = mean(exp(1i*[Aphi1{i} - Aphi2{i}])); % Synchronization Index (SI)
        ASIm{i} = abs(ASI{i}); % SI magnitude (0 = perfect desychronization to 1 = perfect synchronization)
        ASIp{i} = (atan2(imag(ASI{i}),real(ASI{i}))); % preferred phase of the synchronization
                                        
        % SI descriptives
        display([subjectNames{i} ' ASIm: ' num2str(mean(ASIm{i}))]); 
        display([subjectNames{i} ' ASIp: ' num2str(mean(ASIp{i}))]); 
               
    end % of loop over sujects
        
    % Rayleigh test on non-uniform distribution of preferred phase angles over subjects
    [GRp,GRz] = circ_rtest([ASIp{:}]');
    sprintf('%s%0.10f', 'Rayleigh test p-value: ', GRp)
    sprintf('%s%0.10f', 'Rayleigh''s z-value: ', GRz)
    
    % results
    [gcfg.subjects; ASIp{:}; circ_rad2ang([ASIp{:}])]  

    % plot circular stats
    figure;
    circ_plot([ASIp{:}]','pretty','o',true,'linewidth',2,'color','r');
    
    % save figure
    set(gcf, 'Name',[gcfg.conditionLabels{gcfg.thisCondition} '_GA_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_CFCavgwise_' gcfg.timestamp],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
    filename = fullfile(gcfg.paths.results, 'group', [gcfg.conditionLabels{gcfg.thisCondition} '_GA_spindle_' gcfg.spindle_timelockevent '_' gcfg.spindle_tsv_flag '_' gcfg.spindle_timelockchannel '_CFCavgwise_' gcfg.timestamp]);
    saveas(gcf,[filename '.fig']);   
    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [25.4 19.05], 'PaperPosition', [-3 -1 30 19.5]);
    print(gcf, '-dpdf', '-r300', [filename '.pdf']); 
    
end % of phase amplitude coupling per subject




end % of function