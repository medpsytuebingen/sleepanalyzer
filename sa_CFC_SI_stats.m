%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sa_CFC_SI_stats
% by Til Ole Bergmann 2013
% last modified 2017/10/09 by TOB
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

function [results] = sa_CFC_SI_stats(gcfg)
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

%% phase amplitude coupling per trial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    



%%%%%%%%%%%%%%%%%%%%%%%%%
end % of function