%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sa_PAC_SI
% by Til Ole Bergmann 2017
% last modified 2017/12/16 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% performs phase-amplitude coupling (PAC) analysis using the Synchronization Index (cf. Cohen et al., J Neurosci Methods 2008)
%
% output arguments:
% results = PAC_SI results structure
% results.figurehandle =  handle to figure

% cfg = [];
% cfg.subjectNames = 1xN cell array with the names of the subjects 
% cfg.subjects = 1xN cell array with the subject numbers 
% cfg.individuallowfreq = 'yes' or 'no' whether individual frequency ranges are provided for the phase-providing frequencies
% cfg.individualhighfreq = 'yes' or 'no' whether individual frequency ranges are provided for the amplitude-providing frequencies
% cfg.hpfilter = frequency for high-pass filter for phase-providing frequency specified in Hz (no default), in cell array with one cell per subject if cfg.individuallowfreq = 'yes'
% cfg.lpfilter = frequency for low-pass filter for phase-providing frequency specified in Hz (no default), in cell array with one cell per subject if cfg.individuallowfreq = 'yes'
% cfg.timewin = time window around zero to be considered for phase-amplitude coupling, specified as [low high] in s (no default), in cell array with one cell per subject if cfg.individuallowfreq = 'yes' 
% cfg.individualhighfreq = 'yes' or 'no'
% cfg.freqwin = frequency window to be considered for phase-amplitude coupling, specified as [low high] im Hz (no default), in cell array with one cell per subject if cfg.individualhighfreq = 'yes'
% cfg.Vtestvalue = phase angle for which maximum power is expected in power-providing frequency, specified in angular degree (e.g., 0 = peak , 90 = falling flank, 180 = trough, 270 = rising flank) (no default)
% cfg.timelockchannel = string indicating the channel for which PAC_SI shall be calcualted, ERPs time-locked and TFR plotted (no default)
% [results] = sa_PAC_SI(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results] = sa_PAC_SI(gcfg)
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

for iSub = 1:numel(gcfg.subjects) % loop over subjects
    for iChan = 1:length(gcfg.files.eventdata{iSub}); % loop over channels
        
        % load trialwise data of subject i
        temp = load(gcfg.files.eventdata{iSub}{iChan}); % load freq data;
        tempfieldnames = fieldnames(temp);
        eventdata = eval(['temp.' tempfieldnames{~cellfun('isempty',strfind(tempfieldnames,'eventdata'))}]);
        clear temp tempfieldnames;
        temp = load(gcfg.files.freq{iSub}{iChan}); % load freq data
        tempfieldnames = fieldnames(temp);
        freq = eval(['temp.' tempfieldnames{~cellfun('isempty',strfind(tempfieldnames,'freq'))}]);
        clear temp tempfieldnames;
        
        % select trials
        switch gcfg.trials{iSub}{iChan}
            case 'all'
                % don't run trial selection
            otherwise
                cfg = [];
                cfg.trials = gcfg.trials{iSub}{iChan};
                cfg.tolerance = 1e-3;
                eventdata = ft_selectdata(cfg,eventdata);
                freq = ft_selectdata(cfg,freq);
        end
        
        % define data        
        switch gcfg.individualhighfreq 
            case 'no'
                freqindex = freq.freq >= gcfg.freqwin(1) & freq.freq <= gcfg.freqwin(2);
            case 'yes'
                freqindex = freq.freq >= gcfg.freqwin(iSub,1) & freq.freq <= gcfg.freqwin(iSub,2);           
        end        
        switch gcfg.individuallowfreq 
            case 'no'
                timeindex = unique(int16(find(eventdata.time{1} >= gcfg.timewin(1) & eventdata.time{1} <= gcfg.timewin(2))));
            case 'yes'
                timeindex = unique(int16(find(eventdata.time{1} >= gcfg.timewin(iSub,1) & eventdata.time{1} <= gcfg.timewin(iSub,2))));
        end
        
        trialN{iSub}{iChan} = size(eventdata.trial,2); % # trials
        for iTrial = 1:trialN{iSub}{iChan} % loop over trials
            
            % define data for trial j
            LF = eventdata.trial{iTrial};
            HF = squeeze(mean(freq.powspctrm(iTrial,1,freqindex,:),3))';
            
            % preprocess LF data
            switch gcfg.individuallowfreq
                case 'no'
                    LFhp = ft_preproc_bandpassfilter(LF, double(eventdata.fsample), [gcfg.hpfilter(1)  gcfg.lpfilter(1)], 3*fix(eventdata.fsample/gcfg.hpfilter(1)), 'fir', 'twopass'); % Mathilde's formula        
                case 'yes'
                    LFhp = ft_preproc_bandpassfilter(LF, double(eventdata.fsample), [gcfg.hpfilter(iSub)  gcfg.lpfilter(iSub)], 3*fix(eventdata.fsample/gcfg.hpfilter(iSub)), 'fir', 'twopass'); % Mathilde's formula
            end
            
            
            % preprocess HF data
            %             usfactor = double(eventdata.fsample / round(1/mean(diff(freq.time)))); %  upsampling factor to interpolate HF (TFR) signal to LF signal
            HF(isnan(HF)) = 0; % replace NaN by 0
            %             HFus = interp(HF,usfactor); % upsample to fit LF time resolution
            %             HFus = transpose(resample(transpose(HF),eventdata.fsample,round(1/mean(diff(freq.time))))); % upsample to fit LF time resolution
            HFus = transpose(resample(transpose(double(HF)), double(freq.time), double(eventdata.fsample), 'pchip')); % upsample to fit LF time resolution
            switch gcfg.individualhighfreq
                case 'no'
                    HFhp = ft_preproc_highpassfilter(HFus, double(eventdata.fsample), gcfg.hpfilter(1), 3*fix(eventdata.fsample/gcfg.hpfilter(1)), 'fir', 'twopass'); % Mathilde's formula
                case 'yes'
                    HFhp = ft_preproc_highpassfilter(HFus, double(eventdata.fsample), gcfg.hpfilter(iSub) , 3*fix(eventdata.fsample/gcfg.hpfilter(iSub)), 'fir', 'twopass'); % Mathilde's formula
            end
            % LF event-locked HF power change (difference)
            bl = [-2.0 -1.0]; % for SO it would be bl = [-2.5 -1.5]
            blindex = unique(int16(find(eventdata.time{1} >= bl(1) & eventdata.time{1} <= bl(2))-1));
            ERPow{iSub}{iChan}(iTrial) = mean(HFus(:,timeindex),2) - mean(HFus(:,blindex),2);                
            blPow{iSub}{iChan}(iTrial) = mean(HFus(:,blindex),2);
            
            % get phase angles (hilbert transform)
            LFph = angle(hilbert(LFhp(1,:)));
            HFph = angle(hilbert(HFhp(1,:)));
            
            % synchronization index (cf. Cohen et al., J Neurosci Methods 2008)
            phi1{iSub}{iChan}(iTrial,:) = LFph(:,timeindex); % phase values lf modulating signal
            phi2{iSub}{iChan}(iTrial,:) = HFph(:,timeindex); % phase values hf power modulation
            SI{iSub}{iChan}(iTrial) = mean(exp(1i*[phi1{iSub}{iChan}(iTrial,:) - phi2{iSub}{iChan}(iTrial,:)])); % Synchronization Index (SI)
            SIm{iSub}{iChan}(iTrial) = abs(SI{iSub}{iChan}(iTrial)); % SI magnitude (0 = perfect desychronization to 1 = perfect synchronization)
            SIp{iSub}{iChan}(iTrial) = (atan2(imag(SI{iSub}{iChan}(iTrial)),real(SI{iSub}{iChan}(iTrial)))); % preferred phase of the synchronization
            
        end % of loop over trials
        
        
        % calculate mean LF event-locked HF power change (percent change!)
        AERPow{iSub}{iChan} = mean(ERPow{iSub}{iChan})/mean(blPow{iSub}{iChan}) * 100;
        
        % two-sided one-sample t-test of percent LF event-timelocked HF power change against zero
        [AERPowh{iSub}{iChan},AERPowp{iSub}{iChan},AERPowci{iSub}{iChan},AERPowstats{iSub}{iChan}] = ttest([ERPow{iSub}{iChan}(:)], 0, 0.05, 'both');
        AERPowt{iSub}{iChan} = AERPowstats{iSub}{iChan}.tstat;
        
        % calculate mean vector length from preferred phase angles  (measure of consistency of SI preferred phase over LF events)
        r{iSub}{iChan} = circ_r(SIp{iSub}{iChan}');
        
        % calculate mean SI magnitude (measure of average SI magnitude over LF events, independent of their respective phases)
        ASIm{iSub}{iChan} = mean(SIm{iSub}{iChan});
        
        % calculate mean preferred phase (measure of average phase angle of preferred phase over LF events)
        ASIp{iSub}{iChan} = circ_mean(SIp{iSub}{iChan}');
        ASIp_ang{iSub}{iChan} = circ_rad2ang(ASIp{iSub}{iChan});
        
        % generate surrogates (from uniform distribution)
        for k = 1:10000 % iterations for surrogate distribution
            Sur_SIp{iSub}{iChan}(k,:) = rand(1,trialN{iSub}{iChan})*2*pi; % generate as many random 2pi values from uniform distribution as there are trials
            Sur_r(iSub,iChan,k) = circ_r(Sur_SIp{iSub}{iChan}(k,:)'); % calculate mean vector length of surrogates
            Sur_ASIp(iSub,iChan,k) = circ_mean(Sur_SIp{iSub}{iChan}(k,:)'); % calculate mean preferred phase for surrogates
        end
        
        % single-subject permutation test of mean vector length against surrogates from uniform distribution with same trial number
        ixr = find(sort(Sur_r(iSub,iChan,:)) > r{iSub}{iChan}, 1);
        if isempty(ixr), ixr = length(Sur_r(iSub,iChan,:)); end
        Sur_r_p{iSub}{iChan} = 1-(ixr/length(Sur_r(iSub,iChan,:)));
        z{iSub}{iChan} = (r{iSub}{iChan} - mean(Sur_r(iSub,iChan,:))) / std(Sur_r(iSub,iChan,:)); % z-transform actual r-value with respect to surrogate distribution
        z_p{iSub}{iChan} = normcdf(-abs(z{iSub}{iChan}),0,1); % transform z-value to p-value
        
        % Rayleigh test on non-uniform distribution of preferred phase angles within subject
        [Rp{iSub}{iChan},Rz{iSub}{iChan}] = circ_rtest([SIp{iSub}{iChan}]');
        
        % single-subject chi° test on non-uniform distribution
        %         [Chi_h{iSub}{iChan},Chi_p{iSub}{iChan},Chi_st{iSub}{iChan}] = chi2gof(SIp{iSub}{iChan},'edges',linspace(-pi,pi,18),'expected',length(SIp{1})*diff(linspace(-pi,pi,18)));
        
        % single-subject KS test on non-uniform distribution
        %         [KSh{iSub}{iChan},KSp{iSub}{iChan},KSk{iSub}{iChan},KSc{iSub}{iChan}] = kstest([SIp{iSub}{iChan}]',[[SIp{iSub}{iChan}]', unifcdf([SIp{iSub}{iChan}]',-pi,pi)]);
        
        % print SUBJECT results
        display(' '); display(['SUBJECT ' gcfg.subjectNames{iSub} ' RESULTS for channel ' gcfg.channelNames{gcfg.subjects(iSub)}{iChan}]);
        display('LF Event-Locked HF Power Change');
        display([gcfg.subjectNames{iSub} ' AERPow: ' num2str(AERPow{iSub}{iChan})]);
        display(sprintf('%s %s%0.10f', gcfg.subjectNames{iSub}, ['T-test t(' num2str(AERPowstats{iSub}{iChan}.df) ') (AERPowt): '], AERPowt{iSub}{iChan}));
        display(sprintf('%s %s%0.10f', gcfg.subjectNames{iSub}, 'T-test p (AERPowp): ', AERPowp{iSub}{iChan}));
        display('PAC Strength');
        display([gcfg.subjectNames{iSub} ' mean vecotr length (r): ' num2str(r{iSub}{iChan})]);
        display([gcfg.subjectNames{iSub} ' magnitude of synchronization index (ASIm): ' num2str(ASIm{iSub}{iChan})]);
        display(sprintf('%s %s%0.10f', gcfg.subjectNames{iSub} , 'Rayleigh''s z (Rz): ', Rz{iSub}{iChan}));
        display(sprintf('%s %s%0.10f', gcfg.subjectNames{iSub} , 'Rayleigh p (Rp): ', Rp{iSub}{iChan}));
        display(sprintf('%s %s%0.10f', gcfg.subjectNames{iSub} , 'Mean vector length permutation p (Sur_r_p{iSub}{iChan}): ', Sur_r_p{iSub}{iChan}));
        display(sprintf('%s %s%0.10f', gcfg.subjectNames{iSub} , 'Mean vector length permutation z (z{iSub}{iChan}): ', z{iSub}{iChan}));
        display(sprintf('%s %s%0.10f', gcfg.subjectNames{iSub} , 'Mean vector length permutation z-to-p (z_p{iSub}{iChan}): ', z_p{iSub}{iChan}));
        display('PAC Preferred Phase');
        display([gcfg.subjectNames{iSub} ' preferred phase (ASIp): ' num2str(ASIp{iSub}{iChan})]);
    
        % plot circular stats & surrogates
        switch gcfg.plotpersubject
            case 'yes'
                results.figurehandles.PAC_SI{iSub}{iChan} = figure; set(gcf,'name',gcfg.subjectNames{iSub});
                subplot(2,4,1); circ_plot([SIp{iSub}{iChan}]','pretty','o',true,'linewidth',2,'color','r');
                title({['distribution of preferred phase angles of ' num2str(trialN{iSub}{iChan}) ' trials'], ['mean angle: ' num2str(ASIp{iSub}{iChan}) ' rad, ' num2str(circ_rad2ang(ASIp{iSub}{iChan})) 'degree']});
                subplot(2,4,2); circ_plot([SIp{iSub}{iChan}]','hist',[],18,true,true,'linewidth',2,'color','r');
                title({['distribution of preferred phase angles of ' num2str(trialN{iSub}{iChan}) ' trials'], ['mean angle: ' num2str(ASIp{iSub}{iChan}) ' rad, ' num2str(circ_rad2ang(ASIp{iSub}{iChan})) 'degree']});
                subplot(2,4,3); circ_plot(squeeze(Sur_ASIp(iSub,iChan,:)),'pretty','o',true,'linewidth',2,'color','r');
                title({'surrogate distribution of preferred phase angles'});
                subplot(2,4,4); circ_plot(squeeze(Sur_ASIp(iSub,iChan,:)),'hist',[],18,true,true,'linewidth',2,'color','r');
                title({'surrogate distribution of preferred phase angles'});
                subplot(2,4,5:8); hist(squeeze(Sur_r(iSub,iChan,:)),round(k/100));
                hold on; ylim = get(gca,'ylim'); line([r{iSub}{iChan} r{iSub}{iChan}], ylim,'color','r'); hold off;
                title(['surrogate distribution of mean vector lengths with actual r-value ' num2str(r{iSub}{iChan}) ' (red)']);
        end        
    end % of loop over channels
    
    % average across channels per subject
    AERPow{iSub} = mean([AERPow{iSub}{:}]);
    ASIp{iSub} = circ_mean([ASIp{iSub}{:}]');
    ASIp_ang{iSub} = circ_rad2ang(ASIp{iSub});
    r{iSub} = mean([r{iSub}{:}]);
    trialN{iSub} = mean([trialN{iSub}{:}]);
    AERPowt{iSub} = mean([AERPowt{iSub}{:}]);
    AERPowp{iSub} = mean([AERPowp{iSub}{:}]);
    ASIm{iSub} = mean([ASIm{iSub}{:}]);
    Rz{iSub} = mean([Rz{iSub}{:}]);
    Rp{iSub} = mean([Rp{iSub}{:}]);
    z{iSub} = mean([z{iSub}{:}]);
    z_p{iSub} = mean([z_p{iSub}{:}]);
    
end % of loop over subjects

% average across channels per subject
Sur_r = squeeze(mean(Sur_r,2));

    
    

    %% group statistics    
    
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
    Vtestvalue = gcfg.Vtestvalue; % CHANGE!!!!! was 180

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
    results.figurehandles.PAC_SI_GA = figure; 
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

results.cfg = gcfg; % add cfg input info to output

%%%%%%%%%%%%%%%%%%%%%%%%%
end % of function