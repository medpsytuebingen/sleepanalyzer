%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% selectevents
% by Til Ole Bergmann 2015
% last modified 2016/11/25 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% selects events based on eventinfo parameters 
% and gives back event selection vector (esv)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function esv = sa_selectevents(gcfg)
display('Selecting events based on eventinfo parameters...');

    
%% define events

switch gcfg.eventtype
    case 'SO'
        gcfg.event.spec = [gcfg.eventtype '_' gcfg.SOspec '_' gcfg.SO_timelockevent];
    case 'spindle'
         gcfg.event.spec = [gcfg.eventtype '_' gcfg.spindle_timelockevent];
    case 'TW'
         gcfg.event.spec = [gcfg.eventtype '_' gcfg.TWspec '_' gcfg.TW_timelockevent];
end

gcfg.event.eventdata = [gcfg.eventtype '_eventdata'];

switch gcfg.timelockchannel
    case 'HC'
        gcfg.event.chan = gcfg.HCtargetchannel;
    otherwise
        gcfg.event.chan = gcfg.timelockchannel;
end

    
    
%% load individual data
for i = 1:numel(gcfg.subjects) % loop over subjects
    temp = load(fullfile(gcfg.resultsPath,gcfg.subjectNames{i},[gcfg.subjectNames{i} '_' gcfg.event.spec '_' gcfg.event.chan{gcfg.subjects(i)} '_eventdata.mat']), gcfg.event.eventdata);       

    %   1       2           3           4          5        6          7        8       9        10         11        12         13
    % 'stage' 'startTime' 'midTime' 'endTime' 'duration' 'maxTime' 'minTime' 'minAmp' 'maxAmp' 'p2pAmp' 'p2pTime' 'RMSmaxAmp' 'RMSmaxTime'
    event{i} = abs(temp.(gcfg.event.eventdata).trialinfo); % WARNING: negative values become positive! High positive values mean large troughs!!
end
   

%% generate esv depending on logical tests 

for i = 1:numel(event) % loop over subjects

    % prepare criteria 
    
    n = size(event{i},1);
        
    trough_75_percentile = 75;
    peak_75_percentile = 75;      
    trough_25_percentile = 25;
    peak_25_percentile = 25;          
    p2p_75_percentile = 75;    
    p2p_90_percentile = 90;
    p2p_95_percentile = 95;
    
    sorttrough = sort(event{i}(:,8));
    sortpeak = sort(event{i}(:,9));
    sortp2p = sort(event{i}(:,10));

    idxtrough_25 = max(1,round(n/100*trough_25_percentile));   
    idxtrough_75 = max(1,round(n/100*trough_75_percentile));
    idxpeak_25 = max(1,round(n/100*peak_25_percentile));
    idxpeak_75 = max(1,round(n/100*peak_75_percentile));
    idxp2p_75 = max(1,round(n/100*p2p_75_percentile));    
    idxp2p_90 = max(1,round(n/100*p2p_90_percentile));
    idxp2p_95 = max(1,round(n/100*p2p_95_percentile));
    
    trough_25_crit = sorttrough(idxtrough_25);
    trough_75_crit = sorttrough(idxtrough_75);
    peak_25_crit = sortpeak(idxpeak_25);
    peak_75_crit = sortpeak(idxpeak_75);
    p2p_75_crit = sortp2p(idxp2p_75);
    p2p_90_crit = sortp2p(idxp2p_90);    
    p2p_95_crit = sortp2p(idxp2p_95);
    
       
%     for i=1:12
%         p2p_percentile{i} = find(sort(event{i}(:,10)) > 75, 1, 'first') / size(event{i},1);
%         
%         S = sort(event{i}(:,10));
%         p2p_75perc_value{i} = S(round(0.75*size(event{i},1)));
%     end
%     [p2p_percentile{:}];
%     [p2p_75perc_value{:}];


    % test criteria 
    esv{i} = zeros(1,size(event{i},1)); % preset esv    
    for j = 1:size(event{i},1) % loop over events                                     
        ok = 0;
        
        switch gcfg.selectioncriteria
            
            % SOs             
            case 'amp_p2p_>75%' % was default
                if event{i}(j,10) >= p2p_75_crit
                   ok = 1;    
                end
            case 'amp_peak_trough_>75%_>75%'
                if event{i}(j,8) >= trough_75_crit && event{i}(j,9) >= peak_75_crit
                   ok = 1;           
                end        
            case 'amp_trough_peak_p2p_>25µV_>25µV_>75µV'                
                if event{i}(j,8) >= 25 && event{i}(j,9) >= 25 && event{i}(j,10) >= 75 
                   ok = 1;           
                end 
            case 'amp_p2p_>90%' 
                if event{i}(j,10) >= p2p_90_crit
                   ok = 1;           
                end                   
            case 'amp_p2p_<90%'
                if event{i}(j,10) < p2p_90_crit
                    ok = 1;
                end                
            case 'amp_p2p_>75%_<95%'
                if event{i}(j,10) >= p2p_75_crit && event{i}(j,10) < p2p_95_crit
                    ok = 1;
                end                                
            case 'amp_p2p_>75%_dur_400msPerHalfWave'
                if event{i}(j,10) >= p2p_75_crit &&  (event{i}(j,3) - event{i}(j,2))*1000 >= 0.4 && (event{i}(j,4) - event{i}(j,3))*1000 >= 0.4
                    ok = 1;
                end               
            case 'amp_conservative'
                if event{i}(j,10) >= p2p_75_crit && event{i}(j,10) < p2p_95_crit && event{i}(j,8) >= trough_25_crit && event{i}(j,9) >= peak_25_crit && (event{i}(j,3) - event{i}(j,2))*1000 >= 0.4 && (event{i}(j,4) - event{i}(j,3))*1000 >= 0.4
                    ok = 1;
                end               
            case 'dur_400msPerHalfWave'
                if (event{i}(j,3) - event{i}(j,2))*1000 >= 0.4 && (event{i}(j,4) - event{i}(j,3))*1000 >= 0.4
                    ok = 1;
                end                           
                
            % spindles  
            case 'amp_p2p_10µV'
                if event{i}(j,10) >= 10 
                   ok = 1;           
                end 
                
                
        end % of switch
                
        
        esv{i}(j) = ok;
    end % of loop over events
end % of loop over subjects
    


end % of function




