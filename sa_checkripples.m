%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checkripples
% by Til Ole Bergmann 2013
% last modified 2016/11/22 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% check ripple candidates 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sa_checkripples(gcfg, eventdata, freq)

if ismember('all', gcfg.browsechannel)    
    gcfg.browsechannel = eventdata.label;
end

% set parameters
plotchannels = find(ismember(eventdata.label, gcfg.browsechannel))';
plotevents = 1:length(eventdata.trial);
sliderMin = 0;
sliderMax = length(plotevents); 
sliderStep = [1,1] / (sliderMax - sliderMin);

% generate figure with callback function
h = figure('Name','event browser');
hb = uicontrol('Style','slider','Callback',@eventSlider);
set(hb, 'Min', sliderMin, 'Max', sliderMax, 'SliderStep', sliderStep, 'Value', sliderMin); 

% nested function eventSlider
function eventSlider(hObject,eventinput)
    step = round(get(hObject,'Value'));
    if step <= 1
        step = 1;
    elseif step >= length(plotevents); 
        step = length(plotevents);
    end
    event = plotevents(step);
    
    channum = 0;
    for chan = plotchannels
        channum = channum+1;
        
        % RAW
        sphandle = subplot(2,length(plotchannels),channum);
        pos = get(sphandle,'position');
        pos(3) = pos(3)*0.95;  
        set(sphandle, 'position', pos);      
        title(['channel: ' eventdata.label{chan}... 
            ', trial: ' num2str(event)... 
            ', stage: ' num2str(eventdata.trialinfo(event,1))...
            ', duration: ' num2str(eventdata.trialinfo(event,4)) ' s'...
            ', peak-to-peak amplitude: ' num2str(eventdata.trialinfo(event,9)) ' �V'...
            ]);   
        cla;hold on;
        plot(eventdata.time{1}, eventdata.trial{event}(chan,:),'k'); % raw signal
        for i = 1:length(eventdata.label) % loop over channels candidates
            if strfind(eventdata.label{i}, eventdata.label{chan}) & strfind(eventdata.label{i}, '_bp') % bp signal of current raw channel
                plot(eventdata.time{1}, eventdata.trial{event}(i,:),'b'); % bp signal
            elseif strfind(eventdata.label{i}, eventdata.label{chan}) & strfind(eventdata.label{i}, 'rms') % rms signal of current raw channel
                plot(eventdata.time{1}, eventdata.trial{event}(i,:),'g'); % rms signal
            end
        end
        xlim = [-2 2];
%         line([0 length(rmsdata.trial{event}(chan,:))],[gcfg.rmsampthresh(chan) gcfg.rmsampthresh(chan)],'Color','r','LineStyle',':'); % zero line
        hold off;
        
        % TFR
        subplot(2,length(plotchannels),length(plotchannels)+channum);
        cla;
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.channel = gcfg.browsechannel(channum); 
        cfg.trials = event;
%         cfg.baseline = 'no';
        cfg.baseline = 'yes';
        cfg.baseline = [-1 1];
        cfg.baselinetype  = 'relative';
        cfg.zlim = 'maxmin';
%         cfg.zlim = [0 200];
        cfg.hotkeys = 'yes';
        cfg.masknans = 'no';
        
        ft_singleplotTFR(cfg, freq);
        
    end    
end

end % of main function














