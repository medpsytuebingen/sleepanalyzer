%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checkspindles
% by Til Ole Bergmann 2013
% last modified 2016/11/22 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% check spindle candidates 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sa_checkspindles(gcfg, eventdata, freq)

if ismember('all', gcfg.browsechannel)    
    gcfg.browsechannel = eventdata.label;
end

% normalize freq data
normfreq = freq;
[a b c d] = size(freq.powspctrm);
for trial = 1:a
    baseline = repmat(nanmean(freq.powspctrm(trial,:,:,:),4),[1,1,1,d]);
    normfreq.powspctrm(trial,:,:,:) = (freq.powspctrm(trial,:,:,:) - baseline) ./ baseline .* 100;
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

figSize = get(h, 'Position');
hb2 = uicontrol('Style','togglebutton','Position',[figSize(3)-70 20 50 20], 'String', 'reject', 'Callback',@rejectEvent); 

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
        sphandle = subplot(4,length(plotchannels),channum);
%         sphandle = subplot(3,length(plotchannels),channum);
        pos = get(sphandle,'position');
        pos(3) = pos(3)*0.92; 
        set(sphandle, 'position', pos);           
        if channum == 1
            extratext = [', trial: ' num2str(event)... 
            ', stage: ' num2str(eventdata.trialinfo(event,1))...
            ', dur: ' num2str(eventdata.trialinfo(event,4)) ' s'...
            ', amp: ' num2str(eventdata.trialinfo(event,9)) ' �V'...
            ];   
        else
            extratext = '';
        end       
        title(['channel: ' eventdata.label{chan} extratext]); 
        cla; hold on;
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
        subplot(4,length(plotchannels),length(plotchannels)+channum);
        cla;
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.channel = gcfg.browsechannel(channum); 
        cfg.trials = event;
        cfg.baseline = 'no';
%         cfg.baseline = 'yes';
%         cfg.baseline = [-1 -0.5];
%         cfg.baselinetype  = 'relative';
        cfg.zlim = 'maxabs';
        cfg.hotkeys = 'yes';
        cfg.masknans = 'yes';
        ft_singleplotTFR(cfg, normfreq);        
        
        % TFR REDUCED
        subplot(4,length(plotchannels),2*length(plotchannels)+channum);
        cla;
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.channel = gcfg.browsechannel(channum); 
        cfg.trials = event;
        cfg.baseline = 'no';
%         cfg.baseline = 'yes';
%         cfg.baseline = [-1 -0.5];
%         cfg.baselinetype  = 'absolute';
        cfg.zlim = 'maxabs';
        cfg.ylim = [0 30];
        cfg.hotkeys = 'yes';
        cfg.masknans = 'yes';
        ft_singleplotTFR(cfg, normfreq);
        
        % frequency profile
        sphandle = subplot(4,length(plotchannels),3*length(plotchannels)+channum);
        pos = get(sphandle,'position');
        pos(3) = pos(3)*0.92;
        set(sphandle, 'position', pos);
        cla;
        timewin = [-0.05 0.05]; % time window in s to average 
        plotdata_raw = (squeeze(mean(freq.powspctrm(event,chan,:,freq.time > timewin(1) & freq.time < timewin(2)),4)))'; % raw      
        plotdata_log = (log(squeeze(mean(freq.powspctrm(event,chan,:,freq.time > timewin(1) & freq.time < timewin(2)),4))+1).*100)'; % log correction
        plotdata_1overf = (squeeze(mean(freq.powspctrm(event,chan,:,freq.time > timewin(1) & freq.time < timewin(2)),4)).*freq.freq'./10)'; % 1/f correction
        plotdata_norm = (squeeze(mean(normfreq.powspctrm(event,chan,:,normfreq.time > timewin(1) & normfreq.time < timewin(2)),4)))'; % % change      
        hold on;
        plot(freq.freq, plotdata_raw, 'b');
        plot(freq.freq, plotdata_log, 'g');
        plot(freq.freq, plotdata_1overf, 'r');
        plot(normfreq.freq, plotdata_norm, 'm');
        hold off;
        title(['frequency profile averaged from ' num2str(timewin(1)) ' to ' num2str(timewin(2)) ' s: raw (blue), log(raw+1).*100 (green), raw.*freq (red), % change (magenta)']); 
 
    end    
end

% nested function rejectEvent
function varargout = rejectEvent(h, eventdata, handles, varargin)
        button_state = get(h,'Value');
        if button_state == get(h,'Max')
            % toggle button is pressed
        elseif button_state == get(h,'Min')
            % toggle button is not pressed
        end
end % of rejectEvent function

end % of main function