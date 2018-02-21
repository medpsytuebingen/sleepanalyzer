%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epochbrowser
% by Til Ole Bergmann 2013
% last modified 2016/11/22 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% browse through epochs with detected events
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sa_epochbrowser(varargin)
% epochbrowser(gcfg, data, bpdata, rmsdata, supthreshdata)
switch nargin
    case 1
        gcfg = varargin{1};
        display('Not enough input arguments specified!');
        return
    case 2
        gcfg = varargin{1};        
        data = varargin{2};
        bpdata_flag = 0; rmsdata_flag = 0; supthreshdata_flag = 0;
	case 3
        gcfg = varargin{1};        
        data = varargin{2};
        bpdata = varargin{3};
        bpdata_flag = 1; rmsdata_flag = 0; supthreshdata_flag = 0;
	case 4
        gcfg = varargin{1};        
        data = varargin{2};
        bpdata = varargin{3};
        rmsdata = varargin{4};
        bpdata_flag = 1; rmsdata_flag = 1; supthreshdata_flag = 0;
	case 5
        gcfg = varargin{1};        
        data = varargin{2};
        bpdata = varargin{3};
        rmsdata = varargin{4};        
        supthreshdata = varargin{4}; 
        bpdata_flag = 1; rmsdata_flag = 1; supthreshdata_flag = 1;
end


if ismember('all', gcfg.browsechannel)    
    gcfg.browsechannel = data.label;
end

% set parameters
plotchannels = find(ismember(data.label, gcfg.browsechannel))';

if numel(data.trial) == 1 % if continuous data
    windowSize = gcfg.browseWindowLength*data.fsample;
    windowBegins = [1:windowSize:size(data.trial{1},2)];
    stageFilter = ismember(gcfg.scoring(1,:), gcfg.browsestages);
    plotWindowBegins = windowBegins(ismember(gcfg.scoring(1,windowBegins), gcfg.browsestages));
    if bpdata_flag && rmsdata_flag 
        ylimAll = [min([data.trial{1}(plotchannels,stageFilter) bpdata.trial{1}(plotchannels,stageFilter) rmsdata.trial{1}(plotchannels,stageFilter)],[],2) max([data.trial{1}(plotchannels,stageFilter) bpdata.trial{1}(plotchannels,stageFilter) rmsdata.trial{1}(plotchannels,stageFilter)],[],2)];
    elseif bpdata_flag && rmsdata_flag == 0
        ylimAll = [min([data.trial{1}(plotchannels,stageFilter) bpdata.trial{1}(plotchannels,stageFilter)],[],2) max([data.trial{1}(plotchannels,stageFilter) bpdata.trial{1}(plotchannels,stageFilter)],[],2)];        
    elseif bpdata_flag == 0 && rmsdata_flag == 0
        ylimAll = [min([data.trial{1}(plotchannels,stageFilter)],[],2) max([data.trial{1}(plotchannels,stageFilter)],[],2)];                
    end
    sliderMin = 0;
    sliderMax = length(plotWindowBegins); 
    sliderStep = [1,1] / (sliderMax - sliderMin);

elseif numel(data.trial) > 1 % if epoched data
    plotepochs = find(ismember(gcfg.scoring(1,:), gcfg.browsestages))';
    if bpdata_flag && rmsdata_flag     
        ylimAll = [min([data.trial{:} bpdata.trial{:} rmsdata.trial{:}],[],2) max([data.trial{:} bpdata.trial{:} rmsdata.trial{:}],[],2)];
    elseif bpdata_flag && rmsdata_flag == 0
        ylimAll = [min([data.trial{:} bpdata.trial{:}],[],2) max([data.trial{:} bpdata.trial{:}],[],2)];
	elseif bpdata_flag == 0 && rmsdata_flag == 0
        ylimAll = [min([data.trial{:}],[],2) max([data.trial{:}],[],2)];
    end    
    sliderMin = 0;
    sliderMax = length(plotepochs); 
    sliderStep = [1,1] / (sliderMax - sliderMin);
end

% generate figure with callback function
h = figure('Name','epoch browser');
if numel(data.trial) == 1 % if continuous data
    hb = uicontrol('Style','slider','Callback',@windowSlider);
elseif numel(data.trial) > 1 % if epoched data
    hb = uicontrol('Style','slider','Callback',@epochSlider);
end
set(hb, 'Min', sliderMin, 'Max', sliderMax, 'SliderStep', sliderStep, 'Value', sliderMin); 


% nested function windowSlider
function windowSlider(hObject,epochdata)
    step = round(get(hObject,'Value'));
    if step <= 1
        step = 1;
    elseif step >= length(plotWindowBegins); 
        step = length(plotWindowBegins);
    end
    window = [plotWindowBegins(step):plotWindowBegins(step)+windowSize]; 
    windowNum = round(plotWindowBegins(step) / windowSize) + 1;
    switch gcfg.scoring(1,plotWindowBegins(step))
        case 0, stage = 'WAKE';
        case 1, stage = 'NREM S1';
        case 2, stage = 'NREM S2';
        case 3, stage = 'NREM S3';
        case 4, stage = 'NREM S4';
        case 5, stage = 'REM';            
        case 6, stage = 'MT';            
        otherwise stage = 'UNKNOWN';            
    end  
    for chan = plotchannels;
        subplot(length(plotchannels),1,chan);
        cla;hold on;
        plot(data.trial{1}(chan,window),'k');
        if bpdata_flag
            plot(bpdata.trial{1}(chan,window),'b');
        end
        if rmsdata_flag
            plot(rmsdata.trial{1}(chan,window),'g');
            line([0 length(rmsdata.trial{1}(chan,window))],[gcfg.rmsampthresh(chan) gcfg.rmsampthresh(chan)],'Color','r','LineStyle',':'); % zero line
        end        
        if supthreshdata_flag
            plotsupthresh = double(supthreshdata{1}(chan,window));
            plotsupthresh(find(plotsupthresh == 0)) = NaN;
            plot(plotsupthresh .* gcfg.rmsampthresh(chan),'Color','m','LineStyle','-','LineWidth',5);
        end        
        plotartifact = double(~data.artifactFilter{1}(chan,window)); %!!!!!!!!!! invertieren!!!!
        plotartifact(find(plotartifact == 0)) = NaN;
        plot(plotartifact,'Color','r','LineStyle','-','LineWidth',3);        
        
        ylim(ylimAll(chan,:));
        xlim([1 windowSize]);
        title([data.label{chan} ' (bpfilter: ' num2str(gcfg.bpfreq(chan,1)) '-' num2str(gcfg.bpfreq(chan,2)) ' Hz) , trial ' num2str(windowNum) ', ' stage]);   
        hold off;
    end    
end


% nested function epochSlider
function epochSlider(hObject,epochdata)
    step = round(get(hObject,'Value'));
    if step <= 1
        step = 1;
    elseif step >= length(plotepochs); 
        step = length(plotepochs);
    end
    epoch = plotepochs(step);
    switch gcfg.scoring(1,epoch)
        case 0, stage = 'WAKE';
        case 1, stage = 'NREM S1';
        case 2, stage = 'NREM S2';
        case 3, stage = 'NREM S3';
        case 4, stage = 'NREM S4';
        case 5, stage = 'REM';            
        case 6, stage = 'MT';            
        otherwise stage = 'UNKNOWN';            
    end  
    for chan = plotchannels;
        subplot(length(plotchannels),1,chan);
        cla;hold on;
        plot(data.trial{epoch}(chan,:),'k');
        if bpdata_flag
            plot(bpdata.trial{epoch}(chan,:),'b');
        end
        if rmsdata_flag
            plot(rmsdata.trial{epoch}(chan,:),'g');
            line([0 length(rmsdata.trial{epoch}(chan,:))],[gcfg.rmsampthresh(chan) gcfg.rmsampthresh(chan)],'Color','r','LineStyle',':'); % zero line
        end        
        if supthreshdata_flag        
            plotsupthresh = double(supthreshdata{epoch}(chan,:));
            plotsupthresh(find(plotsupthresh == 0)) = NaN;
            plot(plotsupthresh .* gcfg.rmsampthresh(chan),'Color','m','LineStyle','-','LineWidth',5);
        end
        ylim(ylimAll(chan,:));
        title([data.label{chan} ' (bpfilter: ' num2str(gcfg.bpfreq(chan,1)) '-' num2str(gcfg.bpfreq(chan,2)) ' Hz) , trial ' num2str(epoch) ', ' stage]);   
        hold off;
    end    
end

end % of main function