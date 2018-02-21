%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sa_plotevent
% by Til Ole Bergmann 2013
% last modified 2017/10/04 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plots data around events from EEG data
%
% pe = plotevent structure
% pe.figurehandle = handle to figure
%
% cfg = [];
% [pe] = sa_inspectsleep(gcfg,data);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pe] = sa_plotevent(varargin)
display('Plotting events...');

switch nargin
    case 1
        gcfg = varargin{1};
        display('Not enough input arguments specified!');
        return
    case 2
        gcfg = varargin{1};        
        timelock = varargin{2};
        nrows = 1;
	case 3
        gcfg = varargin{1};        
        timelock = varargin{2};
        freq = varargin{3};
        nrows = 2;
end

% define channels to be plotted
if strcmp(gcfg.plotchannel,'all')
    if strcmp(gcfg.channeltype,'raw'); % select raw channels only
        procchannels = {'_bp','_rms'};
        goodchan =[];
        for i = 1:length(timelock.label)
            for j = 1:length(procchannels) % loop over channels postfixes to exclude
                k = strfind(timelock.label{i},procchannels{j});
                goodchan(i,j) = isempty(k);
            end
        end
        goodchan = logical(prod(goodchan,2));
        gcfg.plotchannel = timelock.label(logical(goodchan));
        
    else
        gcfg.plotchannel = timelock.label;
    end
end


% calculate indices
timelock.sem = sqrt(timelock.var)./sqrt(timelock.dof); % calculate SEM

% define colors
switch gcfg.colorset
    case 'allblack' 
        colortable = zeros(length(timelock.label),3);
    case 'mycolor' 
        for i = 1:length(timelock.label)
            switch timelock.label{i}
                case 'Cz'
                   colortable(i,:) = [0 0 1]; 
                case 'HC'
                   colortable(i,:) = [0 1 0];    
                case 'HCA'
                   colortable(i,:) = [0 0.75 0];   
                case 'HCP'
                   colortable(i,:) = [0 0.5 0];                      
                case 'RC'
                   colortable(i,:) = [1 0 0];                                         
%                 case 'Cz_SO_bp'
%                    colortable(i,:) = [0 0 0.5]; 
%                 case 'HC_SO_bp'
%                    colortable(i,:) = [0 0.5 0];    
                otherwise
                   colortable(i,:) = [0 0 0];                    
            end
        end
    case 'varycolor' 
%         colortable = circshift(varycolor(length(gcfg.plotchannel)),2);        
        colortable = circshift(varycolor(length(timelock.label)),2);        
    case 'distinguishable_colors' 
        colortable = distinguishable_colors(length(timelock.label));
    case 'nautic'
        colortable = zeros(length(gcfg.plotchannel),3);
        TL_chanNum = 10; %sum(cell2mat(strfind(gcfg.plotchannel,'TL')));
        TR_chanNum = 10; %sum(cell2mat(strfind(gcfg.plotchannel,'TR')));        
        for i = 1:length(gcfg.plotchannel)
            if ~isempty(strfind(gcfg.plotchannel{i},'Cz'))
                colortable(i,:) = [0 0 1];
            elseif ~isempty(strfind(gcfg.plotchannel{i},'TL'))
                k = strfind(gcfg.plotchannel{i},'TL');
                chan = str2num(gcfg.plotchannel{i}(k+2:k+3));
                colortable(i,:) = [0+(chan/TL_chanNum) 0 0];
            elseif ~isempty(strfind(gcfg.plotchannel{i},'TR'))
                k = strfind(gcfg.plotchannel{i},'TR');
                chan = str2num(gcfg.plotchannel{i}(k+2:k+3));
                colortable(i,:) = [0 0+(chan/TR_chanNum) 0];               
            end
        end
    case 'parula'
        ct = parula;
        colortable = ct(1:length(ct)/length(gcfg.plotchannel):end,:);
    case 'jet'
        ct = jet;
        colortable = ct(1:length(ct)/length(gcfg.plotchannel):end,:);       
    case 'NNS1' % as in Staresina et al., NN 2015, Figure S1
        colortable =   [0,0,0;...
                        246,157,169;...
                        128,29,36;...
                        191,29,42;...
                        200,56,149;...
                        61,25,90;...
                        30,65,153;...
                        20,187,172;...
                        49,179,71;...
                        253,183,18;...
                        237,28,33];
        colortable = colortable/255;
    case 'customcolortable' % custom colortable
        colortable = gcfg.colortable;                        
end % of switch

% !!!Insert to plot some filer related stuff!!!
% keyboard;
% nrows = 3;
% hpfilter = 12;
% freqwin = [80 90];
% timelock_fsample = round(1/mean(diff(timelock.time)));
% freq_fsample = round(1/mean(diff(freq.time)));
% usfactor =  timelock_fsample / freq_fsample;
% 
% % timelocked average of hf power modulation
% sphandle = subplot(nrows,1,3);
% if exist('freq', 'var')
%     pos = get(sphandle,'position');
%     pos(3) = pos(3)*0.85;  
%     set(sphandle, 'position', pos);    
% end
% % hold on;
% % linecount = 0;
% % if isfield(gcfg, 'xlim') % reduce dataset to plot  
% %     HFpow.avg = HFpow.avg(:, timelock.time >= gcfg.xlim(1) & timelock.time <= gcfg.xlim(2));
% %     HFpow.var = HFpow.var(:, timelock.time >= gcfg.xlim(1) & timelock.time <= gcfg.xlim(2));
% %     HFpow.sem = HFpow.sem(:, timelock.time >= gcfg.xlim(1) & timelock.time <= gcfg.xlim(2));
% %     HFpow.dof = HFpow.dof(:, timelock.time >= gcfg.xlim(1) & timelock.time <= gcfg.xlim(2));
% %     HFpow.time = HFpow.time(:, timelock.time >= gcfg.xlim(1) & timelock.time <= gcfg.xlim(2)); % do this last!    
% % end  
% 
% freqindex = freq.freq >= freqwin(1) & freq.freq <= freqwin(2);  
% HF = squeeze(mean(freq.powspctrm(1,freqindex,:),2))';
% HF(isnan(HF)) = 0; % replace NaN by 0
% HFus = interp(HF,usfactor); % upsample to fit LF time resolution
% HFhp = ft_preproc_highpassfilter(HFus, double(timelock_fsample ), hpfilter, 3*fix(timelock_fsample /hpfilter)+1, 'fir', 'twopass'); % Mathilde's formula        
% plot(HFhp, 'r');

pe.figurehandle = figure; 
set(pe.figurehandle,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);

% timelocked average
if exist('freq', 'var')
    sphandle = subplot(nrows,1,2);
else
    sphandle = subplot(nrows,1,1);
end
if exist('freq', 'var')
    pos = get(sphandle,'position');
    pos(3) = pos(3)*0.955;  
    set(sphandle, 'position', pos);    
end
hold on;
linecount = 0;
if isfield(gcfg, 'xlim') % reduce dataset to plot  
    timelock.avg = timelock.avg(:, timelock.time >= gcfg.xlim(1) & timelock.time <= gcfg.xlim(2));
    timelock.var = timelock.var(:, timelock.time >= gcfg.xlim(1) & timelock.time <= gcfg.xlim(2));
    timelock.sem = timelock.sem(:, timelock.time >= gcfg.xlim(1) & timelock.time <= gcfg.xlim(2));
    timelock.dof = timelock.dof(:, timelock.time >= gcfg.xlim(1) & timelock.time <= gcfg.xlim(2));
    timelock.time = timelock.time(:, timelock.time >= gcfg.xlim(1) & timelock.time <= gcfg.xlim(2)); % do this last!    
end  
for chan = find(ismember(timelock.label, gcfg.plotchannel))'
    linecount = linecount+1;
    shadeHandle = fill([timelock.time, fliplr(timelock.time)],[timelock.avg(chan,:)+timelock.sem(chan,:), fliplr(timelock.avg(chan,:)-timelock.sem(chan,:))],colortable(chan,:));
    set(shadeHandle,'edgecolor',colortable(chan,:));
    alpha(0.25);  
        lineHandle(linecount) = plot(timelock.time,timelock.avg(chan,:), 'Color', colortable(chan,:), 'LineWidth', 2);    
    if chan == find(strcmp(gcfg.targetchannel,timelock.label)) % if targetchannel, emphasize with black dashes
        lineHandle(linecount) = plot(timelock.time,timelock.avg(chan,:), 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle','--');
    end
end
line([timelock.time(1) timelock.time(end)], [0 0], 'Color', 'k', 'LineWidth', 1);
% y = get(gca,'ylim');
% line([0 0], y, 'Color', 'k', 'LineWidth', 1, 'LineStyle','--');
hold off;
if isfield(gcfg,'tl_ylim')
    ylim(gcfg.tl_ylim);
end
% legend(lineHandle,timelock.label,'Location','NorthOutside','Orientation','horizontal','Interpreter','none');
legend(lineHandle,gcfg.plotchannel,'Location','NorthOutside','Orientation','horizontal','Interpreter','none');
title(['Average of ' num2str(max(unique(timelock.dof))) ' observations']);

% averaged TFR
if exist('freq', 'var')
    subplot(nrows,1,1);
    cfg = [];
    cfg.channel = gcfg.plotTFRchannel;
    if isfield(gcfg, 'baseline')
        cfg.baseline = gcfg.baseline;
    else
        cfg.baseline = 'yes';    
    end
    if isfield(gcfg, 'baselinetype')
        cfg.baselinetype = gcfg.baselinetype;
    else
        cfg.baselinetype = 'relative';    
    end

    if isfield(gcfg, 'zlim')
        cfg.zlim = gcfg.zlim;
    else
        cfg.zlim = 'maxabs';

    end
    if isfield(gcfg, 'ylim')
        cfg.ylim = gcfg.ylim;
    end
    if isfield(gcfg, 'xlim')
        cfg.xlim = gcfg.xlim;
    end  

    ft_singleplotTFR(cfg, freq);
    
end
pe.cfg = gcfg;
end % of function
