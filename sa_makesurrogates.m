%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makesurrogates
% by Til Ole Bergmann 2013
% last modified 2016/11/22 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% generates surrogate data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [surrogate_freq surrogate_timelock surrogate_fft] = sa_makesurrogates(gcfg, eventdata)
display('Generating surrogates...');

% select raw channels only
procchannels = {'_bp','_rms'};
goodchan =[];
for i = 1:length(eventdata.label)
    for j = 1:length(procchannels) % loop over channels postfixes to exclude
        k = strfind(eventdata.label{i},procchannels{j});  
        goodchan(i,j) = isempty(k);
    end
end
goodchan = logical(prod(goodchan,2));
rawchannels = eventdata.label(logical(goodchan));

eventdata_old = eventdata;

%% generate surrogates
segmentlength = numel(eventdata.time{1});
window = ones(segmentlength,1); 
% window = hann(segmentlength);

for i = 1:numel(eventdata.trial) % over trials        
    for j = ismember(eventdata.label, rawchannels) % over channels
        
        % preprocess data (!!! multo importante !!!)
        eventdata.trial{i}(j,:) = detrend(eventdata.trial{i}(j,:),'linear'); % demean and detrend       
%         eventdata_old.trial{i}(j,:) = detrend(eventdata_old.trial{i}(j,:),'linear'); % demean and detrend       
                
        fourier = fft(window.*eventdata.trial{i}(j,:)');
        fft_ang = angle(fourier);
        fft_abs = abs(fourier);
        
        % add noise to angles
        randval = rand(size(fourier)).*(2*pi);
        ang_new = fft_ang+randval;
        fouriernew = fft_abs.*(cos(ang_new)+1i*sin(ang_new));
        fouriernew([1 end]) = fourier([1 end]); % importance pretty unclear...
        
        % recreate symmetry of fourier transformed data (!!! multo importante !!!)
        for sample = 2:segmentlength/2
            fouriernew(sample) = real(fouriernew(segmentlength-sample+2)) - imag(fouriernew(segmentlength-sample+2))*1i;
        end
        
        % transform back
        trldatnew = ifft(fouriernew);
        
        % replace original trial values
        eventdata.trial{i}(j,:) = real(trldatnew);
        
    end
end

%% FFT of surrogates and actual data to check for differences
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'rectwin'; % 'hanning' produces differences!
cfg.channel = rawchannels; 
cfg.keeptrials = 'no';
cfg.foi = [5:1:200]; 
% cfg.pad = 4;  % padding produces differences!
cfg.output = 'pow';	
cfg.polyremoval = 1; % 0 = mean (default), 1 = linear
eval([gcfg.eventtype '_surrogate_fft = ft_freqanalysis(cfg, eventdata);']);
eval([gcfg.eventtype '_fft = ft_freqanalysis(cfg, eventdata_old);']);


% plot and save figures 
switch gcfg.eventtype
    case 'SO'
        figure;hold on;
            plot(SO_fft.freq, SO_fft.powspctrm); 
            plot(SO_fft.freq, SO_surrogate_fft.powspctrm,'r'); 
        hold off;       
        set(gcf, 'Name',[gcfg.subjectName ' surrogate spectral power check'],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
        filename = fullfile(gcfg.resultspath, [gcfg.subjectName '_SO_surrogate_power_check']);
        saveas(gcf,filename,'fig');
        
	case 'spindle'
        figure;hold on;
            plot(spindle_fft.freq, spindle_fft.powspctrm); 
            plot(spindle_fft.freq, spindle_surrogate_fft.powspctrm,'r'); 
        hold off; 
        set(gcf, 'Name',[gcfg.subjectName ' surrogate spectral power check'],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
        filename = fullfile(gcfg.resultspath, [gcfg.subjectName '_spindle_surrogate_power_check']);
        saveas(gcf,filename,'fig');
        
	case 'ripple'
        figure;hold on;
            plot(ripple_fft.freq, ripple_fft.powspctrm); 
            plot(ripple_fft.freq, ripple_surrogate_fft.powspctrm,'r'); 
        hold off;         
        set(gcf, 'Name',[gcfg.subjectName ' surrogate spectral power check'],'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
        filename = fullfile(gcfg.resultspath, [gcfg.subjectName '_ripple_surrogate_power_check']);
        saveas(gcf,filename,'fig');        
end
close all;


%% TFR of surrogates
windowlength = (size(eventdata.trial{1},2)-1)/eventdata.fsample; % length of time window

cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.channel = rawchannels; 
cfg.keeptrials = 'yes';
switch gcfg.eventtype
    case 'spindle'
        cfg.toi = [-windowlength/2:0.005:windowlength/2];
        cfg.polyremoval = 1; % 0 = mean (default), 1 = linear
    case 'SO'
        cfg.toi = [-windowlength/2:0.01:windowlength/2];
        cfg.polyremoval = 0; % 0 = mean (default), 1 = linear
end
cfg.foi = [5:1:200]; 
for i = 1:200, cycles(i) = floor(100/(1000/i)); end % keep windows about 100ms with integer number of cycles
cycles(cycles < 5) = 5; % but at least 5 cycles!
cfg.t_ftimwin = cycles(cfg.foi(1):cfg.foi(end))./cfg.foi;
% cfg.t_ftimwin = 5./cfg.foi;
cfg.output = 'pow';	
eval([gcfg.eventtype '_surrogate_freq = ft_freqanalysis(cfg, eventdata);']);


%% timelocked averages of surrogates
cfg = [];
cfg.channel = 'all';
eval([gcfg.eventtype '_surrogate_timelock = ft_timelockanalysis(cfg, eventdata);']);


%% rename surrogate data for output of function
switch gcfg.eventtype
    case 'SO'
        surrogate_freq = SO_surrogate_freq; clear SO_surrogate_freq;
        surrogate_timelock = SO_surrogate_timelock; clear SO_surrogate_timelock;
        surrogate_fft = SO_surrogate_fft; clear SO_surrogate_fft; 
    case 'spindle'
        surrogate_freq = spindle_surrogate_freq; clear spindle_surrogate_freq;
        surrogate_timelock = spindle_surrogate_timelock; clear spindle_surrogate_timelock;
        surrogate_fft = spindle_surrogate_fft; clear spindle_surrogate_fft; 
    case 'ripple'
        surrogate_freq = ripple_surrogate_freq; clear ripple_surrogate_freq;
        surrogate_timelock = ripple_surrogate_timelock; clear ripple_surrogate_timelock;
        surrogate_fft = ripple_surrogate_fft; clear ripple_surrogate_fft; 
end


% %% save surrogate data 
% switch gcfg.eventtype
%     case 'SO'
%         save('-v7.3',fullfile(gcfg.resultspath, [gcfg.subjectName '_SO_surrogate_freq.mat']),'SO_surrogate_freq');
%         save('-v7.3',fullfile(gcfg.resultspath, [gcfg.subjectName '_SO_surrogate_timelock.mat']),'SO_surrogate_timelock');
%         save('-v7.3',fullfile(gcfg.resultspath, [gcfg.subjectName '_SO_surrogate_fft.mat']),'SO_surrogate_fft');
%     case 'spindle'
%         save('-v7.3',fullfile(gcfg.resultspath, [gcfg.subjectName '_spindle_surrogate_freq.mat']),'spindle_surrogate_freq');
%         save('-v7.3',fullfile(gcfg.resultspath, [gcfg.subjectName '_spindle_surrogate_timelock.mat']),'spindle_surrogate_timelock');
%         save('-v7.3',fullfile(gcfg.resultspath, [gcfg.subjectName '_spindle_surrogate_fft.mat']),'spindle_surrogate_fft');
%     case 'ripple'
%         save('-v7.3',fullfile(gcfg.resultspath, [gcfg.subjectName '_ripple_' gcfg.ripple_timelockevent '_' timelockchannel '_surrogate_freq.mat']),'ripple_surrogate_freq');
%         save('-v7.3',fullfile(gcfg.resultspath, [gcfg.subjectName '_ripple_surrogate_timelock.mat']),'ripple_surrogate_timelock');            
%         save('-v7.3',fullfile(gcfg.resultspath, [gcfg.subjectName '_ripple_surrogate_fft.mat']),'ripple_surrogate_fft');
% end

end % of function


