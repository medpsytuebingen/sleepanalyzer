%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sa_savefig
% by Til Ole Bergmann 2016
% last modified 2016/12/13 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% inspect sleep architecture in EEG data
%
% cfg = [];
% cfg.fhandle = handle of figuire to save
% cfg.path = full path to folder where figure shall be saved
% cfg.filename = string of filename (is also figure name) 
% cfg.pdf = save pdf in addition to matlab figure, 'yes' or 'no' (default = 'no');
% cfg.close = close figure after saveing 'yes' or 'no' (default = 'no');
% sa_savefig(cfg);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sa_savefig(gcfg)

%% check and adjust config
if ~isfield(gcfg,'pdf'), gcfg.pdf = 'no'; end
if ~isfield(gcfg,'close'), gcfg.close = 'no'; end

%% preparation
display(['Saving figure ' gcfg.filename '...']);
filename = fullfile(gcfg.path,gcfg.filename);

%% save .fig
set(gcfg.fhandle, 'Name',gcfg.filename,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
saveas(gcfg.fhandle,[filename '.fig']);

%% save .pdf
if strcmp(gcfg.pdf,'yes')
    set(gcfg.fhandle, 'PaperUnits', 'centimeters', 'PaperSize', [33.867 19.05], 'PaperPosition', [0 0 33 19]); % [-3 -1 30 19.5]
    print(gcfg.fhandle, '-dpdf', '-r300', [filename '.pdf']);
end
%% close figure
if strcmp(gcfg.close,'yes')
    close(gcfg.fhandle); 
end

end % of function
