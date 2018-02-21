clear all; close all; clc; 
dbstop if error;

% determine operating system 
if strcmp(filesep,'/')
    filesystem='linux';
    my_drive='/gpfs01/born/tbergmann';
    home_drive='/home';
elseif strcmp(filesep,'\')
    filesystem='windows';
    my_drive='X:';
else
    filesystem='unknown';
    display('Unknown filesystem!');
end

% add matlab paths
% filepath = (fullfile(my_drive,'Projects','iEEGsleep2_Bernhard','Data'));
% filepath = (fullfile(my_drive,'Projects','iEEGsleep2_Bernhard','Analyses','Results','S02'));
filepath = (fullfile(my_drive,'Projects','iEEGsleep2_Bernhard','Analyses','Results'));

%% GET DIRECTORY INFORMATION
% dirinfo = dir();
% dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
% dirinfo = dirinfo(3:end); % remove '.' and '..'
dirinfo(1).name = 'S02';
dirinfo(2).name = 'S03';
dirinfo(3).name = 'S05';
dirinfo(4).name = 'S06';
dirinfo(5).name = 'S07';
dirinfo(6).name = 'S09';
dirinfo(7).name = 'S10';
dirinfo(8).name = 'S11';
dirinfo(9).name = 'S12';
dirinfo(10).name = 'S13';
dirinfo(11).name = 'S17';
subdirinfo = cell(1,length(dirinfo));
for i = 1:length(dirinfo)
  thisdir = dirinfo(i).name;
%   subdirinfo{i} = dir(fullfile(thisdir, '*_3sec_travel_ar.mat'));
%   subdirinfo{i} = dir(fullfile(thisdir, [thisdir '_hightheta(16to22Hz)*.mat']));
  subdirinfo{i} = dir(fullfile(thisdir, '*hightheta(16to22Hz)*.mat'));
end

%% CHANGE FILENAMES
for i = 1:length(dirinfo)
    for j = 1:length(subdirinfo{i})
        oldName = subdirinfo{i}(j).name;
%         newName = [subdirinfo{i}(j).name(1:end-13) 'unequated_wake_nREM_REM_ar.mat'];
        newName = [oldName(1:4) 'beta' oldName(14:end)];
        display(['Renaming ' oldName '  to  ' newName]);        
        oldFullPath = fullfile(filepath, dirinfo(i).name, oldName);
        newFullPath = fullfile(filepath, dirinfo(i).name, newName);        
        movefile(oldFullPath, newFullPath);    
    end 
end







