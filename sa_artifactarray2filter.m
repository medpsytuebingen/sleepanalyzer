%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sa_artifactarray2filter
% by Til Ole Bergmann 2016
% last modified 2016/12/15 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% converts artifacts from artifact array format with Nx1 cell array for N 
% channels, each cell containing an (Mx2) matrix with start and end points 
% of M artifact epochs in the 1st and 2nd column, respectively, 
% to the artifactFilter format with NXD matrix , with N channels and D 
% datapoints and 1 and 0 indicating artifact-free and artifact--containing 
% dataoints, respectively. 
% 
% cfg = [];
% cfg.prepad = duation of padding artifact to the left in seconds
% cfg.postpad = duation of padding artifact to the right in seconds
% data = sa_artifactarray2filter(cfg, data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = sa_artifactarray2filter(gcfg, data)
tic;

%% check and adjust config
if ~isfield(gcfg,'prepad'), gcfg.prepad = 0; end
if ~isfield(gcfg,'postpad'), gcfg.postpad = 0; end
if ~isfield(data,'artifact') 
    display('Data structure does not cointain .artifact field!'); 
    return
end

%% preparation
display(['Converting data.artifact to data.artifactFilter structure, padding with ' num2str(cfg.prepad) ' s pre- and ' num2str(cfg.postpad) ' s post-artifact...']);
prepad = gcfg.prepad * data.fsample;
postpad = gcfg.postpad * data.fsample;    

%% convert structures
data.artifactFilter = ones(size(data.trial{1}));
for i = 1:length(data.artifact)
    for j = 1:size(data.artifact{i})
        data.artifactFilter(i,max(1,(data.artifact{i}(j,1) - prepad)):min((data.artifact{i}(j,2) + postpad),size(data.trial{1},2))) = 0;
    end
end

%% finishing
ttoc = toc;
display('Sleep architecture inspected.');
display(['Converting artifact structure took ' num2str(ttoc) ' seconds.']);

end % of function