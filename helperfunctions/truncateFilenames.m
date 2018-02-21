function truncateFilenames(d, nchars)
%TRUNCATEFILENAMES truncates filenames
%   TRUNCATEFILENAMES(D, NCHARS) truncates the file names of all files in
%   directory D so that the filename (not including the extension) has
%   NCHARS characters. Shorter filenames are left unchanged. If the
%   resulting names will not be unique, an error is thrown and no changes
%   are made.
%   NCHARS must be a positive integer.
%
%   Caution! The operation is irreversible.
validateattributes(nchars, {'numeric'}, {'scalar' 'real' 'positive' 'integer'});
if ~exist(d, 'dir')
    error('matlabAnswers:truncateFilenames:noDir', ...
        'Directory not found');
end
fstruct = dir(d);    % examine directory
fnames = {fstruct.name};   % extract full names
% decompose names
[~, fnms, fexts] = cellfun(@fileparts, fnames, 'UniformOutput', false);
% retain only those that are too long initially
needTrunc = cellfun(@length, fnms) > nchars;
fnames = fnames(needTrunc);
fnms = fnms(needTrunc);
fexts = fexts(needTrunc);
% do the truncation, and replace the extension
fnmsout = cellfun(@(sn, se) [sn(1:nchars) se], fnms, fexts, 'UniformOutput', false);
% check for uniqueness
if ~isequal(length(fnmsout), length(unique(fnmsout)))
    error('matlabAnswers:truncateFilenames:nonunique', ...
        'Resulting filenames will not be unique');
end
% do the renames
for k = 1:length(fnms)
    movefile(fullfile(d, fnames{k}), fullfile(d, fnmsout{k}));
end
end