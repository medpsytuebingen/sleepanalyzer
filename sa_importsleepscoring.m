%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sa_importsleepscoring
% by Til Ole Bergmann 2016
% last modified 2016/12/14 by TOB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% loads sleep scorings from SchlafAus Software
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scoring = sa_importsleepscoring(scoringFileName)

fileID = fopen(scoringFileName,'r'); % open file
scoring = [];
inputMat = [];
m = 0;
while (~feof(fileID)) % not end of file
    m = m + 1;
    inputText = textscan(fileID,'%f',2,'Delimiter','\n'); % read all 2 columns of row m at once
    scoring(m,:) = int8(cell2mat(inputText)');
end
fclose(fileID); % close file

end % of function