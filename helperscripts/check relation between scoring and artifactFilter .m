% load('/home/electromag/tilber/Projects/iEEGsleep_Bernhard/Data/CFC/S04/artifactFilter');
% datHC = load('/home/electromag/tilber/Projects/iEEGsleep_Bernhard/Data/CFC/S04/dataone_TRb05.mat');
% datCz = load('/home/electromag/tilber/Projects/iEEGsleep_Bernhard/Data/CFC/S04/dataone_Cz.mat');
% 
% load('/home/electromag/tilber/Projects/iEEGsleep_Bernhard/Data/CFC/S07/artifactFilter');
% datHC = load('/home/electromag/tilber/Projects/iEEGsleep_Bernhard/Data/CFC/S07/dataone_TL06.mat');
% datCz = load('/home/electromag/tilber/Projects/iEEGsleep_Bernhard/Data/CFC/S07/dataone_Cz.mat');
%
load('/home/electromag/tilber/Projects/iEEGsleep_Bernhard/Data/CFC/S17/artifactFilter');
datHC = load('/home/electromag/tilber/Projects/iEEGsleep_Bernhard/Data/CFC/S17/dataone_TL06.mat');
datCz = load('/home/electromag/tilber/Projects/iEEGsleep_Bernhard/Data/CFC/S17/dataone_Cz.mat');

display([num2str(sum(artifactFilter)/numel(artifactFilter))]);

figure;
hold on;
plot(datHC.data.scoring,'k');
plot(artifactFilter,'r');
temp = double(datHC.data.scoring);
temp(~artifactFilter) = NaN;
plot (temp,'Color','g','LineWidth', 5);
ylim([-0.2 5.2]);
hold off;

figure;
hist(temp,[0,1,2,3,4,5]);