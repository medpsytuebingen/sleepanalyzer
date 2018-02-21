a = load('artifactFilter_Cz_12sec_equated_wake_nREM_REM_ar.mat');
b = load('artifactFilter_TL07_12sec_equated_wake_nREM_REM_ar.mat');
c = load('/home/electromag/tilber/Projects/iEEGsleep_Bernhard/Data/CFC/S02/artifactFilter_CzTL07_12sec_equated_wake_nREM_REM_ar.mat');
d.artifactFilter = all([a.artifactFilter;b.artifactFilter]);
e = load('artifactFilter_TL02_12sec_equated_wake_nREM_REM_ar.mat');
figure;
hold on;
% w = 1000000:2000000;
w = 3000000:4000000;
plot(c.artifactFilter(w), 'k');
plot(d.artifactFilter(w)-0.1, 'r');
plot(a.artifactFilter(w)-0.15, 'b');
plot(b.artifactFilter(w)-0.2, 'g')
plot(e.artifactFilter(w)-0.25, 'm');
ylim([-0.5 1.5]);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = load('artifactFilter_Cz_12sec_equated_wake_nREM_REM_ar.mat');
b = load('artifactFilter_TL06_12sec_equated_wake_nREM_REM_ar.mat');
c = load('/home/electromag/tilber/Projects/iEEGsleep_Bernhard/Data/CFC/S17/artifactFilter_CzTL06_12sec_equated_wake_nREM_REM_ar.mat');
d.artifactFilter = all([a.artifactFilter;b.artifactFilter]);
e = load('artifactFilter_TL01_12sec_equated_wake_nREM_REM_ar.mat');
figure;
hold on;
% w = 1000000:2000000;
w = 1000000:4000000;
plot(c.artifactFilter(w), 'k');
plot(d.artifactFilter(w)-0.1, 'r');
plot(a.artifactFilter(w)-0.15, 'b');
plot(b.artifactFilter(w)-0.2, 'g')
plot(e.artifactFilter(w)-0.25, 'm');
ylim([-0.5 1.5]);
hold off;
