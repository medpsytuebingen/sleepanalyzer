addpath('/home/electromag/tilber/Projects/iEEGsleep_Bernhard/Analyses');
filename = '/home/electromag/tilber/Projects/iEEGsleep_Bernhard/Data/task_raw/rawdat.SIG';

[dInfoS, header] = loadSIGcopy_TOB(filename,1);
clear dInfoS;

chan = [1:header.chans];
dLen = header.tt*header.head.rf - 1; % 4032128
dStart = 1;
scalingFactor = 0.0975; % ??? ~0.0975 µV/bit (16bit amp, 65536 possible values, int16)

temp = loadSIGcopy_TOB(filename, chan, dLen, dStart);
data.trial = {temp'.* scalingFactor};
clear temp;

% data.hdr = ;
data.label = header.head.rnames';
data.time = {[0:1/header.head.rf:header.tt]};
data.fsample = header.head.rf;
data.cfg = [];
data.trialinfo = 1;

data.hdr.Fs = header.head.rf;
data.hdr.nChans = size(data.trial,1);
data.hdr.label = header.head.rnames';
data.hdr.nSamples = header.tt*header.head.rf - 1;
data.hdr.nSamplesPre = 0;
data.hdr.nTrials = 1;
data.hdr.orig = [];

% ispect data in fieldtrip
cfg = [];
ft_databrowser(cfg,data);
