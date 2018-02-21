% FFT_surrogate_playground
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all; close all;

% settings
Fs = 1000; % sampling rate
L = 3001; % length of signal
T = 1/Fs; % sample time
t = (0:L-1)*T; % time vector
NFFT = 2^nextpow2(L);

% generate data
offset = 0.5 * ones(1,L);
noise = randn(1,L);
alpha = sin(2*pi*10*t);
spindle_L = 0.5;
spindle_event = sin(2*pi*15*t(1:spindle_L*Fs+1)) .* hanning(spindle_L*Fs+1)';
spindle = [zeros(1,(L-(Fs*spindle_L+1))/2), spindle_event, zeros(1,(L-(Fs*spindle_L+1))/2)];
data = 0.85*spindle + 0.05*alpha + 0.1*noise + offset; % compose fake data

data = eventdata_old.trial{200};

% preprocess data
data = detrend(data,'linear'); % demean and detrend

% FFT
fourier = fft(data);
fabs = abs(fourier);
fphase = angle(fourier);

% phase scramble
randphase = rand(1,L).*(2*pi);
fphasenew = fphase + randphase;
fouriernew = fabs.*(cos(fphasenew)+1i*sin(fphasenew));
fouriernew([1 end]) = fourier([1 end]);

% recreate symmetry of fourier transformed data
for sample = 2:L/2
    fouriernew(sample) = real(fouriernew(L-sample+2)) - imag(fouriernew(L-sample+2))*1i;
end

% iFFT and check
datanew = real(ifft(fouriernew));
fouriercheck = fft(datanew);

% plot stuff
figure;
set(gcf,'WindowStyle','docked')

subplot(3,1,1);
title('data (blue = original, red = after phase-scrambling)');
cla;hold on;
plot(data,'b'); % original data
plot(datanew+1,'r'); % backtransformed data
hold off;

subplot(3,1,2);
title('fourier complex values  (blue = original, red = after phases-crambling)');
cla;hold on;
plot(fourier,'b'); % original data
plot(fouriercheck+1,'r'); % backtransformed data
hold off;

subplot(3,1,3);
title('power spectrum  (blue = original, red = after phase-scrambling)');
cla;hold on;
freqwin = [0 50]; % frequency window to plot
freq = Fs/2 * linspace(0,1,L/2+1);
powspctrm_orig = 2*abs(fourier(1:L/2+1)); % original data
plot(freq(freq>=freqwin(1) & freq<=freqwin(2)), powspctrm_orig(freq>=freqwin(1) & freq<=freqwin(2)),'b'); % original data
powspctrm_new = 2*abs(fouriercheck(1:L/2+1)); % backtransformed data
plot(freq(freq>=freqwin(1) & freq<=freqwin(2)), powspctrm_new(freq>=freqwin(1) & freq<=freqwin(2)),'r'); % original data
hold off;
