t = GA_spindle_freq.time;
f = GA_spindle_freq.freq;
p = squeeze((GA_spindle_freq.powspctrm));

[H,Pval,CI,STATS] = ttest(p);
pcrit = .01;

d = STATS.tstat .* (Pval < pcrit);
% d = nanmean(p);

figure;
imagesc(t,f,squeeze(d))
axis xy
colorbar
ylim([50 200])
xlim([-2 2])
caxis([-max(abs(d(:))) max(abs(d(:)))])
%%