for i = 1:1%length(All_SO_timelock)
    cfg = [];
    cfg.method = 'mtmconvol';
    cfg.taper = 'hanning';
    cfg.keeptrials = 'no';
    cfg.foi = [6:1:200]; 
    cfg.toi = [-2.5:0.001:2.5];
    cfg.t_ftimwin = 5./cfg.foi;
    cfg.output = 'pow';	
    cfg.polyremoval = 1; % 0 = mean (default), 1 = linear
    All_SO_evokedTFR{i} = ft_freqanalysis(cfg, All_SO_timelock{i});

    [a b c] = size(All_SO_evokedTFR{i}.powspctrm);
    bl = [-2.5 -1.5];
    blwin = All_SO_evokedTFR{i}.time >= bl(1) & All_SO_evokedTFR{i}.time <= bl(2);
    baseline = repmat(nanmean(All_SO_evokedTFR{i}.powspctrm(:,:,blwin,:),3),[1,1,c]);
    All_SO_evokedTFR{i}.powspctrm = (All_SO_evokedTFR{i}.powspctrm - baseline) ./ baseline .* 100;
        
    figure;
    sphandle = subplot(2,1,1);
    pos = get(sphandle,'position');
    pos(3) = pos(3)*0.85;  
    set(sphandle, 'position', pos);        
    plot(All_SO_timelock{i}.time,All_SO_timelock{i}.avg);
       
    subplot(2,1,2);
    cfg = [];
    cfg.zlim = 'maxmin';
    ft_singleplotTFR(cfg,All_SO_evokedTFR{i});
    title(['subject ' num2str(subjects(i))]);    
    
end

GA_SO_evokedTFR = ft_freqgrandaverage(cfg, All_SO_evokedTFR{:});

figure;
cfg = [];
ft_singleplotTFR(cfg,GA_SO_evokedTFR);
title(['GA']);

