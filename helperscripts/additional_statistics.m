% additional semi-manual statistics

clear all;
load('GA_SO_[0.50-4.00Hz]_trough_HC_CFCtrialwise_[-0.50 1.00s]_[12.00 16.00Hz]_20140428T104041.mat');
SW2SP.ASIp = cell2mat(ASIp);
SW2SP.ASIp_ang = cell2mat(ASIp_ang);

load('GA_SO_[0.50-4.00Hz]_trough_HC_CFCtrialwise_[-0.50 1.00s]_[80.00 90.00Hz]_20140424T172935.mat');
SW2RP.ASIp = cell2mat(ASIp);
SW2RP.ASIp_ang = cell2mat(ASIp_ang);


diff_rad = [SW2SP.ASIp-SW2RP.ASIp]';
circ_r(diff_rad)


[pval,med,P] = circ_cmtest(SW2SP.ASIp, SW2RP.ASIp) % independent two-sample test

[h mu ul ll] = circ_mtest(diff_rad,0) % one-sample test against zero