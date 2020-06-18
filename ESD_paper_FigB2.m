addpath /Users/meranna/Desktop/linux_home/python/esd_weighting_large_ensembles/old_weightings
% load files
filename = ['CMIP5_neu_9P_DA_fix.nc']; % s2
sigma_i_neu = ncread(filename,'sigma_i');
delta_i_neu = ncread(filename,'delta_i');

% load files
filename = ['CMIP5-ALL_neu_9P_DA_fix.nc'];
sigma_i_all_neu = ncread(filename,'sigma_i');
delta_i_all_neu = ncread(filename,'delta_i');

for jj = 1:88;
delta_i_neu(jj:end,jj) = NaN;
end

for jj = 1:288;
delta_i_all_neu(jj:end,jj) = NaN;
end

delta_i_CESM12 = delta_i_all_neu(14:63,14:63);
delta_i_CanESM2LE = delta_i_all_neu(82:131,82:131);
delta_i_MPILE = delta_i_all_neu(182:281,182:281);

delta_i_ts = delta_i_neu(~isnan(delta_i_neu));
delta_i_m = sort(delta_i_ts);
len = length(delta_i_ts);

delta_i_25 = median(delta_i_m(1:len/2));
delta_i_75 = median(delta_i_m(len/2+1:end));
delta_i_min = min(delta_i_m);
delta_i_max = max(delta_i_m);

delta_i_CESM12_ts = delta_i_CESM12(~isnan(delta_i_CESM12));
delta_i_CESM12_m = sort(delta_i_CESM12_ts);
len = length(delta_i_CESM12_ts);

delta_i_CESM12_25 = delta_i_CESM12_m(306)+0.25*(delta_i_CESM12_m(307)-delta_i_CESM12_m(306));
delta_i_CESM12_75 = delta_i_CESM12_m(918)+0.75*(delta_i_CESM12_m(919)-delta_i_CESM12_m(918));
delta_i_CESM12_med = delta_i_CESM12_m(612)+0.5*(delta_i_CESM12_m(613)-delta_i_CESM12_m(612));
delta_i_CESM12_min = min(delta_i_CESM12_m);
delta_i_CESM12_max = max(delta_i_CESM12_m);

delta_i_CanESM2LE_ts = delta_i_CanESM2LE(~isnan(delta_i_CanESM2LE));
delta_i_CanESM2LE_m = sort(delta_i_CanESM2LE_ts);
len = length(delta_i_CanESM2LE_ts);

delta_i_CanESM2LE_25 = delta_i_CanESM2LE_m(306)+0.25*(delta_i_CanESM2LE_m(307)-delta_i_CanESM2LE_m(306));
delta_i_CanESM2LE_75 = delta_i_CanESM2LE_m(918)+0.75*(delta_i_CanESM2LE_m(919)-delta_i_CanESM2LE_m(918));
delta_i_CanESM2LE_med = delta_i_CanESM2LE_m(612)+0.5*(delta_i_CanESM2LE_m(613)-delta_i_CanESM2LE_m(612));
delta_i_CanESM2LE_min = min(delta_i_CanESM2LE_m);
delta_i_CanESM2LE_max = max(delta_i_CanESM2LE_m);

delta_i_MPILE_ts = delta_i_MPILE(~isnan(delta_i_MPILE));
delta_i_MPILE_m = sort(delta_i_MPILE_ts);
len = length(delta_i_MPILE_ts);

delta_i_MPILE_25 = median(delta_i_MPILE_m(1:len/2));
delta_i_MPILE_75 = median(delta_i_MPILE_m(len/2+1:end));
delta_i_MPILE_med = median(delta_i_MPILE_m);
delta_i_MPILE_min = min(delta_i_MPILE_m);
delta_i_MPILE_max = max(delta_i_MPILE_m);


subplot(311),
plot([delta_i_max delta_i_max],[0.99 1.01],'b','linewidth',3) % max whisker
hold on
plot([delta_i_min delta_i_min],[0.99 1.01],'b','linewidth',3) % min whisker
plot([delta_i_min delta_i_max], [1 1],'b','linewidth',3) % line between
plot([delta_i_75 delta_i_75],[0.99 1.01],'b','linewidth',3) % 75th whisker
plot([delta_i_25 delta_i_25],[0.99 1.01],'b','linewidth',3) % 25th whisker
plot([delta_i_25 delta_i_75], [0.99 0.99],'b','linewidth',3) % line between
plot([delta_i_25 delta_i_75], [1.01 1.01],'b','linewidth',3) % line between
plot([median(delta_i_m) median(delta_i_m)],[0.99 1.01],'b','linewidth',3)

plot([delta_i_CESM12_max delta_i_CESM12_max],[0.05+0.99 0.05+1.01],'r','linewidth',3) % max whisker
hold on
plot([delta_i_CESM12_min delta_i_CESM12_min],[0.05+0.99 0.05+1.01],'r','linewidth',3) % min whisker
plot([delta_i_CESM12_min delta_i_CESM12_max], [0.05+1 0.05+1],'r','linewidth',3) % line between
plot([delta_i_CESM12_75 delta_i_CESM12_75],[0.05+0.99 0.05+1.01],'r','linewidth',3) % 75th whisker
plot([delta_i_CESM12_25 delta_i_CESM12_25],[0.05+0.99 0.05+1.01],'r','linewidth',3) % 25th whisker
plot([delta_i_CESM12_25 delta_i_CESM12_75], [0.05+0.99 0.05+0.99],'r','linewidth',3) % line between
plot([delta_i_CESM12_25 delta_i_CESM12_75], [0.05+1.01 0.05+1.01],'r','linewidth',3) % line between
plot([delta_i_CESM12_med delta_i_CESM12_med],[0.05+0.99 0.05+1.01],'r','linewidth',3)

plot([delta_i_CanESM2LE_max delta_i_CanESM2LE_max],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % max whisker
hold on
plot([delta_i_CanESM2LE_min delta_i_CanESM2LE_min],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % min whisker
plot([delta_i_CanESM2LE_min delta_i_CanESM2LE_max], [0.1+1 0.1+1],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % line between
plot([delta_i_CanESM2LE_75 delta_i_CanESM2LE_75],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % 75th whisker
plot([delta_i_CanESM2LE_25 delta_i_CanESM2LE_25],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % 25th whisker
plot([delta_i_CanESM2LE_25 delta_i_CanESM2LE_75], [0.1+0.99 0.1+0.99],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % line between
plot([delta_i_CanESM2LE_25 delta_i_CanESM2LE_75], [0.1+1.01 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % line between
plot([delta_i_CanESM2LE_med delta_i_CanESM2LE_med],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3)

plot([delta_i_MPILE_max delta_i_MPILE_max],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % max whisker
hold on
plot([delta_i_MPILE_min delta_i_MPILE_min],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % min whisker
plot([delta_i_MPILE_min delta_i_MPILE_max], [0.15+1 0.15+1],'Color',[0 0.5 0],'linewidth',3) % line between
plot([delta_i_MPILE_75 delta_i_MPILE_75],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % 75th whisker
plot([delta_i_MPILE_25 delta_i_MPILE_25],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % 25th whisker
plot([delta_i_MPILE_25 delta_i_MPILE_75], [0.15+0.99 0.15+0.99],'Color',[0 0.5 0],'linewidth',3) % line between
plot([delta_i_MPILE_25 delta_i_MPILE_75], [0.15+1.01 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % line between
plot([delta_i_MPILE_med delta_i_MPILE_med],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3)

plot([sigma_i_neu sigma_i_neu],[0.95, 1.2],'k:','LineWidth',2)
xlim([0, 1.4])
ylim([0.95,1.2])
yticks([])
yticklabels('')
title('DJF NEU, nine predictors')



% load files
filename = ['CMIP5_med_9P_DA_fix.nc']; % s2
sigma_i_med = ncread(filename,'sigma_i');
delta_i_med = ncread(filename,'delta_i');

% load files
filename = ['CMIP5-ALL_med_9P_DA_fix.nc'];
sigma_i_all_med = ncread(filename,'sigma_i');
delta_i_all_med = ncread(filename,'delta_i');

for jj = 1:88;
delta_i_med(jj:end,jj) = NaN;
end

for jj = 1:288;
delta_i_all_med(jj:end,jj) = NaN;
end

delta_i_CESM12 = delta_i_all_med(14:63,14:63);
delta_i_CanESM2LE = delta_i_all_med(82:131,82:131);
delta_i_MPILE = delta_i_all_med(182:281,182:281);

delta_i_ts = delta_i_med(~isnan(delta_i_med));
delta_i_m = sort(delta_i_ts);
len = length(delta_i_ts);

delta_i_25 = median(delta_i_m(1:len/2));
delta_i_75 = median(delta_i_m(len/2+1:end));
delta_i_min = min(delta_i_m);
delta_i_max = max(delta_i_m);

delta_i_CESM12_ts = delta_i_CESM12(~isnan(delta_i_CESM12));
delta_i_CESM12_m = sort(delta_i_CESM12_ts);
len = length(delta_i_CESM12_ts);

delta_i_CESM12_25 = delta_i_CESM12_m(306)+0.25*(delta_i_CESM12_m(307)-delta_i_CESM12_m(306));
delta_i_CESM12_75 = delta_i_CESM12_m(918)+0.75*(delta_i_CESM12_m(919)-delta_i_CESM12_m(918));
delta_i_CESM12_med = delta_i_CESM12_m(612)+0.5*(delta_i_CESM12_m(613)-delta_i_CESM12_m(612));
delta_i_CESM12_min = min(delta_i_CESM12_m);
delta_i_CESM12_max = max(delta_i_CESM12_m);

delta_i_CanESM2LE_ts = delta_i_CanESM2LE(~isnan(delta_i_CanESM2LE));
delta_i_CanESM2LE_m = sort(delta_i_CanESM2LE_ts);
len = length(delta_i_CanESM2LE_ts);

delta_i_CanESM2LE_25 = delta_i_CanESM2LE_m(306)+0.25*(delta_i_CanESM2LE_m(307)-delta_i_CanESM2LE_m(306));
delta_i_CanESM2LE_75 = delta_i_CanESM2LE_m(918)+0.75*(delta_i_CanESM2LE_m(919)-delta_i_CanESM2LE_m(918));
delta_i_CanESM2LE_med = delta_i_CanESM2LE_m(612)+0.5*(delta_i_CanESM2LE_m(613)-delta_i_CanESM2LE_m(612));
delta_i_CanESM2LE_min = min(delta_i_CanESM2LE_m);
delta_i_CanESM2LE_max = max(delta_i_CanESM2LE_m);

delta_i_MPILE_ts = delta_i_MPILE(~isnan(delta_i_MPILE));
delta_i_MPILE_m = sort(delta_i_MPILE_ts);
len = length(delta_i_MPILE_ts);

delta_i_MPILE_25 = median(delta_i_MPILE_m(1:len/2));
delta_i_MPILE_75 = median(delta_i_MPILE_m(len/2+1:end));
delta_i_MPILE_med = median(delta_i_MPILE_m);
delta_i_MPILE_min = min(delta_i_MPILE_m);
delta_i_MPILE_max = max(delta_i_MPILE_m);


subplot(312),
plot([delta_i_max delta_i_max],[0.99 1.01],'b','linewidth',3) % max whisker
hold on
plot([delta_i_min delta_i_min],[0.99 1.01],'b','linewidth',3) % min whisker
plot([delta_i_min delta_i_max], [1 1],'b','linewidth',3) % line between
plot([delta_i_75 delta_i_75],[0.99 1.01],'b','linewidth',3) % 75th whisker
plot([delta_i_25 delta_i_25],[0.99 1.01],'b','linewidth',3) % 25th whisker
plot([delta_i_25 delta_i_75], [0.99 0.99],'b','linewidth',3) % line between
plot([delta_i_25 delta_i_75], [1.01 1.01],'b','linewidth',3) % line between
plot([median(delta_i_m) median(delta_i_m)],[0.99 1.01],'b','linewidth',3)

plot([delta_i_CESM12_max delta_i_CESM12_max],[0.05+0.99 0.05+1.01],'r','linewidth',3) % max whisker
hold on
plot([delta_i_CESM12_min delta_i_CESM12_min],[0.05+0.99 0.05+1.01],'r','linewidth',3) % min whisker
plot([delta_i_CESM12_min delta_i_CESM12_max], [0.05+1 0.05+1],'r','linewidth',3) % line between
plot([delta_i_CESM12_75 delta_i_CESM12_75],[0.05+0.99 0.05+1.01],'r','linewidth',3) % 75th whisker
plot([delta_i_CESM12_25 delta_i_CESM12_25],[0.05+0.99 0.05+1.01],'r','linewidth',3) % 25th whisker
plot([delta_i_CESM12_25 delta_i_CESM12_75], [0.05+0.99 0.05+0.99],'r','linewidth',3) % line between
plot([delta_i_CESM12_25 delta_i_CESM12_75], [0.05+1.01 0.05+1.01],'r','linewidth',3) % line between
plot([delta_i_CESM12_med delta_i_CESM12_med],[0.05+0.99 0.05+1.01],'r','linewidth',3)

plot([delta_i_CanESM2LE_max delta_i_CanESM2LE_max],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % max whisker
hold on
plot([delta_i_CanESM2LE_min delta_i_CanESM2LE_min],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % min whisker
plot([delta_i_CanESM2LE_min delta_i_CanESM2LE_max], [0.1+1 0.1+1],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % line between
plot([delta_i_CanESM2LE_75 delta_i_CanESM2LE_75],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % 75th whisker
plot([delta_i_CanESM2LE_25 delta_i_CanESM2LE_25],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % 25th whisker
plot([delta_i_CanESM2LE_25 delta_i_CanESM2LE_75], [0.1+0.99 0.1+0.99],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % line between
plot([delta_i_CanESM2LE_25 delta_i_CanESM2LE_75], [0.1+1.01 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % line between
plot([delta_i_CanESM2LE_med delta_i_CanESM2LE_med],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3)

plot([delta_i_MPILE_max delta_i_MPILE_max],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % max whisker
hold on
plot([delta_i_MPILE_min delta_i_MPILE_min],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % min whisker
plot([delta_i_MPILE_min delta_i_MPILE_max], [0.15+1 0.15+1],'Color',[0 0.5 0],'linewidth',3) % line between
plot([delta_i_MPILE_75 delta_i_MPILE_75],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % 75th whisker
plot([delta_i_MPILE_25 delta_i_MPILE_25],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % 25th whisker
plot([delta_i_MPILE_25 delta_i_MPILE_75], [0.15+0.99 0.15+0.99],'Color',[0 0.5 0],'linewidth',3) % line between
plot([delta_i_MPILE_25 delta_i_MPILE_75], [0.15+1.01 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % line between
plot([delta_i_MPILE_med delta_i_MPILE_med],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3)

plot([sigma_i_med sigma_i_med],[0.95, 1.2],'k:','LineWidth',2)
xlim([0, 1.4])
ylim([0.95,1.2])
yticks([])
yticklabels('')
title('JJA MED, nine predictors')

% load files
filename = ['CMIP5_ANN_CLIM_v0.nc']; % s2
sigma_i_med = ncread(filename,'sigma_i');
delta_i_med = ncread(filename,'delta_i');

% load files
filename = ['CMIP5-ALL_ANN_CLIM_v0.nc'];
sigma_i_all_med = ncread(filename,'sigma_i');
delta_i_all_med = ncread(filename,'delta_i');

for jj = 1:88;
delta_i_med(jj:end,jj) = NaN;
end

for jj = 1:288;
delta_i_all_med(jj:end,jj) = NaN;
end

delta_i_CESM12 = delta_i_all_med(14:63,14:63);
delta_i_CanESM2LE = delta_i_all_med(82:131,82:131);
delta_i_MPILE = delta_i_all_med(182:281,182:281);

delta_i_ts = delta_i_med(~isnan(delta_i_med));
delta_i_m = sort(delta_i_ts);
len = length(delta_i_ts);

delta_i_25 = median(delta_i_m(1:len/2));
delta_i_75 = median(delta_i_m(len/2+1:end));
delta_i_min = min(delta_i_m);
delta_i_max = max(delta_i_m);

delta_i_CESM12_ts = delta_i_CESM12(~isnan(delta_i_CESM12));
delta_i_CESM12_m = sort(delta_i_CESM12_ts);
len = length(delta_i_CESM12_ts);

delta_i_CESM12_25 = delta_i_CESM12_m(306)+0.25*(delta_i_CESM12_m(307)-delta_i_CESM12_m(306));
delta_i_CESM12_75 = delta_i_CESM12_m(918)+0.75*(delta_i_CESM12_m(919)-delta_i_CESM12_m(918));
delta_i_CESM12_med = delta_i_CESM12_m(612)+0.5*(delta_i_CESM12_m(613)-delta_i_CESM12_m(612));
delta_i_CESM12_min = min(delta_i_CESM12_m);
delta_i_CESM12_max = max(delta_i_CESM12_m);

delta_i_CanESM2LE_ts = delta_i_CanESM2LE(~isnan(delta_i_CanESM2LE));
delta_i_CanESM2LE_m = sort(delta_i_CanESM2LE_ts);
len = length(delta_i_CanESM2LE_ts);

delta_i_CanESM2LE_25 = delta_i_CanESM2LE_m(306)+0.25*(delta_i_CanESM2LE_m(307)-delta_i_CanESM2LE_m(306));
delta_i_CanESM2LE_75 = delta_i_CanESM2LE_m(918)+0.75*(delta_i_CanESM2LE_m(919)-delta_i_CanESM2LE_m(918));
delta_i_CanESM2LE_med = delta_i_CanESM2LE_m(612)+0.5*(delta_i_CanESM2LE_m(613)-delta_i_CanESM2LE_m(612));
delta_i_CanESM2LE_min = min(delta_i_CanESM2LE_m);
delta_i_CanESM2LE_max = max(delta_i_CanESM2LE_m);

delta_i_MPILE_ts = delta_i_MPILE(~isnan(delta_i_MPILE));
delta_i_MPILE_m = sort(delta_i_MPILE_ts);
len = length(delta_i_MPILE_ts);

delta_i_MPILE_25 = median(delta_i_MPILE_m(1:len/2));
delta_i_MPILE_75 = median(delta_i_MPILE_m(len/2+1:end));
delta_i_MPILE_med = median(delta_i_MPILE_m);
delta_i_MPILE_min = min(delta_i_MPILE_m);
delta_i_MPILE_max = max(delta_i_MPILE_m);


subplot(313),
plot([delta_i_max delta_i_max],[0.99 1.01],'b','linewidth',3) % max whisker
hold on
plot([delta_i_min delta_i_min],[0.99 1.01],'b','linewidth',3) % min whisker
plot([delta_i_min delta_i_max], [1 1],'b','linewidth',3) % line between
plot([delta_i_75 delta_i_75],[0.99 1.01],'b','linewidth',3) % 75th whisker
plot([delta_i_25 delta_i_25],[0.99 1.01],'b','linewidth',3) % 25th whisker
plot([delta_i_25 delta_i_75], [0.99 0.99],'b','linewidth',3) % line between
plot([delta_i_25 delta_i_75], [1.01 1.01],'b','linewidth',3) % line between
plot([median(delta_i_m) median(delta_i_m)],[0.99 1.01],'b','linewidth',3)

plot([delta_i_CESM12_max delta_i_CESM12_max],[0.05+0.99 0.05+1.01],'r','linewidth',3) % max whisker
hold on
plot([delta_i_CESM12_min delta_i_CESM12_min],[0.05+0.99 0.05+1.01],'r','linewidth',3) % min whisker
plot([delta_i_CESM12_min delta_i_CESM12_max], [0.05+1 0.05+1],'r','linewidth',3) % line between
plot([delta_i_CESM12_75 delta_i_CESM12_75],[0.05+0.99 0.05+1.01],'r','linewidth',3) % 75th whisker
plot([delta_i_CESM12_25 delta_i_CESM12_25],[0.05+0.99 0.05+1.01],'r','linewidth',3) % 25th whisker
plot([delta_i_CESM12_25 delta_i_CESM12_75], [0.05+0.99 0.05+0.99],'r','linewidth',3) % line between
plot([delta_i_CESM12_25 delta_i_CESM12_75], [0.05+1.01 0.05+1.01],'r','linewidth',3) % line between
plot([delta_i_CESM12_med delta_i_CESM12_med],[0.05+0.99 0.05+1.01],'r','linewidth',3)

plot([delta_i_CanESM2LE_max delta_i_CanESM2LE_max],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % max whisker
hold on
plot([delta_i_CanESM2LE_min delta_i_CanESM2LE_min],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % min whisker
plot([delta_i_CanESM2LE_min delta_i_CanESM2LE_max], [0.1+1 0.1+1],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % line between
plot([delta_i_CanESM2LE_75 delta_i_CanESM2LE_75],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % 75th whisker
plot([delta_i_CanESM2LE_25 delta_i_CanESM2LE_25],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % 25th whisker
plot([delta_i_CanESM2LE_25 delta_i_CanESM2LE_75], [0.1+0.99 0.1+0.99],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % line between
plot([delta_i_CanESM2LE_25 delta_i_CanESM2LE_75], [0.1+1.01 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3) % line between
plot([delta_i_CanESM2LE_med delta_i_CanESM2LE_med],[0.1+0.99 0.1+1.01],'Color',[0.9290 0.6940 0.1250],'linewidth',3)

plot([delta_i_MPILE_max delta_i_MPILE_max],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % max whisker
hold on
plot([delta_i_MPILE_min delta_i_MPILE_min],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % min whisker
plot([delta_i_MPILE_min delta_i_MPILE_max], [0.15+1 0.15+1],'Color',[0 0.5 0],'linewidth',3) % line between
plot([delta_i_MPILE_75 delta_i_MPILE_75],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % 75th whisker
plot([delta_i_MPILE_25 delta_i_MPILE_25],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % 25th whisker
plot([delta_i_MPILE_25 delta_i_MPILE_75], [0.15+0.99 0.15+0.99],'Color',[0 0.5 0],'linewidth',3) % line between
plot([delta_i_MPILE_25 delta_i_MPILE_75], [0.15+1.01 0.15+1.01],'Color',[0 0.5 0],'linewidth',3) % line between
plot([delta_i_MPILE_med delta_i_MPILE_med],[0.15+0.99 0.15+1.01],'Color',[0 0.5 0],'linewidth',3)

plot([sigma_i_all_med sigma_i_all_med],[0.95, 1.2],'k:','LineWidth',2)
xlim([0, 1.4])
ylim([0.95,1.2])
yticks([])
yticklabels('')
title('Annual large-scale CLIM predictors')

