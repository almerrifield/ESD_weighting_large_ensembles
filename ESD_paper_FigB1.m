addpath /Users/meranna/Desktop/linux_home/python/esd_weighting_large_ensembles/old_weightings

% load file
filename = ['CMIP5-ALL_neu_9P_DA_fix_best.nc'];
delta_q_neu = ncread(filename,'delta_q');

filename = ['CMIP5-ALL_med_9P_DA_fix_best.nc'];
delta_q_med = ncread(filename,'delta_q');


subplot(211),
bar(sort(delta_q_neu),'b')
hold on
[a,ii] =sort(delta_q_neu);

for jj = 14:63
    aa = find(ii == jj);
    bar(aa,delta_q_neu(jj),'r','EdgeColor','none')
end

for jj = 82:131
    aa = find(ii == jj);
    bar(aa,delta_q_neu(jj),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none')
end

for jj = 182:281
    aa = find(ii == jj);
    bar(aa,delta_q_neu(jj),'FaceColor',[0,0.5,0],'EdgeColor','none')
end

xticklabels('')
ylim([0, 1.4])
plot([0 288],[0.32 0.32],'w','LineWidth',5)
plot([0 288],[0.32 0.32],'k','LineWidth',3)
ylabel('Normalized RMSE distance from observations')
title('DJF NEU D_i (nine predictors)')

subplot(212),
bar(sort(delta_q_med),'b')
hold on
[a,ii] =sort(delta_q_med);

for jj = 14:63
    aa = find(ii == jj);
    bar(aa,delta_q_med(jj),'r','EdgeColor','none')
end

for jj = 82:131
    aa = find(ii == jj);
    bar(aa,delta_q_med(jj),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none')
end

for jj = 182:281
    aa = find(ii == jj);
    bar(aa,delta_q_med(jj),'FaceColor',[0,0.5,0],'EdgeColor','none')
end

xticklabels('')
ylim([0, 1.4])
plot([0 288],[0.4 0.4],'w','LineWidth',5)
plot([0 288],[0.4 0.4],'k','LineWidth',3)
ylabel('Normalized RMSE distance from observations')
title('JJA MED D_i (nine predictors)')

