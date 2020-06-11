% load file
filename = ['CMIP5-ALL_neu_9P_v2.nc'];
delta_q_neu = ncread(filename,'delta_q');

filename = ['CMIP5-ALL_med_9P_v0.nc'];
delta_q_med = ncread(filename,'delta_q');

figure,
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
ylabel('Normalized RMSE distance from observations')
title('JJA MED D_i (nine predictors)')

