% load files
filename = ['CMIP5_neu_9P_pmt.nc']; 
info = ncinfo(filename);
vars_info = {info.Variables.Name}

% read in variables
delta_i = ncread(filename,'delta_i');

model_ensemble = ["ACCESS1-0_r1i1p1", "ACCESS1-3_r1i1p1",... 
    "BNU-ESM_r1i1p1", "CCSM4_r1i1p1", "CESM1-BGC_r1i1p1",... 
    "CESM1-CAM5_r1i1p1", "CESM12-LE_r00i1p1", "CMCC-CESM_r1i1p1",... 
    "CMCC-CMS_r1i1p1", "CMCC-CM_r1i1p1", "CNRM-CM5_r1i1p1",... 
    "CSIRO-Mk3-6-0_r1i1p1", "CanESM2-LE_r01i1p1", "CanESM2_r1i1p1",... 
    "EC-EARTH_r1i1p1", "FGOALS-g2_r1i1p1", "FIO-ESM_r1i1p1",... 
    "GFDL-CM3_r1i1p1", "GFDL-ESM2G_r1i1p1", "GFDL-ESM2M_r1i1p1",... 
    "GISS-E2-H-CC_r1i1p1", "GISS-E2-H_r1i1p1", "GISS-E2-H_r1i1p2",... 
    "GISS-E2-H_r1i1p3", "GISS-E2-R-CC_r1i1p1", "GISS-E2-R_r1i1p1",... 
    "GISS-E2-R_r1i1p2", "GISS-E2-R_r1i1p3", "HadGEM2-AO_r1i1p1",... 
    "HadGEM2-CC_r1i1p1", "HadGEM2-ES_r1i1p1", "IPSL-CM5A-LR_r1i1p1",... 
    "IPSL-CM5A-MR_r1i1p1", "IPSL-CM5B-LR_r1i1p1",... 
    "MIROC-ESM-CHEM_r1i1p1", "MIROC-ESM_r1i1p1", "MIROC5_r1i1p1",... 
    "MPI-ESM-LR_r1i1p1", "MPI-ESM-MR_r1i1p1", "MPI-ESM_r001i1p3",... 
    "MRI-CGCM3_r1i1p1", "MRI-ESM1_r1i1p1", "NorESM1-ME_r1i1p1",... 
    "NorESM1-M_r1i1p1", "bcc-csm1-1-m_r1i1p1", "bcc-csm1-1_r1i1p1",... 
    "inmcm4_r1i1p1"]


filename = ['tas_neu_djf_perfect_model_ensemble.nc']; % s2
info = ncinfo(filename);
vars_info = {info.Variables.Name}

year = ncread(filename,'year');
tas = ncread(filename,'tas');
tas = tas';

member = ["CESM122-LE_r0i1p1", "MPI-GE_r1i1p3", "CanESM2-LE_r1i1p1",... 
    "ACCESS1-0_r1i1p1", "ACCESS1-3_r1i1p1", "BNU-ESM_r1i1p1", "CCSM4_r1i1p1",... 
    "CESM1-BGC_r1i1p1", "CESM1-CAM5_r1i1p1", "CMCC-CESM_r1i1p1",... 
    "CMCC-CMS_r1i1p1", "CMCC-CM_r1i1p1", "CNRM-CM5_r1i1p1",... 
    "CSIRO-Mk3-6-0_r1i1p1", "CanESM2_r1i1p1", "EC-EARTH_r1i1p1",... 
    "FGOALS-g2_r1i1p1", "FIO-ESM_r1i1p1", "GFDL-CM3_r1i1p1",... 
    "GFDL-ESM2G_r1i1p1", "GFDL-ESM2M_r1i1p1", "GISS-E2-H-CC_r1i1p1",... 
    "GISS-E2-H_r1i1p1", "GISS-E2-H_r1i1p2", "GISS-E2-H_r1i1p3",... 
    "GISS-E2-R-CC_r1i1p1", "GISS-E2-R_r1i1p1", "GISS-E2-R_r1i1p2",... 
    "GISS-E2-R_r1i1p3", "HadGEM2-AO_r1i1p1", "HadGEM2-CC_r1i1p1",... 
    "HadGEM2-ES_r1i1p1", "IPSL-CM5A-LR_r1i1p1", "IPSL-CM5A-MR_r1i1p1",... 
    "IPSL-CM5B-LR_r1i1p1", "MIROC-ESM-CHEM_r1i1p1", "MIROC-ESM_r1i1p1",... 
    "MIROC5_r1i1p1", "MPI-ESM-LR_r1i1p1", "MPI-ESM-MR_r1i1p1",... 
    "MRI-CGCM3_r1i1p1", "MRI-ESM1_r1i1p1", "NorESM1-ME_r1i1p1",... 
    "NorESM1-M_r1i1p1", "bcc-csm1-1-m_r1i1p1", "bcc-csm1-1_r1i1p1",... 
    "inmcm4_r1i1p1"]

member_reorder = [member(4:9),member(1),member(10:14),...
    member(3),member(15:40),member(2),member(41:47)];

tas_fix = [tas(:,4:9),tas(:,1),tas(:,10:14),...
    tas(:,3),tas(:,15:40),tas(:,2),tas(:,41:47)];

% CMIP5 subset for perfect model ensemble

% ind_sub = [2,3,6,10,11,12,14,15,16,17,18,22,26,31,32,37,38,41,43,46,47];

% delta_i = delta_i(ind_sub,ind_sub)
% tas_fix = tas_fix(:,ind_sub)
% member_reorder = member_reorder(ind_sub)

% compute target

aa = find(year >= 1990 & year <= 2009);
bb = find(year >= 2080 & year <= 2099);

tas_diff = nan(47,1);
for ii = 1:47
    tas_diff(ii) = mean(tas_fix(bb,ii)) - mean(tas_fix(aa,ii));
end



% compute target

dis_all = nan(47,201);
weighted_all = nan(47,201);

for jj = 1:47
    obs_ind = jj;
    
    if obs_ind == 1
        pred_ind = [jj+1:47];
    elseif obs_ind == 2
        pred_ind = [1,jj+1:47];
    elseif obs_ind == 46
        pred_ind = [1:jj-1,47];
    elseif obs_ind == 47
        pred_ind = [1:jj-1];
    else
        pred_ind = [1:jj-1,jj+1:47];
    end
   
    tas_o = tas_diff(obs_ind);
    tas_pred = tas_diff(pred_ind);

    sigma_q = [0:0.01:2];
 
    for ii = 1:length(sigma_q)
        perf_metric(:,ii) = exp(-(delta_i(pred_ind,jj)).^2/sigma_q(ii)^2);
        weighted_dist(ii) = sum(perf_metric(:,ii).*tas_pred)./sum(perf_metric(:,ii));
    end
    
    
    weighted_all(obs_ind,:) = weighted_dist;
    dis_all(obs_ind,:) = abs(tas_o-weighted_dist);
    [a,i] = min(abs(tas_o-weighted_dist));
    

    dis_best(obs_ind) = a;
    sigma_best(obs_ind) = sigma_q(i);
    
    
   
    figure,
    plot(sigma_q,weighted_dist)
    hold on
    plot([0 2],[tas_o tas_o])
    plot([sigma_q(i) sigma_q(i)],[0 10])
    title(member_reorder(jj))
end

tas_diff_rep = repmat(tas_diff,1,201);
for ii = 1:length(sigma_q)
    c_all = corrcoef(tas_diff_rep(:,ii),weighted_all(:,ii));
    corr_all(ii) = c_all(1,2);
end

[a,i] = max(corr_all);
figure
plot(sigma_q,corr_all)
hold on
plot([sigma_q(i) sigma_q(i)],[-1 1])


% choices for picking sigma_best (adjust here!)

aa = find(sigma_best == 2); 

for ii = 1:length(aa);
    dy = gradient(dis_all(aa(ii),:),0.01);
    thr = [2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2];
    bb = find(diff(dy)<= thr(ii) & diff(dy)>= 0);
%     figure,
%     plot(dis_all(aa(ii),:))
%     hold on
%     plot(diff(dy))
    bb_choice = [4,9,1,3,1,15,2,2];
%     plot([bb(bb_choice(ii)) bb(bb_choice(ii))],[-1 1])
    sigma_best(aa(ii)) = sigma_q(bb(bb_choice(ii)));
end

mean(sigma_best)
