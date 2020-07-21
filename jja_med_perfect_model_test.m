% load files
filename = ['CMIP5_med_9P_pmt.nc']; 
info = ncinfo(filename);
vars_info = {info.Variables.Name};

% read in variables
delta_i = ncread(filename,'delta_i');

model_ensemble = ["ACCESS1-0_r1i1p1_LE", "ACCESS1-3_r1i1p1_LE",... 
    "BNU-ESM_r1i1p1_LE", "CCSM4_r1i1p1_LE", "CESM1-BGC_r1i1p1_LE",... 
    "CESM1-CAM5_r1i1p1_LE", "CESM12-LE_r00i1p1_LE", "CMCC-CESM_r1i1p1_LE",... 
    "CMCC-CMS_r1i1p1_LE", "CMCC-CM_r1i1p1_LE", "CNRM-CM5_r1i1p1_LE",... 
    "CSIRO-Mk3-6-0_r1i1p1_LE", "CanESM2-LE_r01i1p1_LE", "CanESM2_r1i1p1_LE",... 
    "EC-EARTH_r1i1p1_LE", "FGOALS-g2_r1i1p1_LE", "FIO-ESM_r1i1p1_LE",... 
    "GFDL-CM3_r1i1p1_LE", "GFDL-ESM2G_r1i1p1_LE", "GFDL-ESM2M_r1i1p1_LE",... 
    "GISS-E2-H-CC_r1i1p1_LE", "GISS-E2-H_r1i1p1_LE", "GISS-E2-H_r1i1p2_LE",... 
    "GISS-E2-H_r1i1p3_LE", "GISS-E2-R-CC_r1i1p1_LE", "GISS-E2-R_r1i1p1_LE",... 
    "GISS-E2-R_r1i1p2_LE", "GISS-E2-R_r1i1p3_LE", "HadGEM2-AO_r1i1p1_LE",... 
    "HadGEM2-CC_r1i1p1_LE", "HadGEM2-ES_r1i1p1_LE", "IPSL-CM5A-LR_r1i1p1_LE",... 
    "IPSL-CM5A-MR_r1i1p1_LE", "IPSL-CM5B-LR_r1i1p1_LE",... 
    "MIROC-ESM-CHEM_r1i1p1_LE", "MIROC-ESM_r1i1p1_LE", "MIROC5_r1i1p1_LE",... 
    "MPI-ESM-LR_r1i1p1_LE", "MPI-ESM-MR_r1i1p1_LE", "MPI-ESM_r001i1p3_LE",... 
    "MRI-CGCM3_r1i1p1_LE", "MRI-ESM1_r1i1p1_LE", "NorESM1-ME_r1i1p1_LE",... 
    "NorESM1-M_r1i1p1_LE", "bcc-csm1-1-m_r1i1p1_LE", "bcc-csm1-1_r1i1p1_LE",... 
    "inmcm4_r1i1p1_LE"];


filename = ['tas_med_jja_perfect_model_ensemble.nc']; % s2
info = ncinfo(filename);
vars_info = {info.Variables.Name};

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
    "inmcm4_r1i1p1"];

member_reorder = [member(4:9),member(1),member(10:14),...
    member(3),member(15:40),member(2),member(41:47)];

tas_fix = [tas(:,4:9),tas(:,1),tas(:,10:14),...
    tas(:,3),tas(:,15:40),tas(:,2),tas(:,41:47)];

% compute target

aa = find(year >= 1990 & year <= 2009);
bb = find(year >= 2080 & year <= 2099);

tas_diff = nan(47,1);
for ii = 1:47
    tas_diff(ii) = mean(tas_fix(bb,ii)) - mean(tas_fix(aa,ii));
end

% compute target

dis_all = nan(47,201);

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
    
    dis_all(obs_ind,:) = abs(tas_o-weighted_dist);
    [a,i] = min(abs(tas_o-weighted_dist));

    dis_best(obs_ind) = a;
    sigma_best(obs_ind) = sigma_q(i);
    
   
    figure,
    plot(sigma_q,weighted_dist)
    hold on
    plot([0 2],[tas_o tas_o])
    plot([sigma_q(i) sigma_q(i)],[0 10])
    title(model_ensemble(jj))
end

% choices for picking sigma_best (adjust here)

aa = find(sigma_best <= 0.05); 
% 
for ii = 1:length(aa);
    lwbd = [6,12,21,23,8,3,11,21,15,11]
    figure,
    plot(dis_all(aa(ii),:))
    sigma_best(aa(ii)) = sigma_q(lwbd(ii));
end

mean(sigma_best)
