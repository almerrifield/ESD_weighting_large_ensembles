
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3 - DJF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% for paper version
% addpath /Users/meranna/Desktop/linux_home/python/esd_weighting_large_ensembles/old_weightings
% filename = ['CMIP5-ALL_neu_9P_DA_fix.nc'];

% load file
filename = ['CMIP5-ALL_neu_9P_DA_fix_best.nc'];

weights_all = ncread(filename,'weights'); % RMSE weighting
weights_q_all = ncread(filename,'weights_q'); % Performance weighting
weights_q_all_ic = ncread(filename,'weights_q'); % 1/N IC weighting
weights_q_all_mc = ncread(filename,'weights_q'); % 1/N model weighting

model_ensemble = ["ACCESS1-0_r1i1p1_LE", "ACCESS1-3_r1i1p1_LE",...
    "BNU-ESM_r1i1p1_LE", "CCSM4_r1i1p1_LE", "CCSM4_r2i1p1_LE",...
    "CCSM4_r3i1p1_LE", "CCSM4_r4i1p1_LE", "CCSM4_r5i1p1_LE",...
    "CCSM4_r6i1p1_LE", "CESM1-BGC_r1i1p1_LE", "CESM1-CAM5_r1i1p1_LE",...
    "CESM1-CAM5_r2i1p1_LE", "CESM1-CAM5_r3i1p1_LE", "CESM12-LE_r00i1p1_LE",...
    "CESM12-LE_r01i1p1_LE", "CESM12-LE_r02i1p1_LE", "CESM12-LE_r03i1p1_LE",...
    "CESM12-LE_r04i1p1_LE", "CESM12-LE_r05i1p1_LE", "CESM12-LE_r06i1p1_LE",...
    "CESM12-LE_r07i1p1_LE", "CESM12-LE_r08i1p1_LE", "CESM12-LE_r09i1p1_LE",...
    "CESM12-LE_r10i1p1_LE", "CESM12-LE_r11i1p1_LE", "CESM12-LE_r12i1p1_LE",...
    "CESM12-LE_r13i1p1_LE", "CESM12-LE_r14i1p1_LE", "CESM12-LE_r15i1p1_LE",...
    "CESM12-LE_r16i1p1_LE", "CESM12-LE_r17i1p1_LE", "CESM12-LE_r18i1p1_LE",...
    "CESM12-LE_r19i1p1_LE", "CESM12-LE_r20i1p1_LE", "CESM12-LE_r21i1p1_LE",...
    "CESM12-LE_r22i1p1_LE", "CESM12-LE_r23i1p1_LE", "CESM12-LE_r24i1p1_LE",...
    "CESM12-LE_r25i1p1_LE", "CESM12-LE_r26i1p1_LE", "CESM12-LE_r27i1p1_LE",...
    "CESM12-LE_r28i1p1_LE", "CESM12-LE_r29i1p1_LE", "CESM12-LE_r30i1p1_LE",...
    "CESM12-LE_r31i1p1_LE", "CESM12-LE_r32i1p1_LE", "CESM12-LE_r33i1p1_LE",...
    "CESM12-LE_r34i1p1_LE", "CESM12-LE_r35i1p1_LE", "CESM12-LE_r36i1p1_LE",...
    "CESM12-LE_r37i1p1_LE", "CESM12-LE_r38i1p1_LE", "CESM12-LE_r39i1p1_LE",...
    "CESM12-LE_r40i1p1_LE", "CESM12-LE_r41i1p1_LE", "CESM12-LE_r42i1p1_LE",...
    "CESM12-LE_r43i1p1_LE", "CESM12-LE_r44i1p1_LE", "CESM12-LE_r45i1p1_LE",...
    "CESM12-LE_r46i1p1_LE", "CESM12-LE_r47i1p1_LE", "CESM12-LE_r48i1p1_LE",...
    "CESM12-LE_r49i1p1_LE", "CMCC-CESM_r1i1p1_LE", "CMCC-CMS_r1i1p1_LE",...
    "CMCC-CM_r1i1p1_LE", "CNRM-CM5_r10i1p1_LE", "CNRM-CM5_r1i1p1_LE",...
    "CNRM-CM5_r2i1p1_LE", "CNRM-CM5_r4i1p1_LE", "CNRM-CM5_r6i1p1_LE",...
    "CSIRO-Mk3-6-0_r10i1p1_LE", "CSIRO-Mk3-6-0_r1i1p1_LE",...
    "CSIRO-Mk3-6-0_r2i1p1_LE", "CSIRO-Mk3-6-0_r3i1p1_LE",...
    "CSIRO-Mk3-6-0_r4i1p1_LE", "CSIRO-Mk3-6-0_r5i1p1_LE",...
    "CSIRO-Mk3-6-0_r6i1p1_LE", "CSIRO-Mk3-6-0_r7i1p1_LE",...
    "CSIRO-Mk3-6-0_r8i1p1_LE", "CSIRO-Mk3-6-0_r9i1p1_LE",...
    "CanESM2-LE_r01i1p1_LE", "CanESM2-LE_r02i1p1_LE",...
    "CanESM2-LE_r03i1p1_LE", "CanESM2-LE_r04i1p1_LE",...
    "CanESM2-LE_r05i1p1_LE", "CanESM2-LE_r06i1p1_LE",...
    "CanESM2-LE_r07i1p1_LE", "CanESM2-LE_r08i1p1_LE",...
    "CanESM2-LE_r09i1p1_LE", "CanESM2-LE_r10i1p1_LE",...
    "CanESM2-LE_r11i1p1_LE", "CanESM2-LE_r12i1p1_LE",...
    "CanESM2-LE_r13i1p1_LE", "CanESM2-LE_r14i1p1_LE",...
    "CanESM2-LE_r15i1p1_LE", "CanESM2-LE_r16i1p1_LE",...
    "CanESM2-LE_r17i1p1_LE", "CanESM2-LE_r18i1p1_LE",...
    "CanESM2-LE_r19i1p1_LE", "CanESM2-LE_r20i1p1_LE",...
    "CanESM2-LE_r21i1p1_LE", "CanESM2-LE_r22i1p1_LE",...
    "CanESM2-LE_r23i1p1_LE", "CanESM2-LE_r24i1p1_LE",...
    "CanESM2-LE_r25i1p1_LE", "CanESM2-LE_r26i1p1_LE",...
    "CanESM2-LE_r27i1p1_LE", "CanESM2-LE_r28i1p1_LE",...
    "CanESM2-LE_r29i1p1_LE", "CanESM2-LE_r30i1p1_LE",...
    "CanESM2-LE_r31i1p1_LE", "CanESM2-LE_r32i1p1_LE",...
    "CanESM2-LE_r33i1p1_LE", "CanESM2-LE_r34i1p1_LE",...
    "CanESM2-LE_r35i1p1_LE", "CanESM2-LE_r36i1p1_LE",...
    "CanESM2-LE_r37i1p1_LE", "CanESM2-LE_r38i1p1_LE",...
    "CanESM2-LE_r39i1p1_LE", "CanESM2-LE_r40i1p1_LE",...
    "CanESM2-LE_r41i1p1_LE", "CanESM2-LE_r42i1p1_LE",...
    "CanESM2-LE_r43i1p1_LE", "CanESM2-LE_r44i1p1_LE",...
    "CanESM2-LE_r45i1p1_LE", "CanESM2-LE_r46i1p1_LE",...
    "CanESM2-LE_r47i1p1_LE", "CanESM2-LE_r48i1p1_LE",...
    "CanESM2-LE_r49i1p1_LE", "CanESM2-LE_r50i1p1_LE", "CanESM2_r1i1p1_LE",...
    "CanESM2_r2i1p1_LE", "CanESM2_r3i1p1_LE", "CanESM2_r4i1p1_LE",...
    "CanESM2_r5i1p1_LE", "EC-EARTH_r12i1p1_LE", "EC-EARTH_r1i1p1_LE",...
    "EC-EARTH_r2i1p1_LE", "EC-EARTH_r8i1p1_LE", "EC-EARTH_r9i1p1_LE",...
    "FGOALS-g2_r1i1p1_LE", "FIO-ESM_r1i1p1_LE", "FIO-ESM_r2i1p1_LE",...
    "FIO-ESM_r3i1p1_LE", "GFDL-CM3_r1i1p1_LE", "GFDL-ESM2G_r1i1p1_LE",...
    "GFDL-ESM2M_r1i1p1_LE", "GISS-E2-H-CC_r1i1p1_LE", "GISS-E2-H_r1i1p1_LE",...
    "GISS-E2-H_r1i1p2_LE", "GISS-E2-H_r1i1p3_LE", "GISS-E2-H_r2i1p1_LE",...
    "GISS-E2-H_r2i1p3_LE", "GISS-E2-R-CC_r1i1p1_LE", "GISS-E2-R_r1i1p1_LE",...
    "GISS-E2-R_r1i1p2_LE", "GISS-E2-R_r1i1p3_LE", "GISS-E2-R_r2i1p1_LE",...
    "GISS-E2-R_r2i1p3_LE", "HadGEM2-AO_r1i1p1_LE", "HadGEM2-CC_r1i1p1_LE",...
    "HadGEM2-ES_r1i1p1_LE", "HadGEM2-ES_r2i1p1_LE", "HadGEM2-ES_r3i1p1_LE",...
    "HadGEM2-ES_r4i1p1_LE", "IPSL-CM5A-LR_r1i1p1_LE",...
    "IPSL-CM5A-LR_r2i1p1_LE", "IPSL-CM5A-LR_r3i1p1_LE",...
    "IPSL-CM5A-LR_r4i1p1_LE", "IPSL-CM5A-MR_r1i1p1_LE",...
    "IPSL-CM5B-LR_r1i1p1_LE", "MIROC-ESM-CHEM_r1i1p1_LE",...
    "MIROC-ESM_r1i1p1_LE", "MIROC5_r1i1p1_LE", "MIROC5_r2i1p1_LE",...
    "MIROC5_r3i1p1_LE", "MPI-ESM-LR_r1i1p1_LE", "MPI-ESM-LR_r2i1p1_LE",...
    "MPI-ESM-LR_r3i1p1_LE", "MPI-ESM-MR_r1i1p1_LE", "MPI-ESM_r001i1p3_LE",...
    "MPI-ESM_r002i1p3_LE", "MPI-ESM_r003i1p3_LE", "MPI-ESM_r004i1p3_LE",...
    "MPI-ESM_r005i1p3_LE", "MPI-ESM_r006i1p3_LE", "MPI-ESM_r007i1p3_LE",...
    "MPI-ESM_r008i1p3_LE", "MPI-ESM_r009i1p3_LE", "MPI-ESM_r010i1p3_LE",...
    "MPI-ESM_r011i1p3_LE", "MPI-ESM_r012i1p3_LE", "MPI-ESM_r013i1p3_LE",...
    "MPI-ESM_r014i1p3_LE", "MPI-ESM_r015i1p3_LE", "MPI-ESM_r016i1p3_LE",...
    "MPI-ESM_r017i1p3_LE", "MPI-ESM_r018i1p3_LE", "MPI-ESM_r019i1p3_LE",...
    "MPI-ESM_r020i1p3_LE", "MPI-ESM_r021i1p3_LE", "MPI-ESM_r022i1p3_LE",...
    "MPI-ESM_r023i1p3_LE", "MPI-ESM_r024i1p3_LE", "MPI-ESM_r025i1p3_LE",...
    "MPI-ESM_r026i1p3_LE", "MPI-ESM_r027i1p3_LE", "MPI-ESM_r028i1p3_LE",...
    "MPI-ESM_r029i1p3_LE", "MPI-ESM_r030i1p3_LE", "MPI-ESM_r031i1p3_LE",...
    "MPI-ESM_r032i1p3_LE", "MPI-ESM_r033i1p3_LE", "MPI-ESM_r034i1p3_LE",...
    "MPI-ESM_r035i1p3_LE", "MPI-ESM_r036i1p3_LE", "MPI-ESM_r037i1p3_LE",...
    "MPI-ESM_r038i1p3_LE", "MPI-ESM_r039i1p3_LE", "MPI-ESM_r040i1p3_LE",...
    "MPI-ESM_r041i1p3_LE", "MPI-ESM_r042i1p3_LE", "MPI-ESM_r043i1p3_LE",...
    "MPI-ESM_r044i1p3_LE", "MPI-ESM_r045i1p3_LE", "MPI-ESM_r046i1p3_LE",...
    "MPI-ESM_r047i1p3_LE", "MPI-ESM_r048i1p3_LE", "MPI-ESM_r049i1p3_LE",...
    "MPI-ESM_r050i1p3_LE", "MPI-ESM_r051i1p3_LE", "MPI-ESM_r052i1p3_LE",...
    "MPI-ESM_r053i1p3_LE", "MPI-ESM_r054i1p3_LE", "MPI-ESM_r055i1p3_LE",...
    "MPI-ESM_r056i1p3_LE", "MPI-ESM_r057i1p3_LE", "MPI-ESM_r058i1p3_LE",...
    "MPI-ESM_r059i1p3_LE", "MPI-ESM_r060i1p3_LE", "MPI-ESM_r061i1p3_LE",...
    "MPI-ESM_r062i1p3_LE", "MPI-ESM_r063i1p3_LE", "MPI-ESM_r064i1p3_LE",...
    "MPI-ESM_r065i1p3_LE", "MPI-ESM_r066i1p3_LE", "MPI-ESM_r067i1p3_LE",...
    "MPI-ESM_r068i1p3_LE", "MPI-ESM_r069i1p3_LE", "MPI-ESM_r070i1p3_LE",...
    "MPI-ESM_r071i1p3_LE", "MPI-ESM_r072i1p3_LE", "MPI-ESM_r073i1p3_LE",...
    "MPI-ESM_r074i1p3_LE", "MPI-ESM_r075i1p3_LE", "MPI-ESM_r076i1p3_LE",...
    "MPI-ESM_r077i1p3_LE", "MPI-ESM_r078i1p3_LE", "MPI-ESM_r079i1p3_LE",...
    "MPI-ESM_r080i1p3_LE", "MPI-ESM_r081i1p3_LE", "MPI-ESM_r082i1p3_LE",...
    "MPI-ESM_r083i1p3_LE", "MPI-ESM_r084i1p3_LE", "MPI-ESM_r085i1p3_LE",...
    "MPI-ESM_r086i1p3_LE", "MPI-ESM_r087i1p3_LE", "MPI-ESM_r088i1p3_LE",...
    "MPI-ESM_r089i1p3_LE", "MPI-ESM_r090i1p3_LE", "MPI-ESM_r091i1p3_LE",...
    "MPI-ESM_r092i1p3_LE", "MPI-ESM_r093i1p3_LE", "MPI-ESM_r094i1p3_LE",...
    "MPI-ESM_r095i1p3_LE", "MPI-ESM_r096i1p3_LE", "MPI-ESM_r097i1p3_LE",...
    "MPI-ESM_r098i1p3_LE", "MPI-ESM_r099i1p3_LE", "MPI-ESM_r100i1p3_LE",...
    "MRI-CGCM3_r1i1p1_LE", "MRI-ESM1_r1i1p1_LE", "NorESM1-ME_r1i1p1_LE",...
    "NorESM1-M_r1i1p1_LE", "bcc-csm1-1-m_r1i1p1_LE", "bcc-csm1-1_r1i1p1_LE",...
    "inmcm4_r1i1p1_LE"];

% Equal Weighting
weights_ones = ones(288,1);
weights_ones_norm = weights_ones./sum(weights_ones);

% Performance Weighting
weights_q_all_norm = weights_q_all./sum(weights_q_all);

% 1/N Weighting, IC members
weights_q_all_ic(4:9) = mean(weights_q_all_ic(4:9))./6; %CCSM4
weights_q_all_ic(11:13) = mean(weights_q_all_ic(11:13))./3; %CESM1
weights_q_all_ic(14:63) = mean(weights_q_all_ic(14:63))./50; %CESM12
weights_q_all_ic(67:71) = mean(weights_q_all_ic(67:71))./5; %CNRM
weights_q_all_ic(72:81) = mean(weights_q_all_ic(72:81))./10; %CSIRO
weights_q_all_ic(82:131) = mean(weights_q_all_ic(82:131))./50; %CanESM2LE
weights_q_all_ic(132:136) = mean(weights_q_all_ic(132:136))./5; %CanESM2
weights_q_all_ic(137:141) = mean(weights_q_all_ic(137:141))./5; %EC
weights_q_all_ic(143:145) = mean(weights_q_all_ic(143:145))./3; %FIO
weights_q_all_ic([150,153]) = mean(weights_q_all_ic([150,153]))./2; %GISS-E2-H
weights_q_all_ic([152,154]) = mean(weights_q_all_ic([152,154]))./2; %GISS-E2-H
weights_q_all_ic([156,159]) = mean(weights_q_all_ic([156,159]))./2; %GISS-E2-R
weights_q_all_ic([158,160]) = mean(weights_q_all_ic([158,160]))./2; %GISS-E2-R
weights_q_all_ic(163:166) = mean(weights_q_all_ic(163:166))./4; %Had
weights_q_all_ic(167:170) = mean(weights_q_all_ic(167:170))./4; %IPSL
weights_q_all_ic(175:177) = mean(weights_q_all_ic(175:177))./3; %MIROC
weights_q_all_ic(178:180) = mean(weights_q_all_ic(178:180))./3; %MPI
weights_q_all_ic(182:281) = mean(weights_q_all_ic(182:281))./100; %MPILE

weights_ic_all_norm = weights_q_all_ic./sum(weights_q_all_ic);

% 1/N Weighting, models
weights_q_all_mc(1:2) = mean(weights_q_all_mc(1:2))./2; %ACCESS
weights_q_all_mc(4:63) = mean(weights_q_all_mc(4:63))./60; %CESM
weights_q_all_mc(64:66) = mean(weights_q_all_mc(64:66))./3; %CMCC
weights_q_all_mc(67:71) = mean(weights_q_all_mc(67:71))./5; %CNRM
weights_q_all_mc(72:81) = mean(weights_q_all_mc(72:81))./10; %CSIRO
weights_q_all_mc(82:136) = mean(weights_q_all_mc(82:136))./55; %CanESM2LE
weights_q_all_mc(137:141) = mean(weights_q_all_mc(137:141))./5; %EC
weights_q_all_mc(143:145) = mean(weights_q_all_mc(143:145))./3; %FIO
weights_q_all_mc(146:148) = mean(weights_q_all_mc(146:148))./3; %GFDL
weights_q_all_mc(149:160) = mean(weights_q_all_mc(149:160))./12; %GISS
weights_q_all_mc(161:166) = mean(weights_q_all_mc(161:166))./6; %Had
weights_q_all_mc(167:172) = mean(weights_q_all_mc(167:172))./6; %IPSL
weights_q_all_mc(173:177) = mean(weights_q_all_mc(173:177))./5; %MIROC
weights_q_all_mc(178:281) = mean(weights_q_all_mc(178:281))./104; %MPI
weights_q_all_mc(282:283) = mean(weights_q_all_mc(282:283))./2; %MRI
weights_q_all_mc(284:285) = mean(weights_q_all_mc(284:285))./2; %NorESM1
weights_q_all_mc(286:287) = mean(weights_q_all_mc(286:287))./2; %bcc

weights_mc_all_norm = weights_q_all_mc./sum(weights_q_all_mc);

figure,
bar(1,sum(weights_ones_norm(14:63)),'r','EdgeColor','none')
hold on
bar(2,sum(weights_ones_norm(82:131)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none')
bar(3,sum(weights_ones_norm(182:281)),'FaceColor',[0,0.5,0],'EdgeColor','none')
bar(4,sum(weights_ones_norm([1:13,64:81,132:181,282:288])),'b','EdgeColor','none')

bar(6,sum(weights_q_all_norm(14:63)),'r','EdgeColor','none')
hold on
bar(7,sum(weights_q_all_norm(82:131)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none')
bar(8,sum(weights_q_all_norm(182:281)),'FaceColor',[0,0.5,0],'EdgeColor','none')
bar(9,sum(weights_q_all_norm([1:13,64:81,132:181,282:288])),'b','EdgeColor','none')

bar(11,sum(weights_ic_all_norm(14:63)),'r','EdgeColor','none')
hold on
bar(12,sum(weights_ic_all_norm(82:131)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none')
bar(13,sum(weights_ic_all_norm(182:281)),'FaceColor',[0,0.5,0],'EdgeColor','none')
bar(14,sum(weights_ic_all_norm([1:13,64:81,132:181,282:288])),'b','EdgeColor','none')

bar(16,sum(weights_mc_all_norm(14:63)),'r','EdgeColor','none')
bar(17,sum(weights_mc_all_norm(82:131)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none')
bar(18,sum(weights_mc_all_norm(182:281)),'FaceColor',[0,0.5,0],'EdgeColor','none')
bar(19,sum(weights_mc_all_norm([1:13,64:81,132:181,282:288])),'b','EdgeColor','none')

bar(21,sum(weights_all(14:63)),'r','EdgeColor','none')
bar(22,sum(weights_all(82:131)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none')
bar(23,sum(weights_all(182:281)),'FaceColor',[0,0.5,0],'EdgeColor','none')
bar(24,sum(weights_all([1:13,64:81,132:181,282:288])),'b','EdgeColor','none')

ylabel('Fraction of Total Weight')
xticks([2.5,7.5,12.5,17.5,22.5])
xticklabels({'equal','perf.','1/N ic','1/N mc','RMSE'})
ylim([0,1])
