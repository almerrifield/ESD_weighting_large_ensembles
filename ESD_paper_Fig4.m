%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4 - DJF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /Users/meranna/Desktop/linux_home/python/esd_weighting_large_ensembles/

% load file
filename = ['CMIP5-ALL_neu_9P_DA_fix_best.nc'];
weights_i_all = ncread(filename,'weights_i');

figure,
subplot(131),
bar(0.1, 2,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
hold on
bar(0.14, 1,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.18, 60,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.22, 3,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.26, 5,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.30, 10,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.34, 55,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.38, 5,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.42, 1,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.46, 3,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.50, 3,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.54, 12,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.58, 6,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.62, 6,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.66, 5,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.70, 104,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.74, 2,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.78, 2,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.82, 2,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.86, 1,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
xlim([0.075,0.885])
ylim([0,110])

plot([0.09 0.11],[weights_i_all(1) weights_i_all(1)],'b','linewidth',1) 
hold on
plot([0.09 0.11],[weights_i_all(2) weights_i_all(2)],'b','linewidth',1) 
%BNU
plot([0.13 0.15],[weights_i_all(3) weights_i_all(3)],'b','linewidth',1) 
% NCAR
for ii = 4:13
plot([0.17 0.19],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
for ii = 14:63
plot([0.17 0.19],[weights_i_all(ii) weights_i_all(ii)],'r','linewidth',1) 
end
%CMCC
for ii = 64:66
plot([0.21 0.23],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% CNRM
for ii = 67:71
plot([0.25 0.27],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% CSIRO
for ii = 72:81
plot([0.29 0.31],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% CanESM2
for ii = 82:131
plot([0.33 0.35],[weights_i_all(ii) weights_i_all(ii)],'Color',[0.93, 0.69, 0.13],'linewidth',1) 
end
for ii = 132:136
plot([0.33 0.35],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% EC-Earth
for ii = 137:141
plot([0.37 0.39],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% FGOALS
plot([0.41 0.43],[weights_i_all(142) weights_i_all(142)],'b','linewidth',1) 
% FIO
for ii = 143:145
plot([0.45 0.47],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% GFDL
for ii = 146:148
plot([0.49 0.51],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% GISS
for ii = 149:160
plot([0.53 0.55],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% Had
for ii = 161:166
plot([0.57 0.59],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% IPSL
for ii = 167:172
plot([0.61 0.63],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% MIROC
for ii = 173:177
plot([0.65 0.67],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% MPI
for ii = 178:181
plot([0.69 0.71],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
for ii = 182:281
plot([0.69 0.71],[weights_i_all(ii) weights_i_all(ii)],'Color',[0, 0.5, 0],'linewidth',1) 
end
% MRI
for ii = 282:283
plot([0.73 0.75],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
%NorESM
for ii = 284:285
plot([0.77 0.79],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% bcc
for ii = 286:287
plot([0.81 0.83],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
plot([0.85 0.87],[weights_i_all(288) weights_i_all(288)],'b','linewidth',1) 

ylabel('Independence scaling')
xticks([])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s
% Figure 4 - JJA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load file
filename = ['CMIP5-ALL_med_9P_DA_fix_best.nc'];
weights_i_all = ncread(filename,'weights_i');

subplot(132),
bar(0.1, 2,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
hold on
bar(0.14, 1,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.18, 60,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.22, 3,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.26, 5,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.30, 10,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.34, 55,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.38, 5,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.42, 1,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.46, 3,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.50, 3,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.54, 12,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.58, 6,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.62, 6,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.66, 5,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.70, 104,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.74, 2,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.78, 2,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.82, 2,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.86, 1,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
xlim([0.075,0.885])
ylim([0,110])

plot([0.09 0.11],[weights_i_all(1) weights_i_all(1)],'b','linewidth',1) 
hold on
plot([0.09 0.11],[weights_i_all(2) weights_i_all(2)],'b','linewidth',1) 
%BNU
plot([0.13 0.15],[weights_i_all(3) weights_i_all(3)],'b','linewidth',1) 
% NCAR
for ii = 4:13
plot([0.17 0.19],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
for ii = 14:63
plot([0.17 0.19],[weights_i_all(ii) weights_i_all(ii)],'r','linewidth',1) 
end
%CMCC
for ii = 64:66
plot([0.21 0.23],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% CNRM
for ii = 67:71
plot([0.25 0.27],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% CSIRO
for ii = 72:81
plot([0.29 0.31],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% CanESM2
for ii = 82:131
plot([0.33 0.35],[weights_i_all(ii) weights_i_all(ii)],'Color',[0.93, 0.69, 0.13],'linewidth',1) 
end
for ii = 132:136
plot([0.33 0.35],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% EC-Earth
for ii = 137:141
plot([0.37 0.39],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% FGOALS
plot([0.41 0.43],[weights_i_all(142) weights_i_all(142)],'b','linewidth',1) 
% FIO
for ii = 143:145
plot([0.45 0.47],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% GFDL
for ii = 146:148
plot([0.49 0.51],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% GISS
for ii = 149:160
plot([0.53 0.55],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% Had
for ii = 161:166
plot([0.57 0.59],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% IPSL
for ii = 167:172
plot([0.61 0.63],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% MIROC
for ii = 173:177
plot([0.65 0.67],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% MPI
for ii = 178:181
plot([0.69 0.71],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
for ii = 182:281
plot([0.69 0.71],[weights_i_all(ii) weights_i_all(ii)],'Color',[0, 0.5, 0],'linewidth',1) 
end
% MRI
for ii = 282:283
plot([0.73 0.75],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
%NorESM
for ii = 284:285
plot([0.77 0.79],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% bcc
for ii = 286:287
plot([0.81 0.83],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
plot([0.85 0.87],[weights_i_all(288) weights_i_all(288)],'b','linewidth',1) 

ylabel('Independence scaling')
xticks([])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4 - ANN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load file
filename = ['CMIP5-ALL_ANN_CLIM_v0.nc'];
weights_i_all = ncread(filename,'weights_i');

subplot(133),
bar(0.1, 2,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
hold on
bar(0.14, 1,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.18, 60,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.22, 3,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.26, 5,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.30, 10,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.34, 55,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.38, 5,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.42, 1,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.46, 3,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.50, 3,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.54, 12,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.58, 6,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.62, 6,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.66, 5,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.70, 104,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.74, 2,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.78, 2,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.82, 2,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
bar(0.86, 1,'BarWidth',0.02,'FaceColor',[0.80,0.80,0.80])
xlim([0.075,0.885])
ylim([0,110])

plot([0.09 0.11],[weights_i_all(1) weights_i_all(1)],'b','linewidth',1) 
hold on
plot([0.09 0.11],[weights_i_all(2) weights_i_all(2)],'b','linewidth',1) 
%BNU
plot([0.13 0.15],[weights_i_all(3) weights_i_all(3)],'b','linewidth',1) 
% NCAR
for ii = 4:13
plot([0.17 0.19],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
for ii = 14:63
plot([0.17 0.19],[weights_i_all(ii) weights_i_all(ii)],'r','linewidth',1) 
end
%CMCC
for ii = 64:66
plot([0.21 0.23],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% CNRM
for ii = 67:71
plot([0.25 0.27],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% CSIRO
for ii = 72:81
plot([0.29 0.31],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% CanESM2
for ii = 82:131
plot([0.33 0.35],[weights_i_all(ii) weights_i_all(ii)],'Color',[0.93, 0.69, 0.13],'linewidth',1) 
end
for ii = 132:136
plot([0.33 0.35],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% EC-Earth
for ii = 137:141
plot([0.37 0.39],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% FGOALS
plot([0.41 0.43],[weights_i_all(142) weights_i_all(142)],'b','linewidth',1) 
% FIO
for ii = 143:145
plot([0.45 0.47],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% GFDL
for ii = 146:148
plot([0.49 0.51],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% GISS
for ii = 149:160
plot([0.53 0.55],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% Had
for ii = 161:166
plot([0.57 0.59],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% IPSL
for ii = 167:172
plot([0.61 0.63],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% MIROC
for ii = 173:177
plot([0.65 0.67],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% MPI
for ii = 178:181
plot([0.69 0.71],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
for ii = 182:281
plot([0.69 0.71],[weights_i_all(ii) weights_i_all(ii)],'Color',[0, 0.5, 0],'linewidth',1) 
end
% MRI
for ii = 282:283
plot([0.73 0.75],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
%NorESM
for ii = 284:285
plot([0.77 0.79],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
% bcc
for ii = 286:287
plot([0.81 0.83],[weights_i_all(ii) weights_i_all(ii)],'b','linewidth',1) 
end
plot([0.85 0.87],[weights_i_all(288) weights_i_all(288)],'b','linewidth',1) 

ylabel('Independence scaling')
xticks([])