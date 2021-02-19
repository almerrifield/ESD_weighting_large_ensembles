#######################################################

# ESD submission

#######################################################

** Source files updated after incorrect predictors found (BEST tasCLIM), use *_DA_fix_best.nc **
** Source files updated after incorrect shapefile found (DA.txt), use *_DA_fix.nc **

For basic versions of python figures generated (e.g. ESD_paper_Fig1.py)
Environment: iacpy3_2020

ESD_paper_Fig1.py
ESD_paper_Fig2.py
ESD_paper_Fig5.py

ESD_paper_FigC1_sat.py
ESD_paper_FigC1_slp.py
ESD_paper_FigC2_sat.py
ESD_paper_FigC2_slp.py

For matlab generated figures: Matlab 2019a

ESD_paper_Fig3c.m
ESD_paper_Fig3d.m

ESD_paper_Fig4.m
ESD_paper_FigB1.m
ESD_paper_FigB2.m


For weighting code generated figures: model weighting with patch
Environment: iacpy3_2019

cd model_weighting_paper_version/model_weighting/utils/
python ESD_fig3_neu_box_DA_fix_best.py -p /**path**/ CMIP5_neu_9P_DA_fix_best.nc CMIP5-ALL_neu_9P_DA_fix_best.nc CMIP5_neu_9P_DA_fix_best.nc CMIP5-ALL_neu_9P_DA_fix_best.nc CMIP5_neu_9P_DA_fix_best.nc CMIP5-ALL_neu_9P_DA_fix_best.nc CMIP5_neu_9P_DA_fix_best.nc CMIP5-ALL_neu_9P_DA_fix_best.nc CMIP5_neu_9P_DA_fix_best.nc CMIP5-ALL_neu_9P_DA_fix_best.nc -u -s box_NEU_9P_DA_fix_best.png
python ESD_figB3i.py -p /**path**/ CMIP5_neu_9P_005.nc CMIP5_neu_9P_01.nc CMIP5_neu_9P_02.nc  CMIP5_neu_9P_03.nc CMIP5_neu_9P_04.nc CMIP5_neu_9P_05.nc  CMIP5_neu_9P_06.nc CMIP5_neu_9P_07.nc CMIP5_neu_9P_08.nc -u -s figB3ai.png
python ESD_figB3ii.py -p /**path**/ CMIP5-ALL_neu_9P_005.nc CMIP5-ALL_neu_9P_01.nc CMIP5-ALL_neu_9P_02.nc  CMIP5-ALL_neu_9P_03.nc CMIP5-ALL_neu_9P_04.nc CMIP5-ALL_neu_9P_05.nc  CMIP5-ALL_neu_9P_06.nc CMIP5-ALL_neu_9P_07.nc CMIP5-ALL_neu_9P_08.nc -u -s figB3aii.png

python ESD_fig3_med_box_DA_fix_best.py -p /**path**/ CMIP5_med_9P_DA_fix_best.nc CMIP5-ALL_med_9P_DA_fix_best.nc CMIP5_med_9P_DA_fix_best.nc CMIP5-ALL_med_9P_DA_fix_best.nc CMIP5_med_9P_DA_fix_best.nc CMIP5-ALL_med_9P_DA_fix_best.nc CMIP5_med_9P_DA_fix_best.nc CMIP5-ALL_med_9P_DA_fix_best.nc CMIP5_med_9P_DA_fix_best.nc CMIP5-ALL_med_9P_DA_fix_best.nc -u -s box_MED_9P_DA_fix_best.png
python ESD_figB3i.py -p /**path**/ CMIP5_med_9P_005.nc CMIP5_med_9P_01.nc CMIP5_med_9P_02.nc  CMIP5_med_9P_03.nc CMIP5_med_9P_04.nc CMIP5_med_9P_05.nc  CMIP5_med_9P_06.nc CMIP5_med_9P_07.nc CMIP5_med_9P_08.nc -u -s figB3bi.png
python ESD_figB3ii.py -p /**path**/  CMIP5-ALL_med_9P_005.nc CMIP5-ALL_med_9P_01.nc CMIP5-ALL_med_9P_02.nc  CMIP5-ALL_med_9P_03.nc CMIP5-ALL_med_9P_04.nc CMIP5-ALL_med_9P_05.nc  CMIP5-ALL_med_9P_06.nc CMIP5-ALL_med_9P_07.nc CMIP5-ALL_med_9P_08.nc -u -s figB3bii.png
