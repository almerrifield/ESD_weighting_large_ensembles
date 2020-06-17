#######################################################

# ESD submission

#######################################################

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

ESD_paper_Fig3c.py
ESD_paper_Fig3d.py

ESD_paper_Fig4.py
ESD_paper_FigB1.py
ESD_paper_FigB2.py


For weighting code generated figures: model weighting with patch
Environment: iacpy3_2019

cd model_weighting_paper_version/model_weighting/utils/
python ESD_fig3_neu_box.py -p /**path**/ CMIP5_neu_9P_v2.nc CMIP5-ALL_neu_9P_v2.nc CMIP5_neu_9P_v2.nc CMIP5-ALL_neu_9P_v2.nc CMIP5_neu_9P_v2.nc CMIP5-ALL_neu_9P_v2.nc CMIP5_neu_9P_v2.nc CMIP5-ALL_neu_9P_v2.nc CMIP5_neu_9P_v2.nc CMIP5-ALL_neu_9P_v2.nc -u -s box_NEU_9P.png
python ESD_fig3_med_box.py -p /**path**/ CMIP5_med_9P_v0.nc CMIP5-ALL_med_9P_v0.nc CMIP5_med_9P_v0.nc CMIP5-ALL_med_9P_v0.nc CMIP5_med_9P_v0.nc CMIP5-ALL_med_9P_v0.nc CMIP5_med_9P_v0.nc CMIP5-ALL_med_9P_v0.nc CMIP5_med_9P_v0.nc CMIP5-ALL_med_9P_v0.nc -u -s box_MED_9P.png










