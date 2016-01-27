Control Plots:

python g1_exo_doFit_class.py --control -c el

python g1_exo_doFit_class.py --control -c mu

Simple Try Analysis:

python g1_exo_doFit_class.py  -b -c el

python g1_exo_doFit_class.py  -b -c mu



Do Analysis:
python runLimitsEXO750_mu.py --makeCards -b > xxx.log

mv cards_EXO_mu_HP750_g1/ cards_allCats

python runLimitsEXO750_mu.py --computeLimits

python runLimitsEXO750_mu.py --plotLimits  


