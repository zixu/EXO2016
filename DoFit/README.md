Control Plots:

python g1_exo_doFit_class.py --control -c el -b

python g1_exo_doFit_class.py --control -c mu -b

Simple Try Analysis:

python g1_exo_doFit_class.py  -b -c el

python g1_exo_doFit_class.py  -b -c mu



Do Full Analysis:

python runLimitsEXO750_mu.py --channel mu --makeCards -b > mu.log

python runLimitsEXO750_mu.py --channel el --makeCards -b > el.log

mkdir cards_allCats
cp cards_EXO_*_HP750_g1/* cards_allCats

python runLimitsEXO750_mu.py --channel mu --computeLimits

python runLimitsEXO750_mu.py --channel el --computeLimits

python runLimitsEXO750_mu.py --channel mu --plotLimits  

python runLimitsEXO750_mu.py --channel el --plotLimits  

