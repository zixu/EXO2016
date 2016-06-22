Control Plots:

python g1_exo_doFit_class.py --control -c el -b

python g1_exo_doFit_class.py --control -c mu -b

Simple Try Analysis:

python g1_exo_doFit_class.py  -b -c el

python g1_exo_doFit_class.py  -b -c mu



Do Full Analysis:

python runLimitsEXO_TeV_mu.py --channel mu --makeCards -b 

python runLimitsEXO_TeV_mu.py --channel el --makeCards -b 

mkdir cards_allCats

cp cards_EXO_*/* cards_allCats

python runLimitsEXO_TeV_mu.py --channel mu --computeLimits

python runLimitsEXO_TeV_mu.py --channel el --computeLimits

python runLimitsEXO_TeV_mu.py --channel mu --plotLimits  

python runLimitsEXO_TeV_mu.py --channel el --plotLimits  

