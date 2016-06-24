Control Plots:

python g1_exo_doFit_class.py --control -c el -b

python g1_exo_doFit_class.py --control -c mu -b




Do Full Analysis:

python runLimitsB2GWW.py --channel mu --makeCards -b 

python runLimitsB2GWW.py --channel el --makeCards -b 


python runLimitsB2GWW.py --channel mu --computeLimits

python runLimitsB2GWW.py --channel el --computeLimits

python runLimitsB2GWW.py --channel mu --plotLimits  

python runLimitsB2GWW.py --channel el --plotLimits  

