Control Plots:

python g1_exo_doFit_class.py --control -c el -b

python g1_exo_doFit_class.py --control -c mu -b




Do Full Analysis:

##for WprimeWZ

python runLimitsB2GWW.py --channel mu --makeCards -b --signalmodel WprimeWZ --batchMode --lxbatchCern &
 
python runLimitsB2GWW.py --channel el --makeCards -b --signalmodel WprimeWZ --batchMode --lxbatchCern &
 

rename datacards direcotory "cards_B2GWW_XXXX" to "cards_B2GWW", for example, ln -s cards_B2GWW_closuretest0_HP_ExpN cards_B2GWW

python runLimitsB2GWW.py --channel mu --computeLimits --signalmodel WprimeWZ
 
python runLimitsB2GWW.py --channel el --computeLimits --signalmodel WprimeWZ
 
python runLimitsB2GWW.py --channel em --computeLimits --signalmodel WprimeWZ

 
python runLimitsB2GWW.py --channel mu --plotLimits -b --signalmodel WprimeWZ
 
python runLimitsB2GWW.py --channel el --plotLimits -b --signalmodel WprimeWZ
 
python runLimitsB2GWW.py --channel em --plotLimits -b --signalmodel WprimeWZ
 

##for WprimeWZ-HVT-A

python runLimitsB2GWW.py --channel mu --makeCards -b --signalmodel WprimeWZ-HVT-A --batchMode --lxbatchCern &
 
python runLimitsB2GWW.py --channel el --makeCards -b --signalmodel WprimeWZ-HVT-A --batchMode --lxbatchCern &
 

rename datacards direcotory "cards_B2GWW_XXXX" to "cards_B2GWW", for example, ln -s cards_B2GWW_closuretest0_HP_ExpN cards_B2GWW

python runLimitsB2GWW.py --channel mu --computeLimits --signalmodel WprimeWZ-HVT-A
 
python runLimitsB2GWW.py --channel el --computeLimits --signalmodel WprimeWZ-HVT-A
 
python runLimitsB2GWW.py --channel em --computeLimits --signalmodel WprimeWZ-HVT-A

 
python runLimitsB2GWW.py --channel mu --plotLimits -b --signalmodel WprimeWZ-HVT-A
 
python runLimitsB2GWW.py --channel el --plotLimits -b --signalmodel WprimeWZ-HVT-A
 
python runLimitsB2GWW.py --channel em --plotLimits -b --signalmodel WprimeWZ-HVT-A
 



# for BulkGravWW


python runLimitsB2GWW.py --channel mu --makeCards -b --signalmodel BulkGravWW 

python runLimitsB2GWW.py --channel el --makeCards -b --signalmodel BulkGravWW  

#batchMode at cern lxplus
python runLimitsB2GWW.py --channel mu --makeCards -b --signalmodel BulkGravWW --batchMode --lxbatchCern &
python runLimitsB2GWW.py --channel el --makeCards -b --signalmodel BulkGravWW --batchMode --lxbatchCern &

rename datacards direcotory "cards_B2GWW_XXXX" to "cards_B2GWW", for example, ln -s cards_B2GWW_closuretest0_HP_ExpN cards_B2GWW

python runLimitsB2GWW.py --channel mu --computeLimits --signalmodel BulkGravWW  

python runLimitsB2GWW.py --channel el --computeLimits --signalmodel BulkGravWW  

python runLimitsB2GWW.py --channel em --computeLimits --signalmodel BulkGravWW  


python runLimitsB2GWW.py --channel mu --plotLimits -b --signalmodel BulkGravWW  

python runLimitsB2GWW.py --channel el --plotLimits -b --signalmodel BulkGravWW  

python runLimitsB2GWW.py --channel em --plotLimits -b --signalmodel BulkGravWW  

