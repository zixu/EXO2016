python newstyle_B2GWW_doFit_class.py mu BulkGravWW750 650 850 40 150 600 1500 ExpN Pow -b -m 1  --category HP --closuretest 0  --realdata 1  --keepblind 0 -w 0 &
python newstyle_B2GWW_doFit_class.py el BulkGravWW750 650 850 40 150 600 1500 ExpN Pow -b -m 1  --category HP --closuretest 0  --realdata 1  --keepblind 0 -w 0 &
python newstyle_B2GWW_doFit_class.py mu BulkGravWW4500 4400 4600 40 150 800 5000 ExpN Pow -b -m 1 --category HP --closuretest 0  --realdata 1  --keepblind 0 -w 1  &
python newstyle_B2GWW_doFit_class.py el BulkGravWW4500 4400 4600 40 150 800 5000 ExpN Pow -b -m 1 --category HP --closuretest 0  --realdata 1  --keepblind 0 -w 1  &

python newstyle_B2GWW_doFit_class.py mu WprimeWZ-HVT-A800 700 900 40 150 600 1500 ExpN Pow -b -m 1  --category HP --closuretest 0  --realdata 1  --keepblind 0 -w 0 &
python newstyle_B2GWW_doFit_class.py el WprimeWZ-HVT-A800 700 900 40 150 600 1500 ExpN Pow -b -m 1  --category HP --closuretest 0  --realdata 1  --keepblind 0 -w 0 &
python newstyle_B2GWW_doFit_class.py mu WprimeWZ-HVT-A4500 4400 4600 40 150 800 5000 ExpN Pow -b -m 1 --category HP --closuretest 0  --realdata 1  --keepblind 0 -w 1  &
python newstyle_B2GWW_doFit_class.py el WprimeWZ-HVT-A4500 4400 4600 40 150 800 5000 ExpN Pow -b -m 1 --category HP --closuretest 0  --realdata 1  --keepblind 0 -w 1  &


#python newstyle_B2GWW_doFit_class.py el WprimeWZ-HVT-A4500 4400 4600 40 150 800 5000 ExpN Pow -b -m 1  --category HP --closuretest 0  --realdata 1  --keepblind 0 -w 1 &


#combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m 750 -n _lim_BulkGravWW750_mu_HP -d wwlvj_BulkGravWW750_mu_HP_unbin.txt -v 2 -S 1
#
#combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m 4500 -n _lim_BulkGravWW4500_mu_HP -d wwlvj_BulkGravWW4500_mu_HP_unbin.txt -v 2 -S 1
#
#
#combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m 750 -n _lim_BulkGravWW750_el_HP -d wwlvj_BulkGravWW750_el_HP_unbin.txt -v 2 -S 1 
#
#combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m 4500 -n _lim_BulkGravWW4500_el_HP -d wwlvj_BulkGravWW4500_el_HP_unbin.txt -v 2 -S 1
