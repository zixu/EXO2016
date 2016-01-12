hadd -f mu_PKUTree_SingleTop_xww.root mu_out_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root        mu_out_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root mu_out_ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root  mu_out_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root mu_out_ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root

#ln -s mu_out_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root mu_PKUTree_TTBARpowheg_xww.root 
cp mu_out_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root mu_PKUTree_TTBARpowheg_xww.root 

hadd -f mu_PKUTree_VV_xww.root mu_out_WWToLNuQQ_13TeV-powheg.root mu_out_WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.root mu_out_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root

hadd -f  mu_PKUTree_WJetsPt180_xww.root  mu_out_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root

#hadd -f mu_PKUTree_allBkg_xww.root  mu_PKUTree_SingleTop_xww.root mu_PKUTree_TTBARpowheg_xww.root mu_PKUTree_VV_xww.root mu_PKUTree_WJetsPt180_xww.root

#ln -s mu_out_BulkGravWW600.root mu_PKUTree_BulkGravWW600.root
#ln -s mu_out_BulkGravWW700.root mu_PKUTree_BulkGravWW700.root
#ln -s mu_out_BulkGravWW750.root mu_PKUTree_BulkGravWW750.root
#ln -s mu_out_BulkGravWW800.root mu_PKUTree_BulkGravWW800.root
#ln -s mu_out_BulkGravWW900.root mu_PKUTree_BulkGravWW900.root
#ln -s mu_out_BulkGravWW1000.root mu_PKUTree_BulkGravWW1000.root
cp mu_out_BulkGravWW600.root mu_PKUTree_BulkGravWW600.root
cp mu_out_BulkGravWW700.root mu_PKUTree_BulkGravWW700.root
cp mu_out_BulkGravWW750.root mu_PKUTree_BulkGravWW750.root
cp mu_out_BulkGravWW800.root mu_PKUTree_BulkGravWW800.root
cp mu_out_BulkGravWW900.root mu_PKUTree_BulkGravWW900.root
cp mu_out_BulkGravWW1000.root mu_PKUTree_BulkGravWW1000.root

#root MakePseudoData.C\(\"mu\"\) -q
hadd -f mu_PKUTree_15D.root mu_out_singleMuon15DOct.root  mu_out_singleMuon15Dv4.root 


#hadd -f el_PKUTree_SingleTop_xww.root el_out_TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola.root el_out_TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola.root el_out_TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola.root el_out_TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola.root el_out_T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola.root el_out_Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola.root
#
#hadd -f  el_PKUTree_WJetsPt180_xww.root el_out_WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola.root el_out_WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola.root el_out_WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola.root el_out_WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola.root
#ln -s el_out_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola.root el_PKUTree_TTBARpowheg_xww.root
#ln -s el_out_LOWW-lvjj-PT180.root  el_PKUTree_VV_xww.root
#
#
#hadd -f el_PKUTree_allBkg_xww.root  el_PKUTree_SingleTop_xww.root el_PKUTree_TTBARpowheg_xww.root el_PKUTree_VV_xww.root el_PKUTree_WJetsPt180_xww.root
#
#ln -s el_out_treeEDBR_RSGravitonToWW_2000_bb.root  el_PKUTree_RSGravitonToWW_2000_xww.root
#
#root MakePseudoData.C\(\"el\"\) -q
