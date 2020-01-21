
echo $1
##make clean and make
#make clean;make -j 16
make -j 16
##get workspace
./MakeWorkspace_VR --signal_model=Sig_MS55ct100
##signal_model can be any of the model names listed in data/signal_model_list.txt

##combine data cards
#combineCards.py EleMu=card_elemu.txt EleMuL=card_elemul.txt TwoMuDY=card_twomudy.txt TwoEleDY=card_twoeledy.txt TwoMuZH=card_twomuzh.txt TwoEleZH=card_twoelezh.txt > card.txt
#combineCards.py EleMu=card_elemu.txt EleMuL=card_elemul.txt TwoMuDY=card_twomudy.txt TwoMuZH=card_twomuzh.txt > card.txt
combineCards.py EleMu=card_elemu.txt EleMuL=card_elemul.txt TwoMuDY=card_twomudy.txt TwoMuZH=card_twomuzh_add_sys.txt > card.txt

#############################
####DEFAULT CONFOGURATION####
#############################
##create final workspace
#text2workspace.py card.txt
##get fitDiagnostic command from combine
#combine card.root -M FitDiagnostics --saveNormalizations --plot --saveShapes --saveWithUncertainties --saveWorkspace
#mv fitDiagnostics.root fitDiagnostics_$1.root

#############################
##MASK SEARCH CHANNELS#######
#############################
text2workspace.py card.txt --channel-masks
#combine card.root -M FitDiagnostics --saveNormalizations --plot --saveShapes --saveWithUncertainties --saveWorkspace --setParameters mask_TwoMuZH=1,mask_TwoEleZH=1#original
#combine card.root -M FitDiagnostics --saveNormalizations --plot --saveShapes --saveWithUncertainties --saveWorkspace --setParameters mask_TwoEleZH=1,mask_TwoEleDY=1
combine card.root -M FitDiagnostics --saveNormalizations --plot --saveShapes --saveWithUncertainties --saveWorkspace
mv fitDiagnostics.root fitDiagnostics_$1.root
