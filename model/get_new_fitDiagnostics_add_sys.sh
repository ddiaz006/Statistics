
echo $1
#make clean and make
make clean;make -j 16
#get workspace
./MakeWorkspace_VR --signal_model=Sig_MS55ct100
#signal_model can be any of the model names listed in data/signal_model_list.txt

#combine data cards
combineCards.py EleMu=data/cards/card_elemu.txt EleMuL=data/cards/card_elemul.txt TwoMuDY=data/cards/card_twomudy_add_sys.txt TwoEleDY=data/cards/card_twoeledy_add_sys.txt TwoMuZH=data/cards/card_twomuzh_add_sys.txt TwoEleZH=data/cards/card_twoelezh_add_sys.txt > card.txt
#create final workspace
cp param_ws.root data/cards/. 
text2workspace.py card.txt
#get fitDiagnostic command from combine
combine -M FitDiagnostics card.txt --saveNormalizations --plot --saveShapes --saveWithUncertainties --saveWorkspace
mv fitDiagnostics.root fitDiagnostics_$1.root
