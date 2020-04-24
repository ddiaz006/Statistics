#run from "model" directory ./scripts/limit_test.sh
make -j 12
./MakeWorkspace_VR_SYS --signal_model=Sig_MS55ct100
cp param_ws.root data/cards/exo_20_003/.
combineCards.py EleMu=data/cards/exo_20_003/card_elemu.txt EleMuL=data/cards/exo_20_003/card_elemul.txt TwoMuDY=data/cards/exo_20_003/card_ll_dy.txt TwoMuZH=data/cards/exo_20_003/card_ll_zh.txt > card.txt

#no-systematics
#./MakeWorkspace_VR --signal_model=Sig_MS55ct30
#./MakeWorkspace_VR_SYS --signal_model=Sig_MS55ct30
#combineCards.py EleMu=card_elemu.txt EleMuL=card_elemul.txt TwoMuDY=card_twomudy.txt TwoMuZH=card_twomuzh.txt > card.txt

text2workspace.py card.txt
combine -M FitDiagnostics card.root --saveNormalizations --plot --saveShapes --saveWithUncertainties --saveWorkspace
combine -M AsymptoticLimits card.root 
