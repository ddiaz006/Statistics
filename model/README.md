# Instructions

## Setup
```
#build workspace
root -l 
.L build_ws.C++
build_ws()
.q

#combine data cards
combineCards.py EleMu=card_elemu.txt EleMuL=card_elemul.txt TwoMuDY=card_twomudy.txt TwoEleDY=card_twoeledy.txt TwoMuZH=card_twomuzh.txt TwoEleZH=card_twoelezh.txt > card.txt

#create final workspace
text2workspace.py card.txt
```

## Fits
```
combine -M FitDiagnostics -v 3 card.root
combine -M FitDiagnostics -v 3 card.root -t -1 --expectSignal 0
combine -M FitDiagnostics -v 3 card.root -t -1 --expectSignal 1
combine -M MultiDimFit -v 3 card.root

#For pre and post fit plots
combine -M FitDiagnostics -v 3 card.root --saveShapes 

#Limits
combine -M AsymptoticLimits card.root --run blind
```