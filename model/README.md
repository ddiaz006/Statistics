# Instructions

## Setup
```
root -l extract_ratios.C
#and put ratios into datacards

#build workspace
root -l 
.L build_ws.C++
build_ws()
.q

#combine data cards
combineCards.py EleMu=card_elemu.txt OnePho=card_onepho.txt TwoMuZH=card_twomuzh.txt > card.txt

#create final workspace
text2workspace.py card.txt
```

## Fits
```
combine -M FitDiagnostics -v 3 card.root -t -1
combine -M FitDiagnostics -v 3 card.root
combine -M MultiDimFit -v 3 card.root

#For pre and post fit plots
combine -M FitDiagnostics -v 3 card.root --saveShapes 

```