#Instructions
```
root -l extract_ratios.C
#and put ratios into datacards

#build workspace
root -l build_ws.C

#combine data cards
combineCards.py Name1=card_elemu.txt Name2=card_onepho.txt Name3=card_twomuzh.txt > card.txt

#create final workspace
text2workspace.py card.txt

#fit
combine -M MultiDimFit -v 3 card.root
```