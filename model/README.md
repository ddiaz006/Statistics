# Instructions

## Setup
```
#build workspace
cmsenv
make
./MakeWorkspace --signal_model=Sig_MS15ct1
#signal_model can be any of the model names listed in data/signal_model_list.txt

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

# To run all limits do
```
source  scripts/make_all_signal_models.sh --input_list data/signal_model_list.txt --output_dir <your_directory>
```
