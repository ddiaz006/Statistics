declare -a sig=("Sig_MS40ct1" "Sig_MS40ct10" "Sig_MS40ct100" "Sig_MS40ct1000")

## Loop
for i in "${sig[@]}"
do
    echo -e "\n Starting $i"

    #Prepare output area
    mydir="model_$i"
    rm -rf $mydir
    mkdir $mydir
    
    rm -rf param_ws.root
    root -l -b -q "build_ws.C+(\"$i\")"
    
    cp param_ws.root $mydir/.
    cp card_*.txt $mydir/.
    
    combineCards.py EleMu=$mydir/card_elemu.txt EleMuL=$mydir/card_elemul.txt TwoMuDY=$mydir/card_twomudy.txt TwoEleDY=$mydir/card_twoeledy.txt TwoMuZH=$mydir/card_twomuzh.txt TwoEleZH=$mydir/card_twoelezh.txt > $mydir/card.txt

    text2workspace.py $mydir/card.txt -o $mydir/card.root

    combine -M AsymptoticLimits $mydir/card.root -n ".$i" --run blind

done
