declare -a sig=("Sig_MS40ct1" "Sig_MS40ct10" "Sig_MS40ct100" "Sig_MS40ct1000")

## Loop
for i in "${sig[@]}"
do
    echo -e "\n Starting $i"

    #Prepare output area
    rm -rf "model_$i"
    mkdir "model_$i"
    
    rm -rf param_ws.root
    root -l -b -q "build_ws.C+(\"$i\")"
    
    cp param_ws.root model_$i/.
    cp card_*.txt model_$i/.
    

done
