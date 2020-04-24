#! /bin/bash

while [[ $# -gt 1 ]]
do
    #echo blash $#
    key="$1"
    #echo $key " " $2

    case $key in
	-il|--input_list)
	    input_list="$2"
	        shift #past argument
		    ;;
	-od|--output_dir)
	    output_dir="$2"
	        shift # past argument
		    ;;
	-i2|--inputDir2)
	    indir2="$2"
	        shift # past argument
		    ;;
	-o|--outputFile)
	    outfile="$2"
	        shift # past argument
		    ;;
	-d|--dataFile)
	    datafile="$2"
	        shift # past argument
		    ;;
	-smh|--smhFile)
	    smhfile="$2"
	        shift # past argument
		    ;;
	-s|--signalFile)
	    sfile="$2"
	        shift # past argument
		    ;;
	-b|--binningFile)
	    bfile="$2"
	        shift # past argument
		    ;;
	--default)
	    DEFAULT=YES
	        ;;
	*)
            # unknown option
	    ;;
    esac
    shift # past argument or value
done


echo INPUT LIST  = $input_list
echo OUTFILE DIR = $output_dir

make clean;make -j 12
mkdir -p $output_dir
while read model; do
  echo $model
  mkdir -p $output_dir/$model
  cp data/cards/exo_20_003/card_*.txt $output_dir/$model
  ./MakeWorkspace_VR_SYS --signal_model=$model
  mv param_ws.root $output_dir/$model
  cd $output_dir/$model
  #combineCards.py EleMu=card_elemu.txt EleMuL=card_elemul.txt TwoMuDY=card_twomudy.txt TwoEleDY=card_twoeledy.txt TwoMuZH=card_twomuzh.txt TwoEleZH=card_twoelezh.txt > card.txt
  combineCards.py EleMu=card_elemu.txt EleMuL=card_elemul.txt TwoMuDY=card_ll_dy.txt TwoMuZH=card_ll_zh.txt > card_systematics.txt #only for ee+mumu combined
  text2workspace.py card_systematics.txt
  #combine -M AsymptoticLimits card_systematics.root --run blind
  combine -M AsymptoticLimits card_systematics.root
  cd -
  
done <$input_list

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi
