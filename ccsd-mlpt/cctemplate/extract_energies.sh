#!/bin/bash


echo "#POSCAR        DFT-PBE         HF              MP2          CCSD          CCSD(T) "

for posd in `ls cc_POSCAR_*`
do
pos=$(echo $posd | sed 's/cc_//')
if grep -q "free  en" OUTCAR.DFT.$pos &> /dev/null; then
    pbe_e=`awk <OUTCAR.DFT.$pos "/free  en/ { print \\$5 }"`
else
    pbe_e="N/A          "
fi

if grep -q "free  en" OUTCAR..HFT.$pos &> /dev/null; then
    hf_e=`awk <OUTCAR..HFT.$pos "/free  en/ { print \\$5 }"`

    if grep -q "converged" OUTCAR.MP2full.$pos &> /dev/null; then
        mp2_cor=`awk <OUTCAR.MP2full.$pos "/converged/ { print \\$5 }"`
	mp2=`python -c "e=( $mp2_cor + $hf_e ) ; print(e)" `
    else
       mp2_cor="N/A          "
    fi

    ccfile_5nos=`echo "cc4s.out.latt..enc.400.nbno.624."$pos".1Em5new"`
    if grep -q "corrected Energy:" $ccfile_5nos &> /dev/null; then
        cc_cor_5nos=`awk < $ccfile_5nos "/corrected Energy:/ { print \\$3 }"`
	ccsd_5nos=`python -c "e=( $cc_cor_5nos + $hf_e ) ; print(e)" `

        if grep -q "Atrip: Energy:" $ccfile_5nos &> /dev/null; then
            triples_cor_5nos=`awk < $ccfile_5nos "/Atrip: Energy:/ { print \\$3 }"`
	    triples_5nos=`python -c "e=( $ccsd_5nos + $triples_cor_5nos ) ; print(e)" `
        else
            triples_5nos="N/A          "
        fi
    else
        ccsd_5nos="N/A          "
        triples_5nos="N/A          "
    fi

    ccfile_10nos=`echo "cc4s.out.latt..enc.400.nbno.1144."$pos".1Em5new"`
    if grep -q "corrected Energy:" $ccfile_10nos &> /dev/null; then
        cc_cor_10nos=`awk < $ccfile_10nos "/corrected Energy:/ { print \\$3 }"`
	ccsd_10nos=`python -c "e=( $cc_cor_10nos + $hf_e ) ; print(e)" `


        if grep -q "Atrip: Energy:" $ccfile_5nos &> /dev/null; then
            triples_cor_5nos=`awk < $ccfile_5nos "/Atrip: Energy:/ { print \\$3 }"`
	    triples_5_10_nos=`python -c "e=( $ccsd_10nos + $triples_cor_5nos ) ; print(e)" `
        else
            triples_5_10_nos="N/A          "
        fi

    else
        ccsd_10nos="N/A          "
        triples_5_10_nos="N/A          "
    fi


else
    hf_e="N/A          "
    mp2="N/A          "
    ccsd_5nos="N/A          "
    ccsd_10nos="N/A          "
fi

echo "${pos}  ${pbe_e}  ${hf_e}  ${mp2}  ${ccsd_10nos}  ${triples_5_10_nos}"

done
