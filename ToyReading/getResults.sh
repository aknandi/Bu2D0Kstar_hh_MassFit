#!/bin/bash

INPUTFILES=()

INPUTDIR="/data/lhcb/users/cheungs/B02DKstar_Kshh_toys/"
SUBDIR="2015-06-01_MassFit_modAlexis_relaxFrac010/"
MIN=1
MAX=10000
GAP=25

for (( i=$MIN; i <= $MAX; i+=$GAP ))
do
  INPUTFILES+=("${INPUTDIR}${SUBDIR}toySeed_$i/toy_$i.txt")
done

#INPUTDIR="/data/lhcb/users/cheungs/B02DKstar_Kshh_toys/"
#SUBDIR="2015-05-07_MassFit_noDstKst/"
#MIN=500
#MAX=990
#GAP=10
#
#for (( i=$MIN; i <= $MAX; i+=10 ))
#do
#  INPUTFILES+=("${INPUTDIR}${SUBDIR}toySeed_$i/toy_$i.txt")
#done
#
#INPUTFILES=("../TOYS/test2/toy_31.txt"
#            )

VARS=("bs_signal_mean" 
      "d2kspipi_exp_mix_combs_slope" 
      "frac010_bs_4900" 
      "n_bd_fit_d2kspipi_both_mix_merge" 
      "n_bs_fit_d2kspipi_both_mix_merge" 
      "n_comb_d2kspipi_both_mix_merge" 
      "ratio_bs_drho_both_mix_merge" 
      "ratio_bs_dstkst_both_mix_merge" 
      "ratio_bd_dstkst_both_mix_merge" 
      "ratio_bs_dstkst_010_both_mix_merge"
      "ratio_bs_dstkst_001_both_mix_merge"
      "signal_width" 
      "n_drho_both_mix_merge" 
      "n_bs_dstkst_both_mix_merge"
      )

rm results/*.txt

for VAR in ${VARS[*]}
do
  for FILE in ${INPUTFILES[*]}
  do
    grep ${VAR} ${FILE} >> results/${VAR}.txt
  done
  echo "%s/${VAR} //
  w
  q
  " | ex results/${VAR}.txt
  wc -l results/${VAR}.txt
done
