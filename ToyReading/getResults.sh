#!/bin/bash
INPUTFILES=()

INPUTDIR="/data/lhcb/users/nandia/B2DKstar/ToyStudies/"
SUBDIR="2016-08-01_currentStatus_genscale10/"
MIN=1
MAX=2000
GAP=40

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

VARS=("A_d2kpi" 
      "A_d2kk" 
      "A_d2pipi" 
      #"A_d2pik" 
      "R_d2kk" 
      "R_d2pipi" 
      #"R_d2pik"
      "Rminus_d2pik"
      "Rplus_d2pik"
      #"n_comb_d2pik_plus_DD_all"
      "bu_mean"
      "bu_width"
      "d2kpi_exp_DD_combs_slope"
      "d2kpi_exp_LL_combs_slope"
      "n_comb_d2kk_minus_DD_all"
      "n_comb_d2kk_minus_LL_all"
      "n_comb_d2kk_plus_DD_all"
      "n_comb_d2kk_plus_LL_all"
      "n_comb_d2kpi_minus_DD_all"
      "n_comb_d2kpi_minus_LL_all"
      "n_comb_d2kpi_plus_DD_all"
      "n_comb_d2kpi_plus_LL_all"
      "n_comb_d2pik_minus_DD_all"
      "n_comb_d2pik_minus_LL_all"
      "n_comb_d2pik_plus_DD_all"
      "n_comb_d2pik_plus_LL_all"
      "n_comb_d2pipi_minus_DD_all"
      "n_comb_d2pipi_minus_LL_all"
      "n_comb_d2pipi_plus_DD_all"
      "n_comb_d2pipi_plus_LL_all"
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
