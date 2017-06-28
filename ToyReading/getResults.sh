#!/bin/bash
INPUTFILES=()

if [ "$2" == "toys" ]
then
    echo "Finding toy fits"
    INPUTDIR="/data/lhcb/users/nandia/B2DKstar/ToyStudies/"
elif [ "$2" == "data" ]
then
    echo "Finding data fits"
    INPUTDIR="/data/lhcb/users/nandia/B2DKstar/Systematics/"
else
    echo "Need to choose either data or toys"
    exit 0
fi

SUBDIR="2017-06-27_$1/"
MIN=1
MAX=1000
GAP=20

for (( i=$MIN; i <= $MAX; i+=$GAP ))
do
  INPUTFILES+=("${INPUTDIR}${SUBDIR}toySeed_$i/toy_$i.txt")
done

VARS=("A_d2kpi" 
      "A_d2kk" 
      "A_d2pipi" 
      "R_d2kk" 
      "R_d2pipi" 
      "Rplus_d2pik"
      "Rminus_d2pik"
      "A_d2kpipipi" 
      "A_d2pipipipi" 
      "R_d2pipipipi" 
      "Rplus_d2pikpipi"
      "Rminus_d2pikpipi"
      "bu_mean_kpi"
      "bu_mean_kpipipi"
      "bu_width_kpi"
      "bu_width_kpipipi"
      "exp_kpi_LL_combs_slope"
      "exp_kpi_DD_combs_slope"
      "exp_kpipipi_LL_combs_slope"
      "exp_kpipipi_DD_combs_slope"
#      "n_comb_d2kk_minus_DD_all"
#      "n_comb_d2kk_minus_LL_all"
#      "n_comb_d2kk_plus_DD_all"
#      "n_comb_d2kk_plus_LL_all"
#      "n_comb_d2kpi_minus_DD_all"
#      "n_comb_d2kpi_minus_LL_all"
#      "n_comb_d2kpi_plus_DD_all"
#      "n_comb_d2kpi_plus_LL_all"
#      "n_comb_d2pik_minus_DD_all"
#      "n_comb_d2pik_minus_LL_all"
#      "n_comb_d2pik_plus_DD_all"
#      "n_comb_d2pik_plus_LL_all"
#      "n_comb_d2pipi_minus_DD_all"
#      "n_comb_d2pipi_minus_LL_all"
#      "n_comb_d2pipi_plus_DD_all"
#      "n_comb_d2pipi_plus_LL_all"
       "adsSignificance"
#       "significanceMinus"
      )

rm results/*.txt

for VAR in ${VARS[*]}
do
  for FILE in ${INPUTFILES[*]}
  do
    grep ${VAR}" " ${FILE} >> results/${VAR}.txt
  done
  echo "%s/${VAR} //
  w
  q
  " | ex results/${VAR}.txt
  wc -l results/${VAR}.txt
done
