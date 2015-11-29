#!/bin/bash

INPUTDIR='/home/cheungs/B02DKstar_Kshh/B02DKstar_Kshh_MassFit'
OUTPUTDIR='/home/cheungs/B02DKstar_Kshh/B02DKstar_Kshh_BinnedFit/Settings/PDFShapes'
RANGE='5200'
TAG='KsPiPiKsKK_withBdDstKst'

# Run script that converts range-dependent params to 5200
python convertParamsForBinnedFit.py ${RANGE} ${TAG}
sleep 1.0

# Copy FIXED params files
cd ${INPUTDIR}/Settings/PDFShapes/Gen/${TAG}/${RANGE}/
echo "Reading input from..."
echo "  "${INPUTDIR}/Settings/PDFShapes/Gen/${TAG}/${RANGE}/
echo "Sending output to..."
echo "  "${OUTPUTDIR}/Fit/
cp signalFIXED.txt ${OUTPUTDIR}/Fit/signalFIXED.txt
cp partrecoFIXED.txt ${OUTPUTDIR}/Fit/partrecoFIXED.txt
cp partrecoFIXED_corr.txt ${OUTPUTDIR}/Fit/partrecoFIXED_corr.txt
cp combsFIXED.txt ${OUTPUTDIR}/Fit/combsFIXED.txt
cp ${INPUTDIR}/Settings/PDFShapes/Fit/drho.txt ${OUTPUTDIR}/Fit/drho.txt

# Make an optional file that allows a different combinatoric slope for each bin
cp ${INPUTDIR}/Settings/PDFShapes/Fit/combs.txt ${OUTPUTDIR}/Fit/combs.txt
for i in {1..8}
do
  sed -n "s/mix_combs_slope/mix_binm${i}_combs_slope/gp" ${INPUTDIR}/Settings/PDFShapes/Fit/combs.txt >> ${OUTPUTDIR}/Fit/combs.txt
  sed -n "s/mix_combs_slope/mix_binp${i}_combs_slope/gp" ${INPUTDIR}/Settings/PDFShapes/Fit/combs.txt >> ${OUTPUTDIR}/Fit/combs.txt
done

# Make a copy for Gen/
echo "Making copy in..."
echo "  "${OUTPUTDIR}/Gen
cp ${OUTPUTDIR}/Fit/*.txt $OUTPUTDIR/Gen

##########################
echo "Copying more stuff to..."
echo "  "${OUTPUTDIR}/../../Inputs/${TAG}/

# Give ratios correct label
cp ratiosFIXED_${RANGE}_corr.txt ${OUTPUTDIR}/../../Inputs/${TAG}/ratiosFIXED_corr_${RANGE}.txt
cp ${INPUTDIR}/output/GenTotals_binnedfit_${RANGE}.txt ${OUTPUTDIR}/../../Inputs/${TAG}/GenTotals_${RANGE}.txt

#sed -n 's/ratio_bs_drho_both_mix_merge/Bd_Drho_ratio_to_Bs_5200/gp' ratiosFIXED.txt > ~/B02DKstar_Kshh/B02DKstar_Kshh_BinnedFit/Inputs/Ratio_Bs_to_All.txt
#sed -n 's/ratio_bs_dstkst_both_mix_merge/Bs_DstKst_ratio_to_Bs_5200/gp' ratiosFIXED.txt > ~/B02DKstar_Kshh/B02DKstar_Kshh_BinnedFit/Inputs/Ratio_Bs_to_BsDstKst.txt

