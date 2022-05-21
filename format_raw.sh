#!/bin/bash



WFROOT="/scratch/scw1019/DV/spN"
declare -A WF_FILES
#WF_FILES["WF_2_12_2.40"]="${WFROOT}/sp2/nc2_b2.40_12/WF/out_2_12_2.400000" 
#WF_FILES["WF_2_16_2.475"]="${WFROOT}/sp2/nc2_b2.475_16/WF/out_2_16_2.475000"
WF_FILES["WF_2_20_2.55"]="${WFROOT}/sp2/nc2_b2.55_20/WF/out_2_20_2.550000" 
WF_FILES["WF_2_24_2.60"]="${WFROOT}/sp2/nc2_b2.60_24/WF/out_2_24_2.600000" 
WF_FILES["WF_2_32_2.65"]="${WFROOT}/sp2/nc2_b2.65_32/WF/out_2_32_2.650000" 
WF_FILES["WF_2_32_2.70"]="${WFROOT}/sp2/nc2_b2.7_32/WF/out_2_32_2.700000" 
WF_FILES["WF_4_20_7.7"]="${WFROOT}/sp4/nc4_b7.7_20/WF/out_4_20_7.7" 
WF_FILES["WF_4_20_7.72"]="${WFROOT}/sp4/nc4_b7.72_20/WF/out_4_20_7.72" 
WF_FILES["WF_4_20_7.76"]="${WFROOT}/sp4/nc4_b7.76_20/WF/out_4_20_7.76" 
WF_FILES["WF_4_20_7.78"]="${WFROOT}/sp4/nc4_b7.78_20/WF/out_4_20_7.78" 
WF_FILES["WF_4_20_7.80"]="${WFROOT}/sp4/nc4_b7.80_20/WF/out_4_20_7.80" 
WF_FILES["WF_4_20_7.85"]="${WFROOT}/sp4/nc4_b7.85_20/WF/out_4_20_7.85" 
WF_FILES["WF_4_24_8.2"]="${WFROOT}/sp4/nc4_b8.2_24/WF/out_4_24_8.2" 
WF_FILES["WF_6_18_15.75"]="${WFROOT}/sp6/nc6_b15.75_18/WF/out_6_18_15.750000" 
WF_FILES["WF_6_16_15.9"]="${WFROOT}/sp6/nc6_b15.9_16/WF/out_6_16_15.9" 
WF_FILES["WF_6_16_16.1"]="${WFROOT}/sp6/nc6_b16.1_16/WF/out_6_16_16.1" 
WF_FILES["WF_6_20_16.3"]="${WFROOT}/sp6/nc6_b16.3_20/WF/out_6_20_16.3" 
WF_FILES["WF_8_16_26.5"]="${WFROOT}/sp8/nc8_b26.5_16/WF/out_8_16_26.5" 
WF_FILES["WF_8_16_26.7"]="${WFROOT}/sp8/nc8_b26.7_16/WF/out_8_16_26.7" 
WF_FILES["WF_8_16_27.0"]="${WFROOT}/sp8/nc8_b27.0_16/WF/out_8_16_27.0" 
WF_FILES["WF_8_16_27.2"]="${WFROOT}/sp8/nc8_b27.2_16/WF/out_8_16_27.2" 


TH=0
for f in "${!WF_FILES[@]}"; do
        if [ -f "${WF_FILES[$f]}" ]; then
                echo ${f}
                grep -c "read" ${WF_FILES[$f]}
                #grep "\[WILSONFLOW\]" ${WF_FILES[${f}]} | cut -d"=" -f2 > tmp
                #awk -v TH=${TH} '{if($1==0.0) i++; if( i>TH) print (i-TH),$0;}' tmp > ~/WF_spN/${f}
                grep "Plaquette" ${WF_FILES[${f}]} | cut -d"=" -f2 > ~/WF_spN/${f}_plaq
        else
                echo "file ${WF_FILES[$f]}" does not exist
        fi
done
