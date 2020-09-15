
#  === ============== ============== 
#   1       16228066       16228163  
#   1       44951210       44951269  
#   1       85495021       85495094  
#   1      106171472      106171616  
#  === ============== ============== 
echo "- 1.3) Get the concensus peaks from all groups"
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
consensusPeaksBed="${analysisDir}/${projName}_$(basename ${analysisDir})_consensus_peaks.bed"
cat ${projDir}/results/bwa/mergedLibrary/macs/narrowPeak/consensus/{FKO/FKO.consensus_peaks,PKF1OE/PKF1OE.consensus_peaks}.bed | sort -k1,1V -k2,2g > ${consensusPeaksBed}.tmp
mergeBed -i ${consensusPeaksBed}.tmp > ${consensusPeaksBed}
rm ${consensusPeaksBed}.tmp

# Get the output consensus saf file
# ============ ===== ============== ============== ========= 
#   GeneID      Chr        Start         End         Strand  
# ============ ===== ============== ============== ========= 
# Interval_1     1       16228066       16228163      +     
# Interval_2     1       44951210       44951269      +     
# Interval_3     1      106171472      106171616      +     
# Interval_4     1      156615557      156615645      +     
# ============ ===== ============== ============== ========= 
consensusPeaksSaf="${analysisDir}/${projName}_$(basename ${analysisDir})_consensus_peaks.saf"
# Change float to integer in the python one liner
python -c "import sys; import pandas as pd; import datatable as dt; input_file=sys.argv[1]; outtxt_file=sys.argv[2]; 
peaksDT = dt.fread(input_file, sep='\t', header=False, nthreads=16);peaksDF = peaksDT.to_pandas(); 
peaksDF.rename(columns={ peaksDF.columns[0]: 'PeakChrom' }, inplace = True); peaksDF.rename(columns={ peaksDF.columns[1]: 'PeakStart' }, inplace = True); peaksDF.rename(columns={ peaksDF.columns[2]: 'PeakEnd'   }, inplace = True); peaksDF.insert (0, 'PeakID', ['Interval_{0}'.format(x) for x in range(peaksDF.shape[0])]); peaksDF.insert (4, 'PeakStrand', '+'); peaksDF.to_csv(outtxt_file, header=False, index=False, sep='\t', float_format='%.0f')" ${consensusPeaksBed} ${consensusPeaksSaf}

