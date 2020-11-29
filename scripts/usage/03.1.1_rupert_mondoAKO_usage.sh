cd /home/rad/users/gaurav/projects/workflows/nfchipseq

# Parameters for the script
species="human"
user="rupert"
projName="mondoAKO"
outdir="/media/rad/HDD1/nfchip"
jobdir="/home/rad/users/gaurav/projects/workflows/nfchipseq"


#########################################################################################
# 1) RUN THE PIPELINE
#########################################################################################

# Parameters to run the pipeline
projDir="${outdir}/${user}/${projName}"
fastqDir="${projDir}/fastq"; mkdir -p ${fastqDir}

# # 1.1) Copy the fastq files
# cd /home/rad/media/nas/raw/TUM_NextSeq/201022_NB501802_0292_AHTV53BGXG/Data/Intensities/BaseCalls/A0000576_ChIPseq_AS
# ls *.fastq.gz | parallel --progress --eta -j 16 "rsync --ignore-existing -arzPR {} ${fastqDir}"
# cd -

# 1.2) Run the pipeline
cd ${projDir}
nextflow run /home/rad/users/gaurav/projects/workflows/nfchipseq --input ${projDir}/${projName}_SampleSheet.csv --genome GRCh38 --narrow_peak --single_end -name ${projName} --skip_diff_analysis 0 --deseq2_vst 0

# #########################################################################################
# # 2) PREPROCESS OUTPUT OF THE PIPELINE FOR DOWNSTREAM ANALYSIS
# #########################################################################################
# # 2.1) Parse initial output of the pipeline and generate input file
# echo "bash scripts/generate_rawCount_mergedPeak_files.sh ${projDir} narrowPeak"
# bash scripts/generate_rawCount_mergedPeak_files.sh ${projDir} "narrowPeak"

# #########################################################################################
# # 3) DOWNSTREAM ANALYSIS
# #########################################################################################
# # 3.1) Output paramerts for downstream analysis
# projDir="${outdir}/${user}/${projName}"
# bamDir=${projDir}/results/bwa/mergedLibrary
# analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
# bname="${analysisDir}/${projName}_consensus_peaks"
# rawCountsTxtFile="${analysisDir}/${bname}_rawCounts.txt"
# peaksAnnTxtFile="${analysisDir}/${bname}_annotation.txt"
# origAnnFile="${projDir}/analysis/interimFiles/${projName}_consensus_peaks.mLb.clN.boolean.txt"

# # 3.2) Parse concensus raw matrix and boolean matrix to get annoation files
# echo "bash scripts/parse_nfchip_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${outdir} ${bamDir} ${origAnnFile} ${jobdir}"
# bash scripts/parse_nfchip_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${outdir} ${bamDir} ${origAnnFile} ${jobdir}
# # Output files are:
# echo "- Consensus bed file: ${consensusPeaksBed}"
# echo "- Raw peaks count   : ${rawCountsTxtFile}"
# echo "- Peaks annotation  : ${peaksAnnTxtFile}"


# 3.3) Custom peakcalling with calculated threshold using the --cutoff-analysis parameter
# Output paramerts for downstream analysis
projDir="${outdir}/${user}/${projName}"
bamDir=${projDir}/results/bwa/mergedLibrary
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
bname="${analysisDir}/${projName}_consensus_peaks"

# The output columns are: 
  # 1) p-score (log10(-p)) cutoff, 
  # 2) q-score (log10(-q)) cutoff, 
  # 3) number of peaks called, 
  # 4) total length in basepairs of peaks called, 
  # 5) average length of peak in basepair

customPeaksDir="${analysisDir}/customPeaks"; mkdir -p ${customPeaksDir}

# Cat the following commands in a temporary file and run using gnu parallel
# # modify manually for H3k4me1 marks to run peaks as broad peaks
# for pd in 0_005 0_0001 0_00005 0_000001
for pd in 0_005
do
  pvalue="${pd//_/.}"
  echo ${pvalue}
  peaksOutDir=${customPeaksDir}/p${pd}
  macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t  ${bamDir}/5B4-MYC_R1.mLb.clN.sorted.bam ${bamDir}/5B4-MYC_R2.mLb.clN.sorted.bam -c ${bamDir}/5B4-input_R1.mLb.clN.sorted.bam ${bamDir}/5B4-input_R2.mLb.clN.sorted.bam -n 5B4-MYC_p${pd} --outdir ${peaksOutDir}
done

# After cutoff analysis
for pd in 0_0005
do
  pvalue="${pd//_/.}"
  echo ${pvalue}
  peaksOutDir=${customPeaksDir}/p${pd}
	# myc/input
  macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/5B4-MYC_R1.mLb.clN.sorted.bam ${bamDir}/5B4-MYC_R2.mLb.clN.sorted.bam -c ${bamDir}/5B4-input_R1.mLb.clN.sorted.bam ${bamDir}/5B4-input_R2.mLb.clN.sorted.bam -n 5B4-MYC_over_input_p${pd} --outdir ${peaksOutDir}
	# myc/IgG
  macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/5B4-MYC_R1.mLb.clN.sorted.bam ${bamDir}/5B4-MYC_R2.mLb.clN.sorted.bam -c ${bamDir}/5B4-IgG_R1.mLb.clN.sorted.bam ${bamDir}/5B4-IgG_R2.mLb.clN.sorted.bam -n 5B4-MYC_over_IgG_p${pd} --outdir ${peaksOutDir}
	# 5B4-IgG/input
	macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/5B4-IgG_R1.mLb.clN.sorted.bam ${bamDir}/5B4-IgG_R2.mLb.clN.sorted.bam -c ${bamDir}/5B4-input_R1.mLb.clN.sorted.bam ${bamDir}/5B4-input_R2.mLb.clN.sorted.bam -n 5B4-IgG_over_input_p${pd} --outdir ${peaksOutDir}
	
	# Cas9-MondoA/input
	macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-MondoA_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-MondoA_R1.mLb.clN.sorted.bam -c ${bamDir}/Cas9-input_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-input_R2.mLb.clN.sorted.bam -n Cas9-MondoA_over_input_p${pd} --outdir ${peaksOutDir}
	# Cas9-MYC/input
	macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-MYC_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-MYC_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-input_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-input_R2.mLb.clN.sorted.bam -n Cas9-MYC_over_input_p${pd} --outdir ${peaksOutDir}
	# Cas9-RNAPolymeraseII/input
	macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-RNAPolymeraseII_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-RNAPolymeraseII_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-input_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-input_R2.mLb.clN.sorted.bam -n Cas9-RNAPolymeraseII_over_input_p${pd} --outdir ${peaksOutDir}

	# Cas9-MondoA/IgG
	macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-MondoA_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-MondoA_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-IgG_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-IgG_R2.mLb.clN.sorted.bam -n Cas9-MondoA_over_IgG_p${pd} --outdir ${peaksOutDir}
	# Cas9-MondoA/IgG
	macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-MYC_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-MYC_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-IgG_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-IgG_R2.mLb.clN.sorted.bam -n Cas9-MYC_over_IgG_p${pd} --outdir ${peaksOutDir}
	# Cas9-RNAPolymeraseII/IgG
	macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-RNAPolymeraseII_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-RNAPolymeraseII_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-IgG_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-IgG_R2.mLb.clN.sorted.bam -n Cas9-MYC_over_IgG_p${pd} --outdir ${peaksOutDir}

	# Cas9-IgG/input
	macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-IgG_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-IgG_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-input_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-input_R2.mLb.clN.sorted.bam -n Cas9-MYC_over_input_p${pd} --outdir ${peaksOutDir}
done


####################################################################################
# Final peak calling using the threshold called with Rupert and Christine
####################################################################################
pd="0_008"; pvalue="${pd//_/.}"; peaksOutDir=${customPeaksDir}/customThreshold
# myc/input
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/5B4-MYC_R1.mLb.clN.sorted.bam ${bamDir}/5B4-MYC_R2.mLb.clN.sorted.bam -c ${bamDir}/5B4-input_R1.mLb.clN.sorted.bam ${bamDir}/5B4-input_R2.mLb.clN.sorted.bam -n 5B4-MYC_over_input_p${pd} --outdir ${peaksOutDir}
# myc/IgG
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/5B4-MYC_R1.mLb.clN.sorted.bam ${bamDir}/5B4-MYC_R2.mLb.clN.sorted.bam -c ${bamDir}/5B4-IgG_R1.mLb.clN.sorted.bam ${bamDir}/5B4-IgG_R2.mLb.clN.sorted.bam -n 5B4-MYC_over_IgG_p${pd} --outdir ${peaksOutDir}
# Cas9-Myc/input
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-MYC_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-MYC_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-input_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-input_R2.mLb.clN.sorted.bam -n Cas9-MYC_over_input_p${pd} --outdir ${peaksOutDir}

pd="0_00025"; pvalue="${pd//_/.}"; peaksOutDir=${customPeaksDir}/customThreshold
# Cas9-Myc/IgG
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-MYC_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-MYC_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-IgG_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-IgG_R2.mLb.clN.sorted.bam -n Cas9-MYC_over_IgG_p${pd} --outdir ${peaksOutDir}

pd="0_001"; pvalue="${pd//_/.}"; peaksOutDir=${customPeaksDir}/customThreshold
# Cas9-MondoA/IgG
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-MondoA_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-MondoA_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-IgG_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-IgG_R2.mLb.clN.sorted.bam -n Cas9-MondoA_over_IgG_p${pd} --outdir ${peaksOutDir}

pd="0_004"; pvalue="${pd//_/.}"; peaksOutDir=${customPeaksDir}/customThreshold
# Cas9-MondoA/input
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-MondoA_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-MondoA_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-input_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-input_R2.mLb.clN.sorted.bam -n Cas9-MondoA_over_input_p${pd} --outdir ${peaksOutDir}

pd="0_0005"; pvalue="${pd//_/.}"; peaksOutDir=${customPeaksDir}/customThreshold
# 5B4-IgG/input
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/5B4-IgG_R1.mLb.clN.sorted.bam ${bamDir}/5B4-IgG_R2.mLb.clN.sorted.bam -c ${bamDir}/5B4-input_R1.mLb.clN.sorted.bam ${bamDir}/5B4-input_R2.mLb.clN.sorted.bam -n 5B4-IgG_over_input_p${pd} --outdir ${peaksOutDir}
# Cas9-IgG/input
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-IgG_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-IgG_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-input_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-input_R2.mLb.clN.sorted.bam -n Cas9-IgG_over_input_p${pd} --outdir ${peaksOutDir}
# Cas9-RNAPolymeraseII/IgG
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-RNAPolymeraseII_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-RNAPolymeraseII_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-IgG_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-IgG_R2.mLb.clN.sorted.bam -n Cas9-RNAPolymeraseII_over_IgG_p${pd} --outdir ${peaksOutDir}
# Cas9-RNAPolymeraseII/input
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-RNAPolymeraseII_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-RNAPolymeraseII_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-input_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-input_R2.mLb.clN.sorted.bam -n Cas9-RNAPolymeraseII_over_input_p${pd} --outdir ${peaksOutDir}





pd="0_00025"; pvalue="${pd//_/.}"; peaksOutDir=${customPeaksDir}/customThreshold
# Cas9-MondoA/IgG
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-MondoA_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-MondoA_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-IgG_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-IgG_R2.mLb.clN.sorted.bam -n Cas9-MondoA_over_IgG_p${pd} --outdir ${peaksOutDir}

pd="0_001"; pvalue="${pd//_/.}"; peaksOutDir=${customPeaksDir}/customThreshold
# Cas9-MondoA/input
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-MondoA_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-MondoA_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-input_R1.mLb.clN.sorted.bam ${bamDir}/Cas9-input_R2.mLb.clN.sorted.bam -n Cas9-MondoA_over_input_p${pd} --outdir ${peaksOutDir}
