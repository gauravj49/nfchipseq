cd /home/rad/users/gaurav/projects/workflows/nfchipseq

# Parameters for the script
species="mouse"
user="christine"
projName="pdacBatch1"
outdir="/media/rad/HDD1/nfchip"
jobdir="/home/rad/users/gaurav/projects/workflows/nfchipseq"


#########################################################################################
# 1) RUN THE PIPELINE
#########################################################################################

# Parameters to run the pipeline
projDir="${outdir}/${user}/${projName}"
fastqDir="${projDir}/fastq"; mkdir -p ${fastqDir}

# # 1.1) Copy the fastq files
# cd /home/rad/media/nas/raw/TUM_NextSeq/200519_NB501802_0260_AHTYLNBGXF/Data/Intensities/BaseCalls
# ls A0000504-FKO-372-Foxp1-IP_S5_R1_001.fastq.gz A0000504-FKO-372-Input_S9_R1_001.fastq.gz A0000504-FKO-387-Foxp1-IP_S4_R1_001.fastq.gz A0000504-FKO-387-Input_S8_R1_001.fastq.gz A0000504-PKF1OE-057-Foxp1-IP_S3_R1_001.fastq.gz A0000504-PKF1OE-057-Input_S7_R1_001.fastq.gz A0000504-PKF1OE-062-Foxp1-IP_S2_R1_001.fastq.gz A0000504-PKF1OE-062-Input_S6_R1_001.fastq.gz | parallel --progress --eta -j 16 "rsync --ignore-existing -arzPR {} ${fastqDir}"
# cd -

# 1.2) Run the pipeline
cd ${projDir}
nextflow run /home/rad/users/gaurav/projects/workflows/nfchipseq --input ${projDir}/${projName}_SampleSheet.csv --genome GRCm38 --narrow_peak --single_end -name ${projName}

#########################################################################################
# 1) DOWNSTREAM ANALYSIS
#########################################################################################
# 2.1) Output paramerts for downstream analysis
projDir="${outdir}/${user}/${projName}"
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
consensusPeaksBed="${analysisDir}/${projName}_$(basename ${analysisDir})_consensus_peaks.bed"
rawCountsTxtFile="${analysisDir}/$(basename ${consensusPeaksBed} .bed)_rawCounts.txt"
peaksAnnTxtFile="${analysisDir}/$(basename ${consensusPeaksBed} .bed)_annotation.txt"
origConsFile="${projDir}/results/bwa/mergedLibrary/macs/narrowPeak/consensus/abcam/consensus_peaks.mLb.clN.featureCounts.txt"
origAnnFile="${projDir}/results/bwa/mergedLibrary/macs/narrowPeak/consensus/abcam/abcam.consensus_peaks.boolean.txt"

# 2.2) Parse concensus raw matrix and boolean matrix to get annoation files
echo "bash scripts/parse_nfatac_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${outdir} ${origConsFile} ${origAnnFile} ${jobdir}"
bash scripts/parse_nfatac_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${outdir} ${origConsFile} ${origAnnFile} ${jobdir}

# Output files are:
echo "- Consensus bed file: ${consensusPeaksBed}"
echo "- Raw peaks count   : ${rawCountsTxtFile}"
echo "- Peaks annotation  : ${peaksAnnTxtFile}"
