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

# Annotate using chipseeker
peaksOutDir=${customPeaksDir}/customThreshold
for peaksInputFile in $(ls ${peaksOutDir}/*.narrowPeak)
do
	peaksAnnTxtFile=${peaksOutDir}/$(basename ${peaksInputFile} .narrowPeak)_annotation.txt
	echo Rscript ${jobdir}/scripts/R_annotate_peaks_chipseeker.R -if=${peaksInputFile} -of=${peaksAnnTxtFile} -sp=${species} -mf
	Rscript ${jobdir}/scripts/R_annotate_peaks_chipseeker.R -if=${peaksInputFile} -of=${peaksAnnTxtFile} -sp=${species} -mf
done

# Call individual peaks for replicates individually
pd="0_008"; pvalue="${pd//_/.}"; peaksOutDir=${customPeaksDir}/individualCustomPeaks
# myc/input
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/5B4-MYC_R1.mLb.clN.sorted.bam -c ${bamDir}/5B4-input_R1.mLb.clN.sorted.bam -n 5B4-MYC_R1_over_input_R1_p${pd} --outdir ${peaksOutDir}
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/5B4-MYC_R2.mLb.clN.sorted.bam -c ${bamDir}/5B4-input_R2.mLb.clN.sorted.bam -n 5B4-MYC_R2_over_input_R2_p${pd} --outdir ${peaksOutDir}

# Cas9-Myc/input
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-MYC_R1.mLb.clN.sorted.bam -c ${bamDir}/Cas9-input_R1.mLb.clN.sorted.bam -n Cas9-MYC_R1_over_input_R1_p${pd} --outdir ${peaksOutDir}
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-MYC_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-input_R2.mLb.clN.sorted.bam -n Cas9-MYC_R2_over_input_R2_p${pd} --outdir ${peaksOutDir}

pd="0_001"; pvalue="${pd//_/.}"; peaksOutDir=${customPeaksDir}/individualCustomPeaks
# Cas9-MondoA/input
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-MondoA_R1.mLb.clN.sorted.bam -c ${bamDir}/Cas9-input_R1.mLb.clN.sorted.bam -n Cas9-MondoA_R1_over_input_R1_p${pd} --outdir ${peaksOutDir}
macs2 callpeak -f BAM --seed=39751 -g hs  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t ${bamDir}/Cas9-MondoA_R2.mLb.clN.sorted.bam -c ${bamDir}/Cas9-input_R2.mLb.clN.sorted.bam -n Cas9-MondoA_R2_over_input_R2_p${pd} --outdir ${peaksOutDir}

# Copy the samplesheet to the server
scp -r mondoAKO_SampleSheet_DiffBind.csv iws4:/media/rad/HDD1/nfchip/rupert/mondoAKO/

# Run Diffbind
# https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2020/ChIPSeq/scripts/ChIP_Practical3_DiffBind.html
R

# # In case we do not have these packages
# install.packages("tidyverse")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("DiffBind")

# load required packages:
library(DiffBind)
library(tidyverse)
library(rtracklayer)

# 1) Step 1: Reading a peakset
# 1.1) Read a csv file
samples <- read.csv("/media/rad/HDD1/nfchip/rupert/mondoAKO/mondoAKO_SampleSheet_DiffBind.csv")
# 1.2) Look at the loaded metadata
names(samples)
# 1.3) The peaksets are read in using the following DiffBind function, that will construct a new DBA object from the sample sheet.
mondoAKO <- dba(sampleSheet="/media/rad/HDD1/nfchip/rupert/mondoAKO/mondoAKO_SampleSheet_DiffBind.csv")

# 2) Step 2: Occupancy analysis
dba.plotHeatmap(mondoAKO)

# 3) Step 3: Counting reads
# 3.1) Take the alignment files (BAM) and compute count information for each of the peaks/regions in the consensus set
mondoAKO.counted <- dba.count(mondoAKO, summits=250)
mondoAKO.counted
# 6 Samples, 1621 sites in matrix:
#               ID Condition Caller Intervals
# 1     5B4-MYC_R1       MKO counts      1621
# 2     5B4-MYC_R2       MKO counts      1621
# 3 Cas9-MondoA_R1  MondoAKO counts      1621
# 4 Cas9-MondoA_R2  MondoAKO counts      1621
# 5    Cas9-MYC_R1  MondoAWT counts      1621
# 6    Cas9-MYC_R2  MondoAWT counts      1621

# 3.2) With this matrix, the samples can be re-clustered using affinity, rather than occupancy data:
cat("\n- draw a Sample correlation plot using all 1621 consensus sites\n")
SCRplotFile     <- paste0("/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoAKO_sampleCorrelation.jpg", sep='') 
jpeg(filename=SCRplotFile, height=900, width=1000, bg="white", res=100); par(mar=c(20,4,1,1));
dba.plotHeatmap(mondoAKO.counted)
dev.off()

# 3.3) To see how well the samples cluster with one another, we can draw a PCA plot using all 1621 consensus sites.
cat("\n- draw a PCA plot using all 1621 consensus sites\n")
PCAplotFile     <- paste0("/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoAKO_PCA.jpg", sep='') 
jpeg(filename=PCAplotFile, height=900, width=1000, bg="white", res=100); par(mar=c(20,4,1,1));
dba.plotPCA(mondoAKO.counted,  attributes=DBA_CONDITION, label=DBA_ID)
dev.off()

# 3.4) Open the file from terminal
eog mondoAKO_PCA.jpg
eog mondoAKO_sampleCorrelation.jpg

# 4) Step 4: Differential binding affinity analysis
# 4.1) 	we have to let DiffBind know how we want to group our samples. 
# 			In our case we will group based on condition. 
# 			This is done using the function, as follows:
mondoAKO.counted <- dba.contrast(mondoAKO.counted, categories=DBA_CONDITION, minMembers=2)
# 3 Contrasts:
#     Group1 Members1   Group2 Members2
# 1      MKO        2 MondoAKO        2
# 2      MKO        2 MondoAWT        2
# 3 MondoAKO        2 MondoAWT        2
# 4.2) The main differential analysis function using DESeq2 is invoked as follows:
mondoAKO.analysed <- dba.analyze(mondoAKO.counted)

# 4.3) plots
dba.plotHeatmap(mondoAKO.analysed, contrast=2, th=1)
PCAplotFile     <- paste0("/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoAKO_binding_affinity_heatmap.jpg", sep='') 
jpeg(filename=PCAplotFile, height=900, width=1000, bg="white", res=100); par(mar=c(20,4,1,1));
dba.plotHeatmap(mondoAKO.analysed, ColAttributes = DBA_CONDITION, contrast=2, correlations=FALSE, th=1)
dev.off()

# 4.4) Other plots:
# 4.4.1) MA plots
PCAplotFile     <- paste0("/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoAKO_binding_affinity_MA_plot.jpg", sep='') 
jpeg(filename=PCAplotFile, height=900, width=1000, bg="white", res=100); par(mar=c(20,4,1,1));
dba.plotMA(mondoAKO.analysed)
dev.off()
# 4.4.2) Volcano plots
PCAplotFile     <- paste0("/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoAKO_binding_affinity_volcano_plot.jpg", sep='') 
jpeg(filename=PCAplotFile, height=900, width=1000, bg="white", res=100); par(mar=c(20,4,1,1));
dba.plotVolcano(mondoAKO.analysed)
dev.off()
# 4.4.3) PCA plots for contrast=2
PCAplotFile     <- paste0("/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoAKO_binding_affinity_PCA_MKOvsMondoAWT_plot.jpg", sep='') 
jpeg(filename=PCAplotFile, height=900, width=1000, bg="white", res=100); par(mar=c(20,4,1,1));
dba.plotPCA(mondoAKO.analysed, contrast=2, th=1)
dev.off()
# 4.4.4) # Boxplots can give us an idea about the read distribution differences between the classes 
# - in our case the two conditions. The first two boxes show distribution of reads over all 
# differentially bound sites; the middle two show differences on those sites where the affinity 
# increases in Responsive and the two boxes on the right show differences where the affinity increases in Resistant samples.
PCAplotFile     <- paste0("/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoAKO_binding_affinity_box_plot.jpg", sep='') 
jpeg(filename=PCAplotFile, height=900, width=1000, bg="white", res=100); par(mar=c(20,4,1,1));
dba.plotBox(mondoAKO.analysed)
dev.off()

# 5) Report the differentially bound peak regions, identified by DESeq2. 
# These results files contain the genomic coordinates for all consensus sites and 
# statistics for differential enrichment including fold-change, p-value and FDR.
report <- dba.report(mondoAKO.analysed)
report
# The value columns are described below:
# 	Chr: Chromosome of binding site
# 	Start: Starting base position of binding site
# 	End: End base position of binding site
# 	Conc: mean read concentration over all the samples (the default calculation uses log2 normalized ChIP read counts with control read counts subtracted)
# 	Conc_group1: Group 1 Concentration
# 	Conc_group2: Group 2 Concentration
# 	Fold: Fold difference – mean fold difference of binding affinity of group 1 over group 2 (Conc1 - Conc2)
# 	p-valule and *FDR statistic indicating significance of difference
# 	Before writing to file we need to convert it to a data frame so that genomic coordinates get written as columns and not GRanges.
report.df <- as.data.frame(report)  
write.table(report.df, "/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoAKO_report.csv", sep="\t", quote=F, row.names=F)

# 6) 
# PeakChrom	PeakStart	PeakEnd	PeakID	width	strand	Conc	Fold	p.value	FDR
