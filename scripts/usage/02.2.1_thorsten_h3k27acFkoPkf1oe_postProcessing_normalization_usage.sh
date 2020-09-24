cd /home/rad/users/gaurav/projects/workflows/nfchipseq

# Parameters for the script
species="mouse"
user="thorsten"
projName="h3k27acFkoPkf1oe"
outdir="/media/rad/HDD1/nfchip"
jobdir="/home/rad/users/gaurav/projects/workflows/nfchipseq"

#########################################################################################
# 1) RUN THE PIPELINE
#########################################################################################

# Parameters to run the pipeline
projDir="${outdir}/${user}/${projName}"
fastqDir="${projDir}/fastq"; mkdir -p ${fastqDir}

# # 1.1) Copy the fastq files
# cd /media//nas/raw/TUM_NextSeq/200724_NB501802_0272_AHLGJYBGXF/Data/Intensities/BaseCalls/A0000525_TK-ChIP-seq
# ls *.fastq.gz | parallel --progress --eta -j 16 "rsync --ignore-existing -arzPR {} ${fastqDir}"
# cd -

# 1.2) Run the pipeline
cd ${projDir}
nextflow run /home/rad/users/gaurav/projects/workflows/nfchipseq --input ${projDir}/h3k27acFkoPkf1oe_SampleSheet.csv --genome GRCm38 --narrow_peak --single_end -name ${projName}

#########################################################################################
# 2) PREPROCESS OUTPUT OF THE PIPELINE FOR DOWNSTREAM ANALYSIS
#########################################################################################
# 2.1) Parse initial output of the pipeline and generate input file
echo "bash scripts/generate_rawCount_mergedPeak_files.sh ${projDir} narrowPeak"
bash scripts/generate_rawCount_mergedPeak_files.sh ${projDir} "narrowPeak"

#########################################################################################
# 3) DOWNSTREAM ANALYSIS
#########################################################################################
# 3.1) Output paramerts for downstream analysis
projDir="${outdir}/${user}/${projName}"
bamDir=${projDir}/results/bwa/mergedLibrary
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
bname="${analysisDir}/${projName}_consensus_peaks"
rawCountsTxtFile="${analysisDir}/${bname}_rawCounts.txt"
peaksAnnTxtFile="${analysisDir}/${bname}_annotation.txt"
origAnnFile="${projDir}/analysis/interimFiles/${projName}_consensus_peaks.mLb.clN.boolean.txt"

# 3.2) Parse concensus raw matrix and boolean matrix to get annoation files
echo "bash scripts/parse_nfchip_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${outdir} ${bamDir} ${origAnnFile} ${jobdir}"
bash scripts/parse_nfchip_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${outdir} ${bamDir} ${origAnnFile} ${jobdir}

# Output files are:
echo "- Consensus bed file: ${consensusPeaksBed}"
echo "- Raw peaks count   : ${rawCountsTxtFile}"
echo "- Peaks annotation  : ${peaksAnnTxtFile}"

# 3.3) Custom peakcalling with calculated threshold using the --cutoff-analysis parameter
# The output columns are: 
  # 1) p-score (log10(-p)) cutoff, 
  # 2) q-score (log10(-q)) cutoff, 
  # 3) number of peaks called, 
  # 4) total length in basepairs of peaks called, 
  # 5) average length of peak in basepair

customPeaksDir="${analysisDir}/customPeaks"; mkdir -p ${customPeaksDir}

# Cat the following commands in a temporary file and run using gnu parallel
# # modify manually for H3k4me1 marks to run peaks as broad peaks
for pd in 0_0001 0_00005 0_000001
do
  pvalue="${pd//_/.}"
  echo ${pvalue}
  peaksOutDir=${customPeaksDir}/p${pd}
  macs2 callpeak -f BAM --seed=39751 -g mm  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t  /media/rad/HDD1/nfchip/thorsten/h3k27acFkoPkf1oe/results/bwa/mergedLibrary/A0000525-FKO-372-H3K27ac_R1.mLb.clN.sorted.bam -c/media/rad/HDD1/nfchip/thorsten/h3k27acFkoPkf1oe/results/bwa/mergedLibrary/A0000525-FKO-372-Input_R1.mLb.clN.sorted.bam -n A0000525-FKO-372-H3K27ac_mmu_chipseq_se_R1_${pd} --outdir ${peaksOutDir}
  macs2 callpeak -f BAM --seed=39751 -g mm  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t  /media/rad/HDD1/nfchip/thorsten/h3k27acFkoPkf1oe/results/bwa/mergedLibrary/A0000525-FKO-387-H3K27ac_R1.mLb.clN.sorted.bam -c /media/rad/HDD1/nfchip/thorsten/h3k27acFkoPkf1oe/results/bwa/mergedLibrary/A0000525-FKO-387-Input_R1.mLb.clN.sorted.bam -n A0000525-FKO-387-H3K27ac_mmu_chipseq_se_R1_${pd} --outdir ${peaksOutDir}
  macs2 callpeak -f BAM --seed=39751 -g mm  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t  /media/rad/HDD1/nfchip/thorsten/h3k27acFkoPkf1oe/results/bwa/mergedLibrary/A0000525-PKF1OE-057-H3K27ac_R1.mLb.clN.sorted.bam -c /media/rad/HDD1/nfchip/thorsten/h3k27acFkoPkf1oe/results/bwa/mergedLibrary/A0000525-PKF1OE-057-Input_R1.mLb.clN.sorted.bam -n A0000525-PKF1OE-057-H3K27ac_mmu_chipseq_se_R1_${pd} --outdir ${peaksOutDir}
  macs2 callpeak -f BAM --seed=39751 -g mm  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t  /media/rad/HDD1/nfchip/thorsten/h3k27acFkoPkf1oe/results/bwa/mergedLibrary/A0000525-PKF1OE-062-H3K27ac_R1.mLb.clN.sorted.bam -c /media/rad/HDD1/nfchip/thorsten/h3k27acFkoPkf1oe/results/bwa/mergedLibrary/A0000525-PKF1OE-062-Input_R1.mLb.clN.sorted.bam -n A0000525-PKF1OE-062-H3K27ac_mmu_chipseq_se_R1_${pd} --outdir ${peaksOutDir}
done
# cat > ${customPeaksDir}/${projName}_macs2_callpeak_p005.sh
# # Cltr+C
# chmod 775 ${customPeaksDir}/${projName}_macs2_callpeak_p005.sh
# parallel :::: ${customPeaksDir}/${projName}_macs2_callpeak_p005.sh

# Annotate using chipseeker
for pd in 0_005 0_0001 0_00005 0_000001 
do
 peaksOutDir=${customPeaksDir}/p${pd}
 for peaksInputFile in $(ls ${peaksOutDir}/*.narrowPeak)
 do
  peaksAnnTxtFile=${peaksOutDir}/$(basename ${peaksInputFile} .narrowPeak)_annotation.txt
  Rscript ${jobdir}/scripts/R_annotate_peaks_chipseeker.R -if=${peaksInputFile} -of=${peaksAnnTxtFile} -sp=${species} -mf -ac
 done
done

0.0001 

# 3.3) Run CRCs
# 3.3.1) Get sample folders with respective bam and peaks file
# General pipeline parameters
species="mouse"
user="thorsten"
projName="h3k27acFkoPkf1oe"
outdir="/media/rad/HDD1/nfchip"
jobdir="/home/rad/users/gaurav/projects/workflows/nfchipseq"
# Script specific parameters
# projDir="/media/rad/HDD1/nfchip/${user}/${projName}/gjchip"
projDir="${outdir}/${user}/${projName}"
bamDir=${projDir}/results/bwa/mergedLibrary
analysisDir="${projDir}/analysis"
customPeaksDir="${analysisDir}/customPeaks"
peaksDir=${customPeaksDir}/p0_0001
crcDir="${analysisDir}/crcs"; mkdir -p ${crcDir}

cd /home/rad/users/gaurav/projects/ctrc
bash scripts/create_nf_sample_dirs.sh ${projName} ${crcDir} ${bamDir} ${peaksDir}

rm -rf ${crcDir}/*Input*

# 3.3.2) Call SuperEnhancers using ROSE2 and Core Regulatory Circuits (CRCs) using crc2
speciesGenomeAssembly="mm10" 
bash scripts/get_crcs.sh ${crcDir} ${projName} ${speciesGenomeAssembly} nf

# 3.3.3) Run the CRCs wrapper
scriptsdir="/home/rad/users/gaurav/projects/ctrc/scripts/03_crcs/${projName}"
rm -rf ${scriptsdir}/CRCLogs_crcs.sh
# cmd="parallel ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}
for s in ${scriptsdir}/*.sh; do echo ${s}; chmod 775 ${s}; bash ${s}; done;
cd -

























# 2.2) Generate consensus peaks region file and raw count matrix for:
projDir="${outdir}/${user}/${projName}"
bamDir=${projDir}/results/bwa/mergedLibrary
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
customPeaksDir="${analysisDir}/customPeaks/customConsensus"; mkdir -p ${customPeaksDir}
# 2.2.1) Over expression samples at pvalue 0.005
# 2.2.2) Knockout samples at pvalue 0.0001
# 2.2.3) Overall consensus of OE and FKO for above pvalues

for pd in OE FKO OEFKO;
do 
  peaksOutDir=${customPeaksDir}/${pd}
  consensusDir="${peaksOutDir}/consensus"
  echo -e "- ${projDir}\n- ${peaksOutDir}\n- ${consensusDir}"
  bash /home/rad/users/gaurav/projects/workflows/nfchipseq/scripts/generate_rawCount_mergedPeak_files.sh ${projDir} narrowPeak ${peaksOutDir} ${consensusDir}

  bname="${projName}_p${pvalueDir}_${pd}_consensus_peaks"
  rawCountsTxtFile="${consensusDir}/${bname}_rawCounts.txt"
  peaksAnnTxtFile="${consensusDir}/${bname}_annotation.txt"
  origAnnFile="${consensusDir}/interimFiles/${projName}_consensus_peaks.mLb.clN.boolean.txt"
  scriptDir="/home/rad/users/gaurav/projects/workflows/nfatacseq"
  bash /home/rad/users/gaurav/projects/workflows/nfatacseq/scripts/parse_nfatacseq_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${consensusDir} ${bamDir} ${origAnnFile} ${scriptDir}
  echo ""
done


# bash /home/rad/users/gaurav/projects/workflows/nfatacseq/scripts/parse_nfatacseq_consensus_peaks_annotation.sh mouse thorsten foxp1FkoPkf1oe /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/analysis/customPeaks/customConsensus/OE/consensus /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/results/bwa/mergedLibrary /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/analysis/customPeaks/customConsensus/OE/consensus/interimFiles/foxp1FkoPkf1oe_consensus_peaks.mLb.clN.boolean.txt /home/rad/users/gaurav/projects/workflows/nfatacseq
# bash /home/rad/users/gaurav/projects/workflows/nfatacseq/scripts/parse_nfatacseq_consensus_peaks_annotation.sh mouse thorsten foxp1FkoPkf1oe /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/analysis/customPeaks/customConsensus/FKO/consensus /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/results/bwa/mergedLibrary /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/analysis/customPeaks/customConsensus/FKO/consensus/interimFiles/foxp1FkoPkf1oe_consensus_peaks.mLb.clN.boolean.txt /home/rad/users/gaurav/projects/workflows/nfatacseq
# bash /home/rad/users/gaurav/projects/workflows/nfatacseq/scripts/parse_nfatacseq_consensus_peaks_annotation.sh mouse thorsten foxp1FkoPkf1oe /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/analysis/customPeaks/customConsensus/OEFKO/consensus /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/results/bwa/mergedLibrary /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/analysis/customPeaks/customConsensus/OEFKO/consensus/interimFiles/foxp1FkoPkf1oe_consensus_peaks.mLb.clN.boolean.txt /home/rad/users/gaurav/projects/workflows/nfatacseq