cd /home/rad/users/gaurav/projects/workflows/nfchipseq

# Parameters for the script
species="human"
genome="hs"
user="anja"
projName="nfTALLhs"
outdir="/media/rad/HDD1/atacseq"
jobdir="/home/rad/users/gaurav/projects/workflows/nfchipseq"

#########################################################################################
# 1) RUN THE PIPELINE
#########################################################################################

# Parameters to run the pipeline
projDir="${outdir}/${user}/${projName}"
fastqDir="${projDir}/fastq"; mkdir -p ${fastqDir}

# # 1.1) Copy the fastq files
# cd /media/nas/temporary/PUB_CRCs/tallCRCs/human/fastq
# ls *.gz | parallel --progress --eta -j 32 "rsync --ignore-existing -arzRP {} ${fastqDir}"
# cd -

# 1.2) Run the pipeline
cd ${projDir}
nextflow run /home/rad/users/gaurav/projects/workflows/nfatacseq --input ${projDir}/${projName}_design.csv --genome GRCh38 --single_end --narrow_peak -resume -name ${projName}
cd -

#########################################################################################
# 3) DOWNSTREAM ANALYSIS
#########################################################################################
# 3.1) Output paramerts for downstream analysis
projDir="${outdir}/${user}/${projName}"
bamDir=${projDir}/results/bwa/mergedLibrary
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
bname="${analysisDir}/${projName}_consensus_peaks"

# 3.3) Custom peakcalling with calculated threshold using the --cutoff-analysis parameter
# The output columns are: 
  # 1) p-score (log10(-p)) cutoff, 
  # 2) q-score (log10(-q)) cutoff, 
  # 3) number of peaks called, 
  # 4) total length in basepairs of peaks called, 
  # 5) average length of peak in basepair

origPeaksDir=${projDir}/results/bwa/mergedLibrary/macs/narrowPeak
customPeaksDir="${analysisDir}/customPeaks"; mkdir -p ${customPeaksDir}

# Cat the following commands in a temporary file and run using gnu parallel
# # modify manually for H3k4me1 marks to run peaks as broad peaks
pd=0_005
pvalue="${pd//_/.}"
echo ${pvalue}
peaksOutDir=${customPeaksDir}/p${pd}; mkdir -p ${peaksOutDir}
for b in $(ls ${bamDir}/*.bam);
do 
  echo ${b}
  echo macs2 callpeak -f BAM --seed=39751 -g ${genome} --keep-dup auto -p ${pvalue} --cutoff-analysis -t ${b} -n $(basename ${b} _R1.mLb.clN.sorted.bam)_${pd} --outdir ${peaksOutDir}
done

# Cutoffs using cutoff-analysis:
# +---------------------------------------------------------------+--------+--------+--------+
# |                            Samples                            | pscore | npeaks | pvalue |
# +---------------------------------------------------------------+--------+--------+--------+
# | H3K27ac_DND41_0_005_cutoff_analysis.txt                       |     10 |   4947 | 1E-10  |
# | H3K27ac_DND41_R2.mLb.clN.sorted.bam_0_005_cutoff_analysis.txt |     10 |   5095 | 1E-10  |
# | H3K27ac_DND41_R3.mLb.clN.sorted.bam_0_005_cutoff_analysis.txt |     10 |   4904 | 1E-10  |
# | H3K27ac_DND41_R4.mLb.clN.sorted.bam_0_005_cutoff_analysis.txt |     10 |   5000 | 1E-10  |
# | H3K27ac_RPMI8402_0_005_cutoff_analysis.txt                    |     10 |   6570 | 1E-10  |
# | H3K27ac_CCRFCEM_0_005_cutoff_analysis.txt                     |     15 |        | 1E-15  |
# | H3K27ac_DND41_Broad_0_005_cutoff_analysis.txt                 |     15 |        | 1E-15  |
# | H3K27ac_Jurkat_0_005_cutoff_analysis.txt                      |     20 |        | 1E-20  |
# | H3K27ac_Loucy_0_005_cutoff_analysis.txt                       |     20 |        | 1E-20  |
# | H3K27ac_MOLT4_0_005_cutoff_analysis.txt                       |     20 |   3250 | 1E-20  |
# +---------------------------------------------------------------+--------+--------+--------+

pd=1E-10; pvalue=1E-10; peaksOutDir=${customPeaksDir}/p${pd}; mkdir -p ${peaksOutDir}
for s in  H3K27ac_DND41_ H3K27ac_DND41_R2_ H3K27ac_DND41_R3 H3K27ac_DND41_R4 H3K27ac_RPMI8402_ ; do bam=$(find ${bamDir} -iname *${s}*.bam); macs2 callpeak -f BAM --seed=39751 -g ${genome} --keep-dup auto -p ${pvalue} --cutoff-analysis -t ${bam} -n $(basename ${bam} mLb.clN.sorted.bam)_p${pd} --outdir ${peaksOutDir}; done

pd=1E-15; pvalue=1E-15; peaksOutDir=${customPeaksDir}/p${pd}; mkdir -p ${peaksOutDir}
for s in  H3K27ac_CCRFCEM_ H3K27ac_DND41_Broad_ ; do bam=$(find ${bamDir} -iname *${s}*.bam); macs2 callpeak -f BAM --seed=39751 -g ${genome} --keep-dup auto -p ${pvalue} --cutoff-analysis -t ${bam} -n $(basename ${bam} mLb.clN.sorted.bam)_p${pd} --outdir ${peaksOutDir}; done

pd=1E-20; pvalue=1E-20; peaksOutDir=${customPeaksDir}/p${pd}; mkdir -p ${peaksOutDir}
for s in  H3K27ac_Jurkat_ H3K27ac_Loucy_ H3K27ac_MOLT4_ ; do bam=$(find ${bamDir} -iname *${s}*.bam); macs2 callpeak -f BAM --seed=39751 -g ${genome} --keep-dup auto -p ${pvalue} --cutoff-analysis -t ${bam} -n $(basename ${bam} mLb.clN.sorted.bam)_p${pd} --outdir ${peaksOutDir}; done

# 3.3) Run CRCs
# 3.3.1) Get sample folders with respective bam and peaks file
# General pipeline parameters
species="human"
genome="hs"
user="anja"
projName="nfTALLhs"
outdir="/media/rad/HDD1/atacseq"
jobdir="/home/rad/users/gaurav/projects/workflows/nfchipseq"
projDir="${outdir}/${user}/${projName}"
bamDir=${projDir}/results/bwa/mergedLibrary
analysisDir="${projDir}/analysis"
customPeaksDir="${analysisDir}/customPeaks"
crcDir="${analysisDir}/crcs"; mkdir -p ${crcDir}
peaksDir=${crcDir}/peaks; mkdir -p ${peaksDir}
# Create a softlink for peaks in the peaksdir
for pd in 1E-10 1E-15 1E-20;
do
  peaksInDir=${customPeaksDir}/p${pd}
  for p in $(ls ${peaksInDir}/*);
  do
    echo ln -s ${p} ${peaksDir}/$(basename ${p})
    ln -s ${p} ${peaksDir}/$(basename ${p})
  done
done

# We can only run the crc scripts from the crc project directory
cd /home/rad/users/gaurav/projects/ctrc
bash scripts/create_nf_sample_dirs.sh ${projName} ${crcDir} ${bamDir} ${peaksDir}

rm -rf ${crcDir}/*Input*

# 3.3.2) Call SuperEnhancers using ROSE2 and Core Regulatory Circuits (CRCs) using crc2
speciesGenomeAssembly="hg19" 
bash scripts/get_crcs.sh ${crcDir} ${projName} ${speciesGenomeAssembly}

# 3.3.3) Run the CRCs wrapper
scriptsdir="/home/rad/users/gaurav/projects/ctrc/scripts/03_crcs/${projName}"
# Remove useless files
rm -rf ${scriptsdir}/{CRCLogs,peaks}_crcs.sh
# cmd="parallel ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}
for s in ${scriptsdir}/*.sh; do echo ${s}; chmod 775 ${s}; bash ${s}; done;
cd -

# Zip the crc output excluding certain files and directories
cd ${analysisDir}
zip -r crcs.zip crcs -x *.bam* crcsbamliquidator/**/* crcsFIMO/**/*
cd -
