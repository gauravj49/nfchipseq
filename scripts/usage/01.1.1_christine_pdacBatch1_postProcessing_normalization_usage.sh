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

# 3.3) Run CRCs
# 3.3.1) Get sample folders with respective bam and peaks file
# General pipeline parameters
species="mouse"
user="christine"
projName="pdacBatch1"
jobdir="/home/rad/users/gaurav/projects/workflows/nfchipseq"
# Script specific parameters
projDir="/media/rad/HDD1/nfchip/${user}/${projName}/gjchip"
bamDir="${projDir}/mapping";
peaksDir="${projDir}/peaks"; mkdir -p ${peaksDir}
peaksLogsDir="${projDir}/peaks/logs"; mkdir -p ${peaksLogsDir}
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
crcDir="${analysisDir}/crcs"

cd /home/rad/users/gaurav/projects/ctrc
bash scripts/create_gjchip_sample_dirs.sh ${projName} ${crcDir} ${bamDir} ${peaksDir}

# 3.3.2) Call SuperEnhancers using ROSE2 and Core Regulatory Circuits (CRCs) using crc2
speciesGenomeAssembly="mm10"
bash scripts/get_crcs.sh ${crcDir} ${projName} ${speciesGenomeAssembly}

# 3.3.3) Run the CRCs wrapper
scriptsdir="/home/rad/users/gaurav/projects/ctrc/scripts/03_crcs/${projName}"
rm -rf ${scriptsdir}/CRCLogs_crcs.sh
# cmd="parallel ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}
for s in ${scriptsdir}/*.sh; do echo ${s}; chmod 775 ${s}; bash ${s}; done;
cd -

# # 3.3) Run CRCs
# # 3.3.1) Get sample folders with respective bam and peaks file
# cd /home/rad/users/gaurav/projects/ctrc
# crcDir="${analysisDir}/crcs"
# peaksDir="${projDir}/results/bwa/mergedLibrary/macs/narrowPeak"
# bash scripts/create_nf_sample_dirs.sh ${projName} ${crcDir} ${bamDir} ${peaksDir}

# # 3.3.2) Call SuperEnhancers using ROSE2 and Core Regulatory Circuits (CRCs) using crc2
# speciesGenomeAssembly="mm10"
# bash scripts/get_crcs.sh ${crcDir} ${projName} ${speciesGenomeAssembly}

# # 3.3.3) Run the CRCs wrapper
# scriptsdir="/home/rad/users/gaurav/projects/ctrc/scripts/03_crcs/${projName}"
# rm -rf ${scriptsdir}/CRCLogs_crcs.sh
# # cmd="parallel ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}
# for s in ${scriptsdir}/*.sh; do echo ${s}; chmod 775 ${s}; bash ${s}; done;
# cd -


# 3.3.4) Run Extended network analysis
cd /home/rad/users/gaurav/projects/ctrc/scripts/CLL_TFnetworks_2018
analysisDir="/media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/analysis/crcs/extendedNetworks"
crcDir="/media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/analysis/crcs"
atacPeakDir="/media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/analysis/customPeaks/mergedLibrary/q0_01"
genome="MM10"
# Define the pairs of chip and atac samples
chipSamples=(123abcam-53646-PPT-1 124abcam-53646-LivMet-1 125abcam-53646-LivMet-2 126abcam-53646-LivMet-3 127abcam-5320-PPT-1 128abcam-5320-LungMet-1 129abcam-5320-LivMet-1 130abcam-5320-LivMet-3)
atacSamples=(53646_PPT-1 53646_LivMet-1 53646_LivMet-2 53646_LivMet-3 5320_PPT-1 5320_LungMet-1 5320_LivMet-1 5320_LivMet-3)
# Run the samples
for ((i=0;i<${#chipSamples[@]};i++))
do
	echo python extendedSuperNetwork.py -pn=${chipSamples[$i]} -ad=${analysisDir} -nf=${crcDir}/${chipSamples[$i]}/crcs${chipSamples[$i]}_NODELIST.txt -se=${crcDir}/${chipSamples[$i]}/rose2/${chipSamples[$i]}_summits_SuperEnhancers.table.txt -ap=${atacPeakDir}/${atacSamples[$i]}_005_R1_peaks.broadPeak -gn=${genome}
	python extendedSuperNetwork.py -pn=${chipSamples[$i]} -ad=${analysisDir} -nf=${crcDir}/${chipSamples[$i]}/crcs${chipSamples[$i]}_NODELIST.txt -se=${crcDir}/${chipSamples[$i]}/rose2/${chipSamples[$i]}_summits_SuperEnhancers.table.txt -ap=${atacPeakDir}/${atacSamples[$i]}_005_R1_peaks.broadPeak -gn=${genome}
	echo 
done
cd /home/rad/users/gaurav/projects/workflows/nfchipseq

# 3.3.5) Intersect the extended network file with the super enhancer table
for ((i=0;i<${#chipSamples[@]};i++))
do
	echo ${chipSamples[$i]}
	python -c "import sys; import pandas as pd; import datatable as dt; input_file=sys.argv[1]; input_supenc=sys.argv[2]; outtxt_file=sys.argv[3]; infileDT = dt.fread(input_file, sep='\t', header=True, nthreads=16); supencDT = dt.fread(input_supenc, sep='\t', header=True, nthreads=16); infileDT = dt.fread(input_file, sep='\t', header=True, nthreads=16); supencDT = dt.fread(input_supenc, sep='\t', header=True, nthreads=16); infileDF = infileDT.to_pandas(); infileDF.set_index('REGION_ID', inplace=True); supencDF = supencDT.to_pandas(); supencList = supencDF.columns.tolist(); infileDF = infileDF[infileDF.index.isin(supencList)]; infileDF.reset_index(level=0, inplace=True); infileDFColumns = infileDF.columns.tolist(); infileDFColumns.insert(3,infileDFColumns.pop(0)); infileDF = infileDF[infileDFColumns]; infileDF.to_csv(outtxt_file, header=True, index=False, sep='\t', float_format='%.0f');" ${crcDir}/${chipSamples[$i]}/rose2/${chipSamples[$i]}_summits_SuperEnhancers.table.txt ${crcDir}/extendedNetworks/${chipSamples[$i]}/${chipSamples[$i]}_extendedNetwork.allEnhancers.matrix.txt ${crcDir}/extendedNetworks/${chipSamples[$i]}/${chipSamples[$i]}_allEnhancers_matrix_intersect_superEnhancers_table.txt
	echo 
done
