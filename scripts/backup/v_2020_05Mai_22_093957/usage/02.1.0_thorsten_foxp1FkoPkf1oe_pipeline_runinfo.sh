[rad@imows4 foxp1FkoPkf1oe ] $ nextflow run /home/rad/users/gaurav/projects/workflows/nfchipseq --input ${projDir}/foxp1FkoPkf1oe_SampleSheet.csv --genome GRCm38 --narrow_peak --single_end -name foxp1FkoPkf1oe1
N E X T F L O W  ~  version 20.04.1
Launching `/home/rad/users/gaurav/projects/workflows/nfchipseq/main.nf` [foxp1FkoPkf1oe1] - revision: f362665a84
----------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/chipseq v1.1.0
----------------------------------------------------
Run Name            : foxp1FkoPkf1oe1
Data Type           : Single-End
Design File         : /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/foxp1FkoPkf1oe_SampleSheet.csv
Genome              : GRCm38
Design File         : /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/foxp1FkoPkf1oe_SampleSheet.csv
Genome              : GRCm38
Fasta File          : s3://ngi-igenomes/igenomes//Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa
GTF File            : s3://ngi-igenomes/igenomes//Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf
Gene BED File       : s3://ngi-igenomes/igenomes//Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.bed
BWA Index           : s3://ngi-igenomes/igenomes//Mus_musculus/Ensembl/GRCm38/Sequence/BWAIndex/genome.fa
Blacklist BED       : /home/rad/users/gaurav/projects/workflows/nfchipseq/assets/blacklists/GRCm38-blacklist.bed
MACS2 Genome Size   : 1.87e9
Min Consensus Reps  : 1
MACS2 Narrow Peaks  : Yes
Trim R1             : 0 bp
Trim R2             : 0 bp
Trim 3' R1          : 0 bp
Trim 3' R2          : 0 bp
NextSeq Trim        : 0 bp
Fragment Size       : 200 bp
Fingerprint Bins    : 500000
Save Genome Index   : No
Max Resources       : 128 GB memory, 16 cpus, 10d time per job
Output Dir          : ./results
Launch Dir          : /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe
Working Dir         : /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/work
Script Dir          : /home/rad/users/gaurav/projects/workflows/nfchipseq
User                : rad
Config Profile      : standard
----------------------------------------------------
executor >  local (19)
executor >  local (19)
executor >  local (19)
executor >  local (20)
executor >  local (20)
executor >  local (112)
[2a/17bba7] process > CheckDesign (foxp1FkoPkf1oe_SampleSheet.csv)                  [100%] 1 of 1 ✔
[a3/acf60e] process > MakeTSSBED (genes.bed)                                        [100%] 1 of 1 ✔
[5a/be3a5d] process > MakeGenomeFilter (genome.fa)                                  [100%] 1 of 1 ✔
[78/a38d21] process > FastQC (foxp1ovko_200519062115_A0000504_FKO387Foxp1_mmu_ch... [100%] 8 of 8 ✔
[c3/6a826e] process > TrimGalore (foxp1ovko_200519062115_A0000504_FKO387Foxp1_mm... [100%] 8 of 8 ✔
[77/3b3831] process > BWAMem (foxp1ovko_200519062115_A0000504_FKO387Foxp1_mmu_ch... [100%] 8 of 8 ✔
[ef/4ccdce] process > SortBAM (foxp1ovko_200519062115_A0000504_FKO387Foxp1_mmu_c... [100%] 8 of 8 ✔
[c5/881019] process > MergeBAM (foxp1ovko_200519062115_A0000504_FKO387Foxp1_mmu_... [100%] 8 of 8 ✔
[09/334d42] process > MergeBAMFilter (foxp1ovko_200519062115_A0000504_FKO387Foxp... [100%] 8 of 8 ✔
[d6/65dcd5] process > Preseq (foxp1ovko_200519062115_A0000504_FKO387Foxp1_mmu_ch... [100%] 8 of 8 ✔
[6f/006d54] process > CollectMultipleMetrics (foxp1ovko_200519062115_A0000504_FK... [100%] 8 of 8 ✔
[a8/eb5378] process > BigWig (foxp1ovko_200519062115_A0000504_FKO387Foxp1_mmu_ch... [100%] 8 of 8 ✔
[c5/439bf0] process > PlotProfile (foxp1ovko_200519062115_A0000504_FKO387Foxp1_m... [100%] 8 of 8 ✔
[af/1229de] process > PhantomPeakQualTools (foxp1ovko_200519062115_A0000504_FKO3... [100%] 8 of 8 ✔
[96/b557b2] process > PlotFingerprint (foxp1ovko_200519062115_A0000504_PKF1OE062... [100%] 4 of 4 ✔
[04/936ead] process > MACSCallPeak (foxp1ovko_200519062115_A0000504_FKO387Foxp1_... [100%] 4 of 4 ✔
[0f/b2107f] process > AnnotatePeaks (foxp1ovko_200519062115_A0000504_FKO387Foxp1... [100%] 4 of 4 ✔
[44/15d884] process > PeakQC                                                        [100%] 1 of 1 ✔
[0d/79167a] process > ConsensusPeakSet (PKF1OE)                                     [100%] 2 of 2 ✔
[bf/abd18c] process > ConsensusPeakSetAnnotate (PKF1OE)                             [100%] 2 of 2 ✔
[-        ] process > ConsensusPeakSetDESeq                                         -
[3e/ca40b3] process > IGV                                                           [100%] 1 of 1 ✔
[52/1b1e8a] process > get_software_versions                                         [100%] 1 of 1 ✔
[9e/32282d] process > MultiQC                                                       [100%] 1 of 1 ✔
[9c/8d090f] process > output_documentation                                          [100%] 1 of 1 ✔
[0;35m[nf-core/chipseq] Pipeline completed successfully
Completed at: 21-May-2020 06:35:59
Duration    : 2h 51m 16s
CPU hours   : 51.4
Succeeded   : 112
