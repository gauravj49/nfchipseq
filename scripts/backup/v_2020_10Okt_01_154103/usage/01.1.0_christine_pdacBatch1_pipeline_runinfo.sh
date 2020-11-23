[rad@imows4 pdacBatch1 ] $ nextflow run /home/rad/users/gaurav/projects/workflows/nfchipseq --input /media/rad/HDD1/nfchip/christine/pdacBatch1/pdacBatch1_SampleSheet.csv --genome GRCm38 --narrow_peak --single_end -name pdacBatch1
N E X T F L O W  ~  version 20.04.1
Launching `/home/rad/users/gaurav/projects/workflows/nfchipseq/main.nf` [pdacBatch1] - revision: f362665a84
----------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/chipseq v1.1.0
----------------------------------------------------
----------------------------------------------------
Run Name            : pdacBatch1
Data Type           : Single-End
Design File         : /media/rad/HDD1/nfchip/christine/pdacBatch1/pdacBatch1_SampleSheet.csv
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
Launch Dir          : /media/rad/HDD1/nfchip/christine/pdacBatch1
Working Dir         : /media/rad/HDD1/nfchip/christine/pdacBatch1/work
Script Dir          : /home/rad/users/gaurav/projects/workflows/nfchipseq
User                : rad
Config Profile      : standard
----------------------------------------------------
executor >  local (168)
[28/16cc1c] process > CheckDesign (pdacBatch1_SampleSheet.csv)                                                                                                                                    [100%] 1 of 1 ✔
[98/ed53c2] process > MakeTSSBED (genes.bed)                                                                                                                                                      [100%] 1 of 1 ✔
[83/2aef9d] process > MakeGenomeFilter (genome.fa)                                                                                                                                                [100%] 1 of 1 ✔
[44/f6f60d] process > FastQC (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT1-H3K27ac_mmu_chipseq_se_R1_R1_T1)                                                                             [100%] 12 of 12 ✔
[e5/2c767a] process > TrimGalore (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT1-H3K27ac_mmu_chipseq_se_R1_R1_T1)                                                                         [100%] 12 of 12 ✔
[d6/7b1579] process > BWAMem (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT1-H3K27ac_mmu_chipseq_se_R1_R1_T1)                                                                             [100%] 12 of 12 ✔
[8a/a8ed91] process > SortBAM (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT1-H3K27ac_mmu_chipseq_se_R1_R1_T1)                                                                            [100%] 12 of 12 ✔
[64/935a6c] process > MergeBAM (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT1-H3K27ac_mmu_chipseq_se_R1_R1)                                                                              [100%] 12 of 12 ✔
[54/bffc44] process > MergeBAMFilter (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT1-H3K27ac_mmu_chipseq_se_R1_R1)                                                                        [100%] 12 of 12 ✔
[db/a44dd9] process > Preseq (pdacBatch1_20200424121526_A0000498_127abcam-5320-PPT1-H3K27ac_mmu_chipseq_se_R1_R1)                                                                                 [100%] 12 of 12 ✔
[37/fe5d85] process > CollectMultipleMetrics (pdacBatch1_20200424121526_A0000498_127abcam-5320-PPT1-H3K27ac_mmu_chipseq_se_R1_R1)                                                                 [100%] 12 of 12 ✔
[ae/d6a95b] process > BigWig (pdacBatch1_20200424121526_A0000498_127abcam-5320-PPT1-H3K27ac_mmu_chipseq_se_R1_R1)                                                                                 [100%] 12 of 12 ✔
[ba/08babf] process > PlotProfile (pdacBatch1_20200424121526_A0000498_125abcam-53646-LivMet2-H3K27ac_mmu_chipseq_se_R1_R1)                                                                        [100%] 12 of 12 ✔
[13/191240] process > PhantomPeakQualTools (pdacBatch1_20200424121526_A0000498_127abcam-5320-PPT1-H3K27ac_mmu_chipseq_se_R1_R1)                                                                   [100%] 12 of 12 ✔
[0d/e18159] process > PlotFingerprint (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT1-H3K27ac_mmu_chipseq_se_R1_R1 vs pdacBatch1_20200424121526_A0000498_input53646_mmu_chipseq_se_R1_R1) [100%] 8 of 8 ✔
[5a/2364ef] process > MACSCallPeak (pdacBatch1_20200424121526_A0000498_127abcam-5320-PPT1-H3K27ac_mmu_chipseq_se_R1_R1 vs pdacBatch1_20200424121526_A0000498_input5320_mmu_chipseq_se_R1_R1)      [100%] 8 of 8 ✔
[54/112a6a] process > AnnotatePeaks (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT1-H3K27ac_mmu_chipseq_se_R1_R1 vs pdacBatch1_20200424121526_A0000498_input53646_mmu_chipseq_se_R1_R1)   [100%] 8 of 8 ✔
[3b/259e96] process > PeakQC                                                                                                                                                                      [100%] 1 of 1 ✔
[ad/3fd048] process > ConsensusPeakSet (53646H3k27ac)                                                                                                                                             [100%] 2 of 2 ✔
[13/f7f75a] process > ConsensusPeakSetAnnotate (53646H3k27ac)                                                                                                                                     [100%] 2 of 2 ✔
[-        ] process > ConsensusPeakSetDESeq                                                                                                                                                       -
[9c/ae8737] process > IGV                                                                                                                                                                         [100%] 1 of 1 ✔
[39/dd4a6a] process > get_software_versions                                                                                                                                                       [100%] 1 of 1 ✔
[3d/cb6a8c] process > MultiQC                                                                                                                                                                     [100%] 1 of 1 ✔
[17/08fdbb] process > output_documentation                                                                                                                                                        [100%] 1 of 1 ✔
[0;35m[nf-core/chipseq] Pipeline completed successfully
Completed at: 22-May-2020 05:49:56
Duration    : 2h 28m 46s
CPU hours   : 82.1
Succeeded   : 168
