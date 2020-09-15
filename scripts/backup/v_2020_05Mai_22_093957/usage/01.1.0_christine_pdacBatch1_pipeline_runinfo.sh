[rad@imows4 pdacBatch1 ] $ nextflow run /home/rad/users/gaurav/projects/nfPipelines/nfchipseq --input /media/rad/HDD1/nfchip/pdacBatch1/nfpdacBatch1_SampleSheet.csv --genome GRCm38 --single_end --min_reps_consensus 1 --narrow_peak --save_macs_pileup -name nfpdacBatch1
N E X T F L O W  ~  version 20.04.1
Launching `/home/rad/users/gaurav/projects/nfPipelines/nfchipseq/main.nf` [nfpdacBatch1] - revision: f362665a84
----------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/chipseq v1.1.0
----------------------------------------------------
Run Name            : nfpdacBatch1
Data Type           : Single-End
Design File         : /media/rad/HDD1/nfchip/pdacBatch1/nfpdacBatch1_SampleSheet.csv
Genome              : GRCm38
Design File         : /media/rad/HDD1/nfchip/pdacBatch1/nfpdacBatch1_SampleSheet.csv
Genome              : GRCm38
Fasta File          : s3://ngi-igenomes/igenomes//Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa
GTF File            : s3://ngi-igenomes/igenomes//Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf
Gene BED File       : s3://ngi-igenomes/igenomes//Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.bed
BWA Index           : s3://ngi-igenomes/igenomes//Mus_musculus/Ensembl/GRCm38/Sequence/BWAIndex/genome.fa
Blacklist BED       : /home/rad/users/gaurav/projects/nfPipelines/nfchipseq/assets/blacklists/GRCm38-blacklist.bed
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
Save MACS2 Pileup   : Yes
Max Resources       : 128 GB memory, 16 cpus, 10d time per job
Output Dir          : ./results
Launch Dir          : /media/rad/HDD1/nfchip/pdacBatch1
Working Dir         : /media/rad/HDD1/nfchip/pdacBatch1/work
Script Dir          : /home/rad/users/gaurav/projects/nfPipelines/nfchipseq
User                : rad
Config Profile      : standard
----------------------------------------------------
executor >  local (154)
[6c/7e3a28] process > CheckDesign (nfpdacBatch1_SampleSheet.csv)                                                                                                                                     [100%] 1 of 1 ✔
[93/a41332] process > MakeTSSBED (genes.bed)                                                                                                                                                         [100%] 1 of 1 ✔
[ef/0f8a7c] process > MakeGenomeFilter (genome.fa)                                                                                                                                                   [100%] 1 of 1 ✔
[b5/5b1a3c] process > FastQC (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT-1_mmu_chipseq_se_R1_R1_T1)                                                                                       [100%] 12 of 12 ✔
[b1/8f3608] process > TrimGalore (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT-1_mmu_chipseq_se_R1_R1_T1)                                                                                   [100%] 12 of 12 ✔
[8a/2464e4] process > BWAMem (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT-1_mmu_chipseq_se_R1_R1_T1)                                                                                       [100%] 12 of 12 ✔
[0f/cd1b97] process > SortBAM (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT-1_mmu_chipseq_se_R1_R1_T1)                                                                                      [100%] 12 of 12 ✔
[a6/73a745] process > MergeBAM (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT-1_mmu_chipseq_se_R1_R1)                                                                                        [100%] 12 of 12 ✔
[4a/211cd9] process > MergeBAMFilter (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT-1_mmu_chipseq_se_R1_R1)                                                                                  [100%] 12 of 12 ✔
[ba/86776d] process > Preseq (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT-1_mmu_chipseq_se_R1_R1)                                                                                          [100%] 12 of 12 ✔
[4a/81b82c] process > CollectMultipleMetrics (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT-1_mmu_chipseq_se_R1_R1)                                                                          [100%] 12 of 12 ✔
[65/e67169] process > BigWig (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT-1_mmu_chipseq_se_R1_R1)                                                                                          [100%] 12 of 12 ✔
[76/729af6] process > PlotProfile (pdacBatch1_20200424121526_A0000498_125input-53546-LivMet-2_mmu_chipseq_se_R1_R1)                                                                                  [100%] 12 of 12 ✔
[13/bdb648] process > PhantomPeakQualTools (pdacBatch1_20200424121526_A0000498_127abcam-5320-PPT-1_mmu_chipseq_se_R1_R1)                                                                             [100%] 12 of 12 ✔
[24/88bf7d] process > PlotFingerprint (pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT-1_mmu_chipseq_se_R1_R1 vs pdacBatch1_20200424121526_A0000498_123input-53646-PPT-1_mmu_chipseq_se_R1_R1) [100%] 4 of 4 ✔
[9e/eebadc] process > MACSCallPeak (pdacBatch1_20200424121526_A0000498_127abcam-5320-PPT-1_mmu_chipseq_se_R1_R1 vs pdacBatch1_20200424121526_A0000498_127input-5320-PPT-1_mmu_chipseq_se_R1_R1)      [100%] 4 of 4 ✔
[c9/0a9bf8] process > AnnotatePeaks (pdacBatch1_20200424121526_A0000498_127abcam-5320-PPT-1_mmu_chipseq_se_R1_R1 vs pdacBatch1_20200424121526_A0000498_127input-5320-PPT-1_mmu_chipseq_se_R1_R1)     [100%] 4 of 4 ✔
[ec/8d7667] process > PeakQC                                                                                                                                                                         [100%] 1 of 1 ✔
[54/2981b7] process > ConsensusPeakSet (abcam)                                                                                                                                                       [100%] 1 of 1 ✔
[6d/80d932] process > ConsensusPeakSetAnnotate (abcam)                                                                                                                                               [100%] 1 of 1 ✔
[-        ] process > ConsensusPeakSetDESeq                                                                                                                                                          -
[51/60d9c9] process > IGV                                                                                                                                                                            [100%] 1 of 1 ✔
[1f/09df36] process > get_software_versions                                                                                                                                                          [100%] 1 of 1 ✔
[26/c04635] process > MultiQC                                                                                                                                                                        [100%] 1 of 1 ✔
[53/b687eb] process > output_documentation                                                                                                                                                           [100%] 1 of 1 ✔
[0;35m[nf-core/chipseq] Pipeline completed successfully
Completed at: 12-May-2020 03:18:09
Duration    : 2h 35m 14s
CPU hours   : 68.8
Succeeded   : 154