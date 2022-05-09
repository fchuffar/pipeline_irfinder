rsync -auvP /Volumes/LACIE4TFCH_Storage/datashare/rnaseq_andre_ebi/ cargo:~/projects/datashare/rnaseq_andre_ebi/ --dry-run
ssh cargo
cd ~/projects/datashare/rnaseq_andre_ebi/raw
md5sum *.fastq.gz > md5sum.bettik.txt 
diff md5sum.txt md5sum.bettik.txt
rsync -auvP ~/projects/datashare/rnaseq_andre_ebi/ ~/projects/datashare_epistorage/rnaseq_andre_ebi/
cd ~/projects/datashare_epistorage/rnaseq_andre_ebi/raw/
md5sum *.fastq.gz > md5sum.summer.txt 
diff md5sum.txt md5sum.summer.txt




cd ~/projects/pipeline_irfinder/results/
rsync -auvP ~/projects/pipeline_irfinder/ cargo:~/projects/pipeline_irfinder/

# https://github.com/williamritchie/IRFinder

# wget https://github.com/williamritchie/IRFinder/archive/v1.3.0.tar.gz
# tar -zxvf v1.3.0.tar.gz

# chmod -R a+x IRFinder-1.3.0/bin
# cd IRFinder-1.3.0/

# IRFinder -m BuildRefFromSTARRef -r ReferenceDir -x STARRefDir


rm -Rf ReferenceDir
IRFinder -m BuildRefFromSTARRef -r ReferenceDir -x ~/projects/datashare/genomes/Homo_sapiens/UCSC/hg38/Sequence/StarIndex/


chuffarf@luke61:~/projects/pipeline_irfinder$ IRFinder -m BuildRefFromSTARRef -r ReferenceDir -x ~/projects/datashare/genomes/Homo_sapiens/UCSC/hg38/Seque
nce/StarIndex/
Argument error: -r ReferenceDir. Reference directory must not exist, BuildRef will create it.
chuffarf@luke61:~/projects/pipeline_irfinder$ rm -Rf ReferenceDir                                                                                         
chuffarf@luke61:~/projects/pipeline_irfinder$ IRFinder -m BuildRefFromSTARRef -r ReferenceDir -x ~/projects/datashare/genomes/Homo_sapiens/UCSC/hg38/Seque
nce/StarIndex/
Launching reference build process. The full build might take hours.
<Phase 1: STAR Reference Preparation>
Nov 19 21:16:44 ... copying the genome FASTA file...
cp: failed to close 'ReferenceDir/genome.fa': Input/output error
Nov 19 21:17:02 ... copying the transcriptome GTF file...
Nov 19 21:17:04 ... copying the STAR reference folder...
cp: failed to close 'ReferenceDir/STAR/Genome': Input/output error
<Phase 2: Mapability Calculation>
Nov 19 21:19:02 ... mapping genome fragments back to genome...
Nov 19 21:41:46 ... sorting aligned genome fragments...
[bam_sort_core] merging from 64 files and 32 in-memory blocks...
Nov 19 21:44:56 ... indexing aligned genome fragments...
Nov 19 21:45:16 ... filtering aligned genome fragments by chromosome/scaffold...
Nov 19 21:46:15 ... merging filtered genome fragments...
Nov 19 21:46:59 ... calculating regions for exclusion...
Nov 19 21:51:15 ... cleaning temporary files...
<Phase 3: IRFinder Reference Preparation>
Nov 19 21:51:17 ... building Ref 1...
sort: unknown subpragma '_mergesort' at /summer/epistorage/opt/IRFinder-1.3.0/bin/util/gtf2bed-custom.pl line 33.
BEGIN failed--compilation aborted at /summer/epistorage/opt/IRFinder-1.3.0/bin/util/gtf2bed-custom.pl line 33.
Nov 19 21:51:17 ... building Ref 2...
Nov 19 21:51:20 ... building Ref 3...
Nov 19 21:51:20 ... building Ref 4...
Nov 19 21:51:25 ... building Ref 5...
Nov 19 21:51:30 ... building Ref 6...
Nov 19 21:51:30 ... building Ref 7...
Nov 19 21:51:30 ... building Ref 8...
Nov 19 21:51:30 ... building Ref 9...
Nov 19 21:51:30 ... building Ref 10c...
Nov 19 21:51:30 ... building Ref 11c...
Error: exclude.directional.bed is empty.
Error: introns.unique.bed is empty.
Error: ref-cover.bed is empty.
Error: ref-read-continues.ref is empty.
Error: ref-sj.ref is empty.
Error: IRFinder reference building FAILED.






cd ~/projects/pipeline_irfinder

rsync -auvP cargo:~/projects/pipeline_irfinder/ ~/projects/pipeline_irfinder/ --dry-run --exclude="*.bam" --exclude="*.fa"
rsync -auvP cargo:~/projects/datashare/rnaseq_andre_ebi/ ~/projects/datashare/rnaseq_andre_ebi/ --include="*/" --include="*IRFinder-IR-dir.txt" --exclude="*" --dry-run
rsync -auvP cargo:~/projects/datashare/rnaseq_andre_ebi/ ~/projects/datashare/rnaseq_andre_ebi/ --include="*/" --include="*ResultsByIntron.txt" --exclude="*" --dry-run
rsync -auvP cargo:~/projects/datashare/rnaseq_andre_ebi/ ~/projects/datashare/rnaseq_andre_ebi/ --include="*/" --include="*SIRratio.txt" --exclude="*" --dry-run


 

rsync -auvP cargo:/summer/epistorage/opt/IRFinder-1.3.0/bin/DESeq2Constructor.R /opt/IRFinder-1.3.0/bin/DESeq2Constructor.R --dry-run



mkdir  /bettik/chuffarf/pipeline_irfinder
cd  /bettik/chuffarf/pipeline_irfinder
wget https://github.com/williamritchie/IRFinder/archive/v1.3.0.tar.gz
tar -zxvf v1.3.0.tar.gz

cd IRFinder-1.3.0/
rm -Rf REF/Human-GRCh38-release100
bin/IRFinder -m BuildRef -r REF/Human-GRCh38-release100 \
  -e REF/extra-input-files/RNA.SpikeIn.ERCC.fasta.gz \
  -b REF/extra-input-files/Human_hg38_nonPolyA_ROI.bed \
  ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz








# cd  /bettik/chuffarf/pipeline_irfinder/IRFinder-1.3.0
# bin/IRFinder -r REF/Human-GRCh38-release100 -d test_repro ~/projects/datashare/rnaseq_andre_ebi/raw/DB19_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB19_R2.fastq.gz


cd  /bettik/chuffarf/pipeline_irfinder
cd IRFinder-1.3.0/
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB19_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB19_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB19_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB20_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB20_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB20_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB21_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB21_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB21_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB22_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB22_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB22_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB23_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB23_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB23_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB24_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB24_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB24_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB25_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB25_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB25_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB26_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB26_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB26_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB27_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB27_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB27_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB28_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB28_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB28_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB29_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB29_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB29_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB30_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB30_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB30_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB31_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB31_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB31_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB32_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB32_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB32_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB33_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB33_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB33_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB34_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB34_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB34_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB35_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB35_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB35_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB36_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB36_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB36_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB37_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB37_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB37_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB38_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB38_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB38_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB39_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB39_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB39_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB40_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB40_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB40_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB41_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB41_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB41_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB42_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB42_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB42_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB43_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB43_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB43_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB44_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB44_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB44_R2.fastq.gz
bin/IRFinder -r REF/Human-GRCh38-release100 -d DB45_irfinder ~/projects/datashare/rnaseq_andre_ebi/raw/DB45_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB45_R2.fastq.gz


cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB19_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB19_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB22_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB20_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB25_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB21_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB20_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB22_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB23_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB23_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB26_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB24_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB19_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB25_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB22_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB26_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB25_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB27_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB20_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB28_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB23_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB29_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB26_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB30_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB19_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB31_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB22_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB32_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB25_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB33_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB20_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB34_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB23_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB35_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB26_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB36_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB19_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB37_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB22_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB38_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB25_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB39_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB20_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB40_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB23_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB41_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB26_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB42_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB19_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB43_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB22_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB44_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB25_irfinder/; samtools sort -@ 24 -T /dev/shm/tmp_DB45_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam





cd  /bettik/chuffarf/pipeline_irfinder
cd IRFinder-1.3.0/
ls   DB*_irfinder/*Log.final.out 
multiqc --force -o ~/projects/pipeline_irfinder/IRFinder-1.3.0/ -n multiqc_notrim \
  DB*_irfinder/*Log.final.out \
  DB*_irfinder/Unsorted.bam 

rsync -auvP cargo:~/projects/pipeline_irfinder/ ~/projects/pipeline_irfinder/ --dry-run --exclude="*.bam" --exclude="*.fa"



cd  /bettik/chuffarf/pipeline_irfinder
R



# cd  /bettik/chuffarf/pipeline_irfinder
# cd IRFinder-1.3.0/
# bin/AdaptorDetect.pl ~/projects/datashare/rnaseq_andre_ebi/raw/DB19_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB19_R2.fastq.gz
# bin/AdaptorDetect.pl ~/projects/datashare/rnaseq_andre_ebi/raw/DB20_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB20_R2.fastq.gz
# bin/AdaptorDetect.pl ~/projects/datashare/rnaseq_andre_ebi/raw/DB22_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB22_R2.fastq.gz
# bin/AdaptorDetect.pl ~/projects/datashare/rnaseq_andre_ebi/raw/DB23_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB23_R2.fastq.gz
# bin/AdaptorDetect.pl ~/projects/datashare/rnaseq_andre_ebi/raw/DB25_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB25_R2.fastq.gz
# bin/AdaptorDetect.pl ~/projects/datashare/rnaseq_andre_ebi/raw/DB26_R1.fastq.gz ~/projects/datashare/rnaseq_andre_ebi/raw/DB26_R2.fastq.gz




# rm -Rf pooled_NHS pooled_HS
# bin/IRFinder -m BAM  -r REF/Human-GRCh38-release100 -d pooled_NHS <(samtools cat DB19_irfinder/Unsorted.bam DB22_irfinder/Unsorted.bam)
# bin/IRFinder -m BAM  -r REF/Human-GRCh38-release100 -d pooled_HS <(samtools cat DB20_irfinder/Unsorted.bam DB23_irfinder/Unsorted.bam)
#
# bin/analysisWithLowReplicates.pl \
#   -A pooled_NHS/IRFinder-IR-dir.txt DB19_irfinder/IRFinder-IR-dir.txt DB22_irfinder/IRFinder-IR-dir.txt \
#   -B pooled_HS/IRFinder-IR-dir.txt DB20_irfinder/IRFinder-IR-dir.txt DB23_irfinder/IRFinder-IR-dir.txt  \
#   > NHS-vs-HS.txt
#





cd  /bettik/chuffarf/pipeline_irfinder
cd IRFinder-1.3.0/
echo DB19_irfinder/IRFinder-IR-dir.txt >  filePaths.txt
echo DB22_irfinder/IRFinder-IR-dir.txt >> filePaths.txt
echo DB25_irfinder/IRFinder-IR-dir.txt >> filePaths.txt
echo DB20_irfinder/IRFinder-IR-dir.txt >> filePaths.txt
echo DB23_irfinder/IRFinder-IR-dir.txt >> filePaths.txt
echo DB26_irfinder/IRFinder-IR-dir.txt >> filePaths.txt
cat filePaths.txt


echo -e "SampleNames\tCondition" > experiment.txt
echo -e "DB19_irfinder\tNHS" >>    experiment.txt
echo -e "DB22_irfinder\tNHS" >>    experiment.txt
echo -e "DB25_irfinder\tNHS" >>    experiment.txt
echo -e "DB20_irfinder\tHS"  >>    experiment.txt
echo -e "DB23_irfinder\tHS"  >>    experiment.txt
echo -e "DB26_irfinder\tHS"  >>    experiment.txt
cat experiment.txt  

  





cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB19_irfinder/; samtools sort -@ 32 -T /dev/shm/tmp_DB19_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB22_irfinder/; samtools sort -@ 32 -T /dev/shm/tmp_DB22_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB25_irfinder/; samtools sort -@ 32 -T /dev/shm/tmp_DB25_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB20_irfinder/; samtools sort -@ 32 -T /dev/shm/tmp_DB20_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB23_irfinder/; samtools sort -@ 32 -T /dev/shm/tmp_DB23_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam
cd ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB26_irfinder/; samtools sort -@ 32 -T /dev/shm/tmp_DB26_sortbam Unsorted.bam -o Sorted.bam ; samtools index Sorted.bam

rsync -auvP cargo:~/projects/pipeline_irfinder/IRFinder-1.3.0/DB19_irfinder/Sorted.bam ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB19_irfinder/
rsync -auvP cargo:~/projects/pipeline_irfinder/IRFinder-1.3.0/DB22_irfinder/Sorted.bam ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB22_irfinder/
rsync -auvP cargo:~/projects/pipeline_irfinder/IRFinder-1.3.0/DB25_irfinder/Sorted.bam ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB25_irfinder/
rsync -auvP cargo:~/projects/pipeline_irfinder/IRFinder-1.3.0/DB20_irfinder/Sorted.bam ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB20_irfinder/
rsync -auvP cargo:~/projects/pipeline_irfinder/IRFinder-1.3.0/DB23_irfinder/Sorted.bam ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB23_irfinder/
rsync -auvP cargo:~/projects/pipeline_irfinder/IRFinder-1.3.0/DB26_irfinder/Sorted.bam ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB26_irfinder/
rsync -auvP cargo:~/projects/pipeline_irfinder/IRFinder-1.3.0/DB19_irfinder/Sorted.bam.bai ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB19_irfinder/
rsync -auvP cargo:~/projects/pipeline_irfinder/IRFinder-1.3.0/DB22_irfinder/Sorted.bam.bai ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB22_irfinder/
rsync -auvP cargo:~/projects/pipeline_irfinder/IRFinder-1.3.0/DB25_irfinder/Sorted.bam.bai ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB25_irfinder/
rsync -auvP cargo:~/projects/pipeline_irfinder/IRFinder-1.3.0/DB20_irfinder/Sorted.bam.bai ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB20_irfinder/
rsync -auvP cargo:~/projects/pipeline_irfinder/IRFinder-1.3.0/DB23_irfinder/Sorted.bam.bai ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB23_irfinder/
rsync -auvP cargo:~/projects/pipeline_irfinder/IRFinder-1.3.0/DB26_irfinder/Sorted.bam.bai ~/projects/pipeline_irfinder/IRFinder-1.3.0/DB26_irfinder/





diff volcano_in_WT_NHS_ref_vs_HS_blue_genes.txt ~/projects/pipeline_irfinder.old/results/volcano_WT_NHS_HS_blue_genes.txt 
diff volcano_in_WT_NHS_ref_vs_R3_blue_genes.txt ~/projects/pipeline_irfinder.old/results/volcano_WT_NHS_R3_blue_genes.txt 
diff volcano_in_C1_NHS_ref_vs_HS_blue_genes.txt ~/projects/pipeline_irfinder.old/results/volcano_C1_NHS_HS_blue_genes.txt 
diff volcano_in_C1_NHS_ref_vs_R3_blue_genes.txt ~/projects/pipeline_irfinder.old/results/volcano_C1_NHS_R3_blue_genes.txt 
diff volcano_in_M3_NHS_ref_vs_HS_blue_genes.txt ~/projects/pipeline_irfinder.old/results/volcano_M3_NHS_HS_blue_genes.txt 
diff volcano_in_M3_NHS_ref_vs_R3_blue_genes.txt ~/projects/pipeline_irfinder.old/results/volcano_M3_NHS_R3_blue_genes.txt 
diff volcano_in_WT_NHS_ref_vs_HS_red_genes.txt ~/projects/pipeline_irfinder.old/results/volcano_WT_NHS_HS_red_genes.txt 
diff volcano_in_WT_NHS_ref_vs_R3_red_genes.txt ~/projects/pipeline_irfinder.old/results/volcano_WT_NHS_R3_red_genes.txt 
diff volcano_in_C1_NHS_ref_vs_HS_red_genes.txt ~/projects/pipeline_irfinder.old/results/volcano_C1_NHS_HS_red_genes.txt 
diff volcano_in_C1_NHS_ref_vs_R3_red_genes.txt ~/projects/pipeline_irfinder.old/results/volcano_C1_NHS_R3_red_genes.txt 
diff volcano_in_M3_NHS_ref_vs_HS_red_genes.txt ~/projects/pipeline_irfinder.old/results/volcano_M3_NHS_HS_red_genes.txt 
diff volcano_in_M3_NHS_ref_vs_R3_red_genes.txt ~/projects/pipeline_irfinder.old/results/volcano_M3_NHS_R3_red_genes.txt 


diff volcano_in_WT_NHS_ref_vs_HS_blue_ir.xlsx ~/projects/pipeline_irfinder.old/results/volcano_WT_NHS_HS_blue_ir.xlsx 
diff volcano_in_WT_NHS_ref_vs_R3_blue_ir.xlsx ~/projects/pipeline_irfinder.old/results/volcano_WT_NHS_R3_blue_ir.xlsx 
diff volcano_in_C1_NHS_ref_vs_HS_blue_ir.xlsx ~/projects/pipeline_irfinder.old/results/volcano_C1_NHS_HS_blue_ir.xlsx 
diff volcano_in_C1_NHS_ref_vs_R3_blue_ir.xlsx ~/projects/pipeline_irfinder.old/results/volcano_C1_NHS_R3_blue_ir.xlsx 
diff volcano_in_M3_NHS_ref_vs_HS_blue_ir.xlsx ~/projects/pipeline_irfinder.old/results/volcano_M3_NHS_HS_blue_ir.xlsx 
diff volcano_in_M3_NHS_ref_vs_R3_blue_ir.xlsx ~/projects/pipeline_irfinder.old/results/volcano_M3_NHS_R3_blue_ir.xlsx 
diff volcano_in_WT_NHS_ref_vs_HS_red_ir.xlsx ~/projects/pipeline_irfinder.old/results/volcano_WT_NHS_HS_red_ir.xlsx 
diff volcano_in_WT_NHS_ref_vs_R3_red_ir.xlsx ~/projects/pipeline_irfinder.old/results/volcano_WT_NHS_R3_red_ir.xlsx 
diff volcano_in_C1_NHS_ref_vs_HS_red_ir.xlsx ~/projects/pipeline_irfinder.old/results/volcano_C1_NHS_HS_red_ir.xlsx 
diff volcano_in_C1_NHS_ref_vs_R3_red_ir.xlsx ~/projects/pipeline_irfinder.old/results/volcano_C1_NHS_R3_red_ir.xlsx 
diff volcano_in_M3_NHS_ref_vs_HS_red_ir.xlsx ~/projects/pipeline_irfinder.old/results/volcano_M3_NHS_HS_red_ir.xlsx 
diff volcano_in_M3_NHS_ref_vs_R3_red_ir.xlsx ~/projects/pipeline_irfinder.old/results/volcano_M3_NHS_R3_red_ir.xlsx 


in_WT_NHS_ref_vs_HS
in_WT_NHS_ref_vs_R3
in_C1_NHS_ref_vs_HS
in_C1_NHS_ref_vs_R3
in_M3_NHS_ref_vs_HS
in_M3_NHS_ref_vs_R3



