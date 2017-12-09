Concatenate reads from different lanes (if necessary)

# concatenating reads from the same biological replicate that was sequences on two different runs.  Be sure to keep the F reads together and R reads together

# example: 

$ cat  <sample>.lane1.R1.fastq.gz  <sample>.lane2.R1.fastq.gz  >  <sample>_cat.R1.fastq.gz
$ cat  <sample>.lane1.R2.fastq.gz  <sample>.lane2.R2.fastq.gz  >  <sample>_cat.R2.fastq.gz
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Trimming adapters with NGmerge:

# Running NGmerge on raw sequence reads to trim adapters.  No quality score trimming 

# slurm script is 'sbatch_NGmerge_FINAL.sh 

NGmerge \
-a \
-1 <sample>.R1.fasta.gz \
-2 <sample>.R2.fasta.gz \
-o outputname
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Map reads with Bowtie2

# map trimmed reads (NGmerge output) with bowtie2 

# slurm script is 'sbatch_bowtie2_dovetail_FINAL.sh 

bowtie2 -x hm \
-X 2000 \
-1 <ngmergeoutput>.R1.fastq.gz \
-2 <ngmergeoutput>.R2.fastq.gz \
-p 31 | samtools view -b -S - | samtools sort - outfilename

samtools index outfilename.bam

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Remove mitochondrial genome with removeChrom 

# Remove Hofstenia MT genome using removeChrom program from Harvard RC #

# slurm script is 'sbatch_removeChrom_FINAL.sh' #

# scaffold2257 should be changed to whatever designates the MT genome in your input fasta file 

samtools view -h <inputfiltename>.bam | removeChrom - - scaffold2257 | samtools view -b - | samtools sort - -o outputname_noMt_sorted.bam

samtools index outputname_noMt_sorted.bam

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Remove duplicates with picard

# remove duplicates with Picard 

# slurm script is `sbatch_removeDups.sh' 

# the VALIDATION_STRINGENCY=LENIENT flag is because for some reason, when reads from separate runs have been concatenated together using cat, some of the headers are wrong or missing reads or something.  So I told picard to just ignore these and continue to run, but it would be good to know exactly what its ignoring 

java -jar $PICARD_TOOLS_HOME/picard.jar MarkDuplicates I=<inputfilname_noMT>.bam O=inputfilname_noMT_nodups.bam M=inputnoMTfiledups.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

samtools index inputfilname_noMT_nodups.bam

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Downsample bam file to desired sequwencing depth

# downsampling bam files to reach a similar number of matched reads 

# slurm script is `sbatch_downsample.sh' 

# P= percent of reads that you want to keep.  So in this case, you are telling picard to keep 41% of the reads.  So run flagstat on your final no MT and no DUPS .bam file to get the total # of mapped reads, then downsample using this number as 100% #

java -jar /n/sw/fasrcsw/apps/Core/picard/2.9.0-fasrc01/picard.jar DownsampleSam I=<inputbamfile>_noMT_noDups.bam O=<inputbamfile>_noMT_noDups.bam_down.bam VALIDATION_STRINGENCY=SILENT CREATE_INDEX=TRUE P=0.41

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Call peaks with MACS2

# Peak calling with MACS 

# slurm script is 'sbatch_MACS2_FINAL.sh' 

# There is some debate about the best way to run this program.  Empirically, I found that the parameters below provide a stringency that calls peaks that I was comfortable with (e.g. probably a big lenient).  The parameters were chosen here to aid with footprinting.  For simply calling peaks, it might be better to use other parameters 

# From ATAC literature:  To fine enriched cutting sites such as some DNase-seq datasets, all 5' ends of sequenced reasd should be extended in both directions to smooth the pileup signals.  If the wanted smoothing window is 200 bps, then use --nomodel --shift -100 --extsize 200 

# OmniATAC protocol uses only --nomodel --nolambda as inputs 

# To call nucleosomes, from ATAC literature: For certain nucleosome-seq data, we need to pileup the centers of nucleosomes using a half nucleosomze size for wavelet analysis (e.g. NPS algorith).  Since the DNA wrappted on nucleosome is about 147 bps, this option can be used: --nomodel --shift 37 --extsize 73 

 # keep dup should be all, because they were removed in the previous step 

macs2 callpeak \
 --treatment <inputfilename>_noMT_noDups_down.bam \
 -q 0.05  \
 --gsize  800000000 \
  --nomodel \
  --extsize 200 /
  --shift -100 \
  --nolambda  \
  -n outputfilename \
  --keep-dup all

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Filter peaks that fall withing repetitive regions

# Filter peaks that fall withing repetitive regions

# -v says "report in -a whever is not present in -b".

# -wa says "report the outcome values as they are in -a"

# for Hofstenia genome v1, the repeat file is "hmi_repeats.gff"

$ module load bedtools

$ bedtools intersect -wa -v -a <inputpeakfile.narrowPeak> -b <repeats.bed>  > peakfile_norepeats.narrowPeak

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
