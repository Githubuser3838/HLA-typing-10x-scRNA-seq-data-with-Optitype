# HLA-typing-10x-scRNA-seq-data-with-Optitype
 I wanted to share the procedure I used to successfully type a few samples (around 50) sequenced with 10x Genomics, using Optitype (https://github.com/FRED-2/OptiType). 
 The results from Optitype were verified independently by conventional HLA typing.

The strategy is to filter the bam file from Cell Ranger for chromosome 6, and then, optionally, sub-sample the resulting file (if there are too many reads at the typing step and it crashes). We then make fastq files using 10x bamtofastq, and type with Optitype. Reads from the different cells are collapsed onto one file. We lose single cell information, but because cells come from the same sample, this allows accurate typing: there are collectively many reads to do the typing despite relatively few reads coming from one particular cell.

Sub-sampling may be necessary for single-cell 5’ assays, because the assay selects many reads from the first exons. For 3' assays, sub-sampling should not be used, because the assay doesn't select many reads from the first exons (comparatively), and Optitype relies mostly on exons 2 and 3 to differentiate HLA alleles.
The parameter -s for subsampling can be adjusted, but I got fast and consistent results with 0.001 (after sequencing 5000-10000 cells with 10x v2 chemistry, 5' assay). This depends on the dataset, it should be adjusted so that enough reads are kept after alignment to the reference, but not so many that it crashes at the typing step (typing worked well with about 1000-5000 reads after alignment). Subsampling can be repeated with a different seed to make sure typing is consistent.

I used a conda installation. I also used bwa mem for mapping onto the reference, as described by @armish at hammerlab/biokepi#419, since it solves memory issues that arise with razers3.

Here is the procedure:

### Install Optitype and index the reference:

conda create -n optitype-conda\
conda activate optitype-conda\
conda install -c bioconda optitype\
conda install -c bioconda bwa\
conda install -c bioconda 10x_bamtofastq

REFRNA=\~/miniconda3/envs/optitype-conda/bin/data/hla_reference_rna.fasta\
OPTITYPE_HOME=\~/miniconda3/envs/optitype-conda/bin

bwa index $REFRNA

### Make chr6-filtered bam:

samtools view -t8 -bh sample_alignments.bam chr6 -o chr6.filtered.bam\
or (depending on your bam file):\
samtools view -t8 -bh sample_alignments.bam 6 -o chr6.filtered.bam

### Sub-sample (optional):

samtools view -b -s 0.001 chr6.filtered.bam -o chr6.filtered.sub.bam

(To try a different sub-sampled set of reads, change the seed in the integer part of the -s parameter, e.g. 123.001)

### Make fastq from R2 reads:

bamtofastq chr6.filtered.sub.bam ./bamtofastq\
cd bamtofastq/*/\
cat *_R2_001.fastq.gz > Sample_R2_001.fastq.gz\
(R1 reads correspond to barcodes/UMI)

### Run alignement and Optitype on generated fastq:

bwa mem $REFRNA Sample_R2_001.fastq.gz | samtools fastq -F4 > fished.Sample_R2_001.fastq\
(should yield 1000-5000 reads for the typing step, although accurate results can be obtained with much less)

python $OPTITYPE_HOME/OptiTypePipeline.py -i fished.Sample_R2_001.fastq --rna -v -o ./
