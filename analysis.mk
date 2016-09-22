
all: 

# Pipeline targets:

# Nanopore assembly by canu:

CANU_PREFIX=canu
CANU_DIR=$(WDIR)/canu-assembly
CANU_CONTIGS=$(CANU_DIR)/$(CANU_PREFIX).contigs.fasta

canu_assembly: $(CANU_CONTIGS)
$(CANU_CONTIGS): $(NANOPORE_READS)
	@echo Assembling nanopore reads using canu.
	@canu\
		 -p $(CANU_PREFIX) \
		 -d $(CANU_DIR) \
		 genomeSize=$(CANU_GENOME_SIZE) \
		 -nanopore-raw $(NANOPORE_READS)

# Index contigs, map Illumina reads to contigs by BWA, sorting and indexing using samtools:

BWA_BAM_PREFIX=$(WDIR)/bwa_aligned_reads
BWA_BAM=$(BWA_BAM_PREFIX).bam

bwa_align: $(BWA_BAM)
$(BWA_BAM): $(CANU_CONTIGS) $(ILLUMINA_READS_PAIR1) $(ILLUMINA_READS_PAIR2)
	@echo Indexing contigs.
	@bwa index $(CANU_CONTIGS)
	@echo Aligning Illumina reads using BWA mem.
	@bwa mem -t $(CORES) $(BWA_PARAMETERS) $(CANU_CONTIGS)  $(ILLUMINA_READS_PAIR1) $(ILLUMINA_READS_PAIR2)\
		| samtools sort -@ $(CORES) - $(BWA_BAM_PREFIX)
