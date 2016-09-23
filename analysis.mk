
# Pipeline targets:

# Nanopore assembly by canu:

CANU_PREFIX=canu
CANU_DIR=$(WDIR)/canu-assembly
CANU_CONTIGS=$(CANU_DIR)/$(CANU_PREFIX).contigs.fasta

canu_assembly: $(CANU_CONTIGS) $(WDT)
$(CANU_CONTIGS): $(NANOPORE_READS)
	@echo Assembling nanopore reads using canu.
	@canu\
		 -p $(CANU_PREFIX) \
		 -d $(CANU_DIR) \
		 genomeSize=$(CANU_GENOME_SIZE) \
		 -nanopore-raw $(NANOPORE_READS)

# Map nanopore reads to canu contings and use racon to perform correction based on nanopore reads only:

ifeq ($(USE_RACON),yes)
RACON_CONTIGS=$(WDIR)/racon.contigs.fasta
else
RACON_CONTIGS=$(CANU_CONTIGS)
endif

MINIMAP_OVERLAPS=$(WDIR)/minimap_overlaps.paf

racon_correct: $(RACON_CONTIGS)
$(RACON_CONTIGS): $(NANOPORE_READS) $(CANU_CONTIGS)
ifeq ($(USE_RACON),yes)
	@echo Mapping nanopore reads onto canu contings using minimap.
	@minimap $(CANU_CONTIGS) $(NANOPORE_READS) > $(MINIMAP_OVERLAPS)
	@echo Correcting contigs using racon.
	@racon $(NANOPORE_READS) $(MINIMAP_OVERLAPS) $(CANU_CONTIGS) $(RACON_CONTIGS)
else
	@echo Skipping racon polishing.
endif

# Index contigs, map Illumina reads to contigs by BWA, sorting and indexing using samtools:

BWA_BAM_PREFIX=$(WDIR)/bwa_aligned_reads
BWA_BAM=$(BWA_BAM_PREFIX).bam

bwa_align: $(BWA_BAM)
$(BWA_BAM): $(RACON_CONTIGS) $(ILLUMINA_READS_PAIR1) $(ILLUMINA_READS_PAIR2)
	@echo Indexing contigs.
	@bwa index $(RACON_CONTIGS)
	@echo Aligning Illumina reads using BWA mem.
	@bwa mem -t $(CORES) $(BWA_PARAMETERS) $(RACON_CONTIGS)  $(ILLUMINA_READS_PAIR1) $(ILLUMINA_READS_PAIR2)\
		| samtools view -S -b -u - | samtools sort -@ $(CORES) - $(BWA_BAM_PREFIX)
	@samtools index $(BWA_BAM)

# Correct contigs using pilon based on the Illumina reads:

PILON_CONTIGS=$(WDIR)/pilon.contigs.fasta

pilon_correct: $(PILON_CONTIGS)
$(PILON_CONTIGS): $(RACON_CONTIGS) $(BWA_BAM)
	@echo Correcting contigs using pilon.
	@pilon --threads $(CORES) --genome $(RACON_CONTIGS) --bam $(BWA_BAM) --outdir $(WDIR) --output pilon.contigs $(PILON_PARAMETERS)

all: $(PILON_CONTIGS)
	@echo
	@echo Analysis finishes.
	@echo Contigs assembled by Canu are at: $(CANU_CONTIGS)
	@echo Contigs polished by racon are at: $(RACON_CONTIGS)
	@echo Pilon corrected contigs are at: $(PILON_CONTIGS)

