
all: 

# Pipeline targets:

# Nanopore assembly by canu:

CANU_PREFIX=canu
CANU_DIR=$(WDIR)/canu-assembly
CANU_CONTIGS=$(CANU_DIR)/$(CANU_PREFIX).contigs.fasta

canu-assembly: $(CANU_CONTIGS)
$(CANU_CONTIGS): $(NANOPORE_READS)
	@echo Assembling nanopore reads using canu.
	@canu\
		 -p $(CANU_PREFIX) \
		 -d $(CANU_DIR) \
		 genomeSize=$(CANU_GENOME_SIZE) \
		 -nanopore-raw $(NANOPORE_READS)
