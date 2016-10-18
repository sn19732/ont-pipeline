

# Simulate a synthetic genome:

SIMULATED_GENOME=data/simulated/genome.fas
NR_CHROMS=23
MEAN_CHROM_LENGTH=5000000
CHROM_GAMMA_SHAPE=0.5
CHROM_BASE_FREQS=1,1,1,1

simulate_genome:  $(SIMULATED_GENOME)
$(SIMULATED_GENOME):
	@simulate_genome.py -n $(NR_CHROMS) -m $(MEAN_CHROM_LENGTH) -a $(CHROM_GAMMA_SHAPE) -b $(CHROM_BASE_FREQS) $(SIMULATED_GENOME)

# Download yeast genome from Ensembl:

YEAST_GENOME=data/yeast_genome.fas
fetch_genome: $(YEAST_GENOME)

$(YEAST_GENOME):
	@wget http://ftp.ensembl.org/pub/release-86/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
	@gzip -d Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz; mv Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa $(YEAST_GENOME)
# Simulate long reads:

SIMULATED_LONG_READS=data/simulated/long_reads.fastq
NR_SIM_LONG_READS=50000
LONG_READS_MEAN_LENGTH=6000
LONG_READS_GAMMA_SHAPE=0.5
LONG_READS_MIN_LENGTH=600
LONG_READS_ERROR_RATE=0.15
LONG_READS_ERROR_WEIGHTS=1,1,1

simulate_long_reads: $(SIMULATED_LONG_READS)
$(SIMULATED_LONG_READS): $(YEAST_GENOME)
	@simulate_sequencing_simple.py -n $(NR_SIM_LONG_READS) -m $(LONG_READS_MEAN_LENGTH) -a $(LONG_READS_GAMMA_SHAPE)\
		-l $(LONG_READS_MIN_LENGTH) -e $(LONG_READS_ERROR_RATE) -w $(LONG_READS_ERROR_WEIGHTS) $(YEAST_GENOME) $(SIMULATED_LONG_READS)


# Simulate short reads using simNGS:

RUNFILE=data/s_1_4x.runfile
COVERAGE=120.0
READ_LENGTH=101

simulate_short_reads: $(YEAST_GENOME)
	@(simLibrary -r $(READ_LENGTH) -x $(COVERAGE)  $(YEAST_GENOME)| simNGS -p paired $(RUNFILE) -O data/simulated/short_reads)

# Get rid of simulated data:
.PHONY: clean_simulated
clean_simulated:
	@rm $(SIMULATED_LONG_READS)
	@rm data/simulated/short_reads_end*

evaluate:
	@echo Running canu->racon->pilon pipeline.
	@make all
	@echo $(PILON_CONTIGS) 
	@echo Running canu->pilon pipeline.
	@make all
