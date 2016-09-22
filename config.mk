# General pipeline parameters:

# Parent directory to pipeline workspace:
WORKSPACE=./
# Pipeline name:
PIPELINE_NAME=assembly-polish
# Pipeline working directory:
WDIR=$(WORKSPACE)/$(PIPELINE_NAME)
# Results directory:
RES=$(WDIR)/results
# Pipeline git repo:
REPO=git@git.oxfordnanolabs.local:bsipos/ont-assembly-polish.git

# Custom pipeline parameters:

# Input files:

NANOPORE_READS=
ILLUMINA_READS_PAIR1=
ILLUMINA_READS_PAIR2=


# Number of cores to use for multithreaded applications:
CORES=32
