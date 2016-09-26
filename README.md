ONT assembly and Illumina polishing pipeline
=============================================

This pipeline performs the following steps:
- Assembly of nanopore reads using [Canu](http://canu.readthedocs.io).
- Polish canu contigs using [racon](https://github.com/isovic/racon) (*optional*).
- Map a paired-end Illumina dataset onto the contigs obtained in the previous steps using [BWA](http://bio-bwa.sourceforge.net) mem.
- Perform correction of contigs using [pilon](https://github.com/broadinstitute/pilon/wiki) and the Illumina dataset.

Usage
-----

Edit `config.mk` to set input files and parameters. Specifying the following is mandatory:
- `NANOPORE_READS` - input nanopore reads (note that this **must** be a single valid fastq file, see [here](https://www.biostars.org/p/81924/) how to combine fastq files).
- `ILLUMINA_READS_PAIR1` - fastq with the first reads of the paired-end Illumina dataset.
- `ILLUMINA_READS_PAIR2` - fastq with the second reads of the paired-end Illumina dataset.
- `CANU_GENOME_SIZE` - genome size parameter passed to canu.

The number of cores used can be specified by `CORES` (set this to the number of CPUs in your machine).
Racon corrections can be disabled by setting `USE_RACON=no`.

Then issue issue `make all` to run the pipeline. Issue `make help` for a list of utility make targets. Issue `make clear_wdir` to delete the working directory (including all results!).

Using through docker
--------------------

The easiest way to use the pipeline is through docker. First pull the container:

```bash
docker pull docker-registry.oxfordnanolabs.local:5000/ont-assembly-polish:latest
```

Then run the container:

```bash
docker run -v /path/to/my_data:/data -it ont-assembly-polish:latest
```

You will be dropped into the directory `/home/ont-assembly-polish`, then simply edit `config.mk` and run the pipeline.
The `-v` flag will make the `/path/to/my_data` directory on the host available under `/data` in the container.

Application dependencies
------------------------

- [Canu](http://canu.readthedocs.io)
- samtools
- [BWA](http://bio-bwa.sourceforge.net)
- [racon](https://github.com/isovic/racon) - the pipeline will download and build it
- [minimap](https://github.com/lh3/minimap) and [miniasm](https://github.com/lh3/miniasm)
- [pilon](https://github.com/broadinstitute/pilon/wiki) - the pipeline will download it

TODO
----
- Add the pipeline itself to the docker image.
- Build simulation tool for evaluating the performance of correction.
