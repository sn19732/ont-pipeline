ONT assembly and Illumina polishing pipeline
=============================================

This pipeline performs the following steps:
- Assembly of nanopore reads using [Canu](http://canu.readthedocs.io).
- OPTIONAL: polish canu contigs using [racon](https://github.com/isovic/racon).
- Map a paired-end Illumina dataset onto the contigs obtained in the previous steps using [BWA](http://bio-bwa.sourceforge.net) mem.
- Perform correction of contigs using [pilon](https://github.com/broadinstitute/pilon/wiki) and the Illumina dataset.

Usage
-----

Edit `config.mk`

Then issue issue `make all` to run the pipeline. Issue `make help` for a list of utility make targets.

Application dependencies
------------------------

- Canu [Canu](http://canu.readthedocs.io)
- samtools
- [BWA](http://bio-bwa.sourceforge.net)
- [racon](https://github.com/isovic/racon)
- [pilon](https://github.com/broadinstitute/pilon/wiki)

