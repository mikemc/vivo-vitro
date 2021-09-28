# Data from sequencing center S1 (Diversigen/CoreBiome)

## Source data from Diversigen

The sequence data and analysis report was downloaded by MM from an AWS-hosted website that served an html version of the QC analysis and 'Core Analysis' and provided links to download the following data files.

- fastq.tar: Fastq data as a single tarball with individual gzipped fastq files. Raw and quality-controlled reads are provided.
- fastq.tar.gz: Gzipped archive of fasta files, which were used as the input for Diversigen's BURST-based analyis pipeline.
- report.zip: The html report
- sample_overview.txt: Table of read counts input and mapped in the 'Core Analysis'
<!-- -->

The PicoGreen DNA concentration measurements and AbsQuant 16S qPCR measurements were provided to MM via email.
The samples were split over two plates; hence the two files for each measurement type.
The 'Sequencing performance' pdf gives a sense of the expected relationship between minimum sequencing depth and input DNA concentration.

## Within this folder

The raw fastq reads are within this folder in reads/raw.
The PicoGreen and AbsQuant files provide bulk and 16S DNA concentration measurements, respectively.
The file 'sample_overview.txt' was used for a fast QC check when the data was first received.

## Backups

The source data files are backed up in MM's sequence-data backup drive in /data/vivo-vitro/2021/s1.
