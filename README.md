# CoverMon
A script that does live read mapping of nanopore sequencing. While this can be run standalone, it's main use is to built in as part of other pipelines, see https://github.com/KMA-Aarhus/pappenheim

This is built for bacterial and viral genomes and is not suitable for larger genomes.

Also note that this requires unix exclusive gnome-terminal and browsersync software. While the scripts can be run on windows or macOS, the live updates will not work
## Requirements
- conda installation (mamba recommended
- browser-sync
- samplesheet linking barcodes to samples
- one or more references
## Installation from yaml file
```
mamba env create -f covermon.yaml
```
To run call the script with the following inputs in the following order:
(1) samplesheet
(2) path to run directory
(3) path to reference or reference directory
(4) threshold for minimum coverage
(5) maximum depth displayed in plot
Optional (6) a region file can be specified if all samples use the same reference genome.
```
mamba env create -f covermon.yaml
```
