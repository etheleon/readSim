readSim
=======

A Perl and R pipeline for simulating metagenomic reads from a highly complex community 
with introduced sequencing errors based on empirically derived fastQ files of a 
full HiSeq 2500 Illumina run.


## Prerequisites

### Dependencies 

#### Data
 
* [NCBI Genome Reports](ftp://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS updated daily)
* NCBI taxonomy 
* Sequence database (NCBI RefSeq)

#### Software

* Neo4j database loaded with NCBI taxonomy (setup file is still being made)
* [MetamapsDB R package](https://github.com/etheleon/metamaps)

## User input

* Requires input table consisting of Genera to be simulated and abundance.

**Example**

|Genus                      |Abundance          |
|---------------------------|-------------------|
|Caldilinea                 |25414.868178775    |
|Nitrospira                 |20661.7492551171   |
|Sorangium                  |18241.7170368186   |
|Mycobacterium              |14892.9966161738   |
|Candidatus Accumulibacter  |11906.5856237162   |
Caution: Currently only supports simulating using full name.
This is flawed because there are genera with the same genus epithet. 
eg. *missing*

## Installation
After cloning the repository, 


## Usage
Caution the default for variable `cypherurl` in the `dbquery function` `in MetamapsDB is set to the local server’s address. 
User will have to species otherwise this 

## Description

### Genome selection

**script**: [readSim.0100.chooseGenomes.r](readSim.0100.chooseGenomes.r)

readSIM selects a random genome from an taxon belonging to the provided genus.
The selection can be set to choose only from complete genomes only or 
with scaffolds/contigs taken from WGS data.

#### Complete genomes

In its first run, readSim takes abundance information and selects complete genomes (RefSeq and gapless)
This is done by first searching NCBI taxonomy for the leaf nodes of the given genera and randomly selecting for one.

#### Scaffolds and contigs WGS (Optional)


