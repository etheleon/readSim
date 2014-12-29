readSim
=======

## Introduction 

A Perl and R pipeline for simulating metagenomic reads from a highly complex community based on an abundance profile
with introduced sequencing errors based on empirically derived fastQ files of a 
full HiSeq 2500 Illumina run.

Abundance profile is estimated using RPM values of reads homology mapped to existing database (NCBI RefSeq) using MEGAN.

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
Caution the default for variable `cypherurl` in the `dbquery` `function in MetamapsDB is set to the local server’s address. 
User will have to species otherwise this 

## Description

### 100 Genome selection

**script**: [readSim.0100.chooseGenomes.r](readSim.0100.chooseGenomes.r), [readSim.0101](readSim.0101.chooseContigsXScaffolds.r), 

readSIM selects a random genome from an taxon belonging to the provided genus.

The selection of leaf taxa under given genera with  
1. complete genomes or 
2. with scaffolds/contigs WGS.

#### Complete genomes

In its first run, readSim takes abundance information and selects complete genomes (RefSeq and gapless)
This is done by first searching NCBI taxonomy for the leaf nodes of the given genera and randomly selecting for one.

#### Scaffolds and contigs WGS 

**script** [readSim.0102.chooseContigsXScaffolds.Rmd](readSim.0102.chooseContigsXScaffolds.md)
**R markdown** [readSim.0102.chooseContigsXScaffolds.md](readSim.0102.chooseContigsXScaffolds.md)

We then take the taxa which do not have a completed genome based on the above criteria and select contigs and scaffolds which have a WGS entry
Other filters include choosing only out of all possible leaf taxa under a genus based on the combined concatenated nucleotide length.


Chooses possible leaf taxa, for metagenomic shotgun sequencing recreation
Taxa with the longest nucleotide sequence length are chosen

### 0200

### 0400
Takes empircal fastQ files as template


