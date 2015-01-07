readSim
=======

## Introduction 

A Perl and R pipeline for simulating metagenomic reads from a highly complex community. 
The simulation is based on an \*abundance profile, includes introduced sequencing errors based 
on empirically derived fastQ files.

\*Abundance profile is estimated using RPM values of reads homology mapped to existing database (NCBI RefSeq) using [MEGAN](http://ab.inf.uni-tuebingen.de/software/megan5/).

## Prerequisites

### Dependencies 

#### Data

* [NCBI Genome Reports](ftp://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS updated daily)
* NCBI taxonomy 
* Sequence database 
  * NCBI RefSeq
  * NCBI NR (protein)

#### Software

* Perl 5.18
  * autodie
  * Modern::Perl
  * Parallel::ForkManager
  * Statistics::R
  * File::Basename
  * Bio::SeqIO
  * REST::Neo4p
  * REST::Neo4p::Query
  * Cwd

* R 3.1.1
  * ggplot2
  * dplyr
  * [MetamapsDB R package](https://github.com/etheleon/metamaps) for interfacing with the graphDB
  * getopt

* [Neo4j](http://neo4j.com/download/)
* Installed graphDB, installation pipeline unpublished, msg me if interested otherwise please use existing database on water.bic.nus.edu.sg or from the public url [metamaps.scelse.nus.edusg:7474](http://metamaps.scelse.nus.edu.sg).

## User input
1. [Input table](./example/data/abundanceProfile.txt) of 1. **Genera** and **abundance**
2. [Input table](./example/data/fastqLocations) mapping the fastq pairs 
3. [cypherurl](./cypherurl) which lists neo4j address

|Genus                      |Abundance          |
|---------------------------|-------------------|
|Caldilinea                 |25414.868178775    |
|Nitrospira                 |20661.7492551171   |
|Sorangium                  |18241.7170368186   |
|Mycobacterium              |14892.9966161738   |
|Candidatus Accumulibacter  |11906.5856237162   |

Caution: Currently only supports simulating using full name.
BUG:: This is flawed because there are genera with the same genus epithet. Pipeline throws away some instances where the above cannot be resolved

## Usage

A [bootstrap](bootstrap.pl) perl script automates the whole process. 
An example input will be similar to the following:

```
perl bootstrap.pl config.json 20
```
an example config.json can be found [here](config.json)

## Description

### 100 Genome selection

**script**: [readSim.0100.chooseGenomes.r](readSim.0100.chooseGenomes.r), 

readSIM selects a random genome from a leaf taxon for each genus in list.

The selection of leaf taxa under given genera with  
1. complete genomes or 
2. scaffolds/contigs from WGS projects.

#### Complete genomes 

In its first run, readSim takes abundance information and selects complete genomes (RefSeq and gapless)
This is done by first searching NCBI taxonomy for the leaf nodes of the given genera and randomly selecting for one.

#### Scaffolds and contigs WGS 

**scripts:** [readSim.0101.contigsXscaffolds.pl](readSim.0101.contigsXScaffolds.pl), [readSim.0102.chooseContigsXScaffolds.Rmd](example/readSim.0102.chooseContigsXScaffolds.md)

We then take the taxa which do not have a completed genome based on the above criteria and select contigs and scaffolds which have a WGS entry
Other filters include choosing only out of all possible leaf taxa under a genus based on the combined concatenated nucleotide length.

Chooses possible leaf taxa, for metagenomic shotgun sequencing recreation
Taxa with the longest nucleotide sequence length are chosen

##### Outputs

| Output                              | Description                                           |
| -----                               | -----                                                 |
| readSim.0100.chosen_completeGenomes | list of chosen genomes and their refseqIDs            |
| readSim.0101.output.txt             | all sequences under genera with completed genomes     |
| readSim.0101.output.txt             | all leaf taxa under genera without a completed genome |
| readSim.0102.chosen_scaffolds.txt   | list of chosen scaffolds and their refseqIDs          |

### 0200 Reference Protein Database trimming

**scripts:** [readSim.0201.whoTotrimmed.Rmd](readSim.0201.whoTotrimmed.Rmd)
**scripts:** [readSim.0202.trimDB.pl](readSim.0202.trimDB.pl)

To simulate a closer to reality scenerio where metagenomic reference sequences are not found in the any of the databases, 
we prepare a protein reference database (based on nr) which we will carry out a distant homology search using a blastx like aligner eg. diamond or rapsearch (preference be diamond because of the the computational speed up).

We trim NCBI’s nr protein database of sequences belonging under the same species as the chosen taxa from `0100 Genome selection`.

The first script hunts for all leaf taxa belonging to the same species as the chosen genomes, 
while the second reads in the NR database and removes the former to give a trimmed NR.

NOTE: We leave species in the genus Thauera untrimmed from the NR database.

##### Outputs

| Output                  | Description                                               |
| -----                   | -----                                                     |
| readSim.0201.output.txt | the species of the chosen leaf taka and the exact taxonID |

### 0300 template Sequence extraction

**scripts** [readSim.0300.genome_extraction.pl](readSim.0300.genome_extraction.pl)

Extracts and concatenate sequences from the chosen leaf taxon belonging to genera in original list 
and provides the mapping locations of the individual sequences on the global sequence.

**scripts** [readSim.0301](readSim.0301.getLength.pl), [readSim.0302](readSim.0301.sequenceLengths.Rmd)

* counts the lengths of the sequences and plots and partitions them based on genome completion status.

##### Outputs

| Output                          | Description                                                                                 |
| -----                           | -----                                                                                       |
| readSim.0300 dir                | contains fasta files of the leaf taxa belonging to genera of interest; named by the genusID |
| readSim.0300.out.mapping        | the mapping of the locations                                                                |
| readSim.0300.out.notRemoved     | sequences which were not found in the refseq genomic database                               |
| readSim.0300.out.fna            | collection of sequences used                                                                |
| readSim.0301.lengthtable.txt    | the sequence lengths of any given leaf taxon under a genus                                  |
| readSim.0302.sequenceLengths.md | summary report of the number of sequences                                                   |

### 0400 Generating simulated metagenomic reads

[In progress]
Takes empircal fastQ files from the Illumina platform as template, and generates the simualted reads.
Input fastQ files partitioned by lane.

#### Header of output files 
```
readID|ncbiTaxonomyID|location on concatenated sequence|outputFile/pairID
```

Example:
```
>READ_0|256616|7302-7403|s_1/1
CCATCGCCATCCGCGCCGAGCGCGAAGGCATCACCGTCGAAGTCGCCATGTGGTGGAACGACAGCTACCACGAGAACGTGCTCTGCTTCACCAACAACATC
>READ_0|256616|7302-7403|s_1/2
GAGCGCGCCACGGAAGCCCGCGAGATGCGTGCCGCCATCGCGCTGCGGGATGTTGTTGGTGAAGCAGAGCACGTTCTCGTGGTAGCTGTCGTTCCACCACA
```

The indel rate was set as a percentage at `0.01`.
The insert size is set at 150 bp with SD of 5.
