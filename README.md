readSim
=======

Perl and R pipeline for simulating metagenomic reads from a highly complex community 
with introduced sequencing errors based on empirically derived fastQ files of a 
full HiSeq 2500 Illumina run.


## Prerequisites

### Dependencies 

#### Data
 
* [NCBI Genome Reports](ftp://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS updated daily)
* NCBI taxonomy 
* NR 
* Abundance profile

#### Software

* Neo4j database loaded with NCBI taxonomy (setup file is still being made)


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

## Description

### Genome selection

readSIM selects a random genome from an taxon belonging to the provided genus.
The selection can be set to choose only from complete genomes only or 
with scaffolds/contigs taken from WGS data.

#### Complete genomes

In its first run, readSim takes abundance information and selects complete genomes (RefSeq and gapless)
This is done buy first

#### Scaffolds and contigs WGS (Optional)
