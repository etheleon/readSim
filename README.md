readSim
=======

## Introduction 

A Perl and R pipeline for simulating metagenomic reads from a highly complex community based on an abundance profile
with introduced sequencing errors based on empirically derived fastQ files of a 
full HiSeq 2500 Illumina run.

Abundance profile is estimated using RPM values of reads homology mapped to existing database (NCBI RefSeq) using [MEGAN](http://ab.inf.uni-tuebingen.de/software/megan5/).

## Prerequisites

### Dependencies 

#### Data
 
* [NCBI Genome Reports](ftp://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS updated daily)
* NCBI taxonomy 
* Sequence database (NCBI RefSeq)

#### Software

* Neo4j database loaded with NCBI taxonomy (setup file is still being made)
* [MetamapsDB R package](https://github.com/etheleon/metamaps)
* Installed graphDB, installation pipeline is still in progress, please use existing database on water.bic.nus.edu.sg or from the public url [metamaps.scelse.nus.edusg:7474](http://metamaps.scelse.nus.edu.sg)

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

A [bootstrap](bootstrap.pl) perl script automates the whole process. 
An example input will be similar to the following

```
perl bootstrap.pl --abundance=input --cypherurl=graphDatabaseAddress:port/db/data/cypher --refseqDB_nucl=/home/user/db/refseq/ --taxDB=/home/user/db/NCBItaxonomy/ --threads=#integer
```

| Arguments     | Description                                                                                          |
| ----          | ----                                                                                                 |
| abundance     | input table with genusName and abundance (see above)                                                 |
| cypherurl     | The address of the graph database. Default for `cypherurl` set to head node on water.bic.nus.edu.sg. |
| refseqDB_nucl | NCBI’s refseq database (only archea and bacteria)                                                    |
| taxDB         | NCBI’s taxonomy database                                                                             |
| Threads       | The number of threads/CPU cores available on the machine                                             |

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

### 0200 Database trimming

**scripts:** [readSim.0201.whoTotrimmed.Rmd](readSim.0201.whoTotrimmed.Rmd)
**scripts:** [readSim.0202.trimDB.pl](readSim.0202.trimDB.pl)
To simulate real life metagenomic data (genomes without reference sequences), 
we trim NCBI’s nr database of sequences belonging under the same species as the chosen taxa from `0100 Genome selection`.


The first script hunts for all leaf taxa belonging to the same species as the chosen genomes, 
while the second reads in the NR database and removes the former to give a trimmed NR.

NOTE: We leave species in the genus thauera untrimmed from the NR database

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
readID|taxID|ncbiTaxonomyID|loc|location on concatenated sequence|output|lane/pair
```

Example:
```
>simuREAD_0|taxID|256616|loc|7302-7403|output|s_1/1
CCATCGCCATCCGCGCCGAGCGCGAAGGCATCACCGTCGAAGTCGCCATGTGGTGGAACGACAGCTACCACGAGAACGTGCTCTGCTTCACCAACAACATC
>simuREAD_0|taxID|256616|loc|7302-7403|output|s_1/2
GAGCGCGCCACGGAAGCCCGCGAGATGCGTGCCGCCATCGCGCTGCGGGATGTTGTTGGTGAAGCAGAGCACGTTCTCGTGGTAGCTGTCGTTCCACCACA
```

The indel rate was set as a percentage at `0.01`.
The insert size is set at 150 bp with SD of 5.

#### Breakdown of the script


##### Phred Score calculation

The phred and map functions stores 1.) Phred Quality Score and 2.) base-calling error probabilities 
as key-value pairs in hash `%score`.

```perl
sub phred
{
    my ($score) = @_;
    return 10**($score/(-10));
}

map {$score{$_}=phred($_)} 0..100
```


Function `gaussian_rand` generates a random number from a uniform distribution between 0 and 1
while `random_normal` uses the former to generate a number taken from a normal distribution.

The later is used to generate an insert size before truncation into read pairs.

```perl
sub gaussian_rand 
{
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g2 = $u1 * $w;
    $g1 = $u2 * $w;
    # return both if wanted, else just one
#    return wantarray ? ($g1, $g2) : $g1;
return $g1;
}

sub random_normal
{
    my ($mean, $sd) = @_;
    my $num = gaussian_rand();
    my $newnum = ($num * $sd) + $mean;
}
```

##### Substitution matrix

Substitution of template nucleotide upon failure of accurate base calling (based on phred score from template fastQ) is random (uniform distribution).

```perl
sub rand_nt
{
    my $rand = rand;
    return 'a' if($rand<$globalnt{'a'});
    return 't' if($rand<$globalnt{'a'}+$globalnt{'t'});
    return 'g' if($rand<$globalnt{'a'}+$globalnt{'t'}+$globalnt{'g'});
    return 'c' if($rand<$globalnt{'a'}+$globalnt{'t'}+$globalnt{'g'}+$globalnt{'c'});
    return 'n';
}
```

The mutation is based on ....
```perl
sub mutate
{
    my ($readnt,$qual,$readLength) = @_;
    my @readnt = split(//, $readnt);
    my @qual = split(/\t/, $qual);		#stores the probability
    my $loc = 1;				#this is loc of buffered seqeunce in case of deletion event

    my @outputsequence;	#store sequence

    for(my $i=0; $i < $readLength; $i++) {
        my $qualityscore = $qual[$i];
        if(rand()<$qualityscore){
            #MUTATION
            #INDEL EVENT #####################################
            if(rand()<$indelrate)	{
                ##################################################
                if(int(rand(2))){ 	#INSERTION
                    my $extra = rand_nt();
                    my $extra.= $readnt[$loc];
                    push @outputsequence, $extra;
                    $loc++
                }else{			#DELETION
                    $loc++;
                    push @outputsequence, $readnt[$loc];
                    $loc++
                }
                ##################################################

                #SUBSTITUTION EVENT ##############################
            }else	{
                ##################################################
                my $substitutedNT = rand_nt();
                push @outputsequence, $substitutedNT;
                $loc++
            }
            ##################################################
        }else{
            #NO MUTATION
            push @outputsequence, $readnt[$loc];
            $loc++;
        }
    }
    my $finalseq = join '', @outputsequence;
    return $finalseq;
}
```
