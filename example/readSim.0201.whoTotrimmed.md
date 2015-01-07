



Trimming the database
====

## Sypnosis
For each of the taxa chosen (with and without completed genomes), we 
identity the corresponding species it belongs to. We proceed by removing all
species belonging under these species from the NR database



```r
library(MetamapsDB)
library(dplyr)
library(ggplot2)
library(RJSONIO)
theme_set(theme_bw())
config = fromJSON("config.json")
```

## Part1

Finding the species ID of chosen leaf taxa.


```r
chosen_genomes      = read.table("out/readSim.0100.chosen_completeGenomes",h=T, sep="\t", quote="")
chosen_nogenomes    = read.table("out/readSim.0102.chosen_scaffolds.txt",h=T, sep="\t", quote="")
#Request for species name of the chosen taxa
query = 
"START 
   basetaxa=node:ncbitaxid(taxid={taxaid}) 
MATCH 
   basetaxa-[:childof*0..]->(higher:species)-[:childof*0..]->(highest:genus)
RETURN 
   basetaxa.taxid as taxid, 
   higher.taxid as species, 
   highest.taxid as genus, 
   head(labels(higher)) as rankSpecies, 
   head(labels(highest)) as rankGenus"

##Part1::From the complete genome list 
completeGenomeSelectedSpecies =
do.call(rbind,
        lapply(as.character(unique(chosen_genomes$taxid)), function(x){
               dbquery(query,list(taxaid=x),cypherurl=config$cypherurl)
        })      )
```

```
## Error in fromJSON(content, handler, default.size, depth, allowComments, : invalid JSON input
```

```r
##Part2::From the incomplete genome list 
incompleteGenomeSelectedSpecies  = 
do.call(rbind,
        lapply(as.character(unique(chosen_nogenomes$taxid)), function(x){
               dbquery(query,list(taxaid=x),cypherurl=config$cypherurl)
        })      )

selectedSpecies = rbind(completeGenomeSelectedSpecies, incompleteGenomeSelectedSpecies)
```

```
## Error in rbind(completeGenomeSelectedSpecies, incompleteGenomeSelectedSpecies): object 'completeGenomeSelectedSpecies' not found
```








