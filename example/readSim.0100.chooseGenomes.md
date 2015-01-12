


```r
library(MetamapsDB)
library(dplyr)
library(RJSONIO)
config = fromJSON("config.json")
```

Using NCBI's GENOME_REPORTS `prokaryotes.txt` file we determine which genera have complete genomes. 
Criteria for complete genomes include:
1. Status - Gapless Chromosome
2. A refseq sequence


```r
#Download MappingFile
dir.create("data/genomereports")
```

```
## Warning in dir.create("data/genomereports"): 'data/genomereports' already
## exists
```

```r
if(file.exists("data/genomereports/prokaryotes.txt"))
{
    print("Using existing prokaryotes files")
}else
{
    url = "ftp://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt"
    download.file(url, "data/genomereports/prokaryotes.txt")
}
```

```
## [1] "Using existing prokaryotes files"
```

```r
#Read in mappingfile
data        <- read.table("data/genomereports/prokaryotes.txt", sep ="\t", fill =T, comment.char="",h=T, quote="")
dataAfter   <- filter(data, Chromosomes.RefSeq != "-", Status == 'Gapless Chromosome')
```

Using MetamapsDB we find the taxID of the genera.

```r
#get the genus of these sequences
query0       <-
"START 
    basetaxa=node:ncbitaxid(taxid={taxid})
MATCH 
    basetaxa-[:childof*0..]->(higher)
WHERE 
    higher:genus OR higher:family OR higher:order OR higher:class OR higher:phylum
RETURN 
    head(labels(higher)) AS RANK, higher.taxid AS parentID ,basetaxa.taxid AS taxid"

#match taxID to appropriate rank eg. genus
taxonDF     <- 
make.data.frame(do.call(rbind,
    lapply(unique(dataAfter$TaxID), 
        function(taxID) dbquery(query0, list(taxid=taxID), cypherurl=config$cypherurl)
    )))

taxondf.genus   <-  taxonDF %>%
                    group_by(taxid) %>%
                    filter(RANK == 'genus')
```

Since there exists ambiguities in the string name of the taxa, ie. have 2 genera with the same genus epithet, 
we sort to resolve this by just taking the genus epithet which belonged to the Superkingdom of Bacteria or Archaea.

However, even then some ambiguities will not be resolved and we will remove these from our future analyses.


```r
abu     <-  setNames(read.table(config$abundance,sep="\t",h=T), c("taxon", "total"))
query1  <-
"MATCH
    (genus:genus)
WHERE
    genus.name = {taxaname}
WITH
    genus
MATCH
    path=genus-[:childof*]->(king:superkingdom)
WHERE
    king.name='Bacteria' XOR king.name = 'Archaea'
RETURN
    genus.name as taxon,
    genus.taxid as taxid"

genus2taxID <- do.call(rbind,lapply(unique(abu$taxon), 
function(taxID) 
    dbquery(query1, list(taxaname=taxID), cypherurl=config$cypherurl)
))

abundance = merge(abu, genus2taxID, by="taxon")
head(abundance)
```

```
##            taxon    total  taxid
## 1  Acaryochloris 667.4628 155977
## 2    Acetivibrio 128.3666  35829
## 3    Acetobacter 138.6274    434
## 4  Achromobacter 476.6540    222
## 5 Acidimicrobium 848.2083  53634
## 6   Acidiphilium 332.3117    522
```

We then filter the original list of completed genomes for organisms belonging to our given list.


```r
#Combined avaliable genomes with genera
combined=merge(taxondf.genus, abundance, by.x="parentID", by.y="taxid")
combined=merge(combined, dataAfter, by.x="taxid", by.y="TaxID")
```

We proceed next by randomly choosing a genome under each genus.

```r
#Choosen genome::luckydraw
set.seed(5)
chosen_genomes=do.call(rbind,
lapply(unique(combined$parentID), function(x) { 
       genus=subset(combined, parentID== x)
       genus[sample(1:nrow(genus),1),]
}))

write.table(
    chosen_genomes %>% 
    select(parentID,taxid,taxon,Chromosomes.RefSeq, total) %>%
    arrange(desc(total)), 
    file="out/readSim.0100.chosen_completeGenomes", 
    row.names=F, sep="\t", quote=F)

#Genera without complete genomes
write.table(
    data.frame(
        taxid=abundance$taxid[!abundance$taxid %in% unique(combined$parentID)]
        ), 
    file="out/readSim.0100.woCompleteGenomes_taxidList",
    row.names=F, sep="\t", quote=F)
write.table(abundance, file="out/readSim.0100.abundance_NameTaxid.txt",row.names=F,sep="\t",quote=F)
```
