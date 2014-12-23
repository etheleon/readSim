#!/usr/bin/env Rscript

##################################################
#+------------------------------------------------
#Init
#+------------------------------------------------
##################################################
library(MetamapsDB)
library(dplyr)
library(functional)

args        <- commandArgs(T)
abuFile     <- args[1]
testing     <- args[2]
cypherurl   <- args[3]


if(args[1]){
    dbquery <<- Curry(dbquery, cypherurl=cypherurl)
}
##################################################
#+------------------------------------------------
#Part1: Find map genus (rank)'s taxonID to genomes
#+------------------------------------------------
##################################################

#Download MappingFile
mappingfile="data/genomereports/prokaryotes.txt"
if(file.exists(mappingfile))
{
    print("Using existing prokaryotes files")
}else
{
    #Download and update
    #Unwritten
}

#Read in mappingfile
data        <- read.table(mappingfile, sep ="\t", fill =T, comment.char="",h=T, quote="")
#Choose only genera with refseq sequences
dataAfter   <- filter(data, Chromosomes.RefSeq != "-", Status == 'Gapless Chromosome')

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
        function(taxID) dbquery(query0, list(taxid=taxID))
    )))

taxondf.genus   <-  taxonDF %>%
                    group_by(taxid) %>%
                    filter(RANK == 'genus')


##################################################
#+------------------------------------------------
#Part2: Read in abundance profile
#+------------------------------------------------
##################################################

abu     <-  setNames(read.table(abuFile,sep="\t",h=T), c("taxon", "total")

#Simulating using only genomes from bacteria or archeae, hence 
#genera in the abundance profile NOT of the above 2 superkingdoms will be ignored

#!!!WARNING!!!
#!!Match genus name to taxonID, occurences with ambiguities will occur.!!!

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

genus2taxID<- do.call(rbind,lapply(unique(abu$taxon), 
function(taxID) 
    dbquery(query1, list(taxaname=taxID))
))

abundance = merge(abu, genus2taxID, by="taxon")

#|taxon               |total    |taxid  |
#|--------------------|---------|-------|
#|      Acaryochloris |667.4628 | 155977|
#|      Acetivibrio   |128.3666 |  35829|
#|      Acetobacter   |138.6274 |    434|

#for writing to output as input to sim.0200
#write.table(abundance, file="out/sim.0101.out2.txt",quote=F, row.names=F,sep="\t")

##################################################
#+------------------------------------------------
#Part3
#+------------------------------------------------
##################################################

#Combined avaliable genomes with genera
combined=tbl_df(merge(taxondf.genus, abundance, by.x="parentID", by.y="taxid"))
combined=merge(combined, dataAfter, by.x="taxid", by.y="TaxID")

#Choosen genome::luckydraw
set.seed(5)
chosen_genomes=do.call(rbind,
lapply(unique(combined$parentID), function(x) { 
genus=subset(combined, parentID== x)
genus[sample(1:nrow(genus),1),]
}))

#OUTPUT
#Genera with complete genomes
write.table(chosen_genomes %>% 
            select(parentID,taxid,taxon,Chromosomes.RefSeq, total) %>%
            arrange(desc(total)), 
file="out/readSim.0100.chosen_completeGenomes", row.names=F, sep="\t", quote=F)

#Genera without complete genomes
write.table(
    data.frame(
        taxid=abundance$taxid[!abundance$taxid %in% unique(combined$parentID)]
        ), 
    sep="\t", row.names=F, quote=F,file="out/readSim.0100.woCompleteGenomes_taxidList")
