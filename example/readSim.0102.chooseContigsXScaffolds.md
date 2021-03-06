

Selecting leaf taxa for genera without complete genomes.
====



```r
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
```

## Introduction
For genera without leaf taxa (with completed genome), 
we select a viable leaf taxon based on (1) availability of of WGS records and (2) its total sequence length.

## Exploratory analysis

For each genus belongs numerous leaf taxa. 
We investigate:
1. total number of sequences (within the genus), 
2. the number of taxa within the genus
3. the combined sequence lengths


```r
data=read.table("./out/readSim.0101.output.txt",sep="\t", h=T)
str(data)
```

```
## 'data.frame':	25344 obs. of  5 variables:
##  $ taxid         : int  1304270 1336208 1336208 1336208 1336208 1336208 1336208 1336208 1336208 1336208 ...
##  $ genus         : int  180541 125216 125216 125216 125216 125216 125216 125216 125216 125216 ...
##  $ gi            : int  631253143 653028352 653028789 653029045 653029202 653029333 653029441 653029549 653029663 653029762 ...
##  $ refseq        : Factor w/ 25344 levels "NC_002806.1",..: 406 22898 22899 22900 22901 22902 22903 22904 22905 22906 ...
##  $ combinedLength: int  1462 4612006 4612006 4612006 4612006 4612006 4612006 4612006 4612006 4612006 ...
```

```r
#Exploratory analysis for the distribution of taxa within a given genus with relation to its total sequence length, number of sequences etc.
d2=data         %>% 
group_by(genus) %>%
summarise(numSeq=n(),
          numtaxa=length(unique(taxid)), 
          totLength=sum(as.numeric(combinedLength))) %>%
mutate(ratio = numSeq/numtaxa)
numGenus = length(unique(data$genus))
numTaxa  = length(unique(data$taxid))
```

A total of 80 genera were without completed completed genomes. 
A total of 640 leaf taxa under the former.


```r
ggplot(d2, aes(y=numSeq, x=numtaxa))+
geom_point(aes(size=sqrt(totLength)))+
geom_text(data=subset(d2, numSeq > 20000| numtaxa >=50 ), aes(label=genus), color='red')+
xlab("# taxa")+ylab("# sequences")
```

![plot of chunk unnamed-chunk-4](figures/readSim.0102-unnamed-chunk-4-1.png) <center><p class="caption"><b>Figure:</b> Number of fastA sequences within a genus against the number of taxa within a the same genus. Annotated genus are based on the following: > 20000 sequences or >50 taxa</p></center>

## Data Processing

We choose the taxa based on the availability of WGS data and total sequence length.

```r
#Will have to write the script which produces this wgstaxa, taken from NCBI

#read in wgs annotation list
wgs = setNames(read.table("~/db/refseq/wgstaxa",h=F), c("taxa"))
wgsdata = subset(data, taxid %in% wgs$taxa)

#choose max length taxa
chosen = wgsdata %>% 
group_by(genus) %>% 
filter(combinedLength == max(combinedLength))
str(chosen)
```

```
## Classes 'grouped_df', 'tbl_df', 'tbl' and 'data.frame':	7861 obs. of  5 variables:
##  $ taxid         : int  1120941 1120941 1120941 1120941 1120941 1120941 1120941 1120941 1120941 1120941 ...
##  $ genus         : int  1654 1654 1654 1654 1654 1654 1654 1654 1654 1654 ...
##  $ gi            : int  651276016 651276137 651276214 651276316 651276388 651276541 651276587 651276628 651276679 651276893 ...
##  $ refseq        : Factor w/ 25344 levels "NC_002806.1",..: 14529 14530 14531 14532 14533 14534 14535 14536 14537 14538 ...
##  $ combinedLength: int  3531393 3531393 3531393 3531393 3531393 3531393 3531393 3531393 3531393 3531393 ...
##  - attr(*, "vars")=List of 1
##   ..$ : symbol genus
```

### WGS availability 

404 taxa were removed because they have no WGS projects listed in NCBI, this resulted in 516 sequences trimmed out of a total of 25344.


```r
ggplot() + 
geom_boxplot(data=wgsdata, aes(x=reorder(as.factor(genus), combinedLength, median), y=combinedLength))+
geom_jitter( data=wgsdata, aes(x=reorder(as.factor(genus), combinedLength, median), y=combinedLength))+
geom_point(data=subset(wgsdata, taxid %in% chosen$taxid), aes(x=as.factor(genus), y=combinedLength),color='red')+
theme(axis.text.x = element_text(angle=90, hjust=1))+
xlab("Genera")+ylab("Total Length of sequences in leaf taxa")
```

![plot of chunk unnamed-chunk-6](figures/readSim.0102-unnamed-chunk-6-1.png) <center><p class="caption"><b>Figure:</b>We plot the concatenated lengths of the sequences and across genera. Each data point represents one leaf taxa. Highlighted in red are the selected taxa.</p></center>


```r
write.table(chosen %>% select(-gi), file="out/readSim.0102.chosen_scaffolds.txt", sep="\t", quote=F, row.names=F, 
col.names=c("taxid","parentID","Chromosome.RefSeq", "totalLength"))
```


```r
sessionInfo()
```

```
## R version 3.1.1 (2014-07-10)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggplot2_1.0.0.99 dplyr_0.3.0.2    knitr_1.8       
## 
## loaded via a namespace (and not attached):
##  [1] assertthat_0.1   colorspace_1.2-4 DBI_0.3.1        digest_0.6.8    
##  [5] evaluate_0.5.5   formatR_1.0      grid_3.1.1       gtable_0.1.2    
##  [9] labeling_0.3     lazyeval_0.1.10  magrittr_1.5     MASS_7.3-35     
## [13] munsell_0.4.2    parallel_3.1.1   plyr_1.8.1       proto_0.3-10    
## [17] Rcpp_0.11.3      reshape2_1.4.1   scales_0.2.4     stringr_0.6.2   
## [21] tools_3.1.1
```
