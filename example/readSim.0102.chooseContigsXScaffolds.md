

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
```

```
## Warning in file(file, "rt"): cannot open file
## './out/readSim.0101.output.txt': No such file or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```r
str(data)
```

```
## function (..., list = character(), package = NULL, lib.loc = NULL, 
##     verbose = getOption("verbose"), envir = .GlobalEnv)
```

```r
#Exploratory analysis for the distribution of taxa within a given genus with relation to its total sequence length, number of sequences etc.
d2=data         %>% 
group_by(genus) %>%
summarise(numSeq=n(),
          numtaxa=length(unique(taxid)), 
          totLength=sum(as.numeric(combinedLength))) %>%
mutate(ratio = numSeq/numtaxa)
```

```
## Error in UseMethod("group_by_"): no applicable method for 'group_by_' applied to an object of class "function"
```

```r
numGenus = length(unique(data$genus))
```

```
## Error in data$genus: object of type 'closure' is not subsettable
```

```r
numTaxa  = length(unique(data$taxid))
```

```
## Error in data$taxid: object of type 'closure' is not subsettable
```










