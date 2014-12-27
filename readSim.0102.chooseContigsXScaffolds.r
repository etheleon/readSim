#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

data=read.table("out/sim.0102.out",sep="\t", h=T)
bplength = setNames(read.table("out/sim.0102.out2",h=F, sep="\t"), c("parentGenus", "bplength"))


d2=data %.% 
group_by(parentGenus) %.%
summarise(numSeq=n(),numtaxa=length(unique(taxid)))

d2$ratio=d2$numSeq/d2$numtaxa
d3=merge(d2, bplength, by="parentGenus")

d3$bplength=as.integer(as.character(d3$bplength))
p0 = ggplot(d3, aes(y=numSeq, x=numtaxa))+
geom_point(aes(size=sqrt(bplength)))+
geom_text(data=subset(d3, numSeq > 20000| numtaxa >=50 ), aes(label=parentGenus), color='red')+
ggtitle("Number of Sequences against Number of taxa per Genera")+
xlab("# taxa")+ylab("# sequences")
ggsave(plot=p0, filename="out/sim.0103.summary1.pdf")

#pdf("out/sim.0103.summary.pdf",w=10,h=10)
#lapply(unique(data$parentGenus), function(x) { 
#  df=subset(data, parentGenus == x) 
#  df2=df %.% 
#  group_by(taxid) %.%
#  summarise(numSeq=n())
#  qplot(reorder(as.factor(taxid),numSeq), log10(numSeq), data=df2, geom="bar", stat="identity")+
#  theme(axis.text.x=element_text(angle=45, hjust=1,vjust=1))+ggtitle(x)
#})
#dev.off()


###########
#luckydraw#
###########

#length(unique(data$parentGenus))
#80
#length(unique(subset(data, taxid %in% wgs$taxa)$parentGenus))
#75

wgs = setNames(read.table("~/db/refseq/wgstaxa",h=F), c("taxa"))
data = subset(data, taxid %in% wgs$taxa)

lengthdata = read.table("out/sim.0900.out.txt",sep="\t",h=T)

#set.seed(5)
#chosen = do.call(rbind,lapply(unique(data$parentGenus), function(x) { 
#    chosen=sample(size=1, as.character(unique(subset(data, parentGenus == x)$taxid)))
#    subset(data, taxid == chosen)
#}))

#non-random choosing
data2=unique(merge(data, lengthdata, by.x="taxid",by.y="Taxon",all.x=T)[,c("taxid","Genus","Length")])
chosen = do.call(rbind,lapply(unique(data$parentGenus), function(x) { 
    	    df = subset(data2, Genus == x)
	    subset(data, taxid == df[which(df$Length == max(df$Length)),]$taxid)
}))
straindf=do.call(rbind,lapply(unique(chosen$taxid), function(x) { 
    df = subset(chosen , taxid == x)
    data.frame(taxid=x, strains=
    length(unique(substr(df$refid, 1, 10
))))
    }))

p10=ggplot() + 
geom_boxplot(data=data2, aes(x=reorder(as.factor(Genus), Length, median), y=Length))+
geom_jitter( data=data2, aes(x=reorder(as.factor(Genus), Length, median), y=Length))+
geom_point(data=subset(data2, taxid %in% chosen$taxid), aes(x=as.factor(Genus), y=Length),color='red')+
theme(axis.text.x = element_text(angle=90))
ggsave(plot=p10, file="out/sim.0103.chosenAgainstRest.pdf",w=20)

write.table(chosen, file="out/sim.0103.chosen.txt", sep="\t", quote=F, row.names=F)
