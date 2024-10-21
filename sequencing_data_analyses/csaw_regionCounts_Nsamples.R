#This script uses the function reg.counts 
#from the csaw package to quantify a target bed file
#across many samples 

#Arguments:
#1) comma separated list of bam files
#2) bed file
#3) outfile

cat(date())

library(csaw)

args=commandArgs(trailingOnly = T)

comma_bams=args[1]
bam.paths=strsplit(comma_bams,",")[[1]]

bed.file=args[2]

outfile=args[3]

bed=read.table(bed.file)
my.regions <- GRanges(bed$V1,
                      IRanges(bed$V2, 
                              bed$V3))
param=readParam(pe = "both") #this has to be changed if it is not paired-end data
reg.counts <- regionCounts(bam.files = bam.paths,regions = my.regions,param = param)

#saveRDS(reg.counts,paste0(outfile,".RDS"))
reg.counts=assay(reg.counts)
rownames(reg.counts)=paste(seqnames(my.regions),ranges(my.regions),sep = ":")
colnames(reg.counts)=basename(bam.paths)

write.table(reg.counts,outfile,quote = F,sep="\t")
cat(date())


