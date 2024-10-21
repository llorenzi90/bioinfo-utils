#Usage: du -sh * | Rscript /path/to/scripts/check_directory_sizes.R

cArgs=commandArgs(trailingOnly = T)
stf=file("stdin")
if(length(cArgs)==0){
  print("No command arguments")
  sizes=read.table(file = stf)

} else{
  print("command arguments passed")
   sizes=read.table(cArgs[1])

 }

#sizes=read.table(file = 'stdin')
options(scipen=999)
sizes$V1 <- gsub(",",".",sizes$V1)
values=as.numeric(gsub("([0-9.]+)([a-zA-z]+)","\\1",sizes$V1))
units=sub("([0-9.]+)([a-zA-z]+)","\\2",sizes$V1 ,)
units[grepl("[0-9.]+", units)]="B"
conv=c(T=0.001,G=1,M=1000,K=1000000,B=1000000000)
GBvalues=values/conv[match(units,names(conv))]
names(GBvalues) <- sizes$V2
print(sort(GBvalues))
print(sum(GBvalues))
