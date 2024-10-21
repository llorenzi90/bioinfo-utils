#check md5sum between two files
#usage: Rscript compare_md5sums.R origin_md5sum_file destination_md5sum_file 
#these two files are outputs from md5sum command
file1=commandArgs(trailingOnly=TRUE)[1]
file2=commandArgs(trailingOnly=TRUE)[2]
origin_file=read.table(file1)
dest_file=read.table(file2)

cat("comparison md5sums in destination vs origin:\n")
print(table(dest_file$V1%in%origin_file$V1))

if(any(!dest_file$V1%in%origin_file$V1)){
  cat("\nthese files do not match:\n")
  print(dest_file[!dest_file$V1%in%origin_file$V1,])
}


cat("\ncomparison md5sums in origin vs destination:\n")
print(table(origin_file$V1%in%dest_file$V1))

if(any(!origin_file$V1%in%dest_file$V1)){
  cat("\nthese files do not match:\n")
  print(origin_file[!origin_file$V1%in%dest_file$V1,])
}

