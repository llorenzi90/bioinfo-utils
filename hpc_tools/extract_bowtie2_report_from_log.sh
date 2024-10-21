awk '/reads; of these:/,/overall alignment rate/' $1

#examples of usage: ls |grep .out|while read INPUT;do SUBSTRING=$(echo $INPUT| cut -d'_' -f 1- --output-delimiter="." |cut -d"." -f1,2 );echo $SUBSTRING; /home/llorenzi/scripts/extract_bowtie2_report_from_log.sh $INPUT > $SUBSTRING.bowtie2_report.txt ;done
#grep "Bowtie2 done" *out|cut -d ":" -f1|while read INPUT;do SUBSTRING=$(echo $INPUT| cut -d'_' -f 1- --output-delimiter="." |cut -d"." -f1,2 --output-delimiter="_" );echo $SUBSTRING;/home/llorenzi/scripts/extract_bowtie2_report_from_log.sh $INPUT > $SUBSTRING.bowtie2_report.txt ;done
