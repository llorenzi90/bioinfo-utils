#!/bin/bash
die() {
	printf '%s\n' "$1" >&2
    exit 1
    }

###### ***Parameters that might change with time (because of changes in hpc modules, etc)*** ######
#define variables for loading modules, change if necessary
#load_fastqc="module load FastQC/0.11.9"
load_fastqc="module load FastQC/0.12.1-Java-11"
load_trimgalore="module load TrimGalore/0.6.6"
#load_bowtie2="module load bowtie2/2.4.4 gcc-11.2.0-gcc-4.8.5-pjl2s6u"
load_bowtie2="module load Bowtie2/2.4.4-GCC-11.2.0"
load_samtools="module load SAMtools/1.13-foss-2021b"
load_picard="module load picard/2.26.3-Java-11"
load_deeptools="module load deepTools/3.3.1-foss-2021b-Python-3.8.5"

#indexes and ref files.
bowtie2Index="/mnt/beegfs/public/references/index/bowtie2/GRCm39_UCSC/GRCm39_UCSC"
#bowtie2Index="/mnt/beegfs/public/references/index/bowtie2/mm10_UCSC/mm10"
hg38Index="/mnt/beegfs/public/references/index/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as"
mm10Index="/mnt/beegfs/public/references/index/bowtie2/mm10_UCSC/mm10"
min_frag_len=1000
primary_chrs_bed="/ijc/USERS/llorenzi/references/primary_chrs_nochrM.bed"

show_help() {
cat << EOF
		Usage: ${0##*/} [options] [working_dir]...
		Run ChIP-seq basic data processing in IJC hpc;
		 from QC to alignment to generation of bigWig files for visualization 
		 
		Available options are:
			
			-h			display this help and exit
			-rs runstep		run the analysis starting from this particular step. All steps are run by default (runstep=all) options: all,qc,trim,align,sort_index,markdups,filtstats,bigWig,rmdups_bigWig.
			-only			if this flag is set (no additional arguments) then only the step indicated by -rs will be run (default FALSE)
			-ppn integer		number of proccesors per node (default 10)
			-mm39			use bowtie2 index mouse genome mm39 (default): ${bowtie2Index}
			-hg38			use bowtie2 index human genome hg38: ${hg38Index}
			-mm10			use bowtie2 index mouse genome mm10: ${mm10Index}
			-idx_path		use custom path to bowtie2 index
			-atacseq		set this flag to indicate data is ATAC-seq (default ChIP-seq)	
		   
EOF
}
	
runstep=all
PPN=$SLURM_CPUS_PER_TASK
only=false
atacseq=false
#parse optional arguments	   

while :; do
	case $1 in
		-h|-\?|--help)
			show_help    # Display a usage synopsis.
			exit
			;;
		-rs|--runstep)       # Takes an option argument; ensure it has been specified.
			re='^-'
			if [ "$2" ] &&  ! [[ $2 =~ $re ]]; then
				case $2 in 
				  "all"|"qc"|"trim"|"align"|"sort_index"|"markdups"|"filtstats"|"bigWig"|"rmdups_bigWig")
					  runstep=$2
					  echo "running step $runstep"
				  ;;
				  *)
				    die 'ERROR: non-valid runnning step
					"-rs" requires one of these non-empty option arguments : "all", "qc", "trim", "align", "sort_index", "markdups", "filtstats", "bigWig" or "rmdups_bigWig".'
				  ;;
				esac
	            		shift 2
	        	else
	            		die 'ERROR: "-rs" requires one of these non-empty option arguments : "all", "qc", "trim", "align", "sort_index", "markdups", "filtstats" , "bigWig" or "rmdups_bigWig". '
			fi
			;;
		-ppn)
			if [ "$2" ]; then
				PPN=$2
				re='^[0-9]+$'
				if ! [[ $PPN =~ $re ]] ; then
				   die 'ERROR:  "-ppn" must be followed by an integer' 
				fi
	            	shift 2
	        	else
	           		 die 'ERROR: "-ppn" requires an integer as option argument'
			fi
			;;
		-only)
			only=true
			shift
			echo "Only $runstep will be run"
			;;
		
		-mm39)
			shift			
			;;
		-hg38)	
			bowtie2Index=$hg38Index
			shift			
			;;
		-mm10)	
			bowtie2Index=$mm10Index
			shift			
			;;
		-idx_path)       # Takes a path to custom index; ensure directory exists.
			dn=$(dirname $2)
			if [[ ! -d $dn ]]; then
				die "ERROR: index directory ${dn} does not exist"
			else
				bowtie2Index=$2
			fi
			shift 2
			;;
		-atacseq)
			atacseq=true
			shift
			;;
		
		-?*)
			printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
			show_help
			shift
	    		;;
				  
	    	*)  # Default case: No more options, so break out of the loop.
	    		break
	esac
done	    	





echo "PPN: $PPN"
echo "only: $only"
echo "wdir: $1"



if [ $runstep == all ] && $only
then
	only=false
	echo 'Note that when -rs all then -only must by FALSE, setting -only back to false'
fi

if [[ $atacseq = true ]]
then 
	min_frag_len=2000
	echo "processing ATAC-seq data"
else
	min_frag_len=1000
	echo "processing ChIP-seq data"
fi

#Now the actual analysis starts 
#NOTE: this is a preliminary version, in the future I would like to make a function out of each step for the sake of organisation

sdir=$1
sample=$(basename $sdir) 


cd $sdir

#1) QC
step=qc
if [ $runstep == all ] || [ $runstep == $step ]
then
	read1=$(ls *_1.fq.gz)
	read2=$(ls *_2.fq.gz)
	
	if ! [ -f "$read1" ]; then
		die "input file $read1 not found"
	fi
	read2=$(ls *_2.fq.gz)
	if ! [ -f "$read2" ]; then
		die "input file $read2 not found"
	fi
	echo "read1: $read1"
	echo "read2: $read2"
	
	echo `date`
	echo -e "################################\n\tFASTQC started\n################################\n"
	$load_fastqc
	mkdir fastqc
	fastqc $read1 -o fastqc
	fastqc $read2 -o fastqc
	echo -e "################################\n\tFASTQC done\n################################\n"
	echo `date  +'%r'`
	if $only
	then 
		exit 0
	else 
		runstep=all
	fi
fi


#2)trim_galore
step=trim
if [ $runstep == all ] || [ $runstep == $step ]
then
	read1=$(ls *_1.fq.gz)
	if ! [ -f "$read1" ]; then
		die "input file $read1 not found"
	fi
	read2=$(ls *_2.fq.gz)
	if ! [ -f "$read2" ]; then
		die "input file $read2 not found"
	fi
	echo "read1: $read1"
	echo "read2: $read2"
	
	echo `date  +'%r'`
	echo -e "################################\n\tTrimGalore started\n################################\n"

	$load_trimgalore

	trim_galore --paired --length 35 --fastqc $read1 $read2
	echo -e "################################\n\tTrimGalore done\n################################\n"
	echo `date  +'%r'`
	if $only
	then 
		exit 0
	else 
		runstep=all
	fi
fi
	

#3)alignment
step=align
if [ $runstep == all ] || [ $runstep == $step ]
then
	read1_val=$(ls *val_1.fq.gz)
	if ! [ -f "$read1_val" ]; then
		die "input file $read1_val not found"
	fi
	read2_val=$(ls *val_2.fq.gz)
	if ! [ -f "$read2_val" ]; then
		die "input file $read2_val not found"
	fi
	echo "read1: $read1_val"
	echo "read2: $read2_val"
	echo "bowtie2 index ${bowtie2Index} will be used"
	echo `date  +'%r'`
	echo -e "################################\n\tBowtie2 started\n################################\n"
	
	$load_bowtie2
	
	bowtie2 --verbose --very-sensitive -X $min_frag_len -x $bowtie2Index -1 $read1_val -2 $read2_val -S $sample.bt2.sam
	echo -e "################################\n\tBowtie2 done\n################################\n"
	echo `date  +'%r'`

	
	if $only
	then 
		exit 0
	else 
		runstep=all
	fi
fi
	
#4) convert to bam while sorting and index
step=sort_index
if [ $runstep == all ] || [ $runstep == $step ]
then
	input=$sample.bt2.sam
	if ! [ -f "$input" ]; then
		die "input file $input not found"
	fi
	
	echo `date  +'%r'`
	echo -e "################################\n\tsamtools sort started\n################################\n"
	
	$load_samtools
	
	samtools sort -@ $PPN -O bam -o $sample.sorted.bam $input
	samtools index -@ $PPN $sample.sorted.bam
	
	echo -e "################################\n\tsamtools sort done\n################################\n"
	echo `date  +'%r'`
	
	if $only
	then 
		exit 0
	else 
		runstep=all
	fi
fi



#5) mark duplicates
step=markdups
if [ $runstep == all ] || [ $runstep == $step ]
then 
	input=$sample.sorted.bam
	if ! [ -f "$input" ]; then
		die "input file $input not found"
	fi
	
	echo `date  +'%r'`
	echo -e "################################\n\tpicard Mark duplicates started\n################################\n"
	$load_picard	

	java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$input O=$sample.sorted.markedDups.bam M=$sample.dupmatrix CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false TAGGING_POLICY=All
	
	echo -e "################################\n\tpicard Mark duplicates done\n################################\n"
	echo `date  +'%r'`
	
	if $only
	then 
		exit 0
	else 
		runstep=all
	fi
fi


#6) MAPQ distribution, raw and proper pairs, filtering and stats
step=filtstats
if [ $runstep == all ] || [ $runstep == $step ]
then 
	input=$sample.sorted.markedDups.bam
	if ! [ -f "$input" ]; then
		die "input file $input not found"
	fi
	
	$load_samtools
	
	if [[ $atacseq = true ]]
	then
		echo `date  +'%r'`
		echo -e "################################\n\tsamtools view and stats\n################################\n"
		#A) FILTERS: primary chromosomes only (exclude chrM reads), proper pairs and MAPQ > 2
		#A.1) exclude chrM reads
		samtools view -@ $PPN -b -h -L $primary_chrs_bed $sample.sorted.markedDups.bam > $sample.sorted.markedDups.primary_chr.bam
		#A.2) retain only proper pairs and reads with minimum MAPQ=2
		samtools view -@ $PPN -b -h -f2 -q2 $sample.sorted.markedDups.primary_chr.bam > $sample.sorted.markedDups.primary_chr.proper_pairs.minq2.bam

		#B) MAPQ distribution for raw reads and proper pairs with and without chrM
		for fi in $sample.sorted.markedDups $sample.sorted.markedDups.primary_chr
		do
			samtools view $fi.bam |cut -f5 |sort -n |uniq -c > $fi.MAPQdist
			samtools view -f2 $fi.bam |cut -f5 |sort -n |uniq -c > $fi.proper_pairs.MAPQdist
		done

		#C) STATS and fragment lenght distributions for all bam files (raw, no chrM and no chrM + qulity filter) 
		for fi in $sample.sorted.markedDups $sample.sorted.markedDups.primary_chr $sample.sorted.markedDups.primary_chr.proper_pairs.minq2
		do
			samtools stats $fi.bam > $fi.stats
			samtools view $fi.bam | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > $fi.fragment_length_count.txt
		done
		#D) index the final file
		samtools index -@ $PPN $sample.sorted.markedDups.primary_chr.proper_pairs.minq2.bam

		#E) remove intermediate sam and bam files (only interested in the filtered one for downstream analyses)
		sam=$sample.bt2.sam
		rm $sam $sample.sorted.bam $sample.sorted.bam.bai $sample.sorted.markedDups.bam $sample.sorted.markedDups.bai $sample.sorted.markedDups.primary_chr.bam

		echo -e "################################\n\tsamtools view and stats done\n################################\n"
		echo `date  +'%r'`

		if $only
		then 
			exit 0
		else 
			runstep=all
		fi
	
	else
		echo `date  +'%r'`
		echo -e "################################\n\tsamtools view and stats\n################################\n"
		#A) FILTERS: proper pairs and MAPQ >= 2
	
		#A.1) retain only proper pairs and reads with minimum MAPQ>=2
		samtools view -@ $PPN -b -h -f2 -q2 $sample.sorted.markedDups.bam > $sample.sorted.markedDups.proper_pairs.minq2.bam

		#B) MAPQ distribution for raw reads and proper pairs 
		raw=$sample.sorted.markedDups
		samtools view $raw.bam |cut -f5 |sort -n |uniq -c > $raw.MAPQdist
		samtools view -f2 $raw.bam |cut -f5 |sort -n |uniq -c > $raw.proper_pairs.MAPQdist

		#C) STATS and fragment lenght distributions for both bam files (raw and proper pairs + qulity filter) 
		for fi in $sample.sorted.markedDups $sample.sorted.markedDups.proper_pairs.minq2
		do
			samtools stats $fi.bam > $fi.stats
			samtools view $fi.bam | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > $fi.fragment_length_count.txt
		done

		#D) index the final file
		samtools index -@ $PPN $sample.sorted.markedDups.proper_pairs.minq2.bam

		#E) remove intermediate sam and bam files (only interested in the filtered one for downstream analyses)
		sam=$sample.bt2.sam
		rm $sam $sample.sorted.bam $sample.sorted.bam.bai $sample.sorted.markedDups.bam $sample.sorted.markedDups.bai 
	
		echo -e "################################\n\tsamtools view and stats done\n################################\n"
		echo `date  +'%r'`
		
		if $only
		then 
			exit 0
		else 
			runstep=all
		fi
	fi
fi

#7) bigwig
step=bigWig
if [ $runstep == all ] || [ $runstep == $step ]
then
	if [[ $atacseq = true ]]
	then
		input=$sample.sorted.markedDups.primary_chr.proper_pairs.minq2.bam
	else	
		input=$sample.sorted.markedDups.proper_pairs.minq2.bam
	fi

	if ! [ -f "$input" ]; then
	die "input file $input not found"
	fi
  	
	$load_deeptools

	echo `date  +'%r'`
	echo -e "################################\n\tbamcoverage started\n################################\n"
	bamCoverage --binSize 10 --normalizeUsing CPM --bam $input -o $input.CPM.bw
	echo -e "################################\n\tbamcoverage done\n################################\n"
	echo `date  +'%r'`

	if $only
	then
		exit 0
	else
		runstep=all
	fi
fi

#8) Remove duplicates and generate bigwig
step=rmdups_bigWig
if [ $runstep == all ] || [ $runstep == $step ]
then
	if [[ $atacseq = true ]]
	then
		input=$sample.sorted.markedDups.primary_chr.proper_pairs.minq2.bam
	else	
		input=$sample.sorted.markedDups.proper_pairs.minq2.bam
	fi

	if ! [ -f "$input" ]; then
	die "input file $input not found"
	fi
  	
 
  	$load_samtools
	
	output_rd=${input/markedDups/NoDups}
	
  	echo `date  +'%r'`
	echo -e "################################\n\tsamtools rm duplicates started\n################################\n"

	samtools view -@ $PPN -F 1024 -b -h -o $output_rd $input
	samtools index -@ $PPN $output_rd
	samtools stats $output_rd > $output_rd.stats
	
	echo -e "################################\n\tsamtools rm duplicates done\n################################\n"
	echo `date  +'%r'`
  	
  	
	$load_deeptools
	input=$output_rd
	
	echo `date  +'%r'`
	echo -e "################################\n\tbamcoverage Nodups started\n################################\n"
	bamCoverage --binSize 10 --normalizeUsing CPM --bam $input -o $input.CPM.bw
	echo -e "################################\n\tbamcoverage Nodups done\n################################\n"
	echo `date  +'%r'`

	if $only
	then
		exit 0
	else
		runstep=all
	fi
fi

