## Head -------------------------------------
##
##
## Purpose of script:
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-04-18
##
## Email: lucialorenzi90@gmail.com
##
## Notes ---------------------------

# Inputs: gtf gencode
# 0) data prep: generate gene level bed file and merged exons bed file
#
# 1) compute distance to closest gene as I did before but also keep strand info
#    and also report distance with negative values if upstream relative to lncRNA
#    (-D option). For overlaps calculate the fraction of overlap
#
# 2) rules
#
#  If abs(distance) >500 -> "intergenic"
#       else
#         if abs(distance) > 0 <= 500 (not overlapping but close)
#              if pcg_strand == strand
#                 if distance < 0  -> "sense_downstream"
#                    else -> "sense_upstream"
#                 else if distance < 0 -> "divergent"
#                    else -> "convergent"
#         else (case of overlap)
#                 if frac_olA < 0.33 AND frac_olB < 0.33 (marginal overlap):
#                     if same strand
#                         if strand "+"
#                           if startA < starB -> sense_overlap_upstream
#                           else sense_overlap_downstream
#                         else (strand "-")
#                           if startA < starB -> sense_overlap_downstream
#                           else sense_overlap_upstream
#                     else (different strand)
#                         if strand "+"
#                           if startA < starB -> convergent_overlap
#                           else divergent_overlap
#                         else (strand "-")
#                           if startA < starB -> divergent_overlap
#                           else convergent_overlap
#                 else (non marginal overlap)
#                     if same strand
#                        if intronic (see details below) -> "intronic sense"
#                        else "sense_overlap"
#                     if antisense strand
#                        if intronic -> "intronic_antisense"
#                        else "antisense"

# what to do in case of multiple overlaps?? i.e if a lncRNA overlaps more
# than one PCG?? I'll classify each overlap individually and then
# prioritize: the one that has higher overlap first arrange(abs(dist),-olen)
# but also, If there are both intronic and any kind of overlapping class,
# I want to put always intronic overlaps at the end
#
##
##
##
## Setup ---------------------------

options(scipen = 999)
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
require(rtracklayer)

# load pre-defined funcs ----
source("scripts/drafts/compute_gene_non_redundant_exons_and_exonic_length.functions.R")
#source("scripts/drafts/get_genomic_regions_function.R")

# define funcs ----

# define overlapping classes and intronic classes
overlap_classes=c("sense_overlap","antisense",
                  "sense_overlap_upstream",
                  "sense_overlap_downstream",
                  "divergent_overlap","convergent_overlap")
intronic_classes=c("intronic_antisense", "intronic_sense")



# 1) this function takes a vector of parameters
# needed to classify a gene relative to a PCG.
# The distance to consider something as close
# for close but non-overlapping
# classes such as "divergent", "sense_upstream" etc
# is 2.5 kb by default (taken from previous
# definitions of divergent, e.g DOI:https://doi.org/10.1016/j.cell.2014.05.049)
classify_gene <- function(x,close_dist=2500){
  di=as.numeric(x[1])
  strQ=as.character(x[2])
  strT=as.character(x[3])
  startQ=as.integer(x[4])
  startT=as.integer(x[5])
  fracA=as.numeric(x[6])
  fracB=as.numeric(x[7])
  exon_overlap=as.logical(x[8])

  if(abs(di)>close_dist){
    return("intergenic")
  }  else{
    if(abs(di) > 0 & abs(di) <= close_dist){ #no overlap but close
      if(strQ == strT){
        if(di<0){
          return("sense_downstream")
        }else{
          return("sense_upstream")
        }
      }else{ #different strand
        if(di<0){
          return("divergent")
        }else{
          return("convergent")
        }
      }
    }else{ #overlap
      if(fracA < 0.35 & fracB < 0.35){ #marginal overlap
        if(strQ == strT){
          if(strQ == "+"){
            if(startQ < startT){
              return("sense_overlap_upstream")
            }else{
              return("sense_overlap_downstream")
            }
          }else{ #negative strand
            if(startQ < startT){
              return("sense_overlap_downstream")
            }else{
              return("sense_overlap_upstream")
            }
          }
        } else{ # marginal overlap with different strands
          if(strQ == "+"){
            if(startQ < startT){
              return("convergent_overlap")
            }else{
              return("divergent_overlap")
            }
          }else{ #negative strand
            if(startQ < startT){
              return("divergent_overlap")
            }else{
              return("convergent_overlap")
            }
          }
        }
      } else{ #non-marginal overlap
        if(strQ == strT){
          if(exon_overlap){
            return("sense_overlap")
          }else{
            return("intronic_sense")
          }
        }else{
          if(exon_overlap){
            return("antisense")
          }else{
            return("intronic_antisense")
          }
        }
      }
    }
  }
}


# 2) this function takes the reference gtf and a bed-format table with genomic coordinates for query genes
# and uses bedtools closest to generate a table with all closest/overlapping PCGs for each lncRNA
# each row of this table contains all parameters necessary to classify each query gene relative to
# the target PCG
# it then performs the classification by feeding the necessary parameters
# for each row to the function classify_gene defined above
get_gene_class <- function(ref_gtf,query_bed,N=T,k=4){
  require(bedtoolsr)
  # prepare ref tables ----
  ref_bed=trastools::get_genomic_range_by_gene(ref_gtf)
  ref_bed=ref_bed%>%dplyr::filter(seqnames%in%query_bed$seqnames)
  ref_bed=ref_bed%>%dplyr::arrange(seqnames,start)
  ref_bed_Exons=trastools::get_per_gene_merged_exons(ref_gtf)
  ref_bed_Exons=ref_bed_Exons[,c(3:5,2,6,7)]
  ref_bed_Exons=ref_bed_Exons%>%dplyr::filter(seqnames%in%query_bed$seqnames)
  ref_bed_Exons=ref_bed_Exons%>%dplyr::arrange(seqnames,start)

  query_bed <- query_bed%>%dplyr::arrange(seqnames,start)
  # distance to closest gene ----
  dist_closest=bedtoolsr::bt.closest(a=query_bed,
                          b=ref_bed,
                          D="a", # reports upstream features (relative to a) as negative distances
                          N=N, # if true requires that query and closest hit have different names
                          t = "all", # all ties are reported
                          k=k) # reports the k closest
  dist_closest=bedtoolsr::bt.overlap(i = dist_closest,cols = "2,3,8,9") #reports amount of overlap or distance (negative) for features in each row
  dist_closest$fracA=dist_closest$V14/dist_closest$V5
  dist_closest$fracB=dist_closest$V14/dist_closest$V11
  colnames(dist_closest)=c("chrQ","startQ","endQ","idQ","widthQ","strQ",
                           "chrT","startT","endT","idT","widthT","strT",
                           "dist","olen","fracA","fracB")

  dist_closest$query_target=paste(dist_closest$idQ,
                                  dist_closest$idT,
                                  sep = "_")


  # overlap with ref exons ----
  ol_exons=bedtoolsr::bt.intersect(a=query_bed,
                        b=ref_bed_Exons,wo = T) #only features with overlap are reported
  ol_exons$query_target=paste(ol_exons$V4,
                              ol_exons$V10,
                              sep = "_")

  dist_closest$exon_overlap=dist_closest$query_target%in%ol_exons$query_target

  # select cols ----
  da=dist_closest%>%
    dplyr::select(dist,strQ,strT,startQ,startT,fracA,fracB,exon_overlap)

  classif=apply(da,1,classify_gene)
  dist_closest$class=classif
  return(dist_closest)
}

# 3) this function takes ordered matching vectors of classes and corresponding
# genes for a given query and
# cleans up and returns the sorted classes/genes and the best class/gene
# in the form of a data.frame
get_best_classes_genes <- function(cl,ge) {
  if(any(cl!="intergenic")){
    ind=cl!="intergenic" # remove intergenic if other classes present
  } else ind=1 # else intergenic
  cl=cl[ind]
  ge=ge[ind]

  # in case there are both overlapping classes and intronic, move intronic
  # to the end
  if(any(cl%in%overlap_classes)&any(cl%in%intronic_classes)){
    ind=order(cl%in%intronic_classes)
    cl=cl[ind]
    ge=ge[ind]
  }
  return(data.frame(best_class=cl[1],
                    best_gene=ge[1],
                    class=paste(cl,collapse = ","),
                    gene=paste(ge,collapse = ",")))
}

# 4) this function summarises the table generated with get_gene_class
# in order to get the best class for each query gene and also
# retrieves a comma separated value of all classes relative to all overlapping/closest pcgs
# this info may be useful for record as there are cases in which a lncRNA can
# have several relevant classes that we want to keep track of
# e.g divergent and sense_overlapping relative to two closeby PCGs
# it also takes a conversion table of gene_ids to gene_names
# because gene names are more informative
summarise_classif <- function(closest_class,conversion_table){
  #first if any non-intergenic remove all intergenic
  closest_class$gene_name=conv_table$gene_name[match(closest_class$idT,
                                                     conv_table$gene_id)]
  closest_class_summ=closest_class%>%group_by(idQ)%>%
    arrange(abs(dist),-olen) %>%
    summarise(classes=paste(class,collapse = ","),
              genes=paste(gene_name,collapse = ",")
    )

  classes=strsplit(closest_class_summ$classes,split = ",")
  genes=strsplit(closest_class_summ$genes,split = ",")

  summarised_classes=mapply(function(cl,ge)get_best_classes_genes(cl,ge),
                            cl=classes,
                            ge=genes,SIMPLIFY = F)

  summarised_classes=do.call("rbind",summarised_classes)
  rownames(summarised_classes)=closest_class_summ$idQ
  return(summarised_classes)

}


# load data ----
ref_gtf_path="/home/llorenzi/Downloads/gencode.vM31.primary_assembly.annotation.gtf.gz"
ref_gtf_path="~/Descargas/merged_refs.combined.with_gene_name.gtf"
gene_biot=read.table("~/Descargas/annotated_tracking_file.updated_gene_names.txt",header = T)

# read gtf, assign gene_names as ids and extract the genes of interest (PCGs)
ref_gtf <- readGFF(ref_gtf_path)

tgn=ref_gtf$gene_name
agn=gene_biot$gene_name[match(ref_gtf$transcript_id,gene_biot$V1)]
table(tgn==agn) #not all gene_names in ref_gtf are updated
table(ref_gtf$transcript_id%in%gene_biot$V1)

ref_gtf$gene_name=gene_biot$gene_name[match(ref_gtf$transcript_id,gene_biot$V1)]
ref_gtf <- ref_gtf %>% mutate(gene_id=gene_name)

# remove RefSeq model things (XM_,XR_,XP_)
ref_gtf_tr=ref_gtf%>%filter(type=="transcript")
table(grepl("XM_|XR_",ref_gtf_tr$oId))
Model_transcripts=ref_gtf_tr$transcript_id[grepl("XM_|XR_",ref_gtf_tr$oId)]

ref_gtf <- ref_gtf %>%filter(!transcript_id%in%Model_transcripts)

table(gene_biot$simplified_gene_biotype)
PCG_transcripts <- gene_biot %>% filter(simplified_gene_biotype=="protein_coding") %>% pull(V1)

ref_gtf <- ref_gtf %>% filter(transcript_id%in%PCG_transcripts)



conv_table=ref_gtf%>%dplyr::select(gene_id,gene_name)
conv_table=conv_table[!duplicated(conv_table),]
anyNA(conv_table$gene_name)


# run code ----

# filtered assembled genes
query_data_path="outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.tsv"
query_data <- read.table(query_data_path,header = T)
query_bed <- trastools::get_genomic_range_by_gene(query_data,by = "gene_name")
closest_class_to_PCG=get_gene_class(ref_gtf,query_bed,N = F)
table(closest_class_to_PCG$idT%in%conv_table$gene_id)
#make sure all gene_ids have a gene_name
anyNA(conv_table$gene_name)
class_summary=summarise_classif(closest_class = closest_class_to_PCG,conversion_table = conv_table)

colnames(class_summary) <- c("best_classif_to_PCG","best_closest_PCG","classif_to_PCG","closest_PCG")
class_summary$gene_name <- rownames(class_summary)


# add info to tables ----
query_data <- left_join(query_data, class_summary)


# read data on gene level directly
gene_level_data_path="outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.tsv"
gene_level_data=read.table(gene_level_data_path,header = T)
#source("scripts/source_all_functions.R")

gene_level_data <- gene_level_data %>% filter(pass_filter)
gene_level_data <- left_join(gene_level_data, class_summary)
gene_level_data <- gene_level_data %>% mutate(biotype = ifelse(is.na(biotype),"potNovel",biotype))


query_data$gene_biotype <- gene_level_data$biotype[match(query_data$gene_name,
                                                         gene_level_data$gene_name)]
# write results ----
out_path_transcript_level=gsub(".tsv",".gene_classif.tsv",query_data_path)
write.table(query_data,out_path_transcript_level,quote = F,row.names = F,sep = "\t")

out_path_gene_level=gsub(".tsv",".filtered.gene_classif.tsv",gene_level_data_path)
write.table(gene_level_data,out_path_gene_level,quote = F,row.names = F,sep = "\t")
