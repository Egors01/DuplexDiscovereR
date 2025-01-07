#' This script shows how the the example data used for test and vignettes was generated.
#' It is not intended to be run 'as-is'. To replicate the process, one should
#' start with downloading the raw duplex read read data fro i.e SRA and proceed
#' with the downstream steps:

#' For creating this example data, we used single replicate of the SPLASH ES
#' cells libraries from SRA https://www.ncbi.nlm.nih.gov/sra/?term=SRR3404943.
#' Further details on SPLASH can be found here
#' https://csb5.github.io/splash/
#' Reads were adidtionally pre-processed with fastp using default parameters.
#' SPLASH libraries were sequenced in paied-end mode,
#' but stored as merged (vie seqprep)to single-end in in SRA
#' Following STAR command was used to map raw reads and produce
#' Chimeric.out.junction file and ReadsPerGene.out.tab
##' -------------------------------
##'   ##### Final effective command line:
##'   STAR   --runThreadN 8   --genomeDir ./hg38_masked_index
##' --genomeLoad NoSharedMemory   --readFilesIn ./SRR3404943.filt.fastq
##' --outFileNamePrefix ./SPLASH_ES_1/   --outSAMtype BAM   Unsorted
##' --outSAMattributes All
##' --outSAMunmapped Within
##' --outSJfilterReads All
##' --alignIntronMin 1
##' --alignIntronMax 10
##' --chimSegmentMin 15
##' --chimScoreMin 1
##' --chimScoreDropMax 30
##' --chimScoreJunctionNonGTAG 0
##' --chimJunctionOverhangMin 10
##' --chimOutType Junctions
##' --chimSegmentReadGapMax 3
##' --chimMultimapNmax 5
##' --chimNonchimScoreDropMin 20
##' --chimOutJunctionFormat 1
##' --quantMode GeneCounts
##' ----------------------------------------
##' After calling STAR, one can follow the steps below, modifying the
##' paths to the user location
##'
##' Use STAR Chimeric.out.junction file as the input.
##' For the test data, we use chr 22 and only whole genes as the features.
##'
library(tidyverse)
library(InteractionSet)
stardf = readr::read_tsv("./SPLASH_ES_1/Chimeric.out.junction")
RNADuplexesRawChimSTAR <- stardf %>% dplyr::filter(chr_donorA==chr_acceptorB,
                                                   chr_donorA == 'chr22')

# Import gene read counts
df_counts = read_tsv('./SPLASH_ES_1/ReadsPerGene.out.tab',skip = 4,col_names = F)
colnames(df_counts) = c('id','s_strand','r_strand','unstr')
RNADuplexesGeneCounts = df_counts[c('id','unstr')]

##' Import gencode human refernce annotation
##' requires downloaded GENCODE or other annotation.
##' https://www.gencodegenes.org/human/release_44.html or via Annotation Hub
##' As thsi test data inlcudes only chr 22, we reduce the annotatoin accordingly.

library(rtracklayer)
gtf_gnc = import.gff('./gencode.v44.annotation.gtf')
SampleGeneAnnoGR = gtf_gnc[seqlevels(gtf_gnc) %in% 'chr22']
SampleGeneAnnoGR = SampleGeneAnnoGR[SampleGeneAnnoGR$type=='gene']

##' To create the example splice junction GRanges object, one can use the reference
##' annotation, as we show below.
##' However, for working with the real-life duplex probing data, we
##' recommend to use the custom file containing splice
##' junctions discovered from mapping the RNA-seq reads of the appropriate
##' cell lines. Novel junctions can be found by STAR run and extracted from the
##' output SJ.out.tab

library(GenomicFeatures)
library(SGSeq)
# create txdb from the same gencode annotation
txdb_v44 <- makeTxDbFromGFF("./gencode.v44.annotation.gtf")
saveDb(txdb_v44,"./txdb_v44.txdb")
tx_gtf<-AnnotationDbi::loadDb("./txdb_v44.txdb")
txf <- convertToTxFeatures(tx_gtf)
SampleSpliceJncGR<-GRanges(txf[txf@type=='J'])

# call the pipeline to produce clustered and unclustered reads samples
library(DuplexDiscovereR)
res <- run_duplexdisco(
  df = RNADuplexesRawChimSTAR,
  junctions_gr = SampleSpliceJncGR,
  anno_gr = SampleGeneAnnoGR,
  df_counts = RNADuplexesGeneCounts,
  sample_name = 'test',
  lib_type = "SE",
  table_type = "STAR"
)
# Two-arm chimeric reads, which can be represented in the GInteraction object
RNADuplexSampleGI = res$gi_reads
mcols(RNADuplexSampleGI) = as_tibble(data.frame(mcols(res$gi_reads))) %>%
  dplyr::select(n_reads,n_reads_dg,dg_id,score,readname,map_type,junction_type,
                cigar_alnA ,cigar_alnB) %>% data.frame()
# Reads clustered into duplex groups
RNADuplexSampleClustReads= res$gi_reads
# Reads clustered and collapsed into duplex groups
RNADuplexSampleDGs = res$gi_clusters
# Reads which can be represented as 2-arm chimeras
# We leave only required BEDPE-specific columns in order defined in a format defintion
# see specification at https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format
# We use 'flag field' as the placeholder for the alignment score
RNADuplexesRawBed = DuplexDiscovereR::MakeDfFromGi(RNADuplexSampleGI) %>%
dplyr::select(chromA ,startA,endA,chromB ,startB,endB,readname,score,strandA,strandB)


save(RNADuplexSampleGI,RNADuplexSampleDGs,RNADuplexesRawBed,
     RNADuplexesRawChimSTAR,RNADuplexesGeneCounts,SampleSpliceJncGR,RNADuplexSampleClustReads,
     SampleGeneAnnoGR,file = './data/RNADuplexesSampleData.rda', compress = 'xz')


# save as individual files for tests
write_tsv(RNADuplexesRawBed,'./inst/extdata/test_SPLASH_DuplexesRaw.bedpe',col_names = F)
write_tsv(RNADuplexesRawChimSTAR,'./inst/extdata/test_SPLASH_Chimeric.out.junction',col_names = T)
write_tsv(RNADuplexesGeneCounts,'./inst/extdata/test_ReadsPerGene.out.tab',col_names = T)
export.bed(SampleSpliceJncGR,'./inst/extdata/test_chr22SpliceJunctions.bed')
export.gff(SampleGeneAnnoGR,'./inst/extdata/test_chr22Anno.gtf')



# stardf <- DuplexDiscovereR::run_preproc(stardf,table_type = 'STAR',library_type = 'SE')
# # select
# twoarm <- stardf %>% filter(map_type=='2arm')
# twoarm_gi = MakeGiFromDf(twoarm)
# load('/data/meyer/egor/duplex_workflow/public/DuplexDiscovereR/data/example_data.rda')
# gr_region = GRanges(IRanges(start=17504924,end=17524195),seqnames='chr22',strand = '+')
# add_reads = subsetByOverlaps(twoarm_gi,gr_region+50000)$readname
# add_df = stardf_raw[stardf_raw$read_name %in% add_reads,]
#
# stardf %>% filter(readname %in% example_reads_star$read_name)
# new_removed = example_reads_star[-sample(c(1:5000),size=81),]
# new_test = rbind(new_removed,add_df)
# write_tsv(new_test,'/data/meyer/egor/duplex_workflow/public/DuplexDiscovereR/inst/extdata/sample_Chimeric.out.junction')
# run_duplexdisco()
# setwd("/data/meyer/egor/duplex_workflow/public/DuplexDiscovereR/")
#
#
#
