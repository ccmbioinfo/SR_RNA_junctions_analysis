suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))


option_list <- list(
  make_option(c("-s", "--samplename"), type = "character", help="samplename"),

  # these two need to be less hardcoded
  make_option(c("-d", '--directory'), default = "/hpf/projects/dig2/huayun/processed",help="directory to find results in"),
  make_option(c("-o", "--outname"), type = "character", default = ".qcsum.txt", help = "output filenames after sampleID, default to qcsum.txt"),

  make_option("--star_file", type = "character", help = "STAR Log.final.out file"),
  make_option("--bam_file", type = "character", help = "duplicates marked bam file"),
  make_option("--picard_file", type = "character", help = "picard rnaseq metrics file"),
  make_option("--duplicate_file", type = "character", help = "picard mark duplicates metrics file"),
  make_option("--insert_file", type = "character", help = "picard insert size metrics file" ),
  make_option("--alignment_file" , type = "character", help = "picard alignment metrics file"),
  make_option("--qchist_file" , type = "character", help = "cal_q30 results"),
  make_option("--coverage_file" , type = "character", help = "coverage stat file from get_cov rule"),
  make_option("--expr_file" , type = "character", help = "rsem gene expresion file"),
  make_option("--junc_file",  type = "character", help = "STAR SJ.out.tab file" ),
  make_option("--ercc_file" , type = "character", help = "ercc metrics file"),
  make_option("--sirv_file" , type = "character", help = "sirv metrics files"),
  make_option("--rnaseq_file" , type = "character", help = "rnaseqc results"),
  make_option("--seq_len", type="integer", help="read_length", default=150),
  make_option("--tx_len", type="integer", help="transcriptome_length", default=111579398)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

samplename <- opt$samplename
results_dir <- opt$directory
star_file <- opt$star_file
bam_file <- opt$bam_file
picard_file <- opt$picard_file
duplicate_file <- opt$duplicate_file
insert_file <- opt$insert_file
alignment_file <- opt$alignment_file
qchist_file <- opt$qchist_file
coverage_file <- opt$coverage_file
expr_file <- opt$expr_file
junc_file <- opt$junc_file
ercc_file <- opt$ercc_file
sirv_file <- opt$sirv_file
rnaseq_file <- opt$rnaseq_file
seq_length <- opt$seq_len
tr_length <- opt$tx_len

# Get total and unique mapped reads
# since the star file is a proper dependency in the snakefile if it needs to be present for the pipeline to continue
# I think this is a more robust way of doing things. We should let snakemake make the global decisions about dependencies
message("Getting alignment stats")
star_metric <- read.table(star_file, sep="\t",fill = T, stringsAsFactors = F)
aln_data <- data.frame("SAMPLE_ID" = samplename,
                       "TOTAL_READS" = as.numeric(star_metric[5,2]),
                       "UNIQ_ALIGNED" = as.numeric(star_metric[8,2]),
                       stringsAsFactors = F)

# Get percent MT
# count primary alignments to MT
message("Calculating MT percentage")
MT_cmd <- paste0("samtools view -F 256 -c ", bam_file," MT")
MT_counts <- as.numeric(system(MT_cmd, intern = T))
# calculate percentage of MT reads / total read pairs * 2 
MT_PCT <- data.frame("MT_PCT" = round(100*MT_counts/(aln_data$TOTAL_READS*2),2),stringsAsFactors = F)

# get metrics from picard rnaseq metric output
message("Getting metrices from picard RNAseq metrics")
picard_data <- read.table(picard_file, fill = T, as.is = T, stringsAsFactors = F)
rnaseq_metric <- picard_data[2,]
names(rnaseq_metric) <- picard_data[1,]
used_cols <- c("MEDIAN_5PRIME_TO_3PRIME_BIAS",
               "PCT_RIBOSOMAL_BASES")
used_rnaseq_metric <- rnaseq_metric[1, used_cols] 
used_rnaseq_metric[, grepl("PCT", names(used_rnaseq_metric))] <- round(100*as.numeric(used_rnaseq_metric[, grepl("PCT", names(used_rnaseq_metric))]),2)
used_rnaseq_metric$MEDIAN_5PRIME_TO_3PRIME_BIAS <- round(as.numeric(used_rnaseq_metric$MEDIAN_5PRIME_TO_3PRIME_BIAS), 2)


# rnaseq metrics from rnaseqc
message("Getting metrics from RNAseQC")
rnaseq <- read.table(rnaseq_file, sep = "\t", stringsAsFactors = F, quote = "", header = T, col.names = c("Metrics", "Value")) %>% 
  filter(Metrics %in% c("Exonic Rate", "Intronic Rate", "Intergenic Rate", "Ambiguous Alignment Rate")) %>% 
  mutate(Value = round(100*Value, 2)) %>% 
  spread(Metrics, Value) %>% 
  dplyr::select("PCT_EXON" = "Exonic Rate", 
                "PCT_INTRON" = "Intronic Rate", 
                "PCT_INTERGENIC" = "Intergenic Rate", 
                "AMBIGUOUS_ALN" = "Ambiguous Alignment Rate")


# # get average coverage across gene bodies
# gene_coverage <- picard_data[3:nrow(picard_data), 1:2]
# names(gene_coverage) <- c(as.character(gene_coverage[1,1]), as.character(gene_coverage[1,2]))
# gene_coverage <- gene_coverage[-1,]

# get duplication metrics
message("Getting duplication metrics and insert sizes")
dup_data <- read.table(duplicate_file, fill = T, as.is = T, stringsAsFactors = F, header = T)[1, c("PERCENT_DUPLICATION", "ESTIMATED_LIBRARY_SIZE")]
dup_data[, "PERCENT_DUPLICATION"] <- round(100*dup_data[, "PERCENT_DUPLICATION"],2)

# insert size metrics
insert_data <- read.table(insert_file, fill = T, as.is = T, stringsAsFactors = F, header = T)[1, c("MEAN_INSERT_SIZE", "STANDARD_DEVIATION")]

# # get alignment metrics
# aln_data <- read.table(alignment_file, fill = T, as.is = T, stringsAsFactors = F, header = T)
# 
# aln_data <- data.frame("SAMPLE_ID" = samplename, "TOTAL_READS" = format(round(aln_data[1,"TOTAL_READS"]), nsmall = 0), "PCT_ALIGNED" = round(100*aln_data[3, "PCT_PF_READS_ALIGNED"], 2), stringsAsFactors = F)
# aln_data$COVERAGE <- round(2 * as.numeric(aln_data$TOTAL_READS) * seq_length/tr_length)
# MT_PCT <- data.frame("MT_PCT" = round(100*MT_counts/as.numeric(aln_data["TOTAL_READS"]),2), stringsAsFactors = F)

# q30
message("Getting q30 metrics")
qsum <- read.table(qchist_file, header = T, comment.char = "", as.is = T, stringsAsFactors = F)
qsum <- as.matrix(qsum)
R1_ave <- round(sum(qsum[,"X.Quality"] * qsum[,"count1"])/sum(qsum[,"count1"]), 2)
R2_ave <- round(sum(qsum[,"X.Quality"] * qsum[,"count2"])/sum(qsum[,"count2"]), 2)
R1_q30 <- round(100*sum(qsum[qsum[, "X.Quality"]>30, "fraction1"]), 2)
R2_q30 <- round(100*sum(qsum[qsum[, "X.Quality"]>30, "fraction2"]), 2)

q30_metrics <- data.frame("R1_AVE_Q" = R1_ave, "R2_AVE_Q" = R2_ave, "PCT_R1_Q30" = R1_q30, "PCT_R2_Q30" = R2_q30, stringsAsFactors = F)

# get coverage data
message("Getting coverage data")
coverage <- read.table(coverage_file, as.is = T, stringsAsFactors = F)
names(coverage) <- c("AVE_COV", "PCT_1X_COV", "PCT_10X_COV", "PCT_50X_COV")

# get ERCC
message("Summarize ERCC expression")
ercc <- read.table(ercc_file, as.is = T, stringsAsFactors = F, header = T)

if(sum(ercc$TPM) == 0){
  ercc_result <- matrix(NA, ncol = 4, nrow = 1)
  colnames(ercc_result) <- c("ERCC_TPM_SUM", "ERCC_TPM_0", "ERCC_TPM_1","ERCC_R2")
} else{
  detected <- as.numeric(table(ercc$TPM > 0)[2])
  detected_tpm1 <- as.numeric(table(ercc$TPM >= 1)[2])
  model <- lm(log2(TPM) ~ log2(conc), data = subset(ercc, !is.infinite(log2(ercc$TPM))))
  rsquared <- round(summary(model)$r.squared, 2)
  ercc_result <- data.frame("ERCC_TPM_SUM" = round(sum(ercc$TPM)), "ERCC_TPM_0" = detected, "ERCC_TPM_1" = detected_tpm1, "ERCC_R2" = rsquared, stringsAsFactors = F)
}


# analyze SIRVs
message("Analyze SIVRs expression")
sirv <- read.table(sirv_file, as.is = T, stringsAsFactors = F, header = T)
if(sum(sirv$TPM) == 0){
  sirv_sum <- matrix(NA, ncol = 8, nrow = 1)
  colnames(sirv_sum) <- c("SIRV_TPM_SUM", "SIRV1", "SIRV2", "SIRV3","SIRV4","SIRV5","SIRV6","SIRV7")
} else{
  sirv_sum <- sirv %>%
  group_by(gene_id) %>%
  #summarise(CV = sd(TPM)/mean(TPM))
  mutate(perc_TPM = 100*TPM/sum(TPM)) %>%
  summarise(sd_perc = sd(perc_TPM)) %>%
  mutate(sd_perc = round(sd_perc, 2)) %>% 
  spread(gene_id, sd_perc) %>% 
  as.data.frame() %>% 
  mutate(SIRV_TPM_SUM = round(sum(sirv$TPM))) %>% 
  dplyr::select(SIRV_TPM_SUM, everything())
}

# number of junctions detected 
juncs <- read.table(junc_file, stringsAsFactors = F, as.is = F)
NUM_JUNC <- table(juncs[,7] >= 5)[2]
rm(juncs)

# expression level summary
expr <- read.table(expr_file, header = T, stringsAsFactors = F, as.is = F)
allTPM <- expr[!grepl("ERCC|SIRV", expr$gene_id),"TPM"]
expr_sum <- data.frame("TPM_0" = table(allTPM > 0)[2],
                      "TPM_1" = table(allTPM > 1)[2],
                      "TPM_10" = table(allTPM > 10)[2],
                      "TPM_50" = table(allTPM > 50)[2],
                      "NUM_JUNCS" = NUM_JUNC,
                      stringsAsFactors = F)
rm(expr)
                      


# put together into one line
qc_data_list <- list(aln_data, q30_metrics, dup_data, coverage, rnaseq, used_rnaseq_metric,MT_PCT, expr_sum, ercc_result, sirv_sum)
qc_data_list <- qc_data_list[!sapply(qc_data_list, is.null)]

qc_sum <- do.call("cbind", qc_data_list)
#row.names(qc_sum) <- samplename

# output
outname <- file.path(results_dir, paste0(samplename, opt$outname))
write.table(qc_sum, outname, row.names = F, col.names = T, quote = F, sep = "\t")
