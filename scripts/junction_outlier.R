suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(tibble))

options(dplyr.summarise.inform = FALSE)#this is to avoid excess messaging

# db is a connection not a path
# running an inequality query results in a full table scan makes the junc_id index useless,
# the difference in time 100 mins vs 2 mins is worth the extra memory consumption
get_control_junctions<-function(db, tissue, cutoff=5){
  query<-paste0("select junc_id,sample_id,uniq_map from junc_to_sample where sample_id in (select sample_id from samples where tissue='", tissue, "')")
  tissue_junctions<-dbGetQuery(db, query)
  tissue_junctions<-tissue_junctions %>%
    filter(uniq_map >= cutoff)
  tissue_junctions<-dcast(tissue_junctions, junc_id~sample_id, value.var="uniq_map", fill=0)
  return(tissue_junctions)
}

find_junction_file <- function(sample_id){
  # this function finds .SJ file for control samples. Considering results directory maybe named as "sampleID_xxx"
  #sample_id <- as.character(sample_id)
  junc_file <- file.path(control_dir, sample_id, paste0(sample_id, "SJ.out.tab"))
  if(!file.exists(junc_file)){
    sample_dir <- dir(path = control_dir,
                      pattern = glob2rx(paste0(sample_id, "_*")),
                      full.name = T)
    junc_file <- list.files(path = sample_dir,
                            pattern = glob2rx("*SJ.out.tab"),
                            full.names = T)
  }
  if(length(junc_file) != 1){
    stop(paste0("No or multiple SJ.out.tab file found for sample ", sample_id))
  }else{
    return(junc_file)
  }
}


get_junction<-function(SJ_file, chrs, cutoff=5){
  # this function reads the SJ.out.tab file, filters by read depth and formats
  if(!file.exists(SJ_file)){
    stop(paste0(SJ_file, " does not exist"))
  }
  file<-read.table(SJ_file, header=F, stringsAsFactors = F, sep = "\t",
                   col.names = c("chr", "start", "end", "strand",
                                 "motif", "annotated", "uniq_map", "multi_map",
                                 "max_ohang"))[-6] %>%
    mutate(strand = case_when(strand == 1 ~ "+",
                              strand == 2 ~ "-",
                              TRUE ~ ".")) %>%
    filter(uniq_map >= cutoff, chr %in% chrs) %>%
    unite(chr, start, end, strand, sep = "_", col = "junc_id") %>%
    dplyr::select(junc_id, uniq_map)
  return(file)
}

get_zscores <- function(control_junctions, sample_junctions, sample_id = as.character(opt$samplename), ncores = 1, perc_zscore_cutoff = -1.5){

  # rename sample junction column name to sample_id
  names(sample_junctions)[2] <- sample_id

  # merge control and sample junctions
  all_junctions <- full_join(control_junctions, sample_junctions, by="junc_id") %>%
    separate(junc_id, sep = "_", into = c("chr", "start", "end", "strand"), remove = F) %>%
   # rename(!!sample_id := "uniq_map") %>%
    column_to_rownames("junc_id")


  # change any 0 to NA so to calculate mean counts with only samples in which the junction is detected
  all_junctions[all_junctions == 0] <- NA
  count_means <- rowMeans(all_junctions[, -c(1:4, ncol(all_junctions))], na.rm = T)

  # replace NA values with 0
  all_junctions[is.na(all_junctions)] <- 0

  donor_sum <- suppressMessages(all_junctions %>%
    arrange(chr, start, strand) %>%
    group_by(chr, start, strand) %>%
    mutate_if(is.numeric, sum) %>%
    unite("id", c("chr", "start", "end", "strand")) %>%
    column_to_rownames("id"))


  # calculate read sums of junctions with shared acceptor sites
  acceptor_sum <- suppressMessages(all_junctions %>%
    arrange(chr, end, strand) %>%
    group_by(chr, end, strand) %>%
    mutate_if(is.numeric, sum) %>%
    unite("id", c("chr", "start", "end", "strand")) %>%
    column_to_rownames("id"))


  temp_junctions <- all_junctions %>%
    dplyr::select(-c("chr", "start", "end", "strand"))
 # ratios <- all_junctions/total_counts

  ratios <- temp_junctions/(donor_sum[row.names(temp_junctions),] + acceptor_sum[row.names(temp_junctions),] - temp_junctions)

  # calculated relative to control data (sample is not included in calculation)
  ratio_means <- rowMeans(ratios[,-ncol(ratios)], na.rm = T)
  ratio_sds <- apply(ratios[,-ncol(ratios)], 1, sd, na.rm = T)
  zscores <- abs((ratios[[sample_id]]- ratio_means)/ratio_sds)

  junction_stats <- data.frame(
    junc_id = names(zscores),
    gtex_frac_detected = rowSums(temp_junctions[,-ncol(temp_junctions)] > 0, na.rm = T)/(ncol(temp_junctions) - 1),
    mean_gtex_ratio = ratio_means,
    sd_gtex_ratio = ratio_sds,
    mean_gtex_count = rowMeans(temp_junctions[, -ncol(temp_junctions)], na.rm = T),
    sample_ratio = ratios[ ,sample_id],
    sample_count = temp_junctions[, sample_id],
    zscore = unname(zscores)
  ) %>%
    separate(junc_id, sep = "_", into = c("chrom", "start", "end", "strand")) %>%
    mutate(start = as.integer(start), end = as.integer(end)) %>%
    dplyr::select(chrom, start, end, strand, everything())

  #######################################################
  # calculate additional metrics on junctions with ratio > 0.95
  # associate junctions with genes and calculate percent of reads mapping to each junction
  sample_names_all <- colnames(all_junctions)[-c(1:4)]
  
  # annotate junctions to genes. Very fast step
  junc_gr <- makeGRangesFromDataFrame(all_junctions[, c("chr","start","end","strand")], keep.extra.columns = T)
  olaps <- findOverlaps(junc_gr, genes_gr, type="any", ignore.strand=F, select="all")
  junctions_annot <- all_junctions[olaps@from, c("chr","start","end","strand")]
  gene_names <- genes_gr$gene_id[olaps@to]
  junctions_annot$gene_name <- gene_names
  junctions_annot$gene_width <- width(genes_gr)[olaps@to]
  
  # select only junctions with both gtex and sample ratio of 1
  #intron_juncs <- names(ratio_means)[which(ratio_means >= 0.95 & (ratios[[sample_id]] >= 0.95|is.nan(ratios[[sample_id]])))]
  intron_juncs <- names(ratio_means)[which(ratio_means >= 0.95 & (ratios[[sample_id]] >= 0.95))]
  
  # if there is no junctions annotated to genes or no intron junctions to analyze
  if(nrow(junctions_annot) != 0 & length(intron_juncs) != 0){
    # calculate total number of junction reads mapping to junction-associated genes
    junction_with_genes <- left_join(all_junctions, junctions_annot, by = c("chr", "start", "end", "strand")) %>% 
      unite("id", c("chr", "start", "end", "strand")) %>% 
      filter(!is.na(gene_name) & id %in% intron_juncs) %>% 
      group_by(gene_name, gene_width) %>% 
      mutate_at(sample_names_all, sum) %>% 
      unite("id2", c("id","gene_name","gene_width"), remove = F) 
    
    # calculate percent of reads mapping to this junction/total reads mapping to the gene
    junction_perc <- temp_junctions[junction_with_genes$id,]/junction_with_genes[,sample_names_all]
    row.names(junction_perc) <- junction_with_genes$id2
    
    # calculating z score based on the percentage 
    perc_means <- rowMeans(junction_perc[,-ncol(junction_perc)], na.rm = T)
    perc_sds <- apply(junction_perc[,-ncol(junction_perc)], 1, sd, na.rm = T)
    perc_fc <- junction_perc[[sample_id]]/perc_means
    perc_zscores <- (junction_perc[[sample_id]]- perc_means)/perc_sds
    
    # selecting only junctions to report
    juncs_oi <-names(perc_fc)[which(perc_fc < 0.5 & perc_zscores < perc_zscore_cutoff)]
    if(length(juncs_oi) == 0){
      return(list("ratio" = junction_stats, "percent" = NA))
    } else{
      # assemble output; use the same set of column names 
      mean_gtex_count = rowMeans(temp_junctions[, -ncol(temp_junctions)], na.rm = T)
      gtex_frac_detected = rowSums(temp_junctions[,-ncol(temp_junctions)] > 0, na.rm = T)/(ncol(temp_junctions) - 1)
      
      junction_stats_single <- data.frame(juncs = juncs_oi,
                                          zscore = perc_zscores[juncs_oi],
                                          sample_ratio = junction_perc[juncs_oi, sample_id],
                                          mean_gtex_ratio = perc_means[juncs_oi],
                                          sd_gtex_ratio = perc_sds[juncs_oi],
                                          stringsAsFactors = F) %>% 
        separate(juncs, sep = "_", into = c("chrom","start","end","strand","gene_name","gene_width")) %>% 
        mutate(start = as.integer(start), end = as.integer(end)) %>% 
        unite(id, chrom, start, end, strand, remove = F)
      
      junction_stats_single$sample_count = temp_junctions[junction_stats_single$id, sample_id]
      junction_stats_single$mean_gtex_count = mean_gtex_count[junction_stats_single$id]
      junction_stats_single$gtex_frac_detected = gtex_frac_detected[junction_stats_single$id]
      
      # reorder columns and filter for mean control counts > 10
      junction_stats_single <- junction_stats_single[, c(colnames(junction_stats), "gene_name","gene_width")] %>% 
        filter(mean_gtex_count > 10)
      
      # return the results based on either junction ratios or percent reads mapping to junction as a list 
      return(list("ratio" = junction_stats, "percent" = junction_stats_single))
    }

  } else{
    return(list("ratio" = junction_stats, "percent" = NA))
  }

}


reformat_txdb_exon <- function(txdb_used, get_longest = TRUE, selected_txs = NULL, return_longest = FALSE){
  # function to get exons from a txdb object. Options are;
  # 1. When 'get_longest' is set to TRUE, get the transcript from each gene that has the largest number of exons and name exons based on this transcript
  #    Used for ens37 genes
  # 2. Use selected transcripts when 'selected_txs' is a vector of preferred transcripts
  #    Used for getting exons for hgmd preferred refseq transcripts
  if(get_longest){
    # get num of exons for each transcript
    numExons <- elementNROWS(exonsBy(txdb_used, by = "tx", use.name = T))

    # filtering for the tx with the largest num of exons, if tied, get the longest transcript, if both tied, just get the first one
    txs <- stack(transcriptsBy(txdb_used, by = "gene"), "gene") %>%
      as.data.frame() %>%
      mutate(nExons = numExons[tx_name]) %>%
      group_by(gene) %>%
     # slice_max(order_by = c(nExons, width), n = 1)
      arrange(desc(nExons), desc(width)) %>%
      slice_head(n =1)

    # if return_longest, return a vector with gene_name:number of exons
    if(return_longest){
      return(txs[,c("gene","nExons")] %>% deframe())
    }

    # match between transcript and gene names
    tx_gene_match <- txs[, c("tx_name","gene")] %>%
      deframe()
    used_txs <- txs$tx_name # a list of used transcripts, in this case, longest from each gene
  }

  if(!is.null(selected_txs)){
    # if a list of transcripts are provided, for example, hgmd preferred ones
    used_txs <- selected_txs
    tx_gene_match <- stack(transcriptsBy(txdb_used, by = "gene"), "gene") %>%
      as.data.frame() %>%
      filter(tx_name %in% used_txs) %>%
      dplyr::select(tx_name, gene) %>%
      deframe()
  }

  # get exons per transcript, filter for either the longest transcripts or the preferred transcripts if provided
  exons_by_tx <- exonsBy(txdb_used, "tx", use.name = T)
  exons_by_tx <- stack(exons_by_tx[names(exons_by_tx) %in% used_txs],"tx_name")
  # add "exon" to exon_rank to make exon ids in the formats of exon1, exon2; add gene ids back
  mcols(exons_by_tx)$gene <- tx_gene_match[as.character(exons_by_tx$tx_name)]
  mcols(exons_by_tx)$exon_rank <- paste0("exon",mcols(exons_by_tx)$exon_rank)

  # return grange object of exons from preferred transcripts per gene, with exon ids and gene ids
  return(exons_by_tx)
}

annotate_with_exons <- function(junctions, txdb_used = txdb, exon_col_name = "exon_ids", used_gene_name_col = "gene_name", used_txdb_col = "gene", ncores = 1,...){
  # make a unique list of junctions into granges
  junc_gr <- makeGRangesFromDataFrame(unique(junctions[,c("chrom","start","end","strand")]))

  print(paste(Sys.time(),"formatting exons"))
  # obtain and reformat exons from ensembl txdb object for the longest transcript of the gene
  exon_all_gr <- reformat_txdb_exon(txdb_used,...)

  print(paste(Sys.time(),"overlapping"))
  # perform overlap
  ol <- findOverlaps(junc_gr, exon_all_gr, maxgap = 0)

  # get the data frame version of the junction granges
  junc_df <- as.data.frame(junc_gr, stringsAsFactors = F)
  colnames(junc_df)[1] <- "chrom"

  # reformat to join ori junctions and annotated exons, return this data
  junc_withExon <- cbind(junc_df[from(ol),],
                         mcols(exon_all_gr)[to(ol),c(used_txdb_col,"exon_rank")]) %>%
    as.data.frame() %>%
    group_by_at(vars(-exon_rank)) %>%
    summarise(!!exon_col_name := paste(sort(exon_rank), collapse = ","), .groups = "drop") %>%
    right_join(junctions,  by = c("chrom","start","end","strand", setNames(used_gene_name_col, used_txdb_col))) %>%
    dplyr::select(-width)
  colnames(junc_withExon)[colnames(junc_withExon) == used_txdb_col] <- used_gene_name_col
  return(junc_withExon)
}


annotate_results<-function(junctions, genes_gr, gtf_junctions, db, tissue, ncores = opt$ncores, anno_to_gene = T, add_hgmd_exons = T,HPO = NULL){
  # If the junction is annotated in the gtf file or not
  junctions<-left_join(junctions, gtf_junctions, by=c("chrom", "start", "end", "strand"))
  junctions$annotated[is.na(junctions$annotated)] <- F

  if(anno_to_gene){
    # Annotate junctions with the genes they associated with
    junc_gr <- makeGRangesFromDataFrame(junctions, keep.extra.columns = T)
    olaps <- findOverlaps(junc_gr, genes_gr, type="any", ignore.strand=F, select="all")
    junctions_annot <- junctions[olaps@from,]
    gene_names <- genes_gr$gene_id[olaps@to]
    junctions_annot$gene_name <- gene_names
    junctions_annot$gene_width <- width(genes_gr)[olaps@to]
    print("Associated junctions to genes")
    suppressMessages(junctions<-left_join(junctions, junctions_annot))
  }

  # annotate with disease databases
  annotations<-dbGetQuery(db, paste0("select * from annotations where gene_name in (",
                                      paste(paste0("'", unique(na.omit(junctions$gene_name)), "'"), collapse = ","), ")"))
  #annotations <- readRDS("/home/huayun/dig2/acelik_analysis/clinical_pipeline/genomes/hg37/gene_annotations_withAlias.rds")
  
  if("index" %in% colnames(annotations)){
    annotations <- subset(annotations, select=-c(gene_id,index))#drop the gene id and index columns
  }
  
  print("Add gene annotations")
  junctions <- unique(left_join(junctions, annotations, by="gene_name"))

  #  # annotate with exon ids
  all_genes <- unique(na.omit(junctions$gene_name))

  print("Add exon ids")
  suppressMessages(junctions <- unique(annotate_with_exons(junctions, txdb_used = txdb, ncores = ncores)))

  #suppressMessages(junctions <- left_join(junctions, anno_exons) )

  
  if(add_hgmd_exons){
    # add refseq exon ids from hgmd preferred transcript
    print("Add hgmd transcript exon ids")
    suppressWarnings(junctions <- annotate_with_exons(junctions, txdb_used = ref_txdb, exon_col_name = "hgmd_tx_exon_ids", used_gene_name_col = "hgmd_tx", used_txdb_col = "tx_name", get_longest = FALSE, selected_txs = annotations$hgmd_tx[!is.na(annotations$hgmd_tx)],ncores = ncores))
  } else{
    junctions <- mutate(junctions, hgmd_tx_exon_ids = NA, hgmd_tx = NA)
  }
  

  num_exons <- reformat_txdb_exon(txdb, return_longest = T)
  junctions$num_exons <- num_exons[junctions$gene_name]
  # summarize number of junctions detected per gene
  junctions_gene <- junctions[!is.na(junctions$gene_name), c("gene_name", "chrom")] %>%
    group_by(gene_name) %>%
    summarise(total_juncs = length(chrom), .groups = "drop")
  junctions <- left_join(junctions, junctions_gene, by = "gene_name")

  # annotate with expression
  expressions<-dbGetQuery(db, paste0("select b.gene_name,a.mean_tpm, a.median_tpm, a.sd_tpm from gene_averages a
                                     left join gene_info b on a.gene_id=b.gene_id where tissue='",
                                     tissue, "' and a.gene_id in (",
                                     paste("select gene_id from gene_info where gene_name in (",
                                           paste(paste0("'", all_genes, "'"), collapse = ","), ")"), ")"))

  print("Add gene expression levels")
  suppressMessages(junctions <- left_join(junctions, expressions, by=c("gene_name")))

  print("Checking HPO terms")
  if(is.null(HPO)){
    junctions$relevant_HPO <- NA
  } else{
    junctions$relevant_HPO <- sapply(strsplit(junctions$hpo_id,"\\|"), function(x) ifelse(length(intersect(x, HPO))>0, "yes","no"))
    junctions$relevant_HPO[is.na(junctions$hpo_id) | junctions$hpo_id == "NA"] <- NA
  }

  print("Finished annotating!")
  # reorder columns
  junctions <- junctions %>%
    mutate(width = abs(start - end)) %>%
    unite("all_disease",hpo_term, mim_disease_desc, Orphanet, sep=";", remove = F) %>%
    dplyr::select("hgnc_symbols", "aliases","hgmd_tx", "all_disease","relevant_HPO","num_exons", "exon_ids","hgmd_tx_exon_ids", 
      "gene_width", "total_juncs", "gene_name","chrom", "start", "end", "strand","width", "annotated", "sample_ratio", "sample_count","zscore", "gtex_frac_detected","mean_gtex_ratio", 
      "sd_gtex_ratio", "mean_gtex_count", "mean_tpm", "median_tpm","sd_tpm", 
      'hgnc_ids','hpo_term','hpo_id','Orphanet','mim_gene_acc','mim_gene_desc','mim_disease_acc',
      'mim_disease_desc','gene_constraint_scores','hgmd_tx_version') %>%
    arrange(all_disease, desc(zscore)) %>%
    mutate_if(is.numeric, round, 5)
  return(junctions)
}

filter_junctions<-function(junc_stats, outlier_gtex_frac=0.3, outlier_mean_ratio=1.25, novel_junc_ratio=0.02,
                           novel_frac_gtex=0.02, z_cutoff=2, missing_sample_ratio=0.05, missing_gtex_ratio=0.1, 
                           missing_frac_gtex=0.98, extreme_z_cutoff = 10){
  #junc_stats<-left_join(junc_stats, annotated_junctions, by=c("chrom", "start", "end", "strand"))
  #junc_stats$annotated[is.na(junc_stats$annotated)]<-F
  #TODO stringent annotated z scores >5-6
  # stringent junction fraction 10%
  # stingent delta junction when compared to gtex 30%
  junc_stats_used <- junc_stats %>%
    dplyr::select(chrom, start, end, strand, annotated, gtex_frac_detected, zscore, sample_ratio, mean_gtex_ratio) %>%
    unique()

  outlier<-junc_stats_used %>%
    filter(annotated) %>%
    filter(gtex_frac_detected > outlier_gtex_frac) %>%
    filter(zscore >= z_cutoff) %>%
    filter(sample_ratio/mean_gtex_ratio > outlier_mean_ratio | sample_ratio/mean_gtex_ratio < (1/outlier_mean_ratio)) %>%
    mutate(class="outlier")

  novel<-junc_stats_used %>%
    filter(!annotated) %>%
    filter(gtex_frac_detected < novel_frac_gtex) %>%
    filter(sample_ratio > novel_junc_ratio) %>%
    filter( zscore >= z_cutoff | is.na(zscore)) %>%
  #  filter(mean_gtex_ratio != sample_ratio) %>% # filter out cases when gtex == sample == 1
    mutate(class="novel")

  #I'm hard coding these more might get hardcoded later on
  missing<-junc_stats_used %>%
    filter(sample_ratio < missing_sample_ratio | is.na(sample_ratio)) %>%
    filter(mean_gtex_ratio > missing_gtex_ratio) %>%
    filter(gtex_frac_detected > missing_frac_gtex) %>%
    mutate(class="missing")
  
  # extreme outlier: capturing any other junctions with an extremely high z score
  outlier_juncs <- rbind(outlier[,c("chrom", "start","end","strand")],
                         novel[,c("chrom", "start","end","strand")],
                         missing[,c("chrom", "start","end","strand")]) %>% 
    unite(junc_id, chrom, start, end, strand) %>% 
    pull(junc_id)

  extreme_outlier <- junc_stats_used %>% 
    unite(junc_id, chrom, start, end, strand, remove = F) %>% 
    filter(zscore >= extreme_z_cutoff & zscore != Inf & !junc_id %in% outlier_juncs) %>% 
    mutate(class = "extreme_outlier") %>%
    dplyr::select(-junc_id)

  dfs<-list(outlier, novel, missing, extreme_outlier)
  dfs <- dfs[sapply(dfs, length) > 0]
  results<-do.call("rbind", dfs) %>%
    unique() %>%
    group_by_at(vars(-class)) %>%
    summarise(class = paste(sort(unique(class)), collapse=";"),.groups = "drop") # merge class that are missing + outlier
   results <- left_join(junc_stats, results) %>%
     filter(!is.na(class))
  ##################
  # Summarize numbers/classes of junctions per gene
  #suppressMessages(junc_stats <- left_join(junc_stats, results))
  junctions_gene <- results[!is.na(results$hgnc_symbols), c("hgnc_symbols", "class")] %>%
    group_by(hgnc_symbols, class) %>%
    summarise(num = length(class),.groups = "drop") %>%
    spread(class, num)
  junctions_gene[is.na(junctions_gene)] <- 0
  junctions_gene <- junctions_gene %>%
    rowwise() %>%
    mutate(total_report = sum(missing, novel, outlier))

  # reorder columns
  results <- left_join(results, junctions_gene) %>%
    rowwise() %>%
    mutate(ifdisease = ifelse(hpo_term != "NA"|!is.na(Orphanet)|mim_disease_desc != "", "yes","no")) %>%
    dplyr::select("hgnc_symbols", "aliases","hgmd_tx", "all_disease","relevant_HPO", "num_exons","exon_ids","hgmd_tx_exon_ids",
      "gene_width", "total_juncs","total_report", "missing", "outlier", "novel", "gene_name","chrom", "start", "end", 
      "strand","width", "annotated","class", "sample_ratio", "sample_count","zscore", "gtex_frac_detected","mean_gtex_ratio", 
      "sd_gtex_ratio", "mean_gtex_count", "mean_tpm", "median_tpm","sd_tpm",everything()) %>%
    arrange(desc(relevant_HPO), desc(ifdisease), desc(abs(zscore))) %>%
        dplyr::select(-ifdisease) %>% dplyr::select(-any_of(c('missing;outlier')))

  #rownames(results)<-c(1:nrow(results))
  return(results)
}

option_list = list(
  make_option(c("-d", "--gtex_database"), type="character", help="gtex database path"),
  make_option(c("-t", '--tissue'), type="character",  help="tissue to analyze"),
  make_option(c("-o", '--output'), type="character", help="output name"),
  make_option(c("-g", '--txdb'), type="character", help="path to txDb"),
  make_option(c("-a", '--annotations'), type="character", help="annotated junctions file, part of star index"),
  make_option(c("-s", '--samplename'), type="character", help="samplename"),
  make_option(c("-j", "--junctions"), type="character", help="SJ.out.tab file for the sample"),
  make_option(c("-c", "--junction_cutoff"), type="integer", default=5, help="minimum uniq reads per junction, default=5"),
  make_option("--outlier_frac", default = 0.3, help = "Fraction of GTEx samples with the outlier junction (0.3)"),
  make_option("--outlier_mean_ratio", default = 1.25, help = "the score ratio for outlier boundaries (1.25 and 1/1.25)"),
  make_option("--novel_ratio", default = 0.02, help = "min score for novel junctions"),
  make_option("--novel_frac", default = 0.02, help = "max gtex frac for novel junctions"),
  make_option("--z_cutoff",type="numeric", default = 6, help = "zscore cutoff for novel and outlier junctions 3"),
  make_option("--missing_ratio", default = 0.05, help = "max score for missing junctions"),
  make_option("--missing_gtex_ratio", default = 0.1, help = "min gtex score for missing junctions"),
  make_option("--missing_gtex_frac", default = 0.98, help = "gtex fraction for missing junctions"),
  make_option("--ncores", type = "numeric", default = 4, help = "number of cores to use"),
  make_option("--controls", default = "gtex", help = "which samples to use as controls. Default is GTEx samples of the same tissue. Can provide a text file with processed sample names"),
  make_option("--control_directory", default = NULL, help = "where to find output samples. If not given, then used the root_dir of 'output'"),
  make_option("--refseq_txdb",default = "/hpf/projects/dig2/acelik_analysis/clinical_pipeline/genomes/hg37/hg19.refGene.UCSC.txdb", help = "Txdb of UCSC refseq genes"),
  make_option("--HPO_terms", default = "none", help = "Text file containing relevant HPO terms")
)


opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

# verify required files
if(!file.exists(opt$annotations)){
  stop("Gene annotation file is not found")
}

if(is.null(opt$junctions)){
  # if junction file is not found, try to find it from output dir
  output_dir <- dirname(opt$output)
  opt$junctions <- list.files(path = output_dir,
                          pattern = glob2rx("*SJ.out.tab"),
                          full.names = T)
  if(length(opt$junctions) !=  1){
    stop("No or more than 1 junction file (SJ.out) found")
  }
}

# if not using GTEx as controls, specify directory to find control SJ files
if(opt$controls != "gtex"){
  if(is.null(opt$control_directory)){
    control_dir <- basename(basename(opt$output))
  } else{
    control_dir <- opt$control_directory
  }
}


# read in HPO file
HPO_file <- opt$HPO_terms
HPO_terms <- NULL
if(!HPO_file == "none"){
  if(file.exists(HPO_file)){
    HPO_terms <- scan(HPO_file, what = "character")
    print(paste("Using HPO annotations in", HPO_file))
  }else{
    print(paste("HPO term file ", HPO_file, " does not exist"))
  }
}

# connect to getx database and txdb
db<-dbConnect(drv=SQLite(), opt$gtex_database)
txdb<-loadDb(opt$txdb)
genes_gr<-genes(txdb)
chroms<-unique(as.character(seqnames(genes_gr)@values))

#exons <- exonsBy(txdb, by = "gene")
# load refseq db
ref_txdb <- loadDb(opt$refseq_txdb)


message(Sys.time())
message("Getting control junctions")

if(opt$controls == "gtex"){
  control_junctions<-get_control_junctions(db, opt$tissue, opt$junction_cutoff)
} else{
  # get control junctions (using other processed samples)
  control_samples <- scan(opt$controls, what = "character")


  # if sample is within controls, remove sample from controls
  if(opt$samplename %in% control_samples){
    control_samples <- control_samples[control_samples != opt$samplename]
  }
  print(paste("Using", length(control_samples), "control samples"))

  control_junc_files <- sapply(control_samples, find_junction_file)
  control_junc_list <- lapply(control_junc_files, function(x) get_junction(SJ_file = x, chrs = chroms, opt$junction_cutoff))
  names(control_junc_list) <- control_samples

  control_junctions <- bind_rows(control_junc_list, .id = "sample_id") %>%
    spread("sample_id", "uniq_map", fill = 0)

}


message(Sys.time())
message("Getting sample_junctions")
sample_junctions <- get_junction(opt$junctions, chrs=chroms, opt$junction_cutoff)


message(Sys.time())
message("Calculating z scores")

#juncs <- get_zscores(control_junctions, sample_junctions, opt$samplename, ncores = opt$ncores)

control_junctions$chr <- as.factor(sapply(strsplit(control_junctions$junc_id,"_"),"[[",1))
control_junctions_list <- lapply(split(control_junctions, control_junctions$chr), function(x) x[colnames(x) != "chr"])

sample_junctions$chr <- sapply(strsplit(sample_junctions$junc_id,"_"),"[[",1)
sample_junctions_list <- lapply(split(sample_junctions, sample_junctions$chr), function(x) x[colnames(x) != "chr"])

temp <- mclapply(intersect(names(sample_junctions_list), names(control_junctions_list)), function(x) get_zscores(control_junctions_list[[x]], sample_junctions_list[[x]], opt$samplename, perc_zscore_cutoff = -1.2), mc.cores = opt$ncores)

# z scores based on junction ratios
juncs <- do.call(rbind, lapply(temp, "[[","ratio"))

# z scores based on junction reads percent
juncs_single <- do.call(rbind, lapply(temp, "[[","percent"))
juncs_single <- juncs_single[!is.na(juncs_single$chrom),] # remove entries that are NA

rm(temp)

# read in annotated junctions
star_junctions<-read.table(opt$annotations, header = F, sep = "\t", stringsAsFactors = F)
star_junctions<-star_junctions[,c(1:4)]#only take the first 4 columns
colnames(star_junctions) <- c("chrom", "start", "end", "strand")
star_junctions$annotated <- T


message(Sys.time())
message("Annotating")
annotated_junctions <- annotate_results(juncs, genes_gr, star_junctions, db, opt$tissue, ncores = opt$ncores, HPO = HPO_terms)

#modify columns in junction output to meet the requirement for clinical use
#replace | seperator to ;
#remove the .0 at the end of each accession number
#temp_mimacc <- gsub('\\|',';',annotated_junctions$mim_gene_acc)
#annotated_junctions$mim_gene_acc <- gsub("\\.[0-9]","",temp_mimacc)

#annotated_junctions$entrez_id <- gsub('\\|',';',annotated_junctions$entrez_id)
#annotated_junctions$hpo_id <- gsub('\\|',';',annotated_junctions$hpo_id)
#annotated_junctions$aliases <- gsub('\\|',';',annotated_junctions$aliases)
#annotated_junctions$mim_gene_desc <- gsub('\\|',';',annotated_junctions$mim_gene_desc)

# output all junctions
write.table(annotated_junctions, gzfile(gsub(".tsv", ".all.tsv.gz", opt$output)), col.names = T, row.names = F,quote = F, sep = "\t")
#write.table(annotated_junctions, gsub(".tsv", ".all.tsv", opt$output), col.names = T, row.names = F,quote = F, sep = "\t")


# annotate the junctions identified with percent reads analysis and output
annotated_junctions_single <- annotate_results(juncs_single, genes_gr, star_junctions, db, opt$tissue, 
  ncores = opt$ncores, anno_to_gene = F, HPO = HPO_terms)



#write.table(annotated_junctions_single, gsub(".tsv", ".singleJuncs.tsv", opt$output), col.names = T, row.names = F,
#            quote = F, sep = "\t")


message(Sys.time())
message("Filtering")
sig_junctions<-filter_junctions(junc_stats = annotated_junctions, 
  outlier_gtex_frac=opt$outlier_frac, outlier_mean_ratio=opt$outlier_mean_ratio,
                                novel_junc_ratio=opt$novel_ratio, novel_frac_gtex=opt$novel_frac,
                                z_cutoff=opt$z_cutoff, missing_sample_ratio=opt$missing_ratio,
                                missing_gtex_ratio=opt$missing_gtex_ratio, missing_frac_gtex=opt$missing_gtex_frac)

write.table(sig_junctions, opt$output, col.names = T, row.names = F,
            quote = F, sep = "\t")

#save.image(gsub("tsv","RData",opt$output))

message("Done")