# Focal amplification of the breast cancer genome
# -----------------------------------------------

library(tidyverse)
library(gUtils)
library(cowplot)

# Calculate F-score (focal amplification score)
gr_tiles <- gr.tile(Seqinfo(genome="hg38"), 5e6) %>% as_tibble() %>% filter(seqnames %in% str_c("chr", c(1:22, "X", "Y"))) %>% makeGRangesFromDataFrame()
hg38_seqinfo <- Seqinfo(genome="hg38")

df_centromere <- read_tsv("hg38_centromere_edit.tsv")
gr_centromere <- GRanges(df_centromere$chromosome, IRanges(start=df_centromere$exclude_start, end=df_centromere$exclude_end))

df_focal <- tibble()
for(f in list.files("/home/users/chrono0707/analysis/08_breast/06_seqz", full.names = T, pattern = "*segments.clean.txt$")) {
  print(f)
  df_temp <- read_tsv(f, show_col_types = F) %>% filter((end_pos-start_pos)<5e6) %>% 
    dplyr::slice(-queryHits(findOverlaps(GRanges(seqnames = `#CHROM`, IRanges(start=start_pos, end=end_pos)), resize(gr_centromere, width = width(gr_centromere)+(1e6*2), fix = "center"), type="within"))) # centromere근처 제거
  if(nrow(df_temp)==0) next
  
  gr_temp <- df_temp %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand=TRUE, seqnames.field=c('#CHROM'), start.field=c("start_pos"), end.field=c("end_pos"), seqinfo = hg38_seqinfo)
  
  df_focal <- bind_rows(
    df_focal, 
    gr_tiles %>% as_tibble() %>% transpose() %>% map(function(x) {
      df <- df_temp[queryHits(findOverlaps(gr_temp, GRanges(seqnames=x[1], ranges = IRanges(start=as.numeric(x[2]), as.numeric(x[3]))), ignore.strand=T)),] %>% mutate(width=end_pos-start_pos)
      if(nrow(df)>0) {
        return(tibble(segment=str_c(x[1], ":", x[2], "-", x[3]), narrowest_width=min(df$width), narrowest_majCN=df$majCN[which.min(df$width)], narrowest_minCN=df$minCN[which.min(df$width)], highest_width=df$width[which.max(df$majCN)], highest_majCN=max(df$majCN), highest_minCN=df$minCN[which.max(df$majCN)], n_segments=nrow(df)))
      } 
    }) %>% do.call(bind_rows, .) %>% mutate(tumor_id=str_sub(basename(f), end=19)) %>% dplyr::select(tumor_id, everything())
  )
}

# Load ecDNA results from ampliconarchitect
df_aa <- tibble()
for(f in list.files("/home/users/chrono0707/analysis/08_breast/14_ampliconarchitect", pattern="*.amplicon_classification_profiles.tsv", full.names = T, recursive = T)) {
  df_aa <- bind_rows(
    df_aa,
    read_tsv(f, show_col_types = F) %>% mutate(tumor_id=str_sub(basename(f), end=19)) %>% dplyr::select(tumor_id, everything()) %>% 
      left_join(readLines(str_replace(f, "_amplicon_classification_profiles.tsv", "_summary.txt")) %>% str_extract("^.*OncogenesAmplified.*$") %>% as_tibble() %>% filter(!is.na(value)) %>% mutate(amplicon_number=str_extract(value, 'amplicon\\d+'), genes=str_extract(value, "(?<=OncogenesAmplified = ).*$")) %>% dplyr::select(amplicon_number, genes)) %>% 
      left_join(readLines(str_replace(f, "_amplicon_classification_profiles.tsv", "_summary.txt")) %>% str_extract("^.*TotalIntervalSize.*$")   %>% as_tibble() %>% filter(!is.na(value)) %>% mutate(amplicon_number=str_extract(value, 'amplicon\\d+'), total_interval_size=str_extract(value, "(?<=TotalIntervalSize = ).*$") %>% as.numeric())  %>% dplyr::select(amplicon_number, total_interval_size)) %>% 
      left_join(readLines(str_replace(f, "_amplicon_classification_profiles.tsv", "_summary.txt")) %>% str_extract("^.*AmplifiedIntervalSize.*$")   %>% as_tibble() %>% filter(!is.na(value)) %>% mutate(amplicon_number=str_extract(value, 'amplicon\\d+'), amplified_interval_size=str_extract(value, "(?<=AmplifiedIntervalSize = ).*$") %>% as.numeric())  %>% dplyr::select(amplicon_number, amplified_interval_size)) %>% 
      left_join(readLines(str_replace(f, "_amplicon_classification_profiles.tsv", "_summary.txt")) %>% str_extract("^.*AverageAmplifiedCopyCount.*$")   %>% as_tibble() %>% filter(!is.na(value)) %>% mutate(amplicon_number=str_extract(value, 'amplicon\\d+'), average_amplified_copy_count=str_extract(value, "(?<=AverageAmplifiedCopyCount = ).*$") %>% as.numeric())  %>% dplyr::select(amplicon_number, average_amplified_copy_count)) %>% 
      left_join(readLines(str_replace(f, "_amplicon_classification_profiles.tsv", "_summary.txt")) %>% str_extract("^.*#Intervals.*$")   %>% as_tibble() %>% filter(!is.na(value)) %>% mutate(amplicon_number=str_extract(value, 'amplicon\\d+'), n_intervals=str_extract(value, "(?<=Intervals = ).*$") %>% as.numeric())  %>% dplyr::select(amplicon_number, n_intervals)) %>% 
      left_join(readLines(str_replace(f, "_amplicon_classification_profiles.tsv", "_summary.txt")) %>% str_extract("^.*#Chromosomes.*$")   %>% as_tibble() %>% filter(!is.na(value)) %>% mutate(amplicon_number=str_extract(value, 'amplicon\\d+'), n_chroms=str_extract(value, "(?<=Chromosomes = ).*$") %>% as.numeric())  %>% dplyr::select(amplicon_number, n_chroms)) %>% 
      left_join(readLines(str_replace(f, "_amplicon_classification_profiles.tsv", "_summary.txt")) %>% str_extract("^.*#CoverageShifts =.*$")   %>% as_tibble() %>% filter(!is.na(value)) %>% mutate(amplicon_number=str_extract(value, 'amplicon\\d+'), n_cov_shifts=str_extract(value, "(?<=CoverageShifts = ).*$") %>% as.numeric())  %>% dplyr::select(amplicon_number, n_cov_shifts)) %>% 
      left_join(readLines(str_replace(f, "_amplicon_classification_profiles.tsv", "_summary.txt")) %>% str_extract("^.*#Foldbacks.*$")   %>% as_tibble() %>% filter(!is.na(value)) %>% mutate(amplicon_number=str_extract(value, 'amplicon\\d+'), n_foldbacks=str_extract(value, "(?<=Foldbacks = ).*$") %>% as.numeric())  %>% dplyr::select(amplicon_number, n_foldbacks)) %>% 
      left_join(read_tsv(str_replace(f, "_amplicon_classification_profiles.tsv", "_feature_basic_properties.tsv"), show_col_types = F) %>% mutate(amplicon_number=str_extract(feature_ID, "amplicon\\d+")) %>% group_by(amplicon_number) %>% summarise(maxCN=ifelse(nrow(.), max(max_feature_CN), numeric()))) 
  )
}
df_aa <- df_aa %>% dplyr::rename(ecdna=`ecDNA+`, bfb=`BFB+`) %>% mutate(genes=str_remove(genes, ",$"))

df_ecdna <- tibble()
for(l in df_aa %>% filter(ecdna=="Positive") %>% dplyr::select(tumor_id, amplicon_number) %>% transpose() %>% as.list()) {
  df_ecdna <- bind_rows(
    df_ecdna, 
    read_tsv(str_c("../14_ampliconarchitect/", l[1], "_", str_sub(l[1], end=14), "-10AD/", l[1], "_", str_sub(l[1], end=14), "-10AD_", l[2], "_edges_cnseg.txt"), skip=1, col_names=c("bp_edge", "number_of_read_pairs", "homology_size", "homology_sequence")) %>%
      mutate(
        chrom1=str_extract(bp_edge, "^chr[^:]+(?=:)"),
        pos1  =str_extract(bp_edge, "(?<=:)\\d+(?=[^$])") %>% as.numeric(),
        chrom2=str_extract(bp_edge, "(?<=>)chr[^:]+(?=:)"),
        pos2  =str_extract(bp_edge, "(?<=:)\\d+(?=[-+]$)") %>% as.numeric(),
        tumor_id=l[[1]],
        amplicon_number=l[[2]]
      ) %>% dplyr::select(tumor_id, amplicon_number, chrom1, pos1, chrom2, pos2)
  ) 
}

# Draw figures 
# ------------

# Figure 4A. The distribution of F-score and segments involved in extrachromosomal DNA formation across the genome
df_temp1 <- df_focal %>% 
  filter(tumor_id %in% df_crf$tumor_id) %>% 
  filter(tumor_id %in% (df_decision %>% filter(cellularity>0.20) %>% pull(tumor_id))) %>% 
  filter(tumor_id %in% (df_tmb %>% filter(tmb>1000) %>% pull(tumor_id))) %>% 
  filter(n_segments>1) %>% # If a tile is not divided into multiple segments, it is not a focal amplification smaller than 5e6 in size
  mutate(chr=str_extract(segment, "chr.+(?=:)"), start=as.numeric(str_extract(segment, "(?<=:)\\d+(?=-)")), end=as.numeric(str_extract(segment, "(?<=-)\\d+"))) %>% 
  dplyr::slice(-queryHits(findOverlaps(GRanges(seqnames = chr, IRanges(start=start, end=end)), resize(gr_centromere, width = width(gr_centromere)+(1e6*2), fix = "center")))) %>%  # Remove regions near the centromere
  mutate(ratio=narrowest_majCN*1e6/narrowest_width) %>% 
  group_by(chr, start, end) %>% summarise(ratio=sum(ratio)) %>% 
  filter(chr %in% str_c("chr", c(1:22, "X"))) %>% 
  mutate(segment=str_c(chr, ":", start, "-", as.integer(end))) %>% 
  mutate(chr = factor(chr, levels=str_c("chr", c(1:22, "X")))) %>% 
  arrange(chr, start) %>% mutate(segment=factor(segment, levels=segment)) %>% 
  mutate(bin=ifelse(chr %in% str_c("chr", c(seq(1,22,2), "X")), TRUE, FALSE))

g1 <- df_temp1 %>% 
  ggplot(aes(x=segment, y=ratio)) + 
  geom_bar(aes(fill=bin, y=500), stat='identity', width=1, color=NA) + 
  geom_bar(stat='identity', width=1, color=NA) + 
  scale_fill_manual(values=c("TRUE"="grey90", "FALSE"="white")) + 
  theme_minimal() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid=element_blank(), legend.position='none', axis.title.x=element_blank())

df_temp2 <- df_temp1 %>% as_tibble() %>% transpose() %>% map(function(x) {
  n1 <- df_ecdna %>% filter(chrom1==x[1] & pos1 >= as.numeric(x[2]) & pos1 < as.numeric(x[3])) %>% dplyr::select(tumor_id) %>% distinct() %>% nrow()
  n2 <- df_ecdna %>% filter(chrom2==x[1] & pos2 >= as.numeric(x[2]) & pos2 < as.numeric(x[3])) %>% dplyr::select(tumor_id) %>% distinct() %>% nrow()
  return(tibble(chrom=x[1], start=as.numeric(x[2]), end=as.numeric(x[3]), n=n1+n2))
}) %>% do.call(bind_rows, .) %>% 
  mutate(segment=str_c(chrom, ":", start, "-", as.integer(end))) %>% 
  mutate(chrom = factor(chrom, levels=str_c("chr", c(1:22, "X")))) %>% 
  arrange(chrom, start) %>% mutate(segment=factor(segment, levels=segment)) %>% 
  mutate(bin=ifelse(chrom %in% str_c("chr", c(seq(1,22,2), "X")), TRUE, FALSE))
 
g2 <- df_temp2 %>% ggplot(aes(x=segment, y=n, group=1)) + 
  geom_bar(aes(fill=bin, y=max(n)), stat='identity', width=1, color=NA) + 
  geom_area(fill=NA, color='black') +
  scale_fill_manual(values=c("TRUE"="grey90", "FALSE"="white")) + 
  theme_minimal() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid=element_blank(), legend.position='none', axis.title.x=element_blank())

plot_grid(g1, g2, align="v", rel_heights = c(0.6, 0.4), ncol=1)


# Figure 4B. Mechanisms involved in ERBB2, FGFR1, and CCND1 focal amplification.

# ERBB2
gr_her2   <- GRanges(seqnames = "chr17", ranges = IRanges(start=39700080,     end=39728662)) # ERBB2 location
df_her2 <- tibble()
for(f in list.files("/home/users/chrono0707/analysis/08_breast/06_seqz", full.names = T, pattern = "*segments.clean.txt$")) {
  
  tn_pair <- basename(f) %>% str_remove(".segments.clean.txt")  
  if(!file.exists(str_c("/home/users/chrono0707/analysis/08_breast/06_seqz/", basename(f) %>% str_remove("segments.clean.txt"), "seqz.solution.tsv"))) next
  
  df_seqz <- read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/06_seqz/", basename(f) %>% str_remove("segments.clean.txt"), "seqz.solution.tsv"), show_col_types = FALSE)
  gr_seqz <- GRanges(seqnames = df_seqz$chromosome, ranges = IRanges(start=df_seqz$start.pos, end=df_seqz$end.pos))
  
  idx <- queryHits(findOverlaps(gr_seqz, gr_her2))
  if(length(idx)) {
    df_seqz_her2 <- df_seqz[idx,]  # Keep only the regions overlapping with ERBB2
  } else { # If there is no overlap, take the side that is closer from both directions
    idx <- df_seqz %>% mutate(idx=1:nrow(.)) %>% filter(chromosome=="chr17") %>%  mutate(distance=ifelse(end.pos<start(gr_her2), start(gr_her2)-end.pos, ifelse(start.pos>end(gr_her2), start.pos-end(gr_her2), NA))) %>% arrange(distance) %>% pull(idx) %>% .[1]
    df_seqz_her2 <- df_seqz[idx,]
  }

  df_her2 <- df_her2 %>% bind_rows(tibble(tumor_id         = str_sub(tn_pair, end=19), 
                                            start           = min(df_seqz_her2$start.pos),
                                            end             = max(df_seqz_her2$end.pos),
                                            GCN             = read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/06_seqz/", basename(f) %>% str_remove("segments.clean.txt"), "cancer_gene_cn.tsv"), show_col_types = FALSE) %>% filter(gene_name=="ERBB2") %>% .$GeneCN,
                                            majCN           = max(df_seqz_her2$A),
                                            minCN           = min(df_seqz_her2$B),
                                            mean_CN_chr17   = mean(df_seqz[-idx,] %>% filter(chromosome=='chr17') %>% .$CNt),  # The average CN of regions excluding the ERBB2 region
                                            mean_CN_padding_5mbp = if(nrow(df_seqz[queryHits(findOverlaps(gr_seqz, GRanges(seqnames = "chr17", ranges = IRanges(start=39700080-5e6, end=39728662+5e6)))),])==1) df_seqz[queryHits(findOverlaps(gr_seqz, GRanges(seqnames = "chr17", ranges = IRanges(start=39700080-5e6, end=39728662+5e6)))),]$CNt[1] else mean(df_seqz[-idx,][queryHits(findOverlaps(gr_seqz[-idx], GRanges(seqnames = "chr17", ranges = IRanges(start=39700080-5e6, end=39728662+5e6)))),]$CNt),
                                            mean_CN_padding_1mbp = if(nrow(df_seqz[queryHits(findOverlaps(gr_seqz, GRanges(seqnames = "chr17", ranges = IRanges(start=39700080-1e6, end=39728662+1e6)))),])==1) df_seqz[queryHits(findOverlaps(gr_seqz, GRanges(seqnames = "chr17", ranges = IRanges(start=39700080-1e6, end=39728662+1e6)))),]$CNt[1] else mean(df_seqz[-idx,][queryHits(findOverlaps(gr_seqz[-idx], GRanges(seqnames = "chr17", ranges = IRanges(start=39700080-1e6, end=39728662+1e6)))),]$CNt)
  ))
}

df_her2 <- df_her2 %>% mutate(width=end-start, delta_CN=GCN-mean_CN_chr17, delta_CN_padding_5mbp=GCN-mean_CN_padding_5mbp, delta_CN_padding_1mbp=GCN-mean_CN_padding_1mbp)
df_her2 <- filter(tumor_id %in% df_crf$tumor_id)

df_her2 %>% left_join(df_aa %>% mutate(her2=str_detect(genes, "ERBB2")) %>% filter(her2) %>% dplyr::select(tumor_id, ecdna, bfb, her2, amplicon_decomposition_class)) %>%
  mutate(mechanism=ifelse(ecdna=="Positive", "ecDNA", ifelse(bfb=="Positive", "BFB", amplicon_decomposition_class))) %>% 
  mutate(ifelse(mechanism=="No amp/Invalid", NA, mechanism)) %>% 
  dplyr::count(focal, mechanism, her2) %>% 
  filter(focal) %>% 
  mutate(p=n/sum(n)) %>% 
  mutate(mechanism=factor(mechanism, levels=c("BFB", "ecDNA", "Complex non-cyclic", "Linear amplification"))) %>% 
  ggplot(aes(x=factor(1), y=p, fill=mechanism)) + geom_col(position='fill') + geom_col(aes(x=0, y=0)) +  coord_polar(theta="y") + theme_void() + #scale_fill_(na.value="gray90")
  scale_fill_manual(values=c("BFB"=unikn::uni_mannheim_2[[4]], "ecDNA"=unikn::uni_mannheim_2[[6]], "Complex non-cyclic"=unikn::uni_mannheim_2[[5]], "Linear amplification"=unikn::uni_mannheim_2[[7]]))


# FGFR1
gr_fgfr1 <- GRanges(seqnames = "chr8", ranges = IRanges(start=38411143, end=38468635)) # FGFR1 location
df_fgfr1 <- tibble()
for(f in list.files("/home/users/chrono0707/analysis/08_breast/06_seqz", full.names = T, pattern = "*segments.clean.txt$")) {
  
  tn_pair <- basename(f) %>% str_remove(".segments.clean.txt")
  
  if(!file.exists(str_c("/home/users/chrono0707/analysis/08_breast/06_seqz/", basename(f) %>% str_remove("segments.clean.txt"), "seqz.solution.tsv"))) next
  
  df_seqz <- read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/06_seqz/", basename(f) %>% str_remove("segments.clean.txt"), "seqz.solution.tsv"), show_col_types = FALSE)
  gr_seqz <- GRanges(seqnames = df_seqz$chromosome, ranges = IRanges(start=df_seqz$start.pos, end=df_seqz$end.pos))
  
  idx <- queryHits(findOverlaps(gr_seqz, gr_fgfr1))
  if(length(idx)) {
    df_seqz_fgfr1 <- df_seqz[idx,]  
  } else { 
    idx <- df_seqz %>% mutate(idx=1:nrow(.)) %>% filter(chromosome=="chr8") %>%  mutate(distance=ifelse(end.pos<start(gr_fgfr1), start(gr_fgfr1)-end.pos, ifelse(start.pos>end(gr_fgfr1), start.pos-end(gr_fgfr1), NA))) %>% arrange(distance) %>% pull(idx) %>% .[1]
    df_seqz_fgfr1 <- df_seqz[idx,]
  }
  
  # Check if it has been called as a driver
  tryCatch({
    df_cnv_dv <- read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/03_somatic/", tn_pair, ".cnv_dv.vcf"), comment="##")
    df_cnv_dv <- df_cnv_dv %>% mutate(symbol=str_extract(INFO, "(?<=SYMBOL=)[A-Z0-9]*(?<!;)")) 
    df_cnv_dv <- df_cnv_dv %>% mutate(GCN=str_extract(INFO, "(?<=GCN=)[-.A-Z0-9]*(?<!;)"))
    df_cnv_dv <- df_cnv_dv %>% filter(symbol == "FGFR1")
    is_driver <- if(nrow(df_cnv_dv)) TRUE else FALSE
  },
  error = function(e) {is_driver <- FALSE})  
  
  df_fgfr1 <- df_fgfr1 %>% bind_rows(tibble(tumor_id         = str_sub(tn_pair, end=19), 
                                            start           = min(df_seqz_fgfr1$start.pos),
                                            end             = max(df_seqz_fgfr1$end.pos),
                                            GCN             = read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/06_seqz/", basename(f) %>% str_remove("segments.clean.txt"), "cancer_gene_cn.tsv"), show_col_types = FALSE) %>% filter(gene_name=="FGFR1") %>% .$GeneCN,
                                            majCN           = max(df_seqz_fgfr1$A),
                                            minCN           = min(df_seqz_fgfr1$B),
                                            driver          = is_driver,
                                            mean_CN_chr8    = mean(df_seqz[-idx,] %>% filter(chromosome=='chr8') %>% .$CNt),  
                                            mean_CN_padding_5mbp = if(nrow(df_seqz[queryHits(findOverlaps(gr_seqz, GRanges(seqnames = "chr8", ranges = IRanges(start=38411143-5e6, end=38468635+5e6)))),])==1) df_seqz[queryHits(findOverlaps(gr_seqz, GRanges(seqnames = "chr8", ranges = IRanges(start=38411143-5e6, end=38468635+5e6)))),]$CNt[1] else mean(df_seqz[-idx,][queryHits(findOverlaps(gr_seqz[-idx], GRanges(seqnames = "chr8", ranges = IRanges(start=38411143-5e6, end=38468635+5e6)))),]$CNt),
                                            mean_CN_padding_1mbp = if(nrow(df_seqz[queryHits(findOverlaps(gr_seqz, GRanges(seqnames = "chr8", ranges = IRanges(start=38411143-1e6, end=38468635+1e6)))),])==1) df_seqz[queryHits(findOverlaps(gr_seqz, GRanges(seqnames = "chr8", ranges = IRanges(start=38411143-1e6, end=38468635+1e6)))),]$CNt[1] else mean(df_seqz[-idx,][queryHits(findOverlaps(gr_seqz[-idx], GRanges(seqnames = "chr8", ranges = IRanges(start=38411143-1e6, end=38468635+1e6)))),]$CNt)
  ))
}

df_fgfr1 <- df_fgfr1 %>% mutate(width=end-start, delta_CN=GCN-mean_CN_chr8, delta_CN_padding_5mbp=GCN-mean_CN_padding_5mbp, delta_CN_padding_1mbp=GCN-mean_CN_padding_1mbp)
df_fgfr1 <- df_fgfr1 %>% filter(tumor_id %in% df_crf$tumor_id)

df_fgfr1 %>% left_join(df_aa %>% mutate(fgfr1=str_detect(genes, "FGFR1")) %>% filter(fgfr1) %>% dplyr::select(tumor_id, ecdna, bfb, fgfr1, amplicon_decomposition_class)) %>%
  mutate(mechanism=ifelse(ecdna=="Positive", "ecDNA", ifelse(bfb=="Positive", "BFB", amplicon_decomposition_class))) %>% 
  mutate(ifelse(mechanism=="No amp/Invalid", NA, mechanism)) %>% 
  dplyr::count(driver, mechanism) %>% 
  filter(driver) %>% 
  mutate(p=n/sum(n)) %>% 
  mutate(mechanism=factor(mechanism, levels=c("BFB", "ecDNA", "Complex non-cyclic", "Linear amplification"))) %>% 
  ggplot(aes(x=factor(1), y=p, fill=mechanism)) + geom_col(position='fill') + geom_col(aes(x=0, y=0)) +  coord_polar(theta="y") + theme_void() + #scale_fill_(na.value="gray90")
  scale_fill_manual(values=c("BFB"=unikn::uni_mannheim_2[[4]], "ecDNA"=unikn::uni_mannheim_2[[6]], "Complex non-cyclic"=unikn::uni_mannheim_2[[5]], "Linear amplification"=unikn::uni_mannheim_2[[7]]))

# CCND1
gr_ccnd1 <- GRanges(seqnames = "chr11", ranges = IRanges(start=69641156, end=69654474)) # CCND1 location
df_ccnd1 <- tibble()
for(f in list.files("/home/users/chrono0707/analysis/08_breast/06_seqz", full.names = T, pattern = "*segments.clean.txt$")) {
  
  tn_pair <- basename(f) %>% str_remove(".segments.clean.txt")
  
  if(!file.exists(str_c("/home/users/chrono0707/analysis/08_breast/06_seqz/", basename(f) %>% str_remove("segments.clean.txt"), "seqz.solution.tsv"))) next
  
  df_seqz <- read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/06_seqz/", basename(f) %>% str_remove("segments.clean.txt"), "seqz.solution.tsv"), show_col_types = FALSE)
  gr_seqz <- GRanges(seqnames = df_seqz$chromosome, ranges = IRanges(start=df_seqz$start.pos, end=df_seqz$end.pos))
  
  idx <- queryHits(findOverlaps(gr_seqz, gr_ccnd1))
  if(length(idx)) {
    df_seqz_ccnd1 <- df_seqz[idx,]  
  } else {
    idx <- df_seqz %>% mutate(idx=1:nrow(.)) %>% filter(chromosome=="chr8") %>%  mutate(distance=ifelse(end.pos<start(gr_ccnd1), start(gr_ccnd1)-end.pos, ifelse(start.pos>end(gr_ccnd1), start.pos-end(gr_ccnd1), NA))) %>% arrange(distance) %>% pull(idx) %>% .[1]
    df_seqz_ccnd1 <- df_seqz[idx,]
  }
  
  # Check if it has been called as a driver
  tryCatch({
    df_cnv_dv <- read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/03_somatic/", tn_pair, ".cnv_dv.vcf"), comment="##")
    df_cnv_dv <- df_cnv_dv %>% mutate(symbol=str_extract(INFO, "(?<=SYMBOL=)[A-Z0-9]*(?<!;)")) 
    df_cnv_dv <- df_cnv_dv %>% mutate(GCN=str_extract(INFO, "(?<=GCN=)[-.A-Z0-9]*(?<!;)"))
    df_cnv_dv <- df_cnv_dv %>% filter(symbol == "CCND1")
    is_driver <- if(nrow(df_cnv_dv)) TRUE else FALSE
  },
  error = function(e) {is_driver <- FALSE})  
  
  df_ccnd1 <- df_ccnd1 %>% bind_rows(tibble(tumor_id         = str_sub(tn_pair, end=19), 
                                            start           = min(df_seqz_ccnd1$start.pos),
                                            end             = max(df_seqz_ccnd1$end.pos),
                                            GCN             = read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/06_seqz/", basename(f) %>% str_remove("segments.clean.txt"), "cancer_gene_cn.tsv"), show_col_types = FALSE) %>% filter(gene_name=="CCND1") %>% .$GeneCN,
                                            majCN           = max(df_seqz_ccnd1$A),
                                            minCN           = min(df_seqz_ccnd1$B),
                                            driver          = is_driver,
                                            mean_CN_chr11    = mean(df_seqz[-idx,] %>% filter(chromosome=='chr11') %>% .$CNt), 
                                            mean_CN_padding_5mbp = if(nrow(df_seqz[queryHits(findOverlaps(gr_seqz, GRanges(seqnames = "chr11", ranges = IRanges(start=69641156-5e6, end=69654474+5e6)))),])==1) df_seqz[queryHits(findOverlaps(gr_seqz, GRanges(seqnames = "chr11", ranges = IRanges(start=69641156-5e6, end=69654474+5e6)))),]$CNt[1] else mean(df_seqz[-idx,][queryHits(findOverlaps(gr_seqz[-idx], GRanges(seqnames = "chr11", ranges = IRanges(start=69641156-5e6, end=69654474+5e6)))),]$CNt),
                                            mean_CN_padding_1mbp = if(nrow(df_seqz[queryHits(findOverlaps(gr_seqz, GRanges(seqnames = "chr11", ranges = IRanges(start=69641156-1e6, end=69654474+1e6)))),])==1) df_seqz[queryHits(findOverlaps(gr_seqz, GRanges(seqnames = "chr11", ranges = IRanges(start=69641156-1e6, end=69654474+1e6)))),]$CNt[1] else mean(df_seqz[-idx,][queryHits(findOverlaps(gr_seqz[-idx], GRanges(seqnames = "chr11", ranges = IRanges(start=69641156-1e6, end=69654474+1e6)))),]$CNt)
  ))
}

df_ccnd1 <- df_ccnd1 %>% mutate(width=end-start, delta_CN=GCN-mean_CN_chr11, delta_CN_padding_5mbp=GCN-mean_CN_padding_5mbp, delta_CN_padding_1mbp=GCN-mean_CN_padding_1mbp)
df_ccnd1 <- df_ccnd1 %>% filter(tumor_id %in% df_crf$tumor_id)

df_ccnd1 %>% left_join(df_aa %>% mutate(ccnd1=str_detect(genes, "CCND1")) %>% filter(ccnd1) %>% dplyr::select(tumor_id, ecdna, bfb, ccnd1, amplicon_decomposition_class)) %>%
  mutate(mechanism=ifelse(ecdna=="Positive", "ecDNA", ifelse(bfb=="Positive", "BFB", amplicon_decomposition_class))) %>% 
  mutate(ifelse(mechanism=="No amp/Invalid", NA, mechanism)) %>% 
  dplyr::count(driver, mechanism) %>% 
  filter(driver) %>% 
  mutate(p=n/sum(n)) %>% 
  mutate(mechanism=factor(mechanism, levels=c("BFB", "ecDNA", "Complex non-cyclic", "Linear amplification"))) %>% 
  ggplot(aes(x=factor(1), y=p, fill=mechanism)) + geom_col(position='fill') + geom_col(aes(x=0, y=0)) +  coord_polar(theta="y") + theme_void() + #scale_fill_(na.value="gray90")
  scale_fill_manual(values=c("BFB"=unikn::uni_mannheim_2[[4]], "ecDNA"=unikn::uni_mannheim_2[[6]], "Complex non-cyclic"=unikn::uni_mannheim_2[[5]], "Linear amplification"=unikn::uni_mannheim_2[[7]]))

# Figure 3C. Mutual exclusivity of FGFR1- and CCND1-carrying ecDNA with ERBB2-carrying ecDNA
df_aa %>% filter(ecdna=="Positive") %>% group_by(tumor_id) %>% summarise(erbb2=any(str_detect(genes, "ERBB2")), fgfr1=any(str_detect(genes, "FGFR1"))) %>% dplyr::count(erbb2, fgfr1) %>% mutate(n=n/sum(n)) %>% ggplot(aes(x=erbb2, y=fgfr1, fill=n*100)) + geom_tile() + theme_minimal() + scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 30, l2 = 100, p1 = .9, p2 = 1.2) + theme(panel.background = element_rect(color="black"), panel.grid=element_blank()) 
df_aa %>% filter(ecdna=="Positive") %>% mutate(erbb2=str_detect(genes, "ERBB2"), fgfr1=str_detect(genes, "FGFR1")) %>% dplyr::select(tumor_id, erbb2, fgfr1) %>% distinct() %>% dplyr::count(erbb2, fgfr1) %>% pull(n) %>% matrix(nrow=2) %>% chisq.test()

df_aa %>% filter(ecdna=="Positive") %>% group_by(tumor_id) %>% summarise(erbb2=any(str_detect(genes, "ERBB2")), ccnd1=any(str_detect(genes, "CCND1"))) %>% dplyr::count(erbb2, ccnd1) %>% mutate(n=n/sum(n)) %>% ggplot(aes(x=erbb2, y=ccnd1, fill=n*100)) + geom_tile() + theme_minimal() + scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 30, l2 = 100, p1 = .9, p2 = 1.2) + theme(panel.background = element_rect(color="black"), panel.grid=element_blank())
df_aa %>% filter(ecdna=="Positive") %>% mutate(erbb2=str_detect(genes, "ERBB2"), ccnd1=str_detect(genes, "CCND1")) %>% dplyr::select(tumor_id, erbb2, ccnd1) %>% distinct() %>% dplyr::count(erbb2, ccnd1) %>% pull(n) %>% matrix(nrow=2) %>% chisq.test()
