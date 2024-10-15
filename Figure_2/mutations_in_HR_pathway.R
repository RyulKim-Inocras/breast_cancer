# Germline and somatic mutations in HRD positive samples 
# ------------------------------------------------------


# genes in homologous recombination pathway 
HR_genes <- c("ARID1A", "ATM", "ATRX", "BAP1", "BARD1", "BLM", "BRCA1", "BRCA2", "BRIP1", "CHEK1", "CHECK2", "FANCA", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCL", "MRE11A", "NBN", "PALB2", "RAD50", "RAD51", "RAD51B", "RAD51C", "WRN", "ABL1", "ATR", "CDKN1A", "CDKN2A", "CHEK2", "FANCI", "KDR", "MUTYH", "POLE", "POLQ", "RAD51D", "RAD54L", "TP53")


# Make table for somatic mututations in homologous recombination pathway
df_s <- tibble(SYMBOL=character()) 

for(f in list.files("~/analysis/08_breast/06_seqz/", pattern="*segments.clean.txt$", full.names = TRUE)) {
  
  tn_pair <- basename(f) %>% str_remove(".segments.clean.txt")

  normal_id <- str_split(basename(f) %>% str_remove(".segments.clean.txt"), "_")[[1]][2]
  tumor_id  <- str_split(basename(f) %>% str_remove(".segments.clean.txt"), "_")[[1]][1]

  if(!(tumor_id %in% df_crf$tumor_id)) next # if not in sample list, skip!
  
  df_merged <- tibble()

  df_snv <- read_tsv(str_c("~/analysis/08_breast/03_somatic/", tn_pair, ".snv_dv.vcf"), comment="##", show_col_types = FALSE) %>% rename(CHROM=`#CHROM`) %>% 
    mutate(csq     = str_match(INFO, "(?<=CSQ=)[^;]*(?=;)"),
           consequence = sapply(csq, function(x) str_split(x, "[|]", simplify = T)[2]),
           impact  = sapply(csq, function(x) str_split(x, "[|]", simplify = T)[3]),
           SYMBOL  = sapply(csq, function(x) str_split(x, "[|]", simplify = T)[4]),
           TYPE    = "PM",
           START   = POS,
           END     = POS)
  
  df_snv <- df_snv %>% 
    filter(FILTER=="PASS") %>% 
    filter(!(str_detect(consequence, "intron_variant") & !str_detect(consequence, "splice_region_variant"))) %>%
    filter(!(str_detect(consequence, "upstream_gene_variant"))) %>% 
    filter(SYMBOL %in% HR_genes)
  
  if(nrow(df_snv)) df_merged <- bind_rows(df_merged, df_snv %>% select(CHROM, START, END, SYMBOL, TYPE))
  
  df_ind <- read_tsv(str_c("~/analysis/08_breast/03_somatic/", tn_pair, ".ind_dv.vcf"), comment="##", show_col_types = FALSE) %>% rename(CHROM=`#CHROM`) %>% 
    mutate(csq     = str_match(INFO, "(?<=CSQ=)[^;]*(?=;)"),
           consequence = sapply(csq, function(x) str_split(x, "[|]", simplify = T)[2]),
           impact  = sapply(csq, function(x) str_split(x, "[|]", simplify = T)[3]),
           SYMBOL  = sapply(csq, function(x) str_split(x, "[|]", simplify = T)[4]),
           TYPE    = "PM",
           START   = POS,
           END     = POS)
  
  df_ind <- df_ind %>% 
    filter(FILTER=="PASS") %>% 
    filter(!(str_detect(consequence, "intron_variant") & !str_detect(consequence, "splice_region_variant"))) %>%
    filter(!(str_detect(consequence, "upstream_gene_variant"))) %>% 
    filter(SYMBOL %in% HR_genes)
  
  if(nrow(df_ind)) df_merged <- bind_rows(df_merged, df_ind %>% select(CHROM, START, END, SYMBOL, TYPE))
  
  df_fsn <- read_tsv(str_c("~/analysis/08_breast/03_somatic/", tn_pair, ".fsn_dv.vcf"), comment="##", show_col_types = FALSE) %>% rename(CHROM=`#CHROM`) %>% 
    mutate(csq1    = str_match(INFO, "(?<=CSQ1=)[^;]*(?=;)"),
           csq2    = str_match(INFO, "(?<=CSQ2=)[^;]*(?=;)"),
           symbol1 = sapply(csq1, function(x) str_split(x, "[|]", simplify = T)[2]),
           symbol2 = sapply(csq2, function(x) str_split(x, "[|]", simplify = T)[2]),
           TYPE    = str_match(INFO, "(?<=SVTYPE=)[^;]*(?=;)"),
           POS1    = as.numeric(POS),
           POS2    = as.numeric(ifelse(TYPE=="BND", as.numeric(str_match(INFO, "(?<=POS2=)[^;]*(?=;)")), as.numeric(str_match(INFO, "(?<=END=)[^;]*(?=;)")))))
  
  df_fsn <- df_fsn %>% filter(!(TYPE=="DUP" & symbol1 != symbol2)) %>% filter(FILTER=="PASS") # If it is a DUP, both breakpoints must be within the same gene to pass
  df_fsn <- bind_rows(
    df_fsn %>% select(CHROM, POS1, symbol1, TYPE) %>% rename(START=POS1, SYMBOL=symbol1),
    df_fsn %>% select(CHROM, POS2, symbol2, TYPE) %>% rename(START=POS2, SYMBOL=symbol2)
  )
  
  df_fsn <- df_fsn %>% mutate(END=START) %>% filter(SYMBOL %in% HR_genes)
  
  if(nrow(df_fsn)) df_merged <- bind_rows(df_merged, df_fsn %>% select(CHROM, START, END, SYMBOL, TYPE))
  
  df_tsb <- read_tsv(str_c("~/analysis/08_breast/03_somatic/", tn_pair, ".tsb_dv.vcf"), comment="##", show_col_types = FALSE) %>% rename(CHROM=`#CHROM`) %>% 
    mutate(csq1    = str_match(INFO, "(?<=CSQ1=)[^;]*(?=;)"),
           csq2    = str_match(INFO, "(?<=CSQ2=)[^;]*(?=;)"),
           symbol1 = sapply(csq1, function(x) str_split(x, "[|]", simplify = T)[2]),
           symbol2 = sapply(csq2, function(x) str_split(x, "[|]", simplify = T)[2]),
           TYPE    = str_extract(INFO, "(?<=SVTYPE=)[^;]*(?=;)"),
           POS1    = as.numeric(POS),
           POS2    = as.numeric(ifelse(TYPE=="BND", as.numeric(str_match(INFO, "(?<=POS2=)[^;]*(?=;)")), as.numeric(str_match(INFO, "(?<=END=)[^;]*(?=;)")))))
  
  df_tsb <- df_tsb %>% filter(!(TYPE=="DUP" & symbol1 != symbol2)) %>% filter(FILTER=="PASS")  # If it is a DUP, both breakpoints must be within the same gene to pass
  df_tsb <- bind_rows(
    df_tsb %>% select(CHROM, POS1, symbol1, TYPE) %>% rename(START=POS1, SYMBOL=symbol1),
    df_tsb %>% select(CHROM, POS2, symbol2, TYPE) %>% rename(START=POS2, SYMBOL=symbol2)
  )
  
  df_tsb <- df_tsb %>% mutate(END=START) %>% filter(SYMBOL %in% HR_genes)
  
  if(nrow(df_tsb)) df_merged <- bind_rows(df_merged, df_tsb %>% select(CHROM, START, END, SYMBOL, TYPE))
  
  df_cnv <-  read_tsv(str_c("~/analysis/08_breast/03_somatic/", tn_pair, ".cnv_dv.vcf"), comment="##", show_col_types = FALSE) 
  if(nrow(df_cnv)) { 
    df_cnv <-  df_cnv %>% rename(CHROM=`#CHROM`) %>% 
      mutate(SYMBOL = str_extract(INFO, "(?<=SYMBOL=)[^;]*(?=;)"),
             TYPE   = ifelse(str_detect(str_match(INFO, "(?<=GVM=)[^;]*(?=;)"), "AMP"), "AMP", "BIDEL"),
             START  = POS,
             END    = POS)    
    df_cnv <- df_cnv %>% 
      filter(FILTER=="PASS") %>% 
      filter(!str_detect(TYPE, "AMP")) %>% # Because only bidel is relevant for HR-related genes
      filter(SYMBOL %in% HR_genes)
    
    if(nrow(df_cnv)) df_merged <- bind_rows(df_merged, df_cnv %>% select(CHROM, START, END, SYMBOL, TYPE))  
  }

  if(nrow(df_merged)>0) {
    # Check loss of heterozygosity 
    df_seqz   <- read_tsv(f, show_col_types =FALSE) %>% rename(CHROM=`#CHROM`)
    gr_seqz   <- GRanges(seqnames = df_seqz$CHROM, ranges=IRanges(start=df_seqz$start_pos, end=df_seqz$end_pos), majCN=df_seqz$majCN, minCN=df_seqz$minCN)
    gr_merged <- GRanges(seqnames = df_merged$CHROM, ranges=IRanges(start=df_merged$START, end=df_merged$END))
    
    seqz_idx   <- subjectHits(findOverlaps(gr_merged, gr_seqz))
    merged_idx <- queryHits(findOverlaps(gr_merged, gr_seqz))

    if(length(seqz_idx)) df_merged$TYPE[merged_idx]  <- ifelse(gr_seqz[seqz_idx]$minCN == 0, paste(df_merged$TYPE[merged_idx], "wLOH", sep=" "), df_merged$TYPE[merged_idx]) # wLOH means "with loss of heterozygosity"
    df_merged <- df_merged %>% select(SYMBOL, TYPE) %>% distinct()
    df_merged <- df_merged %>% group_by(SYMBOL) %>% summarise(TYPE=paste0(TYPE, collapse=";")) 
  } else {
    df_merged <- tibble(SYMBOL=character(), TYPE=character())
  }
  
  df_s                       <- df_s %>% full_join(df_merged %>% dplyr::select(SYMBOL, TYPE))
  colnames(df_s)[ncol(df_s)] <- tumor_id
  
}

# Make table for germline mutations in homologous recombination pathway 
df_g <- tibble(SYMBOL=character())

for(f in list.files("/home/users/chrono0707/analysis/08_breast/06_seqz", pattern="*segments.clean.txt$", full.names = TRUE)) {
  
  normal_id <- str_split(basename(f) %>% str_remove(".segments.clean.txt"), "_")[[1]][2]
  tumor_id  <- str_split(basename(f) %>% str_remove(".segments.clean.txt"), "_")[[1]][1]
  
  if(!(tumor_id %in% df_crf$tumor_id)) next
  
  df_pm <- read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/02_germline/", "pm/", normal_id, ".pm.panel.fann.rare5.mvaf.rcfi.acmg.omim.ptchange.vcf"), comment="##", show_col_types = FALSE) %>% dplyr::rename(CHROM=`#CHROM`) %>% 
    mutate(csq     = str_match(INFO, "(?<=CSQ=)[^;]*(?=;)"),
           acmg    = str_match(INFO, "(?<=G_ACMG_class=)[^;]*(?=;)"),
           consequence = sapply(csq, function(x) str_split(x, "[|]", simplify = T)[2]),
           impact  = sapply(csq, function(x) str_split(x, "[|]", simplify = T)[3]),
           clinvar = sapply(csq, function(x) str_split(x, "[|]", simplify = T)[113]),
           SYMBOL  = sapply(csq, function(x) str_split(x, "[|]", simplify = T)[4]),
           TYPE    = "PM",
           START   = POS,
           END     = POS)
  
  df_pm <- df_pm %>% 
    filter(str_detect(acmg, "pathogenic|likely_pathogenic") | str_detect(impact, "HIGH") | str_detect(clinvar, "Pathogenic|Likely_pathogenic")) %>% 
    filter(!(str_detect(consequence, "intron_variant") & !str_detect(consequence, "splice_region_variant"))) %>%
    filter(!(str_detect(consequence, "upstream_gene_variant"))) 
  
  df_merged <- df_pm %>% dplyr::select(CHROM, START, END, SYMBOL, TYPE)
  
  tryCatch({ 
  
  df_delly <- read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/02_germline/", "delly/", normal_id, ".delly_germline_ALL.fann.pon1p.pcg.ptchange.vcf"), comment="##", show_col_types = FALSE) %>% dplyr::rename(CHROM=`#CHROM`) %>% 
    mutate(csq1    = str_match(INFO, "(?<=CSQ1=)[^;]*(?=;)"),
           csq2    = str_match(INFO, "(?<=CSQ2=)[^;]*(?=;)"),
           symbol1 = sapply(csq1, function(x) str_split(x, "[|]", simplify = T)[2]),
           symbol2 = sapply(csq2, function(x) str_split(x, "[|]", simplify = T)[2]),
           TYPE    = str_match(INFO, "(?<=SVTYPE=)[^;]*(?=;)"),
           POS1    = POS,
           POS2    = ifelse(TYPE=="BND", as.numeric(str_match(INFO, "(?<=POS2=)[^;]*(?=;)")), as.numeric(str_match(INFO, "(?<=END=)[^;]*(?=;)"))))

  df_delly <- df_delly %>% filter(!(TYPE=="DUP" & symbol1 != symbol2))
  
  df_delly <- bind_rows(
    df_delly %>% dplyr::select(CHROM, POS1, symbol1, TYPE) %>% dplyr::rename(START=POS1, SYMBOL=symbol1),
    df_delly %>% dplyr::select(CHROM, POS2, symbol2, TYPE) %>% dplyr::rename(START=POS2, SYMBOL=symbol2)
  )
  
  df_delly <- df_delly %>% mutate(END=START) 
  
  df_melt <-  read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/02_germline/", "melt/", normal_id, ".melt.fann.pass.rare1.interx.vcf"), comment="##", show_col_types = FALSE) %>% dplyr::rename(CHROM=`#CHROM`) %>% 
    mutate(csq    = str_match(INFO, "(?<=CSQte=)[^;]*(?=;)"),
           consequence = sapply(csq, function(x) str_split(x, "[|]", simplify = T)[2]),
           SYMBOL = sapply(csq, function(x) str_split(x, "[|]", simplify = T)[4]),
           TYPE   = str_match(INFO, "(?<=SVTYPE=)[^;]*(?=;)"),
           START  = POS,
           END    = POS)    
  df_melt <- df_melt %>% 
    filter(!str_detect(consequence, "intron_variant")) %>%
    filter(!str_detect(consequence, "upstream_gene_variant")) 
  
  df_merged <- bind_rows(
    df_merged,
    df_delly %>% dplyr::select(CHROM, START, END, SYMBOL, TYPE),
    df_melt  %>% dplyr::select(CHROM, START, END, SYMBOL, TYPE))
  
  }, error=function(e) {})
  
  if(nrow(df_merged)>0) {
    # Check loss of heterozygosity
    df_seqz   <- read_tsv(f, show_col_types =FALSE) %>% dplyr::rename(CHROM=`#CHROM`)
    gr_seqz   <- GRanges(seqnames = df_seqz$CHROM, ranges=IRanges(start=df_seqz$start_pos, end=df_seqz$end_pos), majCN=df_seqz$majCN, minCN=df_seqz$minCN)
    gr_merged <- GRanges(seqnames = df_merged$CHROM, ranges=IRanges(start=df_merged$START, end=df_merged$END))
    
    seqz_idx   <- subjectHits(findOverlaps(gr_merged, gr_seqz))
    merged_idx <- queryHits(findOverlaps(gr_merged, gr_seqz))
    
    if(length(seqz_idx)) df_merged$TYPE[merged_idx]  <- ifelse(gr_seqz[seqz_idx]$minCN == 0, paste(df_merged$TYPE[merged_idx], "wLOH", sep=" "), df_merged$TYPE[merged_idx])
    df_merged <- df_merged %>% dplyr::select(SYMBOL, TYPE) %>% distinct()
    df_merged <- df_merged %>% group_by(SYMBOL) %>% summarise(TYPE=paste0(TYPE, collapse=","))     
  } else {
    df_merged <- tibble(SYMBOL=character(), TYPE=character())
  }
  
  df_g                   <- df_g %>% full_join(df_merged %>% dplyr::select(SYMBOL, TYPE))
  colnames(df_g)[ncol(df_g)] <- tumor_id
  
}

df_g <- df_g %>% filter(SYMBOL!="")



# Draw figures 
# ------------

# Figure 2D. 
df_hrd %>% # 
  left_join(df_s %>% gather(tumor_id, s, 2:ncol(.)) %>% filter(!is.na(s)) %>% filter(str_detect(s, "LOH")) %>% filter(!str_detect(s, "BIDEL")) %>% filter(SYMBOL %in% c("BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CHEK2", "NBN", "PALB2", "RAD51", "RAD51B", "RAD51C", "RAD51D")) %>% mutate(s=str_c(SYMBOL, s, sep=" ")) %>% group_by(tumor_id) %>% summarise(s=paste0(s, collapse=","))) %>% 
  left_join(df_g %>% gather(tumor_id, g, 2:ncol(.)) %>% filter(!is.na(g)) %>% filter(str_detect(g, "LOH")) %>% filter(!str_detect(g, "BIDEL")) %>% filter(SYMBOL %in% c("BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CHEK2", "NBN", "PALB2", "RAD51", "RAD51B", "RAD51C", "RAD51D")) %>% mutate(g=str_c(SYMBOL, g, sep=" ")) %>% group_by(tumor_id) %>% summarise(g=paste0(g, collapse=","))) %>% 
  left_join(df_pam50) %>% 
  mutate(hrd=prop>0.20) %>% 
  mutate(g=ifelse(str_detect(g, "BRCA1"), "BRCA1", ifelse(str_detect(g, "BRCA2"), "BRCA2", "Others"))) %>% 
  mutate(s=ifelse(str_detect(s, "BRCA1"), "BRCA1", ifelse(str_detect(s, "BRCA2"), "BRCA2", "Others"))) %>%
  arrange(hrd, g, s) %>%
  mutate(tumor_id=factor(tumor_id, levels=tumor_id)) %>% 
  ggplot() +
  geom_col(aes(x=factor(1), y=factor(0)), fill=NA, color=NA) +
  geom_col(aes(x=factor(2), y=factor(0)), fill=NA, color=NA) + 
  geom_col(aes(x=factor(3), y=factor(0)), fill=NA, color=NA) + 
  geom_tile(aes(x=factor(4), y=tumor_id, fill=hrd)) + 
  geom_tile(aes(x=factor(5), y=tumor_id, fill=hrd)) + 
  geom_tile(aes(x=factor(6), y=tumor_id, fill=hrd)) + 
  geom_col(aes(x=factor(7), y=factor(0)), fill=NA, color=NA) + 
  geom_tile(aes(x=factor(8), y=tumor_id, fill=g)) + 
  geom_tile(aes(x=factor(9), y=tumor_id, fill=s)) +
  scale_fill_manual(values=c("FALSE"=rgb(241/255,241/255,241/255),
                             "TRUE" =rgb(44/255,110/255,179/255),
                             "BRCA1"=pal_locuszoom()(n=2)[1],
                             "BRCA2"=pal_locuszoom()(n=2)[2],
                             "Others"="grey50"),
                    na.value="white") +
  coord_polar(theta="y") + 
  theme_void() 

df_hrd %>% mutate(BIN=prop>0.20) %>% left_join(df_pam50) %>% dplyr::count(pam50, BIN) %>% filter(!is.na(pam50) & pam50!="Normal") %>% 
  ggplot(aes(x=factor(1), y=n, fill=BIN)) + geom_col(position="fill") + geom_col(aes(x=0, y=0)) + coord_polar(theta='y') + scale_fill_manual(values=c("TRUE"=pal_jco()(n=5)[1], "FALSE"=pal_jco()(n=5)[3])) + theme_void() + theme(legend.position="none") +
  facet_grid(rows=vars(pam50)) 

# Figure 2E. 
df_hrd %>% filter(prop>0.2) %>% left_join(df_g %>% filter(SYMBOL %in% c("BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CHEK2", "NBN", "PALB2", "RAD51", "RAD51B", "RAD51C", "RAD51D")) %>% gather(tumor_id, g, 2:ncol(.)) %>% filter(!is.na(g)) %>% filter(str_detect(g, "wLOH")) %>% filter(!str_detect(g, "BIDEL"))) %>% 
  mutate(g=ifelse(str_detect(g, "BND"), "BND", g)) %>% 
  mutate(g=ifelse(str_detect(g, "DUP"), "DUP", g)) %>% 
  mutate(g=ifelse(str_detect(g, "INV"), "INV", g)) %>% 
  mutate(g=ifelse(str_detect(g, "DEL"), "DEL", g)) %>% 
  mutate(g=ifelse(str_detect(g, "PM"), "PM", g)) %>% 
  dplyr::count(SYMBOL, g) %>% group_by(SYMBOL) %>% mutate(t=sum(n)) %>% ungroup() %>% arrange(desc(t)) %>% filter(!is.na(SYMBOL)) %>%  mutate(SYMBOL=factor(SYMBOL, levels=unique(SYMBOL))) %>% ggplot(aes(x=SYMBOL, y=n, fill=g)) + geom_bar(stat='identity') + 
  scale_fill_manual(
    values=c("DUP"= pal_jama()(4)[1],
             "BND"= pal_jama()(4)[2],
             "INV"= pal_jama()(4)[3],
             "DEL"= pal_jama()(4)[4],
             "PM"=pal_jama()(5)[5])
  ) + 
  theme_minimal() + scale_y_continuous(expand=c(0,0)) + theme(axis.text.x=element_text(angle=45, hjust=1.0, vjust=1.0), axis.title.x=element_blank(), panel.grid=element_blank(), axis.line.y=element_line(), axis.ticks.y = element_line(), legend.position='none')

df_hrd %>% filter(prop>0.2) %>%   
  filter(!(tumor_id %in% (df_g %>% filter(SYMBOL %in% c("BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CHEK2", "NBN", "PALB2", "RAD51", "RAD51B", "RAD51C", "RAD51D")) %>% gather(tumor_id, g, 2:ncol(.)) %>% filter(!is.na(g)) %>% filter(str_detect(g, "wLOH")) %>% filter(!str_detect(g, "BIDEL")) %>% dplyr::select(tumor_id) %>% distinct() %>% pull(tumor_id)))) %>%  # Remove samples with germline mutations
  left_join(df_s %>% filter(SYMBOL %in% c("BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CHEK2", "NBN", "PALB2", "RAD51", "RAD51B", "RAD51C", "RAD51D")) %>% gather(tumor_id, s, 2:ncol(.)) %>% filter(!is.na(s)) %>% filter(str_detect(s, "wLOH")) %>% filter(!str_detect(s, "BIDEL"))) %>% 
  mutate(s=ifelse(str_detect(s, "BND"), "BND", s)) %>% 
  mutate(s=ifelse(str_detect(s, "DUP"), "DUP", s)) %>% 
  mutate(s=ifelse(str_detect(s, "INV"), "INV", s)) %>% 
  mutate(s=ifelse(str_detect(s, "DEL"), "DEL", s)) %>% 
  mutate(s=ifelse(str_detect(s, "PM"), "PM", s)) %>% 
  dplyr::count(SYMBOL, s) %>% group_by(SYMBOL) %>% mutate(t=sum(n)) %>% ungroup() %>% arrange(desc(t)) %>% filter(!is.na(SYMBOL)) %>%  mutate(SYMBOL=factor(SYMBOL, levels=unique(SYMBOL))) %>% ggplot(aes(x=SYMBOL, y=n, fill=s)) + geom_bar(stat='identity') + 
  scale_fill_manual(
    values=c("DUP"= pal_jama()(4)[1],
             "BND"= pal_jama()(4)[2],
             "INV"= pal_jama()(4)[3],
             "DEL"= pal_jama()(4)[4],
             "PM"=pal_jama()(5)[5])
  ) + 
  theme_minimal() + scale_y_continuous(expand=c(0,0)) + theme(axis.text.x=element_text(angle=45, hjust=1.0, vjust=1.0), axis.title.x=element_blank(), panel.grid=element_blank(), axis.line.y=element_line(), axis.ticks.y = element_line(), legend.position='none')
