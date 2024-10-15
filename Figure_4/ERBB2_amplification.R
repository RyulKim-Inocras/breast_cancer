# ERBB2 amplification
# -------------------

# SV clusters associagted with ERBB2 amplification
df_cluster <- tibble() 
padding <- 1e6
for(f in list.files("/home/users/chrono0707/analysis/08_breast/03_somatic/", pattern="*cluster$", full.names = T)) {
  print(f)
  df_temp <- read_tsv(f, comment="##", show_col_types = F) %>% dplyr::rename(CHROM=`#CHROM`) %>% 
    mutate(tn_pair=basename(f) %>% str_remove(".SVs.quick_fi.tsv.cluster"),
           tumor_id=str_sub(tn_pair, end=-21L),
           POS2=ifelse(str_detect(ALT,"DUP|DEL|INV"), 
                       as.numeric(str_match(INFO, "(?<=END=)[^;]*(?=;)")),
                       as.numeric(str_match(INFO, "(?<=POS2=)[^;]*(?=;)"))),
           CHROM2=ifelse(str_detect(ALT,"DUP|DEL|INV"), 
                         CHROM,
                         str_match(INFO, "(?<=CHR2=)[^;]*(?=;)")),
           SVTYPE=ifelse(str_detect(ALT,"DUP|DEL|INV"),
                         str_remove(ALT, "<") %>% str_remove(">"),
                         "BND"),
           CSQ      = sapply(str_match(INFO, "(?<=CSQ1=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[1]), # consequence
           CSQ2     = sapply(str_match(INFO, "(?<=CSQ2=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[1]), # consequence
           SYMBOL   = sapply(str_match(INFO, "(?<=CSQ1=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[2]), # symbol (maybe, her2 or not)
           SYMBOL2  = sapply(str_match(INFO, "(?<=CSQ2=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[2]), # symbol 
           EXON     = sapply(str_match(INFO, "(?<=CSQ1=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[6]), # exon
           EXON2    = sapply(str_match(INFO, "(?<=CSQ2=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[6]), # exon
           INTRON   = sapply(str_match(INFO, "(?<=CSQ1=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[7]), # intron
           INTRON2  = sapply(str_match(INFO, "(?<=CSQ2=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[7]), # intron
           CT=as.character(str_match(INFO, "(?<=CT=)[^;]*(?=;)"))
    )
  
  gr_temp_csq1 <- GRanges(seqnames = df_temp$CHROM,  ranges = IRanges(start=df_temp$POS,  end=df_temp$POS ))
  gr_temp_csq2 <- GRanges(seqnames = df_temp$CHROM2, ranges = IRanges(start=df_temp$POS2, end=df_temp$POS2))
  
  suppressWarnings({clusters <- df_temp %>% dplyr::slice(
    unique(
      c(
        queryHits(findOverlaps(gr_temp_csq1, GRanges(seqnames = "chr17", ranges = IRanges(start=39700080 - padding, end=39728662 + padding)))),
        queryHits(findOverlaps(gr_temp_csq2, GRanges(seqnames = "chr17", ranges = IRanges(start=39700080 - padding, end=39728662 + padding))))
      )
    )
  ) %>% .$cluster %>% unique()
  })
  
  df_temp <- df_temp %>% filter(cluster %in% clusters)
  
  if(nrow(df_temp)>0) {
    df_cluster <- bind_rows(
      df_cluster,
      df_temp %>% select(tn_pair, tumor_id, CHROM, POS, CSQ, SYMBOL, EXON, INTRON, SVTYPE, CT, CHROM2, POS2, CSQ2, SYMBOL2, EXON2, INTRON2, cluster)
    )
  }
}
df_cluster <- df_cluster %>% filter(tumor_id %in% df_crf$tumor_id)

# Draw figures
# ------------

# Figure 4D. The segment size and copy number gain of ERBB2 amplicon.
df_her2 %>% 
  mutate(category=ifelse(focal, "focal", ifelse(delta_CN_padding_5mbp>=1, "broad", "no amp"))) %>% 
  ggplot(aes(x=width/1e6, y=delta_CN_padding_5mbp, color=category)) + geom_point() + geom_hline(yintercept = 3) + geom_vline(xintercept = 3) + theme_minimal() + theme(panel.grid=element_blank(), axis.ticks=element_line(), aspect.ratio = 1) + scale_color_jco() 

df_her2 %>% left_join(df_tmb) %>% left_join(df_crf) %>% 
  mutate(category=ifelse(focal, "focal", ifelse(delta_CN_padding_5mbp>=1, "broad", "no amp"))) %>% 
  dplyr::count(HER2, category) %>% filter(!is.na(HER2)) %>% 
  group_by(HER2) %>% mutate(n=n/sum(n)) %>% 
  ggplot(aes(factor(1), y=n, fill=category)) + geom_col(position='fill') + geom_col(aes(x=0, y=0)) + coord_polar(theta="y") + theme_void() + facet_grid(cols=vars(HER2)) + scale_fill_jco()


# Figure 4F. Distribution of structural variations (SVs) associated with ERBB2 focal amplification and breast-cancer-specific super-enhancers
df_temp <- df_cluster %>% left_join(df_her2) %>% mutate(END=POS, END2=POS2) %>% filter(focal)

circos.par("start.degree"=90)
circos.initializeWithIdeogram(species="hg38", plotType=c("ideogram", "labels"), chromosome.index=str_c("chr", c(1:22, "X")))

circos.genomicDensity(bind_rows(
  read_tsv("SE_01_0072_SE_hg38.bed", show_col_types = F) %>% dplyr::select(se_chr, se_start, se_end),   # breast epithelium
  read_tsv("SE_02_0990_SE_hg38.bed", show_col_types = F) %>% dplyr::select(se_chr, se_start, se_end),  # breast cancer cell line from brain metastasis
  read_tsv("SE_02_0991_SE_hg38.bed", show_col_types = F) %>% dplyr::select(se_chr, se_start, se_end) # breast cancer cell line from lung metastasis 
) %>% dplyr::select(se_chr, se_start, se_end), bg.border=NA, col=unikn::pal_signal[[1]], track.height=0.1, bg.lwd=0.5, bg.col="grey90")

circos.genomicTrack(
  bind_rows(
    read_tsv("SE_01_0072_SE_hg38.bed", show_col_types = F) %>% dplyr::select(se_chr, se_start, se_end),   # breast epithelium
    read_tsv("SE_02_0990_SE_hg38.bed", show_col_types = F) %>% dplyr::select(se_chr, se_start, se_end),  # breast cancer cell line from brain metastasis
    read_tsv("SE_02_0991_SE_hg38.bed", show_col_types = F) %>% dplyr::select(se_chr, se_start, se_end) # breast cancer cell line from lung metastasis 
  ) %>% mutate(value1=1) %>% dplyr::select(se_chr, se_start, se_end, value1), 
  ylim=c(0,1), 
  bg.border=NA, 
  bg.lwd = 0.5,
  track.height=0.05,
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, type="h", lwd=0.1, col=unikn::pal_signal[[1]], ...)
  })

circos.genomicDensity(bind_rows(
  df_temp %>% filter(!(CHROM == "chr17" & CHROM2 =="chr17")) %>% filter(SVTYPE=="BND") %>% dplyr::select(CHROM, POS), 
  df_temp %>% filter(!(CHROM == "chr17" & CHROM2 =="chr17")) %>% filter(SVTYPE=="BND") %>% dplyr::select(CHROM2, POS2) %>% dplyr::rename(CHROM=CHROM2, POS=POS2)) %>% mutate(END=POS) %>% filter(CHROM!="chr17"), 
  col = unikn::pal_unikn[[9]], bg.border=NA, bg.col="grey90",
  track.height = 0.1, bg.lwd=0.5)

circos.clear()

# Figure 4G. Permutation test demonstrating the distribution of SVs involved in ERBB2 focal amplification is closely related to that of breast-cancer-specific super-enhancers
library(regioneR)

df_se <- bind_rows(
  read_tsv("SE_01_0072_SE_hg38.bed") %>% dplyr::select(se_chr, se_start, se_end),   # breast epithelium
  read_tsv("SE_02_0990_SE_hg38.bed") %>% dplyr::select(se_chr, se_start, se_end),  # breast cancer cell line from brain metastasis
  read_tsv("SE_02_0991_SE_hg38.bed") %>% dplyr::select(se_chr, se_start, se_end) # breast cancer cell line from lung metastasis 
) %>% filter(!str_detect(se_chr, "alt"))

df_temp <- bind_rows(df_cluster %>% dplyr::select(CHROM, POS), df_cluster %>% dplyr::select(CHROM2, POS2) %>% dplyr::rename(CHROM=CHROM2, POS=POS2)) %>% 
  filter(!((CHROM == "chr17") & (POS>39700080 - padding) & (POS < 39728662 + padding))) %>% 
  filter(CHROM!="chrM")
gr_superenhancer <- GRanges(seqnames = df_se$se_chr,  ranges = IRanges(start = df_se$se_start, end = df_se$se_end))
gr_breakpoints   <- GRanges(seqnames = df_temp$CHROM, ranges = IRanges(start = df_temp$POS, end = df_temp$POS))
permutation_test <- permTest(A=gr_breakpoints, B=gr_superenhancer, ntimes=200, evaluate.function = meanDistance, randomize.function = randomizeRegions)

permutation_test$meanDistance$permuted %>% enframe(value='distance') %>% ggplot(aes(x=distance/1e6)) + geom_histogram(aes(y=after_stat(density)), bins = 200, width=1, fill='grey90', color=NA) + 
  geom_vline(xintercept = permutation_test$meanDistance$observed/1e6) + 
  geom_vline(xintercept = mean(permutation_test$meanDistance$permuted)/1e6, color="red") + 
  stat_function(fun = dnorm, args = list(mean = mean(permutation_test$meanDistance$permuted/1e6), sd = sd(permutation_test$meanDistance$permuted/1e6))) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme_minimal() + theme(axis.line.y=element_line(), axis.ticks.y=element_line(), axis.line.x=element_line(), axis.ticks.x=element_line())

