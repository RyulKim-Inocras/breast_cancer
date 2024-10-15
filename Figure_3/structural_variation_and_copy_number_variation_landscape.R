# Structural variation landscape in breast cancer 
# -----------------------------------------------

library(GenomicRanges)
library(gUtils)

# Make table for structural variations across all samples 
df_sv_all <- tibble() 
for(f in list.files("../03_somatic", pattern="*SVs_quick_fi.vcf", full.names = T)) {
  df_temp <- read_tsv(f, comment="##", show_col_types = F) %>% filter(FILTER=="PASS")
  if(nrow(df_temp)==0) next
  df_temp <- df_temp %>% 
    mutate(
      tumor_id = basename(f) %>% str_sub(end=19),
      
      CHROM1   = `#CHROM`,
      POS1     = POS, 
      SYMBOL1  = sapply(str_match(INFO, "(?<=CSQ1=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[2]),
      EXON1    = sapply(str_match(INFO, "(?<=CSQ1=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[6]),
      INTRON1  = sapply(str_match(INFO, "(?<=CSQ1=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[7]),  
      
      CHROM2   = ifelse(str_detect(ALT, "DUP|DEL|INV"),
                        CHROM1,
                        str_extract(INFO, "(?<=CHR2=)[^;]+(?=;)")),
      POS2     = ifelse(str_detect(ALT, "DUP|DEL|INV"),
                        as.numeric(str_extract(INFO, "(?<=END=)[^;]+(?=;)")),
                        as.numeric(str_extract(INFO, "(?<=POS2=)[^;]+(?=;)"))),
      SYMBOL2   = sapply(str_match(INFO, "(?<=CSQ2=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[2]),
      EXON2     = sapply(str_match(INFO, "(?<=CSQ2=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[6]),
      INTRON2   = sapply(str_match(INFO, "(?<=CSQ2=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[7]),  
      CLASS     = ifelse(str_detect(ALT, "DUP|DEL|INV"),
                        str_remove_all(ALT, "<|>"),
                        "TRA"),
      CT        = as.character(str_extract(INFO, "(?<=CT=)[^;]*(?=;)"))
    )
  df_sv_all <- df_sv_all %>% bind_rows(df_temp %>% dplyr::select(tumor_id, CHROM1, POS1, SYMBOL1, CHROM2, POS2, SYMBOL2, CLASS, CT))
}

df_sv_all <- df_sv_all %>% filter(tumor_id %in% df_crf$tumor_id) # if not in sample list, remove them!

# Make table for SV drivers. 
df_sv <- tibble()
for(f in list.files("/home/users/chrono0707/analysis/08_breast/03_somatic", pattern="*tsb_dv.vcf$", full.names = T)) {
  df_temp <- read_tsv(f, comment="##", show_col_types = F) %>% rename(CHROM=`#CHROM`) %>%  
    filter(FILTER=="PASS") %>% 
    mutate(tn_pair=basename(f) %>% str_remove(".tsb_dv.vcf"),
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
           STRAND   = sapply(str_match(INFO, "(?<=CSQ1=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[12]), # strand 
           STRAND2  = sapply(str_match(INFO, "(?<=CSQ2=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[12]), # strand 
           CT=as.character(str_match(INFO, "(?<=;CT=)[^;]*(?=;)"))
    )
  
  if(nrow(df_temp)>0) {
    df_sv <- bind_rows(
      df_sv,
      df_temp %>% dplyr::select(tn_pair, tumor_id, CHROM, POS, CSQ, SYMBOL, EXON, INTRON, STRAND, SVTYPE, CT, CHROM2, POS2, CSQ2, SYMBOL2, EXON2, INTRON2, STRAND2)
    )
  }
}

for(f in list.files("/home/users/chrono0707/analysis/08_breast/03_somatic", pattern="*fsn_dv.vcf$", full.names = T)) {
  df_temp <- read_tsv(f, comment="##", show_col_types = F) %>% rename(CHROM=`#CHROM`) %>% 
    filter(FILTER=="PASS") %>% 
    mutate(tn_pair=basename(f) %>% str_remove(".fsn_dv.vcf"),
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
           STRAND   = sapply(str_match(INFO, "(?<=CSQ1=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[12]), # strand 
           STRAND2  = sapply(str_match(INFO, "(?<=CSQ2=)[^;]*(?=;)"), function(x) str_split(x, "[|]", simplify = T)[12]), # strand 
           CT=as.character(str_match(INFO, "(?<=;CT=)[^;]*(?=;)"))
    )
  
  if(nrow(df_temp)>0) {
    df_sv <- bind_rows(
      df_sv,
      df_temp %>% dplyr::select(tn_pair, tumor_id, CHROM, POS, CSQ, SYMBOL, EXON, INTRON, STRAND, SVTYPE, CT, CHROM2, POS2, CSQ2, SYMBOL2, EXON2, INTRON2, STRAND2)
    )
  }
}

# Count intra- and interchromosomal interactions between 5Mbp-sized windows across the entire genome
gr_tiles    <- gr.tile(Seqinfo(genome="hg38"), 5e6) %>% as_tibble() %>% filter(seqnames %in% str_c("chr", c(1:22, "X", "Y"))) %>% makeGRangesFromDataFrame()
df_sv_tiles <- gr_tiles %>% as_tibble() %>% transpose() %>% map(function(x) {
  return(
    gr_tiles %>% as_tibble() %>% filter(row_number()>=which(seqnames==x[1] & start==as.numeric(x[2]) & end==as.numeric(x[3]))) %>% transpose() %>% map(function(y) {
      print(x)
      print(y)
      n1 <- df_sv_all %>% filter(CHROM1==x[1] & POS1 >= as.numeric(x[2]) & POS1 < as.numeric(x[3]) & CHROM2==y[1] & POS2 >= as.numeric(y[2]) & POS2 < as.numeric(y[3])) %>% count(tumor_id) %>% nrow()
      n2 <- df_sv_all %>% filter(CHROM1==y[1] & POS1 >= as.numeric(y[2]) & POS1 < as.numeric(y[3]) & CHROM2==x[1] & POS2 >= as.numeric(x[2]) & POS2 < as.numeric(x[3])) %>% count(tumor_id) %>% nrow()
      n <- if(x[1]==y[1] & x[2]==y[2] & x[3]==y[3]) (n1+n2)/2 else (n1+n2)
      print(n)
      return(tibble(CHROM1=x[1], START1=x[2], END1=x[3], CHROM2=y[1], START2=y[2], END2=y[3], n=n))
    }) %>% do.call(bind_rows, .)
  ) 
}) %>% do.call(bind_rows, .)



# Draw figures
# ------------

# Figure 3A. The distribution of SV breakpoints and recurrent CNVs.

# SV pyramid
df_sv_tiles %>%  
  filter(!(CHROM1 %in% c("chrX", "chrY"))) %>% 
  filter(!(CHROM2 %in% c("chrX", "chrY"))) %>% 
  mutate(bp1=str_c(CHROM1, ":", START1, "-", END1), bp2=str_c(CHROM2, ":", START2, "-", END2)) %>% mutate(bp1=factor(bp1, levels=unique(bp1)), bp2=factor(bp2, levels=unique(bp2))) %>% 
  mutate(n=ifelse(n>20, 20, n)) %>% 
  mutate(n=ifelse(n==0, NA, n)) %>% 
  group_by(bp1) %>% filter(row_number() >= which(bp1==bp2)) %>% 
  ggplot(aes(x=bp1, y=bp2, fill=n)) + 
  geom_tile() +
  scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 30, l2 = 100, p1 = .9, p2 = 1.2, na.value=NA) + 
  scale_x_discrete(expand=c(0,0), breaks=gr_tiles %>% as_tibble() %>% filter(!(seqnames %in% c("chrX", "chrY"))) %>% mutate(segment=str_c(seqnames, ":", start, "-", end), segment=factor(segment, levels=segment)) %>% pull(segment)) + 
  scale_y_discrete(expand=c(0,0), breaks=gr_tiles %>% as_tibble() %>% filter(!(seqnames %in% c("chrX", "chrY"))) %>% mutate(segment=str_c(seqnames, ":", start, "-", end), segment=factor(segment, levels=segment)) %>% pull(segment)) + 
  theme_void() + theme(aspect.ratio = 1, panel.background = element_rect(color="black", fill=NA))

# Breakpoints density of SV drivers 
df_temp1 <- gr_tiles %>% as_tibble() %>% transpose() %>% map(function(x) {
  df <- bind_rows(
    df_sv %>% filter(CHROM ==str_c("chr", x[1])  & POS  >= as.integer(x[2])  & POS < as.integer(x[3])) %>% select(tumor_id, SYMBOL),
    df_sv %>% filter(CHROM2==str_c("chr", x[1])  & POS2 >= as.integer(x[2]) & POS2 < as.integer(x[3])) %>% select(tumor_id, SYMBOL2) %>% rename(SYMBOL=SYMBOL2)
  ) %>% distinct() %>% left_join(df_pam50)
  if(nrow(df)) df %>% count(pam50) %>% mutate(CHROM=str_c("chr", x[1]), START=as.integer(x[2]), END=as.integer(x[3]), genes=paste0(unique(df %>% filter(!is.na(SYMBOL)) %>% pull(SYMBOL)), collapse=",")) %>% select(CHROM, START, END, genes, pam50, n) %>% return()
  else tibble(CHROM=str_c("chr", x[1]), START=as.integer(x[2]), END=as.integer(x[3]), genes=NA, pam50=NA, n=0)
}) %>% do.call(bind_rows, .) %>% mutate(segment=str_c(CHROM, ":", START, "-", END)) %>% mutate(segment=factor(segment, levels=unique(segment))) %>% 
  filter(!(CHROM %in% c("chrX", "chrY"))) %>% 
  mutate(pam50=ifelse(is.na(pam50), "Unknown", as.character(pam50)), pam50=factor(pam50, levels=c("Unknown", "LumA", "LumB", "Her2", "Basal")))

ggplot() + 
  geom_bar(data=df_temp1 %>% mutate(odds=ifelse((as.numeric(str_remove(CHROM, "chr")) %% 2)==0, FALSE, TRUE)) %>% dplyr::select(segment, odds) %>% distinct(), aes(x=segment, fill=odds, y=300), width=1.0, color=NA, stat='identity') + 
  geom_bar(data=df_temp1, aes(x=segment, y=n, fill=pam50), width=1, stat='identity') + 
  scale_fill_manual(values=c(c("TRUE"="grey90", "FALSE"="white", "Normal"="gray90", "LumA"=ggsci::pal_nejm()(5)[4], "LumB"=ggsci::pal_nejm()(5)[3], "Her2"=ggsci::pal_nejm()(5)[1], "Basal"=ggsci::pal_nejm()(5)[5]))) +
  scale_x_discrete(expand=c(0,0), breaks=gr_tiles %>% as_tibble() %>% filter(!(seqnames %in% c("chrX", "chrY"))) %>% mutate(segment=str_c(seqnames, ":", start, "-", end), segment=factor(segment, levels=segment)) %>% pull(segment)) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme_minimal() + theme(legend.position='none', panel.grid=element_blank(), axis.text.x=element_blank(), axis.line.y=element_line(), axis.ticks.y=element_line()) 


# GISTIC output
gistic <- readGistic(gisticAllLesionsFile="../10_GISTIC/01_focal/all_lesions.conf_90.txt",
                     gisticAmpGenesFile  ="../10_GISTIC/01_focal/amp_genes.conf_90.txt",
                     gisticDelGenesFile  ="../10_GISTIC/01_focal/del_genes.conf_90.txt",
                     gisticScoresFile    ="../10_GISTIC/01_focal/scores.gistic")


g = getCytobandSummary(gistic)
g = g[qvalues < 0.10]
g[, `:=`(Chromosome, sapply(strsplit(x = g$Wide_Peak_Limits, 
                                     split = ":"), "[", 1))]
g[, `:=`(loc, sapply(strsplit(x = g$Wide_Peak_Limits, split = ":"), 
                     "[", 2))]
g[, `:=`(Start_Position, sapply(strsplit(x = g$loc, split = "-"), 
                                "[", 1))]
g[, `:=`(End_Position, sapply(strsplit(x = g$loc, split = "-"), 
                              "[", 2))]

transformSegments = function(segmentedData, build = 'hg19'){
  
  build.opts = c('hg19', 'hg18', 'hg38')
  
  if(!build %in% build.opts){
    stop('Available reference builds: hg18, hg19, hg38')
  }
  
  if(build == 'hg19'){
    chr.lens = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                 155270560, 59373566)
  } else if(build == 'hg18'){
    chr.lens = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
                 158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
                 114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
                 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  } else if(build == 'hg38'){ #hg38
    chr.lens = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  } else{
    stop('Available reference builds: hg18, hg19, hg38')
  }
  
  segmentedData[,Start_Position := as.numeric(as.character(Start_Position))]
  segmentedData[,End_Position := as.numeric(as.character(End_Position))]
  
  #Replace chr x and y with numeric value (23 and 24) for better ordering
  segmentedData$Chromosome = gsub(pattern = 'chr', replacement = '', x = segmentedData$Chromosome, fixed = TRUE)
  segmentedData$Chromosome = gsub(pattern = 'X', replacement = '23', x = segmentedData$Chromosome, fixed = TRUE)
  segmentedData$Chromosome = gsub(pattern = 'Y', replacement = '24', x = segmentedData$Chromosome, fixed = TRUE)
  
  segmentedData$Chromosome = factor(x = segmentedData$Chromosome, levels = 1:24, labels = 1:24)
  
  segmentedData = segmentedData[order(Chromosome, Start_Position, decreasing = FALSE)]
  
  seg.spl = split(segmentedData, segmentedData$Chromosome)
  
  seg.spl.transformed = seg.spl[[1]]
  if(nrow(seg.spl.transformed) > 0){
    seg.spl.transformed$Start_Position_updated = seg.spl.transformed$Start_Position
    seg.spl.transformed$End_Position_updated = seg.spl.transformed$End_Position
  }
  
  chr.lens.sumsum = cumsum(chr.lens)
  
  for(i in 2:length(seg.spl)){
    
    x.seg = seg.spl[[i]]
    if(nrow(x.seg) > 0){
      x.seg$Start_Position_updated = x.seg$Start_Position + chr.lens.sumsum[i-1]
      x.seg$End_Position_updated = x.seg$End_Position + chr.lens.sumsum[i-1]
    }
    seg.spl.transformed = rbind(seg.spl.transformed, x.seg, fill = TRUE)
  }
  
  return(seg.spl.transformed)
}

gis.scores = transformSegments(segmentedData = gistic@gis.scores, 
                               build = "hg38")
gis.scores$amp = ifelse(test = gis.scores$Variant_Classification == 
                          "Del", yes = -gis.scores$G_Score, no = gis.scores$G_Score)
gis.scores$Variant_Classification = ifelse(test = as.numeric(gis.scores$fdr) > 
                                             -log10(0.10), yes = gis.scores$Variant_Classification, no = "neutral")
gis.scores$Variant_Classification = factor(gis.scores$Variant_Classification, 
                                           levels = c("neutral", "Amp", "Del"))

df_temp2 <- gr_tiles %>% as_tibble() %>% transpose() %>% map(function(x) {
  print(x)
  df <- gis.scores %>% tibble() %>% filter(Chromosome==str_remove(as.character(x[1]), "chr") & Start_Position >= as.numeric(x[2]) & End_Position < as.numeric(x[3])) %>% 
    group_by(Variant_Classification) %>% summarise(amp=mean(amp)) %>% 
    mutate(CHROM=x[1], START=as.numeric(x[2]), END=as.numeric(x[3])) %>% dplyr::select(CHROM, START, END, Variant_Classification, amp) 
  if(nrow(df)) return(df) else return(tibble(CHROM=x[1], START=as.numeric(x[2]), END=as.numeric(x[3]), Variant_Classification=c("neutral", "Amp", "Del"), amp=c(0,0,0)))
}) %>% do.call(bind_rows, .) %>% mutate(segment=str_c(CHROM, ":", START, "-", END)) %>% mutate(segment=factor(segment, levels=unique(segment))) %>% 
  filter(!(CHROM %in% c("chrX", "chrY")))

ggplot() + 
  geom_bar(data=df_temp2 %>% mutate(odds=ifelse((as.numeric(str_remove(CHROM, "chr")) %% 2)==0, FALSE, TRUE)) %>% dplyr::select(segment, odds) %>% distinct(), aes(x=segment, fill=odds, y= 1.0), width=1, color=NA, stat='identity') + 
  geom_bar(data=df_temp2 %>% mutate(odds=ifelse((as.numeric(str_remove(CHROM, "chr")) %% 2)==0, FALSE, TRUE)) %>% dplyr::select(segment, odds) %>% distinct(), aes(x=segment, fill=odds, y=-1.0), width=1, color=NA, stat='identity') + 
  scale_fill_manual(values=c("TRUE"="grey90", "FALSE"="white")) + 
  geom_bar(data=df_temp2 %>% filter(Variant_Classification=="neutral"), aes(x=segment, y=amp), stat='identity', fill="grey70", width=1) + 
  geom_bar(data=df_temp2 %>% filter(Variant_Classification=="Amp"),     aes(x=segment, y=amp), stat='identity', fill=unikn::uni_freiburg_br[[2]], width=1) +  
  geom_bar(data=df_temp2 %>% filter(Variant_Classification=="Del"),     aes(x=segment, y=amp), stat='identity', fill=unikn::uni_freiburg_br[[1]], width=1) + 
  geom_hline(yintercept=0) + 
  scale_x_discrete(expand=c(0,0), breaks=gr_tiles %>% as_tibble() %>% filter(!(seqnames %in% c("chrX", "chrY"))) %>% mutate(segment=str_c(seqnames, ":", start, "-", end), segment=factor(segment, levels=segment)) %>% pull(segment)) + 
  scale_y_continuous(expand=c(0,0), limits = c(-1.0, 1.0)) + 
  theme_minimal() + theme(legend.position='none', panel.grid=element_blank(), axis.text.x=element_blank(), axis.line.y=element_line(), axis.ticks.y=element_line())


# Zoom in the interactions between chr8 and chr11
df_sv_tiles %>% 
  filter(CHROM1 %in% c("chr8", "chr11")) %>% 
  filter(CHROM2 %in% c("chr8", "chr11")) %>% 
  mutate(bp1=str_c(CHROM1, ":", START1, "-", END1), bp2=str_c(CHROM2, ":", START2, "-", END2)) %>% mutate(bp1=factor(bp1, levels=unique(bp1)), bp2=factor(bp2, levels=unique(bp2))) %>% 
  mutate(n=ifelse(n>20, 20, n)) %>% 
  mutate(n=ifelse(n==0, NA, n)) %>% 
  group_by(bp1) %>% filter(row_number() >= which(bp1==bp2)) %>% 
  ggplot(aes(x=bp1, y=bp2, fill=n)) + 
  geom_tile() +
  geom_vline(xintercept="chr8:35000001-40000000",  linetype="dashed", size=0.1) +  # FGFR1
  geom_vline(xintercept="chr11:65000001-70000000", linetype="dashed", size=0.1) +  # CCND1
  geom_hline(yintercept="chr8:35000001-40000000",  linetype="dashed", size=0.1) +  # FGFR1
  geom_hline(yintercept="chr11:65000001-70000000", linetype="dashed", size=0.1) +  # CCND1
  geom_hline(yintercept="chr11:75000001-80000000", linetype="dashed", size=0.1) +  # CCND1
  scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 30, l2 = 100, p1 = .9, p2 = 1.2, na.value=NA) + 
  scale_x_discrete(expand=c(0,0), breaks=gr_tiles %>% as_tibble() %>% filter(!(seqnames %in% c("chrX", "chrY"))) %>% mutate(segment=str_c(seqnames, ":", start, "-", end), segment=factor(segment, levels=segment)) %>% pull(segment)) + 
  scale_y_discrete(expand=c(0,0), breaks=gr_tiles %>% as_tibble() %>% filter(!(seqnames %in% c("chrX", "chrY"))) %>% mutate(segment=str_c(seqnames, ":", start, "-", end), segment=factor(segment, levels=segment)) %>% pull(segment)) + 
  theme_void() + theme(aspect.ratio = 1, panel.background = element_rect(color="black", fill=NA), legend.position="none")


# Zoom in the interactions between chr17 and chr20
df_sv_tiles %>% 
  filter(CHROM1 %in% c("chr17", "chr20")) %>% 
  filter(CHROM2 %in% c("chr17", "chr20")) %>% 
  mutate(bp1=str_c(CHROM1, ":", START1, "-", END1), bp2=str_c(CHROM2, ":", START2, "-", END2)) %>% mutate(bp1=factor(bp1, levels=unique(bp1)), bp2=factor(bp2, levels=unique(bp2))) %>% 
  mutate(n=ifelse(n>20, 20, n)) %>% 
  mutate(n=ifelse(n==0, NA, n)) %>% 
  group_by(bp1) %>% filter(row_number() >= which(bp1==bp2)) %>% 
  ggplot(aes(x=bp1, y=bp2, fill=n)) + 
  geom_tile() +
  geom_vline(xintercept="chr17:35000001-40000000", linetype="dashed", size=0.1) +  # ERBB2
  geom_vline(xintercept="chr17:65000001-70000000", linetype="dashed", size=0.1) +  # ERBB2
  geom_hline(yintercept="chr20:50000001-55000000", linetype="dashed", size=0.1) +  # ERBB2
  geom_hline(yintercept="chr20:60000001-64444167", linetype="dashed", size=0.1) +  # ERBB2
  scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 30, l2 = 100, p1 = .9, p2 = 1.2, na.value=NA) + 
  scale_x_discrete(expand=c(0,0), breaks=gr_tiles %>% as_tibble() %>% filter(!(seqnames %in% c("chrX", "chrY"))) %>% mutate(segment=str_c(seqnames, ":", start, "-", end), segment=factor(segment, levels=segment)) %>% pull(segment)) + 
  scale_y_discrete(expand=c(0,0), breaks=gr_tiles %>% as_tibble() %>% filter(!(seqnames %in% c("chrX", "chrY"))) %>% mutate(segment=str_c(seqnames, ":", start, "-", end), segment=factor(segment, levels=segment)) %>% pull(segment)) + 
  theme_void() + theme(aspect.ratio = 1, panel.background = element_rect(color="black", fill=NA), legend.position="none")

# Superenhancers on chr20
gr_tiles %>% as_tibble() %>% filter(seqnames %in% c("chr17", "chr20")) %>% transpose() %>% map(function(x) {
  n <- bind_rows(
    read_tsv("SE_01_0072_SE_hg38.bed", show_col_types = F) %>% dplyr::select(se_chr, se_start, se_end),   # breast epithelium
    read_tsv("SE_02_0990_SE_hg38.bed", show_col_types = F) %>% dplyr::select(se_chr, se_start, se_end),  # breast cancer cell line from brain metastasis
    read_tsv("SE_02_0991_SE_hg38.bed", show_col_types = F) %>% dplyr::select(se_chr, se_start, se_end) # breast cancer cell line from lung metastasis 
  ) %>% filter(se_chr==str_c("chr", x[1]) & se_start >= as.numeric(x[2]) & se_end < as.numeric(x[3])) %>% nrow()
  tibble(CHROM=x[1], START=as.numeric(x[2]), END=as.numeric(x[3]), n=n) %>% return()
}) %>% do.call(bind_rows, .) %>% 
  mutate(segment=str_c(CHROM, ":", START, '-', format(END, scientific=F, trim=T))) %>% mutate(segment=factor(segment, levels=segment)) %>% 
  ggplot(aes(x=segment, y=n)) + geom_col(width=1) + theme_minimal() + theme(axis.text.x=element_blank(), panel.grid=element_blank(), axis.line.y=element_line(), axis.ticks.y=element_line()) + 
  scale_x_discrete(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))


# Figure 3B. Key cancer driver genes affected by SVs
df_temp <- bind_rows(df_sv %>% dplyr::select(tumor_id, SYMBOL), df_sv %>% dplyr::select(tumor_id, SYMBOL2) %>% dplyr::rename(SYMBOL=SYMBOL2)) %>% distinct() %>% filter(!is.na(SYMBOL)) %>% 
  group_by(SYMBOL) %>% mutate(n=n()) %>% ungroup() %>% 
  filter(n/length(unique(tumor_id))>0.02) %>% 
  arrange(n) %>% mutate(SYMBOL=factor(SYMBOL, levels=unique(SYMBOL))) %>% 
  left_join(df_pam50) %>% dplyr::count(SYMBOL, pam50) %>% 
  mutate(pam50=as.character(pam50), pam50=ifelse(pam50=="Normal", NA, pam50), pam50=factor(pam50, levels=c("LumA", "LumB", "Her2", "Basal"))) 

g1 <- df_temp %>% 
  ggplot(aes(x=SYMBOL, y=n, fill=pam50)) + geom_bar(stat='identity', width=1) + theme_minimal() + 
  theme(legend.position="none", axis.ticks.x=element_line()) + 
  scale_y_continuous(expand=c(0,0)) + 
  scale_fill_manual(values=c("Normal"="gray90", "LumA"=ggsci::pal_nejm()(5)[4], "LumB"=ggsci::pal_nejm()(5)[3], "Her2"=ggsci::pal_nejm()(5)[1], "Basal"=ggsci::pal_nejm()(5)[5]), na.value="grey90") + coord_flip()

g2 <- bind_rows(df_sv %>% dplyr::select(tumor_id, SYMBOL, SVTYPE), df_sv %>% dplyr::select(tumor_id, SYMBOL2, SVTYPE) %>% dplyr::rename(SYMBOL=SYMBOL2)) %>% distinct() %>% filter(!is.na(SYMBOL)) %>% 
  dplyr::count(SYMBOL, SVTYPE) %>% 
  filter(SYMBOL %in% df_temp$SYMBOL) %>% 
  mutate(SYMBOL=factor(SYMBOL, levels=levels(df_temp$SYMBOL))) %>% 
  group_by(SYMBOL) %>% mutate(p=n/sum(n)) %>% 
  ggplot(aes(x=SYMBOL, y=p, fill=SVTYPE)) + geom_col(position='fill', width=1, color=NA) + theme_minimal() + 
  scale_fill_manual(values=c(
    "DUP"                    = pal_jama()(4)[1],
    "BND"                    = pal_jama()(4)[2],
    "INV"                    = pal_jama()(4)[3],
    "DEL"                    = pal_jama()(4)[4]    
  )) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(axis.ticks.x=element_line(), legend.position="none", axis.text.y=element_blank(), axis.title.y=element_blank()) + coord_flip()

plot_grid(g1, g2, align="h", nrow=1, rel_widths = c(0.85, 0.15))  


# Figure 3D. Kaplan-Meier survival curves for overall survival, stratified by TP53 and ETV6 mutation status.
df_temp <- bind_rows(df_sv %>% dplyr::select(tumor_id, SYMBOL), df_sv %>% dplyr::select(tumor_id, SYMBOL2) %>% dplyr::rename(SYMBOL=SYMBOL2)) %>% distinct() %>% filter(!is.na(SYMBOL)) %>% 
  filter(SYMBOL=="ETV6") %>% mutate(BIN=TRUE) %>% 
  right_join(df_crf) %>% mutate(BIN=ifelse(is.na(BIN), FALSE, BIN)) %>% 
  left_join(maf_drivers@data %>% as_tibble() %>% filter(Hugo_Symbol=="TP53") %>% dplyr::select(Tumor_Sample_Barcode) %>% distinct() %>% mutate(tp53=TRUE) %>% dplyr::rename(tumor_id=Tumor_Sample_Barcode)) %>% mutate(tp53=ifelse(is.na(tp53), FALSE, tp53))

df_temp %>% survfit2(Surv(os/365, os_event)~BIN+tp53, data=.) %>% ggsurvfit() +
  scale_color_jco() + 
  theme_minimal() + theme(axis.line = element_line(), axis.ticks = element_line(), panel.grid=element_blank())

df_temp %>% filter(tp53) %>% coxph(Surv(os, os_event)~BIN, data=.) %>% summary()
df_temp %>% filter(!tp53) %>% coxph(Surv(os, os_event)~BIN, data=.)
