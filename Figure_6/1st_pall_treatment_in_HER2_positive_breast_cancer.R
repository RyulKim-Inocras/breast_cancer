# First-line palliative treatment of HER2-positive metastatic breast cancer
# -------------------------------------------------------------------------

# Figure 6A. A swimmer plot of pall 1st anti-her2 treatment 
plot_grid(
  #ER
  df_pall %>% filter(antiher2!="") %>% filter(!is.na(pall_1st_progression)) %>% mutate(pall_1st_progression=factor(pall_1st_progression, levels=unique(pall_1st_progression))) %>% arrange(pfs) %>% mutate(id=str_c("#",id), id=factor(id, levels=id)) %>% 
    ggplot(aes(x=id, y=0, fill=ER_diagnosis)) + geom_tile(color="white") + coord_flip() + scale_fill_manual(values=c("TRUE"=unikn::uni_ulm_2[[1]],  "FALSE"="gray90")) + theme_void() + theme(legend.position="none", axis.text.y=element_text(), axis.title.x=element_text(angle=90, vjust=.5)) + ylab("ER"),
  #PR
  df_pall %>% filter(antiher2!="") %>% filter(!is.na(pall_1st_progression)) %>% mutate(pall_1st_progression=factor(pall_1st_progression, levels=unique(pall_1st_progression))) %>% arrange(pfs) %>% mutate(tumor_id=factor(tumor_id, levels=tumor_id)) %>% 
    ggplot(aes(x=tumor_id, y=0, fill=PR_diagnosis)) + geom_tile(color="white") + coord_flip() + scale_fill_manual(values=c("TRUE"=unikn::uni_ulm_2[[3]],  "FALSE"="gray90")) + theme_void() + theme(legend.position="none", axis.title.x=element_text(angle=90, vjust=0.5)) + ylab("PR"),
  #HER2 TPM
  df_pall %>% filter(antiher2!="") %>% filter(!is.na(pall_1st_progression)) %>% mutate(pall_1st_progression=factor(pall_1st_progression, levels=unique(pall_1st_progression))) %>% arrange(pfs) %>% mutate(tumor_id=factor(tumor_id, levels=tumor_id)) %>% 
    left_join(txi$abundance %>% as_tibble(rownames="symbol") %>% mutate(symbol=mapIds(org.Hs.eg.db, keys=substr(symbol, 1,15), column="SYMBOL", keytype="ENSEMBL", multiVals="first")) %>% filter(symbol=="ERBB2") %>% dplyr::select(-1) %>% gather(tumor_id, tpm) %>% mutate(tumor_id=str_replace(tumor_id, "R", "D"))) %>% 
    ggplot(aes(x=tumor_id, y=0, fill=tpm)) + geom_tile(color="white") + coord_flip() + scale_fill_continuous_sequential(palette="Oranges", na.value="white") + theme_void() + theme(legend.position="none", axis.title.x=element_text(angle=90, vjust=0.5)) + ylab("TPM"),
  #GCN
  df_pall %>% left_join(df_her2) %>% filter(antiher2!="") %>% filter(!is.na(pall_1st_progression)) %>% mutate(pall_1st_progression=factor(pall_1st_progression, levels=unique(pall_1st_progression))) %>% arrange(pfs) %>% mutate(tumor_id=factor(tumor_id, levels=tumor_id)) %>% 
    ggplot(aes(x=tumor_id, y=0, fill=GCN)) + geom_tile(color="white") + coord_flip() + scale_fill_continuous_sequential(palette="Heat") + theme_void() + theme(legend.position="none", axis.title.x=element_text(angle=90, vjust=0.5)) + ylab("HER2 CN"),
  #KI67
  df_pall %>% filter(antiher2!="") %>% filter(!is.na(pall_1st_progression)) %>% mutate(pall_1st_progression=factor(pall_1st_progression, levels=unique(pall_1st_progression))) %>% arrange(pfs) %>% mutate(tumor_id=factor(tumor_id, levels=tumor_id)) %>% 
    ggplot(aes(x=tumor_id, y=0, fill=factor(ki67_diagnosis))) + geom_tile(color="white") + coord_flip() + scale_fill_manual(values=c("0"="gray90", "1"=brewer.pal(5, "GnBu")[2], "2"=brewer.pal(5, "GnBu")[3], "3"=brewer.pal(5, "GnBu")[4], "4"=brewer.pal(5, "GnBu")[5])) + theme_void() + theme(legend.position="none", axis.title.x=element_text(angle=90, vjust=0.5)) + ylab("Ki67"),
  #MATH
  df_pall %>% left_join(df_math) %>% filter(antiher2!="") %>% filter(!is.na(pall_1st_progression)) %>% mutate(pall_1st_progression=factor(pall_1st_progression, levels=unique(pall_1st_progression))) %>% arrange(pfs) %>% mutate(tumor_id=factor(tumor_id, levels=tumor_id)) %>% 
    ggplot(aes(x=tumor_id, y=0, fill=math)) + geom_tile(color="white") + coord_flip() + scale_fill_continuous_sequential(palette="Reds") + theme_void() + theme(legend.position="none", axis.title.x=element_text(angle=90, vjust=0.5)) + ylab("MATH"),
  #HRD
  df_pall %>% left_join(df_hrd) %>% filter(antiher2!="") %>% filter(!is.na(pall_1st_progression)) %>% mutate(pall_1st_progression=factor(pall_1st_progression, levels=unique(pall_1st_progression))) %>% arrange(pfs) %>% mutate(tumor_id=factor(tumor_id, levels=tumor_id)) %>% 
    ggplot(aes(x=tumor_id, y=0, fill=prop)) + geom_tile(color="white") + coord_flip() + scale_fill_continuous_sequential(palette="BuGn") + theme_void() + theme(legend.position="none", axis.title.x=element_text(angle=90, vjust=0.5)) + ylab("HRD score"),
  # ploidy
  df_pall %>% left_join(df_decision) %>% filter(antiher2!="") %>% filter(!is.na(pall_1st_progression)) %>% mutate(pall_1st_progression=factor(pall_1st_progression, levels=unique(pall_1st_progression))) %>% arrange(pfs) %>% mutate(tumor_id=factor(tumor_id, levels=tumor_id)) %>% 
    ggplot(aes(x=tumor_id, y=0, fill=ploidy)) + geom_tile(color="white") + coord_flip() + scale_fill_continuous_sequential(palette="PuBu") + theme_void() + theme(legend.position="none", axis.title.x=element_text(angle=90, vjust=0.5)) + ylab("Ploidy"),
  # regimen
  df_pall %>% left_join(df_decision) %>% filter(antiher2!="") %>% filter(!is.na(pall_1st_progression)) %>% mutate(pall_1st_progression=factor(pall_1st_progression, levels=unique(pall_1st_progression))) %>% arrange(pfs) %>% mutate(tumor_id=factor(tumor_id, levels=tumor_id)) %>% mutate(pall_1st_regimen=ifelse(str_detect(pall_1st_regimen, "trastuzumab|trasutuzmab") & str_detect(pall_1st_regimen, "pertuzumab"), "TP", ifelse(str_detect(pall_1st_regimen, "trastuzumab") & !str_detect(pall_1st_regimen, "pertuzumab"), "T", ifelse(str_detect(pall_1st_regimen, "lapatinib"), "L", pall_1st_regimen)))) %>% 
    ggplot(aes(x=tumor_id, y=0, fill=pall_1st_regimen)) + geom_tile(color="white") + coord_flip() + scale_fill_startrek() + theme_void() + theme(legend.position="none", axis.title.x=element_text(angle=90, vjust=0.5)) + ylab("Regimen"),
  # TP53
  df_pall %>% left_join(df_math) %>% filter(antiher2!="") %>% filter(!is.na(pall_1st_progression)) %>% mutate(pall_1st_progression=factor(pall_1st_progression, levels=unique(pall_1st_progression))) %>% arrange(pfs) %>% mutate(tumor_id=factor(tumor_id, levels=tumor_id)) %>% 
    mutate(pfs=pfs/365) %>% left_join(maf_drivers@data %>% as_tibble() %>% filter(Hugo_Symbol=="TP53") %>% dplyr::select(Tumor_Sample_Barcode) %>% dplyr::rename(tumor_id=Tumor_Sample_Barcode) %>% distinct() %>% mutate(TP53=TRUE)) %>% mutate(TP53=ifelse(is.na(TP53), FALSE, TP53)) %>% 
    ggplot(aes(x=tumor_id, y=0, fill=TP53)) + geom_tile(color="white") + coord_flip() + theme_void() + theme(legend.position="none", axis.title.x=element_text(angle=90, vjust=0.5)) + ylab("TP53") + scale_fill_manual(values=c("TRUE"="red", "FALSE"="white")),
  # PFS
  df_pall %>% left_join(df_math) %>% filter(antiher2!="") %>% filter(!is.na(pall_1st_progression)) %>% mutate(pall_1st_progression=factor(pall_1st_progression, levels=unique(pall_1st_progression))) %>% arrange(pfs) %>% mutate(tumor_id=factor(tumor_id, levels=tumor_id)) %>% 
    mutate(pfs=pfs/365) %>% 
    ggplot(aes(x=tumor_id, y=pfs, fill=pall_1st_progression)) + geom_col(color="white") + coord_flip() + theme_minimal() + theme(legend.position="none", axis.text.y=element_blank(), axis.title.y=element_blank()) + scale_fill_manual(values=c("0"=unikn::uni_bonn_1[[1]], "1"=unikn::uni_bonn_1[[2]])) + scale_y_continuous(expand=c(0,0), limits=c(0,13)) +
    geom_point(aes(x=tumor_id, y=pfs+100/365, shape=factor(pall_1st_reasone_for_discontinuation))) + scale_shape_manual(values=c("1"=4, "3"=5, "0"=1)),
  
  align="h", nrow=1, rel_widths = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 3))


# Figure 6B. survival by math score 
df_pall %>% left_join(df_math) %>% filter(antiher2!="") %>% mutate(bin=math>median(math), pall_1st_progression=ifelse(pall_1st_progression==1 & pfs>900, 0, pall_1st_progression), pfs=ifelse(pfs>900 | is.na(pfs), 900, pfs)) %>% coxph(Surv(pfs, pall_1st_progression)~bin+TNM, data=.) %>% summary() 
df_pall %>% left_join(df_math) %>% filter(antiher2!="") %>% mutate(bin=math>median(math), pall_1st_progression=ifelse(pall_1st_progression==1 & pfs>900, 0, pall_1st_progression), pfs=ifelse(pfs>900 | is.na(pfs), 900, pfs)) %>% survfit2(Surv(pfs/30, pall_1st_progression)~bin, data=.) %>% ggsurvfit() + scale_color_jco() + theme_minimal() + theme(axis.line = element_line(), axis.ticks = element_line(), panel.grid=element_blank())

df_pall %>% left_join(df_math) %>% filter(antiher2!="") %>% mutate(bin=math>median(math), pall_1st_progression=ifelse(pall_1st_progression==1 & pfs>900, 0, pall_1st_progression), pfs=ifelse(pfs>900 | is.na(pfs), 900, pfs)) %>% 
  left_join(maf_drivers@data %>% as_tibble() %>% filter(Hugo_Symbol!="TP53") %>% dplyr::select(Tumor_Sample_Barcode) %>% dplyr::rename(tumor_id=Tumor_Sample_Barcode) %>% distinct() %>% mutate(TP53=TRUE)) %>% mutate(TP53=ifelse(is.na(TP53), FALSE, TP53)) %>% 
  filter(TP53) %>% mutate(bin=math>median(math)) %>% 
  survfit2(Surv(pfs/30, pall_1st_progression)~bin, data=.) %>% ggsurvfit() + scale_color_jco() + theme_minimal() + theme(axis.line = element_line(), axis.ticks = element_line(), panel.grid=element_blank())
