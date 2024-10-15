# Neoadjuvant chemotherapy in HER2-positive breast cancer

# Figure 4H. Pathologic complete remission (pCR) in HER2-positive breast cancer patients who received neoadjuvant docetaxel, carboplatin, trastuzumab, and pertuzumab (TCHP).
library(plotly)

df_temp <- df_crf %>% left_join(df_her2) %>% 
  filter(neoadj==1) %>% 
  filter(neoadj_regimen %in% c(3)) %>%  # neoadj_regimen 3 equals TCHP
  filter(!is.na(pCR))

fig <- plot_ly(
  type="sankey",
  domain = list(
    x=c(0,1),
    y=c(0,1)
  ),
  orientation = "h",
  arrangement="snap",
  
  node=list(
    pad = 10,
    thickness=20,
    list=list(
      color="red", 
      width=0.5
    )
  ),
  
  link=list(
    source=c(0,0,1,1),
    target=c(2,3,2,3),
    value  = c(
      df_temp %>% filter( focal & pCR==1) %>% nrow(),
      df_temp %>% filter( focal & pCR!=1) %>% nrow(),
      df_temp %>% filter(!focal & pCR==1) %>% nrow(),
      df_temp %>% filter(!focal & pCR!=1) %>% nrow()
    )    
    
  )
)

fig <- fig %>% layout(
  font = list(
    size=10
  ),
  xaxis = list(showgrid=F, zeroline=F),
  yaxis = list(showgrid=F, zeroline=F)
)

print(fig) 


# Figure 4I. Oncoplot illustrating clinicopathologic and genomic factors in patients with ERBB2 focal amplification who received neoadjuvant TCHP.

df_temp <- maf_drivers@data %>% as_tibble() %>% 
  filter(FILTER == "PASS") %>% 
  mutate(patient_id=str_sub(Tumor_Sample_Barcode, end=-6L)) %>% 
  dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Consequence) %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% summarise(variant_class=paste0(Variant_Classification, collapse=',')) %>% ungroup() %>% 
  filter(Tumor_Sample_Barcode %in% (df_crf %>% left_join(df_her2) %>% filter(focal & neoadj==1 & neoadj_regimen==3 & !is.na(pCR)) %>% pull(tumor_id))) %>% 
  spread(Tumor_Sample_Barcode, variant_class) %>% 
  rowwise(Hugo_Symbol) %>% mutate(count = sum(!is.na(c_across(everything())))) %>% ungroup()

m <- df_temp %>% 
  filter(count > ncol(.) * 0.03) %>% 
  dplyr::select(-Hugo_Symbol, -count) %>% 
  mutate_all(function(x) ifelse(is.na(x), "", x)) %>% 
  as.matrix()

m <- m[,df_crf %>% left_join(df_her2) %>% filter(tumor_id %in% colnames(m)) %>% filter(focal & neoadj==1 & !is.na(pCR)) %>% arrange(desc(pCR)) %>% pull(tumor_id)] # pCR+/pCR- 비교
rownames(m) <- df_temp %>% filter(count > ncol(.) * 0.03) %>% pull(Hugo_Symbol)

m_cna <- colnames(m) %>% enframe(value="tumor_id") %>% left_join(df_gistic %>% gather(tumor_id, bin, 6:ncol(.))) %>% mutate(bin=ifelse(bin==2, "Y", "N")) %>%  mutate(segment=str_c(chrom, ":", start, "-", end)) %>% dplyr::select(segment, tumor_id, bin) %>% spread(tumor_id, bin) %>% dplyr::select(colnames(m)) %>% as.matrix()
rownames(m_cna) <- colnames(m) %>% enframe(value="tumor_id") %>% left_join(df_gistic %>% gather(tumor_id, bin, 6:ncol(.))) %>% mutate(bin=ifelse(bin==2, "Y", "N")) %>%  mutate(segment=str_c(chrom, ":", start, "-", end)) %>% dplyr::select(segment, tumor_id, bin) %>% spread(tumor_id, bin) %>% pull(segment)

df_temp <- df_sig_snv_proportion %>% filter(tumor_id %in% colnames(m))
m_sig_sbs <- df_temp %>% dplyr::select(starts_with("v3")) %>% replace(is.na(.), 0) %>% as.matrix() 
rownames(m_sig_sbs) <- df_temp$tumor_id

df_temp <- df_sig_ind_proportion %>% filter(tumor_id %in% colnames(m)) 
m_sig_id <- df_temp %>% dplyr::select(starts_with("ID")) %>% replace(is.na(.), 0) %>% as.matrix()
rownames(m_sig_id) <- df_temp$tumor_id

df_temp <- df_sig_sv_proportion %>% filter(tumor_id %in% colnames(m)) 
m_sig_sv <- df_temp %>% dplyr::select(starts_with("RS")) %>% replace(is.na(.), 0) %>%  as.matrix() 
rownames(m_sig_sv) <- df_temp$tumor_id

palette_26colors <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "cyan", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown", "steelblue"
)

draw(
  oncoPrint(m, alter_fun = alter_fun, col=palette_variant, 
            pct_side = "right", row_names_side="left",
            remove_empty_columns = FALSE,
            
            column_split = factor(df_crf %>% left_join(df_her2) %>% filter(tumor_id %in% colnames(m)) %>% filter(focal & neoadj==1 & neoadj_regimen==3 & !is.na(pCR)) %>% arrange(desc(pCR)) %>% pull(pCR)), 
            column_gap=unit(0.3, "mm"),                                                                                                                                                 
            column_title_gp=gpar(fontsize=0),
            
            top_annotation = HeatmapAnnotation(
              tmb         = anno_barplot(df_tmb[match(colnames(m), df_tmb$tumor_id),] %>% dplyr::select(-tumor_id, -tmb) %>% as.matrix(), gp=gpar(col=NA, fill=c(pal_signal[[1]], pal_signal[[2]], pal_signal[[3]])), border = FALSE, bar_width=1),
              er          = df_crf[match(colnames(m), df_crf$tumor_id),] %>% pull(ER_diagnosis),
              pr          = df_crf[match(colnames(m), df_crf$tumor_id),] %>% pull(PR_diagnosis),
              ki67        = df_crf[match(colnames(m), df_crf$tumor_id),] %>% pull(ki67_diagnosis),
              tpm         = colnames(m) %>% enframe() %>% rename(tumor_id=value) %>% left_join(txi$abundance %>% as_tibble(rownames="symbol") %>% mutate(symbol=mapIds(org.Hs.eg.db, keys=substr(symbol, 1,15), column="SYMBOL", keytype="ENSEMBL", multiVals="first")) %>% filter(symbol=="ERBB2") %>% dplyr::select(-1) %>% gather(tumor_id, tpm) %>% mutate(tumor_id=str_replace(tumor_id, "R", "D")) %>% dplyr::slice(match(colnames(m), tumor_id))) %>% pull(tpm), 
              gcn         = df_her2[match(colnames(m), df_her2$tumor_id),] %>% pull(GCN),
              ploidy      = df_decision[match(colnames(m), df_decision$tumor_id),] %>% pull(ploidy),
              tnm         = df_crf[match(colnames(m), df_crf$tumor_id),] %>% pull(TNM),
              na_col="white",
              simple_anno_size=unit(0.2, "cm"),
              annotation_name_side = "left",
              annotation_name_gp = gpar(cex=0.5),
              col =list(
                er          = c("TRUE"=unikn::uni_ulm_2[[1]],  "FALSE"="gray90"),
                pr          = c("TRUE"=unikn::uni_ulm_2[[3]],  "FALSE"="gray90"),
                her2        = c("TRUE"=unikn::uni_ulm_2[[2]],  "FALSE"="gray90"),
                pam50       = c("Normal"="gray90", "LumA"=ggsci::pal_nejm()(5)[4], "LumB"=ggsci::pal_nejm()(5)[3], "Her2"=ggsci::pal_nejm()(5)[1], "Basal"=ggsci::pal_nejm()(5)[5]),
                ki67        = c("0"="gray90", "1"=brewer.pal(5, "GnBu")[2], "2"=brewer.pal(5, "GnBu")[3], "3"=brewer.pal(5, "GnBu")[4], "4"=brewer.pal(5, "GnBu")[5]),
                hrd         = colorRamp2(c(0,1), hcl_palette = "BuGn", reverse = T),
                ploidy      = colorRamp2(c(df_decision[match(colnames(m), df_decision$tumor_id),] %>% pull(ploidy) %>% min(), df_decision[match(colnames(m), df_decision$tumor_id),] %>% pull(ploidy) %>% max()), hcl_palette = "PuBu", reverse = T),
                math        = colorRamp2(c(df_math[match(colnames(m), df_math$tumor_id),] %>% pull(math) %>% min(), df_math[match(colnames(m), df_math$tumor_id),] %>% pull(math) %>% max()), hcl_palette = "Reds", reverse = T),
                tpm         = colorRamp2(c(colnames(m) %>% enframe() %>% rename(tumor_id=value) %>% left_join(txi$abundance %>% as_tibble(rownames="symbol") %>% mutate(symbol=mapIds(org.Hs.eg.db, keys=substr(symbol, 1,15), column="SYMBOL", keytype="ENSEMBL", multiVals="first")) %>% filter(symbol=="ERBB2") %>% dplyr::select(-1) %>% gather(tumor_id, tpm) %>% mutate(tumor_id=str_replace(tumor_id, "R", "D")) %>% dplyr::slice(match(colnames(m), tumor_id))) %>% pull(tpm) %>% min(na.rm=T), colnames(m) %>% enframe() %>% rename(tumor_id=value) %>% left_join(txi$abundance %>% as_tibble(rownames="symbol") %>% mutate(symbol=mapIds(org.Hs.eg.db, keys=substr(symbol, 1,15), column="SYMBOL", keytype="ENSEMBL", multiVals="first")) %>% filter(symbol=="ERBB2") %>% dplyr::select(-1) %>% gather(tumor_id, tpm) %>% mutate(tumor_id=str_replace(tumor_id, "R", "D")) %>% dplyr::slice(match(colnames(m), tumor_id))) %>% pull(tpm) %>% max(na.rm=T)), hcl_palette = "Oranges", reverse = T),
                gcn         = colorRamp2(c(df_her2[match(colnames(m), df_her2$tumor_id),] %>% pull(GCN) %>% min(na.rm=T), df_her2[match(colnames(m), df_her2$tumor_id),] %>% pull(GCN) %>% max(na.rm=T)), hcl_palette = "Heat", reverse = T),
                tnm         = c("IA"=unikn::caltech_pal_3[[1]], "IIA"=unikn::caltech_pal_3[[2]], "IIB"=unikn::caltech_pal_3[[3]], "IIIA"=unikn::caltech_pal_3[[4]], "IIIB"=unikn::caltech_pal_3[[5]], "IIIC"=unikn::caltech_pal_3[[6]])
              )              
            ),           
            right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(axis_param = list(side="top", labels_rot=0, gp=gpar(cex=0.5)))),
            bottom_annotation = HeatmapAnnotation(
              sbs = anno_barplot(m_sig_sbs, gp = gpar(fill=palette_26colors[1:ncol(m_sig_sbs)], col=NA), bar_width=1, height=unit(1.5, "cm"), border=FALSE, axis_param = list(gp=gpar(cex=0.5))),
              id  = anno_barplot(m_sig_id,  gp = gpar(fill=unikn::uni_regensburg_3[1:ncol(m_sig_id)], col=NA), bar_width=1, height=unit(1.5, "cm"), border=FALSE, axis_param = list(gp=gpar(cex=0.5))),
              sv  = anno_barplot(m_sig_sv,  gp = gpar(fill=unikn::uni_koeln_2[1:ncol(m_sig_sv)], col=NA), bar_width=1, height=unit(1.5, "cm"), border=FALSE, axis_param = list(gp=gpar(cex=0.5))),
              annotation_name_gp = gpar(cex=0.5),
              annotation_name_side = "left"
            ),
            show_column_names = F,
            row_names_gp = gpar(cex=0.7),
            pct_gp = gpar(cex=0.7)
  ),
  annotation_legend_list = list(
    Legend(labels=colnames(m_sig_sv), title="SV signature", type="grid", legend_gp = gpar(fill=uni_koeln_2[1:ncol(m_sig_sv)])),
    Legend(labels=colnames(m_sig_sbs), title="SBS signature", type="grid", legend_gp = gpar(fill=palette_26colors[1:ncol(m_sig_sbs)])),
    Legend(labels=colnames(m_sig_id), title="ID signature", type="grid", legend_gp = gpar(fill=unikn::uni_regensburg_3[1:ncol(m_sig_id)]))
  )
)

# calculate statistics 
df_crf %>% filter(tumor_id %in% colnames(m)) %>% filter(!is.na(ER_diagnosis)) %>% dplyr::count(ER_diagnosis, pCR) %>% pull(n) %>% matrix(ncol=2, byrow = T) %>% chisq.test()
df_crf %>% filter(tumor_id %in% colnames(m)) %>% filter(!is.na(PR_diagnosis)) %>% dplyr::count(PR_diagnosis, pCR) %>% pull(n) %>% matrix(ncol=2, byrow = T) %>% chisq.test()
df_crf %>% filter(tumor_id %in% colnames(m)) %>% filter(!is.na(PR_diagnosis) & !is.na(ER_diagnosis)) %>% mutate(hr=ER_diagnosis | PR_diagnosis) %>% dplyr::count(hr, pCR) %>% pull(n) %>% matrix(ncol=2, byrow = T) %>% chisq.test()
df_crf %>% filter(tumor_id %in% colnames(m)) %>% left_join(df_her2) %>% t.test(GCN~pCR, data=.)
df_crf %>% filter(tumor_id %in% colnames(m)) %>% left_join(txi$abundance %>% as_tibble(rownames="symbol") %>% mutate(symbol=mapIds(org.Hs.eg.db, keys=substr(symbol, 1,15), column="SYMBOL", keytype="ENSEMBL", multiVals="first")) %>% filter(symbol=="ERBB2") %>% dplyr::select(-1) %>% gather(tumor_id, tpm) %>% mutate(tumor_id=str_replace(tumor_id, "R", "D")) %>% dplyr::slice(match(colnames(m), tumor_id))) %>% t.test(tpm~pCR, data=.)
df_crf %>% filter(tumor_id %in% colnames(m)) %>% left_join(maf_drivers@data %>% as_tibble() %>% filter(Hugo_Symbol=="PIK3CA") %>% rename(tumor_id=Tumor_Sample_Barcode) %>% mutate(BIN=TRUE) %>% dplyr::select(tumor_id, BIN)) %>% mutate(BIN=ifelse(is.na(BIN), FALSE, BIN)) %>% dplyr::count(pCR, BIN) %>% pull(n) %>% matrix(ncol=2) %>% chisq.test()
