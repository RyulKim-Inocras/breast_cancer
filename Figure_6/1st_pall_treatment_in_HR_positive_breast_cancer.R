# First-line palliative treatment of hormone receptor-positive metastatic breast cancer
# -------------------------------------------------------------------------------------


# Figure 6C. 1st line pall treatment of CDK4/6 inhibitor
df_temp <- maf_drivers@data %>% as_tibble() %>% 
  filter(FILTER == "PASS") %>% 
  mutate(patient_id=str_sub(Tumor_Sample_Barcode, end=-6L)) %>% 
  dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Consequence) %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% summarise(variant_class=paste0(Variant_Classification, collapse=',')) %>% ungroup() %>% 
  filter(Tumor_Sample_Barcode %in% (df_pall %>% filter(!is.na(cdki) & !is.na(pall_1st_progression)) %>% pull(tumor_id))) %>% 
  spread(Tumor_Sample_Barcode, variant_class) %>% 
  rowwise(Hugo_Symbol) %>% mutate(count = sum(!is.na(c_across(everything())))) %>% ungroup()

m <- df_temp %>% 
  filter(count > ncol(.) * 0.03) %>% 
  dplyr::select(-Hugo_Symbol, -count) %>% 
  mutate_all(function(x) ifelse(is.na(x), "", x)) %>% 
  as.matrix()

m <- m[,df_pall %>% filter(tumor_id%in%colnames(m)) %>% filter(!is.na(cdki) & !is.na(pall_1st_progression)) %>% arrange(desc(pall_1st_progression)) %>% pull(tumor_id)] 
rownames(m) <- df_temp %>% filter(count > ncol(.) * 0.03) %>% pull(Hugo_Symbol)

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
            
            column_split = factor(df_pall %>% filter(tumor_id%in%colnames(m)) %>% filter(!is.na(cdki) & !is.na(pall_1st_progression)) %>% arrange(desc(pall_1st_progression)) %>% pull(pall_1st_progression)), 
            column_gap=unit(0.3, "mm"),                                                                                                                                                 
            column_title_gp=gpar(fontsize=0),
            
            top_annotation = HeatmapAnnotation(
              tmb         = anno_barplot(df_tmb[match(colnames(m), df_tmb$tumor_id),] %>% dplyr::select(-tumor_id, -tmb) %>% as.matrix(), gp=gpar(col=NA, fill=c(pal_signal[[1]], pal_signal[[2]], pal_signal[[3]])), border = FALSE, bar_width=1),
              pfs         = anno_barplot(df_pall[match(colnames(m), df_pall$tumor_id),] %>% mutate(pfs=ifelse(is.na(pfs), max(pfs, na.rm=T), pfs)) %>% pull(pfs), gp=gpar(col=NA, fill="steelblue", border = FALSE, bar_width=1)),
              her2        = df_crf[match(colnames(m), df_crf$tumor_id),] %>% pull(HER2_diagnosis),
              ki67        = df_crf[match(colnames(m), df_crf$tumor_id),] %>% pull(ki67_diagnosis),
              pam50       = df_pam50[match(colnames(m), df_pam50$tumor_id),] %>% pull(pam50),
              hrd         = df_hrd[match(colnames(m), df_hrd$tumor_id),] %>% pull(prop),
              ploidy      = df_decision[match(colnames(m), df_decision$tumor_id),] %>% pull(ploidy),
              cdk46i      = df_pall[match(colnames(m), df_pall$tumor_id),] %>% pull(cdki),
              math        = df_math[match(colnames(m), df_math$tumor_id),] %>% pull(math),
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
                cdk46i      = c("abemaciclib"=unikn::uni_mannheim_2[[5]], "palbociclib"=unikn::uni_mannheim_2[[6]], "ribociclib"=unikn::uni_mannheim_2[[7]])
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
  ) 
  ,
  annotation_legend_list = list(
    Legend(labels=colnames(m_sig_sv), title="SV signature", type="grid", legend_gp = gpar(fill=uni_koeln_2[1:ncol(m_sig_sv)])),
    Legend(labels=colnames(m_sig_sbs), title="SBS signature", type="grid", legend_gp = gpar(fill=palette_26colors[1:ncol(m_sig_sbs)])),
    Legend(labels=colnames(m_sig_id), title="ID signature", type="grid", legend_gp = gpar(fill=unikn::uni_regensburg_3[1:ncol(m_sig_id)]))
  )
)

matrix(c(18,26,0,7), nrow=2) %>% fisher.test() # PTEN test
matrix(c(18,30,0,3), nrow=2) %>% fisher.test() # ESR1 test

# Figure 6D. Survival by TMB
df_pall %>% filter(!is.na(cdki) & !is.na(pall_1st_progression)) %>% left_join(df_tmb) %>% mutate(high_tmb=tmb>median(tmb)) %>% coxph(Surv(pfs, pall_1st_progression)~high_tmb, data=.) %>% summary()
df_pall %>% filter(!is.na(cdki) & !is.na(pall_1st_progression)) %>% left_join(df_tmb) %>% mutate(high_tmb=tmb>median(tmb)) %>% survfit2(Surv(pfs/365, pall_1st_progression)~high_tmb, data=.) %>% ggsurvfit() +  scale_color_jco() + theme_minimal() + theme(axis.line = element_line(), axis.ticks = element_line(), panel.grid=element_blank())

# Figure 6E. Survival by HRD
df_pall %>% filter(!is.na(cdki) & !is.na(pall_1st_progression)) %>% left_join(df_hrd) %>% mutate(high_hrd=prop>0.2) %>% coxph(Surv(pfs, pall_1st_progression)~high_hrd, data=.) %>% summary()
df_pall %>% filter(!is.na(cdki) & !is.na(pall_1st_progression)) %>% left_join(df_hrd) %>% mutate(high_hrd=prop>0.2) %>% survfit2(Surv(pfs/365, pall_1st_progression)~high_hrd, data=.) %>% ggsurvfit() +  scale_color_jco() + theme_minimal() + theme(axis.line = element_line(), axis.ticks = element_line(), panel.grid=element_blank())

# Figure 6F. multivariate analysis
df_temp <- df_pall %>% filter(!is.na(cdki) & !is.na(pall_1st_progression)) %>% left_join(df_hrd) %>% left_join(df_tmb) %>% 
  left_join(maf_drivers@data %>% as_tibble() %>% filter(Hugo_Symbol=="PTEN") %>% rename(tumor_id=Tumor_Sample_Barcode) %>% dplyr::select(tumor_id) %>% distinct() %>% mutate(bin=TRUE)) %>% 
  mutate(bin=ifelse(is.na(bin), FALSE, bin))
ggplot() + 
  geom_errorbar(aes(x=v, ymin=`2.5 %`, ymax=`97.5 %`), data=df_temp %>% mutate(hrd=prop>0.20) %>% coxph(Surv(pfs, pall_1st_progression)~hrd+tmb+bin+tmb:hrd, data=.) %>% confint() %>% exp() %>% as_tibble() %>% mutate(v=df_temp %>% mutate(hrd=prop>0.20) %>% coxph(Surv(pfs, pall_1st_progression)~hrd+tmb+bin+tmb:hrd, data=.) %>% confint() %>% exp() %>% rownames()) %>% mutate(v=factor(v, levels=rev(c("hrdTRUE", "tmb", "binTRUE", "hrdTRUE:tmb")))),
                width=0.1) + 
  geom_point(aes(x=v, y=value), data=df_temp %>% mutate(hrd=prop>0.20) %>% coxph(Surv(pfs, pall_1st_progression)~hrd+tmb+bin+tmb:hrd, data=.) %>% coef() %>% exp() %>% as_tibble() %>% mutate(v=df_temp %>% mutate(hrd=prop>0.20) %>% coxph(Surv(pfs, pall_1st_progression)~hrd+tmb+bin+tmb:hrd, data=.) %>% coef() %>% exp() %>% names())%>% mutate(v=factor(v, levels=rev(c("hrdTRUE", "tmb", "binTRUE", "hrdTRUE:tmb")))),
             size=5, shape=15) +
  geom_hline(yintercept = 1) + 
  coord_flip(ylim=c(0,20)) + 
  theme_minimal() 