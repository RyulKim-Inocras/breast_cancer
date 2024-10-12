# Oncoprint 

library(maftools)
library(unikn)
library(circlize)

maf_drivers <- read.maf(maf="/home/users/chrono0707/analysis/08_breast/03_somatic/mafs/drivers/merged.tsv") 

palette_variant = c("Frame_Shift_Del"       = "blue", 
                    "Frame_Shift_Ins"        = "red", 
                    "In_Frame_Del"           = "yellow",
                    "In_Frame_Ins"           = "purple",
                    "Missense_Mutation"      = "#008000",
                    "Nonsense_Mutation"      = "brown",
                    "Nonstop_Mutation"       = "orange",
                    "Splice_Site"            = "cyan",
                    "Translation_Start_Site" = "green",
                    "DUP"                    = pal_jama()(4)[1],
                    "BND"                    = pal_jama()(4)[2],
                    "INV"                    = pal_jama()(4)[3],
                    "DEL"                    = pal_jama()(4)[4]
)

alter_fun = list(
  
  background = function(x,y,w,h) {
    grid.rect(x,y,w-unit(0,'pt'), h-unit(0,'pt'),
              gp=gpar(fill="#CCCCCC", col=NA))
  },
  Frame_Shift_Del = function(x,y,w,h) {
    grid.rect(x,y,(w-unit(0,'pt'))*0.9, (h-unit(0,'pt'))*0.9,
              gp=gpar(fill=palette_variant['Frame_Shift_Del'], col=NA))
  },
  Frame_Shift_Ins = function(x,y,w,h) {
    grid.rect(x,y,(w-unit(0,'pt'))*0.9, (h-unit(0,'pt'))*0.9,
              gp=gpar(fill=palette_variant['Frame_Shift_Ins'], col=NA))
  },
  In_Frame_Del = function(x,y,w,h) {
    grid.rect(x,y,w-unit(0,'pt'), h-unit(0,'pt'),
              gp=gpar(fill=palette_variant['In_Frame_Del'], col=NA))
  },
  In_Frame_Ins = function(x,y,w,h) {
    grid.rect(x,y,w-unit(0,'pt'), h-unit(0,'pt'),
              gp=gpar(fill=palette_variant['In_Frame_Ins'], col=NA))
  },
  Missense_Mutation = function(x,y,w,h) {
    grid.rect(x,y,w-unit(0,'pt'), h*0.33,
              gp=gpar(fill=palette_variant['Missense_Mutation'], col=NA))
  },
  Nonsense_Mutation = function(x,y,w,h) {
    grid.rect(x,y,w-unit(0,'pt'), h*0.33,
              gp=gpar(fill=palette_variant['Nonsense_Mutation'], col=NA))
  },
  Nonstop_Mutation = function(x,y,w,h) {
    grid.rect(x,y,w-unit(0,'pt'), h*0.33,
              gp=gpar(fill=palette_variant['Nonstop_Mutation'], col=NA))
  },
  Splice_Site =function(x,y,w,h) {
    grid.rect(x,y,w-unit(0,'pt'), h*0.33,
              gp=gpar(fill=palette_variant['Splice_Site'], col=NA))
  },
  Translation_Start_Site =function(x,y,w,h) {
    grid.rect(x,y,w-unit(0,'pt'), h*0.33,
              gp=gpar(fill=palette_variant['Translation_Start_Site'], col=NA))
  }
)

df_temp <- maf_drivers@data %>% as_tibble() %>% 
  filter(FILTER == "PASS") %>% 
  mutate(patient_id=str_sub(Tumor_Sample_Barcode, end=-6L)) %>% 
  dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Consequence) %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% 
  summarise(variant_class=paste0(Variant_Classification, collapse=',')) %>% 
  spread(Tumor_Sample_Barcode, variant_class) %>% 
  rowwise(Hugo_Symbol) %>% 
  mutate(count = sum(!is.na(c_across(everything())))) %>% 
  ungroup() 

m <- df_temp %>% 
  filter(count > ncol(.) * 0.01) %>% 
  dplyr::select(-Hugo_Symbol, -count) %>% 
  mutate_all(function(x) ifelse(is.na(x), "", x)) %>% 
  as.matrix()

rownames(m) <- df_temp %>% filter(count > ncol(.) * 0.01) %>% pull(Hugo_Symbol)

df_temp <- df_sig_snv_proportion %>% filter(tumor_id %in% colnames(m))
m_sig_sbs <- df_temp %>% dplyr::select(starts_with("v3")) %>% replace(is.na(.), 0) %>% as.matrix() 
rownames(m_sig_sbs) <- df_temp$tumor_id

df_temp <- df_sig_ind_proportion %>% filter(tumor_id %in% colnames(m)) 
m_sig_id <- df_temp %>% dplyr::select(starts_with("ID")) %>% replace(is.na(.), 0) %>% as.matrix()
rownames(m_sig_id) <- df_temp$tumor_id


df_temp <- df_sig_sv_proportion %>% filter(tumor_id %in% colnames(m)) 
m_sig_sv <- df_temp %>% dplyr::select(starts_with("RS")) %>% replace(is.na(.), 0) %>%  as.matrix() 
rownames(m_sig_sv) <- df_temp$tumor_id


# Load GISTIC results 
df_temp <- read_tsv("~/analysis/08_breast/99_R/output/gistic/all_lesions.conf_90.txt", show_col_types = F) %>% filter(`Amplitude Threshold` == "Actual Copy Change Given") %>% 
  rowwise() %>% mutate(region=str_split(`Region Limits`, "-", simplify = T)[1], CHROM=str_split(region, ":", simplify = T)[1], POS=str_split(region, ":", simplify = T)[2]) %>% ungroup() %>% 
  mutate(CHROM=factor(CHROM, levels=c(str_c("chr", 1:22))), POS=as.numeric(POS)) %>% 
  arrange(CHROM, POS)
m_cna <- df_temp %>% dplyr::select(starts_with("GINS")) %>% as.matrix()
m_cna <- m_cna[, colnames(m)]
rownames(m_cna) <- df_temp$`Wide Peak Limits`

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
            show_heatmap_legend=F,
            heatmap_legend_param = list(legend_direction="horizontal"),
            
            
            #column_split = factor(df_her2 %>% left_join(df_hrd) %>% filter(tumor_id %in% colnames(m)) %>% mutate(hrd=prop>0.7) %>% arrange(hrd, focal) %>% mutate(group=str_c(hrd, focal)) %>% pull(group)), # focal+/hrd+ vs focal-/hrd+ vs hrd- 비교 위해
            #column_gap=unit(0.3, "mm"),                                                                                                                                                                      # focal+/hrd+ vs focal-/hrd+ vs hrd- 비교 위해
            
            
            top_annotation = HeatmapAnnotation(
              tmb         = anno_barplot(df_tmb[match(colnames(m), df_tmb$tumor_id),] %>% dplyr::select(-tumor_id, -tmb) %>% as.matrix(), gp=gpar(col=NA, fill=c(pal_signal[[1]], pal_signal[[2]], pal_signal[[3]])), border = FALSE, bar_width=1),
              #menopause   = df_crf[match(colnames(m), df_crf$tumor_id),] %>% pull(menopausal_status),
              er          = df_crf[match(colnames(m), df_crf$tumor_id),] %>% pull(ER_diagnosis),
              pr          = df_crf[match(colnames(m), df_crf$tumor_id),] %>% pull(PR_diagnosis),
              her2        = df_crf[match(colnames(m), df_crf$tumor_id),] %>% pull(HER2_diagnosis),
              ki67        = df_crf[match(colnames(m), df_crf$tumor_id),] %>% pull(ki67_diagnosis),
              pam50       = df_pam50[match(colnames(m), df_pam50$tumor_id),] %>% pull(pam50),
              hrd         = df_hrd[match(colnames(m), df_hrd$tumor_id),] %>% pull(prop),
              ploidy      = df_decision[match(colnames(m), df_decision$tumor_id),] %>% pull(ploidy),
              math        = df_math[match(colnames(m), df_math$tumor_id),] %>% pull(math),
              #MSI     = df_msi[match(colnames(m), df_msi$tumor_id),] %>% pull(score),
              na_col="white",
              simple_anno_size=unit(0.15, "cm"),
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
                math        = colorRamp2(c(df_math[match(colnames(m), df_math$tumor_id),] %>% pull(math) %>% min(), df_math[match(colnames(m), df_math$tumor_id),] %>% pull(math) %>% max()), hcl_palette = "Reds", reverse = T)
              ),
              show_legend=F,
              annotation_legend_param = list(legend_direction="horizontal")
            ),           
            right_annotation = rowAnnotation(
              subtype = anno_barplot(as_tibble(m) %>% mutate(SYMBOL=rownames(m)) %>% select(SYMBOL, everything()) %>% gather(tumor_id, s, 2:ncol(.)) %>% filter(s!="") %>% left_join(df_pam50) %>% filter(!is.na(pam50)) %>% count(SYMBOL, pam50) %>% group_by(SYMBOL) %>% mutate(p=n/sum(n)*100) %>% ungroup() %>% select(SYMBOL, pam50, p) %>% spread(pam50, p) %>% select(-SYMBOL) %>% replace(is.na(.), 0) %>% as.matrix(),
                                     gp=gpar(col=NA, fill=c("Normal"="gray90", "LumA"=ggsci::pal_nejm()(5)[4], "LumB"=ggsci::pal_nejm()(5)[3], "Her2"=ggsci::pal_nejm()(5)[1], "Basal"=ggsci::pal_nejm()(5)[5])), border = FALSE, bar_width=1),
              rbar = anno_oncoprint_barplot(axis_param = list(side="top", labels_rot=0, gp=gpar(cex=0.5)))),
            bottom_annotation = HeatmapAnnotation(
              sbs = anno_barplot(m_sig_sbs, gp = gpar(fill=palette_26colors[1:ncol(m_sig_sbs)], col=NA), bar_width=1, height=unit(0.7, "cm"), border=FALSE, axis_param = list(gp=gpar(cex=0.5))),
              id  = anno_barplot(m_sig_id,  gp = gpar(fill=unikn::uni_regensburg_3[1:ncol(m_sig_id)], col=NA), bar_width=1, height=unit(0.7, "cm"), border=FALSE, axis_param = list(gp=gpar(cex=0.5))),
              sv  = anno_barplot(m_sig_sv,  gp = gpar(fill=unikn::uni_koeln_2[1:ncol(m_sig_sv)], col=NA), bar_width=1, height=unit(0.7, "cm"), border=FALSE, axis_param = list(gp=gpar(cex=0.5))),
              annotation_name_gp = gpar(cex=0.5),
              annotation_name_side = "left"
            ),
            show_column_names = FALSE,
            row_names_gp = gpar(cex=0.4),
            pct_gp = gpar(cex=0.5)
            
  ), %v%
    Heatmap(
      m_cna, cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = F, height=unit(7, "cm"),
      row_split = factor(str_extract(rownames(m_cna), "chr\\d+"), levels=str_c("chr", 1:22)),
      row_gap = unit(0.3, "mm"),
      row_title_gp=gpar(fontsize=7),
      row_title_rot=0,
      show_heatmap_legend=FALSE
    ),
  annotation_legend_list = list(
    Legend(labels=colnames(m_sig_sv), title="SV signature", type="grid", legend_gp = gpar(fill=uni_koeln_2[1:ncol(m_sig_sv)]), nrow=1),
    Legend(labels=colnames(m_sig_sbs), title="SBS signature", type="grid", legend_gp = gpar(fill=palette_26colors[1:ncol(m_sig_sbs)]), nrow=1),
    Legend(labels=colnames(m_sig_id), title="ID signature", type="grid", legend_gp = gpar(fill=unikn::uni_regensburg_3[1:ncol(m_sig_id)]), nrow=1)
  ),
  gap=unit(0,"mm")
)
