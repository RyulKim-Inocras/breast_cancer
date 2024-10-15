# Timing of copy number amplifications (CNAs) in breast cancer
# ------------------------------------------------------------

# Calculate amplification timing
df_gistic <- read_tsv("output/gistic/all_lesions.conf_90.txt", show_col_types = F) %>% filter(`Amplitude Threshold`=="0: t<0.1; 1: 0.1<t< 0.9; 2: t>0.9") %>% filter(str_detect(`Unique Name`, "Amplification"))

df_gistic <- df_gistic %>% mutate(total=rowSums(df_gistic %>% dplyr::select(starts_with("GINS")) %>% mutate_all(function(x) ifelse(x==2, TRUE, FALSE)))) %>%
  filter(`q values`<0.10) %>% 
  rowwise() %>% mutate(
    `Region Limits` = str_remove(`Region Limits`, "\\(.*\\)"),
    `Wide Peak Limits` = str_remove(`Wide Peak Limits`, "\\(.*\\)"),
    chrom=str_split(`Wide Peak Limits`, ':', simplify = T)[1], 
    start=as.numeric(str_split(str_split(`Wide Peak Limits`, ':', simplify = T)[2], "-", simplify = T)[1]),
    end  =as.numeric(str_split(str_split(`Wide Peak Limits`, ":", simplify = T)[2], "-", simplify = T)[2]),
    gap  =end-start) %>% ungroup() %>% 
  dplyr::select(total, chrom, start, end, gap, starts_with("GINS")) %>% 
  arrange(chrom, start, end, desc(total)) %>% 
  group_by(chrom, start, end) %>% dplyr::slice(1) %>% ungroup() %>% 
  arrange(desc(total))

df_amp <- tibble()
for(i in 1:nrow(df_gistic)) {
  print(str_c(df_gistic$chrom[i], df_gistic$start[i], df_gistic$end[i], sep=" "))
  df_amp <- bind_rows(df_amp,
                      df_gistic %>% dplyr::slice(i) %>% dplyr::select(starts_with("GINS")) %>% gather(tumor_id) %>% filter(value==2) %>% transpose() %>% map(function(x) {
                        print(x$tumor_id)
                        
                        # point mutation
                        df_pm <- bind_rows(
                          read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/03_somatic/", x$tumor_id, "_", str_sub(x$tumor_id, end=14), "-10AD.SNVs.quick_fi.tsv.seqzcn.ccf"),   show_col_types = F, guess_max=1e6) %>% filter(`#CHROM` %in% str_c("chr", 1:22)) %>% mutate(mutCN=as.numeric(mutCN)) %>% dplyr::select(`#CHROM`, POS, REF, ALT, vaf, mutCN, x$tumor_id) %>% dplyr::rename(format=x$tumor_id),
                          read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/03_somatic/", x$tumor_id, "_", str_sub(x$tumor_id, end=14), "-10AD.indels.quick_fi.tsv.seqzcn.ccf"), show_col_types = F, guess_max=1e6) %>% filter(`#CHROM` %in% str_c("chr", 1:22)) %>% mutate(mutCN=as.numeric(mutCN)) %>% dplyr::select(`#CHROM`, POS, REF, ALT, vaf, mutCN, x$tumor_id) %>% dplyr::rename(format=x$tumor_id)
                          ) %>% 
                          filter(`#CHROM`==df_gistic$chrom[i] & POS > df_gistic$start[i] & POS < df_gistic$end[i]) %>% 
                          mutate(chrom=df_gistic$chrom[i], start=df_gistic$start[i], end=df_gistic$end[i], tumor_id=x$tumor_id) %>%
                          dplyr::select(chrom, start, end, tumor_id, everything())
                        
                        df_seqz    <- read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/06_seqz/", x$tumor_id, "_", str_sub(x$tumor_id, end=14), "-10AD.segments.clean.txt"), show_col_types = F) 
                        gr_seqz    <- GRanges(seqnames = df_seqz$`#CHROM`, ranges=IRanges(start=df_seqz$start_pos, end=df_seqz$end_pos), majCN=df_seqz$majCN, minCN=df_seqz$minCN)
                        gr_segment <- GRanges(seqnames = df_gistic$chrom[i], ranges=IRanges(start=df_gistic$start[i], end=df_gistic$end[i]))
                        
                        seqz_idx   <- subjectHits(findOverlaps(gr_segment, gr_seqz))
                        
                        df_pm$maxCN <- max(gr_seqz[seqz_idx]$majCN) 
                        df_pm$minCN <- min(gr_seqz[seqz_idx]$minCN)  
                        df_pm$mean_majCN <- mean(gr_seqz[seqz_idx]$majCN)
                        df_pm$mean_minCN <- mean(gr_seqz[seqz_idx]$minCN)
    
      return(df_pm %>% dplyr::select(chrom, start, end, tumor_id, maxCN, minCN, mean_majCN, mean_minCN, everything()))
    }) %>% do.call(bind_rows, .))
}


# Draw figures
# ------------

# Figure 5A. Relative timing of recurrent CNAs.

df_temp1 <- df_amp %>% 
  left_join(read_tsv("output/gistic/all_lesions.conf_90.txt", show_col_types = F) %>% filter(`Amplitude Threshold`=="0: t<0.1; 1: 0.1<t< 0.9; 2: t>0.9") %>% filter(str_detect(`Unique Name`, "Amplification")) %>% rowwise() %>% mutate(`Wide Peak Limits` = str_remove(`Wide Peak Limits`, "\\(.*\\)"), chrom=str_split(`Wide Peak Limits`, ':', simplify = T)[1], start=as.numeric(str_split(str_split(`Wide Peak Limits`, ':', simplify = T)[2], "-", simplify = T)[1]), end  =as.numeric(str_split(str_split(`Wide Peak Limits`, ":", simplify = T)[2], "-", simplify = T)[2])) %>% dplyr::select(Descriptor, chrom, start, end)) %>% 
  dplyr::slice(-queryHits(findOverlaps(GRanges(seqnames = chrom, IRanges(start=start, end=end)), resize(gr_centromere, width = width(gr_centromere)+(1e6*2), fix = "center"), type="within"))) %>% # Remove regions near the centromore
  group_by(Descriptor, chrom, start, end, tumor_id) %>%
  mutate(cutoff=mean(mean_majCN)*0.75, bin=mutCN>cutoff) %>%  
  summarise(pre=sum(bin), post=n()-pre, meanCN=mean(mean_majCN)) %>% ungroup() %>%
  left_join(df_sig_snv_proportion %>% dplyr::select(tumor_id, v3_1, v3_5)) %>% 
  left_join(df_pam50) %>% 
  mutate(pre_clocklike=pre*(v3_1)/100, post_clocklike=post*(v3_1)/100, ratio=pre_clocklike/(pre_clocklike+post_clocklike)) %>% 
  mutate(chrom=factor(chrom, levels=str_c("chr", c(1:22, "X", "Y")))) %>% arrange(chrom, desc(start)) %>% 
  group_by(Descriptor) %>% mutate(n=n()) %>% ungroup() %>% filter(n>30) %>% # Include only locations with a sufficient number of mutations
  group_by(Descriptor) %>% filter(ratio>quantile(ratio, 0.10) & ratio<quantile(ratio, 0.90)) %>% ungroup() %>%  # remove extreme outliers
  mutate(Descriptor=factor(Descriptor, levels=unique(Descriptor))) 
    
df_temp2 <- df_temp1 %>% 
  mutate(dsc=Descriptor) %>% 
  nest(data=c(-Descriptor, -chrom)) %>% 
  mutate(
    results = sapply(data, function(d) {
      os   <- summary(coxph(Surv(os, os_event)~bin, data=df_gistic %>% filter(Descriptor==d$dsc[1]) %>% gather(tumor_id, bin, 7:ncol(.)) %>% group_by(tumor_id) %>% summarise(bin=any(bin==2)) %>% left_join(df_crf)))
      d <- d %>% filter(ratio>quantile(ratio, 0.10) & ratio<quantile(ratio, 0.90))
      return(
        list(c(
          'n'         = nrow(d), 
          'mean'      = mean(d$ratio),
          'sd'        = 1.96*sd(d$ratio)/sqrt(nrow(d)),
          'meanCN'    = mean(d$meanCN),
          'sdCN'      = 1.96*sd(d$meanCN)/sqrt(nrow(d)),
          'os_hr'     = os$coef[2],
          'os_ll'     = os$conf.int[, "lower .95"],
          'os_ul'     = os$conf.int[, "upper .95"],
          'os_p'      = os$wald["pvalue"],
          'luma'      = sum(d$pam50=="LumA", na.rm = T),
          'lumb'      = sum(d$pam50=="LumB", na.rm = T),
          'her2'      = sum(d$pam50=="Her2", na.rm = T),
          'basal'     = sum(d$pam50=="Basal", na.rm = T),
          'normal'    = sum(d$pam50=="Normal", na.rm = T),
          'unknown'   = sum(is.na(d$pam50))
        )
        )
      )
    })) %>% unnest_wider(results) 

plot_grid(
  df_temp2 %>% ggplot(aes(x=Descriptor, y=n)) + geom_bar(stat='identity', width=1, linewidth=0.1, fill="grey50", color="white") + scale_x_discrete(labels=cytoband) + coord_flip() + ylab("Prevalence") + scale_y_reverse() + theme_minimal() + theme(axis.title.y=element_blank(), axis.text.y=element_text(size=5), panel.grid=element_blank(), axis.line.x=element_line(), axis.ticks.x=element_line(), strip.background = element_blank(), strip.text.y = element_blank()) + facet_grid(rows=vars(chrom), scales="free_y", space = "free_y"),
  df_temp2 %>% dplyr::select(Descriptor, chrom, unknown, normal, luma, lumb, her2, basal) %>% gather(pam50, n, 3:ncol(.)) %>% group_by(Descriptor, chrom) %>% mutate(p=n/sum(n)) %>% ungroup() %>% 
    ggplot(aes(x=Descriptor, y=p, fill=pam50)) + geom_bar(stat='identity', width=1, color="white", linewidth=0.1) + coord_flip() + scale_fill_manual(values=c("unknown"="gray90", "normal"="gray50", "luma"=ggsci::pal_nejm()(5)[4], "lumb"=ggsci::pal_nejm()(5)[3], "her2"=ggsci::pal_nejm()(5)[1], "basal"=ggsci::pal_nejm()(5)[5])) + theme_minimal() + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), legend.position = "none", plot.margin = margin(0,0,0,0, "cm"), panel.grid=element_blank(), axis.line.x=element_line(), axis.ticks.x=element_line(), strip.background = element_blank(), strip.text.y = element_blank()) + facet_grid(rows=vars(chrom), scales="free_y", space="free_y"),
  df_temp1 %>% ggplot(aes(x=Descriptor, y=ratio)) + geom_violin(width=1, linewidth=0.2, fill="grey99") + stat_summary(fun = "mean", geom = "point", color = "black", size=0.5) + coord_flip(ylim=c(0,1)) + theme_minimal() + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_line(), strip.background = element_blank(), strip.text.y = element_blank(), panel.background = element_rect(color="black", fill=NA)) + facet_grid(rows=vars(chrom), scales="free_y", space="free_y"), 
  df_amp %>% 
    left_join(read_tsv("output/gistic/all_lesions.conf_90.txt", show_col_types = F) %>% filter(`Amplitude Threshold`=="0: t<0.1; 1: 0.1<t< 0.9; 2: t>0.9") %>% filter(str_detect(`Unique Name`, "Amplification")) %>% rowwise() %>% mutate(`Wide Peak Limits` = str_remove(`Wide Peak Limits`, "\\(.*\\)"), chrom=str_split(`Wide Peak Limits`, ':', simplify = T)[1], start=as.numeric(str_split(str_split(`Wide Peak Limits`, ':', simplify = T)[2], "-", simplify = T)[1]), end  =as.numeric(str_split(str_split(`Wide Peak Limits`, ":", simplify = T)[2], "-", simplify = T)[2])) %>% dplyr::select(Descriptor, chrom, start, end)) %>% filter(Descriptor%in%df_temp1$Descriptor) %>% mutate(Descriptor=factor(Descriptor, levels=levels(df_temp1$Descriptor)), chrom=factor(chrom, levels=str_c("chr", c(1:22, "X", "Y")))) %>% ggplot(aes(x=Descriptor, y=mean_majCN)) + geom_boxplot(outlier.colour = NA, linewidth=0.2, width=0.51) + coord_flip() + theme_minimal() + theme(axis.text.y=element_blank(), panel.grid=element_line(linewidth=0.1), axis.title.y=element_blank(), axis.line.x=element_line(), axis.ticks.x=element_line(), strip.background = element_blank(), strip.text.y = element_blank()) + facet_grid(rows=vars(chrom), space="free_y", scales="free_y"),
  df_temp2 %>% ggplot(aes(x=Descriptor, y=os_hr, ymin=os_ll, ymax=os_ul)) + geom_pointrange(size=0.3, linewidth=0.1, shape=15) + coord_flip() + geom_hline(yintercept=1, linewidth=0.1) + theme_minimal() + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), panel.grid=element_blank(), axis.line.x=element_line(), axis.ticks.x=element_line(), strip.background = element_blank(), strip.text.y = element_blank()) + facet_grid(rows=vars(chrom), scales="free_y", space="free_y"),
  align = "h", nrow=1, rel_widths = c(0.2, 0.05, 0.3, 0.1, 0.1)
)


# Figure 5B. Kaplan-Meier survival curves for patients with and without 9p23 amplification
df_gistic %>% filter(Descriptor=="9p23") %>% gather(tumor_id, bin, 7:ncol(.)) %>% group_by(tumor_id) %>% summarise(bin=any(bin==2)) %>% left_join(df_crf) %>% left_join(df_pam50) %>% filter(pam50=="Basal") %>% mutate(os_event=ifelse(os_event==1 & os>4000, 0, os_event), os=ifelse(os>4000, 4000, os)) %>% survfit2(Surv(os/365, os_event)~bin, data=.) %>% ggsurvfit() + add_risktable() + scale_color_jco() + theme_minimal() + theme(axis.line = element_line(), axis.ticks = element_line(), panel.grid=element_blank())
df_gistic %>% filter(Descriptor=="9p23") %>% gather(tumor_id, bin, 7:ncol(.)) %>% group_by(tumor_id) %>% summarise(bin=any(bin==2)) %>% left_join(df_crf) %>% left_join(df_pam50) %>% filter(pam50=="Basal") %>% coxph(Surv(os, os_event)~bin+TNM, data=.) %>% summary()


# Figure 5C. The duration of CNAs in patients with homologous recombination DNA repair deficiency (HRD) versus proficiency (HRP).
df_centromere <- read_tsv("hg38_centromere_edit.tsv")
gr_centromere <- GRanges(df_centromere$chromosome, IRanges(start=df_centromere$exclude_start, end=df_centromere$exclude_end))

for(tid in df_crf$tumor_id) {
    df_temp <- read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/06_seqz/", tid, "_", str_sub(tid, end=14), "-10AD.segments.clean.txt"), show_col_types = F) %>% 
        filter(majCN>minCN) %>% 
        dplyr::slice(-queryHits(findOverlaps(GRanges(seqnames = `#CHROM`, IRanges(start=start_pos, end=end_pos)), resize(gr_centromere, width = width(gr_centromere)+(1e6*2), fix = "center"), type="within"))) 

    if(nrow(df_temp)==0) next

    df_temp %>%
        transpose() %>% map(function(x) {
        return(bind_rows(
        read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/03_somatic/", tid, "_", str_sub(tid, end=14), "-10AD.SNVs.quick_fi.tsv.seqzcn.ccf"),   show_col_types = F, guess_max=1e6) %>% filter(`#CHROM` %in% str_c("chr", c(1:22, "X", "Y"))) %>% mutate(mutCN=as.numeric(mutCN)) %>% dplyr::select(`#CHROM`, POS, REF, ALT, vaf, mutCN, tid) %>% dplyr::rename(format=tid),
        read_tsv(str_c("/home/users/chrono0707/analysis/08_breast/03_somatic/", tid, "_", str_sub(tid, end=14), "-10AD.indels.quick_fi.tsv.seqzcn.ccf"), show_col_types = F, guess_max=1e6) %>% filter(`#CHROM` %in% str_c("chr", c(1:22, "X", "Y"))) %>% mutate(mutCN=as.numeric(mutCN)) %>% dplyr::select(`#CHROM`, POS, REF, ALT, vaf, mutCN, tid) %>% dplyr::rename(format=tid)) %>% 
        filter(`#CHROM`==x[1] & POS > as.numeric(x[2]) & POS < as.numeric(x[3])) %>%         
        mutate(chrom=x[1], start=as.numeric(x[2]), end=as.numeric(x[3]), majCN=as.numeric(x[5]), minCN=as.numeric(x[6])) %>%
        dplyr::select(chrom, start, end, majCN, minCN, everything()))
    }) %>% do.call(bind_rows, .) %>% 
        mutate(pre=mutCN>majCN*0.75) %>% 
        group_by(chrom, start, end, majCN, minCN) %>% summarise(n=n(), pre=sum(pre), post=n()-sum(pre)) %>% ungroup() %>% 
        mutate(pre_clock=pre*df_sig_snv_proportion$v3_1[df_sig_snv_proportion$tumor_id==tid]/100, post_clock=post*df_sig_snv_proportion$v3_1[df_sig_snv_proportion$tumor_id==tid]/100, ratio=pre_clock/(pre_clock+post_clock)) %>% 
        mutate(chrom=factor(chrom, levels=c(str_c("chr", c(1:22, "X", "Y"))))) %>% 
        arrange(chrom, start) %>% 
        mutate(segment=str_c(chrom, ":", start, "-", end)) %>% 
        mutate(segment=factor(segment, levels=segment)) %>% write_tsv(str_c("/home/users/chrono0707/analysis/08_breast/99_R/output/amp_timing/", tid, ".tsv"))
}

df_trajectory <- tibble()
for(f in list.files("output/amp_timing/", pattern="*.tsv", full.names = T)) {
  print(f)
  df_trajectory <- bind_rows(
    df_trajectory,
    read_tsv(f, show_col_types = F) %>% mutate(tumor_id=str_sub(basename(f), end=19)) %>% dplyr::select(tumor_id, everything())
  )
}
df_trajectory <- df_trajectory %>% filter(tumor_id %in% df_crf$tumor_id)

df_trajectory %>% filter(majCN>=2) %>% filter(ratio!=0 & ratio!=1) %>% group_by(tumor_id) %>% summarise(gap=max(ratio)-min(ratio)) %>% left_join(df_hrd) %>% filter(!is.na(prop)) %>% mutate(bin=prop>0.2) %>% 
  ggplot(aes(x=gap, fill=bin)) + geom_density(alpha=0.5, linewidth=0.3) + scale_fill_manual(values=c("FALSE"="gray90", "TRUE"=rgb(44/255, 110/255, 179/255))) + theme_minimal() + theme(panel.grid=element_blank(), axis.line.y=element_line(), axis.ticks.y=element_line())

# Figure 5D. Circos plots and CNA trajectories for a triple negative breast cancer patient with HRD, and for a hormone receptor-positive breast cancer patient with HRD carrying a germline BRCA1 frameshift deletion.
read_tsv("output/amp_timing/GINS-0025-0024-01AD.tsv") %>% filter(n>10) %>% filter(majCN>=2) %>% 
  arrange(ratio) %>% mutate(cumratio=cumsum(ratio)) %>% select(ratio, cumratio) %>% bind_rows(tibble(ratio=1, cumratio=max(.$cumratio))) %>% ggplot(aes(x=ratio, y=cumratio, group=1)) + geom_line() + geom_point(shape=15) + ylim(0,50) + theme_minimal() 

read_tsv("output/amp_timing/GINS-0025-1360-01AD.tsv") %>% filter(n>10) %>% filter(majCN>=2) %>% 
  arrange(ratio) %>% mutate(cumratio=cumsum(ratio)) %>% select(ratio, cumratio) %>% bind_rows(tibble(ratio=1, cumratio=max(.$cumratio))) %>% ggplot(aes(x=ratio, y=cumratio, group=1)) + geom_line() + geom_point(shape=15) + ylim(0,50) + theme_minimal()

