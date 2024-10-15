# APOBEC germline deletion
# ------------------------

# Figure 2G. The distribution of germline APOBEC3A/3B deletion in our cohort, European breast cancer patients, and Korean normal population.
df_apobec %>% filter(!normal_control) %>% dplyr::count(gdel) %>% ggplot(aes(x=factor(1), y=n, fill=gdel)) + geom_col(position='fill',width=1,color='white') + coord_polar(theta='y') + geom_col(aes(x=0, y=0)) + scale_fill_nejm() + theme_void()
df_apobec %>% filter( normal_control) %>% dplyr::count(gdel) %>% ggplot(aes(x=factor(1), y=n, fill=gdel)) + geom_col(position='fill',width=1,color='white') + coord_polar(theta='y') + geom_col(aes(x=0, y=0)) + scale_fill_nejm() + theme_void()
tibble(gdel=factor(c("homo", "hetero", "wt"), levels=c("homo", "hetero","wt")), n=c(14, 128, 781)) %>% ggplot(aes(x=factor(1), y=n, fill=gdel)) + geom_col(position='fill',width=1,color='white') + coord_polar(theta='y') + geom_col(aes(x=0, y=0)) + scale_fill_nejm() + theme_void()   # from serena paper

# Figure 2H. TMB vs APOBEC3A/3B germline deletion 
plot_grid(
  df_context %>% left_join(df_apobec) %>% ggplot(aes(x=gdel, y=r_ytca_rtca)) + geom_jitter(color="grey80", width=0.2) + geom_boxplot(outlier.shape = NA) + theme_minimal() + theme(axis.line.y=element_line(), axis.ticks.y=element_line(), axis.text.x=element_blank(), axis.title.x=element_blank()) + geom_signif(comparisons=list(c("homo", "hetero"), c("hetero", "wt"))),
  df_tmb %>% left_join(df_apobec) %>% ggplot(aes(x=gdel, y=log(tmb), fill=gdel)) + 
    geom_violin(draw_quantiles = c(0.5)) + 
    geom_signif(comparisons=list(c("homo", "hetero"), c("hetero", "wt"))) + 
    theme_minimal() + theme(axis.line.y=element_line(), axis.ticks.y=element_line(), legend.position = 'none') + scale_fill_nejm(),
  ncol=1, align="v", rel_heights = c(0.4, 0.6))

# Figure 2I. The proportion of germline APOBEC3A/3B deletion carriers in patients with varying levels of APOBEC signatures.
bin_size=30
plot_grid(
  df_sig_snv_proportion %>% left_join(df_apobec) %>% arrange(desc(v3_2+v3_13)) %>%  mutate(BIN=c(rep(c(1:(nrow(.)%/%bin_size)), each=bin_size), rep((nrow(.)%/%bin_size)+1, each=nrow(.)%%bin_size))) %>% left_join(df_context) %>% ggplot(aes(x=BIN, y=r_ytca_rtca)) + geom_jitter() + scale_x_continuous(expand=c(0,0)) + theme_minimal() + theme(axis.text.x=element_blank(), axis.title.x=element_blank()),
  df_sig_snv_proportion %>% left_join(df_apobec) %>% arrange(desc(v3_2+v3_13)) %>%  mutate(BIN=c(rep(c(1:(nrow(.)%/%bin_size)), each=bin_size), rep((nrow(.)%/%bin_size)+1, each=nrow(.)%%bin_size))) %>% dplyr::count(BIN, gdel) %>% group_by(BIN) %>% mutate(n=n/sum(n), gdel=factor(gdel, levels=rev(c("homo", "hetero", "wt")))) %>% ggplot(aes(x=BIN, y=n*100, fill=gdel)) + geom_bar(stat='identity', width=1) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_minimal() + theme(legend.position = 'none', panel.grid=element_blank(), axis.line.y=element_line(), axis.ticks.y = element_line()) + ylab("Percentage (%)") + 
    scale_fill_manual(values=c("homo"=pal_nejm()(3)[1], "hetero"=pal_nejm()(3)[2], "wt"=pal_nejm()(3)[3])),
  ncol=1, align="v", rel_heights=c(0.4, 0.6))
