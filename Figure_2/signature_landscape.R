# Figure 2A. Mutational signature landscape

df_temp1 <- bind_rows(
  df_sig_snv_exposure %>% gather(sig, n, 2:ncol(.)) %>% filter(sig!="Unexplained") %>% filter(!is.na(n)) %>% filter(n>=1),
  df_sig_ind_exposure %>% gather(sig, n, 2:ncol(.)) %>% filter(!is.na(n)) %>% filter(n>=1),
  df_sig_sv_exposure %>% gather(sig, n, 2:ncol(.)) %>% filter(!is.na(n)) %>% filter(n>=1) 
) %>% mutate(sig=factor(sig, levels=levels(df_temp1$sig))) %>% 
  arrange(sig, n) %>% 
  group_by(sig) %>% mutate(g=cur_group_id()) %>% ungroup() %>% 
  mutate(sig_tid=str_c(sig, tumor_id)) %>% 
  mutate(sig_tid=factor(sig_tid, levels=sig_tid)) 

df_temp2 <- bind_rows(
  df_sig_snv_proportion %>% gather(sig, prop, 2:ncol(.)) %>% filter(sig!="Unexplained") %>% left_join(df_pam50) %>% filter(!is.na(pam50)) %>% mutate(bin=prop>10.0) %>% dplyr::count(sig, pam50, bin) %>% group_by(sig, pam50) %>% mutate(p=n/sum(n)) %>% ungroup() %>% filter(bin),
  df_sig_ind_proportion %>% gather(sig, prop, 2:ncol(.)) %>% left_join(df_pam50) %>% filter(!is.na(pam50)) %>% mutate(bin=prop>10.0) %>% dplyr::count(sig, pam50, bin) %>% group_by(sig, pam50) %>% mutate(p=n/sum(n)) %>% ungroup() %>% filter(bin),
  df_sig_sv_proportion  %>% gather(sig, prop, 2:ncol(.)) %>% left_join(df_pam50) %>% filter(!is.na(pam50)) %>% mutate(bin=prop>10.0) %>% dplyr::count(sig, pam50, bin) %>% group_by(sig, pam50) %>% mutate(p=n/sum(n)) %>% ungroup() %>% filter(bin),
) %>% mutate(sig=factor(sig, levels=levels(df_temp1$sig))) %>% complete(sig, pam50) %>% filter(pam50!="Normal")

plot_grid(
  df_temp1 %>% ggplot(aes(x=sig_tid, y=log10(n))) + geom_bar(aes(y=log10(max(n)), fill=factor(g%%2)), stat='identity', width=1, color=NA) + geom_point(size=0.3) + theme_minimal() + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position="none", axis.ticks.x=element_blank(), panel.grid=element_blank(), axis.line.y=element_line(), axis.ticks.y=element_line()) + 
    scale_y_continuous(expand=c(0,0)) + 
    scale_fill_manual(values=c("1"="grey90", "0"="white")) + facet_grid(cols=vars(sig), scales="free_x"),
  df_temp2 %>% ggplot(aes(x=sig, y=pam50, fill=p)) + geom_tile() + scale_fill_continuous_sequential(palette = "Purple-Yellow", limits=c(0,1)) + scale_x_discrete(expand=c(0.1,0)) + theme_minimal() +  theme(legend.position = "none", axis.text.x=element_text(angle=45, hjust=1.0, vjust=1.0), panel.grid=element_blank(), axis.title.x=element_blank()) + facet_grid(cols=vars(sig), scale="free_x"),
  ncol=1, align="v", rel_heights = c(0.5, 0.5)
)
