# HRD score distribution

df_hrd %>% arrange(prop) %>% mutate(tumor_id=factor(tumor_id, levels=tumor_id)) %>% ggplot(aes(x=prop)) + geom_density(fill='grey70', color=NA) + theme_minimal() + 
  scale_y_continuous(expand=c(0,0)) + 
  scale_x_continuous(expand=c(0,0)) + 
  theme(panel.grid=element_blank(), axis.line.y=element_line(), axis.ticks.y=element_line()) + geom_vline(xintercept = 0.20)
