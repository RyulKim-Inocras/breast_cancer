# HRD for adjuvant chemotherapy in basal-like breast cancer patients 

df_adj <- df_crf %>% 
  mutate(
    date_diagnosis = str_replace_all(date_diagnosis, "-", "."), 
    date_recur     = str_replace_all(date_recur, "-", "."),
    date_adj_start = str_replace_all(adj_start, "-", "."),
    date_recur     = ifelse(recur==0, "2023.11.1", date_recur), # data cutoff as of 2023.11.1.
    recur          = factor(recur, levels=c(0,1))
  ) %>% 
  mutate(
    date_diagnosis = as.POSIXct(date_diagnosis, format="%Y.%m.%d", tz="UTC"), 
    date_recur=as.POSIXct(date_recur, "%Y.%m.%d", tz="UTC"), 
    date_adj_start=as.POSIXct(date_adj_start, "%Y.%m.%d", tz="UTC"),
    rfs = as.numeric(difftime(date_recur, date_adj_start, units = "days"))) 

# KM curve
df_adj %>% 
  filter(tumor_id %in% (df_pam50 %>% filter(pam50=="Basal") %>% pull(tumor_id))) %>% 
  filter(adj==1) %>% 
  filter(adj_regimen %in% c(1,2,3,4,5)) %>% # adj regimens containing anthracycline plus cyclophosphamide 
  left_join(df_hrd) %>% 
  mutate(bin=prop>0.2) %>% 
  mutate(recur=ifelse(rfs>3000 & recur==1, "0", as.character(recur)), recur=as.integer(recur), rfs=ifelse(rfs>3000, 3000, rfs)) %>% 
  survfit2(Surv(rfs/365, recur)~bin, data=.) %>% ggsurvfit() + scale_color_jco() + theme_minimal() + theme(axis.line = element_line(), axis.ticks = element_line(), panel.grid=element_blank(), legend.position='bottom')

# Multivariate cox regression analysis 
df_temp <- 
  df_adj %>% filter(tumor_id %in% (df_pam50 %>% filter(pam50=="Basal") %>% pull(tumor_id))) %>% 
  filter(adj==1) %>% 
  filter(adj_regimen %in% c(1,2,3,4,5)) %>% 
  left_join(df_hrd) %>% mutate(bin=prop>0.2) %>%  
  mutate(recur=ifelse(rfs>3000 & recur==1, "0", as.character(recur)), recur=as.integer(recur), rfs=ifelse(rfs>3000, 3000, rfs)) %>% 
  mutate(stage_bin=ifelse(TNM %in% c("IA", "IIA"), TRUE, FALSE)) 

ggplot() + 
  geom_errorbar(aes(x=v, ymin=`2.5 %`, ymax=`97.5 %`), data=df_temp %>% coxph(Surv(rfs, recur)~age+stage_bin+bin, data=.) %>% confint() %>% exp() %>% as_tibble() %>% mutate(v=df_temp %>% coxph(Surv(rfs, recur)~age+stage_bin+bin, data=.) %>% confint() %>% exp() %>% rownames()) %>% mutate(v=factor(v, levels=rev(c("age", "stage_binTRUE", "binTRUE")))),
                width=0.1) + 
  geom_point(aes(x=v, y=value), data=df_temp %>% coxph(Surv(rfs, recur)~age+stage_bin+bin, data=.) %>% coef() %>% exp() %>% as_tibble() %>% mutate(v=df_temp %>% coxph(Surv(rfs, recur)~age+stage_bin+bin, data=.) %>% coef() %>% exp() %>% names())%>% mutate(v=factor(v, levels=rev(c("age", "stage_binTRUE", "binTRUE")))),
             size=5, shape=15) +
  geom_hline(yintercept = 1) + 
  coord_flip() + 
  theme_minimal() 

