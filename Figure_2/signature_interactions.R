# Signature interactions

ggcorrplot(corr  = df_sig_snv_proportion %>% left_join(df_sig_ind_proportion) %>% left_join(df_sig_sv_proportion) %>% 
             dplyr::select_if(function(x) all(!is.na(x))) %>%dplyr::select(-tumor_id, -Unexplained) %>% 
             dplyr::select(v3_1, v3_5, ID1, ID2, v3_2, v3_13, v3_3, v3_8, ID6, ID8, RS3, RS5, everything()) %>% 
             as.matrix() %>% corr.test(adjust = "none") %>% .$r, 
           type  = "upper",
           insig = "pch",
           sig.lvl = c(0.05, 0.01, 0.001),
           show.diag=F
) + scale_fill_continuous_diverging(palette="Blue-Red", na.value=NA)
