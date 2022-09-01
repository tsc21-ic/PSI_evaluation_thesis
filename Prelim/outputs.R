pacman::p_load(selectiveInference, prospectr, knockoff, mvtnorm, glmnet, stabs, MASS, dplyr,ggthemes, kableExtra, ggplot2, tibble, tidyverse,
               viridis, pbapply, data.table, progress, randomcoloR, cowplot, ggpubr,tcltk, ggh4x, patchwork, profvis, fGarch, mvtnorm, beepr)
source("fun.R")
require(R.utils)

############ Load data ##############
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName, verbose=TRUE)
  get(ls()[ls() != "fileName"])
}

my_results=loadRData('0901_length10502.Rdata')
my_results_DS_UV=loadRData('0901_DS_R_10502.Rdata')

cat_lists <- function(list1, list2) {  
  keys <- unique(c(names(list1), names(list2)))
  map2(list1[keys], list2[keys], c) %>% 
    set_names(keys)
}
my_results = cat_lists(my_results, my_results_DS_UV)
my_results2=data.frame(
  beta=my_results$beta_poly,
  lee=my_results$length_poly_theory,
  mc=my_results$length_poly_ori,
  mc_r=my_results$length_poly_ori_R ) %>%
  group_by(beta)%>%
  pivot_longer(c(lee, mc, mc_r))

my_results_DS=data.frame(beta=my_results$beta_DS,DS=my_results$LENGTH_DS)%>%
  group_by(beta)%>%
  pivot_longer(DS)

my_results_UV=data.frame(beta=my_results$beta_R,UV=my_results$LENGTH_R)%>%
  group_by(beta)%>%
  pivot_longer(UV)

my_results2=rbind(my_results2, my_results_DS, my_results_UV)

method_names=c(expression(paste(Lee['M,s'])),
               expression(paste(MC['M'])),
               expression(paste(MC['M,r+'])),
               expression(DS),
               expression(('U, V')))
beta_names = c(
  expression(paste(beta[i],' = 1')),
  expression(paste(beta[i],' = 0.8')),
  expression(paste(beta[i],' = 0.2')))

my_results2=my_results2 %>% 
  mutate(beta2= factor(beta, levels=c(1,0.8,0.2), labels = beta_names),
         my_methods=factor(name,levels=c("lee","mc","mc_r", 'DS', '(U, V)'),labels=method_names))

my_methods = c("Lee_M,s", "MC_Lee_M", "MC_Lee_M+R", 'DS', '(U, V)')
length_final=ggplot(data = my_results2, mapping = aes(x=my_methods, y=value))+
  geom_boxplot(position=position_dodge(0.7),linetype="dashed")+
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill=my_methods),outlier.shape = 1) +
  stat_boxplot(geom = "errorbar" , aes(ymin= ..ymax..),width=0.2)+
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width=0.2) +
  facet_grid(~beta2, labeller = label_parsed, scales = "free_y")+
  theme_bw()+
  scale_fill_discrete("Method", labels= method_names)+
  scale_x_discrete("",labels=method_names)+
  coord_cartesian(ylim=c(0.2, 2)) +
  ylab("")+
  guides(fill=guide_legend(title="Method"))+
  theme(legend.position = "None",
        strip.background = element_blank(),
        strip.text.x = element_text(size=15),
        axis.text.x = element_text(size=15))
length_final
# ggsave(length_final, file='length_final.pdf', width = 210, height = 297, units = "mm")

######################## Prepare table  #######################

tab_length = data.frame(beta=factor(my_results$beta_poly,levels=c("1","0.8", "0.2")),
                        leeMs=my_results$length_poly_theory,
                        MC_M = my_results$length_poly_ori,
                        MC_MR = my_results$length_poly_ori_R)

tab_length = tab_length %>% 
  dplyr::group_by(beta) %>% 
  summarize(across(everything(), ~format(round(median(.), 3), nsmall=3)))
tab_length = tab_length %>% mutate(layer="Med. length")


tab_length_DS=data.frame(beta=my_results$beta_DS, DS = my_results$LENGTH_DS) %>%
  group_by(beta) %>%
  summarize(across(everything(), ~format(round(median(.), 3), nsmall=3))) %>%
  select(-beta)


tab_length_UV=data.frame(beta=my_results$beta_R, UV = my_results$LENGTH_R) %>%
  group_by(beta) %>%
  summarize(across(everything(), ~format(round(median(.), 3), nsmall=3))) %>%
  select(-beta)

tab_length=cbind(tab_length,tab_length_DS, tab_length_UV)
############

tab_data = data.frame(beta=factor(my_results$beta_poly,levels=c("1","0.8","0.2")),
                      leeMs=my_results$coverage_poly_theory,
                      MC_M = my_results$coverage_poly_ori,
                      MC_MR = my_results$coverage_poly_ori_R)

tab_data = tab_data %>% 
  dplyr::group_by(beta) %>% 
  summarize(across(everything(), ~format(round(mean(.)*100, 2), nsmall=2)))
tab_data = tab_data %>% mutate(layer="Coverage (%)")

tab_cov_DS=data.frame(beta=my_results$beta_DS, DS = my_results$COV_DS_HD) %>%
  dplyr::group_by(beta) %>% 
  summarize(across(everything(), ~format(round(mean(.)*100, 2), nsmall=2))) %>%
  select(-beta)


tab_cov_UV=data.frame(beta=my_results$beta_R, UV = my_results$COV_R_HD) %>%
  dplyr::group_by(beta) %>% 
  summarize(across(everything(), ~format(round(mean(.)*100, 2), nsmall=2))) %>%
  select(-beta)
tab_cov=cbind(tab_data,tab_cov_DS, tab_cov_UV)

######

tab_data_prop = data.frame(beta=factor(my_results$beta_poly,levels=c("1","0.8", "0.2")), 
                           leeMs=my_results$length_poly_theory,
                           MC_M = my_results$length_poly_ori,
                           MC_MR = my_results$length_poly_ori_R)
# MC_MR_more = my_results$length_poly_ori_R_more)
tab_data_prop = tab_data_prop %>% split(.$beta)
# filter(pmap_lgl(across(leeMs), ~is.infinite(...)))

tab_data_prop = lapply(tab_data_prop, function(x){
  sapply(x[, 2:ncol(x)], function(y) format(round(mean(is.infinite(y))*100, 2), nsmall=2))
})

tab_data_prop = data.frame(beta=unique(factor(my_results$beta_poly,levels=c("1","0.8", "0.2"))), 
                           bind_rows(tab_data_prop),
                           layer="Inf. length (%)")


tab_prop_DS= data.frame(beta=factor(my_results$beta_DS,levels=c("1","0.8", "0.2")), DS = my_results$LENGTH_DS) %>%
  group_by(beta) %>%
  summarize(DS=format(round(mean(DS > 1e6), 2),nsmall=2)) %>%
  select(-beta)



tab_prop_UV=data.frame(beta=my_results$beta_R, UV = my_results$LENGTH_R) %>%
  group_by(beta) %>%
  summarize(UV=format(round(mean(UV > 1e6), 2),nsmall=2)) %>%
  select(-beta)

tab_data_prop=cbind(tab_data_prop,tab_prop_DS, tab_prop_UV)

################
tab_data_final=rbind(tab_length, tab_cov, tab_data_prop)

tab_data_final2= tab_data_final %>% 
  pivot_longer(c(leeMs,MC_M,MC_MR,DS,UV), values_to = "label") %>%
  pivot_wider(names_from = c(beta,layer), values_from = label) %>%
  mutate(name = case_when(
    name == "beta" ~ "beta",
    name =="leeMs" ~ my_methods[1],
    name =="MC_M" ~ my_methods[2],
    name =="MC_MR" ~ my_methods[3],
    name =="DS" ~ my_methods[4],
    name =="UV" ~ my_methods[5],
  )) %>% select(name,starts_with('1'),starts_with('0.8'),starts_with('0.2'))

col_names=c("Method", rep(c("Median length", "Coverage", "Infinite length (%)"),3))
header_names = c("",setNames(3, "$\\beta_{i}=1$"),setNames(3, "$\\beta_{i}=0.1$"),setNames(3, "$\\beta_{i}=0$"))


kbl(tab_data_final2,align = 'c',booktabs = T, escape = F, caption = "",
    # format="latex"
    col.names = col_names) %>%
  add_header_above(header_names, bold = T) %>%
  kable_styling(position = "center", latex_options = c("hold_position", "scale_down"), font_size = 12)

#### Selection probability
data.frame(beta=my_results$beta_poly) %>%
  group_by(beta) %>%
  summarise(n = n()) %>%
  mutate(freq = n /300)

