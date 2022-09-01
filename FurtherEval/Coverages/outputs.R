## Load required packages
# install.packages("pacman")
pacman::p_load(selectiveInference, prospectr, knockoff, mvtnorm, glmnet, stabs, MASS, dplyr,ggthemes, kableExtra, ggplot2, tibble, tidyverse,
               viridis, pbapply, data.table, progress, randomcoloR, cowplot, ggpubr,tcltk, ggh4x, patchwork, profvis, fGarch, mvtnorm, beepr)
source("fun.R")
#### Load all data
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName, verbose=TRUE)
  get(ls()[ls() != "fileName"])
}
## Change to B = 200 (first figure) ###################
my_results=loadRData(file="MC_R_1000_B_200.Rdata")

## Change to B = 1000 (second figure) ###################
# my_results=loadRData(file="MC_R_1000_B_1000.Rdata")

################################################################
### Plotting
my_output = list()
times=5
for(i in 1:times){
  my_results2= my_results[[i]]
  my_output[[i]] = data.frame(beta=my_results2$beta_poly,
                              cov_poly=my_results2$coverage_poly_theory,
                              cov_MC=my_results2$coverage_poly_ori,
                              cov_MC_R=my_results2$coverage_poly_ori_R,
                              cov_MC_R_more=my_results2$coverage_poly_ori_R_more,
                              iter=my_results2$iter) %>%
    group_by(beta)%>% 
    mutate(across(1:4, ~(cumsum(.x)/seq_along(.x)),.names = "cummean_{col}")) %>%
    pivot_longer(cols = 7:10, names_to="Methods",values_to = "Coverage") %>%
    mutate(replication=i)
}
my_output=rbindlist(my_output)

beta_names = c(
  expression(paste(beta[i],' = 1')),
  expression(paste(beta[i],' = 0.5')),
  expression(paste(beta[i],' = 0.2')))

method_names=c(expression(paste(Lee['M,s'])),
               expression(paste(MC['M'])),
               expression(paste(MC['M,r'])),
               expression(paste(MC['M,r+'])))

my_output=my_output %>% 
  mutate(replication=factor(replication),
         beta= factor(beta, 
                      levels=c("1", "0.5", "0.2"),
                      labels= beta_names), 
         Methods=factor(Methods,
                        levels=
                          c("cummean_cov_poly",
                            "cummean_cov_MC",
                            "cummean_cov_MC_R",
                            "cummean_cov_MC_R_more"),
                        labels=method_names))
###############################
## Numerical results
agg_results = rbindlist(my_results)
agg_results = data.frame(
  cov_theory=agg_results$coverage_poly_theory,
  cov_ori=agg_results$coverage_poly_ori,
  cov_R=agg_results$coverage_poly_ori_R,
  cov_R_more=agg_results$coverage_poly_ori_R_more,
  beta=agg_results$beta_poly) %>%
  pivot_longer(1:4, names_to="Methods",values_to = "Coverage") %>%
  group_by(beta, Methods) %>%
  summarise(cov=format(round(mean(Coverage),2),nsmall=2),
            MCE=format(round(sd(Coverage),2),nsmall=2))

agg_results=agg_results %>%
  mutate(beta = factor(beta, levels=c("1", "0.5", "0.2"),labels= beta_names),
         Methods=factor(Methods,levels=
                          c("cov_theory",
                            "cov_ori",
                            "cov_R",
                            "cov_R_more"),
                        labels=method_names),
         Cov=as.character(paste0("Mean cov: ", cov, 
                                 "\n",
                                 "MCE: ", MCE)
         ))

## Position the numerical results
agg_results=agg_results %>% mutate(iter=700, Coverage=0.55)

final_p=ggplot(my_output, aes(x = iter, y=Coverage))+
  geom_line(aes(color=replication))+
  geom_hline(yintercept = 0.90, lty=2)+
  facet_grid(rows = vars(beta), cols=vars(Methods), 
             labeller=label_parsed)+
  geom_text(data=agg_results, aes(label = Cov), size=4)+
  guides(color="none")+
  xlab("Replications")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=13),
        strip.text.y = element_text(size=13),
        axis.title = element_text(size=13))
final_p

## Save plots
# ggsave(final_p, filename = "MC_R_1000_B_200_p.pdf",  width = 230, height = 297, units = "mm")
# ggsave(final_p, filename = "MC_R_1000_B_1000_p.pdf",  width = 230, height = 297, units = "mm")