#Data reduction and QAQC

# Prepare the environment -------------------------------------------------
library(tidyverse)
library(here) 
library(openxlsx) 
library(stringr) 
library(tidymodels)
library(viridis)
library(cowplot)
library(corrplot)
library(readxl)
library(zoo)
tidymodels_prefer()

#Clear the workspace
rm(list = ls())

# Read in data ------------------------------------------------------------
#Sr isotope mean roi data
sr_data <- read.csv(here("outputs", "sr_analyses_summary.csv"))
water_data <- read.xlsx(here("data", "EVA_water_samples.xlsx"), sheet="data")

#Mean natal Sr isotope data
sr_data_natal <- sr_data%>%
  filter(roi=="Natal")%>%
  #Manual classification if fish is within the Sr isotope range of Kongsfjordelva
  mutate(classification = case_when(sr87sr86_roi_mean>=0.719294& 
                                      sr87sr86_roi_mean<=0.721690 ~"Local",
                                    TRUE~ "Non-local"))

#Isoscape data
juv_join <-sr_data_natal %>%
  filter(Lifestage=="Juvenile")%>%
  dplyr::select(Sr87Sr86=sr87sr86_roi_mean)%>%
  mutate(River="Kongsfjordelva", sample_material="Juvenile otolith")%>%
  mutate(north=70.654769, east=29.247048)

isoscape_data <- water_data %>%
  filter(rep==1)%>%
  dplyr::select(River,north, east,Sr87Sr86)%>%
  mutate(sample_material="Water")%>%
  rbind(juv_join)%>%
  mutate(sample_material=factor(sample_material, levels=c("Water","Juvenile otolith")))%>%
  arrange((Sr87Sr86))

# Plot data ---------------------------------------------------------------

#Water isoscape
p_isoscape <- ggplot(data=isoscape_data)+
  geom_jitter(aes (x=reorder(River,east, fun=min), y=Sr87Sr86,fill=sample_material, shape= sample_material),
              color="black",size=3, alpha=0.6)+
  xlab("Rivers")+
  scale_shape_manual(name="Sample material", values=c(23,21))+
  scale_fill_manual(name="Sample material", values=c("darkblue","purple"))+
  scale_y_continuous(expression(""^87*Sr*"/"^86*Sr), 
                     breaks = scales::pretty_breaks(n = 5))+
  geom_hline(yintercept = 0.70918, linetype="dotted", color="black")+  #Mean Ocean   
  annotate("text",x =8, y = 0.70938, label = "Mean Global Ocean")+
  geom_hline(yintercept=min(juv_join$Sr87Sr86),  color="steelblue",linetype="dashed")+
  geom_hline(yintercept=max(juv_join$Sr87Sr86),  color="steelblue",linetype="dashed")+   
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position='top',
         axis.text.x = element_text(angle = 90, vjust = 0, hjust=1))
p_isoscape
ggsave(plot=p_isoscape, "outputs/isoscape.png", width = 8, height = 6)

# Natal classification based on Sr range ---------------------------------------

#Plot mean natal adult and juvenile data
p_sr_natal<- ggplot(data=sr_data_natal)+
  geom_hline(yintercept=min(juv_join$Sr87Sr86),  color="steelblue",linetype="dashed")+
  geom_hline(yintercept=max(juv_join$Sr87Sr86),  color="steelblue",linetype="dashed")+
  geom_jitter(aes(x=Lifestage, y=sr87sr86_roi_mean, fill=Lifestage),
              size=3, alpha=0.8, width=0.2, shape=21)+
  scale_y_continuous(expression(""^87*Sr*"/"^86*Sr), 
                     breaks = scales::pretty_breaks(n = 5))+
  scale_fill_manual(values=c("orange","purple"))+
  scale_shape_manual(values=c(21,22))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position='top')
p_sr_natal
ggsave('outputs/p_sr_natal.jpg',plot=p_sr_natal, dpi = 350, height = 10, width = 10, units = 'cm')


#Combined plot
p_natal_comb <-plot_grid(p_sr_natal, p_isoscape, labels = c('A', 'B'),
                                      rel_widths =c(0.4,0.6),nrow = 1)
p_natal_comb 
ggsave(plot=p_natal_comb , "outputs/p_natal_comb.png", width = 20, height = 10,dpi = 350, units = 'cm')

#Natal summary table
write.csv(sr_data_natal, "outputs/sr_data_natal.csv", row.names=F)




#Genetics plot
genetics_data <- read.csv(here("data", "genetics_data.csv"))

figsib<-genetics_data %>%
  mutate(group=factor(group, levels=c("random","observed")))

a<-binom.test(7, 13, p = 0.08, alternative = "two.sided")
figsib$ci_low<-c(NA,NA,NA,a$conf.int[1])
figsib$ci_high<-c(NA,NA,NA,a$conf.int[2])

p_figsib <-ggplot(subset(figsib, comparison=="within"), aes(x=group, y=props))+
  geom_bar(stat="identity")+
  ylab("Proportion of sib-pairs within rivers")+
  xlab("")+
  theme_bw(base_size = 16)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position='none')+
  scale_y_continuous(limits=c(0,1),expand = c(0, 0))+ 
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),width=0.2)
p_figsib
ggsave(plot=p_figsib , "outputs/p_figsib.png", dpi = 350, height = 10, width = 10, units = 'cm') 

