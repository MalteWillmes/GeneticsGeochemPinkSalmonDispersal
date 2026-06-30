#Data reduction and QAQC

# Prepare the environment -------------------------------------------------
library(tidyverse)
library(here) 
library(openxlsx) 
library(stringr) 
library(readxl)
library(zoo)

#Clear the workspace
rm(list = ls())

# Read in data ------------------------------------------------------------
#Fish meta
fish_meta <- read.xlsx(here("data", "sample_overview.xlsx"), sheet="meta")

#Dataframe with all imported strontium isotope data both adult and juvenile
path_data <-"data/20240106_ICP-MS_LA_Sr_Malte Willmes.xlsx"
sr_data_readin <-path_data %>%
  excel_sheets() %>% 
  #Extract only the sheets that have time series data
  str_subset(pattern = "time series", negate = FALSE)%>%
  set_names()%>%
  map_df(~ read_excel(path = path_data, sheet = .x, skip=1), .id = "sheet")
 

# Clean Sr isotope data --------------------------------------------------------------
#Flip profiles run from edge to core (in case)
fish_rev <- c("23-294","23-298","23-300","23-303","23-305",
              "23-306","23-310","23-312","23-318","23-321",
              "23-332","23-333","23-335","23-337","23-338","23-342",
              "23-344")
#Remove fish with failed analyses (too much vaterite, too low V)
fish_exclude <- c("KONG_22_0111","KONG_22_0109","KONG_22_0110",
                  "23-299","23-301a","23-302","23-319","23-320",
                  "23-324","23-326","23-331")


sr_data <-sr_data_readin%>%
  dplyr::select(iolite_id=sheet,duration=`Elapsed Time` , 
                sr87sr86=Sr87_86_CorrRb,totalSr=totalSrBeam)%>%
  mutate(sr87sr86=as.numeric(sr87sr86),
         totalSr=as.numeric(totalSr))%>%
  drop_na()%>%
  mutate(material_type=case_when(grepl("mollusk", iolite_id) ~ "Mollusk",
                                 grepl("Chalk", iolite_id) ~ "Chalk",
                                 grepl("Foram", iolite_id) ~ "Foram",
                                 grepl("BHVO", iolite_id) ~ "BHVO",
                                 grepl("BCR", iolite_id) ~ "BCR",
                                 TRUE~"PinkSalmon"))%>%
  mutate(material_group=case_when(material_type=="PinkSalmon"~"Sample",
                                  TRUE~"Reference"))%>%
  left_join(fish_meta, by="iolite_id")%>%
  group_by(fishID)%>%
  mutate(Distance=duration*scan_speed_ums)%>%
  mutate(Distance=case_when(fishID%in%fish_rev ~(-1)*(Distance-max(Distance)),
                   TRUE ~Distance))%>%
  ungroup()%>%
  #Remove outliers
  #Hard trim
  filter (sr87sr86>0.7 & sr87sr86<0.8)%>%
  filter (totalSr>0.2)%>%
  group_by(fishID)%>%
  #Outlier rejection
  mutate(median = rollapply(sr87sr86, width = 80, FUN = median, partial = T)) %>%
  mutate(upper75prob = rollapply(sr87sr86, width = 80, FUN = quantile, partial = T, probs = 0.75, na.rm=TRUE),
         lower25prob = rollapply(sr87sr86, width = 80, FUN = quantile, partial = T, probs = 0.25, na.rm=TRUE),
         IQR = upper75prob - lower25prob,
         diff_flag = (sr87sr86 - median)/IQR,
         despiked = ifelse(abs(diff_flag) > 2, median, sr87sr86)) %>%
  ungroup()%>%
  #Remove fish that were not correctly analyzed
  filter(!fishID%in%fish_exclude)%>%
  #Remove ref materials
  filter(material_group=="Sample")%>%
  group_by(fishID)%>%
  #Trim to Core
  #adjust for trim
  filter(Distance>=core_start)%>%
  mutate(Distance=Distance-min(Distance))%>%
  #Keep only profile with accurate data (trim of edges of vaterite/low V)
  filter(is.na(trim) | Distance <=trim)%>%
  ungroup()

#Spline data
sr_spline <- sr_data  %>%
  group_by(fishID)%>%
  nest()%>%
  mutate(despiked_spline =map(data, ~mgcv::gam(despiked ~ s(Distance, k=50, bs="tp"), data=.x))) %>%
  mutate(augment_spline= map2(despiked_spline, data, ~broom::augment(.x, newdata = .y, se_fit = TRUE))) %>%
  unnest(augment_spline)%>%
  mutate (sr_spline = .fitted,
          sr_spline_se =.se.fit)%>%
  dplyr::select(fishID, Distance, sr_spline, sr_spline_se) %>%
  mutate (sr_spline_se=sqrt((sr_spline_se^2)+(0.00005^2)))


sr_analyses <- sr_data%>%
  left_join(sr_spline, by=c("fishID", "Distance"))%>%
  #Assign roi
  group_by(fishID)%>%
  mutate(roi=case_when(Distance >= core_start & Distance <= core_stop ~"Core",
                       Distance >= natal_start & Distance <= natal_stop ~"Natal",
                       Distance > natal_stop ~"Adult",
                       TRUE ~"Profile"))%>%
  group_by(fishID, roi)%>%
  mutate(sr87sr86_roi_mean=mean(sr_spline),
         sr87sr86_roi_sd=sd(sr_spline))%>%
  ungroup()%>%
  arrange(fishID, Distance)
write.csv(sr_analyses, "outputs/sr_analyses.csv", row.names=F)

sr_analyses_summary <- sr_analyses %>%
  distinct(fishID, Lifestage, roi,sr87sr86_roi_mean, sr87sr86_roi_sd)
write.csv(sr_analyses_summary, "outputs/sr_analyses_summary.csv", row.names=F)

# Plot data ---------------------------------------------------------------
#Plot adult data
p_sr <- ggplot(data=sr_analyses%>%filter(Lifestage=="Adult"))+
  # geom_point(aes(x=Distance, y=sr87sr86), alpha=0.2, size=1, color="black")+
  geom_ribbon(aes(x=Distance, ymin = sr_spline-(2*sr_spline_se),
                  ymax =sr_spline+(2*sr_spline_se)),
              fill = "black", alpha=0.2)+
  geom_line(aes(x=Distance, y=sr_spline), color="black")+
  geom_hline(yintercept=0.70918,  color="steelblue",linetype="dashed")+
  geom_hline(yintercept=0.72169,  color="red",linetype="dashed")+
  geom_hline(yintercept=0.7173,  color="orange",linetype="dashed")+
  geom_vline(aes(xintercept=core_start), color="black", linetype="dashed")+
  geom_vline(aes(xintercept=core_stop), color="black", linetype="dashed")+
  geom_vline(aes(xintercept=natal_start), color="orange", linetype="dashed")+
  geom_vline(aes(xintercept=natal_stop), color="orange", linetype="dashed")+
  scale_x_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 6),
                     "Distance")+
  scale_y_continuous(expression(""^87*Sr*"/"^86*Sr), 
                     breaks = scales::pretty_breaks(n = 5))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position='none')+
  facet_wrap(~fishID, ncol=4)
p_sr
ggsave('outputs/p_sr.jpg',plot=p_sr, dpi = 350, height = 35, width = 30, units = 'cm')


#Plot juvenile data
p_sr_juv <- ggplot(data=sr_analyses%>%filter(Lifestage=="Juvenile"))+
  geom_point(aes(x=Distance, y=sr87sr86), alpha=0.2, size=1, color="black")+
  geom_ribbon(aes(x=Distance, ymin = sr_spline-(2*sr_spline_se), 
                  ymax =sr_spline+(2*sr_spline_se)),
              fill = "black", alpha=0.2)+
  geom_line(aes(x=Distance, y=sr_spline), color="black")+
  geom_hline(yintercept=0.70918,  color="steelblue",linetype="dashed")+
  geom_hline(yintercept=0.72169,  color="red",linetype="dashed")+
  geom_hline(yintercept=0.7173,  color="orange",linetype="dashed")+
  geom_vline(aes(xintercept=core_start), color="black", linetype="dashed")+
  geom_vline(aes(xintercept=core_stop), color="black", linetype="dashed")+
  # geom_vline(aes(xintercept=natal_start), color="orange", linetype="dashed")+
  # geom_vline(aes(xintercept=natal_stop), color="orange", linetype="dashed")+
  scale_x_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 6),
                     "Distance")+
  scale_y_continuous(expression(""^87*Sr*"/"^86*Sr), 
                     breaks = scales::pretty_breaks(n = 5))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position='none')+
  facet_wrap(~fishID, ncol=4)
p_sr_juv
ggsave('outputs/p_sr_juv.jpg',plot=p_sr_juv, dpi = 350, height = 30, width = 30, units = 'cm')

#Plot strontium isotope profiles
#Show which fish to plot)
fish_ids <- unique(sr_analyses$fishID)
fish_ids #Show list of fish that will be plotted

#Profiles run from core to dorsal edge
pdf("outputs/sr_iso_splined.pdf",16, 8)
for(i in fish_ids){
  plot_chem<- sr_analyses%>% filter (fishID==i)
  p <- ggplot(data=plot_chem)+
    geom_ribbon(aes(x=Distance, ymin = sr_spline-(2*sr_spline_se), ymax =sr_spline+(2*sr_spline_se)),
                fill = "black", alpha=0.2)+
    geom_line(aes(x=Distance, y=sr_spline), color="black")+
    geom_hline(aes(yintercept=0.70918), color="steelblue")+#Ocean¨
    geom_vline(aes(xintercept=natal_start), color="red", linetype="dashed")+
    geom_vline(aes(xintercept=natal_stop), color="red", linetype="dashed")+
    scale_x_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 5),
                       limits=c(0, 1100), "Distance (µm) from Core")+
    scale_y_continuous(expand=c(0,0),expression(""^87*Sr*"/"^86*Sr), limits= c(0.707, 0.730), 
                       breaks = scales::pretty_breaks(n = 8))+
    theme_bw() +
    theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    ggtitle(i)
  print(p)
}
dev.off()

