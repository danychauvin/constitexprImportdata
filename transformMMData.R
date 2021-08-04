# Setting strains
myframes_to_mycells <- myframes_to_mycells %>% 
  mutate(strain='MG1655')
# Converting GFP units
myframes_to_mycells <- myframes_to_mycells %>%
  # append relevant conversion parameters
  mutate(autofluo_predict=NA, fp_per_dn=NA,
         # case of GFPmut2 (from zaslaver library)
         autofluo_predict = ifelse(strain!='asc662', 133.6, autofluo_predict),
         fp_per_dn = ifelse(strain!='asc662', 0.198, fp_per_dn),
         # case of asc662
         autofluo_predict = ifelse(strain=='asc662', 422.8, autofluo_predict),
         fp_per_dn = ifelse(strain=='asc662', 0.0361 * 4, fp_per_dn)) %>% 
  # convert to gfp units (after subtracting autofluorescence)
  mutate(fluogfp_amplitude = fluo_amplitude - autofluo_predict * length_um,
         gfp_nb = fluogfp_amplitude * fp_per_dn ) %>% 
  group_by(cell) %>% 
  arrange(time_sec) %>% 
  mutate(g_birth=first(gfp_nb)) %>% 
  mutate(g_div=last(gfp_nb))
cell_fluo_info <- myframes_to_mycells %>% 
  distinct(cell,g_birth,g_div)


# Propagating fluorescence information to full cell cycle
mycells <- mycells %>% 
  left_join(cell_fluo_info,by=c("cell")) %>% 
  mutate(c_birth=g_birth/l_birth,
         c_div=g_div/l_div,
         c_birth_v=g_birth/volume_birth,
         c_div_v=g_div/volume_div,
         dg=g_div - g_birth,
         dcdt=(g_div/l_div - g_birth/l_birth) / div_time,
         g=log(g_div/g_birth) / div_time, # would be cleaner from a fit but difficult to compare
         gamma=dg/dl,
         q=dg/dl*alpha)
# Change factor values
mycells$condition <- factor(mycells$condition,levels=c("acetate","glycerol","glucose","glucoseaa"))
myframes_to_mycells$condition <- factor(myframes_to_mycells$condition,levels=c("acetate","glycerol","glucose","glucoseaa"))
#myframes_to_mycells$promoter <- factor(myframes_to_mycells$promoter,levels=c("hi1","hi3","med2","med3","rplN","rpmB","rpsB","rrnB","p3F9"))
# Add promoters to mycells
cell_to_promoter <- myframes_to_mycells %>% 
  distinct(cell,promoter,vector)

mycells <- mycells %>% 
  select(-c(vector)) %>% 
  left_join(cell_to_promoter,by=c("cell")) %>% 
  mutate(strain=paste(vector,promoter,sep="-"))

time_interval_df <- myframes_to_mycells %>% 
  mutate(t_interval=t_interval*60) %>% 
  distinct(condition,t_interval)

# Computing pseudo-growth rate
mycells <- mycells %>% 
  group_by(condition,strain) %>% 
  mutate(gr_elongation=mean(alpha*3600/log(2)),
         gr_generation=1/mean(div_time/3600,na.rm = TRUE),
         pred_gr=dichotomic_search_growth_rate(alpha,"alpha")*3600/log(2)) %>% 
  mutate(sd_pred_gr=mean(sd_alpha,na.rm=TRUE)) %>%
  ungroup()

# Computing pseudo-concentration
myconcentrations <- myframes_to_mycells %>% 
  group_by(cell) %>% 
  mutate(mean_concentration=mean(gfp_nb/volume_um),
         sd_concentration=sd(gfp_nb/volume_um)) %>% 
  ungroup() %>% 
  distinct(cell,mean_concentration,sd_concentration)

mycells <- mycells %>% 
  left_join(myconcentrations,by=c("cell"))

mycells <- mycells %>% 
  group_by(condition,promoter,vector) %>% 
  mutate(mean_concentration_condition=mean(mean_concentration)) %>% 
  mutate(sd_concentration_condition=sd(mean_concentration)/n()) %>% 
  ungroup()

