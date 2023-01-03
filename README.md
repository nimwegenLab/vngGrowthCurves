
install with `remotes::install_github('vanNimwegenLab/vngGrowthCurves')`.
NB: this requires authentication (cf `install_github` help; alternatively, clone locally and install from source).

A simple example would look something like:

```
myplates <- bind_rows(
  create_empty_plate(date='20221120') %>% 
    add_col_var(media, rep(c('M9+0.02glu', 'M9+0.2gly', 'M9+0.1ca'), times=c(4, 4, 4)) ) %>% 
    add_col_var(smm, rep(c(0, 0, 20, 20), 3) ) %>% 
    add_col_var(replicate, rep(c('A', 'B'), 6) ) %>% 
    add_row_var(strain, c(A="MG1655", B="DH5a", C="A1", D="A2", E="A4", F='C10', G='D6', H="E6") ),
    
  create_empty_plate(date='20221125') %>% 
    add_col_var(media, rep(c('M9+0.02glu', 'M9+0.2gly', 'M9+0.1ca'), times=c(4, 4, 4)) ) %>% 
    add_col_var(smm, rep(c(0, 20), each=6) ) %>% 
    add_col_var(replicate, rep(c('A', 'B'), 6) ) %>% 
    add_row_var(strain, c(A="MG1655", B="DH5a", C="A1", D="A2", E="A4", F='C10', G='D6', H="E6") ),
) 


myods <- bind_rows(
  read_Biotek_Synergy2_columns(here::here("path/to", "20221120_expt.txt"), date='20221120') %>% 
    filter(! well %in% c('D4') ),
    
  read_Biotek_Synergy2_columns(here::here("path/to", "20221125_expt.txt"), date='20221125'),
  ) %>% 
  group_by(date, well) %>% 
  nest() %>% 
  mutate(b = map(data, ~find_blank_od(.$time, .$value, .tmax=2*3600, .cv_thresh=.02)),
         blank_value=map_dbl(b, ~.$value), 
         blank_time=map_dbl(b, ~.$time), b=NULL) 


myods %>% 
  unnest(data) %>% 
  mutate(od=value-blank_value) %>% 
  left_join(myplates) %>% 
  ggplot(aes(time, od)) +
  facet_grid(date+row~col) +
  geom_line() +
  scale_y_log10() +
  ggCustomTJ::scale_x_hours(5, limits=c(NA, 20*3600)) +
  coord_cartesian(ylim=c(2e-6, NA)) +
  NULL
  
myods %>% 
  unnest(data) %>% 
  mutate(od=value-blank_value) %>% 
  left_join(myplates) %>% 
  ggplot(aes(time, od)) +
  facet_grid(row~media, labeller = labeller(row=c(A="MG1655", B="DH5a", C="A1", D="A2", E="A4", F='C10', G='D6', H="E6")) ) +
  geom_line(aes(col=interaction(date, replicate), lty=factor(smm), group=interaction(date, well))) +
  scale_y_log10() +
  ggCustomTJ::scale_x_hours(20) +
  coord_cartesian(ylim=c(2e-6, NA))
NULL

```