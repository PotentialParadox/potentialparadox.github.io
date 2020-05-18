library(tidyverse)
setwd('~/Documents/paper2/bla')

collect_blas <- function (runtime, infile, state, location, solvent){
  data <- read_csv(infile, col_names = FALSE)
  nsteps <- ncol(data) - 1
  time_labels <- (0:(nsteps-1))*runtime/(nsteps -1)
  names(data) <- c('Bond', time_labels)
  s1_bond_data <- data %>% gather(key = "t", value = "Length", -Bond) %>% 
    group_by(Bond, t) %>% 
    summarise(MeanLength = mean(Length)) %>% 
    spread(Bond, MeanLength) %>% 
    mutate(BLA = (`1`+`3`)/2 - `2`) %>% 
    mutate(t = as.numeric(t)) %>% 
    mutate(State = state) %>% 
    mutate(Location = location) %>% 
    mutate(SolventID = solvent) %>% 
    mutate(Solvent = if_else(location == '', solvent, sprintf("%s-%s", solvent, location))) %>% 
    mutate(Description = if_else(location == '', state, sprintf("%s-%s", state, location)))
  s1_bond_data
}

plot_blas <- function(bond_data, description){
  bond_data %>%
    ggplot(aes_string(x = "t", y="BLA", color = description)) +
    geom_line(alpha = 0.5) +
    geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
    labs(x = "Time (ps)") +
    theme_bw() +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15)
    )
}
######################
# PPV3 Vacuum
######################
# S0 & S1 - Vacuum from ~/backup1/TestRuns/Paper1/ppv-vacuum
# Bonds (24,23), (23,22), (22, 21)
# Bonds (4,7), (7,8), (8,9)
s0_ppv3_vac <- collect_blas(runtime = 4,
                            infile = 'bla-ppv3-s0.csv',
                            state = 'S0',
                            location = '',
                            solvent = "Vacuum") %>% filter(t >=3) %>% mutate(t = t - 3)
s1_ppv3_vac <- collect_blas(runtime = 1,
                            infile = 'bla-ppv3-s1.csv',
                            state = 'S1',
                            location = '',
                            solvent = "Vacuum")
# Sm - Vacuum from ~/backup4/TestRuns/Paper2/ppv3-validation/vacuum-bla
# Bonds (17,16), (16, 15), (15, 12)
# Bonds (4,7), (7, 8), (8, 9)
sm_ppv3_vac <- collect_blas(runtime = 1,
                            infile = 'bla-ppv3-sm-bla.csv',
                            state = 'Sm',
                            location = '',
                            solvent = "Vacuum")
bond_data <- bind_rows(
  s0_ppv3_vac,
  s1_ppv3_vac,
  sm_ppv3_vac
)
plot_blas(bond_data = bond_data, "State")


######################
# PPV3 Sm
######################
# Sm - Vacuum from ~/backup4/TestRuns/Paper2/ppv3-validation/vacuum
# Bonds (17,16), (16, 15), (15, 12)
# Bonds (4,7), (7, 8), (8, 9)
# 254 Trajectories
sm_ppv3_vac <- collect_blas(runtime = 1,
                            infile = 'bla-ppv3-sm-bla.csv',
                            state = 'Sm',
                            location = '',
                            solvent = "Vacuum")

# Sm - CH3OH from ~/backup4/TestRuns/Paper2/ppv3-validation/vacuum
# central benzen rings are atoms 12, 13, 14, 9, 10, 11
bond_data <- bind_rows(
  sm_ppv3_vac
)

plot_blas(bond_data, "Solvent")

######################
# PPV3 - NO2
######################
timestep = 5 #fs
# S0 & S1 - Vacuum from ~/backup1/TestRuns/Paper1/ppvno2-vacuum
# Near to NO2 Bonds: (6, 7), (7, 8), (8, 9)
# Far from NO2 Bonds: (17, 16) (16, 15), (15, 14)
s0_ppv3_no2_vac_near <- collect_blas(runtime = 10,
                                infile = 'bla-ppv3-no2-s0-1.csv',
                                state = 'S0',
                                location = "Near",
                                solvent = "Vacuum") %>% 
  mutate(t = t - 9) %>% 
  filter(t > 0)
s0_ppv3_no2_vac_far <- collect_blas(runtime = 10,
                                     infile = 'bla-ppv3-no2-s0-2.csv',
                                     state = 'S0',
                                     location = "Far",
                                     solvent = "Vacuum") %>% 
  mutate(t = t - 9) %>% 
  filter(t > 0)

s1_ppv3_no2_vac_near <- collect_blas(runtime = 10,
                                infile = 'bla-ppv3-no2-s1-1.csv',
                                state = 'S1',
                                location = 'Near',
                                solvent = "Vacuum") %>% filter(t <= 1)
s1_ppv3_no2_vac_far <- collect_blas(runtime = 10,
                                    infile = 'bla-ppv3-no2-s1-2.csv',
                                    state = 'S1',
                                    location = 'Far',
                                    solvent = 'Vacuum') %>% filter(t <= 1)
# Sm - Vacuum from ~/backup4/TestRuns/Paper2/ppv-no2-validation/vacuum-bla
# Bonds same
sm_ppv3_no2_vac_near <- collect_blas(runtime = 1,
                                     infile = 'bla-ppv3-no2-sm-1.csv',
                                     state = 'Sm-near',
                                     location = "Near",
                                     solvent = "Vacuum")

sm_ppv3_no2_vac_far <- collect_blas(runtime = 1,
                                    infile = 'bla-ppv3-no2-sm-2.csv',
                                    state = 'Sm-far',
                                    location = "Far",
                                    solvent = "Vacuum")

bond_data <- bind_rows(
  s0_ppv3_no2_vac_near,
  s0_ppv3_no2_vac_far,
  s1_ppv3_no2_vac_near,
  s1_ppv3_no2_vac_far,
  sm_ppv3_no2_vac_near,
  sm_ppv3_no2_vac_far
)
plot_blas(bond_data, description = "Description")
