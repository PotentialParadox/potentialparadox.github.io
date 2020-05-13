library(tidyverse)
setwd('~/Documents/paper2/bla')

collect_blas <- function (runtime, timestep, infile, state){
  data <- read_csv(infile, col_names = FALSE)
  nsteps <- ncol(data) - 1
  time_labels <- (0:(nsteps-1))*runtime/(nsteps)
  names(data) <- c('Bond', time_labels)
  s1_bond_data <- data %>% gather(key = "t", value = "Length", -Bond) %>% 
    group_by(Bond, t) %>% 
    summarise(MeanLength = mean(Length)) %>% 
    spread(Bond, MeanLength) %>% 
    mutate(BLA = (`1`+`3`)/2 - `2`) %>% 
    mutate(t = as.numeric(t)) %>% 
    mutate(State = state)
  s1_bond_data
}

plot_blas <- function(bond_data){
  bond_data %>%
    filter(t < 1) %>% 
    ggplot(aes(x = t, y=BLA, color = State)) +
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
# PPV3
######################
timestep = 5 #fs
# S0 & S1 - Vacuum from ~/backup1/TestRuns/Paper1/ppv-vacuum
# Bonds (24,23), (23,22), (22, 21)
# Bonds (4,7), (7,8), (8,9)
s0_ppv3_vac <- collect_blas(runtime = 4,
                            timestep = timestep,
                            infile = 'bla-ppv3-s0.csv',
                            state = 'S0')
s1_ppv3_vac <- collect_blas(runtime = 1,
                            timestep = timestep,
                            infile = 'bla-ppv3-s1.csv',
                            state = 'S1')
# Sm - Vacuum from ~/backup4/TestRuns/Paper2/ppv3-validation/vacuum-bla
# Bonds (17,16), (16, 15), (15, 12)
# Bonds (4,7), (7, 8), (8, 9)
sm_ppv3_vac <- collect_blas(runtime = 1,
                            timestep = timestep,
                            infile = 'bla-ppv3-sm.csv',
                            state = 'Sm')
bond_data <- bind_rows(
  s0_ppv3_vac,
  s1_ppv3_vac,
  sm_ppv3_vac
)
plot_blas(bond_data = bond_data)

######################
# PPV3 - NO2
######################
timestep = 5 #fs
# S0 & S1 - Vacuum from ~/backup1/TestRuns/Paper1/ppvno2-vacuum
# Near to NO2 Bonds: (6, 7), (7, 8), (8, 9)
# Far from NO2 Bonds: (17, 16) (16, 15), (15, 14)
s0_ppv3_no2_vac <- collect_blas(runtime = 10,
                                timestep = timestep,
                                infile = 'bla-ppv3-no2-s0.csv',
                                state = 'S0') %>% 
  mutate(t = t - 4) %>% 
  filter(t > 0)
s1_ppv3_no2_vac_near <- collect_blas(runtime = 10,
                                timestep = timestep,
                                infile = 'bla-ppv3-no2-s1-1.csv',
                                state = 'S1-near')
s1_ppv3_no2_vac_far <- collect_blas(runtime = 10,
                                     timestep = timestep,
                                     infile = 'bla-ppv3-no2-s1-2.csv',
                                     state = 'S1-far')
# Sm - Vacuum from ~/backup4/TestRuns/Paper2/ppv-no2-validation/vacuum-bla
# Bonds same
sm_ppv3_no2_vac_near <- collect_blas(runtime = 1,
                                     timestep = timestep,
                                     infile = 'bla-ppv3-no2-sm-1.csv',
                                     state = 'Sm-near')
sm_ppv3_no2_vac_far <- collect_blas(runtime = 1,
                                    timestep = timestep,
                                    infile = 'bla-ppv3-no2-sm-2.csv',
                                    state = 'Sm-far')

bond_data <- bind_rows(
  s0_ppv3_no2_vac,
  s1_ppv3_no2_vac_near,
  s1_ppv3_no2_vac_far,
  sm_ppv3_no2_vac_near,
  sm_ppv3_no2_vac_far
)
plot_blas(bond_data)
