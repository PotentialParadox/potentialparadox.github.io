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

plot_blas <- function(bond_data, legend_breaks, labels, description, state_colors){
  bond_data %>%
    ggplot(aes_string(x = "t", y="BLA", color = description)) +
    geom_line(alpha = 0.5) +
    geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs"), se=FALSE) +
    scale_color_manual(breaks = legend_breaks, labels = legend_labels, values = state_colors) +
    labs(x = "Time (ps)") +
    facet_wrap(vars(Solute)) +
    theme_bw() +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.text.align = 0,
          legend.position = "top"
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
ppv3_vac_bond_data <- bind_rows(
  s0_ppv3_vac,
  s1_ppv3_vac,
  sm_ppv3_vac
) %>% mutate(Solute = 'PPV3')




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
                                     state = 'Sm',
                                     location = "Near",
                                     solvent = "Vacuum")

sm_ppv3_no2_vac_far <- collect_blas(runtime = 1,
                                    infile = 'bla-ppv3-no2-sm-2.csv',
                                    state = 'Sm',
                                    location = "Far",
                                    solvent = "Vacuum")

ppv3_no2_vac_bond_data <- bind_rows(
  s0_ppv3_no2_vac_near,
  s0_ppv3_no2_vac_far,
  s1_ppv3_no2_vac_near,
  s1_ppv3_no2_vac_far,
  sm_ppv3_no2_vac_near,
  sm_ppv3_no2_vac_far
) %>% mutate(Solute = 'PPV3-NO2')

######################################
# Ploting for BLA Vacuum
######################################
bond_data <- bind_rows(
  ppv3_vac_bond_data,
  ppv3_no2_vac_bond_data
)

legend_breaks <- c(
  "S0",
  "S0-Near",
  "S0-Far",
  "S1",
  "S1-Near",
  "S1-Far",
  "Sm",
  "Sm-Near",
  "Sm-Far"
)
legend_labels <- c(
  expression(paste(S[0])),
  expression(paste(S[0], " Near")),
  expression(paste(S[0], " Far")),
  expression(paste(S[1])),
  expression(paste(S[1], " Near")),
  expression(paste(S[1], " Far")),
  expression(paste(S[m])),
  expression(paste(S[m], " Near")),
  expression(paste(S[m], " Far"))
)
my_colors <- c(
  "S0" = "#EE421D",
  "S0-Near" = "#681F0F",
  "S0-Far" = "#F0755A",
  "S1" = "#27AED6",
  "S1-Near" = "#299999",
  "S1-Far" = "#2ECDFC", 
  "Sm" = "#18A74C",
  "Sm-Near" = "#107836",
  "Sm-Far" = "#22E368"
)

plot_blas(
  bond_data,
  legend_breaks = legend_breaks,
  labels = legend_labels,
  description = "Description",
  state_colors = my_colors
)

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
ppv3_sm_bond_data <- bind_rows(
  sm_ppv3_vac
) %>% mutate(Solute = 'PPV3')

plot_blas(bond_data, "Solvent")