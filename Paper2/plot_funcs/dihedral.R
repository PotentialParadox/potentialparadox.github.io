library(tidyverse)
setwd('~/Documents/paper2/dihedrals/')

collect_dihedrals <- function (runtime, infile, state, location, solvent) {
  data <- read_csv(infile, col_names = FALSE)
  nsteps <- ncol(data)
  time_labels <- (0:(nsteps-1))*runtime/(nsteps-1)
  names(data) <- as.character(time_labels)
  data %>% rowid_to_column("Trajectory") %>% 
    gather(key = "Time", value = "Angle", -Trajectory) %>% 
    group_by(Time) %>% 
    summarise(MeanAngle = mean(Angle)) %>% 
    mutate(State = state) %>%
    mutate(Location = location) %>% 
    mutate(Time = as.numeric(Time)) %>% 
    mutate(SolventID = solvent) %>% 
    mutate(Solvent = if_else(location == '', solvent, sprintf("%s-%s", solvent, location))) %>% 
    mutate(Description = if_else(location == '', state, sprintf("%s-%s", state, location)))
}

plot_dihedrals <- function (dihedral_data, description) {
  dihedral_data %>% ggplot(aes_string(x = "Time", y = "MeanAngle", color = description)) +
    geom_point() +
    geom_smooth(method = 'loess', formula = y ~ x) +
    theme_bw() +
    labs(x = "Time (ps)", y = "Angle (Degree)") +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15)
    )
}


#############################
# PPV3-NO2 Vacuum Comparison
#############################
# S0 & S1 - Vacuum from ~/backup1/TestRuns/Paper1/ppv-vacuum
# 128 Trajectories | 2000 steps
s0_ppv3_no2_vac_near <- collect_dihedrals(runtime = 10,
                                          infile = 'dihedral-ppv3-no2-s0-1.csv',
                                          state = "S0",
                                          solvent = "Vacuum",
                                          location = "Near") %>% 
  filter(Time > 9.5) %>% mutate(Time = Time - 9.5)
s0_ppv3_no2_vac_far <- collect_dihedrals(runtime = 10,
                                          infile = 'dihedral-ppv3-no2-s0-2.csv',
                                          state = "S0",
                                          solvent = "Vacuum",                                        
                                          location = "Far") %>% 
  filter(Time > 9.5) %>% mutate(Time = Time - 9.5)
s1_ppv3_no2_vac_near <- collect_dihedrals(runtime = 10,
                                          infile = 'dihedral-ppv3-no2-s1-1.csv',
                                          state = "S1",
                                          solvent = "Vacuum",     
                                          location = "Near") %>% 
  filter(Time <= 0.5)
s1_ppv3_no2_vac_far <- collect_dihedrals(runtime = 10,
                                         infile = 'dihedral-ppv3-no2-s1-2.csv',
                                         state = "S1",
                                         solvent = "Vacuum",
                                         location = "Far") %>% 
  filter(Time <= 0.5)

# Sm - Vacuum from ~/backup4/TestRuns/Paper2/ppv-no2-validation/vacuum ! Not the bla i needed more data
# Near angles are [[7,8,9,10],[16,15,14,13]
# Far angles are [[18,17,16,15], [16,15,14,13]]
# 260 Trajectories
sm_ppv3_no2_vac_near <- collect_dihedrals(runtime = 0.5,
                                         infile = 'dihedral-ppv3-no2-sm-1.csv',
                                         state = "Sm",
                                         solvent = "Vacuum",
                                         location = 'Near')

sm_ppv3_no2_vac_far <- collect_dihedrals(runtime = 0.5,
                                         infile = 'dihedral-ppv3-no2-sm-2.csv',
                                         state = "Sm",
                                         solvent = "Vacuum",
                                         location = "Far")
dihedral_data <- bind_rows(
  s0_ppv3_no2_vac_near,
  s0_ppv3_no2_vac_far,
  s1_ppv3_no2_vac_near,
  s1_ppv3_no2_vac_far,
  sm_ppv3_no2_vac_near,
  sm_ppv3_no2_vac_far
) %>% 
  mutate(Location = as.factor(Location)) %>%
  mutate(State = as.factor(State)) %>% 
  mutate(Description = as.factor(Description))

plot_dihedrals(dihedral_data = dihedral_data, description = "Description")

#################################
# PPV3-NO2 SM Solvent Comparison
#################################
# Sm - Vacuum from ~/backup4/TestRuns/Paper2/ppv-no2-validation/vacuum ! Not the bla i needed more data
# Near angles are [[7,8,9,10],[16,15,14,13]
# Far angles are [[18,17,16,15], [16,15,14,13]]
# 260 Trajectories
sm_ppv3_no2_vac_near <- collect_dihedrals(runtime = 0.5,
                                          infile = 'dihedral-ppv3-no2-sm-1.csv',
                                          state = "Sm",
                                          solvent = "Vacuum",
                                          location = 'Near')

sm_ppv3_no2_vac_far <- collect_dihedrals(runtime = 0.5,
                                         infile = 'dihedral-ppv3-no2-sm-2.csv',
                                         state = "Sm",
                                         solvent = "Vacuum",
                                         location = "Far")
# Sm - CH3OH from ~/backup4/TestRuns/Paper2/ppv-no2-validation/ch3oh/mm
# Near angles are [[7,8,9,10],[16,15,14,13]
# Far angles are [[18,17,16,15], [16,15,14,13]]
# 560 Trajectories
sm_ppv3_no2_ch3oh_near <- collect_dihedrals(runtime = 0.5,
                                         infile = 'dihedral-ppv3-no2-sm-ch3oh-1.csv',
                                         state = "Sm",
                                         solvent = "Methanol",
                                         location = "Near")
sm_ppv3_no2_ch3oh_far <- collect_dihedrals(runtime = 0.5,
                                           infile = 'dihedral-ppv3-no2-sm-ch3oh-2.csv',
                                           state = "Sm",
                                           solvent = "Methanol",
                                           location = "Far")
# Sm - CH3OH-5s from ~/backup4/TestRuns/Paper2/ppv-no2-validation/ch3oh/5s
# Near angles are [[7,8,9,10],[16,15,14,13]
# Far angles are [[18,17,16,15], [16,15,14,13]]
# 509 Trajectories
sm_ppv3_no2_ch3oh_5s_near <- collect_dihedrals(runtime = 0.5,
                                              infile = 'dihedral-ppv3-no2-sm-ch3oh-5s-1.csv',
                                              state = "Sm",
                                              solvent = "Methanol 5QM",
                                              location = "Near")
sm_ppv3_no2_ch3oh_5s_far <- collect_dihedrals(runtime = 0.5,
                                              infile = 'dihedral-ppv3-no2-sm-ch3oh-5s-2.csv',
                                              state = "Sm",
                                              solvent = "Methanol 5QM",
                                              location = "Far")

dihedral_data <- bind_rows(
  sm_ppv3_no2_vac_near,
  sm_ppv3_no2_vac_far,
  sm_ppv3_no2_ch3oh_near,
  sm_ppv3_no2_ch3oh_far,
  sm_ppv3_no2_ch3oh_5s_near,
  sm_ppv3_no2_ch3oh_5s_far
)
plot_dihedrals(dihedral_data, description = "Solvent")

#############################
# PPV3 Vacuum Comparison
#############################
# Sm Vacuum from ~/home/dustin/backup4/TestRuns/Paper2/ppv3-validation/vacuum
# Note that there were restart_0 (ds = 0.5fs, tf = 500fs) and restart_1 (ds = 0.1fs, ti=500fs, tf=1000fs)
# We only care about restart_0
# We can average all the dihedrals in this one
# [[11, 12, 15, 16], [18, 17, 16, 15], [10, 9, 8, 7], [5, 4, 7, 8]]
# 254 Trajectories
sm_ppv3_vac <- collect_dihedrals(runtime = 3,
                                 infile = 'dihedral-ppv3-sm.csv',
                                 state = "Sm",
                                 solvent = "Vacuum",
                                 location = '') %>% head(11)

# Sm from /home/dustin/backup4/TestRuns/Paper2/ppv3-validation/tammie_reproduction/vacuum-no-trivial
# This has a 0.1fs timestep
sm_ppv3_vac_no_trivial <- collect_dihedrals(runtime = 0.5,
                                            infile = 'dihedral-ppv3-sm-nt.csv',
                                            solvent = "Vacuum",
                                            state = "Sm No Trivial",
                                            location = '')

# From /home/dustin/backup4/TestRuns/Paper2/ppv3-validation/tammie_reproduction/vacuum
# 0.1fs timestep
sm_ppv3_vac_01 <- collect_dihedrals(runtime = 0.5,
                                 infile = 'dihedral-ppv3-sm-0.1fs.csv',
                                 state = "Sm 0.1fs",
                                 solvent = "Vacuum",
                                 location = '') %>% filter(Time <= 0.5)

# S0 & S1 - Vacuum from ~/backup1/TestRuns/Paper1/ppv-vacuum
# We can average all the dihedral in this one
# [[31, 24, 23, 22], [18, 21, 22, 23], [2, 4, 7, 8], [17, 9, 8, 7]]
# 540 Trajectories
s0_ppv3_vac <- collect_dihedrals(runtime = 20,
                                 infile = 'dihedral-ppv3-s0.csv',
                                 state = "S0",
                                 solvent = "Vacuum",
                                 location = '') %>% filter(Time >= 19.5) %>% mutate(Time = Time - 19.5)
s1_ppv3_vac <- collect_dihedrals(runtime = 1,
                                 infile = 'dihedral-ppv3-s1.csv',
                                 state = "S1",
                                 solvent = "Vacuum",
                                 location = '') %>% filter(Time <= 0.5)

dihedral_data <- bind_rows(
  s0_ppv3_vac,
  sm_ppv3_vac,
  sm_ppv3_vac_no_trivial,
  s1_ppv3_vac
)
plot_dihedrals(dihedral_data, description = "Description")
