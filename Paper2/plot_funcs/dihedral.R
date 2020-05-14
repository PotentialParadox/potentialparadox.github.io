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
    mutate(Solvent = sprintf("%s-%s", solvent, location)) %>% 
    mutate(Description = sprintf("%s-%s", state, location))
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
  filter(Time > 9) %>% mutate(Time = Time - 9)
s0_ppv3_no2_vac_far <- collect_dihedrals(runtime = 10,
                                          infile = 'dihedral-ppv3-no2-s0-2.csv',
                                          state = "S0",
                                          solvent = "Vacuum",                                        
                                          location = "Far") %>% 
  filter(Time > 9) %>% mutate(Time = Time - 9)
s1_ppv3_no2_vac_near <- collect_dihedrals(runtime = 10,
                                          infile = 'dihedral-ppv3-no2-s1-1.csv',
                                          state = "S1",
                                          solvent = "Vacuum",     
                                          location = "Near") %>% 
  filter(Time <= 1)
s1_ppv3_no2_vac_far <- collect_dihedrals(runtime = 10,
                                         infile = 'dihedral-ppv3-no2-s1-2.csv',
                                         state = "S1",
                                         solvent = "Vacuum",
                                         location = "Far") %>% 
  filter(Time <= 1)

# Sm - Vacuum from ~/backup4/TestRuns/Paper2/ppv-no2-validation/vacuum ! Not the bla i needed more data
# Near angles are [[7,8,9,10],[16,15,14,13]
# Far angles are [[18,17,16,15], [16,15,14,13]]
# 260 Trajectories
sm_ppv3_no2_vac_near <- collect_dihedrals(runtime = 1,
                                         infile = 'dihedral-ppv3-no2-sm-1.csv',
                                         state = "Sm",
                                         solvent = "Vacuum",
                                         location = 'Near')

sm_ppv3_no2_vac_far <- collect_dihedrals(runtime = 1,
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

plot_dihedrals(dihedral_data = dihedral_data)

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
# Sm - Vacuum from ~/backup4/TestRuns/Paper2/ppv-no2-validation/vacuum ! Not the bla i needed more data
# Near angles are [[7,8,9,10],[16,15,14,13]
# Far angles are [[18,17,16,15], [16,15,14,13]]
# 260 Trajectories
sm_ppv3_no2_vac_near <- collect_dihedrals(runtime = 1,
                                          infile = 'dihedral-ppv3-no2-sm-1.csv',
                                          state = "Sm",
                                          solvent = "Vacuum",
                                          location = 'Near')

sm_ppv3_no2_vac_far <- collect_dihedrals(runtime = 1,
                                         infile = 'dihedral-ppv3-no2-sm-2.csv',
                                         state = "Sm",
                                         solvent = "Vacuum",
                                         location = "Far")