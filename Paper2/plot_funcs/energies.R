library(tidyverse)

setwd('~/Documents/paper2/energies/')
data <- read_csv('ppv3/energies.csv')

data_mean <- data %>% group_by(`Time-fs`, State) %>% 
  summarise(MeanEnergy = mean(`Energy-eV`)) %>% 
  mutate(State = as.factor(State))

#
# Reading Potential Energies
#
read_potential <- function (filein, system_name){
  read_csv(filein) %>% 
    group_by(`Time-fs`) %>% 
    summarise(MeanEnergy = mean(`Potential-Ev`)) %>% 
    mutate(System = system_name)
}

plot_potentials <- function (potentials, legend_breaks, legend_labels){
  potentials %>% 
    ggplot(aes(x = `Time-fs`, y = MeanEnergy, color = System)) +
    geom_point(size = 0.5) +
    labs(x="Time (fs)", y="Energy (eV)")+
    scale_color_discrete(breaks = legend_breaks, labels = legend_labels) +
    theme_bw() +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.text.align = 0
    )
}

#
# PPV3
#
pot_ppv3_no_trivial <- read_potential('ppv3/potential-no-trivial.csv', 'Vacuum (No Trivial)')
pot_ppv3_vacuum <- read_potential('ppv3/potential-vacuum.csv','Vacuum')
pot_ppv3_ch3oh <- read_potential('ppv3/potential-ch3oh.csv', 'Methanol')
pot_ppv3_ch3oh_5s <- read_potential('ppv3/potential-ch3oh-5s.csv', 'Methanol-5s')
pot_ppv3_vacuum_s0 <- read_csv('ppv3/potential-s0-vacuum.csv') %>% 
  filter(State == 0) %>% 
  group_by(`Time-fs`) %>% 
  summarise(MeanEnergy = mean(`Potential-eV`)) %>% 
  filter(`Time-fs` >= 19500) %>% mutate(`Time-fs` = `Time-fs` - 19500) %>% 
  mutate(System = 'S0 Vacuum')
pot_ppv3_vacuum_s1 <- read_csv('ppv3/potential-s1-vacuum.csv') %>% 
  filter(State == 1) %>% 
  group_by(`Time-fs`) %>% 
  summarise(MeanEnergy = mean(`Potential-eV`)) %>% 
  filter(`Time-fs` >= 500) %>% mutate(`Time-fs` = `Time-fs` - 500) %>% 
  mutate(System = 'S1 Vacuum')
  
# This energy is significanly lower because of the additional ch3oh atoms in the
# qm calculation. For comparison, we are going to adjust the values mean energy
# values by the initial value of the solute in mm
adjust_factor <- pot_ppv3_ch3oh_5s[[1,2]] - pot_ppv3_ch3oh[[1,2]]
pot_ppv3_ch3oh_5s <- pot_ppv3_ch3oh_5s %>% mutate(MeanEnergy = MeanEnergy - adjust_factor)

potentials <- bind_rows(pot_ppv3_no_trivial,
                        pot_ppv3_vacuum,
                        pot_ppv3_vacuum_s0,
                        pot_ppv3_vacuum_s1,
                        pot_ppv3_ch3oh,
                        pot_ppv3_ch3oh_5s)

legend_breaks <- c("Vacuum (No Trivial)",
                   "Vacuum",
                   "S0 Vacuum",
                   "S1 Vacuum",
                   "Methanol",
                   "Methanol-5s")
legend_labels <- c(expression(paste(S[m], "Vacuum NT")),
                   expression(paste(S[m], "Vacuum")),
                   expression(paste(S[0], " Vacuum")),
                   expression(paste(S[1], " Vacuum")),
                   expression(paste(S[m], " Methanol")),
                   expression(paste(S[m], " Methanol-5s")))

plot_potentials(potentials, legend_breaks, legend_labels)

#
# PPV3-NO2
#
pot_ppv3_no2_vacuum <- read_potential('ppv3_no2/potential-vacuum.csv','Vacuum')
pot_ppv3_no2_ch3oh <- read_potential('ppv3_no2/potential-ch3oh.csv', 'Methanol')
pot_ppv3_no2_ch3oh_5s <- read_potential('ppv3_no2/potential-ch3oh-5s.csv', 'Methanol-5s')
# The following are taken from 
pot_ppv3_no2_vacuum_s0 <- read_csv('ppv3_no2/potential-s0-vacuum.csv') %>% 
  filter(State == 0) %>% 
  group_by(`Time-fs`) %>% 
  summarise(MeanEnergy = mean(`Potential-eV`)) %>% 
  filter(`Time-fs` >= 9500) %>% mutate(`Time-fs` = `Time-fs` - 9500) %>% 
  mutate(System = 'S0 Vacuum')
pot_ppv3_no2_vacuum_s1 <- read_csv('ppv3_no2/potential-s1-vacuum.csv') %>% 
  filter(State == 1) %>% 
  group_by(`Time-fs`) %>% 
  summarise(MeanEnergy = mean(`Potential-eV`)) %>% 
  filter(`Time-fs` >= 9500) %>% mutate(`Time-fs` = `Time-fs` - 9500) %>% 
  mutate(System = 'S1 Vacuum')
pot_ppv3_no2_ch3oh_s1 <- read_csv('ppv3_no2/potential-s1-methanol.csv') %>% 
  filter(State == 1) %>% 
  group_by(`Time-fs`) %>% 
  summarise(MeanEnergy = mean(`Potential-eV`)) %>% 
  filter(`Time-fs` >= 9500) %>% mutate(`Time-fs` = `Time-fs` - 9500) %>% 
  mutate(System = 'S1 Methanol')

# This energy is significanly lower because of the additional ch3oh atoms in the
# qm calculation. For comparison, we are going to adjust the values mean energy
# values by the initial value of the solute in mm
adjust_factor <- pot_ppv3_no2_ch3oh_5s[[1,2]] - pot_ppv3_no2_ch3oh[[1,2]]
pot_ppv3_no2_ch3oh_5s <- pot_ppv3_no2_ch3oh_5s %>% mutate(MeanEnergy = MeanEnergy - adjust_factor)

potentials <- bind_rows(pot_ppv3_no2_vacuum,
                        pot_ppv3_no2_vacuum_s0,
                        pot_ppv3_no2_vacuum_s1,
                        pot_ppv3_no2_ch3oh,
                        pot_ppv3_no2_ch3oh_s1,
                        pot_ppv3_no2_ch3oh_5s)

legend_breaks <- c("Vacuum",
                   "S0 Vacuum",
                   "S1 Vacuum",
                   "Methanol",
                   "S1 Methanol",
                   "Methanol-5s")
legend_labels <- c(expression(paste(S[m], " Vacuum")),
                   expression(paste(S[0], " Vacuum")),
                   expression(paste(S[1], " Vacuum")),    
                   expression(paste(S[m], " Methanol")),
                   expression(paste(S[1], " Methanol")),
                   expression(paste(S[m], " Methanol-5s")))

plot_potentials(potentials, legend_breaks, legend_labels)


