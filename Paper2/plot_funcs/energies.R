library(tidyverse)

setwd('~/Documents/paper2/energies/')
data <- read_csv('ppv3/energies.csv')

data_mean <- data %>% group_by(`Time-fs`, State) %>% 
  summarise(MeanEnergy = mean(`Energy-eV`)) %>% 
  mutate(State = as.factor(State))

data_mean %>% 
  ggplot(aes(x = `Time-fs`, y = MeanEnergy, group=State, color=State)) +
  geom_line()

#
# Reading Potential Energies
#
read_potential <- function (filein, system_name){
  read_csv(filein) %>% 
    group_by(`Time-fs`) %>% 
    summarise(MeanEnergy = mean(`Potential-Ev`)) %>% 
    mutate(System = system_name)
}

plot_potentials <- function (potentials){
  potentials %>% 
    ggplot(aes(x = `Time-fs`, y = MeanEnergy, color = System)) +
    geom_point() +
    labs(x="Time (fs)", y="Energy (eV)")+
    theme_bw() +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15)
    )
}

plot(1:10, xlab=expression('hi'[5]*'there'[6]^8*'you'[2]))
#
# PPV3
#
pot_ppv3_no_trivial <- read_potential('ppv3/potential-no-trivial.csv', 'Vacuum (No Trivial)')
pot_ppv3_vacuum <- read_potential('ppv3/potential-vacuum.csv','Vacuum')
pot_ppv3_ch3oh <- read_potential('ppv3/potential-ch3oh.csv', expression('CH3OH'))
pot_ppv3_ch3oh_5s <- read_potential('ppv3/potential-ch3oh-5s.csv', 'CH3OH-5s')
# This energy is significanly lower because of the additional ch3oh atoms in the
# qm calculation. For comparison, we are going to adjust the values mean energy
# values by the initial value of the solute in mm
adjust_factor <- pot_ppv3_ch3oh_5s[[1,2]] - pot_ppv3_ch3oh[[1,2]]
pot_ppv3_ch3oh_5s <- pot_ppv3_ch3oh_5s %>% mutate(MeanEnergy = MeanEnergy - adjust_factor)

potentials <- bind_rows(pot_ppv3_no_trivial,
                        pot_ppv3_vacuum,
                        pot_ppv3_ch3oh,
                        pot_ppv3_ch3oh_5s)
plot_potentials(potentials)

#
# PPV3-NO2
#
pot_ppv3_no2_vacuum <- read_potential('ppv3_no2/potential-vacuum.csv','Vacuum')
pot_ppv3_no2_ch3oh <- read_potential('ppv3_no2/potential-ch3oh.csv', 'CH3OH')
pot_ppv3_no2_ch3oh_5s <- read_potential('ppv3_no2/potential-ch3oh-5s.csv', 'CH3OH-5s')
# This energy is significanly lower because of the additional ch3oh atoms in the
# qm calculation. For comparison, we are going to adjust the values mean energy
# values by the initial value of the solute in mm
adjust_factor <- pot_ppv3_no2_ch3oh_5s[[1,2]] - pot_ppv3_no2_ch3oh[[1,2]]
pot_ppv3_no2_ch3oh_5s <- pot_ppv3_no2_ch3oh_5s %>% mutate(MeanEnergy = MeanEnergy - adjust_factor)

potentials <- bind_rows(pot_ppv3_no2_vacuum,
                        pot_ppv3_no2_ch3oh,
                        pot_ppv3_no2_ch3oh_5s)
plot_potentials(potentials)
