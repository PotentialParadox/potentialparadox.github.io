library(tidyverse)

setwd('~/Documents/paper2/energies/ppv3_no2')

#################################################
# Collect Potentials
#################################################
## The goal will be to plot the average potential energies
## In order to this I'm going to have to pair the energy files
## with the coefficient file with the form
## | Trajectory | Time-fs | State

## The energies between different solvents will not be directly
## comparable to one another. However will can plot the differences
## from the initial s0 state, which can get from t=0, s=0.

collect_energy_deltas <- function(filename){
    read_csv(filename) %>%
        group_by(Trajectory) %>%
        mutate(Energy = PotentialeV - first(PotentialeV)) %>%
        rename(Timefs = Timefs) %>%
        mutate(Trajectory = as.numeric(Trajectory)) %>%
        select(Trajectory, Timefs, State, Energy)
}

get_current_state_deltas <- function(state_potential_filename, states_filename){
    all_deltas <- collect_energy_deltas(state_potential_filename)
    states <- read_csv(states_filename)

    states %>%
        left_join(all_deltas, by = c("Trajectory", "Timefs", "State")) %>%
        filter(!is.na(Energy)) %>%
        group_by(Timefs) %>%
        summarize(MeanEnergyeV = mean(Energy), .groups = 'drop')
}

vacuum_energies <- get_current_state_deltas('state-potentials-vacuum.csv', 'states-vacuum.csv') %>%
    mutate(Solvent = "Vacuum")

ch3oh_energies <- get_current_state_deltas('state-potentials-ch3oh.csv', 'states-ch3oh.csv') %>%
    mutate(Solvent = "0")

ch3oh_5s_energies <- get_current_state_deltas('state-potentials-ch3oh-5s.csv', 'states-ch3oh-5s.csv') %>%
    mutate(Solvent = "5")

ch3oh_10s_energies <- get_current_state_deltas('state-potentials-ch3oh-10s.csv', 'states-ch3oh-10s.csv') %>%
    mutate(Solvent = "10")

energies <- bind_rows(
    vacuum_energies,
    ch3oh_energies,
    ch3oh_5s_energies,
    ch3oh_10s_energies
)


energies %>%
    ggplot(aes(x = Timefs, y = MeanEnergyeV, color = Solvent)) +
    geom_line(size=1.5) +
    theme_bw() +
    labs(color = "Number QM Solvents",
         x = "Time (fs)",
         y = "Energy (eV)"
         ) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          legend.text = element_text(size = 18),
          legend.title = element_text(size=20),
          legend.text.align = 0,
          legend.position = "top"
    )
ggsave("~/potentialparadox.github.io/Paper2/Images/potential_energies/solvent_comparison.png", width = 10, height = 10)


#################################################
# Plot S1 Relaxations
#################################################
plot_potentials <- function (potentials, legend_breaks, legend_labels){
  potentials %>% 
    ggplot(aes(x = `Time-fs`, y = MeanEnergy, color = System)) +
    geom_point(size = 0.5) +
    labs(x="Time (fs)", y="Energy (eV)")+
    scale_color_discrete(breaks = legend_breaks, labels = legend_labels) +
    facet_grid(rows = vars(Solute), scales="free") +
    theme_bw() +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.text.align = 0,
          legend.position = c(0.8, 0.7)
    )
}

legend_breaks <- c("S0 Vacuum",
                   "S1 Vacuum",
                   "S1 Methanol"
                   )
legend_labels <- c(expression(paste(S[0], " Vacuum")),
                   expression(paste(S[1], " Vacuum")),
                   expression(paste(S[m], " Methanol")),
                   expression(paste(S[m], " Methanol-5s")))

plot_potentials(potentials, legend_breaks, legend_labels)



