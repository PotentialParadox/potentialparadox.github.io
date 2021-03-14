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
    mutate(Solvent = "Vacuum") %>%
    mutate(Type="Measured")

ch3oh_energies <- get_current_state_deltas('state-potentials-ch3oh.csv', 'states-ch3oh.csv') %>%
    mutate(Solvent = "0") %>%
    mutate(Type="Measured")

ch3oh_5s_energies <- get_current_state_deltas('state-potentials-ch3oh-5s.csv', 'states-ch3oh-5s.csv') %>%
    mutate(Solvent = "5") %>%
    mutate(Type="Measured")

ch3oh_10s_energies <- get_current_state_deltas('state-potentials-ch3oh-10s.csv', 'states-ch3oh-10s.csv') %>%
    mutate(Solvent = "10") %>%
    mutate(Type="Measured")

energies <- bind_rows(
    vacuum_energies,
    ch3oh_energies,
    ch3oh_5s_energies,
    ch3oh_10s_energies
)


energies %>%
    mutate(Solvent = factor(Solvent, levels=c("Vacuum", "0", "5", "10"))) %>%
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



model_decay <- function(params, data){
    A = params[1]
    B = params[2]
    kB = params[3]
    C = params[4]
    kC = params[5]
    A + B * exp(-kB * data$Timefs) + C * exp(-kC * data$Timefs)
}

model <- model_decay

measure_distance <- function(mod_params, data) {
    diff <- data$MeanEnergyeV - model(mod_params, data)
    sqrt(mean(diff^2))
}


best <- optim(c(2.5, 1, 0.005, 0.5, 0.005), measure_distance, data = ch3oh_energies)

s1_vacuum_fit <- purrr::partial(model, c(2.633, 0.9559, 0.0043293, 0.84108, 0.01262))
fit_vacuum_df <- tibble("Timefs" = vacuum_energies$Timefs, Solvent = "Vacuum", Type = "Predicted", "MeanEnergyeV" = s1_vacuum_fit(vacuum_energies))

s1_ch3oh_fit <- purrr::partial(model, c(2.671832, 1.67528, 0.006676))
fit_ch3oh_df <- tibble("Timefs" = ch3oh_energies$Timefs, Solvent = "0", Type = "Predicted", "MeanEnergyeV" = s1_ch3oh_fit(ch3oh_energies))

data <- bind_rows(vacuum_energies, fit_vacuum_df)

data %>%
    ggplot(aes(x=Timefs, y = MeanEnergyeV, color=Type)) +
    geom_line() +
    facet_wrap(~ Solvent)

x = log(ch3oh_energies$MeanEnergyeV)


plot(x)
