library(tidyverse)
library(rio)
library(mosaic)
setwd("~/Documents/paper2/populations/ppv3-no2")

##################################################
## Comparing the populations decay of various solvents
##################################################
read_pops <- function(filename, nQm_Solvent) {
    wide <- as_tibble(import(filename)) %>%
        select(V1, V2, V3, V4)
    names(wide) <- c("t", "Sm", "S1", "S2")
    final <- wide %>%
        pivot_longer(-t, names_to = "State") %>%
        mutate(nQM_Solvent = nQm_Solvent)
}

vacuum <- read_pops("vacuum/all_pops.txt", nQm_Solvent = 'Vacuum')
ch3oh <- read_pops("ch3oh/all_pops.txt", nQm_Solvent = '0')
ch3oh_5s <- read_pops("ch3oh_5s/all_pops.txt", nQm_Solvent = '5')
ch3oh_10s <- read_pops("ch3oh_10s/all_pops.txt", nQm_Solvent = '10')

ch3oh_comparison <- bind_rows(
    vacuum,
    ch3oh,
    ch3oh_5s,
    ch3oh_10s
) %>%
    mutate(State = as.factor(State)) %>%
    mutate(nQM_Solvent = as.factor(nQM_Solvent)) %>%
    rename(Population = value)

ch3oh_comparison %>%
    ggplot(aes(x = t, y = Population, color = nQM_Solvent)) +
    geom_line(size=1.5) +
    facet_grid(~ State) +
    theme_bw() +
    labs(color = "Number QM Solvents",
         x = "Time (fs)"
         ) +
    theme(legend.position = "top",
          legend.text = element_text(size = 12)
          )
ggsave("~/potentialparadox.github.io/Paper2/Images/populations/solvent_comparison.png", width = 16, height = 6)


##################################################
## Fitting Tammie's Equations
##################################################

tammies_fit_func <- function(A, tau, t){
    (A * exp(t / tau)) / (A + exp(t / tau)) - A / (1 + A)
}


vacuum_fit <- purrr::partial(tammies_fit_func, A = 1.62, tau = 100)
vacuum %>%
    filter(State == "S1") %>%
    ggplot(aes(x = t, y = value)) +
    geom_line(color = "blue") +
    stat_function(fun = vacuum_fit) + xlim(-1000,1000)

    
