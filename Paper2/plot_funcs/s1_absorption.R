setwd("~/Documents/paper2/pump_pulse")
library(tidyverse)

get_spectra <- function(filename){
    read_csv(filename) %>%
        select(-EnergyNM) %>%
        pivot_longer(-EnergyEV, names_to="State") %>%
        filter(State != 'S1')
}

vacuum <- get_spectra('vacuum.csv') %>%
    mutate(nSolvent = "Vacuum") %>%
    filter(State %in% c('S9', 'S10')) %>%
    filter(EnergyEV > 3.5 & EnergyEV < 5)
ch3oh <- get_spectra('ch3oh.csv') %>%
    mutate(nSolvent = "0") %>%
    filter(State %in% c('S9', 'S10')) %>%
    filter(EnergyEV > 3.5 & EnergyEV < 5)
ch3oh_5s <- get_spectra('ch3oh-5s.csv') %>%
    mutate(nSolvent = "5") %>%
    filter(State %in% c('S9', 'S10')) %>%
    filter(EnergyEV > 3.5 & EnergyEV < 5)
ch3oh_10s <- get_spectra('ch3oh-10s.csv') %>%
    mutate(nSolvent = "10") %>%
    filter(State %in% c('S9', 'S10')) %>%
    filter(EnergyEV > 3.5 & EnergyEV < 5)

spectras <- bind_rows(
    vacuum,
    ch3oh,
    ch3oh_5s,
    ch3oh_10s
) %>%
    mutate(nSolvent = factor(nSolvent, levels=c("Vacuum", "0", "5", "10")))

s10 <- vacuum %>%
    filter(State == "S10")


gaussian_func <- function(A, sigma, x_0, y_0, x){
    A * exp(-(x_0 - x)**2 / (2*sigma**2)) + y_0}

gaussian <- purrr::partial(gaussian_func, A = 0.1, sigma = .068795, x_0 = 4.3, y_0 = 0.2)
x = seq(4.1, 4.5, by=0.01)
gaussian_plot <- data.frame("x" = x, "y" = gaussian(x), "State" = "S10")


ggplot() +
    geom_line(data = spectras, aes(x=EnergyEV,y=value, color = nSolvent), size=1.5) +
    geom_line(data = gaussian_plot, aes(x = x, y = y), size = 1.5) +
    geom_segment(data = s10, aes(x = 4.3, y = 0.3, xend = 4.3, yend = 0.2), size=1.2, arrow = arrow(length = unit(0.5, "cm"))) +
    facet_wrap(vars(State), ncol=1) +
    labs(x = "Energy (eV)", y = "Intensity", color = "Number QM Solvents") +
    theme_bw() +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20),
          strip.text.x=element_text(size=20),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.title = element_text(size=20),
          legend.text = element_text(size = 18),
          legend.text.align = 0,
          legend.position = c(0.85, 0.9)
    )
ggsave("~/potentialparadox.github.io/Paper2/Images/pulse_pump/spectra.png", width = 12, height = 12)
