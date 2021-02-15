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

model_tammie <- function(params, data){
    A = params[1]
    tau = params[2]
    (A * exp(data$t / tau)) / (A + exp(data$t / tau)) - A / (1 + A)
}


measure_distance <- function(mod_params, data) {
    diff <- data$value - model_tammie(mod_params, data)
    sqrt(mean(diff^2))
}


# Tammies
best <- optim(c(1.62, 100), measure_distance, data = s1)
s1_fit <- purrr::partial(model_tammie, c(1.53, 119.46))
## fit_df <- data.frame("x" = s1$t, "value" = s1_fit(s1))

# Modified Logistic
## logistic_model <- function(params, data){
##     tau = params[1]
##     b = params[2]
##     c = params[3]
##     b / (c + exp(-data$t / tau))
## }
## best <- optim(c(119, 1, 1, 0.5), measure_distance, data = s1)
## s1_fit <- purrr::partial(logistic_model, c(119, 2, 1.5))

s1 <- ch3oh_10s %>%
    filter(State == "S1")

fit_df <- data.frame("x" = s1$t, "value" = s1_fit(s1))
s1 %>%
    ggplot(aes(x = t, y = value)) +
    geom_line(color = "blue") +
    geom_line(data = fit_df, aes(x=x, y=value))

ggsave("~/potentialparadox.github.io/Paper2/Images/populations/s1_fit.png", width = 8, height = 8)

    
