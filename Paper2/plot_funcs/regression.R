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
    mutate(nQM_Solvent = factor(nQM_Solvent, levels=c("Vacuum", "0", "5", "10"))) %>%
    rename(Population = value)

ch3oh_comparison %>%
    ggplot(aes(x = t, y = Population, color = nQM_Solvent)) +
    geom_line(size=1.5) +
    facet_grid(~ State) +
    theme_bw() +
    labs(color = "Number QM Solvents",
         x = "Time (fs)"
         ) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          legend.text = element_text(size = 18),
          strip.text = element_text(size=15),
          legend.title = element_text(size=20),
          legend.text.align = 0,
          legend.position = "top"
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

model_decay <- function(params, data){
    A = params[1]
    B = params[2]
    tau1 = params[3]
    tau2 = params[4]
    B - (A*exp(- data$t / tau1) + (B-A)*exp(- data$t / tau2))
    ## 1 - (A*exp(- data$t / tau1) + B*exp(- data$t / tau2))
    ## s1_ch3oh_10s_fit <- purrr::partial(model, c(0.1571757, 0.9715770, 725.54619, 164.627))
}

model <- model_decay

measure_distance <- function(mod_params, data) {
    diff <- data$value - model(mod_params, data)
    sqrt(mean(diff^2))
}

# Tammies
s1_vacuum <- vacuum %>%
    filter(State == "S1")
s1_ch3oh <- ch3oh %>%
    filter(State == "S1")
s1_ch3oh_5s <- ch3oh_5s %>%
    filter(State == "S1")
s1_ch3oh_10s <- ch3oh_10s %>%
    filter(State == "S1")

s1 <- bind_rows(s1_vacuum, s1_ch3oh, s1_ch3oh_5s, s1_ch3oh_10s) %>%
    mutate(Type="Measured") %>%
    select(t, nQM_Solvent, Type, value)

# Finds the single exponential
test_system <- s1_ch3oh_10s %>%
    mutate(value = ifelse(value == 1, 0.9999, value))
yp = log(1 - test_system$value)
df = data.frame(x = test_system$t, y = yp)
result <- lm(formula = y ~ x, data = df)
1 / result$coefficients[[2]]

# Finds the double exponential
test_system <- s1_vacuum
best <- optim(c(0.0358, 1, 500, 150), measure_distance, data = test_system, method = "L-BFGS-B", lower = c(0, 0, 0, 0), upper=c(4,1,5000,5000))


s1_vacuum_fit <- purrr::partial(model, c(0, 1, 10, 149.20942))
fit_vacuum_df <- tibble("t" = s1_vacuum$t, nQM_Solvent = "Vacuum", Type = "Predicted", "value" = s1_vacuum_fit(s1_vacuum))

s1_ch3oh_fit <- purrr::partial(model, c(0.1491528, 0.98348, 500, 150.00))
fit_ch3oh_df <- tibble("t" = s1_ch3oh$t, nQM_Solvent = "0", Type = "Predicted", "value" = s1_ch3oh_fit(s1_ch3oh))

s1_ch3oh_5s_fit <- purrr::partial(model, c(0.0, 0.9651859, 493.34, 214.152))
fit_ch3oh_5s_df <- tibble("t" = s1_ch3oh_5s$t, nQM_Solvent = "5", Type = "Predicted", "value" = s1_ch3oh_5s_fit(s1_ch3oh_5s))

s1_ch3oh_10s_fit <- purrr::partial(model, c(0.0, 0.97, 492.119, 218.34))
fit_ch3oh_10s_df <- tibble("t" = s1_ch3oh_10s$t, nQM_Solvent = "10", Type = "Predicted", "value" = s1_ch3oh_10s_fit(s1_ch3oh_10s))



s1 <- bind_rows(s1, fit_vacuum_df, fit_ch3oh_df, fit_ch3oh_5s_df, fit_ch3oh_10s_df)

## fit_error <- sqrt(mean((s1$value - fit_df$value)^2))
s1 %>%
    mutate(nQM_Solvent = factor(nQM_Solvent, levels = c("Vacuum", "0", "5", "10"))) %>%
    filter(value >= 0) %>%
    ggplot(aes(x = t, y = value, color = Type)) +
    geom_line() +
    facet_grid(~ nQM_Solvent) +
    labs(x = "Time (fs)", y = "Population", 
         color = "Number QM Solvents") +
    theme_bw() +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          legend.title = element_blank(),
          strip.text = element_text(size=20),
          legend.text = element_text(size = 20),
          legend.text.align = 0,
          legend.position = "top"
    )
ggsave("~/potentialparadox.github.io/Paper2/Images/populations/s1_fits.png", width = 20, height = 6)

    
