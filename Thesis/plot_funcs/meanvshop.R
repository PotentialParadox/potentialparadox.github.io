# Inspired by the 1990 Tully Paper
library(tidyverse)
library(latex2exp)


# Create the manual points for the higher energy state S2
trial_t2s <- c(0,  1,   2, 3, 4,   5, 6,   7, 8,   9, 10)
trial_ehrenfest_e2s <- c(10, 8, 6.5, 6, 6, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75)
trial_e2s <- c(10, 8, 6.5, 6, 6, 6.25, 6.75, 7.25, 7.50, 7.75, 8.00)

trial_df2 <- data.frame('t' = trial_t2s,
                        'Tully' = trial_e2s,
                        'Ehrenfest' = trial_ehrenfest_e2s)

# Fit S2
fit2_tully <- glm(Tully ~ poly(t,3), data = trial_df2)
fit2_ehrenfest <- glm(Ehrenfest ~ poly(t,3), data = trial_df2)


predict_ts <- seq(0, 10, 0.1)
predict_t <- data.frame('t' = predict_ts)

fit2_df <- data.frame('t' = predict_ts,
                      'MF' = predict(fit2_ehrenfest, newdata = predict_t),
                      'SH' = predict(fit2_tully, newdata = predict_t)) %>%
    mutate(SH = case_when(
               t < 6 ~ MF,
               TRUE ~ SH)) %>%
    pivot_longer(-t) %>%
    rename(Simulation = name) %>%
    rename(Energy = value) %>%
    mutate(State = "S2") %>%
    mutate(PES = paste(Simulation, "S2", sep="-"))
    
#
fit2_df %>%
    ggplot(aes(x = t, y = Energy, color = PES)) +
    geom_line(size = 1.5)


# create the manaul points for the lower state S1
trial_t1s <- seq(0, 10, 1)
trial_e1s <- c(5.00, 3.00, 2.00, 3.25, 4.00, 3.50, 3.00, 2.70, 2.50, 2.25, 2.00)
#                 0,    1,    2,    3,    4,    5,    6,    7,    8,    9,  10
trial_ehrenfest_e1s <- c(5.00, 3.00, 2.00, 3.25, 4.00, 3.50, 3.00, 2.70, 2.60, 2.55, 2.50)

trial_df1 <- data.frame('t' = trial_t1s,
                        'MF' = trial_ehrenfest_e1s,
                        'SH' = trial_e1s)

fit1_tully <- glm(Tully ~ poly(t,5), data = trial_df1)
fit1_ehrenfest <- glm(Ehrenfest ~ poly(t,5), data = trial_df1)


fit1_df <- data.frame('t' = predict_ts,
                     'MF' = predict(fit1_ehrenfest, newdata = predict_t),
                     'SH' = predict(fit1_tully, newdata = predict_t)) %>%
    mutate(SH = case_when(
               t < 6.52 ~ MF,
               TRUE ~ SH)) %>%
    pivot_longer(-t) %>%
    rename(Simulation = name) %>%
    rename(Energy = value) %>%
    mutate(State = "S1") %>%
    mutate(PES = paste(Simulation, "S1", sep="-"))
    

fit1_df %>%
    ggplot(aes(x = t, y = Energy, color = Simulation)) +
    geom_line(size = 1.5)

fitted_pes <- bind_rows(fit1_df, fit2_df)

ggplot() +
    geom_line(data = fitted_pes, aes(x=t, y = Energy, color = PES), size = 1.5) +
    scale_x_continuous(limits = c(0, 10)) +
    scale_y_continuous(limits = c(0, 10))


# The probability of each state will be a logistic function from 2.5 to 5.
logistic_model <- function(ts){
    tau = 0.5
    b = 0.5
    c = 1
    b / (c + exp(-(ts-4.5) / tau))
}

logistic_df = data.frame('t' = predict_ts,
                         'P1' = logistic_model(predict_t),
                         'P2' = 1-logistic_model(predict_t)) %>%
    rename(S1 = `t.1`) %>%
    rename(S2 = `t.2`) %>%
    pivot_longer(-t) %>%
    rename(State = name) %>%
    rename(Probability = value)
                        

# We can acually print these probabilities
ggplot() +
    geom_line(data = logistic_df, aes(x = t, y = Probability, color = State), size = 1.5) +
    theme_bw() +
    scale_color_manual(breaks = c("S1", "S2"), values = c("red", "blue")) +
    labs(x = "Time", 
         y = TeX(r'($ | \Psi | ^2$)')) +
    theme(axis.text=element_blank(),
          strip.text = element_text(size=12),
          axis.title=element_text(size=20),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text.align = 0,
          legend.position = "top")
ggsave("~/potentialparadox.github.io/Thesis/Images/probabilities.png", width = 7.5, height = 5)
    
###############################################
# Plot the Ehrenfest Results
###############################################

# First we need to get the mean energy df
mean_energy <- fitted_pes %>%
    filter(Simulation == 'MF') %>%
    left_join(logistic_df, by = c("t", "State")) %>%
    mutate(weightedE = Energy * Probability) %>%
    select(t, State, weightedE) %>%
    pivot_wider(names_from = State, values_from = weightedE) %>%
    mutate(Energy = S1 + S2) %>%
    select(t, Energy) %>%
    mutate(PES = "MF") %>%
    mutate(Simulation = "Ehrenfest") %>%
    mutate(State = "Average")

# Now plot the mean enery on top of the pes plot
ehrenfest_df <- bind_rows(mean_energy, fitted_pes) %>%
    filter(PES != "SH-S1" | t > 4)
head(ehrenfest_df)

ggplot(data = ehrenfest_df, aes(x = t, y = Energy, color = PES, linetype = PES)) +
    geom_line(size = 1.5) +
    geom_segment(aes(x = 4.2, y = 5.95, xend = 4.2, yend = 3.7),
                 arrow = arrow(length = unit(0.5, "cm")),
                 color = 'red',
                 size = 1.5,
                 show.legend = FALSE) +
    scale_color_manual(breaks = c("MF-S1",
                                  "MF-S2",
                                  "MF",
                                  "SH-S1",
                                  "SH-S2"),
                       values = c("red",
                                  "blue",
                                  "black",
                                  "red",
                                  "blue")) +
    scale_linetype_manual(breaks = c("MF-S1",
                                     "MF-S2",
                                     "MF",
                                     "SH-S1",
                                     "SH-S2"),
                          values = c("dashed",
                                     "dashed",
                                     "dashed",
                                     "solid",
                                     "solid")) +
    theme_bw() +
    labs(x = "Time",
         y = "Energy"
         ) +
    theme(axis.text=element_blank(),
        axis.title=element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size =  unit(0.5, "in"),
        legend.position = "top"
    )
ggsave("~/potentialparadox.github.io/Thesis/Images/ehrenfestVsTully.png", width = 7.5, height = 5)
