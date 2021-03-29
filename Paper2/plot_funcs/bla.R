library(tidyverse)
library(latex2exp)
setwd('~/Documents/paper2/bla')



################################
# Data Collection
################################
collect_blas <- function (runtime, infile, state, location, solvent){
  data <- read_csv(infile, col_names = FALSE)
  nsteps <- ncol(data) - 1
  time_labels <- (0:(nsteps-1))*runtime/(nsteps -1)
  names(data) <- c('Bond', time_labels)
  s1_bond_data <- data %>% gather(key = "t", value = "Length", -Bond) %>% 
    group_by(Bond, t) %>% 
    summarise(MeanLength = mean(Length)) %>% 
    spread(Bond, MeanLength) %>% 
    mutate(BLA = (`1`+`3`)/2 - `2`) %>% 
    mutate(t = as.numeric(t)) %>% 
    mutate(State = state) %>% 
    mutate(Location = location) %>% 
    mutate(SolventID = solvent) %>% 
    mutate(Solvent = if_else(location == '', solvent, sprintf("%s-%s", solvent, location))) %>% 
    mutate(Description = if_else(location == '', state, sprintf("%s-%s", state, location)))
  s1_bond_data
}


vacuum_near <- collect_blas(runtime = 1,
                            infile = 'bla-near-vacuum.csv',
                            state = 'S1',
                            location = "Near",
                            solvent = "Vacuum")

ch3oh_near <- collect_blas(runtime = 1,
                            infile = 'bla-near-ch3oh.csv',
                            state = 'S1',
                            location = "Near",
                            solvent = "0")

ch3oh_5s_near <- collect_blas(runtime = 1,
                              infile = 'bla-near-ch3oh-5s.csv',
                              state = 'S1',
                              location = "Near",
                              solvent = "5")

ch3oh_10s_near <- collect_blas(runtime = 1,
                               infile = 'bla-near-ch3oh-10s.csv',
                               state = 'S1',
                               location = "Near",
                               solvent = "10")


bla_near <- bind_rows(
    vacuum_near,
    ch3oh_near,
    ch3oh_5s_near,
    ch3oh_10s_near
)


vacuum_far <- collect_blas(runtime = 1,
                           infile = 'bla-far-vacuum.csv',
                           state = 'S1',
                           location = "Far",
                           solvent = "Vacuum")

ch3oh_far <- collect_blas(runtime = 1,
                          infile = 'bla-far-ch3oh.csv',
                          state = 'S1',
                          location = "Far",
                          solvent = "0")

ch3oh_far_5s <- collect_blas(runtime = 1,
                          infile = 'bla-far-ch3oh-5s.csv',
                          state = 'S1',
                          location = "Far",
                          solvent = "5")

ch3oh_far_10s <- collect_blas(runtime = 1,
                          infile = 'bla-far-ch3oh-10s.csv',
                          state = 'S1',
                          location = "Far",
                          solvent = "10")

bla_far <- bind_rows(
    vacuum_far,
    ch3oh_far,
    ch3oh_far_5s,
    ch3oh_far_10s
)


bla <- bind_rows(
    bla_near,
    bla_far
) %>%
    mutate(Location = factor(Location, levels = c("Near", "Far"))) %>%
    mutate(SolventID = factor(SolventID, levels = c("Vacuum", "0", "5", "10")))

######################################
# Ploting
######################################

bla %>%
    ggplot(aes(x = t, y = BLA, color = SolventID)) +
    geom_line(size=1.5) +
    scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.10, by=0.02)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    labs(x = "Time (ps)",
         y = expression(paste("BLA (", ring(A), ")" )),
         color = "Number QM Solvents") +
    facet_wrap(~ Location) +
    theme_bw() +
    theme(axis.text=element_text(size=15),
          strip.text = element_text(size=12),
          axis.title=element_text(size=20),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.text.align = 0,
          legend.position = "top")
ggsave("~/potentialparadox.github.io/Paper2/Images/bla/solvent_comparison.png", width = 7.5, height = 5)


bondorders %>%
    ggplot(aes(x = Timefs, y = MeanBondOrder, color = NSolvent)) +
    geom_line(size = 1.5) +
    facet_wrap(~ Bond) +
    annotation_custom(diagram, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

