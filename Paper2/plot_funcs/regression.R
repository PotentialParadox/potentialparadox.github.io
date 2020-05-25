library(tidyverse)
library(rio)
library(mosaic)
setwd("~/Documents/paper2/populations")

tammie_fit <- function (df) {
  fitModel(
    Measured~A*exp(t/tau)/(A+exp(t/tau))-(A/(1+A)),
    data=df,
    start=list(tau=671,A=1.15),
    control=list(maxiter=50000,minFactor=1e-15)
)
}

logistic_fit <- fitModel(
  measured ~ A / (1 + exp(-(t-ti)/tau)),
  data = s1_v_time
)

fit <- tammie_fit

get_statenames <- function (df) {
  state_names <- paste("S", as.character(1:(length(df[1,])-2)), sep="")
  state_names <- append(c("t", "Sm"), state_names)
}

tidy <- function (df, solvent) {
   df %>%
    pivot_longer(-t, names_to = "State", values_to = "Measured") %>%
    mutate(Solvent = solvent) %>%
    mutate(State = as.factor(State))
}

fit_s1 <- function (df) {
  s1_v_time <- df %>% filter(State == 'S1') %>% select(Solvent, t, Measured)
  f <- fit(s1_v_time)
  s1_v_time %>% mutate(Fitted = f(t)) %>% gather("Type", "Pop", 3:4)
}

##################################
# PPV3
##################################
tammie_vacuum <- as_tibble(import("ppv3/tammie_reproduction/512-no-trivial/all_pops.txt"))
colnames(tammie_vacuum) <- get_statenames(tammie_vacuum)
tammie_vacuum <- tidy(tammie_vacuum, "Vacuum NT")

vacuum <- as_tibble(import("ppv3/vacuum/all_pops.txt"))
colnames(vacuum) <- get_statenames(vacuum)
vacuum <- tidy(vacuum, 'Vacuum')

ch3oh <- as_tibble(import("ppv3/ch3oh/all_pops.txt"))
colnames(ch3oh) <- get_statenames(ch3oh)
ch3oh <- tidy(ch3oh, 'Methanol')

ch3oh_5s <- as_tibble(import("ppv3/ch3oh_5s/all_pops.txt"))
colnames(ch3oh_5s) <- get_statenames(ch3oh_5s)
ch3oh_5s <- tidy(ch3oh_5s, 'Methanol-5s')

s1_v_time_ppv3 <- map(list(tammie_vacuum, vacuum, ch3oh, ch3oh_5s), fit_s1) %>%
  bind_rows() %>% mutate(Solute = 'PPV3')

vacuum_ppv3 <- vacuum %>% 
  mutate(Solute = 'PPV3')

##################################
# PPV3-NO2
##################################
vacuum <- as_tibble(import("ppv3-no2/vacuum/all_pops.txt"))
colnames(vacuum) <- get_statenames(vacuum)
vacuum <- tidy(vacuum, 'Vacuum')

ch3oh <- as_tibble(import("ppv3-no2/ch3oh/all_pops.txt"))
colnames(ch3oh) <- get_statenames(ch3oh)
ch3oh <- tidy(ch3oh, 'Methanol')

ch3oh_5s <- as_tibble(import("ppv3-no2/ch3oh_5s/all_pops.txt"))
colnames(ch3oh_5s) <- get_statenames(ch3oh_5s)
ch3oh_5s <- tidy(ch3oh_5s, 'Methanol-5s')

s1_v_time_ppv3_no2 <- map(list(vacuum, ch3oh, ch3oh_5s), fit_s1) %>%
  bind_rows() %>% mutate(Solute = 'PPV3-NO2')

vacuum_ppv3_no2 <- vacuum %>% 
  mutate(Solute = 'PPV3-NO2')

##################################
# Plot S1
##################################

s1_v_time <- bind_rows(
  s1_v_time_ppv3,
  s1_v_time_ppv3_no2
)

graph_title = expression(paste(PPV[3]))
graph_caption = expression(paste(S[1], "Population"))

s1_v_time %>% filter(Type == "Measured") %>% 
  filter(t <= 500) %>% 
  ggplot(aes(x=t, y=Pop)) +
  coord_cartesian(ylim = c(0,1))+
  labs(x="Time (fs)", y="Population")+
  geom_line(mapping = aes(color = Solvent), size=1)+
  facet_wrap(~ Solute) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        legend.position = "top",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)
  )

###################################
# Plot all pops in vacuum
###################################
vacuum <- bind_rows(
  vacuum_ppv3 %>% filter(State %in% c('Sm', 'S7', 'S2', 'S1')),
  vacuum_ppv3_no2 %>% filter(State %in% c('Sm', 'S9', 'S2', 'S1'))
)

all_pops <- vacuum %>% 
  filter(t <= 500)
all_pops %>% ggplot(aes(x=t, y=Measured)) +
  geom_line(mapping = aes(color = State), size = 1) +
  labs(x="Time (fs)", y="Population")+
  facet_wrap(~ Solute) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "top"
  )



###################################
# Plot fitted tenua data
###################################
vacuum <- as_tibble(read_table("ppv3/vacuum/ppv3_fitted_sm.txt", col_names = FALSE))
state_names <- c('t', 'Sm Fitted', 'S1 Fitted', 'S2 Fitted', 'SL Fitted', 'S7 Fitted', 'S8 Fitted', 'Sm', 'S1', 'S2', 'SL', 'S7', 'S8')
colnames(vacuum) <-state_names
vac <- vacuum %>% pivot_longer(-t, names_to = "Name", values_to = "Value") %>% 
  mutate(State = str_sub(Name, 1,2)) %>% 
  mutate(Type = ifelse(str_detect(Name, "Fitted"), "Fitted", "Measured")) %>% 
  select(t, State, Type, Value)

svacuum <- vac %>% spread(Type, Value)


svacuum %>% ggplot(aes(x = t, y = Fitted, color = State)) +
  geom_line() +
  geom_line(mapping = aes(y = Measured), size = 1) +
  labs(x="Time (fs)", y="Population")+
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)
  )


