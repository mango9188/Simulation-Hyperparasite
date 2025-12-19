####LV model for project
library(deSolve)
library(tidyverse)
library(patchwork)

times <- seq(0, 3000, by = 0.1)

#state
ini_High_H <- c(H = 0.5, P1H = 0, P2H = 0, P1 = 0.01, P2 = 0.01, S = 0.5) #ini_1 (P1win)
ini_Low_H <- c(H = 0.01, P1H = 0, P2H = 0, P1 = 0.01, P2 = 0.01, S = 0.5) #ini_2 (P2win)

#parms
parms_E1H_E2H = 
  c(epsilon = 1,
    r = 1, K = 10,
    a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
    b1 = 0.2, b2 = 0.45, m1 = 0.03, m2 = 0.001, e1H = 0.5, e2H = 0.5,
    o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)
parms_EC_E2H = 
  c(epsilon = 1,
    r = 1, K = 10,
    a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
    b1 = 0.2, b2 = 0.45, m1 = 0.045, m2 = 0.001, e1H = 0.5, e2H = 0.5,
    o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)


TimeSeries = function(times, state, parms){
  M2 <- function(times, state, parms) {
    with(as.list(c(state, parms)), {
      dH_dt = (h1 * o1 * P1H) + (h2 * o2 * P2H) - (c1 * b1 * P1 + c2 * b2 * P2) * H - (d * H)
      dP1H_dt = (b1 * P1 * H) + DL * (e1H * psi1 * a1 * P1H * S) - (o1 + m1) * P1H
      dP2H_dt = (b2 * P2 * H) + DL * (e2H * psi2 * a2 * P2H * S) - (o2 + m2) * P2H
      dP1_dt = (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * psi1 * a1 * P1H * S) - m1 * P1
      dP2_dt = (e2 * a2 * P2) * S - (b2 * P2) * H + (1 - DL) * (e2H * psi2 * a2 * P2H * S) - m2 * P2
      dS_dt = (r * S * (1-S/K) - (a1 * P1 + a2 * P2 + psi1 * a1 * P1H + psi2 * a2 * P2H) * S)
      return(list(c(dH_dt, dP1H_dt, dP2H_dt, dP1_dt, dP2_dt, dS_dt)))
    })
  }
  pop_size = ode(func = M2, times = times, y = state, parms = parms)
  pop_size %>%
    as.data.frame() #%>%
    # filter(time %% 1 == 0) %>%
    # pivot_longer(cols = c("H", "P1H", "P2H", "P1", "P2", "S"),
    #              names_to = "species", values_to = "biomass") %>%
    # ggplot(mapping = aes(x = time, y = biomass, color = species)) +
    # labs(x = "Time", y = "Abundance") +
    # geom_line(lwd = 1) +
    # scale_y_continuous(limits = c(0, 8))+
    # scale_colour_manual(labels = c("H" = "H", "P1" = expression(P[1]), "P1H" = expression(P[1/H]), "P2" = expression(P[2]), "P2H" = expression(P[2/H]), "S" = "S"),
    #                     values = c("H" = "#C03728", "P1" = "#BCAAA4", "P1H" = "#82491E",
    #                                "P2" = "#B0BEC5", "P2H" = "#546E7A", "S" = "#00AF66"))
}

High_H_E1H_E2H = TimeSeries(times, ini_High_H, parms_E1H_E2H)
Low_H_E1H_E2H = TimeSeries(times, ini_Low_H, parms_E1H_E2H)
High_H_EC_E2H = TimeSeries(times, ini_High_H, parms_EC_E2H)
Low_H_EC_E2H = TimeSeries(times, ini_Low_H, parms_EC_E2H)

High_H_E1H_E2H$ASS = "E1H or E2H"
Low_H_E1H_E2H$ASS = "E1H or E2H"
High_H_EC_E2H$ASS = "EC or E2H"
Low_H_EC_E2H$ASS = "EC or E2H"

High_H_E1H_E2H$m1 = "Low m1"
Low_H_E1H_E2H$m1 = "Low m1"
High_H_EC_E2H$m1 = "High m1"
Low_H_EC_E2H$m1 = "High m1"

High_H_E1H_E2H$ini = "High initial value of H"
Low_H_E1H_E2H$ini = "Low initial value of H"
High_H_EC_E2H$ini = "High initial value of H"
Low_H_EC_E2H$ini = "Low initial value of H"

D = 
  do.call("rbind", list(High_H_E1H_E2H, Low_H_E1H_E2H, High_H_EC_E2H, Low_H_EC_E2H)) %>%
  ##or using bind_rows()
  pivot_longer(cols = c("H", "P1H", "P2H", "P1", "P2", "S"),
               names_to = "species", values_to = "biomass")

#D$m1 = factor(D$m1, levels = c("High m1", "Low m1"))
#D$ini = factor(D$ini, levels = c("High initial value of H", "Low initial value of H"))

TimePlot =
D %>%
  ggplot(data = D, mapping = aes(x = time, y = biomass, color = species))+
  geom_line(lwd = 0.8) +
  labs(x = "Time", y = "Abundance")+
  scale_y_continuous(limits = c(0, 8))+
  scale_colour_manual("Species", labels = c("H" = "H", "P1" = expression(P[1]), "P1H" = expression(P[1/H]), "P2" = expression(P[2]), "P2H" = expression(P[2/H]), "S" = "S"),
                      values = c("H" = "#C03728", "P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", "S" = "#00AF66"))+
  facet_grid(ini ~ fct_rev(m1))

heatmap = heatmap + labs(tag = "(A)") + theme(plot.tag.position = c(0.05, 0.99))


TimePlot = TimePlot + labs(tag = "(B)") + theme(plot.tag.position = c(0, 0.99))

heatmap + plot_spacer() + TimePlot + 
  plot_layout(widths = c(5, 0.2, 8)) & 
  theme(plot.tag = element_text(size = 20, face = "bold"))

heatmap + plot_spacer() + TimePlot + plot_layout(widths = c(5, 1, 10))

ggsave("ESJ Heatmap and ASS time series.png", width = 30, height = 11.25, units = "cm", dpi = 800)
