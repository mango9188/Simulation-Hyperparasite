####For ESJ
#Bifurcation for m (ASS)

comp_out = 
  readRDS("Pre4A1") %>%
  filter(round(psi1, 5) == 1, round(psi2, 5) == 1)
comp_out2 = 
  readRDS("Pre4A1_CLast") %>%
  filter(round(psi1, 5) == 1, round(psi2, 5) == 1)

D = 
  comp_out %>%
  mutate(P1T = P1+P1H) %>%
  mutate(P2T = P2+P2H) %>%
  select(c(m1, m2, P1T, P2T)) %>% #P1, P2, P1H, P2H, H, S
  filter(round(m1, 5) == 0.05) %>%
  filter(round(m2, 5) < 0.057) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1, m2))

D2 = 
  comp_out2 %>%
  mutate(P1T = P1+P1H) %>%
  mutate(P2T = P2+P2H) %>%
  select(c(m1, m2, P1T, P2T)) %>% #P1, P2, P1H, P2H, H, S
  filter(round(m1, 5) == 0.05) %>%
  filter(round(m2, 5) < 0.057) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1, m2))

####################Figure A---------------------------------------------
ggplot() +
  geom_line(D2, mapping = aes(x = m2, y = Abundance, color = Species), lwd = 1) +
  geom_line(D, mapping = aes(x = m2, y = Abundance, color = Species), lwd = 1, linetype = 2) +
  labs(x = expression(m[2]), y = "Abundance", color = "Species")+
  geom_ribbon(mapping = aes(x = c(0, 0.014), ymin = 2.8, ymax = 3), fill = "#C2185B")+#EC or EP2H
  annotate("text", x = 0.007, y = 2.9,
           label = expression(E[C]~"or"~E[P2H]), parse = TRUE, size = 3.5)+ 
  
  geom_ribbon(mapping = aes(x = c(0.015, 0.053), ymin = 2.8, ymax = 3), fill = "#4F4F4F")+ #EC
  annotate("text", x = 0.034, y = 2.9,
           label = expression(E[C]), parse = TRUE, size = 3.5, color = "white")+ 
  
  geom_ribbon(mapping = aes(x = c(0.054, 0.056), ymin = 2.8, ymax = 3), fill = "#546E7A")+ #EP2H
  annotate("text", x = 0.055, y = 2.9,
           label = expression(E[P2H]), parse = TRUE, size = 3.5, color = "white")+ 
  
  scale_y_continuous(limits = c(0, 3), breaks = c(seq(0, 2, by = 0.5))) +
  scale_x_continuous(breaks = c(seq(0.01, 0.053, by = 0.01))) +
  scale_colour_manual(labels = 
                        c("P1" = expression(P[1]), "P1H" = expression(P[1/H]), "P1T" = expression(P[1]~"total"),
                          "P2" = expression(P[2]), "P2H" = expression(P[2/H]), "P2T" = expression(P[2]~"total"),
                          "S" = "S", "H" = "H"),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E", "P1T" = "#BCAAA4",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", "P2T" = "#B0BEC5",
                                 "S" = "#00AF66", "H" = "#C03728",
                                 "total" = "black")) +
  theme(axis.title.y.right = element_text(angle = 90))

####################Figure B---------------------------------------------
D3 = 
  comp_out %>%
  select(c(m1, m2, H, S)) %>% #P1, P2, P1H, P2H, H, S
  filter(round(m1, 5) == 0.05) %>%
  filter(round(m2, 5) < 0.057) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1, m2))

D4 = 
  comp_out2 %>%
  select(c(m1, m2, H, S)) %>% #P1, P2, P1H, P2H, H, S
  filter(round(m1, 5) == 0.05) %>%
  filter(round(m2, 5) < 0.057) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1, m2))

ggplot() +
  geom_line(D3, mapping = aes(x = m2, y = Abundance, color = Species), lwd = 1) +
  geom_line(D4, mapping = aes(x = m2, y = Abundance, color = Species), lwd = 1, linetype = 2) +
  labs(x = expression(m[2]), y = "Abundance", color = "Species")+
  geom_ribbon(mapping = aes(x = c(0, 0.014), ymin = 2.8, ymax = 3), fill = "#C2185B")+#EC or EP2H
  annotate("text", x = 0.007, y = 2.9,
           label = expression(E[C]~"or"~E[P2H]), parse = TRUE, size = 3.5)+ 
  
  geom_ribbon(mapping = aes(x = c(0.015, 0.053), ymin = 2.8, ymax = 3), fill = "#4F4F4F")+ #EC
  annotate("text", x = 0.034, y = 2.9,
           label = expression(E[C]), parse = TRUE, size = 3.5, color = "white")+ 
  
  geom_ribbon(mapping = aes(x = c(0.054, 0.056), ymin = 2.8, ymax = 3), fill = "#546E7A")+ #EP2H
  annotate("text", x = 0.055, y = 2.9,
           label = expression(E[P2H]), parse = TRUE, size = 3.5, color = "white")+ 
  
  scale_y_continuous(limits = c(0, 3), breaks = c(seq(0, 2, by = 0.5))) +
  scale_x_continuous(breaks = c(seq(0.01, 0.053, by = 0.01))) +
  scale_colour_manual(labels = 
                        c("P1" = expression(P[1]), "P1H" = expression(P[1/H]), "P1T" = expression(P[1]~"total"),
                          "P2" = expression(P[2]), "P2H" = expression(P[2/H]), "P2T" = expression(P[2]~"total"),
                          "S" = "S", "H" = "H"),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E", "P1T" = "#BCAAA4",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", "P2T" = "#B0BEC5",
                                 "S" = "#00AF66", "H" = "#C03728",
                                 "total" = "black")) +
  theme(axis.title.y.right = element_text(angle = 90))

ggsave("Test.png", width = 30, height = 11.25, units = "cm", dpi = 800)

A = 
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    legend.text.align = 0,
    legend.position = "bottom")
theme_set(A)
