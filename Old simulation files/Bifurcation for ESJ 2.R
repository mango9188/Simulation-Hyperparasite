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
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1, m2)) #SOLID EP2H

D2 = 
  comp_out2 %>%
  mutate(P1T = P1+P1H) %>%
  mutate(P2T = P2+P2H) %>%
  select(c(m1, m2, P1T, P2T)) %>% #P1, P2, P1H, P2H, H, S
  filter(round(m1, 5) == 0.05) %>%
  filter(round(m2, 5) < 0.057) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1, m2)) #DASH EC


D3 = 
  comp_out %>%
  select(c(m1, m2, H, S)) %>% #P1, P2, P1H, P2H, H, S
  filter(round(m1, 5) == 0.05) %>%
  filter(round(m2, 5) < 0.057) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1, m2)) #SOLID EP2H

D4 = 
  comp_out2 %>%
  select(c(m1, m2, H, S)) %>% #P1, P2, P1H, P2H, H, S
  filter(round(m1, 5) == 0.05) %>%
  filter(round(m2, 5) < 0.057) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1, m2)) #DASH EC

D$Equilibrium = "EP2H"
D3$Equilibrium = "EP2H"
D2$Equilibrium = "EC"
D4$Equilibrium = "EC"

D$Group = "1" #P1T+P2T
D2$Group = "1" #P1T+P2T
D3$Group = "2" #H+S
D4$Group = "2" #H+S

rect_df <- data.frame(
  xmin = c(0, 0.0145, 0.0535),
  xmax = c(0.0145, 0.0535, 0.056),
  ymin = 2.6,
  ymax = 3,
  fill_color = c("#C2185B", "#4F4F4F", "#546E7A"),
  Group = "1"  # 關鍵：指定只出現在 Group 1 (上圖)
)

text_df <- data.frame(
  x = c(0.007, 0.034, 0.055),
  y = 2.8,
  # 這裡用字串格式寫數學式，稍後在 geom_text 開啟 parse = TRUE
  label = c("E[C]~'or'~E[P2H]", "E[C]", "E[P2H]"),
  text_color = c("black", "white", "white"),
  Group = "1"  # 關鍵：指定只出現在 Group 1 (上圖)
)


DataC = rbind(D, D2, D3, D4)

ggplot()+
  geom_rect(data = rect_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color),
            inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_identity() +
  geom_text(data = text_df,
            aes(x = x, y = y, label = label, color = I(text_color)), # 使用 I() 保持原色，不進圖例
            parse = TRUE, size = 4, inherit.aes = FALSE, show.legend = FALSE) +
  
  geom_line(DataC, mapping = aes(x = m2, y = Abundance, linetype = Equilibrium, color = Species), lwd = 1)+
  labs(x = expression(m[2]), y = "Abundance", color = "Species")+
  facet_grid(Group ~ ., scales = "fixed")+
  guides(linetype = FALSE)+
  scale_colour_manual(labels = 
                        c("P1" = expression(P[1]), "P1H" = expression(P[1/H]), "P1T" = expression(P[1]~"total"),
                          "P2" = expression(P[2]), "P2H" = expression(P[2/H]), "P2T" = expression(P[2]~"total"),
                          "S" = "S", "H" = "H"),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E", "P1T" = "#BCAAA4",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", "P2T" = "#B0BEC5",
                                 "S" = "#00AF66", "H" = "#C03728",
                                 "total" = "black")) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y.right = element_text(angle = 90))

# ggplot() +
#   geom_line(DataC, mapping = aes(x = m2, y = Abundance, color = Species, linetype = Equilibrium), lwd = 1) +
#   labs(x = expression(m[2]), y = "Abundance", color = "Species")+
#   geom_ribbon(mapping = aes(x = c(0, 0.014), ymin = 2.8, ymax = 3), fill = "#C2185B")+#EC or EP2H
#   annotate("text", x = 0.007, y = 2.9,
#            label = expression(E[C]~"or"~E[P2H]), parse = TRUE, size = 3.5)+ 
#   
#   geom_ribbon(mapping = aes(x = c(0.015, 0.053), ymin = 2.8, ymax = 3), fill = "#4F4F4F")+ #EC
#   annotate("text", x = 0.034, y = 2.9,
#            label = expression(E[C]), parse = TRUE, size = 3.5, color = "white")+ 
#   
#   geom_ribbon(mapping = aes(x = c(0.054, 0.056), ymin = 2.8, ymax = 3), fill = "#546E7A")+ #EP2H
#   annotate("text", x = 0.055, y = 2.9,
#            label = expression(E[P2H]), parse = TRUE, size = 3.5, color = "white")+ 
#   
#   scale_y_continuous(limits = c(0, 3), breaks = c(seq(0, 2, by = 0.5))) +
#   scale_x_continuous(breaks = c(seq(0.01, 0.053, by = 0.01))) +
#   scale_colour_manual(labels = 
#                         c("P1" = expression(P[1]), "P1H" = expression(P[1/H]), "P1T" = expression(P[1]~"total"),
#                           "P2" = expression(P[2]), "P2H" = expression(P[2/H]), "P2T" = expression(P[2]~"total"),
#                           "S" = "S", "H" = "H"),
#                       values = c("P1" = "#BCAAA4", "P1H" = "#82491E", "P1T" = "#BCAAA4",
#                                  "P2" = "#B0BEC5", "P2H" = "#546E7A", "P2T" = "#B0BEC5",
#                                  "S" = "#00AF66", "H" = "#C03728",
#                                  "total" = "black")) +
#   facet_grid(Group ~ ., scales = "fixed")+
#   guides(linetype = FALSE)+
#   theme(axis.title.y.right = element_text(angle = 90))
#   