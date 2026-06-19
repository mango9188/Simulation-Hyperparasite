####This bifurcation plot shows the alternative stable state
#Bifurcation for m (ASS)
# Method 1 (Old one)----

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

## Figure A---------------------------------------------
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

# values = c("P1" = "#BCAAA4", "P1H" = "#82491E", "P1T" = "#BCAAA4",
# "P2" = "#B0BEC5", "P2H" = "#546E7A", "P2T" = "#B0BEC5",
# "S" = "#00AF66", "H" = "#C03728",
# "total" = "black")


## Figure B---------------------------------------------
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


###############################################################################
# Method 2 (New one)----
comp_out =
  readRDS("Pre4A1") %>%
  filter(round(psi1, 5) == 1, round(psi2, 5) == 1)
comp_out2 =
  readRDS("Pre4A1_CLast") %>%
  filter(round(psi1, 5) == 1, round(psi2, 5) == 1)

comp_out =
  readRDS("Pre4A1_d002_psi08_for")
comp_out2 =
  readRDS("Pre4A1_d002_psi08_rev")


comp_out2 =
  readRDS("Pre4A1_d002_psi08_ForBif_Rev")
comp_out =
  readRDS("Pre4A1_d002_psi08_ForBif_For")





D = 
  comp_out %>%
  mutate(P1T = P1+P1H) %>%
  mutate(P2T = P2+P2H) %>%
  mutate(P.total = P1T+P2T) %>%
  select(c(m1, m2, P.total)) %>% #P1, P2, P1H, P2H, H, S
  #filter(round(m1, 5) == 0.06) %>%
  filter(round(m2, 5) == 0.016) %>% #0.065
  #filter(round(m1 * 1000) %% 2 == 0) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1, m2)) #SOLID EP2H

D2 = 
  comp_out2 %>%
  mutate(P1T = P1+P1H) %>%
  mutate(P2T = P2+P2H) %>%
  mutate(P.total = P1T+P2T) %>%
  select(c(m1, m2, P.total)) %>% #P1, P2, P1H, P2H, H, S
  #filter(round(m1, 5) == 0.06) %>%
  filter(round(m2, 5) == 0.016) %>%
  #filter(round(m1 * 1000) %% 2 == 0) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1, m2)) #DASH EC


D3 = 
  comp_out %>%
  select(c(m1, m2, H, S)) %>% #P1, P2, P1H, P2H, H, S
  #filter(round(m1, 5) == 0.06) %>%
  filter(round(m2, 5) == 0.016) %>%
  #filter(round(m1 * 1000) %% 2 == 0) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1, m2)) #SOLID EP2H

D4 = 
  comp_out2 %>%
  select(c(m1, m2, H, S)) %>% #P1, P2, P1H, P2H, H, S
  #filter(round(m1, 5) == 0.06) %>%
  filter(round(m2, 5) == 0.016) %>%
  #filter(round(m1 * 1000) %% 2 == 0) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1, m2)) #DASH EC

D2$Equilibrium = "Forward"
D4$Equilibrium = "Forward"
D$Equilibrium = "Reverse"
D3$Equilibrium = "Reverse"

D$Group = "1" #P1T+P2T
D2$Group = "1" #P1T+P2T
D3$Group = "2" #H+S
D4$Group = "2" #H+S

# rect_df <- data.frame(
#   xmin = c(0, 0.0335, 0.0475, 0.0745), #From (1,2,3)
#   xmax = c(0.0335, 0.0475, 0.0745, 0.1), #to (A,B,C) -> 1A, 2B, 3C -> Coexist
#   ymin = 3,
#   ymax = 3.4,
#   fill_color = c("#a50f15", "#929292", "#525252", "#2171b5"),
#   Group = "1"
# )

rect_df <- data.frame(
  xmin = c(0, 0.048, 0.0795), #From (1,2,3)
  xmax = c(0.048, 0.0795, 0.1), #to (A,B,C) -> 1A, 2B, 3C -> Coexist
  ymin = 4,
  ymax = 4.68,
  fill_color = c("#a50f15", "#525252", "#2171b5"),
  Group = "1"
)

# text_df <- data.frame(
#   x = c(0.01675, 0.0405, 0.061, 0.08725),
#   y = 3.2,
#   # geom_text -> parse = TRUE
#   label = c("E[P[1]*H]", "E[P[1]*H]~'or'~E[P[2]*H]", "E[C]~'or'~E[P[2]*H]", "E[P[2]*H]"),
#   text_color = c("white", "white", "white", "white"),
#   Group = "1"  # the text will only show on group 1
# )

text_df <- data.frame(
  x = c(0.024, 0.064, 0.08975),
  y = 4.34,
  # geom_text -> parse = TRUE
  #label = c("P[1]*H ~~ 'wins'", "'All coexist' ~~'||'~~ P[2]*H ~~ 'wins'", "P[2]*H ~~ 'wins'"),
  label = c("P[1] + H ", "P[1] + P[2] + H / P[2] + H", "P[2] + H"),
  text_color = c("white", "white", "white"),
  Group = "1"  # the text will only show on group 1
)

DataC = rbind(D, D2, D3, D4)

ggplot()+
  geom_rect(data = rect_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color),
            inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_identity() +
  geom_text(data = text_df,
            aes(x = x, y = y, label = label, color = I(text_color)), # Using I() to hide legend
            parse = TRUE, size = 4, inherit.aes = FALSE, show.legend = FALSE) +
  
  geom_point(filter(DataC), mapping = aes(x = m1, y = Abundance, color = Species))+
  labs(x = expression(P[1]*"'s mortality rate increased by pesticide"), y = "Biomass", color = "Species")+
  ylim(c(0, 4.8))+
  facet_grid(Group ~ ., scales = "fixed")+ #fixed
  #guides(linetype = FALSE)+
  scale_colour_manual(labels = State_labels, values = State_values) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y.right = element_text(angle = 90),
        legend.position = "right")
  #scale_linetype_manual(values = c("Forward" = 1, "Reverse" = 1))
  #scale_shape_manual(values = c("Forward" = 19, "Reverse" = 2))
  #scale_alpha_manual(values = c("Forward" = 1, "Reverse" = 0.5))+
  #scale_size_manual(values = c("Forward" = 1.4, "Reverse" = 2.3))
  
# 1. 建立隱形的邊界資料
# 假設你上圖的 Group 名稱是 "1"，下圖是 "2"，請依照你的 DataC$Group 實際名稱修改
blank_data <- data.frame(
  Group = c("1", "1", "2", "2"), 
  m1 = 0,                         # X 軸隨便塞一個範圍內的值即可
  Abundance = c(1, 4.8, 0, 4.8)   # 上圖要 1~4.8，下圖要 0~4.8
)

ggplot()+
  geom_rect(data = rect_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color),
            inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_identity() +
  geom_text(data = text_df,
            aes(x = x, y = y, label = label, color = I(text_color)), # Using I() to hide legend
            parse = TRUE, size = 4, inherit.aes = FALSE, show.legend = FALSE) +
  geom_point(data = filter(DataC), mapping = aes(x = m1, y = Abundance, color = Species))+
  labs(x = expression("Intrinsic mortality rate of pathogen strain 1"~ (m[1])), y = "Abundance", color = "Species")+
  
  # 2. 關鍵：加入隱形邊界圖層
  geom_blank(data = blank_data, aes(x = m1, y = Abundance), inherit.aes = FALSE) +
  
  # 3. 關鍵：改為 free_y
  facet_grid(Group ~ ., scales = "free_y") + 
  
  scale_colour_manual(labels = State_labels, values = State_values) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y.right = element_text(angle = 90),
        legend.position = "bottom",
        panel.grid.major = element_blank(), # 移除主格線
        panel.grid.minor = element_blank(), # 移除次格線
        panel.background = element_blank(), # 移除灰色背景
        panel.border = element_rect(color = "black", fill = NA))

ggsave("Multi strain bifurcation m same.png", width = 20, height = 13, units = "cm", dpi = 800)
ggsave("Poster multistrain bifur.png", width = 19, height = 13, units = "cm", dpi = 800)
