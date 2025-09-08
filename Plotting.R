library(tidyverse)
library(paletteer)
#library(patchwork)
##Read the simulation result

comp_out = readRDS("Pre4A1")

comp_out = 
  readRDS("Pre4A1") %>%
  filter(round(psi1, 5) == 1, round(psi2, 5) == 1) #%>%
  #filter(round(m2, 5) == 2*m1)
comp_out[comp_out[,"Stable_E"] == "",]$Stable_E = "U"
#comp_out[which(comp_out$Stable_E == ""), "Stable_E"] = "U"

#comp_out = mutate(comp_out, total = P1+P2+P1H+P2H)
#comp_out = filter(comp_out, r<5)
####Theme setting----
A = 
theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.title.y.right = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.text.align = 0)
theme_set(A)

unique_outcomes = unique(comp_out$Outcome)
mycolor = c("TTTTT" = "#BA6338FF",#AC
            "TFTTT" = "#F0E685FF",#C/P1H
            "TTFTT" = "#CC9900FF",#C/P2H
            "FFFTT" = "#CE3D32FF",#C/H
            "TTFTF" = "#466983FF",#P1+H
            "FFFTF" = "#0A47FFFF",#P1
            "TFTFT" = "#749B58FF",#P2+H
            "FFFFT" = "#64DD17",#P2
            "FFFFF" = "#5e0084"#Unstable
)
unspecified_outcomes <- setdiff(unique_outcomes, names(mycolor))
extra_colors <- paletteer_d("ggsci::default_igv")[1:length(unspecified_outcomes)]
final_colors <- c(mycolor, setNames(extra_colors, unspecified_outcomes))

outcome_labels <- c(
  "TTTTT" = "Coexistence",
  "TFTTT" = expression(P[1/H] ~ "excluded"),
  "TTFTT" = expression(P[2/H] ~ "excluded"),
  "FFFTT" = expression(P[1] + P[2] ~ "win"),
  "TTFTF" = expression(P[1] + P[1/H] ~ "win"),
  "FFFTF" = expression(P[1] ~ "wins"),
  "TFTFT" = expression(P[2] + P[2/H] ~ "win"),
  "FFFFT" = expression(P[2] ~ "wins"),
  "FFFFF" = "Unstable"
)

#for ASS
unique_outcomes = unique(comp_out$Stable_E)
mycolor = c("C" = "#BA6338FF",
            "C,P1H" = "#80665DFF",
            "C,P2H" = "#977F48FF",
            "P1H,P2H" = "#5D826D",
            "P1,P2" = "#37928B",
            "P1H" = "#466983FF",
            "P1" = "#0A47FFFF",
            "P2H" = "#749B58FF",
            "P2" = "#64DD17",
            "S" = "#E5E5E5",
            "U" = "#5e0084"#Unstable
)
unspecified_outcomes <- setdiff(unique_outcomes, names(mycolor))
extra_colors <- paletteer_d("ggsci::default_igv")[1:length(unspecified_outcomes)]
final_colors <- c(mycolor, setNames(extra_colors, unspecified_outcomes))

outcome_labels <- c(
  "C" = "Coexistence",
  "C,P1H" = expression(C ~"or"~ P[1] + P[1/H]),
  "C,P2H" = expression(C ~"or"~ P[2] + P[2/H]),
  "P1H,P2H" = expression(P[1] + P[1/H] ~"or"~ P[2] + P[2/H]),
  "P1,P2" = expression(P[1] ~"or"~P[2]),
  "P1H" = expression(P[1] + P[1/H]),
  "P1" = expression(P[1]),
  "P2H" = expression(P[2] + P[2/H]),
  "P2" = expression(P[2]),
  "S" = "S",
  "U" = "Unstable"
)

####Plot the result----

ggplot(comp_out, aes(x = a2, y = b2, fill = Stable_E)) +
  geom_tile() +
  #geom_point(filter(comp_out), mapping = aes(x = a1, y = b1, shape = Stability), color = "black", alpha = 0.2)+
  labs(x = expression(α[2]), y = expression(β[2]))+ #title = "δ = 0 (No vertical transmission)",  ~","~ r == 1.5
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 0.8, by = 0.05)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.8, by = 0.05)) +
  scale_fill_manual(values = final_colors, labels = outcome_labels) +
  #scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+#automatically choose color
  #scale_fill_manual(values = setNames(paletteer_d("ggsci::default_igv")[1:length(all_comb)], all_comb))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  coord_fixed(ratio = 1)


####Plot the result (H as constant)----
#filter(comp_out, a1 > 0, b1 > 0)
unique_outcomes = unique(comp_out$Outcome2)
mycolor = c("TTTT" = "#BA6338FF",#AC
            "FTTT" = "#F0E685FF",#C/P1H
            "TFTT" = "#CC9900FF",#C/P2H
            "TFTF" = "#466983FF",#P1+H
            "FTFT" = "#749B58FF"#P2+H
)
unspecified_outcomes <- setdiff(unique_outcomes, names(mycolor))
extra_colors <- paletteer_d("ggsci::default_igv")[1:length(unspecified_outcomes)]
final_colors <- c(mycolor, setNames(extra_colors, unspecified_outcomes))

outcome_labels <- c(
  "TTTT" = "Coexistence",
  "FTTT" = expression(P[1/H] ~ "excluded"),
  "TFTT" = expression(P[2/H] ~ "excluded"),
  "TFTF" = expression(P[1] + P[1/H] ~ "win"),
  "FTFT" = expression(P[2] + P[2/H] ~ "win")
)

ggplot(comp_out, aes(x = a1, y = b1, z = Outcome, fill = Outcome2)) +
  geom_tile() +
  labs(title = expression(m[1] < m[2] ~","~ "real"), x = expression(α[1]), y = expression(β[1])) + #title = "δ = 0 (No vertical transmission)",  ~","~ r == 1.5
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = final_colors, labels = outcome_labels) +
  #scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+#automatically choose color
  #scale_fill_manual(values = setNames(paletteer_d("ggsci::default_igv")[1:length(all_comb)], all_comb))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  coord_fixed(ratio = 1)


############Plot the heatmap for K------------
#filter(filter(comp_out, round(a1, 3) >= 0.5, round(b1, 3) >= 0.45)
ggplot(comp_out, aes(x = a1, y = K, z = Outcome, fill = Outcome2)) +
  geom_tile() +
  geom_point(filter(comp_out, Cycle == "T"), mapping = aes(x = a1, y = K, shape = Cycle), color = "black", alpha = 0.2)+
  labs(title = expression(β[1] == 0.2), x = expression(α[1]), y = "K")+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0.35, 0.5, by = 0.05)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(10, 100, by = 30)) +
  scale_fill_manual(values = final_colors, labels = outcome_labels) +
  #scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+#automatically choose color
  #scale_fill_manual(values = setNames(paletteer_d("ggsci::default_igv")[1:length(all_comb)], all_comb))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+
  coord_fixed(ratio = 1)


############Plot the heatmap for r------------
ggplot(comp_out, aes(x = a1, y = r, z = Outcome, fill = Outcome2)) +
  geom_tile() +
  geom_point(filter(comp_out, Cycle == "T"), mapping = aes(x = a1, y = r, shape = Cycle), color = "black", alpha = 0.2)+
  labs(title = expression(β[1] == 0.2), x = expression(α[1]), y = "r")+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0.35, 0.5, by = 0.05)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(1, 4, by = 1)) +
  scale_fill_manual(values = final_colors, labels = outcome_labels) +
  #scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+#automatically choose color
  #scale_fill_manual(values = setNames(paletteer_d("ggsci::default_igv")[1:length(all_comb)], all_comb))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))

############Plot the heatmap for m------------
ggplot(filter(comp_out), aes(x = m1, y = m2, z = Stable_E, fill = Stable_E)) +
  geom_raster() +
  #geom_point(filter(comp_out, Stability == "ASS"), mapping = aes(x = m1, y = m2, shape = Stability), color = "black", alpha = 0.2)+
  #geom_point(mapping = aes(x = m1, y = m2), color = "black", alpha = 0.1)+
  labs(title = expression(psi[1] == 1 ~","~ psi[2] == 1), x = expression(m[1]), y = expression(m[2]))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 0.1, by = 0.01)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.1, by = 0.01)) +
  scale_fill_manual(values = final_colors, labels = outcome_labels) +
  #scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+#automatically choose color
  #scale_fill_manual(values = setNames(paletteer_d("ggsci::default_igv")[1:length(all_comb)], all_comb))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+
  coord_fixed(ratio = 1)

ggsave("Heatmap of m1 and m2.png", width = 16, height = 11, units = "cm", dpi = 1600)


ggplot(comp_out, aes(x = factor(m1), y = m2, z = Outcome, fill = Outcome2)) +
  geom_tile() +
  #geom_point(filter(comp_out, Cycle == "T"), mapping = aes(x = m1, y = m2, shape = Cycle), color = "black", alpha = 0.2)+
  labs(title = expression(), x = expression(m[1]), y = expression(m[2]))+
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0.01, 0.1, by = 0.01)) +
  scale_fill_manual(values = final_colors, labels = outcome_labels) +
  #scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+#automatically choose color
  #scale_fill_manual(values = setNames(paletteer_d("ggsci::default_igv")[1:length(all_comb)], all_comb))+
  #theme_minimal()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))



############Plot the heatmap for m + IGR------------
ggplot(filter(comp_out)) +
  geom_raster(mapping = aes(x = m1, y = m2, fill = Stable_E)) +
  geom_contour(filter(comp_out, m2 < 0.057), mapping = aes(x = m1, y = m2, z = IGR2), breaks = 0, color = "grey", linewidth = 1.2)+
  geom_contour(filter(comp_out, m2 < 0.057), mapping = aes(x = m1, y = m2, z = IGR1), breaks = 0, color = "darkred", linewidth = 1.2)+
  #geom_contour(filter(comp_out, m1 > 0.05, m2 < 0.057), mapping = aes(x = m1, y = m2, z = C_State_V), breaks = 0, color = "darkblue", linewidth = 1.2)+
  
  
  geom_point(filter(comp_out, m1 > 0.05, m2 < 0.057, C_State == "Infeasible"), mapping = aes(x = m1, y = m2, shape = C_State), alpha = 0.1)+
  #geom_point(filter(comp_out, IGRH2 > 0), mapping = aes(x = m1, y = m2), alpha = 0.8)+
  #geom_point(filter(comp_out, IGRH1 > 0), mapping = aes(x = m1, y = m2), alpha = 0.8)+
  labs(title = expression(psi[1] == 1 ~","~ psi[2] == 1), x = expression(m[1]), y = expression(m[2]))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 0.1, by = 0.01)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.1, by = 0.01)) +
  scale_fill_manual(values = final_colors, labels = outcome_labels) +
  #scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+#automatically choose color
  #scale_fill_manual(values = setNames(paletteer_d("ggsci::default_igv")[1:length(all_comb)], all_comb))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+
  coord_fixed(ratio = 1)

#ggsave("Heatmap of m1 and m2.png", width = 16, height = 11, units = "cm", dpi = 1600)


ggplot(comp_out, aes(x = factor(m1), y = m2, z = Outcome, fill = Outcome2)) +
  geom_tile() +
  #geom_point(filter(comp_out, Cycle == "T"), mapping = aes(x = m1, y = m2, shape = Cycle), color = "black", alpha = 0.2)+
  labs(title = expression(), x = expression(m[1]), y = expression(m[2]))+
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0.01, 0.1, by = 0.01)) +
  scale_fill_manual(values = final_colors, labels = outcome_labels) +
  #scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+#automatically choose color
  #scale_fill_manual(values = setNames(paletteer_d("ggsci::default_igv")[1:length(all_comb)], all_comb))+
  #theme_minimal()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))



#############Bifurcation diagram--------------
comp_out %>%
  select(c(a1, b1, P1H, P2H, P1, P2)) %>%
  #mutate(Total = P1H + P2H + P1 + P2) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(a1, b1)) %>%
  #gather(key = Species, value = Abundance, -c(a1, b1)) %>% #using gather()
  filter(b1 == 0.2) %>%
  ggplot(aes(x = a1, y = Abundance, color = Species)) +
  geom_line(lwd = 1) + 
  labs(title = expression(β[1] == 0.2), x = expression(α[1]), y = "Abundance", color = "Species")+
  scale_colour_manual(labels = 
                        c("P1" = expression(P[1]), "P1H" = expression(P[1/H]),
                          "P2" = expression(P[2]), "P2H" = expression(P[2/H]),
                          "S" = "Host", "H" = "Hyper", "Total" = "Total"),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", 
                                 "S" = "black", "H" = "red", "Total" = "blue"))


#############Bifurcation for r-------------
comp_out %>%
  select(c(a1, r, P1H, P2H, P1, P2, S, H)) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(a1, r)) %>%
  #gather(key = Species, value = Abundance, -c(a1, r)) %>% #using gather()
  filter(round(a1, 3) == 0.45) %>%
  ggplot(aes(x = r, y = Abundance, color = Species)) +
  geom_line(lwd = 1) + 
  labs(title = expression(β[1] == 0.2 ~","~ α[1] == 0.45), x = "r", y = "Abundance", color = "Species")+
  scale_x_continuous(breaks = seq(1, 22, by = 7)) +
  scale_colour_manual(labels = 
                        c("P1" = expression(P[1]), "P1H" = expression(P[1/H]),
                          "P2" = expression(P[2]), "P2H" = expression(P[2/H]),
                          "S" = "Host", "H" = "Hyper"),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", 
                                 "S" = "black", "H" = "red"))

#############Bifurcation for K-------------
comp_out %>%
  select(c(a1, K, P1H, P2H, P1, P2, H, S)) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(a1, K)) %>%
  #gather(key = Species, value = Abundance, -c(a1, r)) %>% #using gather()
  filter(round(a1, 3) == 0.5) %>%
  ggplot(aes(x = K, y = Abundance, color = Species)) +
  geom_line(lwd = 1) + 
  labs(title = expression(β[1] == 0.2 ~","~ α[1] == 0.5), x = "K", y = "Abundance", color = "Species")+
  scale_y_continuous() +
  scale_x_continuous(breaks = seq(10, 100, by = 30)) +
  scale_colour_manual(labels = 
                        c("P1" = expression(P[1]), "P1H" = expression(P[1/H]),
                          "P2" = expression(P[2]), "P2H" = expression(P[2/H]),
                          "S" = "Host", "H" = "Hyper"),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", 
                                 "S" = "black", "H" = "red"))


#############Bifurcation for h-------------
D = 
comp_out %>%
  select(c(r, h2, P1H, P2H, P1, P2, S, H)) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(r, h2)) %>%
  #gather(key = Species, value = Abundance, -c(a1, r)) %>% #using gather()
  filter(round(r, 3) == 1.5)

  ggplot(filter(D, Species != "H"), aes(x = h2, y = Abundance, color = Species)) +
  geom_line(filter(D, Species == "H"), mapping = aes(x = h2, y = Abundance/50, color = Species), lwd = 1)+
  geom_line(lwd = 1) +
  labs(title = expression(r == 1.5 ~","~ β[1] == 0.2 ~","~ α[1] == 0.35), x = expression(h[2]), y = "Abundance", color = "Species")+
  scale_y_continuous(sec.axis = sec_axis(~.*50, name = "Abundance of H")) +
  scale_x_continuous() +
  scale_colour_manual(labels = 
                        c("P1" = expression(P[1]), "P1H" = expression(P[1/H]),
                          "P2" = expression(P[2]), "P2H" = expression(P[2/H]),
                          "S" = "Host", "H" = "Hyper"),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", 
                                 "S" = "black", "H" = "red")) +
    
  theme(axis.title.y.right = element_text(angle = 90))
    
#############Bifurcation for h-------------
D = 
comp_out %>%
  select(c(r, h2, P1H, P2H, P1, P2, S, H)) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(r, h2)) %>%
  #gather(key = Species, value = Abundance, -c(a1, r)) %>% #using gather()
  filter(round(r, 3) == 1.5)

  ggplot(filter(D, Species != "H"), aes(x = h2, y = Abundance, color = Species)) +
  geom_line(filter(D, Species == "H"), mapping = aes(x = h2, y = Abundance/50, color = Species), lwd = 1)+
  geom_line(lwd = 1) +
  labs(title = expression(r == 1.5 ~","~ β[1] == 0.2 ~","~ α[1] == 0.35), x = expression(h[2]), y = "Abundance", color = "Species")+
  scale_y_continuous(sec.axis = sec_axis(~.*50, name = "Abundance of H")) +
  scale_x_continuous() +
  scale_colour_manual(labels = 
                        c("P1" = expression(P[1]), "P1H" = expression(P[1/H]),
                          "P2" = expression(P[2]), "P2H" = expression(P[2/H]),
                          "S" = "Host", "H" = "Hyper"),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", 
                                 "S" = "black", "H" = "red")) +
    
  theme(axis.title.y.right = element_text(angle = 90))

  
#############Bifurcation for m-------------
D = 
    comp_out %>%
    select(c(m1, m2, P1, P2, P1H, P2H, H, S)) %>% #P1, P2, P1H, P2H, H, S
    #filter(round(m2, 5) == 0.05) %>%
    pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1, m2))
    #gather(key = Species, value = Abundance, -c(a1, r)) %>% #using gather()

ggplot(D, aes(x = m2, y = Abundance, color = Species)) +
  geom_line(filter(D, Species != "total"), mapping = aes(x = m2, y = Abundance, color = Species), lwd = 1) +
  #geom_line(filter(D, Species == "total"), mapping = aes(x = m1, y = Abundance, color = Species), lwd = 0.8, linetype = 2) +
    #scale_linetype_manual(values = c("Stable" = "solid", "Unstable" = "dashed")) +
  labs(title = expression(m[2] == 2*m[1]), x = expression(m[2]), y = "Abundance", color = "Species")+
  scale_y_continuous() +
  scale_x_continuous() + #breaks = c(seq(0.2, 1, by = 0.2))
  scale_colour_manual(labels = 
                        c("P1" = expression(P[1]), "P1H" = expression(P[1/H]),
                          "P2" = expression(P[2]), "P2H" = expression(P[2/H]),
                          "S" = "Host", "H" = "Hyper"),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E",
                                  "P2" = "#B0BEC5", "P2H" = "#546E7A", 
                                  "S" = "#00AF66", "H" = "#C03728",
                                  "total" = "black")) +
    #scale_shape_manual(values = c("Stable" = 16, "Unstable" = 3)) + 
  theme(axis.title.y.right = element_text(angle = 90))
  
D = 
  comp_out %>%
  select(c(m1, m2, P1, P2, P1H, P2H, H, S)) %>% #P1, P2, P1H, P2H, H, S
  filter(round(m1, 5) == 0.05, round(m2, 5) > 0.05, round(m2, 5) < 0.08) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1, m2))
#gather(key = Species, value = Abundance, -c(a1, r)) %>% #using gather()

ggplot(D, aes(x = m2, y = Abundance, color = Species)) +
  geom_line(filter(D, Species != "total"), mapping = aes(x = m2, y = Abundance, color = Species), lwd = 1) +
  #geom_line(filter(D, Species == "total"), mapping = aes(x = m1, y = Abundance, color = Species), lwd = 0.8, linetype = 2) +
  #scale_linetype_manual(values = c("Stable" = "solid", "Unstable" = "dashed")) +
  labs(title = expression(m[1] == 0.05), x = expression(m[2]), y = "Abundance", color = "Species")+
  scale_y_continuous() +
  scale_x_continuous() + #breaks = c(seq(0.2, 1, by = 0.2))
  scale_colour_manual(labels = 
                        c("P1" = expression(P[1]), "P1H" = expression(P[1/H]),
                          "P2" = expression(P[2]), "P2H" = expression(P[2/H]),
                          "S" = "Host", "H" = "Hyper"),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", 
                                 "S" = "#00AF66", "H" = "#C03728",
                                 "total" = "black")) +
  #scale_shape_manual(values = c("Stable" = 16, "Unstable" = 3)) + 
  theme(axis.title.y.right = element_text(angle = 90))

#ggsave("All+total in m2 0057 to 76.png", width = 15, height = 11, units = "cm", dpi = 800)
  