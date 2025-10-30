#Theme setting----
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
  "C" = expression(E[C]),
  "C,P1H" = expression(E[C] ~"or"~ E[1~H]),
  "C,P2H" = expression(E[C] ~"or"~ E[2~H]),
  "P1H,P2H" = expression(E[1~H] ~"or"~ E[2~H]),
  "P1,P2" = expression(E[1] ~"or"~ E[2]),
  "P1H" = expression(E[1~H]),
  "P1" = expression(E[1]),
  "P2H" = expression(E[2~H]),
  "P2" = expression(E[2]),
  "S" = expression(E[S]),
  "U" = "Unstable"
)