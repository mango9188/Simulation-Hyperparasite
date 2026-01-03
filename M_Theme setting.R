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
# mycolor = c("C" = "#BA6338FF",
#             "C,P1H" = "#80665DFF",
#             "C,P2H" = "#977F48FF",
#             "P1H,P2H" = "#5D826D",
#             "P1,P2" = "#37928B",
#             "P1H" = "#466983FF",
#             "P1" = "#0A47FFFF",
#             "P2H" = "#749B58FF",
#             "P2" = "#64DD17",
#             "S" = "#E5E5E5",
#             "U" = "#5e0084"#Unstable
# )
# mycolor <- c(
#   "C"        = "#4D4D4D",
#   "C,P1H"    = "#B38B2D",
#   "C,P2H"    = "#2F5F8F",
#   "P1H,P2H"  = "#5E5E5E",
#   "P1,P2"    = "#B89C4A",
#   "P1H"      = "#E69F00",
#   "P1"       = "#F3C97A",
#   "P2H"      = "#377EB8",
#   "P2"       = "#A6C8E3",
#   "S"        = "#1B9E77",
#   "U"        = "#BDBDBD"
# )
# 
# mycolor <- c(
#   "C"        = "grey50",
#   "C,P1H"    = "#B38B2D",
#   "C,P2H"    = "#2F5F8F",
#   "P1H"      = "#654321", #dark brown
#   "P1"       = "#895129", #brown
#   "P2H"      = "lightblue",
#   "P2"       = "blue",
#   "S"        = "darkgreen",
#   "U"        = "pink"
# )
mycolor <- c(
  "C" = "#4F4F4F",
  "U" = "#C7C7C7",
  "S" = "#00AF66",
  "P1" = "#BCAAA4",
  "P1H" = "#82491E",
  "P2" = "#B0BEC5",
  "P2H" = "#546E7A",
  "P1,P2" = "#B7AEA0",
  "P1H,P2H" = "#6A1B9A",
  "C,P1H" = "pink",
  "C,P2H" = "#C2185B")


unspecified_outcomes <- setdiff(unique_outcomes, names(mycolor))
extra_colors <- paletteer_d("ggsci::default_igv")[1:length(unspecified_outcomes)]
final_colors <- c(mycolor, setNames(extra_colors, unspecified_outcomes))

outcome_labels <- c(
  "C" = expression(E[C]),
  "C,P1H" = expression(E[C] ~"or"~ E[P1H]),
  "C,P2H" = expression(E[C] ~"or"~ E[P2H]),
  "P1H,P2H" = expression(E[P1H] ~"or"~ E[P2H]),
  "P1,P2" = expression(E[P1] ~"or"~ E[P2]),
  "P1H" = expression(E[P1H]),
  "P1" = expression(E[P1]),
  "P2H" = expression(E[P2H]),
  "P2" = expression(E[P2]),
  "S" = expression(E[S]),
  "U" = "Unstable"
)
