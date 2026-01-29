library(paletteer)
#Theme setting----
A = 
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.title.y.right = element_text(size = 15),
    legend.text = element_text(size = 11),
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
#   "P2H"      = "#9ED8F9",
#   "P2"       = "#36B8FA",
#   "S"        = "darkgreen",
#   "U"        = "pink"
# )
# mycolor <- c(
#   "C" = "#C73824", #red
#   "U" = "#C7C7C7", #This will be replaced by other grids
#   "S" = "#65AD62", #Green
#   "P1" = "#663300",
#   "P1H" = "#A37649",
#   "P2"  = "#5D278F",
#   "P2H" = "#8E6FCC",
#   "P1H,P2H" = "#8E0703", #dark red
#   "C,P1H" = "#D6824E", #orange
#   "C,P2H" = "#D66883" #pink
#   )

# mycolor <- c(
#   "C" = "#4F4F4F",
#   "U" = "#C7C7C7",
#   "S" = "#00AF66",
#   "P1" = "#BCAAA4",
#   "P1H" = "#82491E",
#   "P2" = "#B0BEC5",
#   "P2H" = "#546E7A",
#   "P1,P2" = "#B7AEA0",
#   "P1H,P2H" = "#6A1B9A",
#   "C,P1H" = "pink",
#   "C,P2H" = "#C2185B")

# mycolor <- c(
#   "C" = "#5D4D5C",
#   "U" = "#C7C7C7",
#   "S" = "#00AF66",
#   "P1" = "#F99BAF",
#   "P1H" = "#c60c37",
#   "P2" = "#68E1E1",
#   "P2H" = "#0074AE",
#   "P1,P2" = "#B7AEA0",
#   "P1H,P2H" = "#6A1B9A",
#   "C,P1H" = "pink",
#   "C,P2H" = "#C155B8")

mycolor <- c(
  "C" = "#BB7DBE", #BB7DBE #C155B8
  "U" = "#5D4D5C",
  "S" = "#00AF66",
  "P1" = alpha("#a50f15", 0.4),
  "P1H" = "#a50f15",
  "P2" = "#9ecae1",
  "P2H" = "#2171b5",
  "P1,P2" = "#B7AEA0",
  "P1H,P2H" = "#929292",
  "C,P1H" = "pink",
  "C,P2H" = "#525252")

"#D6B701" #H
"#00AF66" #S

unspecified_outcomes <- setdiff(unique_outcomes, names(mycolor))
extra_colors <- paletteer_d("ggsci::default_igv")[1:length(unspecified_outcomes)]
final_colors <- c(mycolor, setNames(extra_colors, unspecified_outcomes))

outcome_labels <- c(
  "C" = expression(E[C]),
  "C,P1H" = expression(E[C] ~"or"~ E[P[1]*H]),
  "C,P2H" = expression(E[C] ~"or"~ E[P[2]*H]),
  "P1H,P2H" = expression(E[P[1]*H] ~"or"~ E[P[2]*H]),
  "P1,P2" = expression(E[P[1]] ~"or"~ E[P[2]]),
  "P1H" = expression(E[P[1]*H]),
  "P1" = expression(E[P[1]]),
  "P2H" = expression(E[P[2]*H]),
  "P2" = expression(E[P[2]]),
  "S" = expression(E[S]),
  "U" = "Unstable"
)

###For bifurcation and time series----
State_labels = c("H" = "H",
                 "P1" = expression(P[1]),
                 "P1H" = expression(P[1/H]),
                 "P2" = expression(P[2]),
                 "P2H" = expression(P[2/H]),
                 "S" = "S")
            
State_values = c("H" = "#D6B701", 
                 "P1" = alpha("#a50f15", 0.4),
                 "P1H" = "#a50f15",
                 "P2" = "#9ecae1",
                 "P2H" = "#2171b5",
                 "S" = "#00AF66")

# State_labels = c("P1" = "Pathogen", "P1H" = "Hyperparasited Pathogen", "S" = "Host", "H" = "Hyperparasite", "P.total" = "Total Pathogen")
# 
# State_values = c("P1" = "#BCAAA4", "P1H" = "#82491E", "S" = "#00AF66", "H" = "#C03728", "P.total" = "black")
