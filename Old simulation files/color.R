generate_combinations <- function() {
  # Generate all combinations using expand.grid
  combinations <- expand.grid(rep(list(c("T", "F")), 5))
  
  # Combine columns into a single character string
  result <- apply(combinations, 1, paste0, collapse = "")
  
  # Return the result as a vector
  return(result)
}
generate_combinations()

all_comb = expand.grid(rep(list(c("T", "F")),5)) %>%
  apply(1, paste0, collapse = "")
all_comb

paletteer_d("ggsci::default_igv")

mycolor = c("TTTTT" = "#BA6338FF",#AC
            "TFTTT" = "#F0E685FF",#C/P1H
            "TTFTT" = "#CC9900FF",#C/P2H
            "FFFTT" = "#CE3D32FF",#C/H
            "TTFTF" =  "#466983FF",#P1+H
            "FFFTF" = "#5050FFFF",#P1
            "TFTFT" = "#749B58FF",#P2+H
            "FFFFT" = "#33CC00FF"#P2
            )

mycolor = c("TTTTT" = "#A2DA5A",#AC
            "TFTTT" = "#BA0000",#C/P1H
            "TTFTT" = "#C49A00FF",#C/P2H
            "FFFTT" = "#3fb5a6",#C/H
            "TTFTF" = "#E6770B",#P1+H
            "FFFTF" = "#F7BB2A",#P1
            "TFTFT" = "#8246AF",#P2+H
            "FFFFT" = "#B37FDBFF"#P2
)
mycolor = c(all_comb = c(paletteer_d("ggsci::default_igv")))

