install.packages("Ryacas")
library(Ryacas)
x = ysym("x")      # define x as a symbolic variable
a = ysym("a") 
b = ysym("b")
lim(x*b + 1 /(x*a + 2), x, 0) # limit of sin(x)/x as x approaches 0
