library(tidyverse)
parms <- list(
           r = 1, K = 10,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

### Leslie matrix for H and PiH's invasion growth rate
Leslie <- matrix(
  data = with(parms, {
    c(
    -(c1 * b1 * P1 + d), h1 * o1,
    b1 * P1, -(m1 + o1)
    )
  }),
  nrow = 2,
  ncol = 2,
  byrow = T)

EIGEN <- eigen(Leslie)
### Dominant eigenvalue -> invasion growth rate
abs(EIGEN$values[1])

####Analytic sol of eigenvalue
EIGEN_A = with(parms, {
  P1 = (r/a1)*(1-(m1)/(e1*a1*K))
  A = 1
  B = (c1*b1*P1 + d + m1 + o1)
  C = -h1*o1*b1*P1 + (c1*b1*P1 +d) * (m1 + o1)
  eigen = (-B + sqrt(B^2-4*A*C))/2*A
  return(eigen)
})

IGR = 
with(parms,{
  P1 = (r/a1)*(1-(m1)/(e1*a1*K))
  A = b1*P1 -c1*b1*P1 - d + h1*o1 - m1*o1
  B = -(h1*o1 - m1*o1) - (c1*b1*P1 + d) - h1 * o1
  C = h1 * o1
  phiH = (-B + sqrt(B^2 - 4*A*C))/2*A 
  IGR = (h1*o1 - m1*o1) * (1-phiH) + (b1*P1- c1*b1*P1 - d) * phiH #it is possible that both roots are positive.
  return(IGR)
})
