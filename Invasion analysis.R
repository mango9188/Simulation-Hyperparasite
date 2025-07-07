library(tidyverse)
parms <- list(
  r = 1, K = 10,
  a1 = 0.35, psi1 = 1, e1 = 0.5,
  b1 = 0.2, m1 = 0.05, e1H = 0.5,
  o1 = 0.8, h1 = 1,  c1 = 0.9, d = 0.03, DL = 0)

#For H and Pi/H's to invade E_SPi----
### Leslie matrix for H and PiH's invasion growth rate-
Leslie <- matrix(
  data = with(parms, {
    P1 = (r/a1)*(1-(m1/(e1*a1*K)))
    M = 
      c(-(c1 * b1 * P1 + d), h1 * o1, b1 * P1, -(m1 + o1))
    return(M)
  }),
  nrow = 2,
  ncol = 2,
  byrow = T)

EIGEN <- eigen(Leslie)
### Dominant eigenvalue -> invasion growth rate
EIGEN$values

####Analytic sol of eigenvalue
with(parms, {
  P1 = (r/a1)*(1-(m1/(e1*a1*K)))
  A = 1
  # B = (c1*b1*P1 + d + m1*o1)
  B = (c1*b1*P1 + d + m1+o1)
  # C = (c1*b1*P1 + d) * (m1*o1) - h1*o1*b1*P1
  C = (c1*b1*P1 + d) * (m1+o1) - h1*o1*b1*P1
  eigen1 = (-B + sqrt(B^2-4*A*C)) / (2*A)
  eigen2 = (-B - sqrt(B^2-4*A*C)) / (2*A)
  return(c(eigen1, eigen2))
})


with(parms,{
  P1 = (r/a1)*(1-(m1/(e1*a1*K)))
  # A = b1*P1 -c1*b1*P1 - d + h1*o1 - m1*o1
  A = -b1*P1 +c1*b1*P1 + d + h1*o1 - m1-o1
  # B = -(h1*o1 - m1*o1) - (c1*b1*P1 + d) - h1 * o1
  B = -(h1*o1 - m1-o1) - (c1*b1*P1 + d) - h1 * o1
  C = h1 * o1
  phiH1 = (-B + sqrt(B^2 - 4*A*C)) / (2*A)
  phiH2 = (-B - sqrt(B^2 - 4*A*C)) / (2*A)
  #note that phiH must be greater than zero and smaller than 1.
  IGR1 = 
    ifelse(phiH1 > 0, 
           ifelse(phiH1 < 1, 
                  # (h1*o1 - m1*o1) * (1-phiH1) + (b1*P1- c1*b1*P1 - d) * phiH1,
                  (h1*o1 - m1-o1) * (1-phiH1) + (b1*P1- c1*b1*P1 - d) * phiH1,
                  "Not meaningful"),
           "Not meaningful")
  IGR2 =
    ifelse(phiH2 > 0,
          ifelse(phiH2 < 1,
                # (h1*o1 - m1*o1) * (1-phiH2) + (b1*P1- c1*b1*P1 - d) * phiH2,
                  (h1*o1 - m1-o1) * (1-phiH2) + (b1*P1- c1*b1*P1 - d) * phiH2,
                  "Not meaningful"),
          "Not meaningful")
  print(c(IGR1, IGR2))
  #print(c(phiH1, phiH2))
})

###Approach2
with(parms,{
  P1 = (r/a1)*(1-(m1/(e1*a1*K)))
  A = h1*o1
  B = m1+o1 - c1*b1*P1 - d
  C = -b1*P1
  Z1 = (-B + sqrt(B^2 - 4*A*C)) / (2*A)
  Z2 = (-B - sqrt(B^2 - 4*A*C)) / (2*A)
  #note that phiH must be greater than zero
  #it is possible that both roots are positive.
  IGR1 = 
    ifelse(Z1 > 0,
           (h1*o1*Z1 - c1*b1*P1 - d),
           "Not meaningful")
  IGR2 =
    ifelse(Z2 > 0,
           (h1*o1*Z2 - c1*b1*P1 - d),
           "Not meaningful")
  print(c(IGR1, IGR2))
  #print(c(Z1, Z2))
})




#For Pi and Pi/H to invade E_SPiPi/H----
### Leslie matrix for Pi and Pi/H's invasion growth rate
Leslie <- matrix(
  data = with(parms, {
    H =
    S = 
    M = 
      c(-(c1 * b1 * P1 + d), h1 * o1, b1 * P1, -(m1 + o1))
    return(M)
  }),
  nrow = 2,
  ncol = 2,
  byrow = T)

EIGEN <- eigen(Leslie)
### Dominant eigenvalue -> invasion growth rate
EIGEN$values