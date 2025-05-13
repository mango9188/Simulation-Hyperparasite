# Load necessary library for root finding
library(nleqslv)

parms <- c(r = 1, K = 10,
           a1 = 0.5, a2 = 0.5, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.45, b2 = 0.45, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)
# Define the system of equations
equilibrium = function(vars) {
  with(parms, {
    H = vars[1]
    P1H = vars[2]
    P2H = vars[3]
    P1 = vars[4]
    P2 = vars[5]
    S = vars[6]
    c(
      (h1 * o1 * P1H) + (h2 * o2 * P2H) - (c1 * b1 * P1 + c2 * b2 * P2) * H - (d * H),
      (b1 * P1 * H) + DL * (e1H * R1 * a1 * P1H * S) - (o1 + m1) * P1H,
      (b2 * P2 * H) + DL * (e2H * R2 * a2 * P2H * S) - (o2 + m2) * P2H,
      (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * R1 * a1 * P1H * S) - m1 * P1,
      (e2 * a2 * P2) * S - (b2 * P2) * H + (1 - DL) * (e2H * R2 * a2 * P2H * S) - m2 * P2,
      r * S * (1-S/K) - (a1 * P1 + a2 * P2 + R1 * a1 * P1H + R2 * a2 * P2H) * S
    )
  })
}

# Initial guesses for x and y
start <- c(H = 0, P1H = 0, P2H = 0, P1 = 0, P2 = 0, S = 0)

# Find equilibrium points
solution <- nleqslv(start, equilibrium)

# Display results
solution$solution
