pop_size = as.data.frame(pop_size)


with(parms, as.list(h1 + o1))
parms

## two strains
pop_size =
  pop_size %>%
  as.data.frame() %>%
  mutate(
    per_H = with(as.list(parms), (h1 * o1 * P1H)/H + (h2 * o2 * P2H)/H - (c1 * b1 * P1 + c2 * b2 * P2) - d),
    per_P1H = with(as.list(parms), (b1 * P1 * H)/P1H + DL * (e1H * psi1 * a1 * S) - (o1 + m1)),
    per_P2H = with(as.list(parms), (b2 * P2 * H)/P2H + DL * (e2H * psi2 * a2 * S) - (o2 + m2)),
    per_P1 = with(as.list(parms), (e1 * a1) * S - (b1 * H) + (1 - DL) * (e1H * psi1 * a1 * P1H * S)/P1 - m1),
    per_P2 = with(as.list(parms), (e2 * a2) * S - (b2 * H) + (1 - DL) * (e2H * psi2 * a2 * P2H * S)/P2 - m2),
    per_S = with(as.list(parms), r * (1-S/K) - (a1 * P1 + a2 * P2 + psi1 * a1 * P1H + psi2 * a2 * P2H)))
#per_H, per_P1H, per_P2H , per_P1, per_P2, per_S

## single strain
pop_size =
  pop_size %>%
  as.data.frame() %>%
  mutate(
    per_H = with(as.list(parms), (h1 * o1 * P1H)/H - (c1 * b1 * P1) - d),
    per_P1H = with(as.list(parms), (b1 * P1 * H)/P1H + DL * (e1H * psi1 * a1 * S) - (o1 + m1)),
    per_P1 = with(as.list(parms), (e1 * a1) * S - (b1 * H) + (1 - DL) * (e1H * psi1 * a1 * P1H * S)/P1 - m1),
    per_S = with(as.list(parms), r * (1-S/K) - (a1 * P1 + psi1 * a1 * P1H)))

#--------------------
with(as.list(c(pop_size, parms)), {
  per_H = (h1 * o1 * P1H)/H + (h2 * o2 * P2H)/H - (c1 * b1 * P1 + c2 * b2 * P2) - d
  per_P1H = (b1 * P1 * H)/P1H + DL * (e1H * psi1 * a1 * S) - (o1 + m1)
  per_P2H = (b2 * P2 * H)/P2H + DL * (e2H * psi2 * a2 * S) - (o2 + m2)
  per_P1 = (e1 * a1) * S - (b1 * H) + (1 - DL) * (e1H * psi1 * a1 * P1H * S)/P1 - m1
  per_P2 = (e2 * a2) * S - (b2 * H) + (1 - DL) * (e2H * psi2 * a2 * P2H * S)/P2 - m2
  per_S = r * (1-S/K) - (a1 * P1 + a2 * P2 + psi1 * a1 * P1H + psi2 * a2 * P2H)
  return(list(c(per_H, per_P1H, per_P2H, per_P1, per_P2, per_S)))
})

pop_size
with(as.list(parms), (h1 * o1 * P1H)/H + (h2 * o2 * P2H)/H - (c1 * b1 * P1 + c2 * b2 * P2) - d))
