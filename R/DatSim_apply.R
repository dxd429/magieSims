### Simulate data for one iteration
DatSim_apply <- function(nreps,
                         nsites,
                         mu,
                         phi,
                         theta,
                         sfx,
                         sfy,
                         res.para,
                         Contrast,
                         sd_multi,
                         seed) {
  res.sim <- SimDat(
    nreps = nreps,
    nsites = nsites,
    mu = mu,
    phi = phi,
    theta = theta,
    sfx = sfx,
    sfy = sfy,
    sd_multi,
    seed = seed
  )
  return(res.sim)
}
