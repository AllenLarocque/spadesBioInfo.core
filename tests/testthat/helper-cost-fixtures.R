# Fixtures come from the model written out here, not from the fitter.
.fixMem <- function(taxa, samples, nw, a = 9.5 / (4000^2 * 200^0.5), b = 2, c = 0.5,
                    k = 0.3, w = 0.5) {
  D <- a * taxa^b * samples^c
  D * (1 + k * (nw - 1)) + w * nw
}
.fixTime <- function(taxa, samples, nw, s = 0.2) {
  1000 * (taxa / 4000)^1.5 * (samples / 200)^0.5 * (s + (1 - s) / nw)
}
.fixRow <- function(taxa, samples, nw, role, memMult = 1, timeMult = 1) {
  data.frame(step = "x", taxa = taxa, samples = samples, nworker = nw, role = role,
             marginalGiB = .fixMem(taxa, samples, nw) * memMult,
             elapsedSec = .fixTime(taxa, samples, nw) * timeMult,
             stringsAsFactors = FALSE)
}
.fixRows <- function(interaction = TRUE, memMult = 1, timeMult = 1) {
  r <- rbind(
    .fixRow(4000, 200, 1, "axis:taxa,axis:samples,axis:nworker"),
    .fixRow(c(1000, 2000), 200, 1, "axis:taxa"),
    .fixRow(4000, c(50, 100), 1, "axis:samples"),
    .fixRow(4000, 200, c(2, 4, 8), "axis:nworker"),
    .fixRow(8000, 400, 4, "corner", memMult, timeMult))
  if (interaction)
    r <- rbind(r, .fixRow(c(1000, 1000, 4000), 200, c(1, 8, 8), "interaction"),
               .fixRow(4000, 200, 1, "interaction"))
  r
}
.fixNote <- function(memory = "", time = "", bd = "", bv = NA_real_) {
  data.frame(step = "x", backend = "", memory = memory, time = time, boundaryDriver = bd,
             boundaryValue = bv, stringsAsFactors = FALSE)
}
near <- function(x, y, tol = 0.05) abs(x - y) <= tol * abs(y)

