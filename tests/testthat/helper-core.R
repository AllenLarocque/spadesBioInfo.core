# Make the project library reachable from PSOCK workers.
#
# iCAMP::pdist.big() (used by subsetPhyloDist) spawns PSOCK workers, which are fresh R
# processes that read ./.Rprofile from THEIR startup working directory. testthat::test_dir()
# runs with the working directory set to THIS folder, not the project root, so the workers
# never see the project .Rprofile, fall back to the system library, and die with
# "there is no package called 'bigmemory'". Exporting the parent's .libPaths() through R_LIBS
# lets workers inherit it via the environment, independent of working directory.
#
# This is a TEST-HARNESS fix, not a production one. Verified: iCAMP::pdist.big() does NOT call
# setwd(), so in a real pipeline run the workers start at the project root, read .Rprofile and
# resolve bigmemory normally -- which is why full runs of calcPhyloDiversity have succeeded.
# Note R_LIBS_USER points at ~/R/..., NOT the standAlone library, so it cannot serve this role.
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))

# Source the package R/ so tests run without installing the package.
local({
  rDir <- normalizePath(file.path("..", "..", "R"), mustWork = FALSE)
  if (dir.exists(rDir))
    for (f in list.files(rDir, pattern = "\\.R$", full.names = TRUE))
      sys.source(f, envir = globalenv())
})


# ---------------------------------------------------------------------------
# Response-routing fixtures, copied (NOT moved) from
# modules/fitResponseModels/tests/testthat/helper-fitResponseModels.R when the routing
# functions moved into core. The module copy STAYS -- fitMLModels' own tests use it.
# ---------------------------------------------------------------------------

#' A core table (nCore rows, 5 treatments x 3 blocks) and a "sequenced" subset of
#' nSeq of them. When biasLevel is given, the OMITTED cores are drawn only from
#' that treatment — reproducing the real defect, where Control never lost a core
#' and every other treatment did (see docs/FINDINGS.md 9.6.1).
#'
#' Every no-op path is an ERROR, not a silent pass: an unmatched biasLevel, a
#' biasLevel with too few cores to drop, or a subset that is not exactly nSeq
#' long all stop. `x[-integer(0)]` returns character(0) in R, so an injection
#' that matched nothing would otherwise "succeed" while dropping everything.
#'
#' `blockEffect` is load-bearing, not decoration. With it zero, `(1 | block)`
#' fits SINGULAR (block variance 0), lmer collapses to OLS cell means, and a
#' contrast between two levels that lost no cores is then mathematically
#' independent of the dropped rows — full and filtered come out identical to
#' 6e-15 and any "the fits differ" assertion fails. A real block effect makes
#' the variance component non-degenerate, which is the channel through which
#' treatment-biased loss actually perturbs every fixed effect.
.coreFixture <- function(nCore = 100, nSeq = 90, biasLevel = "ST", seed = 1L,
                         blockEffect = c(B1 = 0, B2 = 1.5, B3 = -1.5)) {
  stopifnot(nCore >= 1L, nSeq >= 1L, nSeq <= nCore)
  nDrop <- nCore - nSeq
  withr::with_seed(seed, {
    treatment <- rep(c("Control", "60", "30", "ST", "CC"), length.out = nCore)
    block     <- rep(names(blockEffect), length.out = nCore)
    plot      <- paste0(treatment, "_", block)
    # a real treatment effect on p_C so attenuation is detectable, plus a real
    # block effect so the random term is not estimated at zero
    p_C <- 10 + 2 * (treatment == "60") + unname(blockEffect[block]) +
      stats::rnorm(nCore, 0, 1)
    core <- data.frame(unique_name = sprintf("C%03d", seq_len(nCore)),
                       treatment = treatment, block = block, plot = plot,
                       p_C = p_C, stringsAsFactors = FALSE)
    biased <- !is.null(biasLevel) && length(biasLevel) == 1L && !is.na(biasLevel)
    drop <- if (biased) {
      inLevel <- which(core$treatment == biasLevel)
      if (length(inLevel) < nDrop)
        stop(".coreFixture: treatment '", biasLevel, "' has ", length(inLevel),
             " cores, too few to drop ", nDrop, " from.", call. = FALSE)
      utils::head(inLevel, nDrop)
    } else {
      sample(nCore, nDrop)
    }
    sequenced <- if (length(drop)) core$unique_name[-drop] else core$unique_name
    stopifnot(length(sequenced) == nSeq, !anyDuplicated(sequenced))
    list(core = core, sequenced = sequenced)
  })
}
