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
