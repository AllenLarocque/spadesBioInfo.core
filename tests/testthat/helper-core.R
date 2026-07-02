# Source the package R/ so tests run without installing the package.
local({
  rDir <- normalizePath(file.path("..", "..", "R"), mustWork = FALSE)
  if (dir.exists(rDir))
    for (f in list.files(rDir, pattern = "\\.R$", full.names = TRUE))
      sys.source(f, envir = globalenv())
})
