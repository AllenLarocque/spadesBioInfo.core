#' The canonical project plot style: the single source of truth for palettes + theme, consumed by
#' plotting modules ONLY through the accessors below (never by reaching into the list). See the
#' plot-style-foundation design spec.
#' @return list(theme, treatmentColours, qualitative, diverging, treatmentGradient, referenceLine).
#' @export
projectPlotStyle <- function() {
  list(
    theme = ggplot2::theme_classic(base_size = 11) +
            ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"),
                           legend.position = "right"),
    # Exact thesis hues, ordered by disturbance intensity (treatment_continuous): green (intact) ->
    # brown (clearcut). Name-keyed so any subset of levels gets its correct hue.
    treatmentColours = c(Control = "darkolivegreen", ST = "darkolivegreen3",
                         `60` = "tan", `30` = "tan3", CC = "chocolate4"),
    qualitative = c("#E69F00", "#56B4E9", "#009E73", "#F0E442",   # Okabe-Ito, colourblind-safe
                    "#0072B2", "#D55E00", "#CC79A7", "#000000"),
    diverging   = "Blue-Red 3",                    # colorspace diverging: blue = -, red = +
    treatmentGradient = c("darkolivegreen", "chocolate4"),        # continuous green->brown endpoints
    # Disturbance encoding (treatment_continuous). MUST match prepCoreMetadata::deriveTreatment's
    # contMap — lets plots place categorical treatments on the continuous disturbance axis.
    treatmentPositions = c(Control = 0, ST = 10, `60` = 40, `30` = 70, CC = 100),
    referenceLine = "grey60"
  )
}

#' The shared ggplot2 theme.
#' @export
reviewTheme <- function(style) style$theme

#' Name-keyed treatment colours for `levels` (looked up in the disturbance-ordered master vector);
#' a level absent from the master falls back to a qualitative colour + a warning (never errors).
#' @export
treatmentPalette <- function(style, levels) {
  levels <- as.character(levels)
  cols   <- style$treatmentColours[levels]
  miss   <- is.na(cols)
  if (any(miss)) {
    warning("treatmentPalette: level(s) not in the treatment palette: ",
            paste(levels[miss], collapse = ", "), " — using qualitative fallback.")
    cols[miss] <- qualitativePalette(style, levels[miss])
  }
  stats::setNames(unname(cols), levels)
}

#' Colourblind-safe qualitative colours for arbitrary `levels` (Okabe-Ito, position-recycled).
#' @export
qualitativePalette <- function(style, levels) {
  levels <- as.character(levels)
  stats::setNames(rep(style$qualitative, length.out = length(levels)), levels)
}

#' `n` diverging colours (colorspace Blue-Red 3): blue = negative, red = positive.
#' @export
divergingPalette <- function(style, n) colorspace::diverging_hcl(n, palette = style$diverging)

#' The colour for zero/reference lines and other neutral chrome.
#' @export
referenceLineColour <- function(style) style$referenceLine

#' The disturbance-axis position (treatment_continuous encoding) of each treatment level.
#' @export
treatmentPositions <- function(style) style$treatmentPositions
