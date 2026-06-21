# ==============================================================================
# 10_RiskOfBiasFigures.R
# Generates the two ROBINS-E risk-of-bias figures (traffic-light + summary)
# from "Risk of bias assessment - Sheet1.csv" using the robvis package.
#
# The assessment uses a custom 4-level judgement vocabulary:
#   Low  <  Some concern  <  Concern  <  High
# mapped to green / yellow / red / dark-red respectively.
#
# Run from the Code/ directory (paths are relative to it), e.g.:
#   Rscript 10_RiskOfBiasFigures.R
# ==============================================================================

suppressPackageStartupMessages({
  library(robvis)
  library(dplyr)
})

# ----- Locate input/output relative to this script (Code/) --------------------
rob_csv  <- "../Risk of bias assessment - Sheet1.csv"
fig_dir  <- "../Figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ----- Read and reshape the assessment ---------------------------------------
# The CSV stores each cell as "<Label>\n<score>"; we keep only the label and
# coerce it to the judgement vocabulary expected by the figures.
raw <- read.csv(rob_csv, check.names = FALSE, stringsAsFactors = FALSE)

# First token (before the newline) of each cell is the judgement label.
first_line <- function(x) trimws(sub("\n.*$", "", x))

domain_cols <- grep("^Domain", names(raw), value = TRUE)
stopifnot(length(domain_cols) == 7L)

rob <- data.frame(
  Study = first_line(raw[["Author (year)"]]),
  stringsAsFactors = FALSE
)
for (i in seq_along(domain_cols)) {
  rob[[paste0("D", i)]] <- first_line(raw[[domain_cols[i]]])
}
rob[["Overall"]] <- first_line(raw[["Overall risk score"]])

# ----- Map the study's custom vocabulary onto the robvis ROBINS-E template ----
# The ROBINS-E template internally uses five judgement levels
#   (Low, Some concerns, High, Very high, No information).
# This assessment uses a four-level scheme (Low < Some concern < Concern < High)
# which the original figures expressed by relabelling the template levels:
#   Low          -> "Low"            (green)
#   Some concern -> "Some concerns"  (yellow)
#   Concern      -> "High"           (red)       [legend shown as "Concern"]
#   High         -> "Very high"      (dark red)  [legend shown as "High"]
to_template <- c(
  "Low"          = "Low",
  "Some concern" = "Some concerns",
  "Concern"      = "High",
  "High"         = "Very high"
)
map_label <- function(x) {
  x <- gsub("Some concerns", "Some concern", x, fixed = TRUE)  # guard plurals
  unname(to_template[x])
}
rob[-1] <- lapply(rob[-1], map_label)
stopifnot(!anyNA(unlist(rob[-1])))

# ----- Colours for the 5 template levels (5th = "No information", unused) -----
rob_colours <- c(
  "#02C100",  # Low            -> green   (+)
  "#E2DF07",  # Some concerns  -> yellow  (-)
  "#BF0000",  # High           -> red     (x)
  "#820000",  # Very high      -> dark red(!)
  "#CCCCCC"   # No information  (unused here)
)

# Relabel legend text back to the study's vocabulary WITHOUT replacing the
# fill scale (replacing it wipes robvis's internal colour mapping). We edit the
# existing discrete fill scale's labels in place via a relabelling function.
# Map both robvis's full judgement labels (traffic-light) AND its single-letter
# internal codes (summary plot: l/s/h/v/n) to the study's display vocabulary.
legend_relabel <- c(
  "Low"            = "Low",
  "Some concerns"  = "Some concern",
  "High"           = "Concern",
  "Very high"      = "High",
  "No information" = "No information",
  "l" = "Low",
  "s" = "Some concern",
  "h" = "Concern",
  "v" = "High",
  "n" = "No information"
)
# Applies to fill (summary bars) and colour/shape (traffic-light points), which
# is where robvis stores the judgement legend labels for the two figure types.
relabel_legend <- function(p) {
  aes_to_fix <- c("fill", "colour", "color", "shape")
  p$scales$scales <- lapply(p$scales$scales, function(sc) {
    if (!is.null(sc$aesthetics) && any(sc$aesthetics %in% aes_to_fix)) {
      sc$labels <- function(breaks) {
        out <- unname(legend_relabel[as.character(breaks)])
        out[is.na(out)] <- as.character(breaks)[is.na(out)]
        out
      }
    }
    sc
  })
  p
}

# ----- Traffic-light plot -----------------------------------------------------
tl <- rob_traffic_light(
  data    = rob,
  tool    = "ROBINS-E",
  colour  = rob_colours,
  psize   = 13          # larger glyphs / in-cell symbols
)
tl <- relabel_legend(tl)
# In the traffic-light, the in-table text is drawn as facet-strip text:
#   strip.text.x       -> column headers (D1-D7, Overall)
#   strip.text.y.left  -> study names (author (year))
# and the two axis titles ("Risk of bias domains", "Study") come from
# axis.title. Enlarge the strip text and drop the axis titles.
tl <- tl + ggplot2::theme(
  text              = ggplot2::element_text(size = 18),
  axis.title        = ggplot2::element_blank(),  # remove "Study" / "Risk of bias domains"
  strip.text.x      = ggplot2::element_text(size = 18),
  strip.text.y.left = ggplot2::element_text(size = 18, angle = 0),
  legend.text       = ggplot2::element_text(size = 16),
  legend.title      = ggplot2::element_blank(),  # drop the "judgement" legend title
  # Push the "Domains:" caption block down so it clears the bottom table row.
  plot.caption      = ggplot2::element_text(
    size = 14,
    margin = ggplot2::margin(t = 22, unit = "pt")
  )
)
ggplot2::ggsave(
  filename = file.path(fig_dir, "risk_of_bias_traffic_light.png"),
  plot = tl, width = 12, height = 8, dpi = 150, bg = "white"
)

# ----- Summary (stacked bar) plot --------------------------------------------
sm <- rob_summary(
  data     = rob,
  tool     = "ROBINS-E",
  colour   = rob_colours,
  overall  = TRUE,
  weighted = FALSE
)
sm <- relabel_legend(sm)
# Use short domain codes (D1-D7, Overall) for the domain-axis row labels.
# Domains are on the y aesthetic; map by full caption (name-based, order-safe).
domain_short <- c(
  "Bias due to confounding"                                                 = "D1",
  "Bias arising from measurement of the exposure"                           = "D2",
  "Bias in selection of participants into the study (or into the analysis)" = "D3",
  "Bias due to post-exposure interventions"                                 = "D4",
  "Bias due to missing data"                                                = "D5",
  "Bias arising from measurement of the outcome"                            = "D6",
  "Bias in selection of the reported result"                                = "D7",
  "Overall risk of bias"                                                    = "Overall"
)
# Domains are on the x aesthetic (the plot is flipped via coord_flip, so they
# appear as the rows). Named-vector labels match by value (order-safe).
# After the flip: axis.text.x -> the domain row labels (D1-D7, Overall);
#                 axis.text.y -> the percent column labels.
sm <- sm +
  ggplot2::scale_x_discrete(labels = domain_short) +
  ggplot2::theme(
  text         = ggplot2::element_text(size = 18),
  axis.text.x  = ggplot2::element_text(size = 22),  # D1-D7 / Overall row labels
  axis.text.y  = ggplot2::element_text(size = 9),   # 0%-100% column labels
  axis.title   = ggplot2::element_text(size = 18),
  legend.text  = ggplot2::element_text(size = 16),
  legend.title = ggplot2::element_blank()      # drop the "judgement" legend title
)
ggplot2::ggsave(
  filename = file.path(fig_dir, "risk_of_bias_summary.png"),
  plot = sm, width = 11, height = 6, dpi = 150, bg = "white"
)

cat("Risk-of-bias figures written to", normalizePath(fig_dir), "\n")
cat("Studies included:", nrow(rob), "\n")
print(rob[, c("Study", "Overall")])
