# ==============================================================================
# STEP 16 - Figure A0: the pipeline from sources to estimates
# ==============================================================================
#
# One picture of what every source contributes and in which order the steps
# use it: the sources on the left, the three preparation tracks they feed
# (combatants, civilians, and the denominator with its counterfactual), the
# parameter table and the simulation they converge on, and the results. The
# checks that consume the register's linkage but feed no estimate are drawn
# apart, dotted.
#
# The line of each box says how its output enters the simulation: dashed if it
# is drawn in every simulation, solid if it is held at one value and varied
# only in the sensitivity analyses. So the figure is also the map of what the
# intervals cover.
#
# The few numbers it carries come from the run: the register releases from
# ual_releases and the size of the simulation from n_sim (00_setup.R), so the
# figure changes when the pipeline does.
#
# INPUTS   ual_releases, n_sim (00_setup.R)
# OUTPUT   figures/figA0_pipeline.png
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

rel <- ual_releases$release
n_rel <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine")[length(rel)]

# the boxes: column (x), height (y), title, the line below it, the track it
# belongs to, and how its output enters the simulation
col_x <- c(src = 14, prep = 43, est = 72, res = 101)
bx <- tribble(
  ~id,       ~col,   ~y, ~title,                       ~sub,                                         ~track,   ~enters,
  "ual",     "src",  76, "UALosses register",          paste0(str_to_sentence(n_rel), " releases, ", rel[1], " to ", rel[length(rel)]), "source", "none",
  "off",     "src",  66, "Official counts",            "Prisoners held, returned; bodies",            "source", "none",
  "ucdp",    "src",  54, "UCDP GED 26.1",              "Events: low, best, high",                    "source", "none",
  "hrmmu",   "src",  47, "UN HRMMU",                   "Verified civilian casualties",               "source", "none",
  "acled",   "src",  40, "ACLED",                      "Comparison only",                            "check",  "check",
  "sssu",    "src",  33, "SSSU",                       "Base 2022; deaths 1989–2021",           "source", "none",
  "wpp",     "src",  24, "UN WPP 2024",                "Fertility by age",                           "source", "none",
  "stock",   "src",  16, "Stocks abroad",              "Eurostat, UNHCR, Canada, US",                "source", "none",
  "cross",   "src",   9, "Net crossings",              "CES survey, border data",                    "source", "none",
  "p08",     "prep", 77, "08 Register counts",         "Combatants by age, sex, status",             "combat", "fixed",
  "p08b",    "prep", 70, "08b Listing delay",          "Chain ladder to 48 months",                  "combat", "drawn",
  "p09",     "prep", 63, "09 The missing",             "Dead unless alive: prisoners, bound",        "combat", "drawn",
  "p07",     "prep", 49, "07 Civilian deaths",         "Totals with bounds; age profile",            "civil",  "drawn",
  "p0103",   "prep", 33, "01–03 Life tables",     "Kannisto at old ages",                       "denom",  "fixed",
  "p04",     "prep", 26, "04 Lee–Carter",          "Counterfactual and its error",               "denom",  "drawn",
  "p05",     "prep", 19, "05 Fertility",               "2023 carried forward",                       "denom",  "fixed",
  "p06",     "prep", 12, "06 Net migration",           "Two readings bracketed",                     "denom",  "drawn",
  "v09j",    "est",  77, "09j Windows held out",       "Deaths forecast; the alive not",             "check",  "check",
  "v09k",    "est",  71, "09k Cohort test",            "Cohorts differ; none carried",               "check",  "check",
  "v09g",    "est",  65, "09g Multistate model",       "Tested; not adopted",                        "check",  "check",
  "t10",     "est",  45, "10 Parameter table",         "Min, central, max per input",                "est",    "none",
  "t11",     "est",  33, "11 Simulation",              paste0(scales::comma(n_sim), " cohort-component draws"), "est", "none",
  "r13",     "res",  57, "13 Sensitivities",           "One choice at a time",                       "res",    "none",
  "r14e",    "res",  48, "14e Second interval",        "Structural choices sampled",                 "res",    "none",
  "r14",     "res",  39, "14 Decomposition",           "Loss by age and component",                  "res",    "none",
  "r14bd",   "res",  30, "14b–14d Further measures", "YLL, 45q15; two tiers",                  "res",    "none"
) |>
  mutate(x = col_x[col], w = 24)
bx <- bind_rows(bx, tibble(id = "r15", col = "res", y = 13, title = "15 Tables and figures",
                           sub = "Every number in the papers", track = "res", enters = "none",
                           x = mean(col_x[c("est", "res")]), w = 53))
h <- 5.2
bx <- bx |> mutate(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2)
b <- function(i) bx[bx$id == i, ]

# The base on 1 January 2022 enters the projection itself, not only the
# counterfactual fitted to its series; and the second interval (14e) reads
# every sensitivity step's output.
# the edges, as paths: straight from the right edge of one box to the left edge
# of the next, or down a column from the bottom of one to the top of the next
right <- function(i, dy = 0) c(b(i)$xmax, b(i)$y + dy)
left <- function(i, dy = 0) c(b(i)$xmin, b(i)$y + dy)
down <- function(from, to) rbind(c(b(from)$x, b(from)$ymin), c(b(to)$x, b(to)$ymax))
across <- function(from, to, dy_to = 0, dy_from = 0) rbind(right(from, dy_from), left(to, dy_to))
ed <- list(
  ual_p08 = across("ual", "p08"), ual_p09 = across("ual", "p09", 1.2, -1.2),
  off_p09 = across("off", "p09", -0.6),
  ucdp_p07 = across("ucdp", "p07", 1.2), hrmmu_p07 = across("hrmmu", "p07", 0),
  acled_p07 = across("acled", "p07", -1.2),
  sssu_p0103 = across("sssu", "p0103"), wpp_p05 = across("wpp", "p05"),
  stock_p06 = across("stock", "p06", 0.8), cross_p06 = across("cross", "p06", -0.8),
  p08_p08b = down("p08", "p08b"), p08b_p09 = down("p08b", "p09"), p0103_p04 = down("p0103", "p04"),
  p09_v09j = across("p09", "v09j"), p09_v09k = across("p09", "v09k"), p09_v09g = across("p09", "v09g"),
  p09_t10 = across("p09", "t10", 2.0), p07_t10 = across("p07", "t10", 1.0),
  p0103_t10 = across("p0103", "t10", 0.4), p04_t10 = across("p04", "t10", -0.2),
  p05_t10 = across("p05", "t10", -1.1), p06_t10 = across("p06", "t10", -2.0),
  t10_t11 = down("t10", "t11"),
  t10_r13 = across("t10", "r13"), r13_r14e = down("r13", "r14e"),
  t11_r14e = across("t11", "r14e", -1.2, 1.2), t11_r14 = across("t11", "r14"), t11_r14bd = across("t11", "r14bd", 0, -1.2)
)
check_edges <- c("acled_p07", "p09_v09j", "p09_v09k", "p09_v09g")
edges <- imap_dfr(ed, function(m, nm) tibble(edge = nm, px = m[, 1], py = m[, 2])) |>
  mutate(kind = if_else(edge %in% check_edges, "check", "flow"))

# the results collect on a bus down the right and enter the tables and figures
bus_x <- max(bx$xmax) + 2.5
into15 <- c("r13", "r14e", "r14", "r14bd")
stubs <- map_dfr(into15, function(i) tibble(edge = paste0(i, "_bus"), px = c(b(i)$xmax, bus_x), py = b(i)$y))
bus <- tibble(edge = "bus", px = bus_x, py = c(b("r13")$y, b("r15")$y))
bus_in <- tibble(edge = "bus_r15", px = c(bus_x, b("r15")$xmax), py = b("r15")$y)

fill_of <- c(source = "#F6E9CC", combat = "#DCE6F2", civil = "#F3DDD8", denom = "#D9EEE9",
             est = "#E7E0F0", res = "#DDEEDA", check = "#F2F2F2")
lt_of <- c(none = "solid", fixed = "solid", drawn = "22", check = "12")
bx <- bx |> mutate(fill = fill_of[track], lt = lt_of[enters],
                   border = if_else(enters == "check", "grey55", "grey25"))

headers <- tibble(x = col_x, y = 84,
                  label = c("Sources", "Preparation", "Estimation and checks", "Results"))
# the key: how an output enters the simulation, and the three tracks
key_y <- 3
key <- tibble(
  x = c(2, 22, 44, 66, 82, 98), w = 3.2,
  lt = c("22", "solid", "12", "solid", "solid", "solid"),
  fill = c("white", "white", "#F2F2F2", fill_of[["combat"]], fill_of[["civil"]], fill_of[["denom"]]),
  border = c("grey25", "grey25", "grey55", "grey25", "grey25", "grey25"),
  label = c("Drawn in every simulation", "Held; varied one at a time", "Check or comparison, no estimate",
            "Combatants", "Civilians", "Denominator and counterfactual"))


# every layer carries its own constant line type, so the boxes' borders and
# the edges need no shared scale
rect_layer <- function(d, lt) {
  geom_rect(data = d, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, colour = border),
            linetype = lt, linewidth = 0.35)
}
arrow_head <- arrow(length = unit(1.5, "mm"), type = "closed")
key_boxes <- key |> mutate(xmin = x, xmax = x + w, ymin = key_y - 1.1, ymax = key_y + 1.1)

p <-
  ggplot() +
  # the edges first, so the boxes sit on top of their ends
  geom_path(data = bind_rows(stubs, bus), aes(px, py, group = edge), colour = "grey45", linewidth = 0.3) +
  geom_path(data = bus_in, aes(px, py, group = edge), colour = "grey45", linewidth = 0.3, arrow = arrow_head) +
  geom_path(data = filter(edges, kind == "flow"), aes(px, py, group = edge),
            colour = "grey45", linewidth = 0.3, arrow = arrow_head) +
  geom_path(data = filter(edges, kind == "check"), aes(px, py, group = edge),
            colour = "grey60", linewidth = 0.3, linetype = "12", arrow = arrow_head) +
  # the boxes, by how their output enters the simulation
  rect_layer(filter(bx, lt == "solid"), "solid") +
  rect_layer(filter(bx, lt == "22"), "22") +
  rect_layer(filter(bx, lt == "12"), "12") +
  geom_text(data = bx, aes(x, y + 0.95, label = title), fontface = "bold", size = 2.9, colour = "grey10") +
  geom_text(data = bx, aes(x, y - 1.15, label = sub), size = 2.45, colour = "grey25") +
  # column headers
  geom_text(data = headers, aes(x, y, label = label), fontface = "bold", size = 3.4, colour = "grey15") +
  # the key
  rect_layer(filter(key_boxes, lt == "solid"), "solid") +
  rect_layer(filter(key_boxes, lt == "22"), "22") +
  rect_layer(filter(key_boxes, lt == "12"), "12") +
  geom_text(data = key_boxes, aes(xmax + 1, key_y, label = label), hjust = 0, size = 2.45, colour = "grey20") +
  scale_fill_identity() +
  scale_colour_identity() +
  coord_cartesian(xlim = c(0, bus_x + 2), ylim = c(0, 87), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

out <- "figures/figA0_pipeline.png"
ggsave(out, p, width = 12, height = 8.6, dpi = 300, device = ragg::agg_png, bg = "white")
message("  figure: ", basename(out))
cat("Done. ", out, " written.\n", sep = "")
