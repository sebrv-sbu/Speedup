library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)
source("specific.R")

hellinger_data_frame <- read.csv("csvs/h_dist.csv")
work_energy          <- read.csv("csvs/work_energy.csv") |>
  select(-total, -count, -X)
energy_density_bias <- read.csv("csvs/energy_density_bias.csv") |>
  select(-X)

energy <- read.csv("csvs/energy.csv")

create_label <- function(p, scheme, lambda) {
  if (p == 1) {
    "Serial"
  } else {
    paste(p, scheme)
  }
}


m_0_labeller <- function(string) { paste ("m_0 =", string) }
m_1_labeller <- function(string) { paste ("m_1 =", string) }

#hCDR Heat Maps

bias_demo <- energy_density_bias |>
  filter(P == 1 | (P == 200 & m_0 == 1000 & m_1 == 1000))

energy_density_bias <- energy_density_bias |>
  mutate(CleanLabel = mapply(create_label, P, scheme, lambda))



energy_density_bias$CleanLabel <- factor(
  energy_density_bias$CleanLabel,
  levels = c("Serial", "200 hcdr")
)

energy_plot_base <- ggplot(
  data = energy_density_bias,
  aes(x = CleanLabel, y = E / min_energy)
) +
  xlab("Scheme and Processors") +
  ylab("E/E_min") +
  ggtitle("lin105 Best Serial Run and Best hCDR run 
          Energy Distribution") +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text       = element_text(face = "bold"),
    panel.spacing    = unit(1.5, "lines"),
    panel.border     = element_rect(colour = "black", fill = NA),
    plot.title       = element_text(hjust = 0.5)
  )

energy_violin_plot <- energy_plot_base + geom_violin(draw_quantiles=0.5)

ggsave(
  energy_violin_plot,
  filename = "pdfs/BiasViolinplot.pdf", 
  height = 5,
  width = 5
)
ggsave(energy_violin_plot, filename = "pdfs/BiasViolinplot.svg")

lims_work   <- seq(
  min(energy_density_bias$work_bin),
  max(energy_density_bias$work_bin),
  by = 3
)
lims_energy <- seq(
  min(energy_density_bias$energy_bin),
  max(energy_density_bias$energy_bin),
  by = 2
)

work_comparison_plot <- ggplot(
  data = energy_density_bias,
  aes(
    x = energy_bin,
    y = work_bin,
    color = CleanLabel,
    shape = scheme,
    group = scheme,
  )
) +
  geom_point(stroke = 1) +
  scale_shape_manual(values = c(5, 0)) +
  scale_color_manual(
    values = c(
      "Serial" = "orange",
      "200 hcdr" = "magenta"
    )
  ) +
  scale_y_continuous(
    name = "log(I*P)",
    breaks =  lims_work,
    expand = c(0.05, 0.05),
    labels = lapply(lims_work, max_work_of_bin)
  ) +
  scale_x_continuous(
    name = "E/E_min",
    breaks = lims_energy,
    expand = c(0.05, 0.05),
    labels = lapply(lims_energy, max_frac_score_of_bin)
  ) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text       = element_text(face = "bold"),
    panel.spacing    = unit(1.5, "lines"),
    panel.border     = element_rect(colour = "black", fill = NA),
    plot.title       = element_text(hjust = 0.5)
  ) +
  ggtitle("lin105 Work Done by Best Serial and hCDR runs.")

ggsave(
  work_comparison_plot,
  filename = "pdfs/bias_work_comparison.pdf",
  height = 4,
  width = 4 * 1.5
)
