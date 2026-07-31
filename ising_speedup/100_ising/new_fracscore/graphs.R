library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)
source("specific.R")

hellinger_data_frame <- read.csv("csvs/h_dist.csv")
work_energy          <- read.csv("csvs/work_energy.csv") |>
  select(-total, -count, -X)
energy_density_opt <- read.csv("csvs/energy_density_opt.csv") |>
  select(-X)


m_0_labeller <- function(string) { paste ("m_0 =", string) }
m_1_labeller <- function(string) { paste ("m_1 =", string) }

#hCDR Heat Maps

hcdr_work_energy <- work_energy |> filter("hcdr" == scheme)

procs <- unique(hcdr_work_energy$P)
for (proc in procs) {
  hcdr_work_energy_subset <- hcdr_work_energy |>
    filter(P == proc) |>
    complete(work_bin, energy_bin, m_0, m_1, P, lambda, scheme,
             fill = list(density = 0))

  work_bin_scale   <- seq(0, max(hcdr_work_energy_subset$work_bin), by = 3)
  energy_bin_scale <- seq(0, max(hcdr_work_energy_subset$energy_bin), by = 3)

  heat_map_plot <- ggplot(data = hcdr_work_energy_subset,
    aes(x = work_bin, y = energy_bin, fill = density)
  ) +
    scale_fill_gradientn(colors = c("white", "orange", "purple"),
      values = rescale(c(0, 0.05, max(hcdr_work_energy_subset$density)))
    ) +
    geom_tile() +
    coord_fixed() +
    scale_x_continuous(
      name   = "log(I*P)",
      breaks = work_bin_scale,
      expand = c(0, 0),
      labels = lapply(work_bin_scale, max_work_of_bin)
    ) +
    scale_y_continuous(
      name   = "E/E_min",
      breaks = energy_bin_scale,
      expand = c(0, 0),
      labels = lapply(energy_bin_scale, max_energy_of_bin)
    ) +
    facet_grid(
      rows = vars(m_0),
      cols = vars(m_1),
      scales = "fixed",
      space = "fixed",
      labeller = labeller(.rows = m_0_labeller, .cols = m_1_labeller)
    ) +
    theme(
      panel.spacing = unit(.01, "lines"),
      panel.border  = element_rect(color = "black",
                                   fill  = NA,
                                   linewidth = 0.2)
    ) +
    ggtitle(paste(
      paste("heat map ", proc, sep = ""),
      paste(" processors, ", name, sep = "")
    ))
  ggsave(
    filename = paste(
      paste("pdfs/heat_map_fracscore_", proc, sep = ""),
      paste(name, ".pdf", sep = ""), sep = ""
    ),
    plot = heat_map_plot,
    width = 14, height = 14
  )

  ggsave(
    filename = paste(
      paste("pdfs/heat_map_fracscore_", proc, sep = ""),
      paste(name, ".svg", sep = ""), sep = ""
    ),
    plot = heat_map_plot,
    width = 14, height = 14
  )
}

energy_density_opt <- unite(energy_density_opt,
                            SchemeAndProcessors,
                            P,
                            scheme,
                            lambda,
                            sep = " ",
                            remove = FALSE)

lambda_values <- energy_density_opt |>
  filter(P == 1) |>
  pull(lambda) |>
  unique() |>
  sort()

min_lambda <- min(lambda_values)

create_label <- function(p, scheme, lambda) {
 if (p == 1 && lambda == min_lambda) {
    "Serial"
  } else {
    paste(p, scheme)
  }
}

energy_density_opt <- energy_density_opt |>
  mutate(CleanLabel = mapply(create_label, P, scheme, lambda))


vpbp_levels <- c(
  "Serial",
  energy_density_opt |>
    filter(scheme == "hcdr") |>
    arrange(P) |>
    pull(CleanLabel) |>
    unique()
)


energy_density_opt$CleanLabel <- factor(
  energy_density_opt$CleanLabel,
  levels = vpbp_levels
)


energy_plot_base <- ggplot(
  data = energy_density_opt,
  aes(x = CleanLabel, y = E / min_energy)
) +
  xlab("Scheme and Processors") +
  ylab("E/E_min")+
  ggtitle(paste(name, " Energy Density", sep = "")) +
  theme(
    panel.border     = element_rect(colour = "black", fill = NA),
    plot.title       = element_text(hjust = 0.5)
  )

energy_box_plot <- energy_plot_base + geom_boxplot()
energy_violin_plot <- energy_plot_base +
  geom_violin(draw_quantiles = 0.5)

ggsave(energy_box_plot, filename = "pdfs/EnergyBoxplot.pdf")
ggsave(energy_box_plot, filename = "pdfs/EnergyBoxplot.svg")

ggsave(energy_violin_plot, filename = "pdfs/EnergyViolinplot.pdf", width=5, height=5)
ggsave(energy_violin_plot, filename = "pdfs/EnergyViolinplot.svg", width=5, height=5)


lims_work   <- seq(
  min(energy_density_opt$work_bin),
  max(energy_density_opt$work_bin),
  by = 3
)
lims_energy <- seq(
  min(energy_density_opt$energy_bin),
  max(energy_density_opt$energy_bin),
  by = 2
)

work_comparison_plot <- ggplot(
  data = energy_density_opt |>
    arrange(P) |>
    distinct(
      energy_bin,
      work_bin,
      P,
      scheme,
      lambda,
      .keep_all = TRUE
    ) |>
    group_by(energy_bin, work_bin, scheme) |>
    mutate(
      n_in_group = n(),
      group_index = row_number(),
      work_jitter = work_bin + (group_index - (n_in_group + 1) / 2) * 0.05
    ) |>
    ungroup(),
  aes(
    x = energy_bin,
    y = work_jitter,
    color = CleanLabel,
    group = CleanLabel,
    shape = CleanLabel
  )
) +
  geom_point(stroke = 1) +
  scale_shape_manual(values = c(1, 5, 0), name = NULL) +
  scale_color_manual(
    name = "Scheme and P",
    values = c(
      "Serial" = "black",
      "100 hcdr"    = "green",
      "200 hcdr"    = "red"
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
  ggtitle(paste(name, " Work Comparison", sep = "")) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text       = element_text(face = "bold"),
    panel.spacing    = unit(1.5, "lines"),
    panel.border     = element_rect(colour = "black", fill = NA),
    plot.title       = element_text(hjust = 0.5)
  )

ggsave(
  work_comparison_plot,
  filename = "pdfs/work_comparison.pdf",
  height = 5,
  width = 5
)

ggsave(
  work_comparison_plot,
  filename = "pdfs/work_comparison.svg",
  height = 5,
  width = 5
)
