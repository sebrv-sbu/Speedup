library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)
library(patchwork)
source("specific.R")

hellinger_data_frame <- read.csv("csvs/h_dist.csv")

rank_info <- hellinger_data_frame |>
  filter(P > 1) |>
  select(scheme, P, m_0, m_1, rank, bin_width) |>
  group_by(scheme) |>
  arrange(P, m_0, m_1, .group_by = TRUE)


print(rank_info)

work_energy          <- read.csv("csvs/work_energy_1.csv") |>
  select(-total, -count, -X)
energy_density_1_opt <- read.csv("csvs/energy_density_1_opt.csv") |>
  select(-X)

energy_density_2_opt <- read.csv("csvs/energy_density_2_opt.csv") |>
  select(-X)

energy_density_05_opt <- read.csv("csvs/energy_density_05_opt.csv") |>
  select(-X)






m_0_labeller <- function(string) { paste ("m_0 =", string) }
m_1_labeller <- function(string) { paste ("m_1 =", string) }

#hCDR Heat Maps

hcdr_work_energy_1 <- work_energy |> filter("hcdr" == scheme)

energy_density_1_opt <- unite(energy_density_1_opt,
                              SchemeAndProcessors,
                              P,
                              scheme,
                              lambda,
                              sep = " ",
                              remove = FALSE)

lambda_values <- energy_density_1_opt |>
  filter(P == 1) |>
  pull(lambda) |>
  unique() |>
  sort()

max_lambda <- max(lambda_values)
min_lambda <- min(lambda_values)

create_label <- function(p, scheme, lambda) {
  if (p == 1 && lambda == max_lambda) {
    "Serial High"
  } else if (p == 1 && lambda == min_lambda) {
    "Serial Low"
  } else {
    paste(p, scheme)
  }
}

energy_density_1_opt <- energy_density_1_opt |>
  mutate(CleanLabel = mapply(create_label, P, scheme, lambda))

cdr_levels <- energy_density_1_opt |>
  filter(scheme == "cdr") |>
  arrange(P) |>
  pull(CleanLabel) |>
  unique()

hcdr_levels <- energy_density_1_opt |>
  filter(scheme == "hcdr") |>
  arrange(P) |>
  pull(CleanLabel) |>
  unique()

serial_cdr <- energy_density_1_opt |>
  filter(P == 1, lambda == max_lambda) |>
  pull(CleanLabel) |>
  unique()

serial_hcdr <- energy_density_1_opt |>
  filter(P == 1, lambda == min_lambda) |>
  pull(CleanLabel) |>
  unique()


vpbp_levels <- c(serial_cdr, cdr_levels, serial_hcdr, hcdr_levels)

energy_density_1_opt$CleanLabel <- factor(
  energy_density_1_opt$CleanLabel,
  levels = vpbp_levels
)

energy_density_1_opt$Group <- ifelse(
  energy_density_1_opt$CleanLabel %in% c(serial_cdr, cdr_levels),
  "CDR",
  "hCDR"
)

energy_plot_base <- ggplot(
  data = energy_density_1_opt,
  aes(x = CleanLabel, y = E / min_energy)
) +
  xlab("Scheme and Processors") +
  ylab("E/E_min")+
  ggtitle(paste(name, " Energy Density", sep = "")) +
  facet_grid(
    cols   = vars(Group),
    scales = "free_x",
    space  = "free_x"
  ) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text       = element_text(face = "bold"),
    panel.spacing    = unit(1.5, "lines"),
    panel.border     = element_rect(colour = "black", fill = NA),
    plot.title       = element_text(hjust = 0.5)
  )

energy_box_plot <- energy_plot_base + geom_boxplot()
energy_violin_plot <- energy_plot_base +
  geom_violin(draw_quantiles = 0.5)

ggsave(energy_box_plot, filename = "pdfs/EnergyBoxplot.pdf")
ggsave(energy_box_plot, filename = "pdfs/EnergyBoxplot.svg")

ggsave(energy_violin_plot, filename = "pdfs/EnergyViolinplot.pdf", width = 5,
       height = 5)
ggsave(energy_violin_plot, filename = "pdfs/EnergyViolinplot.svg", width = 5,
       height = 5)


lims_work   <- seq(
  min(energy_density_1_opt$work_bin),
  max(energy_density_1_opt$work_bin),
  by = 3
)
lims_energy <- seq(
  min(energy_density_1_opt$energy_bin_1),
  max(energy_density_1_opt$energy_bin_1),
  by = 2
)

work_comparison_plot <- ggplot(
  data = energy_density_1_opt |>
    arrange(P) |>
    distinct(
      energy_bin_1,
      work_bin,
      P,
      scheme,
      lambda,
      .keep_all = TRUE
    ) |>
    group_by(energy_bin_1, work_bin, scheme) |>
    mutate(
      n_in_group = n(),
      group_index = row_number(),
      work_jitter = work_bin + (group_index - (n_in_group + 1) / 2) * 0.05
    ) |>
    ungroup(),
  aes(
    x = energy_bin_1,
    y = work_jitter,
    color = CleanLabel,
    shape = scheme,
    group = scheme,
  )
) +
  geom_point(stroke = 1) +
  scale_shape_manual(values = c(1, 5, 0), name = NULL) +
  scale_color_manual(
    name = "Scheme and P",
    values = c(
      "Serial High" = "brown",
      "Serial Low"  = "orange",
      "5 cdr"       = "green",
      "10 cdr"      = "red",
      "20 cdr"      = "blue",
      "20 hcdr"     = "black",
      "100 hcdr"    = "purple",
      "200 hcdr"    = "magenta"
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

energy_density_2_opt <- unite(energy_density_2_opt,
                              SchemeAndProcessors,
                              P,
                              scheme,
                              lambda,
                              sep = " ",
                              remove = FALSE)

lambda_values <- energy_density_2_opt |>
  filter(P == 1) |>
  pull(lambda) |>
  unique() |>
  sort()

max_lambda <- max(lambda_values)
min_lambda <- min(lambda_values)

create_label <- function(p, scheme, lambda) {
  if (p == 1 && lambda == max_lambda) {
    "Serial High"
  } else if (p == 1 && lambda == min_lambda) {
    "Serial Low"
  } else {
    paste(p, scheme)
  }
}

energy_density_2_opt <- energy_density_2_opt |>
  mutate(CleanLabel = mapply(create_label, P, scheme, lambda))

cdr_levels <- energy_density_2_opt |>
  filter(scheme == "cdr") |>
  arrange(P) |>
  pull(CleanLabel) |>
  unique()

hcdr_levels <- energy_density_2_opt |>
  filter(scheme == "hcdr") |>
  arrange(P) |>
  pull(CleanLabel) |>
  unique()

serial_cdr <- energy_density_2_opt |>
  filter(P == 1, lambda == max_lambda) |>
  pull(CleanLabel) |>
  unique()

serial_hcdr <- energy_density_2_opt |>
  filter(P == 1, lambda == min_lambda) |>
  pull(CleanLabel) |>
  unique()


vpbp_levels <- c(serial_cdr, cdr_levels, serial_hcdr, hcdr_levels)

energy_density_2_opt$CleanLabel <- factor(
  energy_density_2_opt$CleanLabel,
  levels = vpbp_levels
)

energy_density_2_opt$Group <- ifelse(
  energy_density_2_opt$CleanLabel %in% c(serial_cdr, cdr_levels),
  "CDR",
  "hCDR"
)

energy_plot_base <- ggplot(
  data = energy_density_2_opt,
  aes(x = CleanLabel, y = E / min_energy)
) +
  xlab("Scheme and Processors") +
  ylab("E/E_min")+
  ggtitle(paste(name, " Energy Density", sep = "")) +
  facet_grid(
    cols   = vars(Group),
    scales = "free_x",
    space  = "free_x"
  ) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text       = element_text(face = "bold"),
    panel.spacing    = unit(1.5, "lines"),
    panel.border     = element_rect(colour = "black", fill = NA),
    plot.title       = element_text(hjust = 0.5)
  )

energy_box_plot <- energy_plot_base + geom_boxplot()
energy_violin_plot <- energy_plot_base +
  geom_violin(draw_quantiles = 0.5)

ggsave(energy_box_plot, filename = "pdfs/EnergyBoxplot_2.pdf")
ggsave(energy_box_plot, filename = "pdfs/EnergyBoxplot_2.svg")

ggsave(energy_violin_plot, filename = "pdfs/EnergyViolinplot_2.pdf", width=5, height=5)
ggsave(energy_violin_plot, filename = "pdfs/EnergyViolinplot_2.svg", width=5, height=5)


lims_work   <- seq(
  min(energy_density_2_opt$work_bin),
  max(energy_density_2_opt$work_bin),
  by = 3
)
lims_energy <- seq(
  min(energy_density_2_opt$energy_bin_2),
  max(energy_density_2_opt$energy_bin_2),
  by = 2
)

work_comparison_plot <- ggplot(
  data = energy_density_2_opt |>
    arrange(P) |>
    distinct(
      energy_bin_1,
      work_bin,
      P,
      scheme,
      lambda,
      .keep_all = TRUE
    ) |>
    group_by(energy_bin_2, work_bin, scheme) |>
    mutate(
      n_in_group = n(),
      group_index = row_number(),
      work_jitter = work_bin + (group_index - (n_in_group + 1) / 2) * 0.05
    ) |>
    ungroup(),
  aes(
    x = energy_bin_2,
    y = work_jitter,
    color = CleanLabel,
    shape = scheme,
    group = scheme,
  )
) +
  geom_point(stroke = 1) +
  scale_shape_manual(values = c(1, 5, 0), name = NULL) +
  scale_color_manual(
    name = "Scheme and P",
    values = c(
      "Serial High" = "brown",
      "Serial Low"  = "orange",
      "5 cdr"       = "green",
      "10 cdr"      = "red",
      "20 cdr"      = "blue",
      "20 hcdr"     = "black",
      "100 hcdr"    = "purple",
      "200 hcdr"    = "magenta"
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
  filename = "pdfs/work_comparison_2.pdf",
  height = 5,
  width = 5
)

ggsave(
  work_comparison_plot,
  filename = "pdfs/work_comparison_2.svg",
  height = 5,
  width = 5
)

energy_density_05_opt <- unite(energy_density_05_opt,
                              SchemeAndProcessors,
                              P,
                              scheme,
                              lambda,
                              sep = " ",
                              remove = FALSE)

lambda_values <- energy_density_05_opt |>
  filter(P == 1) |>
  pull(lambda) |>
  unique() |>
  sort()

max_lambda <- max(lambda_values)
min_lambda <- min(lambda_values)

create_label <- function(p, scheme, lambda) {
  if (p == 1 && lambda == max_lambda) {
    "Serial High"
  } else if (p == 1 && lambda == min_lambda) {
    "Serial Low"
  } else {
    paste(p, scheme)
  }
}

energy_density_05_opt <- energy_density_05_opt |>
  mutate(CleanLabel = mapply(create_label, P, scheme, lambda))

cdr_levels <- energy_density_05_opt |>
  filter(scheme == "cdr") |>
  arrange(P) |>
  pull(CleanLabel) |>
  unique()

hcdr_levels <- energy_density_05_opt |>
  filter(scheme == "hcdr") |>
  arrange(P) |>
  pull(CleanLabel) |>
  unique()

serial_cdr <- energy_density_05_opt |>
  filter(P == 1, lambda == max_lambda) |>
  pull(CleanLabel) |>
  unique()

serial_hcdr <- energy_density_05_opt |>
  filter(P == 1, lambda == min_lambda) |>
  pull(CleanLabel) |>
  unique()


vpbp_levels <- c(serial_cdr, cdr_levels, serial_hcdr, hcdr_levels)

energy_density_05_opt$CleanLabel <- factor(
  energy_density_05_opt$CleanLabel,
  levels = vpbp_levels
)

energy_density_05_opt$Group <- ifelse(
  energy_density_05_opt$CleanLabel %in% c(serial_cdr, cdr_levels),
  "CDR",
  "hCDR"
)

energy_plot_base <- ggplot(
  data = energy_density_05_opt,
  aes(x = CleanLabel, y = E / min_energy)
) +
  xlab("Scheme and Processors") +
  ylab("E/E_min")+
  ggtitle(paste(name, " Energy Density", sep = "")) +
  facet_grid(
    cols   = vars(Group),
    scales = "free_x",
    space  = "free_x"
  ) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text       = element_text(face = "bold"),
    panel.spacing    = unit(1.5, "lines"),
    panel.border     = element_rect(colour = "black", fill = NA),
    plot.title       = element_text(hjust = 0.5)
  )

energy_box_plot <- energy_plot_base + geom_boxplot()
energy_violin_plot <- energy_plot_base +
  geom_violin(draw_quantiles = 0.5)

ggsave(energy_box_plot, filename = "pdfs/EnergyBoxplot_05.pdf")
ggsave(energy_box_plot, filename = "pdfs/EnergyBoxplot_05.svg")

ggsave(energy_violin_plot, filename = "pdfs/EnergyViolinplot_05.pdf", width=5, height=5)
ggsave(energy_violin_plot, filename = "pdfs/EnergyViolinplot_05.svg", width=5, height=5)


lims_work   <- seq(
  min(energy_density_05_opt$work_bin),
  max(energy_density_05_opt$work_bin),
  by = 3
)
lims_energy <- seq(
  min(energy_density_05_opt$energy_bin_05),
  max(energy_density_05_opt$energy_bin_05),
  by = 2
)

work_comparison_plot <- ggplot(
  data = energy_density_05_opt |>
    arrange(P) |>
    distinct(
      energy_bin_1,
      work_bin,
      P,
      scheme,
      lambda,
      .keep_all = TRUE
    ) |>
    group_by(energy_bin_05, work_bin, scheme) |>
    mutate(
      n_in_group = n(),
      group_index = row_number(),
      work_jitter = work_bin + (group_index - (n_in_group + 1) / 2) * 0.05
    ) |>
    ungroup(),
  aes(
    x = energy_bin_05,
    y = work_jitter,
    color = CleanLabel,
    shape = scheme,
    group = scheme,
  )
) +
  geom_point(stroke = 1) +
  scale_shape_manual(values = c(1, 5, 0), name = NULL) +
  scale_color_manual(
    name = "Scheme and P",
    values = c(
      "Serial High" = "brown",
      "Serial Low"  = "orange",
      "5 cdr"       = "green",
      "10 cdr"      = "red",
      "20 cdr"      = "blue",
      "20 hcdr"     = "black",
      "100 hcdr"    = "purple",
      "200 hcdr"    = "magenta"
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
  filename = "pdfs/work_comparison_05.pdf",
  height = 5,
  width = 5
)

ggsave(
  work_comparison_plot,
  filename = "pdfs/work_comparison_05.svg",
  height = 5,
  width = 5
)

lambda_proc_05 <- energy_density_05_opt |>
  select(P, scheme, m_0, m_1, lambda) |>
  distinct(P, scheme, m_0, m_1, lambda)

lambda_proc_1 <- energy_density_1_opt |>
  select(P, scheme, m_0, m_1, lambda) |>
  distinct(P, scheme, m_0, m_1, lambda)

lambda_proc_2 <- energy_density_2_opt |>
  select(P, scheme, m_0, m_1, lambda) |>
  distinct(P, scheme, m_0, m_1, lambda)

print(lambda_proc_05)
print(lambda_proc_1)
print(lambda_proc_2)

write.csv(rank_info, "./csvs/rank_info.csv")
rank_info <- rank_info |>
  pivot_wider(names_from = bin_width, values_from = rank) |>
  pivot_longer(c("0.005", "0.01", "0.02"), names_to = "bin_width",
               values_to = "rank") |>
  rename("wasserstein_rank" = "0")

scheme_and_ps <- rank_info |>
  select(scheme, P) |>
  unique() |>
  arrange(scheme, P)

plots <- c()
print(scheme_and_ps)

for (i in seq_len(nrow(scheme_and_ps))) {
  size_of_point = 2
  scheme_and_p <- scheme_and_ps[i, ]
  print(scheme_and_p)
  p <- ggplot(
    data = rank_info |>
      filter(scheme == scheme_and_p$scheme & P == scheme_and_p$P),
    aes(
      x = wasserstein_rank,
      y = rank,
      color = bin_width,
      shape = bin_width,
      size = bin_width
    )
  ) +
    geom_point(stroke = 1) +
    scale_color_manual(
      name = "Bin Width",
      values = c(
        "wasserstein" = "black",
        "0.005" = "red",
        "0.01"  = "darkgreen",
        "0.02"  = "darkblue"
      )
    ) +
    scale_shape_manual(values = c(1, 5, 0), guide = "none") +
    scale_size_manual(
      values = c(
        "0.005" = 0.7 * size_of_point,
        "0.01"  = size_of_point,
        "0.02"  = 1.3 * size_of_point
      ),
      guide = "none"
    ) +
    labs(
      y = "Hellinger Rank",
      x = "Wasserstein Rank"
    ) +
    ggtitle(paste(scheme_and_p$scheme, scheme_and_p$P, sep = ", "))
  plots <- append(plots, p)
}

final_p <- plots[[1]]

for (i in seq_len(length(plots) - 1)){
  final_p <- final_p + plots[[i + 1]]
}

final_p <- final_p + plot_layout(nrow = 3)
ggsave("pdfs/ranks.pdf", width = 12, height = 12)
