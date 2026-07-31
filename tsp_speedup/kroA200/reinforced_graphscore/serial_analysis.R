source("specific.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(scales)

data_serial <- read.csv("./csvs/energy.csv") |>
  filter(1 == P)
data_serial_avg <- read.csv("./csvs/work_per_run.csv") |>
  filter(1 == P)

data_serial$log_I <- log(data_serial$tot_iter_bin)
lambda_int_model <- lm(
  log(data_serial_avg |> filter(lambda < 0.1) |> pull(mean_iter)) ~
    log(data_serial_avg |> filter(lambda < 0.1) |> pull(lambda))
)
summary(lambda_int_model)
lambda_int_model_log_10 <- lambda_int_model <- lm(
  log(data_serial_avg |> filter(lambda < 0.1) |> pull(mean_iter)) ~
    log10(data_serial_avg |> filter(lambda < 0.1) |> pull(lambda))
)

min_density <- min(data_serial$density)
max_density <- max(data_serial$density)
spread_density <- max_density - min_density

c_bin_e <- min(data_serial$log_I)
m_bin_e <- max(data_serial$log_I) - c_bin_e
transform_bin_e <- function(x) {
  x <- x - min(data_serial$energy_bin)
  x <- x / (max(data_serial$energy_bin) - min(data_serial$energy_bin))
  m_bin_e * x + c_bin_e
}
inverse_bin_e <- function(x) {
  x <- (x - c_bin_e) / m_bin_e
  x <- x * (max(data_serial$energy_bin) - min(data_serial$energy_bin))
  x + min(data_serial$energy_bin)
}



serial_plot <- ggplot(data = data_serial, aes(x = log10(lambda))) +
  geom_point(
    aes(y = transform_bin_e(energy_bin), fill = density),
    shape = 22,
    size = 1,
    stroke = NA
  ) +
  geom_point(
    size = 0.5,
    aes(y = log(mean_iter), color = mean_frac_score),
    data = data_serial_avg
  ) +
  scale_color_gradientn(
    name = "E/E_min", colors = c("orange", "purple", "green"),
    values = rescale(c(1, 1.07, max(data_serial$energy_bin)))
  ) +
  scale_fill_gradientn(
    colors = c("#66cc00", "#cc3399", "#43123c"),
    values = rescale(c(0, 0.2, 0.6))
  ) +
  scale_y_continuous(
    name = "log iterations",
    sec.axis = sec_axis(
      transform = ~inverse_bin_e(.),
      name = "energy score"
    )
  ) +
  scale_x_continuous(name = "log_10 lambda") +
  geom_abline(
    intercept = coef(lambda_int_model_log_10)[1],
    slope = coef(lambda_int_model_log_10)[2]
  ) +
  ggtitle(paste(name, "summary of serial runs")) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text       = element_text(face = "bold"),
    panel.spacing    = unit(1.5, "lines"),
    panel.border     = element_rect(colour = "black", fill = NA),
    plot.title       = element_text(hjust = 0.5)
  )
ggsave(serial_plot, filename = "pdfs/serial_summary.pdf", width=5, height=5)
ggsave(serial_plot, filename = "pdfs/serial_summary.svg", width=5, height=5)
