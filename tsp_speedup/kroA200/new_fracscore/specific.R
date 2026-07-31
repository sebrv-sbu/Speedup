library(dplyr)
library(tidyr)
library(reshape2)


#Constants
min_energy            <- 29368
init_moves            <- 100000
energy_bin_resolution <- 100
work_bin_resolution   <- 10
name                  <- "kroA200"

#Helper Functions
max_energy_of_bin <- function(bin) {
  ((bin / energy_bin_resolution) + 1) * min_energy
}

max_work_of_bin <- function(bin) {
  bin / work_bin_resolution
}

max_frac_score_of_bin <- function(bin) {
  max_energy_of_bin(bin) / min_energy
}

hellinger_distance <- function(v_1, v_2) {
  (1 / sqrt(2)) * sqrt(sum((sqrt(v_1) - sqrt(v_2))^2))
}

min_hellinger_distance <- function(row_data, serial_energy_density) {
  if ("serial" == row_data$scheme) {
    list(distance = NA, min = NA)
  } else {
    row_e_density <- as.vector(t(row_data[1, -(1:5)]))
    result <- serial_energy_density |>
      rowwise() |>
      mutate(h_distance =
          hellinger_distance(
            as.vector(
              t(pick(everything())[1, -(1:5)])
            ), row_e_density
          )
      ) |>
      ungroup() |>
      slice_min(h_distance, n = 1) |>
      select(h_distance, lambda)

    list(distance = result$h_distance[1], min = result$lambda[1])
  }
}

estimate_par_eff <- function(min, n_proc, mean_iter, work_per_run) {
  if (is.na(min)) {
    NA
  } else {
    work_per_run[which(work_per_run$lambda == min), ]$mean_iter[1] /
      (mean_iter * n_proc)
  }
}
