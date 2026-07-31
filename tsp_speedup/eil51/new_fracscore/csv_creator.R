library(dplyr)
library(tidyr)
library(reshape2)
source("specific.R")

data_s     <- read.table("../s_results", header = TRUE)
data_5_c   <- read.table("../5_c_results", header = TRUE)
data_10_c  <- read.table("../10_c_results", header = TRUE)
data_20_c  <- read.table("../20_c_results", header = TRUE)
data_20_h  <- read.table("../20_h_results", header = TRUE)
data_100_h <- read.table("../100_h_results", header = TRUE)
data_200_h <- read.table("../200_h_results", header = TRUE)

#Format Data

data_s    $P  <- rep(1,   nrow(data_s))
data_5_c  $P  <- rep(5,   nrow(data_5_c))
data_10_c $P  <- rep(10,  nrow(data_10_c))
data_20_c $P  <- rep(20,  nrow(data_20_c))
data_20_h $P  <- rep(20,  nrow(data_20_h))
data_100_h$P  <- rep(100, nrow(data_100_h))
data_200_h$P  <- rep(200, nrow(data_200_h))

cdr_temp  <- bind_rows(data_5_c, data_10_c, data_20_c)
hcdr_temp <- bind_rows(data_20_h, data_100_h, data_200_h)

names(cdr_temp)["Interval" == names(cdr_temp)] <- "m_0"

cdr_temp $scheme <- rep("cdr",    nrow(cdr_temp))
hcdr_temp$scheme <- rep("hcdr",   nrow(hcdr_temp))
data_s   $scheme <- rep("serial", nrow(data_s))

cdr_temp $lambda <- rep(0.000005, nrow(cdr_temp))
hcdr_temp$lambda <- rep(0.000005, nrow(hcdr_temp))

#remove initial moves

cdr_temp $I <- cdr_temp $I - init_moves - (init_moves / cdr_temp$P)
hcdr_temp$I <- hcdr_temp$I - init_moves - (init_moves /
  case_when(
    as.integer(hcdr_temp$P) == 20 ~ 4,
    TRUE                          ~ 20
  )
)
#In the case of 20 processors, we need to divide by 4, and
#in the case of 100 or 200, we need to divide by 20.

data_s$I <- data_s$I - init_moves * 2

data_univ <- bind_rows(data_s, cdr_temp, hcdr_temp)

data_univ$log_I      <- log(data_univ$I)
data_univ$log_work   <- data_univ$log_I + log(data_univ$P)
data_univ$frac_score <- data_univ$E / min_energy
data_univ$energy_bin <- ceiling(energy_bin_resolution * (data_univ$frac_score - 1))
data_univ$work_bin   <- ceiling(work_bin_resolution * data_univ$log_work)

#Energy Only.



energy_density <- data_univ |>
  group_by(lambda, energy_bin, P, m_0, m_1, scheme) |>
  summarise(count = n(), tot_iter_bin = sum(I), .groups = "drop") |>
  group_by(lambda, P, m_0, m_1, scheme) |>
  mutate(total = sum(count),
    density = count / total
  )

work_per_run <- energy_density |>
  group_by(lambda, P, m_0, m_1, scheme) |>
  summarise(
    mean_iter = sum(tot_iter_bin) / first(total),
  ) |>
  left_join(data_univ |>
      group_by(lambda, P, m_0, m_1, scheme) |>
      summarise(mean_frac_score = mean(frac_score), .groups = "drop")
  )



hellinger_dist_dataframe <- energy_density |>
  select(-tot_iter_bin) |>
  melt(
    id = c(
      "lambda", "energy_bin", "P", "m_0", "m_1", "scheme"
    ), na.rm = FALSE
  ) |>
  filter("density" == variable) |>
  complete(
    nesting(
      lambda, P, m_0, m_1, scheme
    ),
    energy_bin,
    fill = list(value = 0)
  ) |>
  select(-variable) |>
  pivot_wider(names_from = energy_bin, values_from = value)

serial_e_d <- hellinger_dist_dataframe |> filter("serial" == scheme)


hellinger_dist_dataframe <- hellinger_dist_dataframe |>
  rowwise() |>
  mutate(min_and_distance =
      list(min_hellinger_distance(pick(everything()), serial_e_d)
      )
  ) |>
  unnest_wider(min_and_distance) |>
  select(lambda, P, m_0, m_1, scheme, min, distance) |>
  full_join(work_per_run) |>
  rowwise() |>
  mutate(par_efficiency = estimate_par_eff(min, P, mean_iter, work_per_run))


#Now we approximate the parallel efficiency.
#Fuck.




#Energy and Iterations

work_energy_density <- data_univ |>
  group_by(lambda, energy_bin, work_bin, P, m_0, m_1, scheme) |>
  summarise(count = n(), .groups = "drop") |>
  group_by(lambda, P, m_0, m_1, scheme) |>
  mutate(total = sum(count),
    density = count / total
  )


#Graph Data Frames
max_h <- hellinger_dist_dataframe |>
  filter(P != 1) |>
  group_by(P, lambda, scheme) |>
  mutate(max_par_efficiency = max(par_efficiency)) |>
  filter(par_efficiency == max_par_efficiency)

energy_density_opt <- data_univ |>
  semi_join(max_h, by = c("m_0", "m_1", "P", "scheme"))

max_log_work <- energy_density_opt |>
  filter("hcdr" == scheme) |>
  group_by(P) |>
  summarise(mean_log_work = mean(log_work)) |>
  pull(mean_log_work) |>
  max()

max_log_work_cdr <- energy_density_opt |>
  filter("cdr" == scheme) |>
  group_by(P) |>
  summarise(mean_log_work = mean(log_work)) |>
  pull(mean_log_work) |>
  max()



serial_run_comp <- data_univ |>
  filter(1 == P) |>
  semi_join(
    work_per_run |>
      ungroup() |>
      filter(1 == P & log(mean_iter) > max_log_work) |>
      arrange(mean_iter) |>
      filter(lambda == first(lambda)),
    by = c("lambda", "P")
  )
serial_run_comp_cdr <- data_univ |>
  filter(1 == P) |>
  semi_join(
    work_per_run |>
      ungroup() |>
      filter(1 == P & log(mean_iter) > max_log_work_cdr) |>
      arrange(mean_iter) |>
      filter(lambda == first(lambda)),
    by = c("lambda", "P")
  )

serial_run_bias <- data_univ |>
  filter(1 == P) |>
  semi_join(
    work_per_run |>
      ungroup() |>
      arrange(mean_iter) |>
      filter(lambda == last(lambda)),
    by = c("lambda", "P")
  )

hcdr_run_bias <- data_univ |> 
  filter(200 == P & m_0 == 1000 & m_1 == 1000)

energy_density_bias <- hcdr_run_bias |>
  bind_rows(serial_run_bias)

energy_density_opt <- energy_density_opt |>
  bind_rows(serial_run_comp) |>
  bind_rows(serial_run_comp_cdr)

write.csv(hellinger_dist_dataframe, "./csvs/h_dist.csv")
write.csv(work_energy_density, "./csvs/work_energy.csv")
write.csv(data_univ, "./csvs/data_univ.csv")
write.csv(energy_density, "./csvs/energy.csv")
write.csv(energy_density_opt, "./csvs/energy_density_opt.csv")
write.csv(work_per_run, "./csvs/work_per_run.csv")
write.csv(energy_density_bias, "./csvs/energy_density_bias.csv")
