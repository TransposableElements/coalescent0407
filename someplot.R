#### Change in Ht over Generations for Different Population Sizes ####

library(tidyverse)
options(scipen = 1)
calculate_Ht <- function(H0, N, t) {
  H0 * ((1 - (1 / (2 * N))) ^ t)
}


H0 <- 1
generations <- 1:1000
population_sizes <- c(1, 10, 100, 1000000)


data <- expand.grid(Generation = generations, Population_Size = population_sizes) %>%
  mutate(Ht = map2_dbl(Population_Size, Generation, ~calculate_Ht(H0, .x, .y)))


ggplot(data, aes(x = Generation, y = Ht, color = as.factor(Population_Size), group = Population_Size)) +
  geom_line(size = 1.2) +  
  scale_y_log10() +  
  labs(x = "Generation",
       y = expression(H[t]),
       color = "Population Size",
       title = "Change in Ht over Generations for Different Population Sizes") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 14, face = "bold"), 
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

ggsave("Ht_change_over_generations.png", width = 10, height = 6, dpi = 900)

#### Exact vs Approximate Half-Life for Varying Population Size ####

exact_half_life <- function(N) {
  -log(2) / log(1 - 1/(2*N))
}

approx_half_life <- function(N) {
  2 * N * log(2)
}

N_values <- tibble(N = 1:1000)

half_lives <- N_values %>%
  mutate(
    Exact = map_dbl(N, exact_half_life),
    Approximate = map_dbl(N, approx_half_life)
  )

long_half_lives <- half_lives %>%
  pivot_longer(
    cols = c(Exact, Approximate),
    names_to = "Formula",
    values_to = "Half_life"
  )

ggplot(long_half_lives, aes(x = N, y = Half_life, color = Formula, group = Formula)) +
  geom_line(size = 1.4) +
  scale_color_manual(values = c("Exact" = "blue", "Approximate" = "red")) +
  theme_minimal() +
  labs(
    x = "Population Size (N)",
    y = "Half-Life (t1/2)",
    title = "Exact vs Approximate Half-Life for Varying Population Size",
    color = "Formula"
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggsave("half_life_plot.png", width = 10, height = 6, dpi = 900)

#### 4Nμ ####

x_values <- seq(1e-10, 1000, by = 0.01)
y_values <- x_values / (1 + x_values)

data <- data.frame(x = x_values, y = y_values)

ggplot(data, aes(x = x, y = y)) +
  geom_line(size = 1.4) +
  theme_minimal() +
  labs(x = "4Nμ", y = "H", title = "Plot of H = 4Nμ / (1 + 4Nμ)") +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggsave("plot_H_4Nμ.png", width = 8, height = 6, dpi = 900)

#### ΔH ####

N <- 10^4
u <- 5 * 10^-5

delta_H <- function(H) {
  (-1 / (2 * N)) * H + 2 * u * (1 - H)
}

H_values <- seq(0, 1, length.out = 1000)
delta_H_values <- delta_H(H_values)
data <- tibble(H = H_values, Delta_H = delta_H_values)

ggplot(data, aes(x = H, y = Delta_H)) +
  geom_line(size = 1.4) +
  labs(title = "Plot of ΔH as a Function of H",
       x = "H",
       y = "ΔH") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggsave("plot_ΔH_H.png", width = 8, height = 6, dpi = 900)


#### WF ####


library(viridis)

wf_df <- data.frame()
sizes <- c(50, 100, 1000, 5000)
starting_p <- c(.01, .1, .5, .8)

n_gen <- 100
n_reps <- 50

for(N in sizes){
  for(p in starting_p){
    p0 <- p
    for(j in 1:n_gen){
      X <- rbinom(n_reps, 2*N, p)
      p <- X / (2*N)
      rows <- data.frame(replicate = 1:n_reps, N = rep(N, n_reps), 
                         gen = rep(j, n_reps), p0 = rep(p0, n_reps), 
                         p = p)
      wf_df <- bind_rows(wf_df, rows)
    }
  }
}

p <- ggplot(wf_df, aes(x = gen, y = p, group = replicate, color = as.factor(gen))) +
  geom_path(alpha = .5) + facet_grid(N ~ p0) +
  scale_color_viridis_d() +  
  guides(color="none")

p
ggsave("plot_W_F.png", width = 8, height = 6, dpi = 900)


#### tajima ####

require(ape)
data(woodmouse)
tajima.test <- function(x)
{
  n <- if (is.list(x)) length(x) else dim(x)[1]
  if (n < 4) {
    warning("Tajima test requires at least 4 sequences")
    return(list(D = NaN, Pval.normal = NaN, Pval.beta = NaN))
  }
  khat <- mean(dist.dna(x, "N"))
  S <- length(seg.sites(x))
  if (!S) {
    warning("no segregating sites")
    return(list(D = NaN, Pval.normal = NaN, Pval.beta = NaN))
  }
  tmp <- 1:(n - 1)
  a1 <- sum(1/tmp)
  a2 <- sum(1/tmp^2)
  b1 <- (n + 1)/(3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
  c1 <- b1 - 1/a1
  c2 <- b2 - (n + 2)/(a1 * n) + a2/a1^2
  e1 <- c1/a1
  e2 <- c2/(a1^2 + a2)
  D <- (khat - S/a1)/sqrt(e1 * S + e2 * S * (S - 1))
  Dmin <- (2/n - 1/a1)/sqrt(e2)
  Dmax <- ((n + 1)/(2 * n) - 1/a1)/sqrt(e2)
  tmp1 <- 1 + Dmin * Dmax
  tmp2 <- Dmax - Dmin
  a <- -tmp1 * Dmax/tmp2
  b <- tmp1 * Dmin/tmp2
  p <- pbeta((D - Dmin)/tmp2, b, a)
  p <- if (p < 0.5) 2 * p else 2 * (1 - p)
  list(D = D, Pval.normal = 2 * pnorm(-abs(D)), Pval.beta = p)
}

tajima.test(woodmouse)




