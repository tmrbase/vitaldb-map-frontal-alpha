# Purpose: create the residual diagnostics reported in Table S2 and Figure S2.
source("R/common.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

objects <- readRDS(file.path(output_root, "r_random_slope_primary_objects.rds"))
fit <- objects$nonlinear
residual_values <- as.numeric(resid(fit, type = "normalized"))
fitted_values <- as.numeric(fitted(fit, level = 1))
acf_table <- as.data.table(nlme::ACF(fit, resType = "normalized", maxLag = 8))
setnames(acf_table, c("lag", "ACF"), c("Lag", "Residual_ACF"))

skewness <- mean((residual_values - mean(residual_values))^3) / sd(residual_values)^3
excess_kurtosis <- mean((residual_values - mean(residual_values))^4) / sd(residual_values)^4 - 3
summary_table <- data.table(
  Metric = c(
    "Normalized residual mean", "Normalized residual SD", "Skewness",
    "Excess kurtosis", "Residual 2.5th percentile", "Residual median",
    "Residual 97.5th percentile"
  ),
  Value = c(
    mean(residual_values), sd(residual_values), skewness, excess_kurtosis,
    quantile(residual_values, 0.025), median(residual_values), quantile(residual_values, 0.975)
  )
)
acf_export <- copy(acf_table)
acf_export[, Metric := sprintf("Normalized residual ACF, lag %d", Lag)]
acf_export <- acf_export[, .(Metric, Value = Residual_ACF)]
diagnostics <- rbindlist(list(summary_table, acf_export))
fwrite(diagnostics, file.path(output_root, "r_model_residual_diagnostics.csv"))

set.seed(202608273L)
plot_sample <- data.table(Fitted = fitted_values, Residual = residual_values)
if (nrow(plot_sample) > 10000L) plot_sample <- plot_sample[sample(.N, 10000L)]

theme_diagnostic <- theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#EAF2F8", color = "#6C8EA3"),
    strip.text = element_text(face = "bold")
  )

p1 <- ggplot(plot_sample, aes(Fitted, Residual)) +
  geom_hline(yintercept = 0, color = "#555555", linewidth = 0.4) +
  geom_point(alpha = 0.13, size = 0.65, color = "#0072B2") +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "#D55E00", linewidth = 0.7) +
  labs(title = "A. Residuals versus conditional fitted values", x = "Conditional fitted logit(RA)", y = "Normalized residual") +
  theme_diagnostic

qq <- qqnorm(residual_values, plot.it = FALSE)
qq_data <- data.table(Theoretical = qq$x, Observed = qq$y)
p2 <- ggplot(qq_data, aes(Theoretical, Observed)) +
  geom_abline(slope = 1, intercept = 0, color = "#555555", linewidth = 0.5) +
  geom_point(alpha = 0.15, size = 0.55, color = "#0072B2") +
  labs(title = "B. Normal Q-Q plot", x = "Theoretical quantile", y = "Observed normalized residual") +
  theme_diagnostic

p3 <- ggplot(acf_table[Lag > 0], aes(Lag, Residual_ACF)) +
  geom_hline(yintercept = 0, color = "#555555", linewidth = 0.4) +
  geom_segment(aes(xend = Lag, y = 0, yend = Residual_ACF), color = "#009E73", linewidth = 0.8) +
  geom_point(color = "#009E73", size = 2) +
  scale_x_continuous(breaks = acf_table$Lag) +
  labs(title = "C. Residual autocorrelation after CAR(1)", x = "Observation lag", y = "Residual ACF") +
  theme_diagnostic

layout <- grid::grid.layout(2, 2)
save_plot <- function(path, device, width, height, units, res = NULL) {
  if (is.null(res)) device(path, width = width, height = height, units = units) else device(path, width = width, height = height, units = units, res = res)
  grid.newpage()
  pushViewport(viewport(layout = layout))
  print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p3, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  dev.off()
}

save_plot(file.path(release_asset_root, "FigureS2_residual_diagnostics.png"), png, 2100, 1800, "px", 300)
save_plot(file.path(release_asset_root, "FigureS2_residual_diagnostics.tiff"), tiff, 7, 6, "in", 300)
pdf(file.path(release_asset_root, "FigureS2_residual_diagnostics.pdf"), width = 7, height = 6, useDingbats = FALSE)
grid.newpage(); pushViewport(viewport(layout = layout))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p3, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
dev.off()

message("Residual diagnostics complete")
