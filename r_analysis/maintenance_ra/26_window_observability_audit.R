# Purpose: audit EEG-window observability and create Table S1 and Figure S3 inputs.
source("R/common.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

input_root <- file.path(repo_root, "outputs", "maintenance_ra_full_v1")

message("[1/4] Reading scheduled windows and defining the audit denominator")
windows <- fread(file.path(input_root, "full_windows.csv"))
cases <- fread(file.path(input_root, "eligible_cases.csv"))

primary_counts <- windows[logical_true(primary_complete_window), .N, by = caseid]
analytic_caseids <- primary_counts[N >= 8L, caseid]

case_fields <- cases[, .(
  caseid,
  operation_duration_minutes = (opend - opstart) / 60,
  emergency_case = as.integer(emop > 0)
)]
windows <- merge(windows, case_fields, by = "caseid", all.x = TRUE, sort = FALSE)

clinical_opportunity <-
  windows$caseid %in% analytic_caseids &
  logical_true(windows$art_qc) &
  logical_true(windows$broad_no_positive_volatile) &
  is.finite(windows$propofol_ce_reconstructed) &
  is.finite(windows$remifentanil_ce_reconstructed) &
  is.finite(windows$etco2) &
  is.finite(windows$body_temperature) &
  is.finite(windows$propofol_rate_observed_fraction) &
  windows$propofol_rate_observed_fraction >= 0.80 &
  is.finite(windows$remifentanil_rate_observed_fraction) &
  windows$remifentanil_rate_observed_fraction >= 0.80

opportunities <- copy(windows[clinical_opportunity])
setorder(opportunities, subjectid, caseid, time_seconds)
opportunities[, eligible := as.integer(logical_true(primary_ra_eligible))]
opportunities[, map_category := cut(
  art_map,
  breaks = c(-Inf, 60, 65, 70, 80, Inf),
  right = FALSE,
  labels = c("<60", "60-64", "65-69", "70-79", ">=80")
)]
opportunities[, time_category := cut(
  minutes_since_opstart,
  breaks = c(-Inf, 30, 60, 120, 180, 240, Inf),
  right = FALSE,
  labels = c("<30", "30-59", "60-119", "120-179", "180-239", ">=240")
)]

rate_table <- function(frame, grouping, label_name) {
  output <- frame[, .(
    Scheduled_opportunities = .N,
    Eligible_RA_windows = sum(eligible),
    Eligibility_percent = 100 * mean(eligible)
  ), by = grouping]
  setnames(output, grouping, label_name)
  output[, Ineligible_windows := Scheduled_opportunities - Eligible_RA_windows]
  setcolorder(output, c(
    label_name, "Scheduled_opportunities", "Eligible_RA_windows",
    "Ineligible_windows", "Eligibility_percent"
  ))
  output
}

by_map <- rate_table(opportunities, "map_category", "MAP_category_mmhg")
by_time <- rate_table(opportunities, "time_category", "Minutes_since_operation_start")
fwrite(by_map, file.path(output_root, "r_window_observability_by_map.csv"))
fwrite(by_time, file.path(output_root, "r_window_observability_by_time.csv"))

message("[2/4] Classifying mutually exclusive EEG exclusion reasons")
fail <- !logical_true(opportunities$primary_ra_eligible)
reason <- rep("Eligible", nrow(opportunities))
unresolved <- fail

assign_reason <- function(mask, label) {
  selected <- unresolved & mask
  reason[selected] <<- label
  unresolved[selected] <<- FALSE
}

assign_reason(!is.finite(opportunities$eeg_raw_coverage) |
  opportunities$eeg_raw_coverage < 0.80, "Raw coverage <80%")
assign_reason(
  !is.finite(opportunities$eeg_flatline_fraction) |
    !is.finite(opportunities$eeg_flatline_longest_seconds) |
    opportunities$eeg_flatline_fraction >= 0.05 |
    opportunities$eeg_flatline_longest_seconds >= 1.0,
  "Flat waveform criterion"
)
assign_reason(
  !is.finite(opportunities$eeg_saturation_fraction) |
    !is.finite(opportunities$eeg_saturation_longest_seconds) |
    opportunities$eeg_saturation_fraction >= 0.0001 |
    opportunities$eeg_saturation_longest_seconds >= 0.01,
  "Saturation criterion"
)
assign_reason(!is.finite(opportunities$welch_valid_segments) |
  opportunities$welch_valid_segments < 23, "Valid Welch segments <23")
assign_reason(
  !is.finite(opportunities$sqi_coverage) |
    !is.finite(opportunities$sqi_ge80_fraction) |
    opportunities$sqi_coverage < 0.80 |
    opportunities$sqi_ge80_fraction < 0.80,
  "BIS SQI criterion"
)
assign_reason(
  !is.finite(opportunities$sr_coverage) |
    !is.finite(opportunities$sr_lt20_fraction) |
    opportunities$sr_coverage < 0.80 |
    opportunities$sr_lt20_fraction < 0.80,
  "Suppression-ratio criterion"
)
reason[unresolved] <- "Other EEG eligibility failure"
opportunities[, exclusion_reason := factor(
  reason,
  levels = c(
    "Eligible", "Raw coverage <80%", "Flat waveform criterion",
    "Saturation criterion", "Valid Welch segments <23",
    "BIS SQI criterion", "Suppression-ratio criterion",
    "Other EEG eligibility failure"
  )
)]

reason_by_map <- opportunities[, .N, by = .(map_category, exclusion_reason)]
reason_by_map[, Percent_within_MAP_category := 100 * N / sum(N), by = map_category]
setnames(reason_by_map, c("map_category", "N"), c("MAP_category_mmhg", "Windows"))
fwrite(reason_by_map, file.path(output_root, "r_window_exclusion_reasons_by_map.csv"))

message("[3/4] Auditing event-window observability")
event_subwindows <- fread(file.path(input_root, "event_subwindows.csv"))
event_intervals <- fread(file.path(input_root, "event_intervals.csv"))
event_subwindow_summary <- event_subwindows[, .(
  Planned_30s_subwindows = .N,
  Eligible_30s_subwindows = sum(logical_true(primary_ra_eligible)),
  Eligible_30s_subwindows_percent = 100 * mean(logical_true(primary_ra_eligible))
), by = .(event_type, interval)]
event_interval_summary <- event_intervals[, .(
  Events = uniqueN(event_id),
  Eligible_intervals = sum(logical_true(ra_interval_eligible)),
  Eligible_intervals_percent = 100 * mean(logical_true(ra_interval_eligible))
), by = .(event_type, interval)]
event_observability <- merge(
  event_subwindow_summary,
  event_interval_summary,
  by = c("event_type", "interval"),
  all = TRUE,
  sort = FALSE
)
setnames(event_observability, c("event_type", "interval"), c("Event_type", "Interval"))
fwrite(event_observability, file.path(output_root, "r_event_window_observability.csv"))

message("[4/4] Creating Figure S3 and writing the audit summary")
map_plot_data <- copy(by_map)
map_plot_data[, Panel := "MAP category"]
setnames(map_plot_data, "MAP_category_mmhg", "Category")
time_plot_data <- copy(by_time)
time_plot_data[, Panel := "Time since operation start"]
setnames(time_plot_data, "Minutes_since_operation_start", "Category")
observability_plot_data <- rbindlist(list(map_plot_data, time_plot_data), fill = TRUE)

observability_plot <- ggplot(
  observability_plot_data,
  aes(x = Category, y = Eligibility_percent, group = 1)
) +
  geom_line(linewidth = 0.8, colour = "#1F4E79") +
  geom_point(size = 2.4, colour = "#1F4E79") +
  geom_text(aes(label = sprintf("%.1f%%", Eligibility_percent)), vjust = -0.7, size = 3.2) +
  facet_wrap(~Panel, scales = "free_x", ncol = 2) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL, y = "Eligible RA windows (%)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(
  file.path(release_asset_root, "FigureS3_window_observability.png"),
  observability_plot,
  width = 8.4,
  height = 4.2,
  dpi = 300
)
ggsave(
  file.path(release_asset_root, "FigureS3_window_observability.pdf"),
  observability_plot,
  width = 8.4,
  height = 4.2
)

payload <- list(
  status = "window-observability audit",
  subjects = uniqueN(opportunities$subjectid),
  surgeries = uniqueN(opportunities$caseid),
  scheduled_opportunities = nrow(opportunities),
  eligible_windows = sum(opportunities$eligible),
  eligibility_fraction = mean(opportunities$eligible),
  by_map = by_map,
  by_time = by_time,
  event_observability = event_observability
)
write_json_pretty(payload, file.path(output_root, "r_window_observability_audit_results.json"))

print(by_map)
print(by_time)
print(event_observability)
