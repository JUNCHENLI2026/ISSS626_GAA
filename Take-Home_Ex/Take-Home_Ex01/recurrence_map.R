# Presentation of verified results; no records or tests are recomputed.
recurrence_location_map <- function(result, presentation = FALSE) {
  stopifnot(requireNamespace("sf", quietly = TRUE),
    requireNamespace("ggplot2", quietly = TRUE))
  cells <- result$main_grid
  selected <- cells[(cells$gx * 5 == 400 & cells$gy * 5 == 9845) |
    (cells$gx * 5 == 405 & cells$gy * 5 == 9785), ]
  selected <- selected[order(-selected$records), ]
  stopifnot(nrow(selected) == 2L,
    identical(as.integer(selected$records), c(152L, 70L)),
    identical(as.integer(selected$recorded_dates), c(22L, 26L)))
  labels <- data.frame(x = (selected$gx + 0.5) * 5000,
    y = (selected$gy + 0.5) * 5000, label_x = 545000,
    label_y = c(9875000, 9760000),
    text = c("A: Highest record count\nMatan Hilir Utara\n152 records / 22 dates",
      "B: Most recorded dates\nMatan Hilir Selatan\n70 records / 26 dates"))
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = result$districts, fill = "#eef2f4", colour = "#64748b", linewidth = 0.25) +
    ggplot2::geom_sf(data = cells, ggplot2::aes(fill = recorded_dates), colour = NA) +
    ggplot2::geom_sf(data = result$districts, fill = NA, colour = "#64748b", linewidth = 0.2) +
    ggplot2::geom_sf(data = selected, fill = NA, colour = "#111827", linewidth = 1) +
    ggplot2::geom_segment(data = labels,
      ggplot2::aes(x = x, y = y, xend = label_x, yend = label_y), colour = "#111827", linewidth = 0.5) +
    ggplot2::geom_label(data = labels,
      ggplot2::aes(x = label_x, y = label_y, label = text), size = if (presentation) 5 else 3.7,
      fill = "white", colour = "#111827", lineheight = 1.1) +
    ggplot2::scale_fill_viridis_c(option = "C", trans = "sqrt", name = "Recorded\nUTC dates") +
    ggplot2::coord_sf(crs = sf::st_crs(32749), datum = sf::st_crs(4326),
      default_crs = sf::st_crs(32749), xlim = if (presentation) c(360000, 640000) else c(350000, 710000)) +
    ggplot2::labs(title = "Record volume and recurrence identify different cells",
      subtitle = "Ketapang, January-August 2026, confidence >=30, 5 km grid",
      caption = paste("District names describe the modal assignment of records within each highlighted cell.",
        "Grey = no retained records. Observation availability is unknown. Counts do not measure risk.", sep = "\n"),
      x = "Longitude", y = "Latitude") +
    ggplot2::theme_minimal(base_size = if (presentation) 17 else 13) +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 10, hjust = 0),
      panel.grid.minor = ggplot2::element_blank())
}
