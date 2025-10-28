#!/usr/bin/env Rscript
# Improved Sunburst - All layers properly visible
# Key: Aggregate data at each level instead of showing individual targets

set.seed(42)
n_total <- 1400
classifications <- c("Exonic", "Intronic", "Intergenic", "Promoter", "UTR")

# Generate data
data <- data.frame(
  id = 1:n_total,
  classification = sample(classifications, n_total, replace = TRUE,
                          prob = c(0.3, 0.25, 0.2, 0.15, 0.1)),
  risk_score = runif(n_total, 0, 100),
  conservation = runif(n_total, 0, 1),
  has_motif = sample(c("Has_Motif", "No_Motif"), n_total, replace = TRUE, prob = c(0.4, 0.6))
)

data$risk_level <- cut(data$risk_score, breaks = c(0, 40, 70, 100),
                       labels = c("Low_Risk", "Med_Risk", "High_Risk"))
data$selected <- rank(-data$risk_score) <= 350
data$validated <- data$selected & runif(n_total) < 0.8
data$interesting <- data$validated & data$risk_score > 60 &
  data$conservation > 0.5 & runif(n_total) < 0.35

# Print summary
cat("\n=== PIPELINE SUMMARY ===\n")
cat("Total Predictions:", nrow(data), "\n")
cat("Selected:", sum(data$selected), "\n")
cat("Validated:", sum(data$validated), "\n")
cat("Interesting Hits:", sum(data$interesting), "\n\n")

# Function to draw ring segment
draw_segment <- function(r_inner, r_outer, theta_start, theta_end, col,
                         border_col = "white", lwd = 2) {
  if (theta_end <= theta_start) return(invisible())

  n_points <- max(50, ceiling((theta_end - theta_start) * 50))
  theta_seq <- seq(theta_start, theta_end, length.out = n_points)

  x_outer <- r_outer * cos(theta_seq)
  y_outer <- r_outer * sin(theta_seq)
  x_inner <- r_inner * cos(rev(theta_seq))
  y_inner <- r_inner * sin(rev(theta_seq))

  x_all <- c(x_outer, x_inner)
  y_all <- c(y_outer, y_inner)

  polygon(x_all, y_all, col = col, border = border_col, lwd = lwd)
}

# Add text with rotation
add_text <- function(r, theta_start, theta_end, label, cex = 0.7, col = "black") {
  if ((theta_end - theta_start) < 0.15) return(invisible())

  theta_mid <- (theta_start + theta_end) / 2
  x <- r * cos(theta_mid)
  y <- r * sin(theta_mid)

  angle <- theta_mid * 180 / pi - 90
  if (angle < -90) angle <- angle + 180
  if (angle > 90) angle <- angle - 180

  text(x, y, label, cex = cex, srt = angle, col = col, font = 2)
}

# Create comprehensive plot
png("/home/claude/sunburst_improved.png", width = 1600, height = 1600, res = 150)
par(mar = c(1, 1, 4, 1), bg = "white")
plot(0, 0, type = "n", xlim = c(-8, 8), ylim = c(-8, 8),
     asp = 1, axes = FALSE, xlab = "", ylab = "")

title("Off-Target Prediction Multi-Layer Pipeline",
      cex.main = 2, font.main = 2, line = 2)
mtext(sprintf("%d predictions → %d selected → %d validated → %d interesting hits",
              n_total, sum(data$selected), sum(data$validated), sum(data$interesting)),
      cex = 1.3, line = 0.5)

# Define colors
class_colors <- c("Exonic" = "#8DD3C7", "Intronic" = "#FFFFB3",
                  "Intergenic" = "#BEBADA", "Promoter" = "#FB8072",
                  "UTR" = "#80B1D3")
risk_colors <- c("Low_Risk" = "#4DAF4A", "Med_Risk" = "#FF7F00",
                 "High_Risk" = "#E41A1C")
motif_colors <- c("Has_Motif" = "#FDB462", "No_Motif" = "#B3B3B3")

# Layer radii - make them substantial
layers <- list(
  classification = list(r_in = 1, r_out = 2.5),
  risk = list(r_in = 2.7, r_out = 4.0),
  motif = list(r_in = 4.2, r_out = 5.3),
  selected = list(r_in = 5.5, r_out = 6.4),
  validated = list(r_in = 6.6, r_out = 7.5)
)

# LAYER 1: Classification
theta <- 0
for (cls in names(class_colors)) {
  n <- sum(data$classification == cls)
  theta_end <- theta + (n / n_total) * 2 * pi
  draw_segment(layers$classification$r_in, layers$classification$r_out,
               theta, theta_end, class_colors[cls])
  add_text((layers$classification$r_in + layers$classification$r_out) / 2,
           theta, theta_end, paste0(cls, "\n", n), cex = 0.8)
  theta <- theta_end
}

# LAYER 2: Risk Level (within each classification)
for (cls in names(class_colors)) {
  cls_data <- data[data$classification == cls, ]
  n_cls <- nrow(cls_data)

  # Find angular position for this classification
  cls_start_idx <- min(which(data$classification == cls))
  cls_theta_start <- (cls_start_idx - 1) / n_total * 2 * pi

  risk_theta <- cls_theta_start
  for (risk in c("Low_Risk", "Med_Risk", "High_Risk")) {
    n_risk <- sum(cls_data$risk_level == risk)
    if (n_risk > 0) {
      risk_theta_end <- risk_theta + (n_risk / n_total) * 2 * pi
      draw_segment(layers$risk$r_in, layers$risk$r_out,
                   risk_theta, risk_theta_end, risk_colors[risk])
      if (n_risk > 30) {
        add_text((layers$risk$r_in + layers$risk$r_out) / 2,
                 risk_theta, risk_theta_end, n_risk, cex = 0.7)
      }
      risk_theta <- risk_theta_end
    }
  }
}

# LAYER 3: Motif presence (aggregated within classification + risk)
for (cls in names(class_colors)) {
  cls_data <- data[data$classification == cls, ]
  cls_start_idx <- min(which(data$classification == cls))
  cls_theta_start <- (cls_start_idx - 1) / n_total * 2 * pi

  for (risk in c("Low_Risk", "Med_Risk", "High_Risk")) {
    risk_data <- cls_data[cls_data$risk_level == risk, ]
    if (nrow(risk_data) == 0) next

    # Find position
    risk_start_idx <- min(which(data$classification == cls & data$risk_level == risk))
    risk_theta_start <- (risk_start_idx - 1) / n_total * 2 * pi

    motif_theta <- risk_theta_start
    for (motif in c("Has_Motif", "No_Motif")) {
      n_motif <- sum(risk_data$has_motif == motif)
      if (n_motif > 0) {
        motif_theta_end <- motif_theta + (n_motif / n_total) * 2 * pi
        draw_segment(layers$motif$r_in, layers$motif$r_out,
                     motif_theta, motif_theta_end, motif_colors[motif])
        motif_theta <- motif_theta_end
      }
    }
  }
}

# LAYER 4: Selected vs Not Selected (aggregated)
for (cls in names(class_colors)) {
  cls_data <- data[data$classification == cls, ]
  cls_start_idx <- min(which(data$classification == cls))
  cls_theta_start <- (cls_start_idx - 1) / n_total * 2 * pi

  for (risk in c("Low_Risk", "Med_Risk", "High_Risk")) {
    risk_data <- cls_data[cls_data$risk_level == risk, ]
    if (nrow(risk_data) == 0) next

    risk_start_idx <- min(which(data$classification == cls & data$risk_level == risk))
    risk_theta_start <- (risk_start_idx - 1) / n_total * 2 * pi

    for (motif in c("Has_Motif", "No_Motif")) {
      motif_data <- risk_data[risk_data$has_motif == motif, ]
      if (nrow(motif_data) == 0) next

      motif_start_idx <- min(which(data$classification == cls &
                                     data$risk_level == risk &
                                     data$has_motif == motif))
      motif_theta_start <- (motif_start_idx - 1) / n_total * 2 * pi

      sel_theta <- motif_theta_start
      for (selected in c(TRUE, FALSE)) {
        n_sel <- sum(motif_data$selected == selected)
        if (n_sel > 0) {
          sel_theta_end <- sel_theta + (n_sel / n_total) * 2 * pi
          sel_col <- if (selected) "#3498DB" else "#DADADA"
          draw_segment(layers$selected$r_in, layers$selected$r_out,
                       sel_theta, sel_theta_end, sel_col)
          sel_theta <- sel_theta_end
        }
      }
    }
  }
}

# LAYER 5: Validated and Interesting (only for selected)
selected_data <- data[data$selected, ]
for (i in 1:nrow(selected_data)) {
  idx <- which(data$id == selected_data$id[i])
  theta_start <- (idx - 1) / n_total * 2 * pi
  theta_end <- idx / n_total * 2 * pi

  if (selected_data$validated[i]) {
    if (selected_data$interesting[i]) {
      col <- "#FFD700"  # Gold for interesting
    } else {
      col <- "#9B59B6"  # Purple for validated but not interesting
    }
  } else {
    col <- "#95A5A6"  # Gray for not validated
  }

  draw_segment(layers$validated$r_in, layers$validated$r_out,
               theta_start, theta_end, col, border_col = "white", lwd = 0.5)
}

# Add center text
text(0, 0, paste0("All\nPredictions\n", n_total), cex = 1.5, font = 2)

# Add layer labels on the left
labels_x <- -7.5
label_data <- data.frame(
  y = c(7.5, 6.5, 5.5, 4.5, 3.5),
  text = c("Layer 1: Classification",
           "Layer 2: Risk Level",
           "Layer 3: Known Motif",
           "Layer 4: Selected",
           "Layer 5: Validated → Interesting"),
  col = c(NA, NA, NA, "#3498DB", "#FFD700")
)

for (i in 1:nrow(label_data)) {
  if (!is.na(label_data$col[i])) {
    points(labels_x - 0.5, label_data$y[i], pch = 15,
           col = label_data$col[i], cex = 2.5)
  }
  text(labels_x, label_data$y[i], label_data$text[i],
       adj = 0, cex = 0.95, font = 2)
}

# Add legend for colors
legend("bottomright",
       legend = c("Classification:", names(class_colors), "",
                  "Risk Level:", names(risk_colors)),
       fill = c(NA, class_colors, NA, NA, risk_colors),
       border = c(NA, rep("white", 5), NA, NA, rep("white", 3)),
       bty = "n", cex = 0.85)

dev.off()

cat("✅ Improved sunburst saved to: /home/claude/sunburst_improved.png\n")

# Create an alternative: exploded sunburst showing progression
png("/home/claude/sunburst_exploded.png", width = 1800, height = 1200, res = 150)
par(mfrow = c(2, 3), mar = c(2, 2, 3, 2))

# Panel 1: Classification
plot(0, 0, type = "n", xlim = c(-2, 2), ylim = c(-2, 2),
     asp = 1, axes = FALSE, xlab = "", ylab = "")
title("1. Classification", font.main = 2, cex.main = 1.5)
theta <- 0
for (cls in names(class_colors)) {
  n <- sum(data$classification == cls)
  theta_end <- theta + (n / n_total) * 2 * pi
  draw_segment(0.5, 1.8, theta, theta_end, class_colors[cls])
  add_text(1.15, theta, theta_end, paste0(cls, "\n", n))
  theta <- theta_end
}
text(0, 0, paste0(n_total, "\ntargets"), cex = 1.2, font = 2)

# Panel 2: Risk Level
plot(0, 0, type = "n", xlim = c(-2, 2), ylim = c(-2, 2),
     asp = 1, axes = FALSE, xlab = "", ylab = "")
title("2. Risk Scoring", font.main = 2, cex.main = 1.5)
risk_counts <- table(data$risk_level)
theta <- 0
for (risk in names(risk_colors)) {
  n <- risk_counts[risk]
  theta_end <- theta + (n / n_total) * 2 * pi
  draw_segment(0.5, 1.8, theta, theta_end, risk_colors[risk])
  add_text(1.15, theta, theta_end, paste0(risk, "\n", n))
  theta <- theta_end
}

# Panel 3: Motif
plot(0, 0, type = "n", xlim = c(-2, 2), ylim = c(-2, 2),
     asp = 1, axes = FALSE, xlab = "", ylab = "")
title("3. Known Motif", font.main = 2, cex.main = 1.5)
motif_counts <- table(data$has_motif)
theta <- 0
for (motif in names(motif_colors)) {
  n <- motif_counts[motif]
  theta_end <- theta + (n / n_total) * 2 * pi
  draw_segment(0.5, 1.8, theta, theta_end, motif_colors[motif])
  add_text(1.15, theta, theta_end, paste0(motif, "\n", n))
  theta <- theta_end
}

# Panel 4: Selection
plot(0, 0, type = "n", xlim = c(-2, 2), ylim = c(-2, 2),
     asp = 1, axes = FALSE, xlab = "", ylab = "")
title("4. Selection", font.main = 2, cex.main = 1.5)
sel_data <- data.frame(
  label = c("Selected", "Not Selected"),
  n = c(sum(data$selected), sum(!data$selected)),
  col = c("#3498DB", "#DADADA")
)
theta <- 0
for (i in 1:nrow(sel_data)) {
  theta_end <- theta + (sel_data$n[i] / n_total) * 2 * pi
  draw_segment(0.5, 1.8, theta, theta_end, sel_data$col[i])
  add_text(1.15, theta, theta_end, paste0(sel_data$label[i], "\n", sel_data$n[i]))
  theta <- theta_end
}

# Panel 5: Validation
plot(0, 0, type = "n", xlim = c(-2, 2), ylim = c(-2, 2),
     asp = 1, axes = FALSE, xlab = "", ylab = "")
title("5. Validation", font.main = 2, cex.main = 1.5)
val_data <- data.frame(
  label = c("Validated", "Not Validated", "Not Selected"),
  n = c(sum(data$validated), sum(data$selected & !data$validated),
        sum(!data$selected)),
  col = c("#9B59B6", "#95A5A6", "#E8E8E8")
)
theta <- 0
for (i in 1:nrow(val_data)) {
  theta_end <- theta + (val_data$n[i] / n_total) * 2 * pi
  draw_segment(0.5, 1.8, theta, theta_end, val_data$col[i])
  if (val_data$n[i] > 50) {
    add_text(1.15, theta, theta_end, paste0(val_data$label[i], "\n", val_data$n[i]), cex = 0.6)
  }
  theta <- theta_end
}

# Panel 6: Interesting
plot(0, 0, type = "n", xlim = c(-2, 2), ylim = c(-2, 2),
     asp = 1, axes = FALSE, xlab = "", ylab = "")
title("6. Interesting Hits", font.main = 2, cex.main = 1.5)
int_data <- data.frame(
  label = c("⭐ Interesting", "Validated", "Not Validated", "Not Selected"),
  n = c(sum(data$interesting),
        sum(data$validated & !data$interesting),
        sum(data$selected & !data$validated),
        sum(!data$selected)),
  col = c("#FFD700", "#9B59B6", "#95A5A6", "#E8E8E8")
)
theta <- 0
for (i in 1:nrow(int_data)) {
  theta_end <- theta + (int_data$n[i] / n_total) * 2 * pi
  draw_segment(0.5, 1.8, theta, theta_end, int_data$col[i])
  if (int_data$n[i] > 20) {
    add_text(1.15, theta, theta_end, paste0(int_data$label[i], "\n", int_data$n[i]), cex = 0.6)
  }
  theta <- theta_end
}

dev.off()

cat("✅ Exploded view saved to: /home/claude/sunburst_exploded.png\n")
cat("\nAll layers are now clearly visible!\n")
