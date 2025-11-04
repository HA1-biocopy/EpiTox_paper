library(UpSetR)
library(ggplot2)
library(tidyr)
library(qpdf)  # install.packages("qpdf")

biocopy_colors = c("#A2C510", "#99CFE9", "#FBB800", "#939597", "#C61E19", "#438D99", "#958BB2", "#6B7B88",
                   "#338232", "#F08000", "#3373A1", "#64686A", "#D14B47", "#98D0BC", "#4F3D7F", "#2C4255")


# Read data
data <- openxlsx::read.xlsx("/Users/hoor.alhasani/Documents/Projects/MAGEA3/P048 binders.xlsx")

# Convert to list format (one vector per cell line)
PEPs_lists <- lapply(data, function(x) x[!is.na(x)])

# Create UpSet plot

cairo_pdf("/Users/hoor.alhasani/Documents/Projects/MAGEA3/AB_OT_overlap.pdf",
          width = 14,
          height = 9,
          family = "Arial",
          onefile = TRUE)

upset(fromList(PEPs_lists),
      nsets = 7,
      order.by = "freq",
      main.bar.color = "#3373A1",
      sets.bar.color = "#A2C510",
      point.size = 4,
      line.size = 2,
      text.scale = c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5))
dev.off()
pdf_subset("/Users/hoor.alhasani/Documents/Projects/MAGEA3/AB_OT_overlap.pdf", pages = 2, output = "/Users/hoor.alhasani/Documents/Projects/MAGEA3/AB_OT_overlap_1.pdf")


# Get all unique peptides
all_peptides <- unique(unlist(lapply(data, function(x) x[!is.na(x)])))

# Calculate statistics
stats <- data.frame(
  CellLine = names(data),
  Total = sapply(data, function(x) sum(!is.na(x))),
  Unique = sapply(1:ncol(data), function(i) {
    others <- unlist(data[-i])
    sum(!data[[i]] %in% others, na.rm = TRUE)
  })
)

stats_long <- pivot_longer(stats, cols = c(Total, Unique))

g = ggplot(stats_long, aes(x = CellLine, y = value, fill = name)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(y = "Number of Peptides",
       x = "Antibody",
       fill = "Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = biocopy_colors)

ggsave(plot = g, "/Users/hoor.alhasani/Documents/Projects/MAGEA3/AB_OT_counts.pdf", width = 12, height = 5, device = cairo_pdf)

####################
# peptide level
####################
# Create a data frame with peptide and cell line information
peptide_antibody <- data.frame()

for(i in 1:ncol(data)) {
  ab <- names(data)[i]
  peptides <- data[[i]][!is.na(data[[i]])]

  if(length(peptides) > 0) {
    temp_df <- data.frame(
      peptide = peptides,
      antibody = ab
    )
    peptide_antibody <- rbind(peptide_antibody, temp_df)
  }
}

# Aggregate to get counts and cell line lists
peptide_summary <- peptide_antibody %>%
  group_by(peptide) %>%
  summarise(
    Count = n(),
    antibodys = paste(antibody, collapse = ", ")
  ) %>%
  arrange(desc(Count))

# View the results
print(head(peptide_summary, 20))

# Save to file if needed
write.csv(peptide_summary, "peptide_summary.csv", row.names = FALSE)

# Plot top 20 with the additional information
top_n <- 20
top_peptides <- peptide_summary[1:top_n,]

ggplot(top_peptides, aes(x = reorder(peptide, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#3373A1") +
  coord_flip() +
  labs(x = "peptide ID",
       y = "Number of antibody",
       title = "Most binding peptides Across antibody") +
  theme_minimal() +
  geom_text(aes(label = Count), hjust = -0.3, size = 3) +
  theme(axis.text.y = element_text(size = 8))
