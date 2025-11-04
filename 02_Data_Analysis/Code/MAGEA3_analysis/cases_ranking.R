library(dplyr)

epitox = openxlsx::read.xlsx("../../data/OT_SCORE_full_annotation.xlsx") %>%
  arrange(mismatch)

case_1 = epitox %>%
  filter(Bi_feature_rank == "Low", multibiophys_category == "High") %>%
  arrange(-multibiophys_similarity)

epitox %>%
  filter(Bi_feature_rank == "High", multibiophys_category == "Low") %>%
  nrow()

case_2 = epitox %>%
  filter(Bi_feature_rank == "High", multibiophys_category == "Low")
