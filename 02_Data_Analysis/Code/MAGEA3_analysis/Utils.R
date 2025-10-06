calculate_enrichment <- function(n, predictions_df, offtarget_column, K, N) {
  # k = number of known off-targets in top n predictions
  # Since the dataframe is already sorted, we just look at the first n rows
  k <- sum(predictions_df[[offtarget_column]][1:n])

  # Expected by chance
  expected <- n * (K / N)

  # Enrichment factor
  enrichment <- (k / n) / (K / N)

  # Create contingency table for Fisher's exact test
  # |                | Known OT | Not Known OT |
  # |----------------|----------|--------------|
  # | Top n          | k        | n - k        |
  # | Remaining      | K - k    | N - n - K + k|

  contingency_table <- matrix(c(k, n - k,
                                K - k, N - n - (K - k)),
                              nrow = 2, byrow = TRUE)

  # Fisher's exact test
  fisher_result <- fisher.test(contingency_table, alternative = "greater")

  # Hypergeometric test (alternative approach)
  # P(X >= k) where X ~ Hypergeometric(N, K, n)
  hypergeo_pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)

  return(list(
    n = n,
    k = k,
    expected = expected,
    enrichment = enrichment,
    fisher_pval = fisher_result$p.value,
    hypergeo_pval = hypergeo_pval,
    odds_ratio = fisher_result$estimate
  ))
}
