test_that("isRankAbstention recognises every absence token and nothing else", {
  expect_true(all(isRankAbstention(
    c(NA, "", "unassigned", "unclassified_Root", "unclassified_Fungi"))))
  expect_false(any(isRankAbstention(
    c("Fungi", "Ascomycota", "Eukaryota_kgd_Incertae_sedis", "Didymella_sp"))))
})

test_that("mutual abstention is NOT agreement", {
  ours <- data.frame(Phylum = c("Ascomycota", "Basidiomycota", "unassigned"),
                     row.names = c("m1", "m2", "m3"), stringsAsFactors = FALSE)
  oracle <- data.frame(Phylum = c("Ascomycota", "Ascomycota", "unassigned"),
                       row.names = c("m1", "m2", "m3"), stringsAsFactors = FALSE)

  out <- compareTaxonomies(list(ours = ours, oracle = oracle), ranks = "Phylum")
  ag  <- out$pairwiseAgreement

  expect_equal(ag$both_assigned, 2L)   # m3/m3 excluded
  expect_equal(ag$n_agree, 1L)         # only m1
  expect_equal(ag$agreement_given_both_assigned, 0.5)
  expect_equal(ag$both_abstain, 1L)
})

test_that("a vocabulary mismatch between abstention tokens does not count as disagreement", {
  ours   <- data.frame(Kingdom = "unclassified_Root", row.names = "m1",
                       stringsAsFactors = FALSE)
  oracle <- data.frame(Kingdom = "unassigned", row.names = "m1",
                       stringsAsFactors = FALSE)

  ag <- compareTaxonomies(list(ours = ours, oracle = oracle),
                          ranks = "Kingdom")$pairwiseAgreement

  expect_equal(ag$both_abstain, 1L)
  expect_equal(ag$both_assigned, 0L)
  expect_true(is.na(ag$agreement_given_both_assigned))
})

test_that("one-sided abstention is attributed to the correct side and drives recall", {
  ours   <- data.frame(Phylum = c("unassigned", "Ascomycota"),
                       row.names = c("m1", "m2"), stringsAsFactors = FALSE)
  oracle <- data.frame(Phylum = c("Ascomycota", "unassigned"),
                       row.names = c("m1", "m2"), stringsAsFactors = FALSE)

  ag <- compareTaxonomies(list(ours = ours, oracle = oracle),
                          ranks = "Phylum")$pairwiseAgreement

  expect_equal(ag$a_abstain_only, 1L)          # we abstain, oracle calls
  expect_equal(ag$b_abstain_only, 1L)
  expect_equal(ag$recall_of_b_calls, 0)        # oracle called m1; we agreed on 0 of 1
})

test_that("Incertae_sedis is a call, so a confident mis-assignment shows as disagreement", {
  ours   <- data.frame(Kingdom = "Eukaryota_kgd_Incertae_sedis", row.names = "m1",
                       stringsAsFactors = FALSE)
  oracle <- data.frame(Kingdom = "Fungi", row.names = "m1", stringsAsFactors = FALSE)

  out <- compareTaxonomies(list(ours = ours, oracle = oracle), ranks = "Kingdom")

  expect_equal(out$pairwiseAgreement$both_assigned, 1L)
  expect_equal(out$pairwiseAgreement$n_agree, 0L)
  expect_equal(out$disagreements$a_value, "Eukaryota_kgd_Incertae_sedis")
  expect_equal(out$disagreements$b_value, "Fungi")
  expect_equal(out$disagreements$n, 1L)
})

test_that("assignmentRates uses the same predicate", {
  ours <- data.frame(Kingdom = c("Fungi", "unclassified_Root", "unassigned"),
                     row.names = c("m1", "m2", "m3"), stringsAsFactors = FALSE)

  rates <- compareTaxonomies(list(ours = ours), ranks = "Kingdom")$assignmentRates

  expect_equal(rates$n_assigned, 1L)
  expect_equal(rates$n_features, 3L)
  expect_equal(rates$pct_assigned, round(100 / 3, 1))
})

test_that("compareTaxonomies intersects rownames when ids is NULL and honours ids when given", {
  a <- data.frame(Phylum = c("X", "Y"), row.names = c("m1", "m2"),
                  stringsAsFactors = FALSE)
  b <- data.frame(Phylum = c("X", "Z"), row.names = c("m1", "m3"),
                  stringsAsFactors = FALSE)

  ag <- compareTaxonomies(list(a = a, b = b), ranks = "Phylum")$pairwiseAgreement
  expect_equal(ag$n_compared, 1L)          # only m1 is shared

  ag2 <- compareTaxonomies(list(a = a, b = b), ranks = "Phylum",
                           ids = "m1")$pairwiseAgreement
  expect_equal(ag2$n_compared, 1L)
})
