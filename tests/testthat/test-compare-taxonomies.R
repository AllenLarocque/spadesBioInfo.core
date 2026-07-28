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

test_that("isRankPlaceholder flags rank-uninformative CALLS", {
  expect_true(isRankPlaceholder("Fusarium_sp"))
  expect_true(isRankPlaceholder("Alveolata_sp"))
  expect_true(isRankPlaceholder("Eukaryota_kgd_Incertae_sedis"))
  expect_true(isRankPlaceholder("Fungi_phy_Incertae_sedis"))
  expect_true(isRankPlaceholder("unidentified"))
})

test_that("isRankPlaceholder does NOT flag real taxa", {
  expect_false(isRankPlaceholder("Fungi"))
  expect_false(isRankPlaceholder("Ascomycota"))
  expect_false(isRankPlaceholder("Didymella_pinodes"))   # a real binomial
  expect_false(isRankPlaceholder("Sanchytriomycota"))
})

test_that("placeholder and abstention are DISJOINT -- they mean different things", {
  # abstention = no call was made. placeholder = a call was made that says nothing at this rank.
  # Conflating them destroys information; that is why this is a sibling, not an extension.
  absent <- c(NA, "", "unassigned", "Unassigned", "unclassified_Fungi")
  place  <- c("Fusarium_sp", "Eukaryota_kgd_Incertae_sedis", "unidentified")

  expect_true(all(isRankAbstention(absent)))
  expect_false(any(isRankPlaceholder(absent)))
  expect_true(all(isRankPlaceholder(place)))
  expect_false(any(isRankAbstention(place)))
})

test_that("isRankPlaceholder handles a SINGLE element and NA without collapsing", {
  expect_length(isRankPlaceholder("Fusarium_sp"), 1L)
  expect_false(isRankPlaceholder(NA_character_))
  expect_length(isRankPlaceholder(character(0)), 0L)
})

test_that("compareTaxonomies reports pct_informative beside pct_assigned", {
  tt <- data.frame(row.names = c("a", "b", "c"),
                   Kingdom = c("Fungi", "Fungi", "unassigned"),
                   Phylum  = c("Ascomycota", "Fungi_phy_Incertae_sedis", "unassigned"),
                   stringsAsFactors = FALSE)

  out <- compareTaxonomies(list(one = tt), ranks = c("Kingdom", "Phylum"))$assignmentRates

  expect_true("pct_informative" %in% names(out))
  # Phylum: 2 of 3 assigned, but only 1 of 3 informative -- the Incertae_sedis says nothing.
  expect_equal(out$pct_assigned[out$rank == "Phylum"], 66.7)
  expect_equal(out$pct_informative[out$rank == "Phylum"], 33.3)
  # pct_assigned is UNCHANGED, so every previously published figure stays comparable.
  expect_equal(out$pct_assigned[out$rank == "Kingdom"], 66.7)
})
