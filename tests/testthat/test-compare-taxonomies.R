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

# ---------------------------------------------------------------------------
# SILVA. The predicate is SHARED across references, so every UNITE case above is a
# REGRESSION test for what follows: measured 2026-08-13 on the SILVA 138.2 V4 release
# (313,734 records), 78.6% of Species values carry no information and the UNITE-only
# predicate matched 0.3% of them.
#
# The twelve records below are the real fixture, verbatim (prefix-stripped) from
# modules/prepReferenceDb/tests/testthat/fixtures/refdb/silva-taxonomy-sample.tsv.
# They are inlined rather than read from that path so this package's suite stays
# runnable standalone after scripts/publish-standalone.sh.
# ---------------------------------------------------------------------------
.silvaFixture <- function() {
  data.frame(
    row.names = c("CP012870.431955.433510", "CP014004.4881510.4883063", "JQ221335.1.914",
                  "AB905621.1.1268", "KC486078.1.1647", "KF130316.1.1760",
                  "GQ132430.1.1360", "AAQK01003598.434.1912", "HG975445.39219252.39220730",
                  "U29214.1.1374", "AP017652.14792.16335", "AB906274.1.1260"),
    Domain = c("Bacteria", "Bacteria", "Archaea", "Archaea", "Eukaryota", "Eukaryota",
               "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Archaea"),
    Phylum = c("Proteobacteria", "Proteobacteria", "Thermoplasmatota", "Euryarchaeota",
               "Protalveolata", "Protalveolata", "Bacteroidota", "Firmicutes",
               "Cyanobacteria", "Proteobacteria", "Firmicutes", "Thermoplasmatota"),
    Class = c("Gammaproteobacteria", "Gammaproteobacteria", "Thermoplasmata",
              "Methanobacteria", "Colpodellida", "Syndiniales", "Bacteroidia", "Clostridia",
              "Cyanobacteriia", "Alphaproteobacteria", "Bacilli", "Thermoplasmata"),
    Order = c("Enterobacterales", "Enterobacterales", "Thermoplasmata",
              "Methanobacteriales", "Colpodellida", "Syndiniales", "Bacteroidales",
              "Oscillospirales", "Chloroplast", "Rickettsiales", "Lactobacillales",
              "Methanomassiliicoccales"),
    Family = c("Enterobacteriaceae", "Enterobacteriaceae", "Thermoplasmata",
               "Methanobacteriaceae", "Colpodellida", "Syndiniales_Group_II",
               "Rikenellaceae", "Ruminococcaceae", "Chloroplast", "Mitochondria",
               "Streptococcaceae", "Methanomethylophilaceae"),
    Genus = c("Escherichia-Shigella", "Klebsiella", "Marine_Group_III",
              "Methanobrevibacter", "Colpodella", "Syndiniales_Group_II", "uncultured",
              "Subdoligranulum", "Chloroplast", "Mitochondria", "Streptococcus",
              "uncultured"),
    Species = c("Escherichia_coli", "Klebsiella_pneumoniae", "uncultured_marine",
                "uncultured_archaeon", "uncultured_Colpodellidae", "uncultured_eukaryote",
                "uncultured_bacterium", "metagenome", "Solanum_pennellii",
                "Aspergillus_flavus", "Streptococcus_sp.", "uncultured_archaeon"),
    stringsAsFactors = FALSE)
}

.silvaRanks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

test_that("BOTH `_sp` spellings are placeholders -- SILVA's has a trailing period", {
  # ⛔ THE test. SILVA writes `_sp.` in 24,640 V4 records and `_sp` in ZERO; UNITE is the
  # exact opposite. `grepl("_sp$", v)` matched none of SILVA's, silently -- no error, no
  # warning, all 24,640 reported as genuine species calls. Verified to FAIL against that
  # old regex before this line was written.
  expect_true(isRankPlaceholder("Streptococcus_sp."))   # AP017652.14792.16335, verbatim
  expect_true(isRankPlaceholder("Fusarium_sp"))         # UNITE spelling, unchanged
  expect_true(all(isRankPlaceholder(c("Streptococcus_sp.", "Fusarium_sp"))))

  # ...and a real binomial is still a call. The period must not turn `_sp` into "anything
  # ending in an abbreviation".
  expect_false(isRankPlaceholder("Escherichia_coli"))
  expect_false(isRankPlaceholder("Klebsiella_pneumoniae"))
  expect_false(isRankPlaceholder("Didymella_pinodes"))
})

test_that("isRankPlaceholder knows SILVA's environmental vocabulary", {
  expect_true(all(isRankPlaceholder(c(
    "uncultured",                # the bare Genus form: 36,896 records
    "uncultured_bacterium",      # 134,031 records on its own
    "uncultured_archaeon", "uncultured_organism", "uncultured_eukaryote",
    "uncultured_marine",         # JQ221335.1.914
    "uncultured_Colpodellidae",  # KC486078.1.1647
    "metagenome",                # AAQK01003598.434.1912
    "soil_metagenome"))))
})

test_that("`uncultured*` does not swallow genuine values", {
  # Every one of these is a real name from the fixture. A prefix test without a token
  # boundary, or a case-blind `metagenome` substring test, would be the same class of
  # silent error in the OTHER direction -- the one that hid 1,442 of 1,445 Kingdom
  # disagreements.
  expect_false(any(isRankPlaceholder(c(
    "Escherichia_coli", "Klebsiella_pneumoniae", "Marine_Group_III",
    "Syndiniales_Group_II", "Methanobrevibacter", "Solanum_pennellii",
    "Aspergillus_flavus", "Colpodella", "Chloroplast", "Mitochondria",
    "Escherichia-Shigella", "Subdoligranulum"))))
})

test_that("every existing UNITE case still classifies exactly as it did (REGRESSION)", {
  place  <- c("Fusarium_sp", "Alveolata_sp", "Eukaryota_kgd_Incertae_sedis",
              "Fungi_phy_Incertae_sedis", "unidentified")
  call   <- c("Fungi", "Ascomycota", "Didymella_pinodes", "Sanchytriomycota")
  absent <- c(NA, "", "unassigned", "Unassigned", "unclassified_Fungi")

  expect_true(all(isRankPlaceholder(place)))
  expect_false(any(isRankPlaceholder(call)))
  # placeholder and abstention stay DISJOINT after the SILVA vocabulary is added
  expect_false(any(isRankPlaceholder(absent)))
  expect_true(all(isRankAbstention(absent)))
  expect_false(any(isRankAbstention(c(place, call))))
})

test_that("isRankRedundant finds SILVA's propagated parent names", {
  tax <- .silvaFixture()

  # f__Syndiniales_Group_II -> g__Syndiniales_Group_II, plus the two organelle rows.
  genus <- isRankRedundant(tax, "Genus", .silvaRanks)
  expect_equal(rownames(tax)[genus],
               c("KF130316.1.1760", "HG975445.39219252.39220730", "U29214.1.1374"))

  # the c=o=f Thermoplasmata chain: JQ221335 repeats Thermoplasmata at Order AND Family.
  order <- isRankRedundant(tax, "Order", .silvaRanks)
  family <- isRankRedundant(tax, "Family", .silvaRanks)
  jq <- rownames(tax) == "JQ221335.1.914"
  expect_true(order[jq] && family[jq])
  expect_equal(rownames(tax)[order],
               c("JQ221335.1.914", "KC486078.1.1647", "KF130316.1.1760"))
  expect_equal(rownames(tax)[family],
               c("JQ221335.1.914", "KC486078.1.1647", "HG975445.39219252.39220730"))

  # Class is NOT redundant anywhere here: Thermoplasmata != Thermoplasmatota. A predicate
  # doing prefix or partial matching would report these, and be wrong.
  expect_false(any(isRankRedundant(tax, "Class", .silvaRanks)))
  expect_length(genus, nrow(tax))
})

test_that("the root rank is never redundant", {
  tax <- .silvaFixture()
  root <- isRankRedundant(tax, "Domain", .silvaRanks)
  expect_length(root, nrow(tax))
  expect_false(any(root))

  # ...and "root" means first in `ranks`, not a hardcoded name: with Phylum as the root
  # of a 6-rank order, Phylum is all-FALSE and Class is judged against Phylum.
  sub <- .silvaRanks[-1]
  expect_false(any(isRankRedundant(tax, "Phylum", sub)))
  expect_equal(isRankRedundant(tax, "Class", sub), isRankRedundant(tax, "Class", .silvaRanks))
})

test_that("isRankRedundant survives a SINGLE row without collapsing to a vector", {
  # R drops dimensions on single-row subsetting -- a bug class that has produced five
  # defects on this project. One row in, one logical out.
  one <- .silvaFixture()["KF130316.1.1760", , drop = FALSE]
  expect_length(isRankRedundant(one, "Genus", .silvaRanks), 1L)
  expect_true(isRankRedundant(one, "Genus", .silvaRanks))
  expect_length(isRankRedundant(one, "Domain", .silvaRanks), 1L)
  expect_length(isRankUninformative(one, "Genus", .silvaRanks), 1L)

  # a matrix behaves identically, including the 1-row matrix
  m <- as.matrix(.silvaFixture())
  expect_equal(isRankRedundant(m, "Genus", .silvaRanks),
               isRankRedundant(.silvaFixture(), "Genus", .silvaRanks))
  expect_length(isRankRedundant(m["U29214.1.1374", , drop = FALSE], "Genus", .silvaRanks), 1L)

  # and a dimensionless input is an ERROR, not a silent wrong answer
  expect_error(isRankRedundant(.silvaFixture()["KF130316.1.1760", ][["Genus"]], "Genus",
                               .silvaRanks), "drop = FALSE")
})

test_that("absence is not redundancy, and a missing parent column is an error", {
  tax <- data.frame(row.names = c("m1", "m2", "m3"),
                    Family = c("unassigned", NA, "Rikenellaceae"),
                    Genus  = c("unassigned", NA, "Rikenellaceae"),
                    stringsAsFactors = FALSE)
  out <- isRankRedundant(tax, "Genus", c("Family", "Genus"))
  # two abstentions in a row is one classifier declining twice, not a propagated name
  expect_equal(out, c(FALSE, FALSE, TRUE))

  # a parent named in `ranks` but ABSENT from the table stops, rather than quietly
  # returning all-FALSE -- a silent all-FALSE is exactly what `_sp$` did to SILVA.
  expect_error(isRankRedundant(tax, "Family", c("Order", "Family", "Genus")), "no column")
  expect_error(isRankRedundant(tax, "Species", c("Family", "Genus")), "not in .ranks.")
})

test_that("the two mechanisms stay separately callable, and the convenience is their OR", {
  tax <- .silvaFixture()
  # At Genus the fixture has BOTH mechanisms and they are disjoint here: 2 rows say
  # `uncultured` (vocabulary), 3 repeat the Family name (propagation). A caller asking
  # "how much of this reference is genuinely resolved at Genus" must be able to tell
  # those apart -- which is why these are two predicates, not one.
  vocab <- isRankPlaceholder(tax$Genus)
  prop  <- isRankRedundant(tax, "Genus", .silvaRanks)
  expect_equal(sum(vocab), 2L)
  expect_equal(sum(prop), 3L)
  expect_false(any(vocab & prop))

  expect_equal(isRankUninformative(tax, "Genus", .silvaRanks), vocab | prop)
  expect_equal(sum(isRankUninformative(tax, "Genus", .silvaRanks)), 5L)
  # the root rank falls through to the vocabulary alone
  expect_equal(isRankUninformative(tax, "Domain", .silvaRanks),
               isRankPlaceholder(tax$Domain))
})

test_that("the fixture's Species values are dominated by placeholders", {
  # 8 of 12 = 66.7%. The measured figure on the full 313,734-record V4 release is 78.6%.
  # ⚠️ These are NOT the same measurement: the fixture is 12 hand-picked records chosen to
  # cover every FORM, not a sample of the release, so it cannot reproduce the release's
  # proportion and this test does not claim it does. What it pins is that every form the
  # release uses is recognised.
  tax <- .silvaFixture()
  sp <- isRankPlaceholder(tax$Species)
  expect_equal(sum(sp), 8L)
  expect_equal(rownames(tax)[!sp],
               c("CP012870.431955.433510", "CP014004.4881510.4883063",
                 "HG975445.39219252.39220730", "U29214.1.1374"))
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
