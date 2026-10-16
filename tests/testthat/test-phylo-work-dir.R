# The phylo-distance working directory is keyed on the ASV SET, not on the subset label.
#
# ⭐ THE REGRESSION THIS PINS. Both assembly modules used to build the dir as
# `file.path(cachePath(sim), "pd.wd", subsetName)`. `cache/pd.wd/whole/pd.desc` was written on
# 2026-07-12 when `ps_raw` held 25,070 ASVs; Increment 4 replaced them with 16,778 DIFFERENT
# ASVs (different md5 ids). The label "whole" was unchanged, so iCAMP silently reused the
# seventeen-day-old distance matrix ("the pd.wd already has a file named pd.desc, which is
# directly used"), then NST died 24 h into a production run with
# "pd.spname has some OTUs not in community matrix" and "comm.colname has 3148 mismatched names".
#
# A label cannot express which ASVs it stood for, so the fingerprint of the taxa set is what
# names the directory. Same taxa (in any order) -> same dir -> the expensive distance is reused;
# one taxon different -> a different dir, and the stale matrix is unreachable rather than wrong.

taxaSetA <- paste0("ASV", 1:10)
# 0.1.39: a phylo distance dir is named from the taxa AND the tree, so every call passes one.
treeFor <- function(taxa, seed = 1) {
  set.seed(seed); tr <- ape::rtree(length(taxa)); tr$tip.label <- as.character(taxa); tr
}
trA <- treeFor(taxaSetA)
pwd <- function(root, name, taxa, tree = treeFor(taxa)) phyloWorkDir(root, name, taxa, tree = tree)

test_that("phyloWorkDir is deterministic and order-insensitive", {
  expect_identical(pwd("/cache", "whole", taxaSetA),
                   pwd("/cache", "whole", taxaSetA))
  # Order must not matter: phyloseq::taxa_names() order follows the otu_table, which
  # prevalence-filtering and re-assembly can permute without changing the ASV SET.
  expect_identical(pwd("/cache", "whole", taxaSetA),
                   pwd("/cache", "whole", rev(taxaSetA), trA))
  expect_identical(pwd("/cache", "whole", taxaSetA),
                   pwd("/cache", "whole", sample(taxaSetA), trA))
})

test_that("phyloWorkDir discriminates: ONE different taxon means a different directory", {
  changed <- c(taxaSetA[-1], "ASV999")
  expect_false(identical(pwd("/cache", "whole", taxaSetA),
                         pwd("/cache", "whole", changed)))
  # ... and a subset of the same taxa is also a different set.
  expect_false(identical(pwd("/cache", "whole", taxaSetA),
                         pwd("/cache", "whole", taxaSetA[1:9])))
})

test_that("phyloWorkDir keeps subsets apart, as the subset-scoped path always did", {
  expect_false(identical(pwd("/cache", "whole", taxaSetA),
                         pwd("/cache", "saprotroph", taxaSetA)))
  expect_false(identical(pwd("/cacheA", "whole", taxaSetA),
                         pwd("/cacheB", "whole", taxaSetA)))
})

test_that("phyloWorkDir returns <cacheRoot>/pd.wd/<subsetName>-<8 char hash>", {
  p <- pwd("/cache", "whole", taxaSetA)
  expect_identical(dirname(p), file.path("/cache", "pd.wd"))
  expect_match(basename(p), "^whole-[0-9a-f]{8}$")
  # The bare label is exactly the path that broke; it must no longer be produced.
  expect_false(identical(p, file.path("/cache", "pd.wd", "whole")))
})

test_that("phyloWorkDir accepts a phyloseq object and uses its taxa AND its tree", {
  skip_if_not_installed("phyloseq")
  set.seed(42)
  otu <- matrix(sample(0:20, 10 * 4, replace = TRUE), nrow = 10,
                dimnames = list(taxaSetA, paste0("S", 1:4)))
  ps <- phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE), trA)
  expect_identical(phyloWorkDir("/cache", "whole", ps),
                   phyloWorkDir("/cache", "whole", taxaSetA, tree = trA))
  # a phyloseq with no tree has nothing to fingerprint a DISTANCE dir from
  ps0 <- phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE))
  expect_error(phyloWorkDir("/cache", "whole", ps0), "tree")
})

test_that("phyloWorkDir refuses inputs that would silently share a directory", {
  # A zero-taxa set would hash to one constant, so every degenerate subset would collide.
  expect_error(phyloWorkDir("/cache", "whole", character(0), tree = trA), "at least one")
  expect_error(phyloWorkDir("/cache", "whole", NULL, tree = trA), "at least one")
  # A subset name becomes a directory name (the referenceId() lesson).
  expect_error(phyloWorkDir("/cache", "../escape", taxaSetA, tree = trA), "unsafe")
  expect_error(phyloWorkDir("/cache", "a/b", taxaSetA, tree = trA), "unsafe")
  expect_error(phyloWorkDir("/cache", "", taxaSetA, tree = trA), "subsetName")
  # no tree, or a tree over other taxa
  expect_error(phyloWorkDir("/cache", "whole", taxaSetA), "tree")
  expect_error(phyloWorkDir("/cache", "whole", taxaSetA, tree = treeFor(paste0("X", 1:10))),
               "tip labels")
})

test_that("phyloWorkDir's hash reflects CONTENT, not the R session that computed it", {
  # Pinned as LITERALS computed outside R (printf | md5sum), not by re-running the encoder: a hash
  # that changed with the R/digest version, or with an encoder edit, would silently orphan every
  # directory earned by an earlier run.
  #   inner = md5("(4#ASV1:0x1p-1,4#ASV2:0x1p-2)")                        = cfafa42aec6e3c4cdefcd85094189dbf
  #   root  = md5("(4#ASV3:0x1p+0,cfafa42aec6e3c4cdefcd85094189dbf:0x1p-3)") = 659dfec9ab8957e69bfbb4eff01da08a
  #   final = md5("root=659dfec9ab8957e69bfbb4eff01da08a\nroot.edge=none")
  tr <- ape::read.tree(text = "((ASV1:0.5,ASV2:0.25):0.125,ASV3:1);")
  expect_identical(treeFingerprint(tr), "10aca9925a1885c7715d93bd37a79cc0")
  # artifact = md5("tree=10aca9925a1885c7715d93bd37a79cc0\ntipOrder=4#ASV1,4#ASV2,4#ASV3")
  expect_identical(phyloArtifactId(tr), "385d96d80ec91df6eedd8e33262d62fc")
  expect_identical(phyloWorkDir("/c", "whole", c("ASV1", "ASV2", "ASV3"), tree = tr),
                   contentWorkDir("/c", "whole", c("ASV1", "ASV2", "ASV3",
                                                   "tree=385d96d80ec91df6eedd8e33262d62fc")))
})

# ── The tree is part of the identity (0.1.39, Codex review finding 1) ─────────────────────────
test_that("same taxa, different branch lengths -> a different directory", {
  tr2 <- trA; tr2$edge.length[3] <- tr2$edge.length[3] * 1.5
  expect_false(identical(phyloWorkDir("/cache", "whole", taxaSetA, tree = trA),
                         phyloWorkDir("/cache", "whole", taxaSetA, tree = tr2)))
  # a change below write.tree()'s 10 printed digits still separates them
  tr3 <- trA; tr3$edge.length[1] <- tr3$edge.length[1] + 1e-13
  expect_false(identical(treeFingerprint(trA), treeFingerprint(tr3)))
})

test_that("same taxa, different topology -> a different directory", {
  trB <- treeFor(taxaSetA, seed = 2)
  expect_false(identical(phyloWorkDir("/cache", "whole", taxaSetA, tree = trA),
                         phyloWorkDir("/cache", "whole", taxaSetA, tree = trB)))
})

test_that("the fingerprint does not depend on the tree's edge ORDER", {
  pw <- ape::reorder.phylo(trA, "postorder")
  expect_identical(treeFingerprint(pw), treeFingerprint(trA))
})

test_that("the fingerprint is canonical under sibling order; real changes still change it", {
  a <- ape::read.tree(text = "((ASV1:0.5,ASV2:0.25):0.125,(ASV3:1,ASV4:2):0.75);")
  b <- ape::read.tree(text = "((ASV4:2,ASV3:1):0.75,(ASV2:0.25,ASV1:0.5):0.125);")  # siblings swapped
  expect_identical(treeFingerprint(b), treeFingerprint(a))
  # the same swap done on the edge rows (with their lengths) of one stored tree
  sw <- a; i <- which(sw$edge[, 1] == sw$edge[1, 1])        # the root's two child edges
  sw$edge[i, ] <- sw$edge[rev(i), ]; sw$edge.length[i] <- sw$edge.length[rev(i)]
  expect_identical(treeFingerprint(sw), treeFingerprint(a))
  # a topology change, a length change, a label change
  expect_false(identical(treeFingerprint(ape::read.tree(
    text = "((ASV1:0.5,ASV3:0.25):0.125,(ASV2:1,ASV4:2):0.75);")), treeFingerprint(a)))
  expect_false(identical(treeFingerprint(ape::read.tree(
    text = "((ASV1:0.5,ASV2:0.25):0.125,(ASV3:1,ASV4:2.5):0.75);")), treeFingerprint(a)))
  expect_false(identical(treeFingerprint(ape::read.tree(
    text = "((ASV1:0.5,ASV2:0.25):0.125,(ASV3:1,ASV5:2):0.75);")), treeFingerprint(a)))
})

test_that("treeFingerprint scales linearly on a large polytomy", {
  # A star tree puts every tip under the root. Growing the root's child list with c() made this
  # quadratic (measured: 20k tips 1.0 s, 100k tips ~30 s). After the fix: 50k-tip star 0.14 s,
  # 100k-tip star ~0.3 s, 25k-tip random binary tree 0.9 s.
  # The bound is deliberately generous so it is stable on a loaded machine.
  star <- ape::stree(100000); star$edge.length <- rep(1, 100000)
  el <- system.time(fp <- treeFingerprint(star))[["elapsed"]]
  expect_match(fp, "^[0-9a-f]{32}$")
  expect_lt(el, 10)
  set.seed(1); bin <- ape::rtree(25000)
  expect_lt(system.time(treeFingerprint(bin))[["elapsed"]], 20)
})

test_that("phyloArtifactId cannot be fooled by a separator inside a label", {
  # t2 is t1 with tips 1 and 2 RENUMBERED: the same tree (same canonical fingerprint), its
  # labels stored as c("a\ta", "a", "b") instead of c("a", "a\ta", "b"). Tab-joined, both
  # orders read "a<TAB>a<TAB>a<TAB>b".
  t1 <- ape::read.tree(text = "((x:1,y:2):3,z:4);"); t1$tip.label <- c("a", "a\ta", "b")
  t2 <- t1; t2$tip.label <- t1$tip.label[c(2L, 1L, 3L)]
  e <- t2$edge[, 2]; t2$edge[e == 1L, 2] <- 2L; t2$edge[e == 2L, 2] <- 1L
  expect_identical(treeFingerprint(t2), treeFingerprint(t1))
  expect_false(identical(phyloArtifactId(t1), phyloArtifactId(t2)))
})
