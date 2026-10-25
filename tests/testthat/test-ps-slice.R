mkPs <- function(n = 6, s = 5) {
  set.seed(1)
  m <- matrix(rpois(n * s, 20) + 1L, n, s,
              dimnames = list(paste0("t", seq_len(n)), paste0("s", seq_len(s))))
  sd <- data.frame(block = rep(c("a", "b"), length.out = s), trt = seq_len(s),
                   other = letters[seq_len(s)], row.names = colnames(m))
  tt <- matrix(paste0(rep(c("K", "P", "C"), each = n), seq_len(n)), n, 3,
               dimnames = list(rownames(m), c("Kingdom", "Phylum", "Class")))
  tr <- ape::rtree(n, tip.label = rownames(m))
  phyloseq::phyloseq(phyloseq::otu_table(m, taxa_are_rows = TRUE),
                     phyloseq::sample_data(sd), phyloseq::tax_table(tt), tr)
}

test_that("psSlice keeps exactly the declared parts and columns", {
  ps <- mkPs()
  x <- psSlice(ps, "otu")
  expect_null(x@sam_data); expect_null(x@tax_table); expect_null(x@phy_tree)
  x <- psSlice(ps, c("otu", "tree", "tax"), sampleCols = "block", taxCols = "Phylum")
  expect_equal(colnames(phyloseq::sample_data(x)), "block")
  expect_equal(colnames(phyloseq::tax_table(x)), "Phylum")
  expect_false(is.null(x@phy_tree))
  expect_null(psSlice(ps, c("otu", "tax"))@phy_tree)
  expect_true(methods::validObject(psSlice(ps, c("otu", "tree"), sampleCols = c("block", "trt"))))
})

test_that("psSlice stops on an absent column, naming it", {
  ps <- mkPs()
  expect_error(psSlice(ps, "otu", sampleCols = "nope"), "nope")
  expect_error(psSlice(ps, c("otu", "tax"), taxCols = "Order"), "Order")
})

test_that("psSlice does not copy the otu matrix", {
  ps <- mkPs()
  tracemem(ps@otu_table)
  out <- capture.output(x <- psSlice(ps, "otu"))
  untracemem(ps@otu_table)
  expect_false(any(grepl("tracemem", out)))
  expect_identical(x@otu_table, ps@otu_table)
})

test_that("psSlice on lists keeps names, order and NULL", {
  ps <- mkPs()
  x <- psSlice(list(b = ps, a = NULL, c = ps), "otu")
  expect_identical(names(x), c("b", "a", "c"))
  expect_null(x$a)
  expect_null(x$b@sam_data)
})

test_that("taxCols order follows the table, not the request", {
  ps <- mkPs()
  a <- psSlice(ps, c("otu", "tax"), taxCols = c("Class", "Kingdom"))
  b <- psSlice(ps, c("otu", "tax"), taxCols = c("Kingdom", "Class"))
  expect_identical(a, b)
  expect_equal(colnames(phyloseq::tax_table(a)), c("Kingdom", "Class"))
})

test_that("only declared data enters the digest", {
  ps <- mkPs(); h <- function(p, ...) digest::digest(psSlice(p, "otu", sampleCols = "block"))
  ps2 <- ps; ps2@sam_data$other <- rev(ps2@sam_data$other)
  expect_identical(h(ps), h(ps2))
  ps3 <- ps; ps3@sam_data$block <- c("z", "b", "a", "b", "a")
  expect_false(identical(h(ps), h(ps3)))
})

test_that("psRestore equals the unsliced computation", {
  ps <- mkPs(n = 8, s = 6)
  fn1 <- function(p) phyloseq::rarefy_even_depth(p, rngseed = 3, verbose = FALSE)
  fn2 <- function(p) phyloseq::transform_sample_counts(
    phyloseq::prune_samples(c("s4", "s2"), p), function(x) x / sum(x))
  for (fn in list(fn1, fn2)) {
    expect_identical(psRestore(psSlice(ps, "otu"), ps), ps)
    expect_equal(psRestore(fn(psSlice(ps, "otu")), ps), fn(ps))
  }
})

test_that("psRestore re-annotates when result's order differs from full's", {
  ps <- mkPs()
  m <- as(phyloseq::otu_table(ps), "matrix")
  ord <- rev(rownames(m)); so <- c("s3", "s1", "s5")
  res <- phyloseq::otu_table(m[ord, so], taxa_are_rows = TRUE)
  got <- psRestore(res, ps)
  expect_identical(phyloseq::taxa_names(got), ord)
  expect_identical(rownames(phyloseq::tax_table(got)), ord)
  expect_identical(rownames(phyloseq::sample_data(got)), so)
  tt <- phyloseq::tax_table(ps)
  expect_identical(unname(as(phyloseq::tax_table(got), "matrix")), unname(as(tt, "matrix")[ord, ]))
  expect_identical(as.character(phyloseq::sample_data(got)$trt), as.character(phyloseq::sample_data(ps)[so, ]$trt))
})

test_that("psRestore keeps .Cache and works on otu-only full and lists", {
  ps <- mkPs(); cp <- tempfile()
  f <- function(p) p
  r <- reproducible::Cache(f, psSlice(ps, "otu"), cachePath = cp)
  expect_false(is.null(attr(r, ".Cache")))
  expect_false(is.null(attr(psRestore(r, ps), ".Cache")))
  o <- psSlice(ps, "otu")
  expect_s4_class(psRestore(o, o), "phyloseq")
  l <- psRestore(list(a = o, b = NULL), list(a = ps, b = ps))
  expect_identical(names(l), c("a", "b")); expect_null(l$b)
  expect_identical(l$a, ps)
})

test_that("a Cache()d result restores identical to the uncached one", {
  ps <- mkPs(n = 8, s = 6)
  f <- function(p) phyloseq::rarefy_even_depth(p, rngseed = 3, verbose = FALSE)
  r1 <- reproducible::Cache(f, psSlice(ps, "otu"), cachePath = tempfile())
  cp <- tempfile()
  reproducible::Cache(f, psSlice(ps, "otu"), cachePath = cp)
  r2 <- reproducible::Cache(f, psSlice(ps, "otu"), cachePath = cp)
  want <- f(ps)
  expect_identical(dropCacheAttr(psRestore(r1, ps)), dropCacheAttr(psRestore(r2, ps)))
  expect_equal(dropCacheAttr(psRestore(r1, ps)), want)
})

test_that("dropCacheAttr removes .Cache, recursing into plain lists", {
  x <- structure(list(a = structure(1, .Cache = "z")), .Cache = "y")
  y <- dropCacheAttr(x)
  expect_null(attr(y, ".Cache")); expect_null(attr(y$a, ".Cache"))
  expect_identical(y$a, 1)
})

test_that("sliceContract records errors per dropped item", {
  ps <- mkPs()
  spec <- list(parts = c("otu", "tree"), sampleCols = "block")
  reads <- function(p) { stopifnot("block" %in% colnames(phyloseq::sample_data(p)), !is.null(p@phy_tree)); 1 }
  tol <- function(p) 1
  r <- sliceContract(reads, ps, spec)
  expect_setequal(r$dropped, c("part:tree", "sample:block"))
  expect_true(all(r$errored))
  expect_false(any(sliceContract(tol, ps, spec)$errored))
  expect_error(sliceContract(function(p) stop("x"), ps, spec))
})

test_that("sliceContract emits no phantom rows", {
  ps <- mkPs(); id <- function(p) 1
  expect_identical(sliceContract(id, ps, list(parts = c("otu", "tree")))$dropped, "part:tree")
  expect_equal(nrow(sliceContract(id, ps, list(parts = "otu"))), 0L)
})

test_that("wholeScopeResponses is character(0) on an empty design", {
  expect_identical(wholeScopeResponses(list(responses = list())), character(0))
})

test_that("psRestore keeps Cache attributes off the otu_table slot", {
  ps <- mkPs(); cp <- withr::local_tempdir()
  f <- function(p) phyloseq::otu_table(p)
  miss <- psRestore(reproducible::Cache(f, psSlice(ps, "otu"), cachePath = cp), ps)
  hit  <- psRestore(reproducible::Cache(f, psSlice(ps, "otu"), cachePath = cp), ps)
  for (a in c(".Cache", "tags", "callInCache")) {
    expect_null(attr(miss@otu_table, a)); expect_null(attr(hit@otu_table, a))
  }
  expect_identical(digest::digest(miss@otu_table), digest::digest(hit@otu_table))
  expect_false(is.null(attr(hit, ".Cache")))
})
