# designSpec(roles=), designRole(), designGroup() -- study roles instead of column names.

st <- list(whole = ~ x)
mk <- function(roles = NULL) designSpec("richness", st, roles = roles)

fullRoles <- list(
  id = "sample",
  group = list(column = "site_type", label = "site type (disturbance)",
               levels = c("a", "b", "c"), reference = "a",
               positions = c(a = 0, b = 50, c = 100), continuous = "site_type_num"),
  strata = "blk", unit = "plt", covariates = c("p_C", "p_N"))

test_that("string group is shorthand; defaults are normalised", {
  ds <- mk(list(id = "sample", group = "site_type"))
  expect_equal(designRole(ds, "id"), "sample")
  expect_equal(designRole(ds, "group"), "site_type")
  expect_identical(designRole(ds, "strata"), NA_character_)
  expect_equal(designGroup(ds)$label, "site_type")
  expect_null(designGroup(ds)$levels)
  expect_equal(names(designGroup(ds)),
               c("column", "label", "levels", "reference", "positions", "continuous"))
})

test_that("full roles resolve, unnamed", {
  ds <- mk(fullRoles)
  expect_identical(designRole(ds, "covariates"), c("p_C", "p_N"))
  expect_identical(designRole(ds, "continuous"), "site_type_num")
  expect_identical(designRole(ds, "strata"), "blk")
  expect_identical(designRole(ds, "unit"), "plt")
  expect_equal(designGroup(ds)$reference, "a")
})

test_that("override wins; NA override defers to the role", {
  ds <- mk(fullRoles)
  expect_identical(designRole(ds, "group", override = "x"), "x")
  expect_identical(designRole(ds, "group", override = NA), "site_type")
  expect_identical(designRole(ds, "group", override = NA_character_), "site_type")
  expect_identical(designRole(ds, "covariates", override = c(a = "u", b = "v")), c("u", "v"))
})

test_that("no roles: NA, override still honoured, NULL ds is fine", {
  ds <- designSpec("richness", list(whole = ~ x))
  expect_false("roles" %in% names(ds))
  expect_identical(designRole(ds, "group"), NA_character_)
  expect_identical(designRole(ds, "group", override = "treatment"), "treatment")
  expect_identical(designRole(NULL, "id"), NA_character_)
  expect_null(designGroup(ds))
  expect_null(designGroup(NULL))
})

test_that("roles = NULL leaves the object exactly as before", {
  expect_identical(mk(NULL), designSpec("richness", st))
})

test_that("invalid roles error, naming the role", {
  expect_error(mk(list(group = "g")), "id")
  expect_error(mk(list(id = "s", group = list(label = "g"))), "group\\$column")
  expect_error(mk(list(id = "s", group = list(column = "g", levels = c("a", "b"),
                                              reference = "z"))), "reference")
  expect_error(mk(list(id = "s", group = list(column = "g", levels = c("a", "b"),
                                              positions = c(a = 1, c = 2)))), "positions")
  expect_error(mk(list(id = "s", group = list(column = "g", positions = c(a = 1)))),
               "positions.*levels")
})

test_that("attachSequencingDepth has no idCol default", {
  core <- data.frame(unique_name = "a")
  expect_error(attachSequencingDepth(core, NULL), "idCol")
})
