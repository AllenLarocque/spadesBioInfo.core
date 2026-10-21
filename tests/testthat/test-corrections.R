mkTables <- function() {
  list(veg = data.frame(Plot = c("A", "B", "C"), Block = c(1L, 2L, 2L),
                        UNIQUE_NAME = c("u1", "u2", "u3"), stringsAsFactors = FALSE),
       ea  = data.frame(Plot = c("A", "B"), Block = c(1L, 2L), stringsAsFactors = FALSE))
}
rule <- function(table = "veg", match_column = "Plot", match_value = "B", action = "set",
                 field = "Block", from = "2", to = "3") {
  data.frame(table = table, match_column = match_column, match_value = match_value,
             action = action, field = field, from = from, to = to,
             evidence = "ev", added = "2026-10-20", stringsAsFactors = FALSE)
}

test_that("set changes only the matched rows and keeps the column class", {
  out <- suppressMessages(applyCorrections(mkTables(), rule()))
  expect_equal(out$tables$veg$Block, c(1L, 3L, 2L))
  expect_type(out$tables$veg$Block, "integer")
  expect_equal(out$tables$ea$Block, c(1L, 2L))
  expect_equal(nrow(out$log), 1L)
  expect_equal(out$log$n_rows, 1L)
  expect_named(out$log, c("table", "match_column", "match_value", "action", "field",
                          "from", "to", "n_rows", "evidence", "added"))
  r2 <- rule(field = "UNIQUE_NAME", from = "u2", to = "u2x")
  out2 <- suppressMessages(applyCorrections(mkTables(), r2))
  expect_equal(out2$tables$veg$UNIQUE_NAME, c("u1", "u2x", "u3"))
})

test_that("drop removes the matched rows", {
  out <- suppressMessages(applyCorrections(mkTables(),
           rule(action = "drop", field = NA, from = NA, to = NA)))
  expect_equal(out$tables$veg$Plot, c("A", "C"))
  expect_true(is.na(out$log$field) && is.na(out$log$from) && is.na(out$log$to))
  expect_equal(out$log$n_rows, 1L)
})

test_that("applying the same set twice stops", {
  once <- suppressMessages(applyCorrections(mkTables(), rule()))$tables
  expect_error(suppressMessages(applyCorrections(once, rule())), "expected '2' but found '3'")
})

test_that("a rule matching no rows stops", {
  expect_error(suppressMessages(applyCorrections(mkTables(), rule(match_value = "Z"))),
               "matches no rows")
})

test_that("a missing column or table stops", {
  expect_error(suppressMessages(applyCorrections(mkTables(), rule(match_column = "plot"))),
               "column 'plot' not in table 'veg'")
  expect_error(suppressMessages(applyCorrections(mkTables(), rule(table = "nope"))),
               "table 'nope' not in tables")
})

test_that("an NA current value stops", {
  tb <- mkTables(); tb$veg$Block[2] <- NA
  expect_error(suppressMessages(applyCorrections(tb, rule())), "found 'NA'")
})

test_that("rules apply in order", {
  rules <- rbind(rule(action = "drop", field = NA, from = NA, to = NA),
                 rule(match_value = "C", from = "2", to = "5"))
  out <- suppressMessages(applyCorrections(mkTables(), rules))
  expect_equal(out$tables$veg$Block, c(1L, 5L))
  expect_equal(nrow(out$log), 2L)
})

test_that("readCorrections keeps leading zeros and checks columns", {
  f <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(rule(match_value = "02"), f, row.names = FALSE)
  expect_equal(readCorrections(f)$match_value, "02")
  utils::write.csv(rule()[, -1], f, row.names = FALSE)
  expect_error(readCorrections(f), "lacks column\\(s\\): table")
})

test_that("an unknown action stops", {
  expect_error(suppressMessages(applyCorrections(mkTables(), rule(action = "dorp"))),
               "unknown action 'dorp'")
})

test_that("set can match several rows", {
  out <- suppressMessages(applyCorrections(mkTables(), rule(match_column = "Block", match_value = "2")))
  expect_equal(out$tables$veg$Block, c(1L, 3L, 3L))
  expect_equal(out$log$n_rows, 2L)
})
