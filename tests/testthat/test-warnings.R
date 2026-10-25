f <- function() {
  warning("odd thing"); warning("odd thing"); warning("other thing")
  42
}

test_that("collectWarnings returns the value and summarises each distinct warning once", {
  msgs <- character()
  got <- withCallingHandlers(
    expect_no_warning(collectWarnings(f(), "engineX")),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_identical(got, suppressWarnings(f()))
  expect_length(msgs, 2L)
  expect_match(msgs[1], "^engineX: 2 x warning: odd thing")
  expect_match(msgs[2], "^engineX: 1 x warning: other thing")
})

test_that("collectWarnings is silent without warnings", {
  expect_silent(v <- collectWarnings(list(a = 1), "x"))
  expect_identical(v, list(a = 1))
})

test_that("collectWarnings lets errors through unchanged", {
  err <- tryCatch(collectWarnings(stop("boom"), "x"), error = function(e) e)
  expect_identical(conditionMessage(err), "boom")
})

test_that("collectWarnings handles an empty warning message", {
  msgs <- character()
  v <- withCallingHandlers(collectWarnings({ warning(""); 1 }, "L"),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_identical(v, 1)
  expect_identical(msgs, "L: 1 x warning: \n")
})
