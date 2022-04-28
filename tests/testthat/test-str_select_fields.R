test_that("str_select_fields selects fields from a string", {
  expect_equal(str_select_fields("apple-banana-pear-mango-peach", "-", c(1,5)), c("apple-peach"))
})
