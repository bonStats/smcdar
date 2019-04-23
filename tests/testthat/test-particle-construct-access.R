context("Particles object")

# setup
alpha_vec <- 10:12

beta_mat <- matrix(1:9, nrow = 3, ncol = 3)

kappa_array <- array(13:24, dim = c(3,2,2))
kappa_mat <- matrix(kappa_array, nrow = 3)

pt <- particles(alpha = alpha_vec, beta = beta_mat, kappa = kappa_array)


test_that("Construction correct",{

  expect_true( all( pt[] == cbind(alpha_vec, beta_mat, kappa_mat) ) )

  expect_equal(pt[["alpha"]], alpha_vec)

  expect_equal(pt[["beta"]], beta_mat)

  expect_equal(pt[["kappa"]], kappa_array)

  })


test_that("Reassignment correct",{


  pt[["alpha"]] <- pt[["alpha"]] + 100

  pt[["beta"]] <- pt[["beta"]] + 100

  pt[["kappa"]] <- pt[["kappa"]] + 100

  expect_equal(pt[["alpha"]], alpha_vec + 100)

  expect_equal(pt[["beta"]], beta_mat + 100)

  expect_equal(pt[["kappa"]], kappa_array + 100)


})

