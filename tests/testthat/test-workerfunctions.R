test_that(".is works", {
  expect_equal(dim(.is(rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))),c(380,420))
  expect_equal(class(.is(rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))),c("matrix","array"))
  expect_equal(dim(.is(c(rast(vals=c(1:(180*360))),rast(vals=c(1:(180*360)))))),c(180,360,2))
  expect_equal(class(.is(c(rast(vals=c(1:(180*360))),rast(vals=c(1:(180*360)))))),c("array"))
  expect_equal(.is(c(1,2,3)),c(1,2,3))
  expect_equal(.is(1),1)
  expect_equal(.is(array(1:3, c(2,4))),array(1:3, c(2,4)))
  expect_equal(.is(array(1:3, c(2,4,5))),array(1:3, c(2,4,5)))
  expect_equal(class(.is(array(1:3, c(2,4)))),c("matrix","array"))
  expect_equal(class(.is(array(1:3, c(2,4,5)))),c("array"))
})


test_that(".rta works", {
  expect_equal(dim(.rta(rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")),1)),c(380,420,1))
  expect_equal(dim(.rta(rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")),2)),c(380,420,2))
  expect_equal(class(.rta(rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")),1)),c("array"))

  expect_equal(dim(.rta(c(rast(vals=c(1:(180*360))),rast(vals=c(1:(180*360)))),6)),c(180,360,6))

  expect_equal(.rta(1,3), array(c(1,1,1)))
  expect_equal(class(.rta(1,3)),c("array"))

  expect_equal(.rta(c(1,2,3),2),array(c(1,2)))
  expect_equal(class(.rta(c(1,2,3),2)),c("array"))

  expect_equal(dim(.rta(array(1:3, c(2,4)),3)), c(2,4,3))
  expect_equal(dim(.rta(matrix(1:6, c(2,3)),3)), c(2,3,3))
  expect_equal(dim(.rta(matrix(1:6, c(1,6)),3)),c(1,6,3))
})
# could do with testing for wrong parameter types
test_that(".vta works", {
  expect_equal(class(.vta(c(1,2,3),rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))),c("array"))
  expect_equal(dim(.vta(c(1,2,3),rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))),c(380,420,3))

  expect_equal(dim(.vta(c(1,2,3),matrix(1:6, c(2,3)))), c(2,3,3))
  expect_equal(dim(.vta(c(1,2,3),matrix(1:6, c(1,6)))),c(1,6,3))
})

test_that(".rast works", {
  #  throw error if m provided as vector even if of length ncell(r)...
  expect_error(.rast(seq(1:(380*420)),rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")) ))

  # ...unless configured as matrix or array
  expect_equal(class(.rast( matrix(seq(1:(380*420)),c(380,420)), rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim"))))[1],"SpatRaster" )
  expect_equal(class(.rast( array(seq(1:(380*420)),c(380,420)), rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim"))))[1],"SpatRaster" )

  expect_equal(dim(.rast( matrix(seq(1:(380*420)),c(380,420)), rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))),c(380,420,1) )
  expect_equal(dim(.rast( array(seq(1:(380*420)),c(380,420)), rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))),c(380,420,1) )

  # works if m = multiple of template cells - creating multilayer rast
  expect_equal(dim(.rast( array(seq(1:(380*420*2)),c(380,420,2)), rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))),c(380, 420, 2 ))

  # Warning if dim of matrix does not match dim of template layer
  expect_warning(.rast( matrix(seq(1:(380*410)),c(380,410)), rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim"))))

  # Works if template has multiple layers as long as dim 1 2 match
  expect_equal(dim(.rast( array(seq(1:(380*420)),c(380,420)), c(rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")),rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim"))))),c(380, 420, 1))
  expect_equal(dim(.rast( array(seq(1:(380*420*3)),c(380,420,3)), c(rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")),rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim"))))),c(380, 420, 3))

})

test_that(".ehr works", {
  # Throws error if vector or 1 or 2D array used as parameter
  expect_error(.ehr(c(1,2,3)))
  expect_error(.ehr( array(1:6,c(2,3)) ))
  # Works for 3D numeric incl NAs
  expect_equal(dim(.ehr( array(1:24,c(2,3,4)) )),c(2,3,96))
  expect_equal(dim(.ehr( array(rep(NA,24),c(2,3,4)) )),c(2,3,96))

  # Works if "vector" expressed as 3D array with first 2 dimensions = 1
  expect_equal(dim(.ehr( array(1:24,c(1,6,4)) )),c(1,6,96))
  expect_equal(dim(.ehr( array(1:24,c(1,1,24)) )),c(1,1,576))

  expect_equal(.ehr( array(c(1,2),c(1,1,2)) ),array(c(rep(1,24),rep(2,24)),c(1,1,48))  )
  # Simply replicates values of daily array for each hour in hourly array
  expect_equal( .ehr( array(1:24,c(2,3,4)) )[,,12], array(rep(1:24,24) ,c(2,3,4))[,,1] )
})

test_that(".mav works", {
  # ts class
  expect_equal(class(.mav(c(1:10),1)),"ts")
  # values extracted as.numeric
  expect_equal(as.numeric(.mav(c(1:10),1)),c(1:10))
})

test_that(".latsfromr works", {
  expect_equal(dim(.latsfromr(rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))),c(380,420))
  expect_equal(class(.latsfromr(rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))),c("matrix","array"))
  })

test_that(".lonsfromr works", {
  expect_equal(dim(.lonsfromr(rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))),c(380,420))
  expect_equal(class(.lonsfromr(rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))),c("matrix","array"))
})

test_that(".latslonsfromr works", {
  expect_equal(class(.latslonsfromr(rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))),c("list"))
  expect_equal(length(.latslonsfromr(rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))[[1]]),380*420)
  expect_equal(length(.latslonsfromr(rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))[[2]]),380*420)
  # Works for multilayer rasts - returns lat lon for 1 layer
  expect_equal( length(.latslonsfromr(c(rast(vals=c(1:(180*360))),rast(vals=c(1:(180*360)))))[[1]]),64800)
})

test_that(".resample works", {
  .resample(rast(system.file("extdata/dtms/era5dtm.tif",package="mesoclim")),rast(system.file("extdata/haduk/rainfall1km.tif",package="mesoclim"))[[1]])
  expect_equal(dim(.resample(rast(system.file("extdata/haduk/rainfall1km.tif",package="mesoclim"))[[1]],rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))),c(380,420,1))
  expect_equal(dim(.resample(rast(system.file("extdata/haduk/rainfall1km.tif",package="mesoclim"))[[1:2]],rast(system.file("extdata/dtms/dtmf.tif",package="mesoclim")))),c(380,420,2))
  expect_equal(NA %in% values(.resample(rast(system.file("extdata/dtms/era5dtm.tif",package="mesoclim")),rast(system.file("extdata/haduk/rainfall1km.tif",package="mesoclim"))[[1]],msk=TRUE)),TRUE)
  expect_equal(NA %in% values(.resample(rast(system.file("extdata/dtms/era5dtm.tif",package="mesoclim")),rast(system.file("extdata/haduk/rainfall1km.tif",package="mesoclim"))[[1]],msk=FALSE)),FALSE)
})

test_that(".sea_to_coast works", {

})
