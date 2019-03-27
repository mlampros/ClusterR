#============================================================

# data

data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

X = center_scale(dat)


# parameters for the predict_GMM function

CENTROIDS = matrix(runif(84), nrow = 2, ncol = 42)

COVARIANCE = matrix(runif(84), nrow = 2, ncol = 42)

WEIGHTS = c(0.5, 0.5)

#===================================================================================================

# data for the determinant error

gmm_object = structure(list(
  centroids = structure(c(0, 0, 0, 0.03766484839859, 
      0.0171257548711382, 0, 0, 0, 0.03766484839859, 0.0171257548711382, 
      0, 0, 0, 0.03766484839859, 0.0171257548711382, 0, 0, 0, 0.03766484839859, 
      0.0171257548711382, 0, 0, 0, 0.03766484839859, 0.0171257548711382, 
      0, 0, 0, 0.03766484839859, 0.0171257548711382, 0, 0, 0, 0.03766484839859, 
      0.0171257548711382, 0, 0, 0, 0.03766484839859, 0.0171257548711382, 
      0, 0, 0, 0.03766484839859, 0.0171257548711382, 0, 0, 0, 0.036542908854579, 
      0.0156329382385248, 0, 0, 0, 0.03766484839859, 0.0171257548711382, 
      0, 0, 0, 0.03766484839859, 0.0171257548711382, 0, 0, 0, 0.0366401725944694, 
      0.0171257548711382, 0, 0, 0, 0.03766484839859, 0.0171257548711382, 
      0, 0, 0, 0.03766484839859, 0.0171257548711382, 0, 0.0087719298245614, 
      0, 5.96350150593148e-05, 0, 0.00639393939393939, 0.0112781954887218, 
      0, 0.000298388155823344, 0.00284489027558641, 0.00639393939393939, 
      0, 0, 1.1085194972591e-05, 9.60576174798685e-05, 0, 0.216400363768785, 
      0.0204191886018415, 0.0254409452055727, 0.0157856255758543, 0, 
      0.00445558340295182, 0.000977517106549365, 0.00459638760195956, 
      0, 0, 0, 0, 0.03766484839859, 0.0171257548711382, 0, 0, 0, 0.03766484839859, 
      0.0171257548711382, 0, 0, 0, 0.0354457535821118, 0.0168852189682521, 
      0.113824908424908, 0.0805074555074555, 0.0332980115839228, 0.0231918686668483, 
      0.0373579465727526, 0.071493228993229, 0.0210850468287539, 0.00988639117888608, 
      0.00566624349575226, 0.00486928620512207, 0, 0, 0.00134408602150538, 
      3.43353117008176e-05, 0, 0.257077125798178, 0.0704223397587242, 
      0.0347048788962759, 0.00899501647368337, 0.021075630862602, 0.0361256854256854, 
      0.010484544695071, 0.00512524282842558, 0.000639572791353075, 
      0.00161862342059266, 0, 0.0025062656641604, 0.00112581564194467, 
      0.00135010652194564, 0, 0.0875955296165822, 0.0279063987268322, 
      0.020282970359359, 0.00660531341299909, 0.0126483154351555, 0, 
      0.131543279362636, 0.0211046275517559, 0.0103013791992964, 0.0179922835865542, 
      0, 0, 0, 4.33412950977534e-05, 0, 0.00125252525252525, 0, 0, 
      0.00021238564795241, 0.00269837634511686, 0.0420601889338731, 
      0.0825860501254229, 0.0377652869045812, 0.0142784136906197, 0.0254359501799185, 
      0.00943333333333333, 0.00258980785296575, 0.00297238428168219, 
      0.0015681735505349, 0.00128433891082907, 0, 0.107829714842772, 
      0.0862212861111289, 0.0158256727358109, 0.0633946221384378, 0, 
      0, 0.0619933279032652, 0.00680476826267524, 0.0377461116567145, 
      0, 0, 0.0644791399416024, 0.00719747365860426, 0.040396820916943, 
      0, 0, 0.0644791399416024, 0.00719747365860426, 0.040396820916943, 
      0, 0, 0.0644791399416024, 0.00719747365860426, 0.040396820916943, 
      0, 0, 0.0656755679659875, 0.00731168889201231, 0.0407283655482459, 
      0, 0, 0.0656755679659875, 0.00731168889201231, 0.0407283655482459, 
      0, 0, 0.0656755679659875, 0.00731168889201231, 0.0407283655482459, 
      0, 0, 0.0656755679659875, 0.00731168889201231, 0.0407283655482459, 
      0, 0, 0.0192448629877526, 0.00228448046040164, 0.00189443586721685, 
      0, 0, 0.0129785783839116, 0.00075251742753932, 0.000404593248796, 
      0.279785378948537, 0.165990234868954, 0.151650313620455, 0.0598668292019649, 
      0.115740758633808, 0.0885642164852691, 0.0556427892812332, 0.0227655383480004, 
      0.0867324031265648, 0.0464779933786618), .Dim = c(5L, 48L)), 
covariance_matrices = structure(c(1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 1e-10, 1e-10, 8.56801175955926e-05, 
                6.85916841193681e-05, 1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 1e-10, 1e-10, 7.87425302761094e-05, 
                5.83334519802252e-05, 1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 0.00430901815943367, 1e-10, 
                3.40451446706317e-07, 1e-10, 0.00153912213039486, 0.00461680517082179, 
                1e-10, 1.934210624387e-06, 0.000136202995188529, 0.00153912213039486, 
                1e-10, 1e-10, 2.31165209106817e-08, 5.52513387222738e-07, 
                1e-10, 0.0915722402454508, 0.00353033024844606, 0.00177445794711472, 
                0.00122184680336181, 1e-10, 0.000554776588427968, 8.79096518108911e-05, 
                8.73731251795789e-05, 1e-10, 1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 1e-10, 1e-10, 7.11576312759953e-05, 
                5.83334519802253e-05, 1e-10, 1e-10, 1e-10, 9.98516455655598e-05, 
                6.33716901289178e-05, 0.0606842874648943, 0.0198048351145863, 
                0.00260820014558815, 0.00147600317201864, 0.00258713475341651, 
                0.0231256610158258, 0.00301445192491625, 0.000707989274195841, 
                0.000148469777154507, 7.86853943120892e-05, 1e-10, 1e-10, 
                0.000166204185454966, 2.21777655856236e-07, 1e-10, 0.126754365352372, 
                0.0231623519497956, 0.00650227534743119, 0.00067382585501657, 
                0.00135997185429724, 0.016319331293225, 0.00241036313981836, 
                0.000429173195992455, 1.00643085564529e-05, 3.47367283584679e-05, 
                1e-10, 0.000351756584443565, 7.06547402326668e-05, 2.0391125106759e-05, 
                1e-10, 0.0433031467773546, 0.00458975434596521, 0.0014870895804591, 
                0.000159367025544011, 0.000300963796618169, 1e-10, 0.0261732144969538, 
                0.00212468505944598, 0.000376888144590918, 0.000552107366953691, 
                1e-10, 1e-10, 1e-10, 3.53378049333949e-07, 1e-10, 0.000113725538210387, 
                1e-10, 1e-10, 1.97013589531563e-06, 0.000132931480827695, 
                0.00857645455340338, 0.0110432507545232, 0.00431970659624033, 
                0.000362024312361024, 0.000846134496067155, 0.00220670666666667, 
                0.000208513486445158, 0.000197957748910284, 3.73651659263766e-05, 
                2.09003402993198e-05, 1e-10, 0.0390417582094696, 0.00221104565828973, 
                0.000371452120822612, 0.000658514745557353, 1e-10, 1e-10, 
                0.00147390817205996, 0.000124789232851087, 0.000468642963090161, 
                1e-10, 1e-10, 0.00143831981625308, 0.000133032527285851, 
                0.000509941760180745, 1e-10, 1e-10, 0.00143831981625308, 
                0.000133032527285851, 0.000509941760180745, 1e-10, 1e-10, 
                0.00143831981625308, 0.000133032527285851, 0.000509941760180745, 
                1e-10, 1e-10, 0.00146916898823009, 0.000139942152066178, 
                0.000525826304180549, 1e-10, 1e-10, 0.00146916898823009, 
                0.000139942152066178, 0.000525826304180549, 1e-10, 1e-10, 
                0.00146916898823009, 0.000139942152066178, 0.000525826304180549, 
                1e-10, 1e-10, 0.00146916898823009, 0.000139942152066178, 
                0.000525826304180549, 1e-10, 1e-10, 0.00188492104327745, 
                3.64566999222642e-05, 3.12732792872787e-05, 1e-10, 1e-10, 
                0.0013644566462401, 1.45350793180041e-05, 3.94423371612239e-06, 
                0.122895123916256, 0.0467274268024201, 0.029418200125597, 
                0.00481865631102736, 0.012114743213734, 0.0296403461479179, 
                0.0147350697035418, 0.00191844746326634, 0.0022029130747748, 
                0.000891714437826205), .Dim = c(5L, 48L)), weights = c(0.238095238095238, 
                                                                       0.108571428571429, 0.177142857142857, 0.360229256251683, 
                                                                       0.115961219938794)), .Names = c("centroids", "covariance_matrices", 
                                                                                                       "weights"))



data_determinant = matrix(
  c(0, 0.0210526315789474, 0, 0, 0, 0, 0, 
0.0384615384615385, 0.037037037037037, 0.0272727272727273,0, 
0.0210526315789474, 0, 0, 0, 0, 0, 0.0384615384615385, 0.037037037037037, 
0.0272727272727273,0, 0.0210526315789474, 0, 0, 0, 
0, 0, 0.0384615384615385, 0.037037037037037, 0.0272727272727273
,0, 0.0210526315789474, 0, 0, 0, 0, 0, 0.0384615384615385, 
0.037037037037037, 0.0272727272727273,0, 0.0210526315789474, 
0, 0, 0, 0, 0, 0.0384615384615385, 0.037037037037037, 0.0272727272727273
,0, 0.0210526315789474, 0, 0, 0, 0, 0, 0.0384615384615385, 
0.037037037037037, 0.0272727272727273,0, 0.0210526315789474, 
0, 0, 0, 0, 0, 0.0384615384615385, 0.037037037037037, 0.0272727272727273
,0, 0.0210526315789474, 0, 0, 0, 0, 0, 0.0384615384615385, 
0.037037037037037, 0.0272727272727273,0, 0.0210526315789474, 
0, 0, 0, 0, 0, 0.0384615384615385, 0.037037037037037, 0.0272727272727273
,0, 0.0210526315789474, 0, 0, 0, 0, 0, 0.0384615384615385, 
0.037037037037037, 0.0272727272727273,0, 0.0210526315789474, 
0, 0, 0, 0, 0, 0.0384615384615385, 0.037037037037037, 0.0272727272727273,0, 0.0210526315789474, 0, 0, 0, 0, 0, 0.0384615384615385, 
0.037037037037037, 0.0272727272727273,0, 0.0210526315789474, 
0, 0, 0, 0, 0, 0.0384615384615385, 0.037037037037037, 0.0272727272727273,0, 0.0210526315789474, 0, 0, 0, 0, 0, 0.0384615384615385, 
0.037037037037037, 0.0272727272727273,0, 0.0210526315789474, 
0, 0, 0, 0, 0, 0.0384615384615385, 0.037037037037037, 0.0272727272727273
,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 
0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 
0, 0, 0, 0, 0, 0, 0.0384615384615385, 0.037037037037037, 0,0, 0, 0, 0, 0, 0, 0, 0.0192307692307692, 0, 0,0, 
0.0210526315789474, 0, 0, 0, 0, 0, 0.0384615384615385, 0.037037037037037, 
0.0272727272727273,0, 0.0210526315789474, 0, 0, 
0, 0, 0, 0.0384615384615385, 0.037037037037037, 0.0272727272727273
,0, 0.0210526315789474, 0, 0, 0, 0, 0, 0.0384615384615385, 
0.037037037037037, 0.0272727272727273,0, 
0.0210526315789474, 0, 0, 0, 0, 0, 0, 0.0740740740740741, 
0.109090909090909,0, 0, 0, 0, 0, 0.0769230769230769, 
0.0344827586206897, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 
0, 0, 0,0.833333333333333, 0, 0, 0, 0, 0, 
0.275862068965517, 0, 0.0555555555555556, 0,0, 
0, 0, 0, 0, 0, 0.0344827586206897, 0, 0, 0,0, 0, 
0, 0, 0, 0, 0, 0, 0, 0.00909090909090909,0.166666666666667, 
0, 0.0285714285714286, 1, 0.0769230769230769, 0, 0, 0, 0, 
0,0, 0.0105263157894737, 0, 0, 0.461538461538462, 
0, 0.0344827586206897, 0.0384615384615385, 0, 0.00909090909090909,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 
0, 0, 0, 0, 0, 0, 0, 0.00909090909090909,0, 0.0105263157894737, 
0, 0, 0.230769230769231, 0, 0.0344827586206897, 0.0192307692307692, 
0.037037037037037, 0.0454545454545455,0, 0.0105263157894737, 
0, 0, 0, 0, 0, 0, 0, 0.00909090909090909,0, 0.0842105263157895, 
0.114285714285714, 0, 0.230769230769231, 0, 0.0689655172413793, 
0.0384615384615385, 0, 0.00909090909090909,0, 
0.0526315789473684, 0.0857142857142857, 0, 0, 0, 0.0344827586206897, 
0, 0, 0,0, 0.0526315789473684, 0.0857142857142857, 
0, 0, 0, 0.0344827586206897, 0, 0, 0,0, 0.0526315789473684, 
0.0857142857142857, 0, 0, 0, 0.0344827586206897, 0, 0, 0,0, 0.0526315789473684, 0.0857142857142857, 0, 0, 
0, 0.0344827586206897, 0, 0, 0,0, 0.0526315789473684, 
0.0857142857142857, 0, 0, 0, 0.0344827586206897, 0, 0, 0,0, 0.0526315789473684, 0.0857142857142857, 0, 0, 
0, 0.0344827586206897, 0, 0, 0,0, 0.0526315789473684, 
0.0857142857142857, 0, 0, 0, 0.0344827586206897, 0, 0, 0,0, 0.0526315789473684, 0.0857142857142857, 0, 0, 
0, 0.0344827586206897, 0, 0, 0,0, 0, 0, 0, 0, 
0, 0, 0, 0, 0.0181818181818182, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0.0105263157894737, 0.171428571428571, 
0, 0, 0.846153846153846, 0.137931034482759, 0, 0.037037037037037, 
0.190909090909091, 0, 0.0526315789473684, 
0, 0, 0, 0.0769230769230769, 0.103448275862069, 0.153846153846154, 
0.0925925925925926, 0.1), nrow = 10, ncol = 48, byrow = T)

#==========================================================================================================================


context('gaussian mixture models')


##############################
# error handling GMM function
##############################

testthat::test_that("in case that the data is not a matrix or data frame, it returns an error", {
  
  tmp_x = list(X)
  
  testthat::expect_error( GMM(tmp_x, 2, "maha_dist", "random_subset", 10, 10)  )
})


testthat::test_that("in case that the number of gaussian mixture components is less than or equal to 0, it returns an error", {
  
  testthat::expect_error( GMM(X, 0, "maha_dist", "random_subset", 10, 10)  )
})


testthat::test_that("in case that the dist_mode is not one of 'eucl_dist', 'maha_dist', it returns an error", {
  
  testthat::expect_error( GMM(X, 2, "unknown_dist", "random_subset", 10, 10)  )
})


testthat::test_that("in case that the seed_mode is not one of 'static_subset','random_subset','static_spread','random_spread', it returns an error", {
  
  testthat::expect_error( GMM(X, 2, "maha_dist", "unknown_subset", 10, 10)  )
})


testthat::test_that("in case that the km_iter is negative, it returns an error", {
  
  testthat::expect_error( GMM(X, 2, "maha_dist", "random_subset", -1, 10)  )
})


testthat::test_that("in case that the em_iter is negative, it returns an error", {
  
  testthat::expect_error( GMM(X, 2, "maha_dist", "random_subset", 10, -1)  )
})


testthat::test_that("in case that the verbose parameter is not logical, it returns an error", {
  
  testthat::expect_error( GMM(X, 2, "maha_dist", "random_subset", 10, 10, verbose = 'NA')  )
})


testthat::test_that("in case that the var_floor parameter is negative, it returns an error", {
  
  testthat::expect_error( GMM(X, 2, "maha_dist", "random_subset", 10, 10, var_floor = -1)  )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {
  
  tmp_dat = X
  
  tmp_dat[1,1] = NaN

  testthat::expect_error( GMM(tmp_dat, 2, "maha_dist", "random_subset", 10, 10)  )
})


#################
# GMM function
#################


testthat::test_that("in case that the data is a matrix the result is a list of length 4 and the class is 'Gaussian Mixture Models' ", {
  
  res = GMM(X, 2, "maha_dist", "random_subset", 10, 10)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == 4 && sum(names(res) %in% c("centroids", "covariance_matrices", "weights", "Log_likelihood")) == 4 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 4 && class(res) == "Gaussian Mixture Models" )
  }
})


testthat::test_that("in case that the data is a data frame the result is a list of length 4 and the class is 'Gaussian Mixture Models' ", {
  
  res = GMM(dat, 2, "eucl_dist", "static_subset", 10, 10)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == 4 && sum(names(res) %in% c("centroids", "covariance_matrices", "weights", "Log_likelihood")) == 4 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 4 && class(res) == "Gaussian Mixture Models" )
  }
})


testthat::test_that("in case that the data is a matrix the result is a list of length 4 and the class is 'Gaussian Mixture Models' ", {
  
  res = GMM(X, 2, "maha_dist", "random_subset", 10, 10)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == 4 && sum(names(res) %in% c("centroids", "covariance_matrices", "weights", "Log_likelihood")) == 4 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 4 && class(res) == "Gaussian Mixture Models" )
  }
})


testthat::test_that("in case that the data is a data frame the result is a list of length 4 and the class is 'Gaussian Mixture Models' ", {
  
  res = GMM(dat, 2, "eucl_dist", "static_subset", 10, 10)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == 4 && sum(names(res) %in% c("centroids", "covariance_matrices", "weights", "Log_likelihood")) == 4 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 4 && class(res) == "Gaussian Mixture Models" )
  }
})



testthat::test_that("in case that the data is a matrix the result is a list of length 4 and the class is 'Gaussian Mixture Models' ", {
  
  res = GMM(X, 2, "maha_dist", "random_spread", 10, 10)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == 4 && sum(names(res) %in% c("centroids", "covariance_matrices", "weights", "Log_likelihood")) == 4 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 4 && class(res) == "Gaussian Mixture Models" )
  }
})


testthat::test_that("in case that the data is a data frame the result is a list of length 4 and the class is 'Gaussian Mixture Models' ", {
  
  res = GMM(dat, 2, "eucl_dist", "static_spread", 10, 10)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == 4 && sum(names(res) %in% c("centroids", "covariance_matrices", "weights", "Log_likelihood")) == 4 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 4 && class(res) == "Gaussian Mixture Models" )
  }
})



######################################
# error handling predict_GMM function
######################################


testthat::test_that("in case that the data is not a matrix or data frame, it returns an error", {
  
  tmp_x = list(X)
  
  testthat::expect_error( predict_GMM(tmp_x, CENTROIDS, COVARIANCE, WEIGHTS)  )
})


testthat::test_that("in case that the data is a matrix AND the CENTROIDS is NOT a matrix or data frame, it returns an error", {
  
  tmp_c = list(CENTROIDS)

  testthat::expect_error( predict_GMM(X, tmp_c, COVARIANCE, WEIGHTS) )
})


testthat::test_that("in case that the data is a matrix AND the COVARIANCE is NOT a matrix or data frame, it returns an error", {
  
  tmp_c = list(COVARIANCE)
  
  testthat::expect_error( predict_GMM(X, CENTROIDS, tmp_c, WEIGHTS) )
})


testthat::test_that("in case that the number of columns of the data do not equal the number of columns of the COVARIANCE or CENTROIDS matrices, it returns an error", {

  testthat::expect_error( predict_GMM(X[, -1], CENTROIDS, COVARIANCE, WEIGHTS) )
})


testthat::test_that("in case that the length of WEIGHTS vector does not equal the number of rows of the COVARIANCE or CENTROIDS matrices, it returns an error", {
  
  testthat::expect_error( predict_GMM(X, CENTROIDS, COVARIANCE, WEIGHTS[-1]) )
})


testthat::test_that("in case that the class of the WEIGHTS vector is not numeric, it returns an error", {
  
  tmp_w = matrix(WEIGHTS, nrow = 1)
  
  testthat::expect_error( predict_GMM(X, CENTROIDS, COVARIANCE, tmp_w) )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {
  
  tmp_dat = X
  
  tmp_dat[1,1] = Inf
  
  testthat::expect_error( predict_GMM(tmp_dat, CENTROIDS, COVARIANCE, WEIGHTS) )
})



testthat::test_that("in case that the determinant is zero the function returns an error", {

  testthat::expect_error( predict_GMM(data_determinant, gmm_object$centroids, gmm_object$covariance_matrices, gmm_object$weights) )
})


#######################
# predict_GMM function
#######################


testthat::test_that("in case that the data is a matrix the result is a list of length 3 and the class is 'Gaussian Mixture Models' ", {
  
  res = predict_GMM(X, CENTROIDS, COVARIANCE, WEIGHTS)
  
  testthat::expect_true( length(res) == 3 && sum(names(res) %in% c("log_likelihood", "cluster_proba", "cluster_labels")) == 3 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 3 && class(res) == "Gaussian Mixture Models")
})


testthat::test_that("in case that the data is a data frame the result is a list of length 3 and the class is 'Gaussian Mixture Models' ", {
  
  res = predict_GMM(dat, CENTROIDS, COVARIANCE, WEIGHTS)
  
  testthat::expect_true( length(res) == 3 && sum(names(res) %in% c("log_likelihood", "cluster_proba", "cluster_labels")) == 3 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 3 && class(res) == "Gaussian Mixture Models" )
})


testthat::test_that("in case that the data is a matrix AND the CENTROIDS is a data frame the result is a list of length 3 and the class is 'Gaussian Mixture Models' ", {
  
  tmp_c = data.frame(CENTROIDS)
  
  res = predict_GMM(X, tmp_c, COVARIANCE, WEIGHTS)
  
  testthat::expect_true( length(res) == 3 && sum(names(res) %in% c("log_likelihood", "cluster_proba", "cluster_labels")) == 3 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 3 && class(res) == "Gaussian Mixture Models")
})


testthat::test_that("in case that the data is a matrix AND the COVARIANCE is a data frame the result is a list of length 3 and the class is 'Gaussian Mixture Models' ", {
  
  tmp_c = data.frame(COVARIANCE)
  
  res = predict_GMM(X, CENTROIDS, tmp_c, WEIGHTS)
  
  testthat::expect_true( length(res) == 3 && sum(names(res) %in% c("log_likelihood", "cluster_proba", "cluster_labels")) == 3 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 3 && class(res) == "Gaussian Mixture Models")
})



###############################################
# error handling Optimal_Clusters_GMM function
###############################################


testthat::test_that("in case that the data is not a matrix or a data frame, it returns an error", {
  
  tmp_x = list(X)
  
  testthat::expect_error( Optimal_Clusters_GMM(tmp_x, Nr_clusters, criterion = "BIC", plot_data = FALSE) )
})


testthat::test_that("in case that the max_clusters parameter is less than 2 and plot_data is TRUE, it returns an error", {

  testthat::expect_error( Optimal_Clusters_GMM(X, 1, criterion = "BIC", plot_data = T) )
})


testthat::test_that("in case that the max_clusters parameter is not numeric, it returns an error", {
  
  tmp_m = data.frame(1)
  
  testthat::expect_error( Optimal_Clusters_GMM(X, tmp_m, criterion = "BIC", plot_data = F) )
})


testthat::test_that("in case that the max_clusters parameter is not a numeric vector (of length 1 or greater), it returns an error", {
  
  tmp_m = list(1,2)
  
  testthat::expect_error( Optimal_Clusters_GMM(X, tmp_m, criterion = "BIC", plot_data = F) )
})


testthat::test_that("in case that the max_clusters parameter is less than 2 AND plot_data is TRUE, it returns an error", {
  
  tmp_m = 1
  
  testthat::expect_error( Optimal_Clusters_GMM(X, tmp_m, criterion = "BIC", plot_data = T) )
})


testthat::test_that("in case that the criterion parameter is not 'AIC', 'BIC', it returns an error", {

  testthat::expect_error( Optimal_Clusters_GMM(X, 5, criterion = "invalid", plot_data = F) )
})



testthat::test_that("in case that the dist_mode parameter is not 'eucl_dist', 'maha_dist', it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_GMM(X, 5, dist_mode = 'invalid', criterion = "BIC", plot_data = F) )
})



testthat::test_that("in case that the seed_mode parameter is not 'static_subset','random_subset','static_spread','random_spread', it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_GMM(X, 5, dist_mode = 'eucl_dist', criterion = "BIC", seed_mode = 'invalid', plot_data = F) )
})



testthat::test_that("in case that the km_iter is less than 0, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_GMM(X, 5, dist_mode = 'eucl_dist', criterion = "BIC", seed_mode = 'random_subset', km_iter = -1, plot_data = F) )
})



testthat::test_that("in case that the em_iter is less than 0, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_GMM(X, 5, dist_mode = 'eucl_dist', criterion = "BIC", seed_mode = 'random_subset', km_iter = 5, 
                                               
                                               em_iter = -1, plot_data = F) )
})



testthat::test_that("in case that the verbose parameter is not logical, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_GMM(X, 5, dist_mode = 'eucl_dist', criterion = "BIC", seed_mode = 'random_subset', km_iter = 5, 
                                               
                                               em_iter = 5, plot_data = F, verbose = 'invalid') )
})


testthat::test_that("in case that the verbose parameter is not logical, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_GMM(X, 5, dist_mode = 'eucl_dist', criterion = "BIC", seed_mode = 'random_subset', km_iter = 5, 
                                               
                                               em_iter = 5, plot_data = F, var_floor = -1) )
})


# testthat::test_that("in case that max_clusters is greater than the number of the columns and verbose is TRUE, it returns a warning", {
#   
#   testthat::expect_warning( Optimal_Clusters_GMM(X, ncol(X) + 1, dist_mode = 'eucl_dist', criterion = "BIC", seed_mode = 'random_subset', km_iter = 5, 
#                                                  
#                                                  em_iter = 5, plot_data = F, verbose = T) )
# })


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {
  
  tmp_dat = X
  
  tmp_dat[1,1] = -Inf
  
  Nr_clusters = 5
  
  testthat::expect_error( Optimal_Clusters_GMM(tmp_dat, Nr_clusters, criterion = "BIC", plot_data = FALSE) )
})


testthat::test_that("in case that 0 is in the vector of clusters, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_GMM(X, 0:3, criterion = "BIC", plot_data = F) )
})



################################
# Optimal_Clusters_GMM function
################################


testthat::test_that("in case that the data is a matrix the result is a vector and the class is 'Gaussian Mixture Models' ", {

  Nr_clusters = 5

  res = Optimal_Clusters_GMM(X, Nr_clusters, criterion = "BIC", plot_data = FALSE)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == Nr_clusters && class(res) == "Gaussian Mixture Models" )
  }
})



testthat::test_that("in case that the data is a data frame the result is a vector and the class is 'Gaussian Mixture Models' ", {

  Nr_clusters = 5

  res = Optimal_Clusters_GMM(dat, Nr_clusters, criterion = "BIC", plot_data = T)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == Nr_clusters && class(res) == "Gaussian Mixture Models" )
  }
})



testthat::test_that("in case of different parameters the result is a vector and the class is 'Gaussian Mixture Models' ", {

  Nr_clusters = 5

  res = Optimal_Clusters_GMM(dat, Nr_clusters, criterion = "AIC", dist_mode = 'maha_dist', seed_mode = 'static_spread', plot_data = T)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == Nr_clusters && class(res) == "Gaussian Mixture Models" )
  }
})



testthat::test_that("in case of a contiguous-vector it returns the correct output ", {
  
  Nr_clusters = 2:6
  
  res = Optimal_Clusters_GMM(dat, Nr_clusters, criterion = "AIC", dist_mode = 'maha_dist', seed_mode = 'static_spread', plot_data = T)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == length(Nr_clusters) && class(res) == "Gaussian Mixture Models" )
  }
})



testthat::test_that("in case of a non-contiguous-vector it returns the correct output ", {
  
  Nr_clusters = c(2,4,6)
  
  res = Optimal_Clusters_GMM(dat, Nr_clusters, criterion = "AIC", dist_mode = 'maha_dist', seed_mode = 'static_spread', plot_data = T)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == length(Nr_clusters) && class(res) == "Gaussian Mixture Models" )
  }
})
