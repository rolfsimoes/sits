% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sits_accuracy.R
\name{sits_accuracy}
\alias{sits_accuracy}
\alias{sits_accuracy.sits}
\alias{sits_accuracy.classified_image}
\title{Area-weighted classification accuracy assessment}
\usage{
sits_accuracy(data, ...)

\method{sits_accuracy}{sits}(data, ...)

\method{sits_accuracy}{classified_image}(data, ..., validation_csv)
}
\arguments{
\item{data}{Either a data cube with classified images or
a set of time series}

\item{\dots}{Specific parameters}

\item{validation_csv}{A CSV file path with validation data}
}
\value{
A list of lists: The error_matrix, the class_areas, the unbiased
estimated areas, the standard error areas, confidence interval 95% areas,
and the accuracy (user, producer, and overall), or NULL if the data is empty.
A confusion matrix assessment produced by the caret package.
}
\description{
This function calculates the accuracy of the classification
result. For a set of time series, it creates a confusion matrix and then
calculates the resulting statistics using the R packge "caret". The time
series needs to be classified using \code{\link[sits]{sits_classify}}.

Classified images are generated using \code{\link[sits]{sits_classify}}
followed by \code{\link[sits]{sits_label_classification}}.
For a classified image, the function uses an area-weighted technique
proposed by Olofsson et al. according to [1-3] to produce more reliable
accuracy estimates at 95% confidence level.

In both cases, it provides an accuracy assessment of the classified,
including Overall Accuracy, Kappa, User's Accuracy, Producer's Accuracy
and error matrix (confusion matrix)
}
\examples{
\donttest{
# Case (1) - Accuracy for classification of time series
# read a tibble with 400 time series of Cerrado and 346 of Pasture
data(cerrado_2classes)
# create a model for classification of time series
svm_model <- sits_train(cerrado_2classes, sits_svm())
# classify the time series
predicted <- sits_classify(cerrado_2classes, svm_model)
# calculate the classification accuracy
acc <- sits_accuracy(predicted)

# Case (2) - Accuracy for classification of raster data
# select a training set with two bands
samples_modis_2bands <- sits_select(samples_modis_4bands,
                                    bands = c("NDVI", "EVI"))

# build an extreme gradient boosting model
xgb_model <- sits_train(
  samples_modis_2bands,
  sits_xgboost(nrounds = 50, verbose = FALSE)
)

# create a data cube based on files
data_dir <- system.file("extdata/raster/mod13q1", package = "sits")
cube <- sits_cube(
    source = "LOCAL",
    name = "sinop-2014",
    satellite = "TERRA",
    sensor = "MODIS",
    data_dir = data_dir,
    delim = "_",
    parse_info = c("X1", "X2", "tile", "band", "date")
)


# classify the data cube with xgb model
probs_cube <- sits_classify(cube,
  xgb_model,
  output_dir = tempdir(),
  memsize = 4,
  multicores = 2
)
# label the classification
label_cube <- sits_label_classification(probs_cube,
  output_dir = tempdir()
)
# get ground truth points
ground_truth <- system.file("extdata/samples/samples_sinop_crop.csv",
  package = "sits"
)
# calculate accuracy according to Olofsson's method
accuracy <- suppressWarnings(sits_accuracy(label_cube,
            validation_csv = ground_truth))
}
}
\references{
[1] Olofsson, P., Foody, G.M., Stehman, S.V., Woodcock, C.E. (2013).
Making better use of accuracy data in land change studies: Estimating
accuracy and area and quantifying uncertainty using stratified estimation.
Remote Sensing of Environment, 129, pp.122-131.

[2] Olofsson, P., Foody G.M., Herold M., Stehman, S.V.,
Woodcock, C.E., Wulder, M.A. (2014)
Good practices for estimating area and assessing accuracy of land change.
Remote Sensing of Environment, 148, pp. 42-57.

[3] FAO, Map Accuracy Assessment and Area Estimation: A Practical Guide.
National forest monitoring assessment working paper No.46/E, 2016.
}
\author{
Rolf Simoes, \email{rolf.simoes@inpe.br}

Alber Sanchez, \email{alber.ipia@inpe.br}
}
