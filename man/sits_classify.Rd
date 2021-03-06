% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sits_classification.R
\name{sits_classify}
\alias{sits_classify}
\alias{sits_classify.sits}
\alias{sits_classify.raster_cube}
\title{Classify time series or data cube using machine learning models}
\usage{
sits_classify(data, ml_model, ...)

\method{sits_classify}{sits}(data, ml_model, ..., filter_fn = NULL, multicores = 2)

\method{sits_classify}{raster_cube}(
  data,
  ml_model,
  ...,
  roi = NULL,
  filter_fn = NULL,
  impute_fn = sits_impute_linear(),
  interp_fn = NULL,
  compose_fn = NULL,
  start_date = NULL,
  end_date = NULL,
  memsize = 8,
  multicores = 2,
  output_dir = tempdir(),
  version = "v1",
  verbose = FALSE
)
}
\arguments{
\item{data}{data cube}

\item{ml_model}{R model trained by \code{\link[sits]{sits_train}}.}

\item{...}{other parameters to be passed to specific functions}

\item{filter_fn}{smoothing filter to be applied (if desired).}

\item{multicores}{number of cores to be used for classification.}

\item{roi}{a region of interest (see above)}

\item{impute_fn}{impute function to replace NA}

\item{interp_fn}{function to interpolate points from cube to match
samples}

\item{compose_fn}{function to compose points from cube to match
samples}

\item{start_date}{starting date for the classification}

\item{end_date}{end date for the classification}

\item{memsize}{memory available for classification (in GB).}

\item{output_dir}{directory for output file}

\item{version}{version of the output (for multiple
classifications)}

\item{verbose}{print information about processing time?}
}
\value{
Predicted data (classified time series)
                          or a data cube with probabilities for each class.
}
\description{
This function classifies a set of time series or data cube given
a trained model prediction model created by \code{\link[sits]{sits_train}}.

SITS supports the following models:
\itemize{
 \item{support vector machines: } {see \code{\link[sits]{sits_svm}}}
 \item{random forests: }          {see \code{\link[sits]{sits_rfor}}}
 \item{linear discriminant analysis: } {see \code{\link[sits]{sits_lda}}}
 \item{quadratic discriminant analysis: } {see \code{\link[sits]{sits_qda}}}
 \item{multinomial logit: }       {see \code{\link[sits]{sits_mlr}}}
 \item{extreme gradient boosting: } {see \code{\link[sits]{sits_xgboost}}}
 \item{multi-layer perceptrons: } {see \code{\link[sits]{sits_mlp}}}
 \item{mixed 1D CNN and MLP networks: } {see \code{\link[sits]{sits_TempCNN}}}
 \item{deep residual netwroks:}{see \code{\link[sits]{sits_ResNet}}}
 }
}
\note{
The "roi" parameter defines a region of interest. It can be
   an sf_object, a shapefile, or a bounding box vector with
   named XY values ("xmin", "xmax", "ymin", "ymax") or
   named lat/long values ("lat_min", "lat_max", "lon_min", "lon_max")

   The "filter_fn" parameter specifies a smoothing filter to be applied to
   time series for reducing noise. Currently, options include
   Savtizky-Golay (see \code{\link[sits]{sits_sgolay}}) and Whittaker
   (see \code{\link[sits]{sits_whittaker}}).

   The "impute_fn" function is used to remove invalid or cloudy pixels
   from time series. The default is a linear interpolator, available
   in \code{\link[sits]{sits_impute_linear}}. Users can add their custom
   functions.

   The "interp_fn" function is used when the training samples which
   were used to generate the classification model have a larger number of
   time instances than the data cube. In this case, pixel time series
   of the data cube will have to be interpolated to fit that of the samples.

   The "compose_fn" function is used when the training samples which
   were used to generate the classification model have a small number of
   time instances than the data cube. In this case, pixel time series
   of the data cube will have to be composed and/or merged
   to fit that of the samples.

   The "memsize" and "multicores" parameters are used for multiprocessing.
   The "multicores" parameter defines the number of cores used for
   processing. The "memsize" parameter  controls the amount of memory
   available for classification.
}
\examples{
\donttest{
# Example of classification of a time series
# Retrieve the samples for Mato Grosso
# select an extreme gradient boosting model
samples_2bands <- sits_select(samples_modis_4bands,
                            bands = c("NDVI", "EVI"))
xgb_model <- sits_train(samples_2bands,
    ml_method = sits_xgboost(verbose = FALSE)
)
# classify the point
point_2bands <- sits_select(point_mt_6bands,
                            bands = c("NDVI", "EVI"))
point_class <- sits_classify(point_2bands, xgb_model)
plot(point_class)

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

# classify the raster image
probs_cube <- sits_classify(cube,
    ml_model = xgb_model,
    output_dir = tempdir(),
    memsize = 4, multicores = 2
)

# label the classified image
label_cube <- sits_label_classification(probs_cube, output_dir = tempdir())
plot(label_cube)
}

}
\author{
Rolf Simoes, \email{rolf.simoes@inpe.br}

Gilberto Camara, \email{gilberto.camara@inpe.br}
}
