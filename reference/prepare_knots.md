# Prepare knot locations for Hmsc GPP models

Prepare the locations of knots for use in Gaussian Predictive Process
(GPP) models within the HMSC framework. It ensures that knots are spaced
at a minimum specified distance and applies jitter to any identical
coordinates to avoid overlap.

## Usage

``` r
prepare_knots(
  coordinates = NULL,
  min_distance = NULL,
  jitter_distance = 100,
  min_lf = NULL,
  max_lf = NULL,
  alphapw = list(Prior = NULL, Min = 10, Max = 1500, Samples = 101)
)
```

## Arguments

- coordinates:

  Numeric matrix or data frame containing the (x, y) coordinates of
  sampling units.

- min_distance:

  Numeric. Minimum distance between knots in meters. This distance is
  used for both `knotDist` and `minKnotDist` parameters of the
  [Hmsc::constructKnots](https://rdrr.io/pkg/Hmsc/man/constructKnots.html)
  function.

- jitter_distance:

  Numeric. The jitter distance applied to overlapping coordinates to
  avoid exact duplicates. Defaults to 100 meters.

- min_lf, max_lf:

  Integer. Minimum and maximum number of latent factors to be used. Both
  default to `NULL` which means that the number of latent factors will
  be estimated from the data. If either is provided, the respective
  values will be used as arguments to
  [Hmsc::setPriors](https://rdrr.io/pkg/Hmsc/man/setPriors.html).

- alphapw:

  Prior for the alpha parameter. Defaults to a list with `Prior = NULL`,
  `Min = 10`, `Max = 1500`, and `Samples = 101`. If `alphapw` is `NULL`
  or a list with all `NULL` list items, the default prior will be used.
  If `Prior` is a matrix, it will be used as the prior. If
  `Prior = NULL`, the prior will be generated using the minimum and
  maximum values of the alpha parameter (`min` and `max`, respectively;
  in kilometre) and the number of samples (`Samples`). Defaults to a
  prior with 101 samples ranging from 10 to 1500 km, with the first
  value in the second column set to 0.5.

## Value

An object of class `HmscRandomLevel`, suitable for specifying the random
level in HMSC GPP models. This object contains the prepared knot
locations as a data frame with columns `var_1` and `var_2` (numeric
coordinates), after possible jittering and conversion to avoid overlap.

## Author

Ahmed El-Gabbas
