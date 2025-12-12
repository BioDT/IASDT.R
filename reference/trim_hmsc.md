# Trim an Hmsc Model Object by Removing Specified Components

Removes specified components from an `Hmsc` model object one at a time.
This is used to keep only a smaller `Hmsc` object, containing only
needed list items. This is useful in situations in which the fitted Hmsc
model is large (e.g. in GBs) and downstream computations are implemented
on parallel.

## Usage

``` r
trim_hmsc(model, names_to_remove = NULL)
```

## Arguments

- model:

  An object of class `Hmsc`, containing fitted Hmsc model. Must not be
  `NULL`.

- names_to_remove:

  A character vector specifying the names of components to remove from
  the model (e.g., `"postList"`, `"Y"`). If `NULL`, no trimming is
  implemented. Must be non-empty and match names in the model.

## Value

An `Hmsc` object with the specified components removed.

## Details

This function is used to reduce the memory footprint of `Hmsc` models by
removing unnecessary components before further processing. It converts
the model to a plain list to avoid S3 method overhead, removes
components iteratively, and restores the `Hmsc` class. The iterative
approach is slower than vectorized subsetting but may be preferred for
specific use cases requiring step-by-step removal. The simple trimming
of list items; e.g. `model[c("postList", X")] <- NULL` took comparably
too much time for trimming done inside functions for jobs submitted via
SLURM

## Author

Ahmed El-Gabbas

## Examples

``` r
require(Hmsc)

(model <- Hmsc::TD$m)
#> Hmsc object with 50 sampling units, 4 species, 3 covariates, 3 traits and 2 random levels
#> Posterior MCMC sampling with 2 chains each with 100 samples, thin 1 and transient 50 

(trimmed_model <- trim_hmsc(
  model, names_to_remove = c("postList", "rL", "ranLevels")))
#> Hmsc object with 50 sampling units, 4 species, 3 covariates, 3 traits and 2 random levels
#> Posterior MCMC sampling with 100 samples, thin 1 and transient 50 

setdiff(names(model), names(trimmed_model))
#> [1] "ranLevels" "rL"        "postList" 

lobstr::obj_size(model)
#> 937.01 kB
lobstr::obj_size(trimmed_model)
#> 44.06 kB
```
