# Function to scale and/or center the data

Function to scale and/or center the data

## Usage

``` r
center_scale(data, mean_center = TRUE, sd_scale = TRUE)
```

## Arguments

- data:

  matrix or data frame

- mean_center:

  either TRUE or FALSE. If mean_center is TRUE then the mean of each
  column will be subtracted

- sd_scale:

  either TRUE or FALSE. See the details section for more information

## Value

a matrix

## Details

If sd_scale is TRUE and mean_center is TRUE then each column will be
divided by the standard deviation. If sd_scale is TRUE and mean_center
is FALSE then each column will be divided by sqrt( sum(x^2) / (n-1) ).
In case of missing values the function raises an error. In case that the
standard deviation equals zero then the standard deviation will be
replaced with 1.0, so that NaN's can be avoided by division

## Examples

``` r
data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

dat = center_scale(dat, mean_center = TRUE, sd_scale = TRUE)
```
