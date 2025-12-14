# Affinity propagation clustering

Affinity propagation clustering

## Usage

``` r
AP_affinity_propagation(
  data,
  p,
  maxits = 1000,
  convits = 100,
  dampfact = 0.9,
  details = FALSE,
  nonoise = 0,
  time = FALSE
)
```

## Arguments

- data:

  a matrix. Either a similarity matrix (where number of rows equal to
  number of columns) or a 3-dimensional matrix where the 1st, 2nd and
  3rd column correspond to (i-index, j-index, value) triplet of a
  similarity matrix.

- p:

  a numeric vector of size 1 or size equal to the number of rows of the
  input matrix. See the details section for more information.

- maxits:

  a numeric value specifying the maximum number of iterations (defaults
  to 1000)

- convits:

  a numeric value. If the estimated exemplars stay fixed for convits
  iterations, the affinity propagation algorithm terminates early
  (defaults to 100)

- dampfact:

  a float number specifying the update equation damping level in \[0.5,
  1). Higher values correspond to heavy damping, which may be needed if
  oscillations occur (defaults to 0.9)

- details:

  a boolean specifying if details should be printed in the console

- nonoise:

  a float number. The affinity propagation algorithm adds a small amount
  of noise to *data* to prevent degenerate cases; this disables that.

- time:

  a boolean. If TRUE then the elapsed time will be printed in the
  console.

## Details

The *affinity propagation* algorithm automatically determines the number
of clusters based on the input preference *p*, a real-valued N-vector.
p(i) indicates the preference that data point i be chosen as an
exemplar. Often a good choice is to set all preferences to median(data).
The number of clusters identified can be adjusted by changing this value
accordingly. If *p* is a scalar, assumes all preferences are that shared
value.

The number of clusters eventually emerges by iteratively passing
messages between data points to update two matrices, A and R (Frey and
Dueck 2007). The "responsibility" matrix R has values r(i, k) that
quantify how well suited point k is to serve as the exemplar for point i
relative to other candidate exemplars for point i. The "availability"
matrix A contains values a(i, k) representing how "appropriate" point k
would be as an exemplar for point i, taking into account other points'
preferences for point k as an exemplar. Both matrices R and A are
initialized with all zeros. The AP algorithm then performs updates
iteratively over the two matrices. First, "Responsibilities" r(i, k) are
sent from data points to candidate exemplars to indicate how strongly
each data point favors the candidate exemplar over other candidate
exemplars. "Availabilities" a(i, k) then are sent from candidate
exemplars to data points to indicate the degree to which each candidate
exemplar is available to be a cluster center for the data point. In this
case, the responsibilities and availabilities are messages that provide
evidence about whether each data point should be an exemplar and, if
not, to what exemplar that data point should be assigned. For each
iteration in the message-passing procedure, the sum of r(k; k) + a(k; k)
can be used to identify exemplars. After the messages have converged,
two ways exist to identify exemplars. In the first approach, for data
point i, if r(i, i) + a(i, i) \> 0, then data point i is an exemplar. In
the second approach, for data point i, if r(i, i) + a(i, i) \> r(i, j) +
a(i, j) for all i not equal to j, then data point i is an exemplar. The
entire procedure terminates after it reaches a predefined number of
iterations or if the determined clusters have remained constant for a
certain number of iterations... (
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5650075/ – See chapter 2 )

Excluding the main diagonal of the similarity matrix when calculating
the median as preference ('p') value can be considered as another option
too.

## References

https://www.psi.toronto.edu/index.php?q=affinity

https://www.psi.toronto.edu/affinitypropagation/faq.html

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5650075/ ( SEE chapter 2 )

## Examples

``` r
set.seed(1)
dat = matrix(sample(1:255, 2500, replace = TRUE), 100, 25)

smt = 1.0 - distance_matrix(dat, method = 'euclidean', upper = TRUE, diagonal = TRUE)
diag(smt) = 0.0

ap = AP_affinity_propagation(smt, p = median(as.vector(smt)))

str(ap)
#> List of 10
#>  $ K                  : int 13
#>  $ N                  : num 100
#>  $ netsim             : num -41124
#>  $ dpsim              : num -34361
#>  $ expref             : num -6763
#>  $ iterations         : int 200
#>  $ exemplars          : int [1:13] 0 8 18 26 36 42 64 66 73 75 ...
#>  $ idx                : num [1, 1:100] 0 73 75 42 66 64 64 18 8 75 ...
#>  $ clusters           :List of 13
#>   ..$ 18: int [1:7] 7 18 51 53 68 79 93
#>   ..$ 64: int [1:7] 5 6 37 44 46 64 78
#>   ..$ 92: int [1:6] 11 34 40 77 84 92
#>   ..$ 66: int [1:9] 4 24 31 58 61 66 72 81 83
#>   ..$ 94: int [1:13] 20 28 32 41 48 63 65 69 82 91 ...
#>   ..$ 42: int [1:8] 3 16 35 42 49 60 62 99
#>   ..$ 36: int [1:10] 13 17 19 23 29 30 33 36 85 87
#>   ..$ 88: int [1:6] 10 55 70 80 88 89
#>   ..$ 75: int [1:11] 2 9 47 56 57 71 75 86 90 97 ...
#>   ..$ 8 : int [1:6] 8 14 21 22 45 50
#>   ..$ 73: int [1:5] 1 25 52 59 73
#>   ..$ 26: int [1:6] 26 38 39 54 67 76
#>   ..$ 0 : int [1:6] 0 12 15 27 43 74
#>  $ clusters_vectorized: int [1:100] 0 73 75 42 66 64 64 18 8 75 ...
```
