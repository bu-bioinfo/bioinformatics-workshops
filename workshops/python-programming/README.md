
<!-- README.md is generated from README.Rmd. Please edit that file -->

``` r
library(reticulate)
```

# Example

You can use Jupyter but R Markdown is sort of nice too…

``` python
print("Hello world")
```

    Hello world

## List comprehension

``` python
DNA = "ATCGACGCTAGCATCAG"
GC = [x for x in DNA if x == "C" or x == "G"]
print(GC)
```

    ['C', 'G', 'C', 'G', 'C', 'G', 'C', 'C', 'G']

``` python
GC_content = round(len(GC)/len(DNA)*100, 2)
print("GC Content {}%".format(GC_content))
```

    GC Content 52.94%
