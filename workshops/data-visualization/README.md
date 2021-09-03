
<!-- README.md is generated from README.Rmd. Please edit that file -->

``` r
library(ggplot2)
library(scales)
library(ggpubr)
library(janitor)
library(dplyr)
```

# Data melting

When using `ggplot` your data needs to be in canonical or melted format.
The format has only one independently observed value(s) in each row. All
of the other columns are metadata that explain what was measured.

``` r
# Bad
df <- data.frame(gene=LETTERS[1:5], group_a=rnorm(5), group_b=rnorm(5))
df
```

``` 
  gene    group_a     group_b
1    A  1.0238116  1.16843270
2    B -1.5122120 -1.57720685
3    C  0.7889155 -0.39208111
4    D  0.5386369  0.15199257
5    E -0.6234121 -0.06627217
```

``` r
# Good
df.melt <- reshape2::melt(df, value.name="expression")
df.melt
```

``` 
   gene variable  expression
1     A  group_a  1.02381156
2     B  group_a -1.51221196
3     C  group_a  0.78891552
4     D  group_a  0.53863690
5     E  group_a -0.62341213
6     A  group_b  1.16843270
7     B  group_b -1.57720685
8     C  group_b -0.39208111
9     D  group_b  0.15199257
10    E  group_b -0.06627217
```

``` r
df.melt %>%
  ggplot(aes(x=gene, y=expression, group=variable, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge())
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

A useful trick for when your data is in a nested list. If you need to
compute things in a for loop, it’s better to use an apply function first
then convert to a dataframe all at once (rather than trying to insert a
new row into an existing dataframe).

``` r
# Bad
data <- mapply(function(g) {
  data.frame(gene=g, group_a=rnorm(1), group_b=rnorm(1))
}, LETTERS[1:5], SIMPLIFY=FALSE, USE.NAMES=FALSE)
data
```

    [[1]]
      gene  group_a  group_b
    1    A 2.538771 1.218381
    
    [[2]]
      gene   group_a    group_b
    1    B -1.205954 -0.1796525
    
    [[3]]
      gene    group_a    group_b
    1    C -0.7756864 -0.6308236
    
    [[4]]
      gene  group_a    group_b
    1    D 1.086061 -0.9887432
    
    [[5]]
      gene    group_a     group_b
    1    E -0.1319474 -0.04987285

``` r
# Fine
df <- do.call(rbind, data)
df
```

``` 
  gene    group_a     group_b
1    A  2.5387708  1.21838068
2    B -1.2059539 -0.17965252
3    C -0.7756864 -0.63082359
4    D  1.0860610 -0.98874325
5    E -0.1319474 -0.04987285
```

``` r
# Good
df.melt <- reshape2::melt(df, value.name="expression")
df.melt
```

``` 
   gene variable  expression
1     A  group_a  2.53877079
2     B  group_a -1.20595394
3     C  group_a -0.77568638
4     D  group_a  1.08606097
5     E  group_a -0.13194741
6     A  group_b  1.21838068
7     B  group_b -0.17965252
8     C  group_b -0.63082359
9     D  group_b -0.98874325
10    E  group_b -0.04987285
```

# Style and Rendering

## Fonts

*Some of these might require a license or are difficult to get for free*

    - Avenir
    - Bell Gothic
    - Courier
    - Frankin Gothic
    - Frutiger
    - Futura
    - Garamond
    - Helvetica
    - Futura

R does not necessarily have access to your fonts library

``` r
library(extrafont)
font_import() # Searches for system-level fonts
fonts() # Here's what it found
loadfonts() # Register the fonts with R
```

## Style

``` r
data(iris)
df <- janitor::clean_names(iris)
head(df)
```

``` 
  sepal_length sepal_width petal_length petal_width species
1          5.1         3.5          1.4         0.2  setosa
2          4.9         3.0          1.4         0.2  setosa
3          4.7         3.2          1.3         0.2  setosa
4          4.6         3.1          1.5         0.2  setosa
5          5.0         3.6          1.4         0.2  setosa
6          5.4         3.9          1.7         0.4  setosa
```

``` r
df %>%
  ggplot(aes(x=sepal_length, y=petal_length, color=species)) +
  geom_point()
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

The key ggplot parameters you will constantly modify. There are [many
more](https://ggplot2.tidyverse.org/reference/theme.html) if you’re
interested.

``` r
ggstyle <- function(font="Helvetica", scale=1) {
  fs <- function(x) x*scale # Dynamic font scaling
  ggplot2::theme(
    plot.title = ggplot2::element_text(family=font, size=fs(26), face="bold", color="#222222"),
    plot.subtitle = ggplot2::element_text(family=font, size=fs(18), margin=ggplot2::margin(0,0,5,0)),
    plot.caption = ggplot2::element_blank(),
    legend.position = "right",
    legend.text.align = 0,
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(family=font, size=fs(18), color="#222222"),
    axis.title =  ggplot2::element_text(family=font, size=fs(18), color="#222222"),
    axis.text = ggplot2::element_text(family=font, size=fs(18), color="#222222"),
    axis.text.x = ggplot2::element_text(margin=ggplot2::margin(5, b=10)),
    axis.line = ggplot2::element_line(color="#222222"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill="white"),
    strip.text = ggplot2::element_text(size=fs(22), hjust=0)
  )
}
```

``` r
df %>%
  ggplot(aes(x=sepal_length, y=petal_length, color=species)) +
  labs(title="Iris Data", subtitle="Sepal Length vs. Petal Length", x="Sepal Length", y="Petal Length") +
  geom_point() +
  ggstyle()
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
# You can use a scaling factor to increase/decreasing font size while maintaining ratios
# this is useful when you need to boost legibility of a small figure for grant proposals with limited space
df %>%
  ggplot(aes(x=sepal_length, y=petal_length, color=species)) +
  labs(title="Iris Data", subtitle="Sepal Length vs. Petal Length", x="Sepal Length", y="Petal Length") +
  geom_point() +
  ggstyle(font="Frutiger", scale=0.75)
```

![](README_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

## Color

Use the HCL color space. The HCL space is particularly useful for
specifying individual colors and color palettes, as its three axes match
those of the human visual system very well.

``` r
library(colorspace)
hcl_palettes("Sequential (multi-hue)", n=8, plot=TRUE)
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
df %>%
  ggplot(aes(x=species, y=sepal_width, color=sepal_length)) +
  labs(title="Iris Data", subtitle="Sepal Dimensions Across Species", x="Species", y="Sepal Width" , color="Sepal Length") +
  geom_jitter(width=0.3, size=3) +
  ggstyle(font="Frutiger", scale=0.75) +
  colorspace::scale_color_continuous_sequential(palette="Plasma", begin=0.1, end=0.9) +
  # Overwrite the ggtheme to change the legend title and position
  theme(legend.position="top", 
        legend.title=element_text(family="Frutiger", size=18, color="#222222"))
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
hcl_palettes("Qualitative", n=8, plot = TRUE)
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

## Publication Quality Rendering

When you’re making figures, save a high quality raw TIFF that can be
imported into Photoshop for alignment with other figures. If you want to
guarantee that a journal will accept your image quality, make the raw
images close to the size of A4 standard sheet (8.3 × 11.7 inches) with a
DPI \>= 300 and then resize them later. Some journals require 600 DPI
for certain figure types so always check with the requirements first.

``` r
ggplot2::ggsave(
    filename="Figure_01.tiff",
    plot=p,
    path="/Figures",
    scale=1,
    width=8,
    height=10,
    units="in",
    dpi=300
)
```
