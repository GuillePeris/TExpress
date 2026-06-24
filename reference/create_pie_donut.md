# Create Donut Pie Chart for Single TE Type

Creates a donut/pie chart for a specific TE type showing the breakdown
of expression types. Requires the `webr` package.

## Usage

``` r
create_pie_donut(
  res,
  TE_feature,
  plot.title = NULL,
  output_folder = ".",
  height = 7,
  width = 7,
  device = "png",
  prefix = NULL
)
```

## Arguments

- res:

  Data frame with classification results

- TE_feature:

  Character string. Specific TE to plot (e.g., "LINE", "L1")

- plot.title:

  Character string. Plot title

- output_folder:

  Character string. Output directory

- height:

  Numeric. Plot height in inches. Default is 7.

- width:

  Numeric. Plot width in inches. Default is 7.

- device:

  Character vector. Output format(s). Default is "png".

- prefix:

  Character string. Optional prefix for filename

## Value

Invisible NULL. Called for side effects (creating plot file).

## Details

Please, note that `webr` package is not maintained and this function
throws several warnings of deprecated uses.

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires webr package
if (requireNamespace("webr", quietly = TRUE)) {
  create_pie_donut(
    res = classified_TEs,
    TE_feature = "LINE",
    output_folder = "plots",
    prefix = "detailed"
  )
}
} # }
```
