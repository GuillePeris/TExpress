# Save ggplot to Multiple File Formats

Internal function to save a ggplot object to one or more file formats.
Supports common vector (SVG, EPS) and raster (PNG, TIFF, JPEG) formats.

## Usage

``` r
printDevice(plot, basename, device, width = 7, height = 7)
```

## Arguments

- plot:

  A ggplot object to save.

- basename:

  Character string. Base name for output files (without extension).
  Extensions will be added automatically based on requested formats.

- device:

  Character vector. File format(s) to save. Supported formats: "svg",
  "eps", "png", "tiff", "jpeg". Case-insensitive. Dots are automatically
  removed (e.g., ".png" becomes "png").

- width:

  Numeric. Plot width in inches. Default is 7.

- height:

  Numeric. Plot height in inches. Default is 7.
