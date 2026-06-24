# Download Test Data for Package Examples and Testing

Downloads sample data files from Dropbox for testing and running package
examples. The data includes human genome annotations (chromosome 1 only)
and TE count data, suitable for demonstrating differential expression
analysis workflows.

## Usage

``` r
downloadTestData(folder = NULL, timeout = 300)
```

## Arguments

- folder:

  Character string. Directory where test data will be downloaded. If
  NULL (default), uses a temporary directory via
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html). The directory
  will be created if it doesn't exist.

- timeout:

  Numeric. Download timeout in seconds. Default is 300 (5 minutes).
  Increase if downloads fail due to slow connections.

## Value

A list with four elements:

- folder:

  Character string with path to download directory

- metafile:

  Full path to metadata CSV file (data.csv)

- gtf.gene.file:

  Full path to gene GTF file

- gtf.TE.file:

  Full path to TE GTF file

## Details

The function downloads three datasets:

1.  Gene GTF (**~78 MB**): Human gene annotations for chromosome 1 from
    Ensembl GRCh38.115

2.  TE GTF (**~10 MB**): Transposable element annotations for chromosome
    1 from RepeatMasker (hg38)

3.  Count Data (**~2 MB compressed**): Example TEtranscripts count files
    in a ZIP archive. Contains 6 sample count files (3 control, 3
    treatment) and a csv file for TExpress analysis

Files are downloaded from Dropbox using direct download links. The count
data ZIP file is automatically extracted and the archive is removed
after extraction.

## Downloaded Files

The following files will be created in the specified folder:

- `Homo_sapiens.GRCh38.115.chr22.gtf` - Gene annotations

- `hg38_rmsk_TE_chr22.gtf` - TE annotations

- `*.cntTable` - Multiple count files (extracted from ZIP)

- `data.csv` - csv file including info on count files

## File Sizes

Total download size is approximately 90 MB. Ensure adequate disk space
and internet bandwidth. On slow connections, increase the `timeout`
parameter.

## Examples

``` r
if (FALSE) { # \dontrun{
# Download to temporary directory (default)
my.data <- downloadTestData()

# Download to specific directory
my.data <- downloadTestData(folder = "test_data")

# Use downloaded data in analysis
results <- TE_DEA(
  metafile = my.data$metafile,
  folder = my.data$folder,
  gtf.TE.file = my.data$gtf.TE.file
)
} # }
```
