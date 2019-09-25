# SCAP
An Rshiny portal for analysis of scRNA-seq data

## Setup

```R
install.packages("devtools")
devtools::install_github("JoelPHoward/SCAP")
```

## Running the app
```R
SCAP::SCAP()
```
## Uploading your data
Currently SCAP can only process Seurat objects.
For better performance, especially with large datasets, the `.rds` files are first converted to `.loom`.

* Navitgate to the `File Conversion` tab.
* Select the Seurat object you wish to analyze and the folder where you want to save the `.loom` file(s) (note, initially this folder should be empty, such that only the `.loom` files will be present in the folder).
* Click `Convert`.

When the conversion is complete, navigate back to the main tab and upload your data.