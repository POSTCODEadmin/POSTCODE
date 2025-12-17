# POSTCODE

**POSTCODE** is a computational framework designed to characterize and interpret
post-transcriptional regulatory mechanisms by integrating transcriptomic,
translational, and protein-related features.

The associated **Shiny web application** allows users to upload gene-level data
(e.g. log2 fold-changes) and explore how intrinsic and extrinsic transcript
features relate to protein expression and post-transcriptional regulation.

🔗 **Live Shiny application**:  
👉 https://postcode-lab.shinyapps.io/POSTCODE/

---

##  Features

- Annotation of gene-level input with post-transcriptional proxies
- Rolling-median and spline-based trend visualization
- Integrated fGSEA analysis across multiple regulatory layers
- Downloadable tables and plots
- Robust error handling and deployment-ready architecture

---

##  Input format

The application expects a **tab-delimited `.txt` file** with at least two columns:

| Column | Description |
|------|-------------|
| Gene name | Gene symbol (HGNC, uppercase recommended) |
| log2FC | Log2 fold-change value |

Example:
```text
TP53    -1.25
MYC      2.10

---

##  More informations

POSTCODE : A tool for deciphering post-transcriptional regulations in omics studies

Several studies have described the lack of correlation between expression at transcriptional and protein levels. This discrepancy can be explained by a multitude of post-transcriptional regulatory mechanisms, such as transcript sequestration, RNA decay, translational efficiency, and protein degradation.
Due to the lack of existing tools, these post-transcriptional mechanisms have remained in the blind spot of omics studies. To bridge this gap, we developed POSTCODE, a bioinformatics pipeline for annotating and reanalyzing omics datasets with PTR-related proxies. By incorporating PTR’s proxies into differential expression analyses, POSTCODE enables the detection of post-transcriptional deregulations, refining functional analysis interpretations and revealing regulatory mechanisms that standard transcriptome or proteome analyses might miss.

1 - Installation
To run the POSTCODE Shiny app, it is essential to have R and RStudio installed. Additionally, the application requires several R packages to function properly. The necessary packages and their specific versions are as follows:  

Core dependencies:  
- 'shiny' (1.8.1.1) and 'shinyjs' (2.1.0) for building the interactive user interface.  
- 'DT' (0.33) for rendering interactive tables.  
- 'dplyr' (1.1.4) and 'tidyverse' (2.0.0) for data manipulation.  
- 'readr' (2.1.5) for reading structured data files.  
- 'DBI' (1.2.3) and 'RSQLite' (2.3.7) for database interactions.  
- 'writexl' (1.5.1) and 'openxlsx' (4.2.8) for exporting data to Excel format.  

Data processing and visualization: 
- 'gridExtra' (2.3) and 'ggplot2' (3.5.1) for advanced plotting.  
- 'tidyr' (1.3.1) for data reshaping.  
- 'plotly' (4.10.4) for interactive visualizations.  
- 'ggrepel' (0.9.6) for enhanced label placement in plots.  
- 'data.table' (1.17.0) for efficient data handling.
- 'pheatmap' (1.0.13) for static heatmap generation.  
- 'stringr' (1.5.1) for string manipulation.  
- 'RColorBrewer' (1.1.3) for generating color palettes.  

Statistical analysis and functional enrichment: 
- 'mgcv' (1.9.0) for generalized additive models.  
- 'fgsea' (1.28.0) for functional gene set enrichment analysis.  

If a package is missing or outdated, it should be updated accordingly. Once all dependencies are installed, the application can be launched within RStudio.

2 - Usage  

The input file must be formatted as a tab-separated (.txt) file with two columns:  
First column named "Genename" – list of gene names.  
Second column named depending on the experiment – corresponding fold change (FC) values.  
Ensure that the file itslef is named 'Input.txt' and follows this structure for proper functioning of the application. A sample input file is provided in the repository for testing.  

Application workflow  

Once the input file is uploaded, the app provides several interactive panels:  
- Annotations Panel: Displays genes annotations, their associated FC values, and various PTR-related proxies. Users can search for specific genes using the search bar.  
- Plots Panel: Visualizes data through three types of graphs.  
- fGSEA Panel: Allows users to perform fGSEA analysis by selecting a dataset and running the analysis. The "Download Volcano Plot Data" button enables users to download a table to recreate the volcano plot.  
- About Section: Provides descriptions of different PTR proxies used in the analysis.  
- Each time postcode is run with a dataset, a new time stamped folder is automatically created inside the "postcode_results" folder. This folder contains the results of statistical tests that identify which proxies are significantly associated with gene expression changes, visualized through heatmaps.
Several download buttons are available across different tabs, allowing users to save tables and graphs for further use.  

Requirements  

To ensure proper functionality, all files from the GitHub repository must be downloaded and placed in the same directory. 
If the input file is located elsewhere, users must specify its correct path in the code before proceeding.

3 - Contact & further inquiries  
For any questions or assistance, please contact:  
ismael.boussaid@aphp.fr  

Link to the study: [we should add the link here later]
