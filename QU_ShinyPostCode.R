# Load required libraries for building and running the Shiny app
library(shiny)
library(shinyjs)
library(DT)
library(dplyr)
library(readr)
library(DBI)
library(RSQLite)
library(writexl)
library(gridExtra)
library(tidyr)
library(ggplot2)
library(mgcv)
library(fgsea)
library(data.table)
library(ggrepel)
library(plotly)
library(tidyverse)
library(openxlsx)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(viridis)

# Ensures the GMT file ends with a blank line by appending one if needed
fix_gmt_file <- function(gmt_file) {
  lines <- readLines(gmt_file)
  if (lines[length(lines)] != "") {
    writeLines(c(lines, ""), gmt_file)
  }
}

# Function to process a GMT file for fGSEA analysis and save the results
process_gmt_file <- function(gmt_file, ranks, output_file) {
  fix_gmt_file(gmt_file)  # Ensure proper GMT format
  pathways <- gmtPathways(gmt_file)
  fgsea_results <- fgseaMultilevel(
    pathways = pathways, stats = ranks,
    nPermSimple = 20000, minSize = 15, maxSize = 5000
  )
  # Convert to data frame and ensure 'pathway' column exists
  fgsea_results <- as.data.frame(fgsea_results)
  if (!"pathway" %in% colnames(fgsea_results)) {
    fgsea_results$pathway <- rownames(fgsea_results)
  }
  fwrite(fgsea_results, file = output_file, sep = "\t", sep2 = c("", " ", ""))
  return(fgsea_results)
}

# Lists of GMT files and corresponding output file names for FGSEA analysis
gmt_files <- list(
  motif3UTR = "motif3UTRfgsea.gmt",
  motif5UTR = "motif5UTRfgsea.gmt",
  mirna = "mirnafgsea.gmt",
  locPROT = "locPROTfgsea.gmt",
  locRNA = "locRNAfgsea.gmt",
  protinfo = "protinfofgsea.gmt",
  RBP = "RBPfgsea.gmt",
  proteins = "proteinsfgsea.gmt",
  motifRNAepigenetic = "motifRNAepigenetfgsea.gmt"
)

output_files <- list(
  motif3UTR = "fgseamotif3UTR_multilevel.txt",
  motif5UTR = "fgseamotif5UTR_multilevel.txt",
  mirna = "fgseamirna_multilevel.txt",
  locPROT = "fgsealocPROT_multilevel.txt",
  locRNA = "fgsealocRNA_multilevel.txt",
  protinfo = "fgseaprotinfo_multilevel.txt",
  RBP = "fgseaRBP_multilevel.txt",
  proteins = "fgseaproteins_multilevel.txt",
  motifRNAepigenetic = "fgseamotifRNAepigenetic_multilevel.txt"
)

# Define the user interface (UI) for the Shiny app
ui <- fluidPage(
  titlePanel("POSTCODE"),
  useShinyjs(),  # To enable shinyjs for modal visibility
  
  # Add styles for the modal dialog
  tags$head(
    tags$style(HTML("
      .shiny-notification {
        position: fixed;
        top: 93%;
        left: 50%;
        transform: translate(-50%, -50%);
        font-size: 16px;
        font-weight: bold;
        color: black;
        background-color: rgba(211, 211, 211, 0.9);
        border-radius: 8px;
        padding: 20px;
        box-shadow: 0px 4px 6px rgba(0, 0, 0, 0.2);
        z-index: 1050;
      }
      .modal {
        display: none; 
        background-color: rgba(0, 0, 0, 0.5); 
        position: fixed; 
        top: 0; 
        left: 0; 
        width: 100%; 
        height: 100%; 
        z-index: 1050;
      }
      .modal-dialog {
        margin: auto; 
        display: flex; 
        align-items: center; 
        justify-content: center; 
        height: 100%;
      }
      .modal-content {
        width: 50%; 
        padding: 20px; 
      }
      .modal-header {
        justify-content: center; 
        position: relative;
      }
      .modal-body {
        text-align: center; 
        font-size: 16px; 
        font-weight: bold; 
        margin-top: 20px; 
        color: red;
      }
    "))
  ),
  
  # Modal dialog displayed during fgsea analysis
  div(
    id = "loading-modal",  
    class = "modal",  
    tabindex = "-1",  
    role = "dialog",
    div(
      class = "modal-dialog modal-dialog-centered",  
      role = "document",
      div(
        class = "modal-content",  
        div(
          class = "modal-header",  
          h5("fGSEA analysis in progress", style = "font-size: 18px; font-weight: bold; text-align: center;")
        ),
        div(
          class = "modal-body",  
          tags$p("Please wait !")
        )
      )
    )
  ),
  # Layout with sidebar for controls and main panel for outputs
  sidebarLayout(
    sidebarPanel(
      # File input for uploading a TXT file
      fileInput("file1", "Choose a TXT file", accept = c(".txt")),
      # Button to submit the uploaded file
      actionButton("enter", "Submit"),
      # Download buttons for data and plot exports
      downloadButton("downloadTable", "Download Data"),
      downloadButton("downloadPlots", "Download Graphs")
    ),
    
    mainPanel(
      # Tabs for organizing app content
      tabsetPanel(
        # Tab for viewing data table annotations
        tabPanel("Annotations", DTOutput("table")),
        # Tab for displaying plots with selectable options
        tabPanel("Plots",
                 selectInput("curve_selection", "Choose curves to display :",
                             choices = list(
                               "Raw Data" = "raw",
                               "Raw Data and Median" = "both",
                               "Spline" = "spline"
                             ),
                             
                 ),
                 plotOutput("plots", width = "1100px", height = "2500px")),
        # Tab for fGSEA analysis options and volcano plot display
        tabPanel("fGSEA", 
                 selectInput("dataset", "Select Dataset", choices = names(gmt_files)),
                 actionButton("run_fgsea", "Run fGSEA analysis"), 
                 downloadButton("downloadPlot", "Download Volcano Plot"),
                 downloadButton("downloadData", "Download Volcano Plot Data"),
                 plotlyOutput("volcanoPlot")
        ),
        
        # Tab for project information and metric descriptions
        tabPanel("About",
                 HTML("
            <h3>About</h3>
            <p>Our tool focuses directly on intrinsic and extrinsic characteristics of transcripts that regulate their post-transcriptional fate leading to greater or lesser protein expression. We have shown that all these characteristics form a grammar and punctuation in the language of gene expression in the cell. Decoding the latter enables us to measure a set of dynamic post-transcriptional regulations that can occur between two conditions and completely remodel gene expression. </p>
            <p>If you use this tool, please cite (Khadra et al 20XX).</p>
            <p>For details on how to determine proxies, please refer to the materials and methods section of the article.</p>
            
            <b>Metrics:</b>
            <ul>
              <li><b>PTR_AML</b> (PTR = Protein to mRNA ratio)
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> The ratio of protein levels (quantified by LFQ, Label-Free Quantification) to mRNA levels (quantified by TPM). Our findings show that PTRs are conserved between AML samples and human tissues, underscoring the universality of underlying mechanisms.</li>
                </ul>
              </li>
              <br>
              <li><b>UTR5_RNA_folding_per_base</b>: 5' energy, enthalpy per base of secondary structures
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Measures the energy of secondary structures in the 5' Untranslated Region (UTR) per base. Secondary structures can influence translation by affecting ribosome access to mRNA. Higher energy may indicate more stable and potentially repressive structures.</li>
                </ul>
              </li>
              <br>
              <li><b>UTR5_RNA_length</b>: Average 5' length
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Average length of the 5' UTR among transcript isoforms.</li>
                </ul>
              </li>
              <br>
              <li><b>UTR5_RNA_folding</b>: Enthalpy of secondary structures in the 5' UTR
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Total energies associated with secondary structures in the 5' UTR among transcript isoforms.</li>
                </ul>
              </li>
              <br>
              <li><b>UTR5_GC_percent</b>: Percentage of G and C bases in the 5' UTR
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Percentage of guanine (G) and cytosine (C) bases in the 5' UTR. GC-rich regions can form more stable secondary structures, influencing translation.</li>
                </ul>
              </li>
              <br>
              <li><b>Mean_MRL</b>: Mean Ribosome Load
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Rate of translation initiation influenced by 5’UTR sequence  We used data from a massively parallel translation assay that determined a proxy for initiation rate called “Mean Ribosome Load” (MRL). This was assessed by inserting a 50-nt preSTART 5'UTR sequence into plasmid constructs for 20,530 human transcripts corresponding to 12,503 genes. </li>
                  <li>If you use this metric, also cite (Sample et al 2019). </li>
                  <li>MRL indicates the average ribosome load on a given mRNA, reflecting the efficiency of translation initiation driven by the 5' UTR. </li>
                </ul>
              </li>
              <br>
              <li><b>UTR3_RNA_folding_per_base</b>: 3' energy, enthalpy per base of secondary structures
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Measures the energy of secondary structures in the 3' UTR per base .</li>
                </ul>
              </li>
              <br>
              <li><b>UTR3_RNA_length</b>: Average length in 3'
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Average length of the 3' UTR downstream of the stop codon.</li>
                  <li> The length of the 3' UTR can affect mRNA stability and the likelihood of interactions with regulatory proteins (RBP = RNA Binding Proteins) or repressive microRNAs.</li>
                </ul>
              </li>
              <br>
              <li><b>UTR3_RNA_folding</b>: Enthalpy of secondary structures in the 3' UTR
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Total energies associated with secondary structures in the 5' UTR. Stable secondary structures in the 3' UTR could influence interactions with regulatory proteins (RBP = RNA Binding Proteins) or repressive microRNAs.</li>
                </ul>
              </li>
              <br>
              <li><b>UTR3_GC_percent</b>: Percentage of G and C bases in the 3' UTR
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Percentage of guanine (G) and cytosine (C) bases in the 3' UTR.</li>
                </ul>
              </li>
              <br>
              <li><b>Poly_A_length</b>: PolyA Tail
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Length of the PolyA tail added to the mRNA.</li>
                  <li>If you use this metric, also cite (Chang et al 2014).</li>
                </ul>
              </li>
              <br>
              <li><b>CAI</b>: Codon Adaptation Index
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> The Codon Adaptation Index (CAI) is the most widespread technique for analyzing codon usage bias. A high CAI indicates a high use of the preferential codon usage, which could enhance translation efficiency or transcript stability. </li>
                  <li>If you use this metric, also cite (Puigbo et al 2008).</li>
                </ul>
              </li>
              <br>
              <li><b>CDS_GC_percent</b>: Percentage of G and C bases in the CDS
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Percentage of G and C bases in the coding region (CDS).</li>
                </ul>
              </li>
              <br>
              <li><b>EJC_density</b>: EJC per nucleotide
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Density of exon-junction complexes (EJC) per nucleotide.</li>
                </ul>
              </li>
              <br>
              <li><b>TR</b>: Translation Rate (HeLa)
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Translation rate measured in HeLa cells.</li>
                  <li>If you use this metric, also cite (Lian et al 2016).</li>
                </ul>
              </li>
              <br>
              <li><b>EVI</b>: Elongation Velocity Index (HeLa)
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Index of ribosome elongation speed on mRNA in HeLa cells.</li>
                  <li>If you use this metric, also cite (Lian et al 2016).</li>
                </ul>
              </li>
              <br>
              <li><b>Ribosome_density</b>: Density (TR/EVI)
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Ratio of translation rate (TR) to elongation velocity index (EVI).</li>
                  <li>If you use this metric, also cite (Lian et al 2016).</li>
                </ul>
              </li>
              <br>
              <li><b>RNA_Halflife</b>: RNA Half-life
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Lifetime of the mRNA before decay.</li>
                  <li> If you use this metric, also cite (Agarwal et al 2022).</li>
                </ul>
              </li>
              <br>
              <li><b>MirDB_sites_number</b>: Number of miRNA sites per gene
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Number of sites potentially targeted by microRNAs (miRNA).</li>
                  <li>If you use this metric, also cite (Chen et al 2020).</li>
                </ul>
              </li>
              <br>
              <li><b>Protein_halflife</b>: Protein Half-life
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Lifetime of the protein before decay.</li>
                  <li>If you use this metric, also cite (Zecha et al 2018).</li>
                </ul>
              </li>
              <br>
              <li><b>Translation_Efficiency</b>: Log2(Translation Efficiency)
                <ul>
                  <li><FONT color= 'blue'>Description:</FONT> Logarithm (base 2) of translation efficiency.</li>
                  <li>If you use this metric, also cite (Khajuria et al 2018).</li>
                </ul>
              </li>
              <br>
              <li><b>Protein_length</b>
                <ul>
                  <li><FONT color='blue'>Description:</FONT> The number of amino acids in the protein sequence.</li>
                  <li> The protein length affects its structure and function; longer proteins often have more complex structures and roles.</li>
                </ul>
              </li>
              <br>
              <li><b>mW</b> (Molecular Weight)
                <ul>
                  <li><FONT color='blue'>Description:</FONT> The molecular weight of the protein, measured in Daltons.</li>
                  <li>It affects the protein's behavior in biological systems, such as how it folds and interacts with other molecules.</li>
                </ul>
              </li>
              <br>
              <li><b>pI</b> (Isoelectric Point)
                <ul>
                  <li><FONT color='blue'>Description:</FONT> The pH at which the protein has no net electrical charge.</li>
                  <li><It determines protein solubility and behavior in different pH environments, impacting its stability and interaction with other molecules.</li>
                </ul>
              </li>
              <br>
              <li><b>Instability_index</b>
                <ul>
                  <li><FONT color='blue'>Description:</FONT> Predicts the stability of the protein in a test tube environment.</li>
                  <li>It helps to determine if the protein is likely to be stable under experimental conditions, impacting its usability in studies.</li>
                </ul>
              </li>
              <br>
              <li><b>Hydropathy</b>
                <ul>
                  <li><FONT color='blue'>Description:</FONT> Indicates the protein’s hydrophobic or hydrophilic nature .</li>
                  <li>It affects how the protein folds and interacts with other molecules, influencing its function and localization within the cell.</li>
                </ul>
              </li>
              <br>
              <li><b>Aromaticity</b>
                <ul>
                  <li><FONT color='blue'>Description:</FONT> Percentage of aromatic amino acids (e.g., phenylalanine, tyrosine, tryptophan) in the protein sequence.</li>
                  <li>Aromatic residues contribute to protein stability and are often involved in molecular interactions.</li>
                </ul>
              </li>
              <br>
              <li><b>Amino Acid Composition</b>
                <ul>
                  <li><FONT color='blue'>Description:</FONT> The proportion of each amino acid type in the protein sequence .</li>
                  <li>It determines the protein structure and function, influencing stability, folding, and interaction with other molecules.</li>
                </ul>
              </li>
            </ul>
      	")
        )
      )
    )
  )
)

# Function to calculate a custom rolling median for a numeric vector
calculate_custom_rolling_median <- function(x) {
  result <- numeric(length(x))
  # Loop through each element in the vector
  for (i in 1:length(x)) {
    if (i <= 101) {
      result[i] <- median(x[1:(100 + i)], na.rm = TRUE)
    } else if (i > (length(x) - 100)) {
      result[i] <- median(x[(i - 100):length(x)], na.rm = TRUE)
    } else {
      result[i] <- median(x[(i - 100):(i + 100)], na.rm = TRUE)
    }
  }
  return(result)
}

# Add statistical tests
run_proxy_comparison <- function(data, log2fc_colname = "log2FC", genename_col = "Genename") {
  
  # Define colors for palettes
  palette_combined <- colorRampPalette(c("blue", "white", "red"))(100)
  palette_direction <- colorRampPalette(c("blue", "white", "red"))(3) 
  palette_fdr <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # Create timestamped output directory
  timestamp <- format(Sys.time(), "%d_%m_%Y_%H%M%S")
  outdir <- file.path("postcode_results", timestamp)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Identify numeric proxy columns (excluding _median)
  proxy_cols <- setdiff(colnames(data), c(genename_col, log2fc_colname))
  proxy_cols <- proxy_cols[sapply(data[proxy_cols], is.numeric)]
  proxy_cols <- proxy_cols[!grepl("_median$", proxy_cols)]
  
  # Filter over and under-expressed genes
  over_expr <- data %>% filter(!!sym(log2fc_colname) > 0.58)
  under_expr <- data %>% filter(!!sym(log2fc_colname) < -0.58)
  
  # Mann–Whitney tests 
  results <- lapply(proxy_cols, function(proxy) {
    # Extract proxy values for over- and under-expressed gene sets
    a <- over_expr[[proxy]]
    b <- under_expr[[proxy]]
    
    # Count non-missing values in each group
    n_over <- sum(!is.na(a))
    n_under <- sum(!is.na(b))
    
    # Proceed only if both groups have at least 3 values
    if (n_over >= 3 && n_under >= 3) {
      # Try running the Mann–Whitney U test 
      pval <- tryCatch({
        wilcox.test(a, b, exact = FALSE)$p.value
      }, error = function(e) NA)
      
      # Compute direction: +1 if median is higher in over group, -1 otherwise
      direction <- ifelse(median(a, na.rm = TRUE) > median(b, na.rm = TRUE), 1, -1)
      
      # Compute group means and fold-change
      mean_over <- mean(a, na.rm = TRUE)
      mean_under <- mean(b, na.rm = TRUE)
      fc_ratio <- ifelse(mean_under == 0, NA, mean_over / mean_under)
      
      # Return one-row tibble with test results
      return(tibble(
        Proxy = proxy,
        P_value = pval,
        Direction = direction,
        N_over = n_over,
        N_under = n_under,
        Mean_Over = mean_over,
        Mean_Under = mean_under,
        FoldChange_Over_Under = fc_ratio
      ))
      
    } else {
      # If not enough data, log a message and exclude this proxy
      message(paste0("Proxy '", proxy, "' excluded: not enough data (N_over = ", n_over, ", N_under = ", n_under, ")"))
      return(NULL)
    }
  }) %>% bind_rows()
  
  # Adjust p-values (FDR)
  results$FDR <- p.adjust(results$P_value, method = "BH")
  
  # Define filtering thresholds and subfolders
  thresholds <- list(
    FDR_0.05 = function(df) df$FDR < 0.05,
    FDR_0.01 = function(df) df$FDR < 0.01,
    FDR_ALL = function(df) rep(TRUE, nrow(df))
  )
  
  # Prevents cut() error caused by non-unique breaks when matrix has a single constant value
  safe_breaks <- function(mat, n = 100) {
    val <- unique(na.omit(as.vector(mat)))
    if (length(val) == 1) {
      delta <- ifelse(val == 0, 0.01, abs(val) * 0.001)
      return(seq(val - delta, val + delta, length.out = n + 1))
    } else {
      return(seq(min(val), max(val), length.out = n + 1))
    }
  }
  
  # Loop over thresholds to create subfolders and plots
  for (folder in names(thresholds)) {
    subdir <- file.path(outdir, folder)
    dir.create(subdir, showWarnings = FALSE)
    
    filtered <- results[thresholds[[folder]](results), ]
    
    if (nrow(filtered) == 0) {
      message(paste("No significant values for", folder, "- skipping heatmap generation."))
      if (dir.exists(subdir)) {
        unlink(subdir, recursive = TRUE)
      }
      next
    }
    
    # Save filtered Excel file
    write.xlsx(filtered, file.path(subdir, "proxy_comparisons.xlsx"))
    
    # Create heatmap matrices
    df_fdr <- filtered %>% select(Proxy, FDR) %>% column_to_rownames("Proxy")
    df_dir <- filtered %>% select(Proxy, Direction) %>% column_to_rownames("Proxy")
    df_combined <- filtered %>% mutate(Combined = -log10(FDR) * Direction) %>%
      select(Proxy, Combined) %>% column_to_rownames("Proxy")
    
    # Format for display
    formatted_fdr <- format(df_fdr, digits = 5, scientific = TRUE)
    n_proxies <- nrow(df_fdr)
    height_px <- max(1400, 40 * n_proxies)
    
    # Heatmap 1: FDR
    df_fdr <- as.matrix(df_fdr)                     
    png(file.path(subdir, "heatmap_fdr.png"), width = 1600, height = height_px, res = 150)
    pheatmap(df_fdr,
             main = paste("Heatmap of adjusted p-values", folder),
             cluster_rows = FALSE, cluster_cols = FALSE,
             display_numbers = formatted_fdr,
             color = palette_fdr,
             breaks = safe_breaks(df_fdr))         
    dev.off()
    
    # Heatmap 2: Direction
    df_dir <- as.matrix(df_dir)
    png(file.path(subdir, "heatmap_direction.png"), width = 1600, height = height_px, res = 150)
    pheatmap(df_dir,
             main = paste("Direction of proxy variation", folder),
             cluster_rows = FALSE, cluster_cols = FALSE,
             display_numbers = TRUE,
             color = palette_direction,
             breaks = c(-1.5, -0.5, 0.5, 1.5))
    dev.off()
    
    # Heatmap 3: Combined score
    df_combined <- as.matrix(df_combined)         
    png(file.path(subdir, "heatmap_combined.png"), width = 1600, height = height_px, res = 150)
    pheatmap(df_combined,
             main = paste("Combined score: direction × -log10(FDR)", folder),
             cluster_rows = FALSE, cluster_cols = FALSE,
             display_numbers = TRUE,
             color = palette_combined,
             breaks = safe_breaks(df_combined))   
    dev.off()
  }
}

# Server function for handling data processing and interactions in the Shiny app
server <- function(input, output, session) {
  # Reactive storage for FGSEA results
  fgsea_results_list <- reactiveValues()
  
  # Observe the 'run_fgsea' button click
  observeEvent(input$run_fgsea, {
    req(input$dataset) # Ensure a dataset is selected
    
    # Show the loading modal
    shinyjs::show("loading-modal")  
    
    # Function to create .rnk file from uploaded input with automatic header detection
    create_rnk_file <- function(input_file = "Input.txt", output_file = "ranked_data.rnk") {
      # Try reading with header = TRUE
      data_try <- tryCatch(
        read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      
      if (!is.null(data_try) && ncol(data_try) >= 2) {
        # Assume first two columns are ID and score
        colnames(data_try)[1:2] <- c("ID", "t")
        data <- data_try[, 1:2]
      } else {
        # Fallback: read without header
        data <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        colnames(data)[1:2] <- c("ID", "t")
        data <- data[, 1:2]
      }
      
      # Force numeric conversion and remove invalid rows
      data$t <- as.numeric(data$t)
      data <- data[!is.na(data$t), ]
      
      # Sort descending
      data <- data[order(-data$t), ]
      
      # Save to .rnk with header
      write.table(data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    
    # Run the function to create the .rnk file
    create_rnk_file()
    
    # Load the ranks
    ranks <- read.table("ranked_data.rnk",
                        header = TRUE, colClasses = c("character", "numeric"))
    ranks <- setNames(ranks$t, ranks$ID)
    
    duplicates  <- names(ranks)[duplicated(names(ranks))]
    print(duplicates)
    
    ranks <- ranks[!duplicated(names(ranks))]
    
    # Process selected dataset 
    fgsea_results <- process_gmt_file(gmt_files[[input$dataset]], ranks, output_files[[input$dataset]])
    fgsea_results_list[[input$dataset]] <- fgsea_results
    stopifnot("pathway" %in% colnames(fgsea_results)) # Ensure 'pathway' is present
    
    # Hide the loading modal after completion
    shinyjs::hide("loading-modal")
    
    # Show a success notification
    showNotification("fGSEA analysis completed successfully!", type = "message")
  })
  
  final_data <- reactiveVal(data.frame())
  # Observe the 'enter' button event to process the user-uploaded file
  observeEvent(input$enter, {
    # Ensure a file is uploaded
    req(input$file1)
    # Read and preprocess the uploaded genes data
    genes_data <- read.table(input$file1$datapath, sep = "\t", header = FALSE, 
                             col.names = c("Genename", "log2FC")) %>%
      mutate(
        log2FC = suppressWarnings(as.numeric(as.character(log2FC))),     # Convert log2FC to numeric
        Genename = toupper(Genename)                   # Normalize gene names to uppercase
      ) %>%
      filter(!is.na(log2FC)) %>%                       # remove genes with missing FC values
      arrange(log2FC)                                  # sort by log2FC ascending before any calculation
    
    # Connect to the SQLite database and retrieve necessary tables
    con <- dbConnect(SQLite(), "Postcode.sqlite")
    CDS <- dbReadTable(con, "CDS")
    Gene <- dbReadTable(con, "Gene")
    Protein <- dbReadTable(con, "Protein")
    UTR3 <- dbReadTable(con, "UTR3")
    UTR5 <- dbReadTable(con, "UTR5")
    # Disconnect from the database
    dbDisconnect(con)
    
    # Merge genes data with database tables to include additional gene metrics
    merged_data <- genes_data %>%
      left_join(Gene, by = c("Genename" = "Genename")) %>%
      left_join(CDS, by = c("ID" = "GeneID")) %>%
      left_join(Protein, by = c("ID" = "GeneID")) %>%
      left_join(UTR5, by = c("ID" = "GeneID")) %>%
      left_join(UTR3, by = c("ID" = "GeneID"))
    
    # Select relevant columns for final data processing
    final_data_raw <- merged_data %>%
      select(
        Genename, log2FC, PTR_AML, RNA_Halflife, Ribosome_density, CDS_GC_percent, EJC_density, CAI, TR, EVI, Translation_Efficiency,
        UTR5_RNA_folding_per_base, UTR5_RNA_length, UTR5_RNA_folding, UTR5_GC_percent, Mean_MRL,
        UTR3_RNA_folding_per_base, UTR3_RNA_length, UTR3_RNA_folding, UTR3_GC_percent, Poly_A_length, MirDB_sites_number, MirDB_Density, Pb_enr,
        SG_enr, Targetscan_sites_number, Context_score, Targetscan_Density, Targetscan_Kd_predicted, Protein_halflife,
        Protein_length, mW, pI, Instability_index, Hydropathy, Aromaticity, 
        Composition_A, Composition_C, Composition_D, Composition_E, Composition_F, 
        Composition_G, Composition_H, Composition_I, Composition_K, Composition_L, 
        Composition_M, Composition_N, Composition_P, Composition_Q, Composition_R, 
        Composition_S, Composition_T, Composition_V, Composition_W, Composition_Y
      )
    
    final_data_raw[,-1] <- lapply(final_data_raw[,-1], as.numeric)
    
    # Calculate rolling median for selected metrics
    metrics <- c("PTR_AML", "RNA_Halflife", "Ribosome_density", "CDS_GC_percent", "EJC_density", "CAI", "TR", "EVI", "Translation_Efficiency",
                 "UTR5_RNA_folding_per_base", "UTR5_RNA_length", "UTR5_RNA_folding", "UTR5_GC_percent", "Mean_MRL",
                 "UTR3_RNA_folding_per_base", "UTR3_RNA_length", "UTR3_RNA_folding", "UTR3_GC_percent", "Poly_A_length", "MirDB_sites_number", "MirDB_Density", "Pb_enr",
                 "SG_enr", "Targetscan_sites_number", "Context_score", "Targetscan_Density", "Targetscan_Kd_predicted", "Protein_halflife",
                 "Protein_length", "mW", "pI", "Instability_index", "Hydropathy", "Aromaticity", 
                 "Composition_A", "Composition_C", "Composition_D", "Composition_E", "Composition_F", 
                 "Composition_G", "Composition_H", "Composition_I", "Composition_K", "Composition_L", 
                 "Composition_M", "Composition_N", "Composition_P", "Composition_Q", "Composition_R", 
                 "Composition_S", "Composition_T", "Composition_V", "Composition_W", "Composition_Y"
    )
    
    for (metric in metrics) {
      # Check if metric exists in data before calculating the median
      if (metric %in% names(final_data_raw)) {
        final_data_raw[[paste0(metric, "_median")]] <- calculate_custom_rolling_median(final_data_raw[[metric]])
      } else {
        warning(paste("Metric", metric, "not found in the dataframe"))
      }
    }
    
    # Extract columns with calculated medians for the final dataset
    median_columns <- grep("_median$", names(final_data_raw), value = TRUE)
    final_median_data <- final_data_raw %>% select(Genename, log2FC, all_of(median_columns))
    final_data(final_median_data)
  
  # Run the statistical analysis
  run_proxy_comparison(final_data_raw)
  })
  
  # Render data table in the UI
  output$table <- renderDT({
    final_data()
  }, options = list(pageLength = 10))
  
  # Render plot based on selected metrics
  output$plots <- renderPlot({
    data <- final_data() %>%
      filter(!is.na(log2FC)) %>%   
      arrange(log2FC)             
    
    data_long <- data %>%
      pivot_longer(cols = starts_with("PTR_AML_median") | starts_with("RNA_Halflife_median") | starts_with("Ribosome_density_median") | starts_with("CDS_GC_percent_median") | starts_with("EJC_density_median") | starts_with("CAI_median") | starts_with("TR_median") | starts_with("EVI_median") | starts_with("Translation_Efficiency_median") | starts_with("UTR5_RNA_folding_per_base_median") | starts_with("UTR5_RNA_length_median") | starts_with("UTR5_RNA_folding_median") | starts_with("UTR5_GC_percent_median") | starts_with("Mean_MRL_median") | starts_with("UTR3_RNA_folding_per_base_median") | starts_with("UTR3_RNA_length_median") | starts_with("UTR3_RNA_folding_median") | starts_with("UTR3_GC_percent_median") | starts_with("Poly_A_length_median") | starts_with("MirDB_sites_number_median") | starts_with("MirDB_Density_median") | starts_with("Pb_enr_median") |
                     starts_with("SG_enr_median") | starts_with("Targetscan_sites_number_median") | starts_with("Context_score_median") | starts_with("Targetscan_Density_median") | starts_with("Targetscan_Kd_predicted_median") | starts_with("Protein_halflife_median") |
                     starts_with("Protein_length_median") | starts_with("mW_median") | starts_with("pI_median") | 
                     starts_with("Instability_index_median") | starts_with("Hydropathy_median") | starts_with("Aromaticity_median") |
                     starts_with("Composition_A_median") | starts_with("Composition_C_median") | starts_with("Composition_D_median") |
                     starts_with("Composition_E_median") | starts_with("Composition_F_median") | starts_with("Composition_G_median") |
                     starts_with("Composition_H_median") | starts_with("Composition_I_median") | starts_with("Composition_K_median") |
                     starts_with("Composition_L_median") | starts_with("Composition_M_median") | starts_with("Composition_N_median") |
                     starts_with("Composition_P_median") | starts_with("Composition_Q_median") | starts_with("Composition_R_median") |
                     starts_with("Composition_S_median") | starts_with("Composition_T_median") | starts_with("Composition_V_median") |
                     starts_with("Composition_W_median") | starts_with("Composition_Y_median"),
                   names_to = "Metric",
                   values_to = "MedianValue")
    
    # Calculate median values for each metric
    medians <- data_long %>%
      group_by(Metric) %>%
      summarise(median_value = median(MedianValue))
    
    # Calculate the difference between the top 200 maximum and minimum values for each metric
    diffs <- data_long %>%
      group_by(Metric) %>%
      summarise(
        median_min_200 = median(MedianValue[order(data$log2FC)][1:min(200, n())], na.rm = TRUE),  # Median of 200 smallest log2FC values
        median_max_200 = median(MedianValue[order(data$log2FC, decreasing = TRUE)][1:min(200, n())], na.rm = TRUE),  # Median of 200 largest log2FC values
        diff_200 = median_max_200 - median_min_200  # Difference between the two median values
      )
    # Add a column for text labels to display the calculated difference
    medians <- medians %>%
      left_join(diffs, by = "Metric")
    
    # Initialize an empty list to store plots based on selected curve types
    plots <- list()
    
    # Generate raw data plot
    if ("raw" %in% input$curve_selection) {
      p_raw <- ggplot(data_long, aes(x = log2FC, y = MedianValue)) +
        geom_line(aes(y = scales::squish(MedianValue, range = c(-1e5, 1e5))), color = "black", alpha = 0.5, na.rm = TRUE) +
        
        # Center and adjust the main title, remove "_median" from facet labels
        labs(title = "Metric Analysis by log2FC", x = "log2FC", y = "Median Value") +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 16),  # Center title
          strip.text = element_text(size = 14),               # Increase facet label size
          panel.spacing = unit(1, "lines"),
          axis.title = element_text(face = "bold", size = 14)
        ) +
        
        # Customize facets, removing "_median" from facet labels and setting rows and columns
        facet_wrap(~ sub("_median$", "", Metric), scales = "free_y", ncol = 3)
      
      # Store the plot in the list
      plots[["raw"]] <- p_raw
    }
    
    # Plot the Spline curve if selected in curve options
    if ("spline" %in% input$curve_selection) {
      p_spline <- ggplot(data_long, aes(x = log2FC, y = MedianValue)) +
        
        # Spline curve with color purple and no confidence interval
        geom_smooth(aes(y = scales::squish(MedianValue, range = c(-1e5, 1e5))), method = "gam", formula = y ~ s(x), color = "purple", se = FALSE, na.rm = TRUE) +
        
        # Add text for the difference in top 200 max and min values, positioned at the top right
        geom_text(data = medians, aes(x = Inf, y = Inf, label = round(diff_200, 2)), 
                  hjust = 1.1, vjust = 1.1, color = "red") +
        
        # Customize facets: remove "_median" from facet labels, set 3 columns per row, and increase space between facets
        facet_wrap(~ sub("_median$", "", Metric), scales = "free_y", ncol = 3) +
        
        # Center the main title, adjust text sizes, and increase space between facets
        labs(title = "Metric Analysis by log2FC (Spline Fit)", x = "log2FC", y = "Median Value") +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 16),  # Center title
          strip.text = element_text(size = 14),               # Increase facet label size
          panel.spacing = unit(1, "lines"),
          axis.title = element_text(face = "bold", size = 14)
        )
      
      # Store the spline plot in the plots list
      plots[["spline"]] <- p_spline
    }
    
    
    # Plot Raw + Median curve if selected in curve options
    if ("both" %in% input$curve_selection) {
      p_both <- ggplot(data_long, aes(x = log2FC, y = MedianValue)) +
        
        # Raw data line in black with transparency
        geom_line(aes(y = scales::squish(MedianValue, range = c(-1e5, 1e5))), color = "black", alpha = 0.5, na.rm = TRUE) +
        
        # Add median horizontal line for each metric in orange, dashed
        geom_hline(data = medians, aes(yintercept = median_value), color = "orange", linetype = "dashed") +
        
        # Add text for difference in top 200 max and min values, positioned at top right
        geom_text(data = medians, aes(x = Inf, y = Inf, label = round(diff_200, 2)), 
                  hjust = 1.1, vjust = 1.1, color = "red") +
        
        # Customize facets: remove "_median" from facet labels, set 3 columns per row, and increase space between facets
        facet_wrap(~ sub("_median$", "", Metric), scales = "free_y", ncol = 3) +
        
        # Center the main title, adjust text sizes, and increase space between facets
        labs(title = "Metric Analysis by log2FC (Raw + Median)", x = "log2FC", y = "Median Value") +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 16),  # Center title
          strip.text = element_text(size = 14),               # Increase facet label size
          panel.spacing = unit(1, "lines"),
          axis.title = element_text(face = "bold", size = 14)
        )
      
      # Store the combined plot in the plots list
      plots[["both"]] <- p_both
    }
    
    # Display selected plots or show a message if none are selected
    if (length(plots) > 0) {
      do.call(grid.arrange, c(plots, ncol = 1))
    } else {
      ggplot() + labs(title = "No Curve Selected") + theme_minimal()
    }
  })
  
  # Download data table as Excel file
  output$downloadTable <- downloadHandler(
    filename = function() { "final_data.xlsx" },
    content = function(file) {
      write_xlsx(final_data(), file)
    }
  )
  
  # Download selected plots as a PDF file
  output$downloadPlots <- downloadHandler(
    filename = function() {
      paste("Plots", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 13, height = 30)
      # Prepare data for plotting by pivoting relevant columns
      data <- final_data() %>%
        filter(!is.na(log2FC)) %>%
        arrange(log2FC)
      data_long <- data %>%
        pivot_longer(cols = starts_with("PTR_AML_median") | starts_with("RNA_Halflife_median") | starts_with("Ribosome_density_median") | starts_with("CDS_GC_percent_median") | starts_with("EJC_density_median") | starts_with("CAI_median") | starts_with("TR_median") | starts_with("EVI_median") | starts_with("Translation_Efficiency_median") | starts_with("UTR5_RNA_folding_per_base_median") | starts_with("UTR5_RNA_length_median") | starts_with("UTR5_RNA_folding_median") | starts_with("UTR5_GC_percent_median") | starts_with("Mean_MRL_median") | starts_with("UTR3_RNA_folding_per_base_median") | starts_with("UTR3_RNA_length_median") | starts_with("UTR3_RNA_folding_median") | starts_with("UTR3_GC_percent_median") | starts_with("Poly_A_length_median") | starts_with("MirDB_sites_number_median") | starts_with("MirDB_Density_median") | starts_with("Pb_enr_median") |
                       starts_with("SG_enr_median") | starts_with("Targetscan_sites_number_median") | starts_with("Context_score_median") | starts_with("Targetscan_Density_median") | starts_with("Targetscan_Kd_predicted_median") | starts_with("Protein_halflife_median") |
                       starts_with("Protein_length_median") | starts_with("mW_median") | starts_with("pI_median") | 
                       starts_with("Instability_index_median") | starts_with("Hydropathy_median") | starts_with("Aromaticity_median") |
                       starts_with("Composition_A_median") | starts_with("Composition_C_median") | starts_with("Composition_D_median") |
                       starts_with("Composition_E_median") | starts_with("Composition_F_median") | starts_with("Composition_G_median") |
                       starts_with("Composition_H_median") | starts_with("Composition_I_median") | starts_with("Composition_K_median") |
                       starts_with("Composition_L_median") | starts_with("Composition_M_median") | starts_with("Composition_N_median") |
                       starts_with("Composition_P_median") | starts_with("Composition_Q_median") | starts_with("Composition_R_median") |
                       starts_with("Composition_S_median") | starts_with("Composition_T_median") | starts_with("Composition_V_median") |
                       starts_with("Composition_W_median") | starts_with("Composition_Y_median"),
                     names_to = "Metric",
                     values_to = "MedianValue")
      
      # Calculate medians and differences for each metric
      medians <- data_long %>% 
        group_by(Metric) %>% 
        summarise(
          median_value = median(MedianValue, na.rm = TRUE),  # Overall median
          diff_200 = median(MedianValue[order(data$log2FC, decreasing = TRUE)][1:min(200, n())], na.rm = TRUE) -  
            median(MedianValue[order(data$log2FC)][1:min(200, n())], na.rm = TRUE)  # Difference using log2FC sorting
        )
      
      # Initialize an empty list to store plots based on selected curve types
      plots <- list()
      
      # Generate raw data plot
      if ("raw" %in% input$curve_selection) {
        p_raw <- ggplot(data_long, aes(x = log2FC, y = MedianValue)) +
          geom_line(aes(y = scales::squish(MedianValue, range = c(-1e5, 1e5))), color = "black", alpha = 0.5, na.rm = TRUE) +
          
          # Center and adjust the main title, remove "_median" from facet labels
          labs(title = "Metric Analysis by log2FC", x = "log2FC", y = "Median Value") +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 16),  # Center title
            strip.text = element_text(size = 11),               # Increase facet label size
            panel.spacing = unit(1, "lines"),
            axis.title = element_text(face = "bold", size = 14)
          ) +
          
          # Customize facets, removing "_median" from facet labels and setting rows and columns
          facet_wrap(~ sub("_median$", "", Metric), scales = "free_y", ncol = 3)
        
        # Store the plot in the list
        plots[["raw"]] <- p_raw
      }
      
      # Plot the Spline curve if selected in curve options
      if ("spline" %in% input$curve_selection) {
        p_spline <- ggplot(data_long, aes(x = log2FC, y = MedianValue)) +
          
          # Spline curve with color purple and no confidence interval
          geom_smooth(aes(y = scales::squish(MedianValue, range = c(-1e5, 1e5))), method = "gam", formula = y ~ s(x), color = "purple", se = FALSE, na.rm = TRUE) +
          
          # Add text for the difference in top 200 max and min values, positioned at the top right
          geom_text(data = medians, aes(x = Inf, y = Inf, label = round(diff_200, 2)), 
                    hjust = 1.1, vjust = 1.1, color = "red") +
          
          # Customize facets: remove "_median" from facet labels, set 3 columns per row, and increase space between facets
          facet_wrap(~ sub("_median$", "", Metric), scales = "free_y", ncol = 3) +
          
          # Center the main title, adjust text sizes, and increase space between facets
          labs(title = "Metric Analysis by log2FC (Spline Fit)", x = "log2FC", y = "Median Value") +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 16),  # Center title
            strip.text = element_text(size = 14),               # Increase facet label size
            panel.spacing = unit(1, "lines"),
            axis.title = element_text(face = "bold", size = 14)
          )
        
        # Store the spline plot in the plots list
        plots[["spline"]] <- p_spline
      }
      
      
      # Plot Raw + Median curve if selected in curve options
      if ("both" %in% input$curve_selection) {
        p_both <- ggplot(data_long, aes(x = log2FC, y = MedianValue)) +
          
          # Raw data line in black with transparency
          geom_line(aes(y = scales::squish(MedianValue, range = c(-1e5, 1e5))), color = "black", alpha = 0.5, na.rm = TRUE) +
          
          # Add median horizontal line for each metric in orange, dashed
          geom_hline(data = medians, aes(yintercept = median_value), color = "orange", linetype = "dashed") +
          
          # Add text for difference in top 200 max and min values, positioned at top right
          geom_text(data = medians, aes(x = Inf, y = Inf, label = round(diff_200, 2)), 
                    hjust = 1.1, vjust = 1.1, color = "red") +
          
          # Customize facets: remove "_median" from facet labels, set 3 columns per row, and increase space between facets
          facet_wrap(~ sub("_median$", "", Metric), scales = "free_y", ncol = 3) +
          
          # Center the main title, adjust text sizes, and increase space between facets
          labs(title = "Metric Analysis by log2FC (Raw + Median)", x = "log2FC", y = "Median Value") +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 16),  # Center title
            strip.text = element_text(size = 14),               # Increase facet label size
            panel.spacing = unit(1, "lines"),
            axis.title = element_text(face = "bold", size = 14)
          )
        
        # Store the combined plot in the plots list
        plots[["both"]] <- p_both
      }
      
      # Print each plot type in the PDF file
      for (plot in plots) {
        print(plot)
      }
      
      dev.off()
    }
  )
  
  # Data to store the state of the highlighted points
  highlighted <- reactiveVal(NULL)
  
  # Observer to update highlighted points when the dataset changes
  observeEvent(input$dataset, {
    req(fgsea_results_list[[input$dataset]])
    fgsea_results <- fgsea_results_list[[input$dataset]]
    highlighted(rep(FALSE, nrow(fgsea_results)))
  })
  
  # Function to toggle labels on a clicked data point in the volcano plot
  toggle_label <- function(index) {
    current_state <- highlighted()
    current_state[index] <- !current_state[index]
    highlighted(current_state)
  }
  
  # Reactive function to retrieve the currently selected dataset
  selected_data <- reactive({
    req(input$dataset, fgsea_results_list[[input$dataset]])
    fgsea_results_list[[input$dataset]]
  })
  
  # Create the volcano plot with ggplot2
  # Reactive data for the volcano plot
  volcano_plot <- reactive({
    data <- selected_data()
    validate(need("pathway" %in% colnames(data), "Error: Missing 'pathway' column."))
    data$significant <- data$padj < 0.05
    ggplot(data, aes(x = NES, y = -log10(padj), text = pathway)) +
      geom_point(aes(color = significant)) +
      geom_text_repel(data = data[highlighted(), ], aes(label = pathway),
                      size = 3, segment.color = NA, box.padding = 0.5) +
      labs(title = "Volcano Plot of fGSEA Results",
           x = "Normalized Enrichment Score (NES)",
           y = "-log10(adjusted p-value)") +
      theme_minimal() +
      scale_color_manual(values = c("black", "turquoise"), labels = c("padj>0.05", "padj<0.05"), name = "padj<0.05")
  })
  
  # Convert the ggplot to an interactive Plotly plot
  output$volcanoPlot <- renderPlotly({
    req(volcano_plot())
    ggplotly(volcano_plot(), tooltip = "text") %>%
      plotly::layout(dragmode = "select") %>%
      event_register("plotly_click")
  })
  
  # Observer to handle clicks on the volcano plot and toggle point labels
  observeEvent(event_data("plotly_click"), {
    click_data <- event_data("plotly_click")
    if (!is.null(click_data)) {
      clicked_point <- click_data$pointIndex + 1  # Ensure correct indexing
      toggle_label(clicked_point)
    }
  })
  
  # Function to download the volcano plot as a PNG image
  output$downloadPlot <- downloadHandler(
    filename = function() { paste("volcano_plot_", input$dataset, ".png", sep = "") },
    content = function(file) {
      req(volcano_plot())
      ggsave(file, plot = volcano_plot(), device = "png", width = 10, height = 6)
    }
  )
  
  # Function to download the volcano plot data
  output$downloadData <- downloadHandler(
    filename = function() { paste("volcano_plot_data_", input$dataset, ".xlsx", sep = "") },
    content = function(file) {
      # Get the data currently displayed in the volcano plot
      data <- selected_data()
      
      # Ensure the data has the relevant columns
      validate(need("pathway" %in% colnames(data), "Error: Missing 'pathway' column."))
      
      # Select relevant columns for reproduction
      volcano_data <- data %>% 
        select(pathway, NES, padj) %>%
        mutate(neg_log10_padj = -log10(padj))  # Add - log10 transformed p-values for plotting
      
      # Write to an Excel file
      write.xlsx(volcano_data, file = file, rowNames = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)