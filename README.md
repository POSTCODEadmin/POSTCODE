# POSTCODE

**POSTCODE** is a computational framework designed to characterize and interpret
post-transcriptional regulatory mechanisms by integrating transcriptomic,
translational, and protein-related features.

The associated **Shiny web application** allows users to upload gene-level data
(e.g. log2 fold-changes) and explore how intrinsic and extrinsic transcript
features relate to protein expression and post-transcriptional regulation.

🔗 **Live Shiny application (recommended)**  
👉 https://postcode-lab.shinyapps.io/POSTCODE/

---

## 🚀 Features

- Annotation of gene-level input with post-transcriptional regulatory proxies
- Rolling-median, raw, and spline-based trend visualization
- Integrated fGSEA analysis across multiple regulatory layers
- Downloadable tables and figures (Excel, PDF, PNG)
- Robust, deployment-ready Shiny architecture

---

## 📥 Input format

The application expects a **tab-delimited `.txt` file** with two columns:

| Column | Description |
|------|-------------|
| Genename | Gene symbol (HGNC, uppercase recommended) |
| log2FC | Log2 fold-change value |

Example:

```text
TP53    -1.25
MYC      2.10
