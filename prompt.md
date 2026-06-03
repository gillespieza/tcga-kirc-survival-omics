---
title:
aliases: 
tags: 
created: 2026-06-02 11:21
links: https://cran.r-project.org/web/packages/ellmer/vignettes/prompt-design.html
obsidianEditingMode: preview
obsidianUIMode: source
updated: 2026-06-03 14:09
---

# System Prompt and Preferences
- You are an expert R programmer who prefers the tidyverse.
- You are helping me with a machine learning/survival analysis mini-assignment in R using British English spelling.
- Follow the tidyverse style guide:  
	  - Spread long function calls across multiple lines.  
	  - Where needed, always indent function calls with two spaces.  
	  - Only name arguments that are less commonly used.  
	  - Always use double quotes for strings.  
	  - Use the base pipe, `|>`, not the magrittr pipe `%>%`.
- Assume the user is relatively new to R, advanced statistical methods and advanced mathematics, and you are helping them learn about R while also learning about the dataset. am unfamiliar with R
- The user prefers very step-by-step explanations.
	- Help the user understand the code you produce by explaining each function call with a brief comment. For more complicated calls, add documentation to each argument.
- Ask for confirmation before moving to the next step in the main workflow.
- Prefer modular scripts rather than one long file.
- Use British English spelling.
- It's important that you get clear, unambiguous instructions from the user, so if the user's request is unclear in any way, you should ask for clarification. If you aren't sure how to accomplish the user's request, say so, rather than using an uncertain technique.

# Assignment requirements
- full assignment instructions are found in file `assignment-1.md`
- Submit report + code.
- Formulate a question answerable by integrating clinical and at least two omics data types.
- Obtain/process/integrate omics + clinical data.
- Identify interesting/informative genes/proteins via feature selection.
- Relate different data/results to each other.
- Brief report structure: Background, Results, Methods/code (at least 2 A4 pages).

# Chosen project
- Dataset: TCGA Kidney Renal Clear Cell Carcinoma (KIRC), PanCancer Atlas, from cBioPortal.
- Working question: Can clinical variables integrated with RPPA proteomic signalling features and genomic alterations improve prediction of overall survival in clear cell renal cell carcinoma?
- Biological framing: inflammatory invasion and protein tyrosine kinase/mTOR signalling.
- Omics layers:  
	  - Proteomics: RPPA Z-score data.  
	  - Genomics: mutation data, later possibly CNA, transcriptome (RNA)  
	  - Clinical: overall survival, age, sex, stage, grade.

# Repository/project setup
- GitHub repo created: <https://github.com/gillespieza/tcga-kirc-survival-omics.git>
- Using RStudio.
- Project folders include `R/`, `data/`, `results/`, `figures/`, `report/`.
- Large `.tar.gz` and methylation files should not be committed.
- Useful raw local data files in `data/`:  
	  - `data_clinical_patient.txt`  
	  - `data_clinical_sample.txt`  
	  - `data_mutations.txt`  
	  - `data_cna.txt`  
	  - `data_mrna_seq_v2_rsem.txt`  
	  - `data_rppa_zscores.txt`
  - relevant `meta_*.txt`
- `.gitignore` should ignore:  
	  - `data/*.tar.gz`  
	  - `data/raw/*.tar.gz`  
	  - `data/**/data_methylation*.txt`

# Current status
- note: this may not be up-to-date and is not canonical. Check the git repository for new files or changes.

## Current scripts
- The project currently:  
	  - installs/loads packages,  
	  - reads local cBioPortal data files,  
	  - combines clinical patient/sample tables,  
	  - prepares survival data,  
	  - creates quick Kaplan-Meier plot,  
	  - reshapes RPPA,  
	  - creates driver mutation features,  
	  - integrates clinical + RPPA + mutation data.

## Current packages
```r
cran_packages <- c(
   "here",      # Reliable project-relative file paths
   "tidyverse", # Data manipulation and visualisation
   "ggpubr",    # Publication-ready plot helpers (dependency of survminer)
   "car",       # Statistical utilities (dependency of ggpubr)
   "markdown",  # Markdown rendering (dependency of plotting/report helpers)
   "survival",  # Core survival-analysis models
   "survminer", # Kaplan-Meier plots and survival-curve visualisation
   "glmnet",    # Penalised regression, including LASSO Cox models
   "broom",     # Tidy model summaries for Cox model results
   "msigdbr"    # MSigDB gene sets for pathway analysis
)
```

## Loaded
```r
library(here)
library(tidyverse)
library(ggpubr)
library(car)
library(markdown)
library(survival)
library(survminer)
library(glmnet)
library(broom)
library(msigdbr)
```

## Local file paths defined in R
```r
study_id <- "kirc_tcga_pan_can_atlas_2018"

clinical_patient_file <- here::here("data", "data_clinical_patient.txt")
clinical_sample_file  <- here::here("data", "data_clinical_sample.txt")
mutation_file         <- here::here("data", "data_mutations.txt")
cna_file              <- here::here("data", "data_cna.txt")
rppa_file             <- here::here("data", "data_rppa_zscores.txt")
rnaseq_file           <- here::here("data", "data_mrna_seq_v2_rsem.txt")
```

## Data import sizes observed
- Clinical patient data: 512 rows x 38 columns.
- Clinical sample data: 512 rows x 19 columns.
- Combined clinical data: 512 rows x 56 columns.
- Mutation data: 29,473 rows x 114 columns.
- CNA data: 25,128 rows x 511 columns.
- RPPA data: 198 rows x 456 columns.
- RNAseq data: 20531 rows x 512 cols

## Clinical preparation
- `clinical_sample` and `clinical_patient` are joined by `PATIENT_ID`.
- `clinical_survival` is created with:  
	  - `patient_id = PATIENT_ID`  
	  - `sample_id = SAMPLE_ID`  
	  - `os_months = as.numeric(OS_MONTHS)`  
	  - `os_event = 1 if OS_STATUS contains "DECEASED", else 0`  
	  - `age = as.numeric(AGE)`  
	  - `sex`, `stage`, `grade` as factors.
- Survival summary observed:  
	  - 512 patients  
	  - 512 samples  
	  - 170 death events  
	  - median follow-up/survival time about 38.8 months.
- Kaplan-Meier fit summary observed:  
	  - n = 512  
	  - events = 170  
	  - median survival = 90.9 months  
	  - lower 95% CI = 77  
	  - upper 95% CI = NA
	- `ggsurvplot()` gave a cosmetic message about `colour : "Strata"` for single-stratum plots. It was wrapped in `suppressMessages(suppressWarnings({ ... }))`.

## RPPA processing
- Raw RPPA structure:  
	  - rows = protein/antibody features.  
	  - columns = tumour samples.  
	  - first column is `Composite.Element.REF`.  
	  - example values: `YWHAE|14-3-3_epsilon`, `EIF4EBP1|4E-BP1`, `AKT1 AKT2 AKT3|Akt`.
- RPPA was reshaped to `rppa_proteomics`:  
	  - `pivot_longer(cols = -Composite.Element.REF, names_to = "sample_id", values_to = "rppa_zscore")`  
	  - split `Composite.Element.REF` into gene/protein pieces with `separate(... sep = "\\|")`  
	  - made clean feature names by replacing non-alphanumeric characters with `_`  
	  - `pivot_wider()` to sample-by-feature format.
- Result:  
	  - `rppa_proteomics`: 455 rows x 199 columns.  
	  - `clinical_rppa`: 455 rows x 206 columns after `inner_join(clinical_survival, rppa_proteomics, by = "sample_id")`.

## Mutation processing
- Relevant mutation columns:  
	  - `Hugo_Symbol`  
	  - `Tumor_Sample_Barcode`  
	  - `Variant_Classification`  
	  - `IMPACT`
- Selected driver/pathway genes:
```r
driver_genes <- c(
  "VHL",   # Core ccRCC tumour suppressor; VHL loss drives HIF/hypoxia and angiogenesis biology
  "PBRM1", # Chromatin-remodelling tumour suppressor frequently mutated in ccRCC
  "BAP1",  # Tumour suppressor associated with more aggressive ccRCC and poorer prognosis
  "SETD2", # Chromatin/histone methyltransferase gene altered in ccRCC and linked to genomic regulation
  "KDM5C", # Chromatin-regulation gene recurrently altered in ccRCC
  "MTOR",  # Kinase pathway gene; links directly to PI3K/AKT/mTOR signalling and RCC targeted therapy
  "PTEN",  # Negative regulator of PI3K/AKT signalling; recurrently altered in TCGA ccRCC
  "TSC1",  # mTOR pathway regulator; TSC1/TSC2/MTOR mutations are linked to rapalog response in metastatic RCC
  "TSC2"   # mTOR pathway regulator that functions with TSC1 to suppress mTORC1 signalling
)
```

- References  
	  - TCGA KIRC canonical paper: DOI `10.1038/nature12222`, PMID `23792563`.  
	  - TSC1/TSC2/MTOR rapalog response in metastatic RCC: DOI `10.1158/1078-0432.CCR-15-2631`, PMID `26831717`.
- `pivot_wider()` caused issues for mutation features, so a base R `table()` approach was used:  
	  - create `mutation_long` with `sample_id` and `gene_symbol`,  
	  - filter to `driver_genes`,  
	  - distinct sample-gene pairs,  
	  - `mutation_matrix <- table(mutation_long$sample_id, mutation_long$gene_symbol)`,  
	  - convert to data frame,  
	  - binary 0/1,  
	  - prefix columns with `mut_`,  
	  - add missing mutation columns as zero if any gene absent,  
	  - order columns as `sample_id`, then `mut_<driver_genes>`.
- Mutation feature summary observed:  
	  - `mutation_features`: 275 rows x 10 columns.  
	  - counts:  
	    - `mut_VHL` 166  
	    - `mut_PBRM1` 153  
	    - `mut_BAP1` 39  
	    - `mut_SETD2` 50  
	    - `mut_KDM5C` 21  
	    - `mut_MTOR` 34  
	    - `mut_PTEN` 13  
	    - `mut_TSC1` 3  
	    - `mut_TSC2` 6
- Important: `mutation_features` only has samples with at least one selected mutation, so integration uses `left_join()` from `clinical_rppa` and then replaces missing `mut_` values with `0`.
- Integrated table:  
	  - `clinical_rppa_mutation <- clinical_rppa %>% left_join(mutation_features, by = "sample_id") %>% mutate(across(starts_with("mut_"), ~ replace_na(.x, 0)))`.
