# Microbiome–Chemotherapy Toxicity Analyzer

## Overview
This project provides a Python-based command-line tool for exploratory analysis of associations between gut microbiome features and chemotherapy-related toxicity severity in pediatric leukemia patients.

The tool integrates microbiome abundance data with clinical metadata in order to identify bacterial taxa (or other molecular features, such as metabolites) that show differential abundance between patients experiencing mild versus severe treatment-related symptoms.

The analysis pipeline is fully automated, reproducible, and configurable via command-line arguments, making it suitable for repeated use across different datasets and research questions.

> **Important note:**  
> This tool is intended for exploratory data analysis and hypothesis generation only.  
> It does not provide clinical recommendations and should not be used for diagnostic or therapeutic decision-making.

---

## Scientific Motivation
Chemotherapy-induced toxicity is a major challenge in pediatric leukemia treatment. While some patients experience relatively mild side effects, others suffer from severe toxicity that can limit treatment intensity or negatively impact quality of life.

Recent research suggests that the gut microbiome may modulate host responses to chemotherapy, influencing inflammation, drug metabolism, and immune function. However, initial microbiome analyses in many labs are often performed using ad hoc scripts or manual workflows, which can be error-prone and difficult to reproduce.

This project aims to address this gap by providing a standardized, transparent, and reusable analysis pipeline that:
- Performs consistent statistical comparisons
- Corrects for multiple hypothesis testing
- Produces interpretable visual summaries
- Can be easily reused on new datasets without modifying the source code

---

## Input Data

The tool requires **two input files**, provided as CSV or TSV.

### 1. Microbiome Abundance Table
A table where:
- Each row corresponds to a biological sample
- Each column (except `SampleID`) represents a bacterial taxon or molecular feature
- Values represent relative abundances or normalized counts

**Required column:**
- `SampleID`

Example:
```text
SampleID,Bacteroides,Faecalibacterium,Escherichia
S01,0.30,0.12,0.01
S02,0.05,0.35,0.08
S03,0.28,0.10,0.02
2. Clinical Metadata Table
A table containing sample-level clinical or experimental information.

Required columns:

SampleID

A grouping column defining symptom severity (default: Severity)

Example:

SampleID,Severity,Age
S01,Mild,7
S02,Severe,9
S03,Mild,6
The grouping column can be customized (e.g. toxicity grade, response category) using command-line arguments.

Analysis Workflow
The analysis pipeline consists of the following steps:

Input validation

Verification of required columns

Detection of mismatched or missing sample identifiers

Data integration

Merging of microbiome abundance data with clinical metadata

Group definition

Separation of samples into two user-defined groups (e.g. Mild vs Severe)

Statistical testing

Feature-wise comparison using the Mann–Whitney U test

Calculation of mean abundance per group

Computation of log2 fold change

Multiple hypothesis correction

Adjustment of p-values using the Benjamini–Hochberg False Discovery Rate (FDR) method

Generation of q-values for downstream interpretation

Visualization

Volcano plot summarizing effect size and statistical significance

Boxplots for the most significant features

Output
All results are written to a user-specified output directory.

Generated files include:
results_table.csv
A summary table containing one row per feature, with the following columns:

Mean abundance in each group

log2 fold change (group2 vs group1)

Raw p-value

FDR-adjusted q-value

volcano_plot.png
A volcano plot where:

The x-axis represents log2 fold change

The y-axis represents −log10(p-value)

Vertical lines indicate fold-change thresholds

Feature labeling can be controlled via command-line options

Statistical significance is determined using FDR (q-values)

boxplot_<feature>.png
Boxplots showing the distribution of abundance values across groups for the top-ranked features.

Installation
Clone the repository and install the required dependencies:

git clone https://github.com/Danielahuppert/Microbiome-Chemotherapy-Toxicity-Analyzer.git
cd Microbiome-Chemotherapy-Toxicity-Analyzer
pip install -r requirements.txt
Usage
Basic usage (example data)
python analyze_microbiome.py
Custom input files
python analyze_microbiome.py \
  --abundance data/my_microbiome.csv \
  --metadata data/my_metadata.csv
TSV input files
python analyze_microbiome.py \
  --abundance data/my_microbiome.tsv \
  --metadata data/my_metadata.tsv \
  --sep "\t"
Volcano plot labeling options
Label top N features by statistical significance:

python analyze_microbiome.py --label_mode top --top_n 10
Label only features passing FDR and fold-change thresholds:

python analyze_microbiome.py --use_q --q_thresh 0.1 --fc_thresh 1
Interpretation Guidelines
Features with low q-values represent statistically significant associations after correcting for multiple testing

Effect size (log2 fold change) should be interpreted alongside significance

Results should be validated in independent cohorts or follow-up experiments

Scope and Limitations
Designed for exploratory analysis

Limited statistical power for small sample sizes

Does not model confounding variables

Not suitable for clinical decision-making

Future Extensions
Possible future improvements include:

Feature prevalence filtering

Support for longitudinal or repeated-measures data

Automated HTML or PDF report generation

Integration of diversity metrics and ordination analyses

Course Information
This project was developed as part of a Python programming course.


