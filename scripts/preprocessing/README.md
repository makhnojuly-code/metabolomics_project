# Preprocessing scripts (LC-MS, neg mode)

This directory contains scripts for raw metabolomics preprocessing.

## Order of execution

1) 00_setup.R  
   loads libraries, configures project paths, sets color palettes.

2) 01_create_sample_sheet_neg.R  
   reads filenames from data/raw/neg_mode and generates sample sheet
   (sample / blank / QC detection based on file naming)

3) 01_preprocessing.R  
   - import raw .mzML (onDisk mode)
   - peak picking (CentWave)
   - RT alignment (Obiwarp)
   - feature grouping
   - filling missing intensities
   - export matrices → data/processed/*.csv
   - export QC summary tables → results/qc/

4) 01a_qc_plots_fast.R  
   quick QC plots: RT shift / BPC subset / TIC boxplot

5) 01b_qc_plots_full.R  
   slow but best visualization: full BPC traces + pastel colors

---
outputs:
- intermediate RDS files: data/intermediate/
- final feature tables: data/processed/
- QC plots + reports: results/qc/