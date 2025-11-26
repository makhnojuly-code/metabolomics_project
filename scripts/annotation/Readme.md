NEG Annotation Pipeline

Script 01 combines PLS-DA VIP values, limma differential statistics, and feature information (mz and rt).  
It produces a full table of NEG-mode candidate features and a simplified file prepared for MetaboAnalyst.  
No annotation is performed at this stage. The script only prepares and structures the input data.

Script 02 converts the simplified file into the text format required by the Mummichog module  
(columns `m.z`, `p.value`, `t.score`, `r.t`).  
The generated txt file is uploaded manually to MetaboAnalyst.

After running Mummichog, MetaboAnalyst produces two key outputs:

- `mummichog_matched_compound_all.csv` (matched compounds)
- `mummichog_pathway_enrichment_integ.csv` (pathway enrichment)

These files are placed into `results/annotation`.  
The annotation process starts from this point.

Script 03 integrates the Mummichog results with the candidate table created in Script 01.  
It uses rounded mz and rt values to reliably match LC–MS features with KEGG compounds.  
The script produces a full combined table and a smaller table with only the features that received KEGG matches.  
This is the first step where LC–MS features are connected to specific metabolites.

Script 04 uses KEGGREST to download chemical information for each KEGG ID detected by Mummichog:  
compound name, molecular formula, and exact mass.  
This information is merged into the feature table.  
Now the table contains statistical results (limma), VIP values, mz/rt, KEGG IDs and KEGG chemical properties.  
This creates a complete dataset for biological interpretation and pathway analysis.

Script 05 identifies the core set of stress-related metabolic signals.  
It filters features using thresholds for `adj.P.Val`, `|logFC|`, and `VIP`, producing `core_features_neg`.  
Then it aggregates data by KEGG ID to get metabolite-level results (`core_metabolites_neg`).  
For each metabolite, it calculates:

- number of supporting features  
- minimum FDR  
- average and maximum logFC  
- maximum VIP  

Each metabolite is assigned a direction (“Up in stress” or “Down in stress”).  
Short names and simple stress-related groups (e.g., JA/oxylipins, flavonoids, aromatic amino acids) are also added.

Script 06 performs a quality check of the final tables.  
It verifies:

- required columns are present  
- no incorrect KEGG duplicates  
- NA values are reasonable  
- direction labels match `mean_logFC`  
- pathway information is consistent (if available)

The results are written to `validation_report_neg.txt`.  
This confirms the tables are ready for visualization and biological interpretation.