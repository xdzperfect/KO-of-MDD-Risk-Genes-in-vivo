# In vivo Perturb-seq for MDD risk genes

These notebooks contain code reproducing the single-cell analyses from:

> Liansheng Zhang, Xiangrui Kong, Qi Ma, Xinde Hu, Shicheng Cai, Bo Wang, Weijuan Zou, Tao Bai, Meimei Zhang, Liu Fan, Runlin Tan, Ziji Dai, Zhiheng Jia, Tianwen Li, Xingyu Liu, Huatai Xu, Jianrong Wu, Yuanyi Zheng, Zhengzheng Xu, and Haibo Zhou, "Linking GWAS risk genes to transcriptional features of major depressive disorder via *in vivo* Perturb-seq", Nature Genetics, 2026

## Defining Cells with Multiplicity of Infection (MOI) = 1

Please refer to the `sg_assignment` directory for details.

In this directory, we demonstrated the MOI calculation workflow for one of the reactions. Other reactions can be processed similarly. Briefly, since our project adopts a dual-sgRNA mode for KO Perturb-seq, we define MOI=1 cells primarily based on cells defined by sg2, and then remove cells with multiple sg1 assignments. For details, refer to **Step 4: "Identify confident MOI=1 cells and assign guides"** in the `sg_assignment.Rmd` within that directory.

> **Note:** It is important to emphasize that this part of the analysis is largely based on the analytical code from the AAV-Perturb-seq project published in Nature by Santinha, A.J, et al. in 2023. More specifically, refer to https://github.com/plattlab/AAV-Perturb-seq. We strongly recommend familiarizing yourself with their scripts and rationale first.

After defining MOI=1 cells, all subsequent analyses are built upon this foundation.

## Application 1: Gene co-expression network analysis

**To reproduce Fig. 2a-2c and Extended Data Fig. 3b-3m**

This part of the analysis was conducted following the online protocol for hdWGCNA. Specifically, we used only MOI=1 cells as the input data, focusing primarily on neuron cells. The analysis strictly adhered to the parameters and example code provided in the hdWGCNA protocol.

Upon completion of this analysis, the results generated will correspond to Figure 2a, 2b, and Extended Data Fig. 3b-3c. Furthermore, the analysis will yield the enrichR annotation results shown in Extended Data Fig. 3d-3m. These enrichment results for GO pathways can then be visualized using ggplot2.

**Implementation details:**
- Script: `hdWGCNA.analysis.r` located in the `hdWGCNA_MDD_risk_genes_perturbation_with_Module` directory.

To investigate the perturbation effect of each perturbed gene relative to the control across different gene co-expression modules, we performed a deeper analysis of how these perturbations affect the activity of individual modules.

**Implementation details:**
- Script: `MDD_risk_gene_knockout_on_module_activity.r` located in the `hdWGCNA_MDD_risk_genes_perturbation_with_Module` directory.
- The results from this analysis correspond to Figure 2c.

## Application 2: Specific pathway similarity analysis of MDD risk gene perturbations

**To reproduce Fig. 2d-2f**

This analysis is primarily based on the DEGs (Differentially Expressed Genes) derived from comparing MDD risk gene perturbations versus control within the MOI=1 cell population. The methodological framework is built upon the DEG analysis code from the AAV-Perturb-seq project published by Santinha, A.J, et al. in Nature (2023).

**Reference demonstration:** Code example provided in the directory: `specific_pathway_cluster_analysis/DEGs_demo_data_and_code/`

Upon obtaining the DEG files, the analysis proceeds with the following scripts:

- **Shell Script:** `2025-6-18.specific_pathway.correlation.sh`
- **R Script:** `specific-pathway.r`

### Key Steps in the Analysis:

1. **Data Preparation:** The DEG files for each MDD risk gene perturbation (compared to control) are loaded and consolidated.
2. **Pathway Extraction:** Specific pathways of interest are extracted from databases such as KEGG, GO, or Reactome, focusing on pathways relevant to MDD pathophysiology and neuronal function.
3. **Similarity Calculation:** For each MDD risk gene perturbation, compute enrichment scores or regulation patterns across these specific pathways.
4. **Correlation Matrix Construction:** Generate a similarity matrix by comparing the pathway regulation profiles of every pair of MDD risk genes.
5. **Clustering and Visualization:** Apply hierarchical clustering to the correlation matrix and visualize results using clustered heatmaps (Fig. 2d-2f).
6. **Validation and Interpretation:** Assess biological coherence through functional enrichment analysis.

## Application 3: Identify neuron cell subtype-specific DEGs

**To reproduce Figure 2g-2h**

The analysis begins by identifying neuronal cell subtypes—specifically excitatory and inhibitory neuronal subpopulations—based on established marker genes (visually detailed in Extended Data Fig. 6a-6b).

Once neuronal subtypes are annotated, calculate subtype-specific Differential Expressed Genes (DEGs) by comparing MOI=1 cells with MDD risk gene knockouts against control cells.

**Implementation details:**
- Script: `get_specific_DEGs.r` located in the `subcell_specific_DEGs_in_Ex_neuron_and_In_neuron/` directory.
- This script processes DEG outputs to identify subtype-specific transcriptional responses and generates Figure 2g-2h.

## Application 4: Cluster and reduce dimensionality of 44 MDD risk genes' perturbation in neurons

**To reproduce Fig. 3a and Extended Data Fig. 4f**

The `perturbation_clustering` directory provides the core scripts:

- `direct_cluster.rmd`
- `umap_cluster.rmd`

**Analysis workflow:**
1. **Data Integration:** Compile transcriptional perturbation profiles for all 44 MDD risk gene knockouts.
2. **Dimensionality Reduction:** Apply UMAP for visualization of perturbation relationships.
3. **Clustering Analysis:** Identify groups of MDD risk genes with similar transcriptional responses.
4. **Visualization:** Generate plots for Figure 3a and Extended Data Fig. 4f.

## Application 5: Map mouse perturbation data back to MDD patients' single-cell data

**To reproduce Fig. 3g**

This cross-species integration analysis uses scripts in the `Seurat_v5mapQuery_MDD_and_mouse_MDD_risk_gene_KO` directory:

- `Seurat_v5mapQuery.analysis.sh`
- `get_mapped_cell.r`

**Analysis workflow:**
1. **Reference Construction:** Use human MDD patient snRNA-seq data as reference.
2. **Query Mapping:** Project mouse perturbation profiles onto human reference using Seurat v5 `mapQuery()` function.
3. **Analysis & Interpretation:** Analyze mapping results to generate Figure 3g.

## Application 6: Correlate six MDD risk genes in mouse neurons to MDD patients

**To reproduce Fig. 6e-f and Extended Data Fig. 10b-g, and j**

The correlation analysis is implemented in the `MDD_patients_correlation_analysis` directory:

- Script: `MDD_correlation_with_mouse_pool_screening.sh`

**Analysis workflow:**
1. **Focus Selection:** Concentrate on six MDD risk genes (cluster1) based on prior results.
2. **Human Gene Expression Profiling:** Extract expression profiles of human orthologs from MDD patient data.
3. **GSEA Enrichment Analysis:** Evaluate enrichment in Oxytocin signaling pathway.
4. **Visualization:** Generate GSEA plots for Figure 6e-6f and Extended Data Fig. 10b-g, and j.

## Special Note

This project did not develop new algorithms and primarily utilizes previously published methods:

- **snRNA-seq analysis:** Based on Hao et al., 2023 (Seurat)
- **gRNA expression libraries:** Count tables generated as shown in Hill et al., 2018
- **MOI calling and DEGs analysis:** Based on Santinha, A.J., et al., 2023
- **hdWGCNA analysis:** Based on Morabito et al., 2023

The analysis process for other array data is consistent with that of the pools data and is not provided here. For further inquiries, please contact us.

## Data Availability

Download the R objects from NGDC to use the notebooks: **OMIX010695** (https://ngdc.cncb.ac.cn/omix/preview/AUGqe6Cm)

As the project has not yet been formally published, the raw data will be automatically released upon publication and can be retrieved from the NGDC database under project **PRJCA031709** (GSA ID: CRA020051).
