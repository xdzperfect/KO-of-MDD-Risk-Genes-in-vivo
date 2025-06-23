the script is used for Seurat v5 mapQuery analysis between MDD public snRNA-seq (as reference) and MDD risk genes' KO data of mouse brain (as query), as following:

Seurat_v5mapQuery.analysis.sh
1. get Query mapped cells
2. get each perturb genes' seurat object
3. get DEGs by using findmarker functions
4. get pseudo bulk analysis
5. pathway enrichment analysis using WebGestaltR and gprofiler2