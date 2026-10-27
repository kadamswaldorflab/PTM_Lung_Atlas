# Pigtail Macaque Lung Single-Cell Atlas: Analysis Pipeline

Code and parameters for processing, integrating, annotating and comparing the pigtail macaque lung scRNA-seq atlas (61,323 cells, 37 cell types, 5 animals). This repository accompanies the Methods section of [manuscript citation].

## Repository contents

| File or folder | What it does |
|---|---|
| `Lung_Atlas.ipynb` | QC, doublet removal, scVI and scANVI integration, UMAP and Leiden clustering. Writes `scanvi_output_v2.h5ad` |
| `Lung_Atlas_Analysis.Rmd` | Applies the annotation, exports the atlas for the IM step, adds the cross-species labels and saves `lung_final.RDS` |
| `IM_Subclustering.ipynb` | Splits interstitial macrophages (IM) into Mo-IM, IM and cDC2 and removes two subclusters |
| `Cross_Species_Reference_HLCA.ipynb` | Maps the atlas to the Human Lung Cell Atlas (HLCA) with scArches |
| `Cross_Species_Reference_RIRA.Rmd` | Transfers labels from the Rhesus Immune Reference Atlas (RIRA) to the immune cells with Seurat |
| `Lung_Atlas_Figures.Rmd` | Makes all figures from `lung_final.RDS` |
| `final_cell_labels.csv` | Final cell type of each of the 61,323 cells |
| `removed_cells.csv` | Reason for removing each of the 7,541 removed cells |
| `IM/` | Original IM subcluster of each IM cell |
| `IM_subclustering_output/` | IM labels and the scANVI latent space from `IM_Subclustering.ipynb` |
| `HLCA/` | Mapped cells, the ortholog table and the HLCA label transfer results |
| `RIRA/` | Reference and query cell lists and the RIRA label transfer results |

Large intermediate objects (`.h5ad`, `.RDS`) and trained models are not included. Folder paths are set at the top of each notebook and Rmd.

## Run order

1. `Lung_Atlas.ipynb`
2. `Lung_Atlas_Analysis.Rmd`, part 1
3. `IM_Subclustering.ipynb`
4. `Cross_Species_Reference_HLCA.ipynb` and `Cross_Species_Reference_RIRA.Rmd`. These only need `scanvi_output_v2.h5ad` and the cell lists in `HLCA/` and `RIRA/`, so they can run any time after step 1.
5. `Lung_Atlas_Analysis.Rmd`, part 2
6. `Lung_Atlas_Figures.Rmd`

## Software versions

| Tool | Version |
|---|---|
| Cell Ranger | 9.0.1 (reference built with mkref 7.2.0) |
| CellBender | 0.3.0 |
| Python | 3.10.19 |
| Scanpy | 1.11.5 |
| anndata | 0.11.4 |
| Scrublet | 0.2.3 |
| scvi-tools | 1.3.3 |
| scArches | 0.6.1 |
| leidenalg | 0.10.2 |
| umap-learn | 0.5.9.post2 |
| R | 4.4.0 |
| Seurat | 5.5.0 |

Random seeds are set in every notebook and Rmd (42 unless stated otherwise).

## Pipeline overview

```
FASTQ
  -> 1. Cell Ranger alignment (Mmul10, Ensembl 113 protein-coding)
  -> 2. CellBender ambient RNA removal (per sample)
  -> 3. QC filtering + Scrublet doublet removal (per sample)
  -> 4. HVG selection -> scVI integration -> scANVI refinement
  -> 5. UMAP + Leiden clustering
  -> 6. Conversion to Seurat, marker-based annotation and subclustering (35 cell types)
  -> 7. IM subclustering (37 cell types)
  -> 8. Cross-species mapping (HLCA via scArches, RIRA via Seurat label transfer)
  -> 9. Figures
```

## 1. Reference and alignment

- Genome is Mmul10 (GenBank GCA_003339765.3)
- Annotation is Ensembl release 113, filtered to protein-coding genes
- Reference built with `cellranger mkref` (v7.2.0), reads aligned with `cellranger count` (v9.0.1)

Rhesus macaque was used as the reference for pigtail macaque samples because of its more complete assembly and annotation.

```bash
cellranger mkgtf Macaca_mulatta.Mmul_10.113.gtf Mmul10.113.pc.gtf \
  --attribute=gene_biotype:protein_coding

cellranger mkref --genome=Mmul10_ens113_pc \
  --fasta=Macaca_mulatta.Mmul_10.dna.toplevel.fa \
  --genes=Mmul10.113.pc.gtf

cellranger count --id=<sample> --transcriptome=Mmul10_ens113_pc \
  --fastqs=<fastq_dir> --sample=<sample> --create-bam false
```

## 2. Ambient RNA removal

CellBender `remove-background` on each raw feature-barcode matrix with default parameters.

```bash
cellbender remove-background --input raw_feature_bc_matrix.h5 \
  --output cellbender.h5 --cuda
```

## 3. Quality control

Mitochondrial genes in the Mmul_10 annotation lack an `MT-` prefix, so the mitochondrial fraction is computed from the 13 mitochondrially encoded protein-coding genes.

| Filter | Threshold |
|---|---|
| Detected genes | > 200 and < 7,500 |
| UMIs | > 500 |
| Mitochondrial fraction | < 2.5% |
| Doublets | Scrublet per sample (expected rate 0.05, threshold 0.25) |

```python
mt_genes = ["ND1","ND2","ND3","ND4","ND4L","ND5","ND6",
            "ATP6","ATP8","COX1","COX2","COX3","CYTB"]
adata.var["mt"] = adata.var_names.isin(mt_genes)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

adata = adata[(adata.obs.n_genes_by_counts > 200) &
              (adata.obs.n_genes_by_counts < 7500) &
              (adata.obs.total_counts > 500) &
              (adata.obs.pct_counts_mt < 2.5)].copy()

import scrublet as scr
scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.05)
scores, _ = scrub.scrub_doublets()
predicted = scrub.call_doublets(threshold=0.25)
adata = adata[~predicted].copy()
```

68,864 cells from the 5 animals passed QC (21,591 genes).

## 4. Integration

**Highly variable genes** from raw counts.

```python
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000,
                            batch_key="sample", layer="counts")
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
```

**scVI** on raw counts of the HVGs, with sample as batch.

| Parameter | Value |
|---|---|
| n_layers | 2 |
| n_hidden | 128 |
| n_latent | 50 |
| gene_likelihood | nb |
| max_epochs | 116 |

```python
scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts", batch_key="sample")
vae = scvi.model.SCVI(adata_hvg, n_layers=2, n_hidden=128,
                      n_latent=50, gene_likelihood="nb")
vae.train(max_epochs=116, train_size=0.9, early_stopping=True,
          early_stopping_patience=20,
          plan_kwargs={"lr": 1e-3, "weight_decay": 1e-4})
adata_hvg.obsm["X_scvi"] = vae.get_latent_representation()
```

**scANVI** seeded with preliminary Leiden clusters (resolution 0.5, 30-NN graph on the scVI latent space).

```python
sc.pp.neighbors(adata_hvg, use_rep="X_scvi", n_neighbors=30)
sc.tl.leiden(adata_hvg, resolution=0.5, key_added="leiden_initial", random_state=42)

lvae = scvi.model.SCANVI.from_scvi_model(vae, labels_key="leiden_initial",
                                         unlabeled_category="Unknown")
lvae.train(max_epochs=10, early_stopping=True, early_stopping_patience=15,
           plan_kwargs={"lr": 5e-3, "weight_decay": 1e-4})
adata_hvg.obsm["X_scanvi"] = lvae.get_latent_representation()
```

## 5. Embedding and clustering

```python
sc.pp.neighbors(adata_hvg, use_rep="X_scanvi", n_neighbors=30)
sc.tl.umap(adata_hvg, min_dist=0.3, random_state=42)
sc.tl.leiden(adata_hvg, resolution=0.5, random_state=42)
```

**Integration metric.** Sample mixing is the inverse Simpson index of sample identity among each cell's 30 nearest neighbors. It is computed in `Lung_Atlas_Figures.Rmd` on the first 50 principal components of the unintegrated data (4,000 variable genes) and on the scANVI latent space. The median increased from 1.40 (unintegrated) to 2.27 (integrated).

## 6. Annotation (R/Seurat)

The AnnData object is converted to Seurat with `anndata` and `reticulate`. The 28 Leiden clusters were annotated iteratively using cluster markers and canonical lineage markers. Heterogeneous clusters were subclustered on the scVI latent space, with the number of neighbors (10 to 30) and the resolutions (0.02 to 1.0) chosen for each cluster.

```r
markers <- FindAllMarkers(obj, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25)

# Subclustering of a heterogeneous cluster
sub <- FindNeighbors(sub, reduction = "scvi", dims = 1:50, k.param = 30)
sub <- FindClusters(sub, resolution = c(0.02, 0.05, 0.1, 0.2))
```

The outcome of the annotation is recorded in `final_cell_labels.csv` and `removed_cells.csv`. `Lung_Atlas_Analysis.Rmd` applies these files to the integrated object to rebuild the annotated atlas.

Low-quality, ambient-dominated and doublet cells were removed in two stages. 7,397 cells were removed during annotation, which gave 61,467 cells across 35 cell types. 144 more were removed in the IM step (section 7).

| Reason in `removed_cells.csv` | Cells | Evidence |
|---|---|---|
| Ambient RNA | 3,172 | Initial clusters 9 and 12. Very few genes, 0% mitochondrial reads and no specific markers |
| Neutrophil/T-NK doublets | 378 | Initial cluster 26. Highest doublet score among the initial clusters, with neutrophil and T/NK markers in the same cells |
| AM doublets/damaged | 1,490 | Alveolar macrophage subclusters with endothelial genes (doublets), mesothelial and lymphocyte genes (contamination), or a loss of abundant cytoplasmic transcripts (damaged cells) |
| NK/B low quality | 1,265 | NK and B cells with about 420 median genes, expressing only their most abundant transcripts |
| Sparse endothelium | 522 | Endothelial cluster with the lowest gene counts in the atlas and only broadly expressed endothelial markers |
| CD8 low quality/RBC | 425 | CD8 T subcluster with half the genes of the main subcluster, and a small subcluster carrying hemoglobin transcripts |
| AT2 low quality | 145 | AT2 subcluster with low gene counts |
| IM low quality | 80 | IM subcluster 7 (section 7) |
| IM possible doublets | 64 | IM subcluster 8 (section 7) |
| **Total** | **7,541** | |

## 7. IM subclustering

The IM cluster (1,680 cells) was split in `IM_Subclustering.ipynb`. A scANVI model was trained on the full annotated atlas, and the IM cells were subclustered on its latent space.

| Parameter | Value |
|---|---|
| Input | 61,467 cells, raw counts of all 21,591 genes |
| Labels | 35 annotated cell types |
| Batch | sample |
| n_layers | 2 |
| n_latent | 30 |
| max_epochs | 130 (scvi-tools default) |
| Neighbors (IM cells) | 15 |
| Leiden resolution | 0.6 (compared with 1.0) |

Of the nine subclusters, four were named Mo-IM (monocyte-derived IM, 1,066 cells), two IM (335 cells) and one cDC2 (135 cells). Subcluster 7 (80 cells) was removed as low quality. It had the lowest gene counts, a higher mitochondrial fraction than the other subclusters, and 88% of its cells came from one animal. Subcluster 8 (64 cells) was removed as possible doublets. It had the highest gene and UMI counts and co-expressed monocyte and macrophage genes.

The neighbor graph depends on the software environment, so rerunning Leiden can give slightly different subclusters (94% of cells received the same label in a rerun). The labels therefore come from the original subcluster assignment in `IM/im_subclusters_original.csv`, and the rerun is kept as a check. The trained model is not included. With `LOAD_MODEL = False` the notebook retrains it, which changes the latent space slightly but not the labels.

## 8. Cross-species comparison

### Human (HLCA, scArches)

Following the HLCA mapping tutorial with the core scANVI reference model from scvi-hub.

1. Convert macaque genes to human Ensembl IDs using one-to-one orthologs (Ensembl 113).
2. Map the query as a single unlabeled batch by scArches surgery.
3. Transfer labels at all five HLCA annotation levels.

| Parameter | Value |
|---|---|
| freeze_dropout | True |
| max_epochs | 500 |
| Early stopping | training ELBO, patience = 10 |
| kNN neighbors (label transfer) | 50 |
| Uncertainty cutoff | > 0.2 labeled "Unknown" |

```python
query = scvi.model.SCANVI.load_query_data(adata_query, hlca_ref_model,
                                          freeze_dropout=True)
query.train(max_epochs=500, early_stopping=True,
            early_stopping_monitor="elbo_train", early_stopping_patience=10,
            plan_kwargs={"weight_decay": 0.0})

knn = weighted_knn_trainer(train_adata=ref_emb, train_adata_emb="X",
                           n_neighbors=50)
labels, uncert = weighted_knn_transfer(query_adata=query_emb,
                                       query_adata_emb="X",
                                       label_keys="Level",
                                       knn_model=knn,
                                       ref_adata_obs=ref_emb.obs)
```

The mapping was run on 61,612 cells from an earlier version of the annotation (`HLCA/final_cells_for_hlca.csv`), and its labels are used for the 61,323 final cells. `HLCA/hlca_label_transfer.csv` has both the unfiltered labels and the labels after the uncertainty cutoff. The figures show the unfiltered labels with the median uncertainty of each cell type.

### Rhesus (RIRA, Seurat)

Pigtail macaque immune cells were compared with the Rhesus Immune Reference Atlas (GSE277821).

- The query is the 51,364 immune cells of the atlas, matched to RIRA by gene symbol (15,888 shared genes).
- The reference (about 19,000 cells) uses lung cells from the three RIRA objects. Cell types with fewer than 150 lung cells were topped up with up to 300 cells from other tissues. Non-lung macrophages (Mac. CD169+, Mac. Liver), promyelocytes and the CD45-High category were removed.
- The cell lists in `RIRA/` record the exact reference and query, so the transfer can be rerun without repeating the random sampling.
- A sensitivity analysis capped each RIRA cell type at 400 cells (`RIRA/rira_label_transfer_balanced_sensitivity.csv`).

```r
rira_ref <- NormalizeData(rira_ref)
rira_ref <- FindVariableFeatures(rira_ref, nfeatures = 2000)
rira_ref <- ScaleData(rira_ref)
rira_ref <- RunPCA(rira_ref, npcs = 30)
immune   <- NormalizeData(immune)

anchors <- FindTransferAnchors(reference = rira_ref, query = immune,
                               dims = 1:30, reference.reduction = "pca")
preds <- list()
for (lv in c("rira_fine", "rira_lineage", "rira_effector")) {
  refdata <- as.character(rira_ref[[lv, drop = TRUE]])
  refdata[is.na(refdata)] <- "not T/NK"
  names(refdata) <- colnames(rira_ref)
  preds[[lv]] <- TransferData(anchorset = anchors, refdata = refdata, dims = 1:30)
}
```

## 9. Figures

`Lung_Atlas_Figures.Rmd` loads `lung_final.RDS` and writes every figure panel to `FINAL_FIGURES/`. The main UMAP is the integration UMAP from section 5. The compartment UMAPs are recomputed from the scANVI latent space of the IM step (30 dimensions, 30 neighbors).

## Data availability

Raw and processed data are available from GEO under accession [GSE348023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348023).
