# Macaque Lung Single-Cell Atlas: Analysis Pipeline

Code and parameters for processing, integrating, annotating, and comparing the macaque lung scRNA-seq atlas. This repository accompanies the Methods section of [manuscript citation].

## Software versions

| Tool | Version |
|---|---|
| Cell Ranger | 9.0.1 (reference built with mkref 7.2.0) |
| CellBender | 0.3.0 |
| Python | 3.10.19 |
| Scanpy | 1.11.5 |
| Scrublet | 0.2.3 |
| scvi-tools | 1.3.3 |
| scArches | 0.6.1 |
| R | 4.4.0 |
| Seurat | 5.5.0 |

A random seed of 42 is used throughout.

## Pipeline overview

```
FASTQ
  -> 1. Cell Ranger alignment (Mmul10, Ensembl 113 protein-coding)
  -> 2. CellBender ambient RNA removal (per sample)
  -> 3. QC filtering + Scrublet doublet removal (per sample)
  -> 4. HVG selection -> scVI integration -> scANVI refinement
  -> 5. UMAP + Leiden clustering
  -> 6. Conversion to Seurat, marker-based annotation and subclustering
  -> 7. Cross-species mapping (HLCA via scArches; RIRA via Seurat label transfer)
```

## 1. Reference and alignment

- Genome: Mmul10 (GenBank GCA_003339765.3)
- Annotation: Ensembl release 113, filtered to protein-coding genes
- Reference built with `cellranger mkref` (v7.2.0); reads aligned with `cellranger count` (v9.0.1)

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

## 4. Integration

**Highly variable genes** (raw counts):

```python
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000,
                            batch_key="sample", layer="counts")
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
```

**scVI** on raw counts of HVGs, sample as batch:

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

**scANVI** seeded with preliminary Leiden clusters (resolution 0.5, 30-NN graph on scVI latent):

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

**Integration metric:** inverse Simpson index of sample identity among each cell's 30 nearest neighbors, computed in R on the first 50 principal components of the unintegrated data (4,000 variable genes) and on the scANVI latent space. Median increased from 1.40 (unintegrated) to 2.27 (integrated).

## 6. Annotation (R/Seurat)

The AnnData object is converted to Seurat with `anndata` and `reticulate`.

```r
markers <- FindAllMarkers(obj, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25)

# Subclustering of heterogeneous clusters
sub <- FindNeighbors(sub, reduction = "scvi", dims = 1:50, k.param = 30)
sub <- FindClusters(sub, resolution = c(0.02, 0.05, 0.1, 0.2))
```

Clusters were annotated iteratively using cluster markers and canonical lineage markers. Low-quality, ambient-dominated and doublet clusters were removed (7,397 cells), yielding 61,467 cells across 35 cell types. 
Final labels and removed cells are provided in `final_cell_labels.csv` and `removed_cells.csv`.

## 7. Cross-species comparison

### Human (HLCA, scArches)

Following the HLCA mapping tutorial with the core scANVI reference model from scvi-hub.

1. Convert macaque genes to human Ensembl IDs using one-to-one orthologs (Ensembl 113).
2. Map query as a single unlabeled batch by scArches surgery.
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

### Rhesus (RIRA, Seurat)

Macaque immune cells compared with the Rhesus Immune Reference Atlas (GSE277821).

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

## Data availability

Raw and processed data: GSE348023[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348023]
