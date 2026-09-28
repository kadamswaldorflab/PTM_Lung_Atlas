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
  -> 6. Conversion to Seurat, marker-based annotation (HLCA-aligned)
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
  --fastqs=<fastq_dir> --sample=<sample>
```

## 2. Ambient RNA removal

CellBender `remove-background` on each raw feature-barcode matrix with default parameters.

```bash
cellbender remove-background --input raw_feature_bc_matrix.h5 \
  --output cellbender.h5 --fpr 0.01 --cuda
```

## 3. Quality control

Mitochondrial genes in the Mmul_10 annotation lack an `MT-` prefix, so the mitochondrial fraction is computed from the 13 mitochondrially encoded protein-coding genes.

| Filter | Threshold |
|---|---|
| Detected genes | 200 to 7,500 |
| UMIs | > 500 |
| Mitochondrial fraction | < 2.5% |
| Doublets | Removed with Scrublet (per sample) |

```python
mt_genes = ["ND1","ND2","ND3","ND4","ND4L","ND5","ND6",
            "ATP6","ATP8","COX1","COX2","COX3","CYTB"]
adata.var["mt"] = adata.var_names.isin(mt_genes)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

adata = adata[(adata.obs.n_genes_by_counts >= 200) &
              (adata.obs.n_genes_by_counts <= 7500) &
              (adata.obs.total_counts > 500) &
              (adata.obs.pct_counts_mt < 2.5)].copy()

sc.pp.scrublet(adata)  # or scrublet v0.2.3 directly
adata = adata[~adata.obs.predicted_doublet].copy()
```

## 4. Integration

**Highly variable genes** (raw counts):

```python
sc.pp.highly_variable_genes(adata, flavor="seurat_v3",
                            n_top_genes=4000, batch_key="sample")
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
vae.train(max_epochs=116)
adata.obsm["X_scVI"] = vae.get_latent_representation()
```

**scANVI** seeded with preliminary Leiden clusters (resolution 0.5, 30-NN graph on scVI latent):

```python
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=30)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_scvi")

lvae = scvi.model.SCANVI.from_scvi_model(vae, labels_key="leiden_scvi",
                                         unlabeled_category="Unknown")
lvae.train(max_epochs=10)
adata.obsm["X_scANVI"] = lvae.get_latent_representation()
```

## 5. Embedding and clustering

```python
sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, resolution=0.5)
```

**Integration metric:** inverse Simpson index of sample identity among each cell's 30 nearest neighbors. Median increased from 1.40 (unintegrated) to 2.27 (integrated).

## 6. Annotation (R/Seurat)

The AnnData object is converted to Seurat with `anndata` and `reticulate`.

```r
markers <- FindAllMarkers(obj, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25)

# Subclustering of heterogeneous clusters
sub <- FindNeighbors(sub, reduction = "scvi", dims = 1:50, k.param = 30)
sub <- FindClusters(sub, resolution = seq(0.02, 1.0, by = 0.02))
```

Clusters were annotated iteratively using cluster markers and canonical lineage markers, with nomenclature aligned to the Human Lung Cell Atlas (HLCA).

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
| Uncertainty cutoff | > 0.2 labeled "Unassigned" |

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
anchors <- FindTransferAnchors(reference = rira, query = immune, dims = 1:30)
pred <- TransferData(anchorset = anchors, refdata = rira$cell_type, dims = 1:30)
```

## Data availability

Raw and processed data: [GEO accession]
