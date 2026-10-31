# Bioinformatic-Scripts

A collection of bioinformatics scripts, mostly Google Colab notebooks, for working with *Xenopus* data from [Xenbase](https://www.xenbase.org/). Most of them take RNA-Seq data that Xenbase has reprocessed from GEO and turn it into gene-level tables, expression timecourse plots or heatmaps. One notebook runs Gene Ontology (GO) enrichment on *Xenopus* gene lists.

## Contents

| File | What it does | Species / data |
|------|--------------|----------------|
| [`Parsing_RNA_SEQ_data.ipynb`](#parsing_rna_seq_dataipynb) | Converts any Xenbase-processed GEO RNA-Seq TPM matrix into a readable table | Any GSE, *X. laevis* or *X. tropicalis* |
| [`Heatmap_RNA_SEQ_data.ipynb`](#heatmap_rna_seq_dataipynb) | Builds that table for any GSE and draws heatmaps of selected genes | Any GSE, *X. laevis* or *X. tropicalis* |
| [`Tissue_heatmaps.ipynb`](#tissue_heatmapsipynb) | Heatmaps of selected genes across adult tissues | GSE73419 (*X. laevis*), GSE43520 (*X. tropicalis*) |
| [`Owens_Expression_data_graph.ipynb`](#owens_expression_data_graphipynb) | Plots a high-resolution developmental timecourse for selected genes | GSE65785 (*X. tropicalis*) |
| [`Session_Expression_data_graph.ipynb`](#session_expression_data_graphipynb) | Plots an oocyte-to-tadpole timecourse for selected genes | GSE73430 (*X. laevis*) |
| [`Isoform_expression_chart.ipynb`](#isoform_expression_chartipynb) | Plots the same timecourse per transcript isoform, labelled with RefSeq accessions | GSE65785 (*X. tropicalis*) |
| [`Xenopus_GO_enrichment.ipynb`](#xenopus_go_enrichmentipynb) | GO enrichment of a gene list against Xenbase GO annotations | *X. tropicalis* (default) or *X. laevis* |
| [`streamlit_app.py`](#streamlit_apppy) | Plain-Python version of the isoform timecourse notebook | GSE65785 (*X. tropicalis*) |
| `requirements.txt` | Pinned Python dependencies for running locally | — |

## Quick start

Each notebook has an **Open in Colab** badge at the top. Open it, change the values in the form fields (cells marked `# @param`), then choose **Runtime → Run all**. Plots appear inline. Files the notebook writes (for example `output.txt` or `<GSE>_TPM.tsv`) can be downloaded from the Colab file browser.

To run locally instead:

```bash
pip install -r requirements.txt
```

`requirements.txt` covers pandas, NumPy, SciPy, matplotlib, seaborn and requests. The heatmap notebooks also need `scikit-learn`, and the GO notebook installs `goenrich` and Graphviz itself (its install commands use `apt-get`, so they expect a Colab or Debian/Ubuntu environment).

## How the RNA-Seq notebooks work

Every notebook except the GO one follows the same pipeline.

1. **Download the TPM matrix.** Xenbase reprocesses public GEO/SRA RNA-Seq through its [RNA-Seq pipeline](https://gitlab.com/Xenbase/bioinformatics/RNAseq-Pipeline) ([Fortriede et al., 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7145613/) describes the general process; the pipeline has since been updated). The results are published per GSE and genome build:
   `https://bigfrog.xenbase.org/xenbase/genomics/GEO/<GSE>/<BUILD>/RNA-Seq/ExpressionFiles/Genes_TPM_Matrix.txt`
   In this matrix the rows are gene model IDs and the columns are SRR run accessions.
2. **Swap gene model IDs for gene symbols.** The notebooks use Xenbase's `models_gene_info` report. A gene model with no Xenbase gene ID keeps its model ID.
3. **Swap SRR columns for sample names.** SRR runs are mapped to GSM samples using `geo_srr_metadata_chd.txt`, and GSM samples to Xenbase track names using the dataset's `gsm_to_track.txt`.
4. **Average replicates.** Runs that share a track name are averaged into one column.
5. **Plot or export.** Depending on the notebook, the result is written to a TSV, plotted as a timecourse, or drawn as a heatmap.

Builds currently supported are `XENTR_10.0` (*X. tropicalis* v10) and `XENLA_10.1` (*X. laevis* v10.1). Not every GSE is processed against every build. Check the [Xenbase GEO download site](https://bigfrog.xenbase.org/pub/xenbase/genomics/GEO/) to see which builds exist for a given GSE.

### Gene naming

- *X. tropicalis* symbols have no suffix: `shh`, `wnt7a`, `tbxt`.
- *X. laevis* symbols carry a homeolog suffix: `shh.L`, `shh.S`. `Tissue_heatmaps.ipynb` adds `.L` and `.S` for you automatically. In the timecourse notebooks you must give the suffix yourself.
- Gene models with no assigned symbol stay as IDs such as `LOC108716104`.

---

## Notebooks

### `Parsing_RNA_SEQ_data.ipynb`
**Xenbase GEO RNA-Seq matrix mapper.** Takes a GSE ID and genome build and produces a gene-by-sample TPM table with gene symbols for rows, sample names for columns and replicates averaged.

- **Inputs:** `GSE_ID` (e.g. `GSE81458`), `Genome_build` (`XENTR_10.0` or `XENLA_10.1`)
- **Output:** `<GSE_ID>_TPM.tsv`, which is downloaded automatically from Colab
- Retries the gene-info download up to 5 times with a long timeout, because that report can be slow.

### `Heatmap_RNA_SEQ_data.ipynb`
Runs the same mapping for any GSE and build, then draws heatmaps for a list of genes you choose.

- **Inputs:** `GSE_ID`, `Genome_build`, `selected_genes` (a Python list)
- **Outputs:** `<GSE_ID>_TPM.tsv` and two heatmaps, samples × genes:
  - **min–max normalised** per gene, which shows *where* each gene is highest;
  - **log2(TPM + 1)**, which shows absolute expression levels.
- Gene matching ignores case, and any genes not found are listed.
- Downloaded files are cached. If a server returns an HTML error page instead of data, the notebook detects it, deletes the file and stops with an error.

### `Tissue_heatmaps.ipynb`
**Xenbase tissue heatmaps.** Same heatmap output, preset to adult tissue panels:

| Species | GSE | Build |
|---------|-----|-------|
| *X. laevis* | [GSE73419](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73419) | `XENLA_10.1` |
| *X. tropicalis* | [GSE43520](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43520) | `XENTR_10.0` |

- **Inputs:** `Species`, `selected_genes` (e.g. `['PAX8','tbxt','shh']`)
- For *X. laevis*, both homeologs (`.L` and `.S`) are matched for each gene you enter.
- The two species' datasets cover somewhat different sets of tissues.

### `Owens_Expression_data_graph.ipynb`
Plots *X. tropicalis* developmental expression from the high-temporal-resolution dataset [GSE65785](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65785) (Owens et al., 2016). The data run from fertilisation (NF1, 0 hpf) to NF43 (66 hpf), sampled every 30 minutes early on and hourly later.

- **Input:** `genes_to_plot`, e.g. `["nadsyn1", "wnt7a", "shh"]`
- **Outputs:** `output.txt` (the mapped TPM table) and two line plots of TPM over development, each showing the data points and a cubic-spline smoothed curve:
  - **Plot Graph**, with samples evenly spaced along the x-axis;
  - **Graph with proper timing**, with the x-axis scaled to hours post-fertilisation.
- Only the Clutch A series is plotted. Clutch B and the ribo-depleted (`ClutchA rd`) samples are dropped.

### `Session_Expression_data_graph.ipynb`
Plots *X. laevis* expression from oocyte to tadpole using [GSE73430](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73430) (Session et al., 2016, the *X. laevis* genome paper). There are 14 stages: oocyte I–II, III–IV and V–VI, unfertilised egg, then NF8, 9, 10–10.5, 12, 15, 20, 25, 29/30, 35/36 and 40.

- **Input:** `genes_to_plot`, using *X. laevis* homeolog names, e.g. `["shh.L", "tbxt.L", "wnt7a.L"]`
- **Outputs:** `output.txt` and two smoothed line plots, one with evenly spaced stages and one on an approximate time axis.

### `Isoform_expression_chart.ipynb`
An isoform-level version of the Owens timecourse. It uses `Isoforms_TPM_Matrix.txt` for GSE65785 and the *X. tropicalis* v10 GFF3 to label each transcript model with its RefSeq accession, so you can compare how a gene's isoforms change over development.

- **Input:** `genes_to_plot`, e.g. `["ldb1"]`
- **Outputs:** `output.txt` (isoform × sample, rows labelled `gene<RefSeq>`) and the same two plots as the Owens notebook, with one line per isoform.
- It downloads the full XENTR_10.0 GFF3 (gzipped), so the first run takes longer.

### `Xenopus_GO_enrichment.ipynb`
Runs GO term enrichment on a list of *Xenopus* gene symbols using [goenrich](https://github.com/jdrudolph/goenrich).

- **Input:** `gene_set`, a list of gene symbols (the default is a list of mesoderm and endoderm genes)
- **Background:** GO annotations from Xenbase (`xenbase.gaf.gz`) for one species, using the `go-basic.obo` ontology. The species is set by `desired_taxon_ids`: `'8364'` is *X. tropicalis* (the default) and `'8355'` is *X. laevis*. Commented-out lines let you switch to human, yeast or NCBI `gene2go` annotations instead.
- **Outputs:** an enrichment table (`example2.html`, top rows) and a GO graph of the enriched terms (`test.dot`, rendered to `test.png` with Graphviz).

### `streamlit_app.py`
A plain-Python script version of `Isoform_expression_chart.ipynb`. Running it writes `output.txt` and shows the isoform timecourse plot for the genes in `genes_to_plot` (default `["tbxt.2"]`). Despite its name it does not import Streamlit yet, so run it with `python streamlit_app.py` rather than `streamlit run`.

---

## Data sources

| Resource | URL |
|----------|-----|
| Xenbase-processed GEO RNA-Seq | `https://bigfrog.xenbase.org/xenbase/genomics/GEO/` |
| Gene model → symbol | `https://tools.xenbase.org/cgi-bin/reports/models_gene_info.cgi` |
| GSM ↔ SRR metadata | `https://download.xenbase.org/xenbase/GenePageReports/geo_srr_metadata_chd.txt` |
| *X. tropicalis* v10 GFF3 | `https://download.xenbase.org/xenbase/Genomics/JGI/Xentr10.0/XENTR_10.0_Xenbase.gff3.gz` |
| *X. laevis* v10.1 GFF3 | `https://download.xenbase.org/xenbase/Genomics/JGI/Xenla10.1/latest/XENLA_10.1_Xenbase.gff3.gz` |
| Xenbase GO annotations | `https://download.xenbase.org/xenbase/GenePageReports/xenbase.gaf.gz` |
| GO ontology | `http://purl.obolibrary.org/obo/go/go-basic.obo` |

## Notes and limitations

- TPM values are averaged across replicates. No differential-expression statistics are calculated.
- The smoothed curves are cubic-spline interpolations for display only. They are clipped at 0 and are not a statistical model.
- Gene models with no Xenbase gene symbol keep their original IDs, so a gene you search for may be missing from the output.
- Every notebook downloads its data from Xenbase servers when it runs. If a download times out, re-run the cell.

## Citation

If you use these scripts, please cite [Xenbase](https://www.xenbase.org/xenbase/static-xenbase/citingMOD.jsp) and the original GEO dataset(s).
