#!/usr/bin/env bash
# Raw inputs of the pipeline, one block per source (accessions and citations in data_sources.tsv)
# files kept with the data release are marked as such
# Large downloads (du -sh --apparent-size): Handsaker 17G, Mesecar 9.6G, Zemke 4.6G, Sloan 2.3G; GTEx gene counts 168M
set -e
cd "$(dirname "$0")/.."

# get <path under data/raw> <url>
get() {  # a finished file is kept; a transfer goes to <file>.part and is renamed when complete (wget -nc fails with -O)
  mkdir -p "$(dirname "data/raw/$1")"
  [[ -s data/raw/$1 ]] || { wget -O "data/raw/$1.part" "$2" && mv "data/raw/$1.part" "data/raw/$1"; }
}
# geo <folder> <GSE or GSM accession> <file>: supplementary file of a GEO series or sample
geo() {
  local kind=series
  [[ $2 == GSM* ]] && kind=samples
  get "$1/$3" "https://ftp.ncbi.nlm.nih.gov/geo/$kind/${2:0:${#2}-3}nnn/$2/suppl/$3"
}

# Proteomes ----
# Takasugi 2024: supplementary workbook and GEO counts
get takasugi2024_pxd040722/supp/41467_2024_52845_MOESM5_ESM.xlsb "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-52845-x/MediaObjects/41467_2024_52845_MOESM5_ESM.xlsb"
geo takasugi2024_pxd040722/geo GSE225576 GSE225576_Brain.txt.gz
# Langfelder 2016 / CHDI: HDinHD workbooks, GEO counts and the MaxQuant proteinGroups table
get langfelder2016_pxd006302/hdinhd/KIQ_144062_Evotec_Proteome_Part2_WP10_reprocessed_Striatum_Results.xlsx.gz "http://repository.hdinhd.org/data/hdinhd/proteomics/KIQ_144062_Evotec_Proteome_Part2_WP10_reprocessed_Striatum_Results.xlsx.gz"
get langfelder2016_pxd006302/hdinhd/KIQ-Brain-Samples-160725-Animal-to-Sample-mapping.xlsx.gz "http://repository.hdinhd.org/data/hdinhd/proteomics/KIQ-Brain-Samples-160725-Animal-to-Sample-mapping.xlsx.gz"
get langfelder2016_pxd006302/hdinhd/Allelic_Series_Decoder_Ring-1.5.xlsx.gz "http://repository.hdinhd.org/data/hdinhd/decoder-ring/Allelic_Series_Decoder_Ring-1.5.xlsx.gz"
geo langfelder2016_pxd006302/geo GSE65774 GSE65774_Striatum_mRNA_counts_processedData.txt.gz
# deposited with the data release: data/raw/langfelder2016_pxd006302/maxquant/striatum_PXD006302/proteinGroups.txt.gz (member proteinGroups.txt of maxQuantData.zip)
# Luo 2024: supplementary workbook
get luo2024_syn32136022/supp/41380_2024_2576_MOESM2_ESM.xlsx "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41380-024-02576-8/MediaObjects/41380_2024_2576_MOESM2_ESM.xlsx"
# Koopmans 2026: DIA-NN report and sample sheet
# deposited with the data release: data/raw/koopmans2026_pxd058441/diann/20220120_rklaas_SUN-schizo-cortex-layer1-3+4-6_MBR.pg_matrix.tsv (member 20220120_rklaas_SUN-schizo-cortex-layer1-3+4-6_MBR.pg_matrix.tsv.zst of diann.zip)
get koopmans2026_pxd058441/supp/MOESM3_supplementary_data_1_samples.xlsx "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-68950-y/MediaObjects/41467_2026_68950_MOESM3_ESM.xlsx"
# Eggers 2025: normalised proteins and sample sheet
get eggers2025_pxd051145/pride/norm_proteins.xlsx "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/05/PXD051145/norm_proteins.xlsx"
get eggers2025_pxd051145/pride/sample_information.xlsx "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/05/PXD051145/sample_information.xlsx"
# Anderton 2026: Spectronaut table and metadata (GitHub, pinned commit)
get anderton2026_biorxiv714902/github/data_raw/EA6_lysate_SN17_2pep_Analysis_V02.xlsx "https://raw.githubusercontent.com/Eanderton91/human_hippocampus_insoluble_proteome/b0c19de372b2e220ba0dfd5bd698c64e46eb2e90/data_raw/EA6_lysate_SN17_2pep_Analysis_V02.xlsx"
get anderton2026_biorxiv714902/github/data_raw/final_metadata.csv "https://raw.githubusercontent.com/Eanderton91/human_hippocampus_insoluble_proteome/b0c19de372b2e220ba0dfd5bd698c64e46eb2e90/data_raw/final_metadata.csv"
# Lu 2025: PulseDIA matrix and sample sheet (Google Drive)
get lu2025_pxd053154/gdrive/Aging_mouse_organ_PulseDIA_Matrix.tsv.xz "https://drive.usercontent.google.com/download?export=download&confirm=t&id=1E22Ye-Ov0jnbolS9qphiQR2vyISP_9FC"
get lu2025_pxd053154/gdrive/Aging_mouse_organ_SampleInformation.xlsx "https://drive.usercontent.google.com/download?export=download&confirm=t&id=1LQCiad0cMsSeCpCpow1WjcDcmxMlIBE-"
# Tsumagari 2023: supplementary workbooks
get tsumagari2023_pxd041485/supp/41598_2023_45570_MOESM3_ESM.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-023-45570-w/MediaObjects/41598_2023_45570_MOESM3_ESM.xlsx"
get tsumagari2023_pxd041485/supp/41598_2023_45570_MOESM4_ESM.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-023-45570-w/MediaObjects/41598_2023_45570_MOESM4_ESM.xlsx"
# Wang 2025: DIA-NN protein groups
get wang2025_pxd047296/pride/BrainAging_report.pg_matrix.tsv "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD047296/BrainAging_report.pg_matrix.tsv"
# Keele 2025: figshare tables S1 and S3
get keele2025_pxd058160/figshare/table_s1.xlsx "https://ndownloader.figshare.com/files/51129356"
get keele2025_pxd058160/figshare/table_s3.xlsx "https://ndownloader.figshare.com/files/51129338"
# Marino 2025: ubiquitin-enriched TMT10 PSMs (MassIVE)
get marino2025_pxd047951/massive/200904_AM_Brain_Mouse_UbEnr_TMT10_fr1_TargetPeptideSpectrumMatch.txt "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?forceDownload=true&file=f.MSV000093690/quant/200904_AM_Brain_Mouse_UbEnr_TMT10_fr1_TargetPeptideSpectrumMatch.txt"
# Lee 2023: Spectronaut protein groups (PRIDE), source data and sample information of the paper
get lee2023_pxd043880/pride/JL_EXP3_300819_Upload.xls "https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/09/PXD043880/JL_EXP3_300819_Upload.xls"
get lee2023_pxd043880/supp/13024_2023_650_MOESM1_ESM.xlsx "https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13024-023-00650-3/MediaObjects/13024_2023_650_MOESM1_ESM.xlsx"
get lee2023_pxd043880/supp/13024_2023_650_MOESM2_ESM.xlsx "https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13024-023-00650-3/MediaObjects/13024_2023_650_MOESM2_ESM.xlsx"
# Teng 2023: PSM tables of the 12 TMTpro plexes, channel map and Supplementary Table 1
pride=https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/10/PXD040884
for run in QEHFX3:1 QEHFX1:2 QEHFX3:3 QEHFX2:4 QEHFX1:5 QEHFX2:6 QEHFX3:7 QEHFX4:8 QEHFX1:9 QEHFX2:10 QEHFX5:11 QEHFX1:12; do
  IFS=: read instrument plex <<< "$run"
  file=NIAAA_TMTpro_${instrument}_020521_P${plex}_Global_PSMs.txt
  get teng2023_pxd040884/pride/$file "$pride/$file"
done
get teng2023_pxd040884/pride/Brain_AUD_Sample_Reference_File.xlsx "$pride/Brain_AUD_Sample_Reference_File.xlsx"
get teng2023_pxd040884/supp/41398_2023_2605_MOESM3_ESM.xlsx "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41398-023-02605-0/MediaObjects/41398_2023_2605_MOESM3_ESM.xlsx"
# Merrihew 2023: batch-adjusted protein group tables and sample sheets of four brain regions (Panorama Public)
panorama="https://panoramaweb.org/_webdav/Panorama%20Public/2022/MacCoss%20-%20Human%20AD%20Clean%20Diagnosis%20DIA%20Data"
for region in SMTG:SMTG:B1-B5 Hippocampus:Hipp:B1-B3 IPL:IPL:B1-B5 Caudate:Caudate:B1-B4; do
  IFS=: read folder prefix batches <<< "$region"
  get merrihew2023_pxd034525/panorama/$folder/level3b/$prefix-protGrp.batchadj.csv "$panorama/$folder/%40files/Level%203B/$prefix-protGrp.batchadj.csv"
  get merrihew2023_pxd034525/panorama/$folder/metadata/HuAD_Clean_Diagnosis_Samples-$prefix-$batches-Meta.csv "$panorama/$folder/%40files/MetaData/HuAD%20Clean%20Diagnosis%20Samples-$prefix-$batches-Meta.csv"
done
# Blumenreich 2024: PRIDE donor table, Supplementary Data 3 and the MaxQuant protein groups of four searches
get blumenreich2024_pxd047134/pride/data_tissue_key.xlsx "https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/02/PXD047134/data_tissue_key.xlsx"
get blumenreich2024_pxd047134/supp/41531_2024_645_MOESM4_ESM.xlsx "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41531-024-00645-x/MediaObjects/41531_2024_645_MOESM4_ESM.xlsx"
# deposited with the data release: data/raw/blumenreich2024_pxd047134/maxquant/{9619_STR,10401_SN,10017_CG,9105_MTG}/proteinGroups.txt (member combined/txt/proteinGroups.txt of combined_{9619,10401,10017,9105}.zip)
# deposited with the data release: data/raw/blumenreich2024_pxd047134/manifests/raw_file_map.tsv (run set and donor of each raw file, parsed from the PRIDE file names)
# Ping 2020: Proteome Discoverer protein table and sample traits (PRIDE)
get ping2020_pxd020296/pride/totalonly_multiconsensus__Proteins.txt "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020296/totalonly_multiconsensus__Proteins.txt"
get ping2020_pxd020296/pride/sampletraits.xlsx "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020296/sampletraits.xlsx"
# Umoh 2018 (ALS frontal cortex): Datasets EV1 (samples) and EV2 (log2 LFQ) of the EMBO Mol Med supplement
# deposited with the data release: data/raw/umoh2018_emmm201708202/emmm201708202-sup-000{2,3}-DatasetEV{1,2}.xlsb (supplement of doi 10.15252/emmm.201708202)
# Bulk RNA-seq ----
# BrainSeq III caudate (Benjamin 2022): gene counts
get benjamin2022_brainseq3_caudate/stage2/caudate_brainseq_phase3_hg38_rseGene_merged_n464.rda "https://caudate-eqtl.s3.us-west-2.amazonaws.com/caudate_brainseq_phase3_hg38_rseGene_merged_n464.rda"
# BrainSeq II (Collado-Torres 2019): unfiltered gene counts
get colladotorres2019_phs000979/libd/rse_gene_unfiltered.Rdata "https://s3.us-east-2.amazonaws.com/libd-brainseq2/rse_gene_unfiltered.Rdata"
# GTEx v10 brain gene counts and the open-access annotations
for tissue in amygdala anterior_cingulate_cortex_ba24 caudate_basal_ganglia cerebellar_hemisphere cerebellum cortex frontal_cortex_ba9 hippocampus hypothalamus nucleus_accumbens_basal_ganglia putamen_basal_ganglia spinal_cord_cervical_c-1 substantia_nigra; do
  get gtexv10_phs000424/gene_reads/gene_reads_v10_brain_$tissue.gct.gz "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/counts-by-tissue/gene_reads_v10_brain_$tissue.gct.gz"
done
get gtexv10_phs000424/annotations/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
get gtexv10_phs000424/annotations/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
# Markunas 2021 nucleus accumbens: Salmon counts, TPM and series matrix
geo markunas2020_gse171936/stage2 GSE171936 GSE171936_nac_rnaseq_salmon_dge_gene_data_counts_filters.csv.gz
geo markunas2020_gse171936/stage2 GSE171936 GSE171936_nac_rnaseq_salmon_dge_gene_data_tpm_filters.csv.gz
get markunas2020_gse171936/stage2/GSE171936_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171936/matrix/GSE171936_series_matrix.txt.gz"
# Matlik 2024 sorted nuclei: one count file per sample (GSM7106325-GSM7106548) and the series matrix
for n in $(seq 7106325 7106548); do
  wget -nc -q -r -np -nd -e robots=off -A 'GSM*_Gene.Counts.txt.gz' -P data/raw/matlik2024_gse227729/stage2 "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7106nnn/GSM$n/suppl/"
done
get matlik2024_gse227729/stage2/GSE227729_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227729/matrix/GSE227729_series_matrix.txt.gz"
# Hahn 2023 bulk RNA-seq of 15 brain regions: count table and series matrix (GEO)
geo hahn2023_gse212336/geo GSE212336 GSE212336_BulkSeq_Counttable.txt.gz
get hahn2023_gse212336/geo/GSE212336_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212336/matrix/GSE212336_series_matrix.txt.gz"
# Shavlakadze 2026 mouse multi-tissue RNA-seq: raw counts and series matrix (GEO)
geo shavlakadze2026_gse314817/geo GSE314817 GSE314817_mouse_raw_counts_wo_muscle.txt.gz
get shavlakadze2026_gse314817/geo/GSE314817_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE314nnn/GSE314817/matrix/GSE314817_series_matrix.txt.gz"
# Ramasamy 2014 UKBEC exon arrays: GPL5175 transcript-cluster series matrix (GEO)
get ramasamy2014_gse46706/geo/GSE46706-GPL5175_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46706/matrix/GSE46706-GPL5175_series_matrix.txt.gz"
# deposited with the data release: data/raw/ramasamy2014_gse46706/platform/GPL5175_data.txt.gz (GEO text view of the GPL5175 platform table, gzipped; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL5175&targ=self&form=text&view=data)
# Single-nucleus data ----
# Sloan 2026 prefrontal cortex: raw-count h5ad (Zenodo) and the donor sheet (Google Sheet)
get sloan2026_snt456vhrv853/zenodo/adata_snrna_rawcounts_0807_filtered.h5ad "https://zenodo.org/records/17467925/files/adata_snrna_rawcounts_0807_filtered.h5ad?download=1"
get sloan2026_snt456vhrv853/supp/zenodo_donor_metadata_sheet.csv "https://docs.google.com/spreadsheets/d/1rmVE_ZtojFaL0ChXuxNjC9N2wWLPpffSWVh0U2asZP4/export?format=csv"
# Zemke 2026 hippocampus: one 10x h5 per sample and the labelled nuclei table
for id in hc1134 hc1153 hc11 hc1203 hc1216 hc1265 hc1271 hc12 hc13344 hc13394 hc13414 hc1745 hc19 hc212191 hc26 hc29 hc35 hc40 hc46426 hc4781 hc5021 hc5087 hc5551 hc5579 hc5610 hc5614 hc6021 hc6052 hc69984 hc73787 hc73 hc76 hc77 hc78 hc81 hc8 hc935 hc937 hc98 hc9; do
  geo zemke2026_gse278576/geo GSE278576 GSE278576_${id}_raw_feature_bc_matrix.h5
done
geo zemke2026_gse278576/geo GSE278576 GSE278576_hippocampus_RNA_seurat_object_filtered_cells_metadata.tsv.gz
# Handsaker 2025 caudate (NeMO): donor-level h5 per village reaction, nuclei labels and donor sheet
nemo=https://data.nemoarchive.org/other/grant/broad_mccarroll_HD_2024/mccarroll/transcriptome/sncell/10x_v3.1/human/processed
for village in 1:0.550:12345678 14:0.550:12345678 20:0.530:12345678 4:0.600:1234578 5:0.470:12345678 8:0.600:12345678; do
  IFS=: read number cut reactions <<< "$village"
  for ((i = 0; i < ${#reactions}; i++)); do
    file=CHDI_Village${number}_rxn${reactions:i:1}.umi_500-Inf_intronic_${cut}-1.000.10X.exonic+intronic.donors.h5
    get handsaker2025_nemo_broadhdcaudate/stage2/$file "$nemo/counts/$file"
  done
done
get handsaker2025_nemo_broadhdcaudate/stage2/10x_village_cell_metadata_with_SPN_subtypes.txt "https://data.nemoarchive.org/other/grant/broad_mccarroll_HD_2024/mccarroll/transcriptome/sncell/10x_v3.1/human/processed/other/10x_village_cell_metadata_with_SPN_subtypes.txt"
get handsaker2025_nemo_broadhdcaudate/stage2/donor_metadata.txt "https://data.nemoarchive.org/other/grant/broad_mccarroll_HD_2024/mccarroll/transcriptome/sncell/10x_v3.1/human/processed/other/donor_metadata.txt"
# Mesecar 2026 putamen and other regions: raw-count h5ad (Zenodo, 9.6G)
get mesecar2026_zenodo7803697/stage2/aging.pegasus.leiden_085.subclustered.h5ad "https://zenodo.org/api/records/17784771/files/aging.pegasus.leiden_085.subclustered.h5ad/content"
# Hahn 2023 (GSE212576), Green 2025 (GSE277313) and the Allen mouse ageing atlas (Jin 2025; adult and aged pseudobulks, script 06 uses the adult animals): pseudobulk tables kept with the data release
# deposited with the data release: data/raw/{hahn2023_gse212576,green2025_gse277313,jin2025_zengagingmouse10xv3}/pseudobulk/{matrix.parquet,samples.tsv,genes.tsv}
# Modification site maps ----
# Lu 2024: m5C-TAC-seq sites (supplementary table S4, hg19)
get lu2024_gse242724/1-s2.0-S1097276524005288-mmc5.xlsx "https://ars.els-cdn.com/content/image/1-s2.0-S1097276524005288-mmc5.xlsx"
# Dai 2024: UBS-seq m5C sites of HEK293T and HeLa (GEO)
geo dai2024_gse225614 GSE225614 GSE225614_HEK293T-WT_sites.tsv.gz
geo dai2024_gse225614 GSE225614 GSE225614_HeLa-WT_sites.tsv.gz
# Huang 2019: RNA-BisSeq m5C sites of human tissues (GEO, hg19)
geo modmaps/GSE122254 GSE122254 GSE122254_human_sites.xls.gz
# Zhang 2019: base-resolution m7G sites curated by Sci-ModoM
get zhang2019_m7g_scimodom/hela_m7G.bedrmod "https://scimodom.dieterichlab.org/api/v0/transfer/dataset/JDLWzBfFGvGJ"
get zhang2019_m7g_scimodom/hepg2_m7G.bedrmod "https://scimodom.dieterichlab.org/api/v0/transfer/dataset/GFGhkyJXxJSx"
# Zhang 2019 and Zhao 2023: m7G-MeRIP peaks and METTL1 knock-down differences (GEO)
for peaks in HeLa_HighFC HepG2_HighFC HEK293T_HighFC HeLa_OverallFC HepG2_OverallFC HeLa_siMETTL1_ConsistentDiff_peak HepG2_shMETTL1_ConsistentDiff_peak; do
  geo zhang2019_gse112276 GSE112276 GSE112276_$peaks.txt.gz
done
geo zhao2023_gse193039 GSM5766794 GSM5766794_r1.mrna.peak.txt.gz
geo zhao2023_gse193039 GSM5766795 GSM5766795_r2.mrna.peak.txt.gz
geo zhao2023_gse193039 GSM5766798 GSM5766798_m7g.1.mrna.peak.txt.gz
geo zhao2023_gse193039 GSM5766799 GSM5766799_m7g.2.mrna.peak.txt.gz
# Pseudouridine maps: BID-seq, BACS, PRAISE, McCormick nanopore calls and PRAISE-seq
get dai2023_bidseq/41587_2022_1505_MOESM3_ESM.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-022-01505-w/MediaObjects/41587_2022_1505_MOESM3_ESM.xlsx"
get xu2024_bacs/41592_2024_2439_MOESM3_ESM.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-024-02439-8/MediaObjects/41592_2024_2439_MOESM3_ESM.xlsx"
get zhang2023_praise/41589_2023_1304_MOESM3_ESM.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs41589-023-01304-7/MediaObjects/41589_2023_1304_MOESM3_ESM.xlsx"
get mccormick2026_nar/gkag353_supplemental_files.zip "https://pmc-oa-opendata.s3.amazonaws.com/PMC13096809.1/gkag353_supplemental_files.zip"
geo praiseseq2026_gse326442 GSE326442 GSE326442_PRAISE-seq-HEK293T.xlsx
# Sun 2025: GLORI 3.0 of mouse hippocampus (GEO)
geo modmaps/GSE270643 GSE270643 GSE270643_RAW.tar
# Dai 2023: BID-seq of 12 mouse tissues (GEO)
geo modmaps/GSE179798 GSE179798 GSE179798_Mouse_mRNA_12-Tissues_BID-seq.xlsx
# Zhang 2024: BID-seq of mouse embryonic stem cells (GEO, GRCm39 coordinates)
geo modmaps/GSE238245 GSE238245 GSE238245_Pseudouridine_mESC-WT_polyA-RNA.xlsx
# He 2025: Uli-epic BID-seq of mouse neural stem cells (GEO)
geo modmaps/GSE277358 GSM8520916 GSM8520916_WT-NSC-Uli-epic-BIDSeq-rep1.xls.gz
geo modmaps/GSE277358 GSM8520917 GSM8520917_WT-NSC-Uli-epic-BIDSeq-rep2.xls.gz
# Yang 2017: RNA-BisSeq m5C sites of mouse tissues (GEO)
geo modmaps/GSE93749 GSE93749 GSE93749_Mouse_m5C_sites_information.xls.gz
# Amort 2017: RNA-BisSeq m5C sites of mouse brain (GEO)
geo modmaps/GSE83432 GSE83432 GSE83432_Brain_total.xlsx
# Wang 2024: m7GHub miCLIP regions of mouse brain and embryonic stem cells
get modmaps/m7ghub2/m7G-miCLIP-seq.csv "http://www.rnamd.org/m7GHub2/download/m7G-miCLIP-seq.csv"
# Shafik 2026: m6A-seq peaks of five human brain regions (supplementary table, hg38)
get modmaps/shafik2026/41593_2025_2112_MOESM2_ESM.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-025-02112-z/MediaObjects/41593_2025_2112_MOESM2_ESM.xlsx"
# Annotation ----
# UCSC liftOver chains hg19 to hg38 and mm39 to mm10
get ucsc_liftover/hg19ToHg38.over.chain.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
get ucsc_liftover/mm39ToMm10.over.chain.gz "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/liftOver/mm39ToMm10.over.chain.gz"
# GENCODE human v44 (basic GTF, GRCh38 genome)
get gencode_v44/gencode.v44.basic.annotation.gtf.gz "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz"
get gencode_v44/GRCh38.primary_assembly.genome.fa.gz "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"
# GENCODE mouse vM25 (GTF, GRCm38 genome) and the vM33 GTF for canonical tags
get gencode_vM25/gencode.vM25.basic.annotation.gtf.gz "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gtf.gz"
get gencode_vM25/GRCm38.primary_assembly.genome.fa.gz "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz"
get gencode_vM33/gencode.vM33.basic.annotation.gtf.gz "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.basic.annotation.gtf.gz"
# Kluever 2022: protein half-lives of mouse brain (PMC open-access tables)
get kluever2022_sciadv/sciadv.abn4437_tables_s1_to_s7.zip "https://pmc-oa-opendata.s3.amazonaws.com/PMC9122331.1/sciadv.abn4437_tables_s1_to_s7.zip"
# Agarwal and Kelley 2022: mRNA half-life PC1 of mouse and human (GitHub, pinned commit)
get saluki_halflife/all_HLs_mouse_PC1.txt "https://raw.githubusercontent.com/vagarwal87/saluki_paper/6aba76beca92f64a760f252db80d1363855a9bd2/Fig2_S3/all_HLs_mouse_PC1.txt"
get saluki_halflife/all_HLs_human_PC1.txt "https://raw.githubusercontent.com/vagarwal87/saluki_paper/6aba76beca92f64a760f252db80d1363855a9bd2/Fig2_S3/all_HLs_human_PC1.txt"
# Gene Ontology release: ontology and annotations of mouse and human
get go_release/go-basic.obo "https://release.geneontology.org/2026-08-05/ontology/go-basic.obo"
get go_release/mgi.gaf.gz "https://release.geneontology.org/2026-05-19/annotations/mgi.gaf.gz"
get go_release/goa_human.gaf.gz "https://release.geneontology.org/2026-05-19/annotations/goa_human.gaf.gz"
# HGNC (2026-09-26), MGI reports, UniProt gene tables and the Orthrus and RiboNN predictions: snapshots kept with the data release
# deposited with the data release: data/raw/{hgnc,mgi,uniprot,gpu_precomputed}/
# Reader CLIP peaks (data/raw/clip_standard): public eCLIP, iCLIP and irCLIP libraries (accessions in data_sources.tsv)
# processed with nf-core/clipseq 1.0.0 (PureCLIP) on the UCR HPCC cluster: pipeline/extra/select_clip_datasets.py writes
# selection.tsv, clip_peaks.sbatch runs one row per array task and returns the binding regions of each IP library;
# deposited with the data release
# Reader binding and perturbation compendium (script 37) ----
# ENCODE eCLIP IDR peaks (data/raw/encode2020_eclip_idr) and the clipflow peaks (data/raw/clipflow_peaks): deposited with the data release
# gene lists of the reader binding sources and of the Mettl14 knockout (ref/target_sets/*.tsv, stage2/shi2026_ythdf1_clip.tsv and stage2/koranda2018_mettl14.tsv: Darnell 2011, Maurin 2018, Wu 2019, Shi 2026, Ascano 2012, Liu 2015, Liu 2017, Wang 2014, Koranda 2018)
# and the Wang 2014 supplementary table (ref/target_sets/raw/wang2014): deposited with the data release, not tracked in git
# IGF2BP2 binding: Conway 2016 (H9 eCLIP), Degrauwe 2016 (MGG8 PAR-CLIP), Han 2022 (A172 eCLIP), Preitner 2016 (E14 brain CLIP; mm9 to mm10 chain), Huang 2018 (HEK293T RIP), GSE133242 (MEF RIP), GSE192705 (adipocyte RIP)
geo conway2016_gse78507 GSM2071743 GSM2071743_IGF2BP2_H9ES_Rep1_eCLIP.InputNormalizedPeaks.bed.gz
geo conway2016_gse78507 GSM2071746 GSM2071746_IGF2BP2_H9ES_Rep2_eCLIP.InputNormalizedPeaks.bed.gz
for file in PARCLIP_IMP2a_r1.bed PARCLIP_IMP2a_r2.bed PARCLIP_IMP2s_r1.bed PARCLIP_IMP2s_r2.bed RNA_IMP2s_KD.txt RNA_IMP2a_OE.txt RNA_brain_IMP2_KO.txt; do
  geo degrauwe2016_gse73845 GSE73845 GSE73845_$file.gz
done
geo han2022_gse192792 GSM5765415 GSM5765415_IMP2_IP_1.peaks.bed.annot.txt.gz
geo han2022_gse192792 GSM5765416 GSM5765416_IMP2_IP_2.peaks.bed.annot.txt.gz
geo preitner2016_gse83822 GSE83822 GSE83822_IMP2_peaks_RefSeq_mRNA.txt.gz
get preitner2016_gse83822/mm9ToMm10.over.chain.gz "https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz"
f=data/raw/huang2018_gse90639/GSE90639_IGF2BPs_RIP_bingding_gene_and_exp_info.xlsx  # GEO serves the workbook gzipped
[[ -s $f ]] || { geo huang2018_gse90639 GSE90639 ${f##*/}.gz; gunzip $f.gz; }
geo huang2018_gse90642 GSE90642 GSE90642_shM14_vs_shCont_diff_gene_and_diff_m6A.xlsx
geo huang2018_gse90684 GSE90684 GSE90684_expression_matrix_of_genes_rpkm_average_result.txt.gz
for sample in IgG.1:GSM3903891 IgG.2:GSM3903892 IgG.3:GSM3903893 IMP2.1:GSM3903897 IMP2.2:GSM3903898 IMP2.3:GSM3903899; do
  geo gse133242_imp2_mef_rip ${sample#*:} ${sample#*:}_IP.${sample%:*}.txt.gz
done
geo gse192705_igf2bp2_rip_adipocyte GSE192705 GSE192705_DESeq.input_noAb_IgG.Igfbp2.over10_either_over2-fold_enrichment.underp05_Adjpval.txt.gz
# IGF2BP2 perturbation: ENCODE shRNA RNA-seq (HepG2, K562), GSE146726, GSE292558, Regue 2019 (GSE119801)
get encode_regulator_rnaseq/ENCSR478FJK/ENCFF302MAE.tsv "https://www.encodeproject.org/files/ENCFF302MAE/@@download/ENCFF302MAE.tsv"
get encode_igf2bp2_kd_k562/ENCFF407VNW.tsv "https://www.encodeproject.org/files/ENCFF407VNW/@@download/ENCFF407VNW.tsv"
for cell in J82 UMUC3; do geo gse146726_igf2bp2_kd GSE146726 GSE146726_${cell}_DESeq2_IGF2BP2.txt.gz; done
geo gse292558_igf2bp2_ko_bmm GSE292558 GSE292558_DESeq2_IGF2BP2_KO_vs_WT.txt.gz
for sample in WT.total.1:GSM3384028 WT.total.2:GSM3384029 WT.total.3:GSM3384030 KO.Total.1:GSM3384031 KO.Total.2:GSM3384032 KO.Total.3:GSM3384033; do
  geo regue2019_gse119801 ${sample#*:} ${sample#*:}_${sample%:*}.count.txt.gz
done
# YTHDF binding: Shi 2018 (hippocampus CLIP), Niu 2022 (retina and cortex RIP), Flamand 2022 (neuron RIP), Yu 2021 (cerebellar granule RIP; Adv Sci supplement, Europe PMC), GSE266096 (NSC RIP), Zhang 2024 (HEK293T CLIP), Wang 2015 (HeLa PAR-CLIP)
# the YTHDF1 striatal CLIP (Shi 2026) and the Wang 2014 YTHDF2 targets are gene lists in ref/target_sets, deposited with the data release
geo shi2018_gse106607 GSE106607 GSE106607_A_Ythdf1_CLIP_seq_peak_list.xlsx
geo niu2022_gse199962 GSE199962 GSE199962_gene_count.txt.gz
for i in 1 2 3; do geo flamand2022_gse194207 GSE194207 GSE194207_DF${i}_RIP_DESeq.txt.gz; done
get yu2021_advsci_supp_pmc8596126/supp.zip "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8596126/supplementaryFiles"
[[ -s data/raw/yu2021_advsci_supp_pmc8596126/ADVS-8-2101329-s008.xlsx ]] || unzip -oq data/raw/yu2021_advsci_supp_pmc8596126/supp.zip -d data/raw/yu2021_advsci_supp_pmc8596126
geo gse266096_ythdf2_nsc GSE266096 GSE266096_Raw_count_table_RNA-seq.txt.gz
geo gse266096_ythdf2_nsc GSE266096 GSE266096_Raw_count_table_YTHDF2-RIP-seq.txt.gz
geo zhang2024_gse264498 GSE264498 GSE264498_YTHDF1_peaks_WT.bed.gz
geo wang2015_gse63591 GSE63591 GSE63591_A-Y1-PARCLIP.xlsx
# YTHDF perturbation: Yu 2021 (GSE153688), Zhuang 2023 (GSE171790), GSE306850, Shi 2018 (RNA-seq), Li 2018, GSE266096, Zou 2023 (neurons, HeLa), Lasman 2020 (mESC Ribo-seq and RNA), Wang 2015 and Shi 2017 (HeLa Ribo-seq), Heck 2020, GSE270417
geo yu2021_gse153688 GSE153688 GSE153688_Raw_gene_counts_matrix_KD.txt.gz
geo huang2021_gse171790 GSE171790 GSE171790_cKO-RNAseq_gene_count.xls.gz
geo gse306850_ythdf1 GSE306850 GSE306850_gene_count_matrix.xlsx
geo shi2018_gse106607 GSE106607 GSE106607_C_Hippocampal_RNA_seq_WT_KO_processed.xlsx
geo li2018_gse104867 GSE104867 GSE104867_RAW.tar
[[ -s data/raw/li2018_gse104867/GSM2809105_input_KO3.raw.counts.txt.gz ]] || tar xf data/raw/li2018_gse104867/GSE104867_RAW.tar -C data/raw/li2018_gse104867
geo zou2023_gse214882 GSE214882 GSE214882_posterior_mean_estimate.xlsx
geo zou2023_gse201070 GSE201070 GSE201070_HeLa_DF_Knockdown.csv.gz
geo lasman2020_gse147849 GSE147849 GSE147849_DESeq_mESCs.txt.gz
geo lasman2020_gse148039 GSE148039 GSE148039_RAW.tar
[[ -s data/raw/lasman2020_gse148039/GSM4473278_Ythdf321_KO_B.ribo.txt.gz ]] || tar xf data/raw/lasman2020_gse148039/GSE148039_RAW.tar -C data/raw/lasman2020_gse148039
geo wang2015_gse63591 GSE63591 GSE63591_C-Y1-Ribosome_profiling.xlsx
geo shi2017_gse86214 GSE86214 GSE86214_siYTHDF3_RiboProfile.xlsx
geo heck2020_gse133898 GSE133898 GSE133898_Control_YTH2KD_HTseq_raw_counts.txt.gz
geo gse270417_ythdf_gastric GSE270417 GSE270417_raw_counts.txt.gz
# Perturbation tables of the sign check (script 34) ----
# Fmr1 knockout: Shu 2020 (GSE140565, cortex RNA and ribosome footprints), Aryal 2021 (GSE143659, folder seo2022_gse143659), Jung 2023 (GSE207145, cortex, hippocampus, cerebellum), Suresh 2026 (GSE325852, RiboTag of Pvalb and Camk2a neurons)
geo shu2019_gse140565 GSE140565 GSE140565_RNADEseq_WTxFK.txt.gz
geo shu2019_gse140565 GSE140565 GSE140565_RPFDEseq_FKxWT.txt.gz
for kind in rna ribo te; do geo seo2022_gse143659 GSE143659 GSE143659_${kind}_fxsVwt_all.csv.gz; done
geo gse207145_fmr1_tissues GSE207145 GSE207145_RNAseq.DE.log2FC_Brain.tsv.gz
geo gse325852_fmr1_translatome GSE325852 GSE325852_ALL_Counts_Pvalb_Camk2.csv.gz
# Wu 2019 (GSE100491): brain RNA-seq counts of the Prrc2a conditional knockouts and the Fto transgene, one sample per group
for sample in GSM2684967_WT GSM2684968_cKO GSM2803005_prrc2a-WT GSM2803006_prrc2a-olig2-cre GSM2803007_Fto-WT GSM2803008_Fto-Tg; do
  geo wu2019_gse100491 ${sample%%_*} ${sample}-reads-count.txt.gz
done
# Engel 2018 (GSE113793): Mettl3 and Fto conditional knockouts of adult neurons; RNA_analysis.csv is the fourth member of the 2.9 GB processed-data archive, so only its first 8 MB are read
d=data/raw/gse113793_mettl3_fto_cortex
[[ -s $d/RNA_analysis.csv ]] || { mkdir -p $d; wget -qO- --header "Range: bytes=0-8388607" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113793/suppl/GSE113793_processed_data_files.tar.gz" | tar xzO --occurrence=1 RNA_analysis.csv > $d/RNA_analysis.csv.part && mv $d/RNA_analysis.csv.part $d/RNA_analysis.csv; }
# Du 2021 (GSE154992): Mettl3 and Fto conditional knockouts of the embryonic cortex, RNA-seq and Ribo-seq read counts
for kind in RNA-Seq Ribo-seq; do geo gse154992_mettl3_ribo GSE154992 GSE154992_${kind}_read_count.txt.gz; done
# Mettl14 knockout of the striatum (Koranda 2018): ref/target_sets/stage2/koranda2018_mettl14.tsv, deposited with the data release
# Curated RNA modification regulators: ref/gene_sets/rna_mod_regulators.tsv, tracked in the repository

# External validation proteomes (mouse brain ageing, control animals; prepared by 01) ----
# Keele 2023: Table S1 (supp_tables.zip), the peptide tables of striatum, cerebellum and hippocampus (data.zip, 399 MB) and the UniProt mouse ID-mapping table (release 2026_03, relnotes.txt saved with it)
d=keele2023_pxd034029
get $d/figshare/supp_tables.zip "https://ndownloader.figshare.com/files/41032439"
get $d/figshare/data.zip "https://ndownloader.figshare.com/files/39688207"
[[ -s data/raw/$d/figshare/supp_tables/table_s1.xlsx ]] || unzip -oq data/raw/$d/figshare/supp_tables.zip -d data/raw/$d/figshare
[[ -d data/raw/$d/figshare/data_zip/data/PeptideQuantitation ]] || unzip -oq data/raw/$d/figshare/data.zip "data/PeptideQuantitation/*_Striatum.tsv" "data/PeptideQuantitation/*_Cerebellum.tsv" "data/PeptideQuantitation/*_Hippocampus.tsv" -d data/raw/$d/figshare/data_zip
uniprot=https://ftp.uniprot.org/pub/databases/uniprot/current_release
get $d/uniprot/MOUSE_10090_idmapping_selected.tab.gz "$uniprot/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping_selected.tab.gz"
get $d/uniprot/relnotes.txt "$uniprot/relnotes.txt"
# Scifo 2025: the first 48.5 MB of the brain mzTab (metadata, protein and peptide tables; the PSM section follows) from the jPOST FTP
f=data/raw/scifo2025_pxd058684/jpost/Mar032021_AgeingMouse_Brain_LFQ_120min_60K_TopS.MTD-PRT-PEP.mzTab.gz
[[ -s $f ]] || { mkdir -p "$(dirname "$f")"; curl -sS --range 0-48499999 "ftp://ftp.jpostdb.org/JPST003472/Mar032021_AgeingMouse_Brain_LFQ_120min_60K_TopS.mzTab" | gzip > "$f.part" && mv "$f.part" "$f"; }
# Yu 2020: Dataset S5, member of the supplementary zip of Europe PMC
d=yu2020_pxd017385
get $d/supp/PMC7211924_supplementaryFiles.zip "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7211924/supplementaryFiles"
[[ -s data/raw/$d/supp/pnas.1919410117.sd05.xlsx ]] || unzip -oq data/raw/$d/supp/PMC7211924_supplementaryFiles.zip pnas.1919410117.sd05.xlsx -d data/raw/$d/supp
# Huang 2024: MaxQuant protein groups and the PRIDE record (sample-to-channel scheme)
pride=https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/12/PXD043352
get huang2024_pxd043352/pride/proteinGroups.txt "$pride/proteinGroups.txt"
get huang2024_pxd043352/pride/PXD043352.project.json "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD043352"
# Jha 2026: MaxQuant protein groups and Supplementary Table 2 (sample labels of the reporter channels)
get jha2026_pxd060735/pride/proteinGroups.txt "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/12/PXD060735/proteinGroups.txt"
get jha2026_pxd060735/supp/MOESM3_supplementary_table_2.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs44323-025-00066-w/MediaObjects/44323_2025_66_MOESM3_ESM.xlsx"
# Oliviero 2022: Supplementary file 2 (whole-brain total lysate)
get oliviero2022_elife/elife/elife-73524-supp2-v2.xlsx "https://cdn.elifesciences.org/articles/73524/elife-73524-supp2-v2.xlsx"
# Mohallem 2024: Supplemental Table S1 (MaxQuant protein groups)
get mohallem2024_msv000094441/supp/mmc2.xlsx "https://ars.els-cdn.com/content/image/1-s2.0-S1535947624001099-mmc2.xlsx"
# Marchionini 2025: run overview (age, condition, sample number)
get marchionini2025_pxd054495/pride/sample_overview.xlsx "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/04/PXD054495/sample_overview.xlsx"
# deposited with the data release: data/raw/marchionini2025_pxd054495/pride_stream/precursors.parquet (precursor rows of the 9.5 GB DIA-NN report.tsv of PXD054495 with run, protein group, precursor, quantity and q-values, streamed from PRIDE)
# Vignane 2026: published log2 intensities of the total proteome
get filipovic2026_pxd045785/pride/Results_Totalproteome.xlsx "https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/05/PXD045785/Results_Totalproteome.xlsx"
# deposited with the data release: data/raw/gebru2024_pxd049439/derived/{Hippocampus,Amygdala}_pg_maxlfq.tsv.gz (PG.MaxLFQ of the protein groups passing the DIA-NN pg_matrix filter per run, streamed from DIANN_LB_Niat_Hippo_M_F_results_Rerun.tsv and Amyg_WTvKO_MvF_rerun.tsv of PXD049439, 3.5 and 3.7 GB)
# Yarbro 2025: Supplementary Data 1-24
get yarbro2025_nat_commun/supp/MOESM3_supplementary_data_1-24.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-56853-3/MediaObjects/41467_2025_56853_MOESM3_ESM.xlsx"
# Chen 2025: Level 3B protein groups and run metadata (Panorama Public)
panorama="https://panoramaweb.org/_webdav/Panorama%20Public/2024/MacCoss%20-%20AD-BXD%20Mouse%20PreFrontal%20Cortex%20Proteomics/%40files/Supplementary%20Files"
get chen2025_pxd045403/panorama/level3b/MouAD_PFC_AXBDX_Protein-BatchAdj-v40125.csv "$panorama/Level%203B/MouAD_PFC_AXBDX_Protein-BatchAdj-v40125.csv"
get chen2025_pxd045403/panorama/metadata/MouAD-PFC-meta-SAGE-2.0.xlsx "$panorama/Metadata/MouAD-PFC-meta-SAGE-2.0.xlsx"

# Published proteomic ageing clocks: linear weights, benchmarks of the zero-shot clock (read by 06) ----
# Lehallier 2019: Supplementary Tables (ST1 SomaScan annotation, ST7 LASSO coefficients of the plasma clock)
get lehallier2019_clock/41591_2019_673_MOESM3_ESM.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-019-0673-2/MediaObjects/41591_2019_673_MOESM3_ESM.xlsx"
# Oh 2023: Supplementary Tables (ST3 SomaScan annotation, ST7 coefficients of the bagged LASSO organ models)
get oh2023_clock/41586_2023_6802_MOESM5_ESM.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06802-1/MediaObjects/41586_2023_6802_MOESM5_ESM.xlsx"
# Melendez 2024: Appendix S1 (elastic-net weights of the CSF clock), member of the supplementary zip of Europe PMC
d=melendez2024_clock
get $d/PMC11488306_supplementaryFiles.zip "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11488306/supplementaryFiles"
[[ -s data/raw/$d/ACEL-23-e14230-s002.xlsx ]] || unzip -oq data/raw/$d/PMC11488306_supplementaryFiles.zip ACEL-23-e14230-s002.xlsx -d data/raw/$d
# Oh 2025: Supplementary Tables (ST2 Olink assays, ST3 LASSO weights of the organ models)
get oh2025_clock/41591_2025_3798_MOESM3_ESM.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-025-03798-1/MediaObjects/41591_2025_3798_MOESM3_ESM.xlsx"
# Goeminne 2025: Table S1 (elastic-net coefficients of the organ models) and Table S3 (standard deviations of Olink NPX)
for m in mmc2 mmc4; do
  get goeminne2025_clock/1-s2.0-S1550413124004017-$m.xlsx "https://ars.els-cdn.com/content/image/1-s2.0-S1550413124004017-$m.xlsx"
done
