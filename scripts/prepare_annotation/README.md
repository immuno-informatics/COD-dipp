# Create the annotation files

## Create the gene annotation files

```
conda env create -f generate_annotation.yml

conda activate cruzdb
python generate_annotation.py
conda deactivate
```

The generate_annotation.py script creates 3 files:
    - UCSC_knownGene_hg38.tsv: The UCSC hg38 gene table
    - UCSC_knownGene_hg38_features.tsv: necessary to locate on which features (UTR, Intron, Exon) denovo peptides are located
    - df_features_inframes.tsv: necessary for checking if denovo peptides are in frame with upstream exons

## Create the 3 frame translation (3FT) database

```
conda env create -f genes3FT_generator.yml

conda activate database
Rscript --vanilla genes3FT_generator.R
conda deactivate
```

The genes3FT_generator.R script creates:
    - 3FTgenes_coding.fasta: 3FT translation of all the coding pre-spliced mRNA

## Create PTMiner singularity image
 - check the how to prepare [PTMiner_singularity_image](./PTMiner_singularity_image)

## download TITER (Translation Initiation Site Predictor) model
The data can be downloade from https://github.com/zhangsaithu/titer
the data must be organized like the following
```
git clone https://github.com/zhangsaithu/titer
cd titer
mv codes/dict_piror_front_Gaotrain.npy ./
unzip data
rm -rf README.md codes __MACOSX

# how the titer dir should look like
titer
  ├── data
  │   ├── neg_seq_test_all_upstream.npy
  │   └── pos_seq_test.npy
  ├── dict_piror_front_Gaotrain.npy
  └── model
      ├── bestmodel_0.hdf5
      ├── bestmodel_10.hdf5
      ├── bestmodel_11.hdf5
      ├── bestmodel_12.hdf5
      ├── bestmodel_13.hdf5
      ├── bestmodel_14.hdf5
      ├── bestmodel_15.hdf5
      ├── bestmodel_16.hdf5
      ├── bestmodel_17.hdf5
      ├── bestmodel_18.hdf5
      ├── bestmodel_19.hdf5
      ├── bestmodel_1.hdf5
      ├── bestmodel_20.hdf5
      ├── bestmodel_21.hdf5
      ├── bestmodel_22.hdf5
      ├── bestmodel_23.hdf5
      ├── bestmodel_24.hdf5
      ├── bestmodel_25.hdf5
      ├── bestmodel_26.hdf5
      ├── bestmodel_27.hdf5
      ├── bestmodel_28.hdf5
      ├── bestmodel_29.hdf5
      ├── bestmodel_2.hdf5
      ├── bestmodel_30.hdf5
      ├── bestmodel_31.hdf5
      ├── bestmodel_3.hdf5
      ├── bestmodel_4.hdf5
      ├── bestmodel_5.hdf5
      ├── bestmodel_6.hdf5
      ├── bestmodel_7.hdf5
      ├── bestmodel_8.hdf5
      └── bestmodel_9.hdf5
```
