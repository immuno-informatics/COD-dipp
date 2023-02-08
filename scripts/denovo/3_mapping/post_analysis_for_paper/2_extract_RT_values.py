"""
    extracts retention time values from mgf files and saves it to TSV file
"""

import glob
import pandas as pd
import os

root = "/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis"
output = os.path.join(root, "data_all/tables/all_specra_rt_exp.tsv")

studies= [
"PXD000394",
"MSV000080527-84172-84442",
"PXD001898-PASS00270",
"PXD003790",
"PXD004023",
"PXD004233",
"PXD004233-ETHCD",
"PXD004746",
"PXD004894",
"PXD005231",
"PXD006939",
"PXD007203",
"PXD007596",
"PXD007860",
"PXD008937",
"PXD009531",
"PXD009738",
"PXD009749-9753-7935-9750-9751-9752-9754-9755",
"PXD010808",
"PXD011257",
"PXD011628",
"PXD011723",
"PXD011766",
"PXD012083",
"PXD012308",
"PXD013057",
"PXD014017"]

spec = []
rts = []
sample_names = []
study_names = []

for study in studies:
    print("Processing study: " + study)
    samples = glob.glob(os.path.join(root, study, "sample_*"))
    for sample in samples:
        sample_name = os.path.basename(sample)
        print("Processing sample: " + sample_name)
        mgf_files = glob.glob(os.path.join(root, study, sample, "denovo/*.mgf"))
        for file in mgf_files:
            if not file.endswith("spectrums.mgf"):
                with open(file) as fh:
                    for line in fh:
                        if line.startswith("TITLE="):
                            spec.append(line[6:].split(" ")[0])
                            sample_names.append(sample_name)
                            study_names.append(study)
                        if line.startswith("RTINSECONDS="):
                            rts.append(line[12:-1])
    print("\n")


df = pd.DataFrame({"spectrum": spec, "RT exp": rts, "sample": sample_names, "pride id": study_names})
df.to_csv(output, sep="\t", header=True, index=False)