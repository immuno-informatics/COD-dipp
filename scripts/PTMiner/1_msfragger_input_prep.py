import argparse
import os
import pickle as pkl
import re

from Bio import SeqIO


def args():
    parser = argparse.ArgumentParser(
        description=
        'Makes MSFragger output (pepXML) compatible with PTMiner 1.1.2'\
        'by converting protein accessions to numerical IDs')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '-mgf',
        type=str,
        nargs="+",
        help='Path to mgf files',
        required=True)
    requiredNamed.add_argument(
        '-pepXML',
        type=str,
        nargs="+",
        help='Path to MSFragger pepXML output',
        required=True)

    requiredNamed.add_argument(
        '-fasta', type=str, help='Path to input fasta', required=True)

    requiredNamed.add_argument(
        '-decoytag',
        type=str,
        help='Decoy tag used in the database (i.e rev_)',
        required=True)

    requiredNamed.add_argument(
        '-output_dir',
        type=str,
        help='Path to output directory location. '\
        'Outputs: \n- pepXML files with protein accessions as numerical IDs\n'\
        '- database___numerical___.fasta: database with numerical IDs\n'\
        '- db_headers___numerical___.pkl: database {"Numerical ID": "Protein Accession"} dict'\
        '- spectra.pkl: spectra {"modified spectrum name": "original spectrum name"} dict'\
        '- pram.parameters: PTMiner parameters file',
        required=True)

    return parser


def remove_hashtag(line, pattern):
    pat = re.findall(pattern, line)[0]
    new_pat = pat.replace('#', '')
    pat_el = pattern.split('"')
    line = re.sub(pattern, pat_el[0] + '"' + new_pat + '"' + pat_el[2], line)
    return line


if __name__ == "__main__":
    arguments = args().parse_args()
    files_input = arguments.pepXML
    fasta_input = arguments.fasta
    output_dir = arguments.output_dir
    decoytag = arguments.decoytag
    mgf_files = arguments.mgf

    suffix = "___numerical___"
    fasta_output = os.path.join(output_dir, "database" + suffix + ".fasta")
    anno_file = os.path.join(output_dir, "db_headers" + suffix + ".pkl")

    dict_conv = {}
    dict_conv_rev = {}
    pepXML_files = [] # for the param file
    with open(fasta_output, "w") as fh:
        for i, hit in enumerate(SeqIO.parse(fasta_input, format="fasta")):
            header = hit.description
            if header.startswith(decoytag):
                continue
            dict_conv[header] = f"sp|{i}|X"
            dict_conv_rev[f"sp|{i}|X"] = header
        for i, hit in enumerate(SeqIO.parse(fasta_input, format="fasta")):
            header = hit.description
            if not header.startswith(decoytag):
                continue
            dict_conv[header] = decoytag + dict_conv[header[4:]]
            dict_conv_rev[decoytag + dict_conv[header[4:]]] = header
        for i, hit in enumerate(SeqIO.parse(fasta_input, format="fasta")):
            header = hit.description
            fh.write(f">{dict_conv[header]}\n{str(hit.seq)}\n")

    pkl.dump(dict_conv_rev, open(anno_file, "wb"))

    for file_input in files_input:
        if os.stat(file_input).st_size == 0:
            continue
        replace_hashtag = False
        basename = os.path.basename(file_input)
        if '#' in basename:
            basename = basename.replace('#', '')
            replace_hashtag = True
        file_output = basename.rsplit('.', 1)
        file_output[0] = file_output[0] + suffix
        file_output = ".".join(file_output)
        pepXML_files.append(file_output)
        file_output = os.path.join(output_dir, file_output)
        with open(file_input) as fh:
            with open(file_output, "w") as out:
                for line in fh:
                    if replace_hashtag and 'spectrum="' in line:
                        pattern = 'spectrum="(.+)" end_scan'
                        line = remove_hashtag(line, pattern)
                        out.write(line)
                    elif replace_hashtag and "msms_run_summary base_name=" in line:
                        pattern = 'msms_run_summary base_name="(.+)" raw_data_type'
                        line = remove_hashtag(line, pattern)
                        out.write(line)
                    elif replace_hashtag and "summary_xml=" in line:
                        pattern = 'summary_xml="(.+)" xsi:schemaLocation'
                        line = remove_hashtag(line, pattern)
                        out.write(line)
                    elif replace_hashtag and "base_name=" in line:
                        pattern = 'base_name="(.+)" precursor_mass_type'
                        line = remove_hashtag(line, pattern)
                        out.write(line)
                    elif 'protein="' in line:
                        line = line.replace("&apos;", "'")
                        prot = re.findall('protein="(.+)" peptide_prev_aa', line)[0]
                        line = re.sub('protein="(.+)" peptide_prev_aa', 'protein="' + dict_conv[prot] + '" peptide_prev_aa', line)
                        out.write(line)
                    else:
                        out.write(line)


line_return = """
"""
spectrum_conv = {}
mgf_files_string = f"dataset_number = {len(mgf_files)}{line_return}"
for i, file_ in enumerate(mgf_files):
    i += 1
    original_name = file_.replace("_uncalibrated", "")
    original_name = original_name.replace("_calibrated", "")
    if '#' in original_name:
        original_name = original_name.replace("#", "")
        with open(original_name, 'w') as out:
            with open(file_) as fh:
                for line in fh:
                    if line.startswith("TITLE"):
                        old_title = line.split("=")[1]
                        line = line.replace('#', '')
                        new_title = line.split("=")[1]
                        spectrum_conv[new_title.strip("\n")] = old_title.strip("\n")
                    out.write(line)
    elif "_uncalibrated" in file_ or "_calibrated" in file_:
        with open(file_) as fh:
            for line in fh:
                if line.startswith("TITLE"):
                    old_title = line.split("=")[1].strip("\n")
                    spectrum_conv[old_title] = old_title
        os.system(f" mv {file_} {original_name}")
    spectrum_conv[os.path.basename(original_name)] = os.path.basename(file_)
    file_ = original_name
    file_path = os.path.basename(file_)
    mgf_files_string = mgf_files_string + rf"dataset_filename_{i} = Z:\data\{file_path}{line_return}"

pkl.dump(spectrum_conv, open("spectra.pkl", "wb"))

if len(pepXML_files) == 0:
    print("No pepXML files were found. exiting ...")
    exit(1)

pepXML_files_string = f"pep_ident_number = {len(pepXML_files)}{line_return}"
for i, file_ in enumerate(pepXML_files):
    i += 1
    pepXML_files_string = pepXML_files_string + rf"pep_ident_filename_{i} = Z:\data\{file_}{line_return}"


parm_file = rf"""
[Search]
# pFind 3 (.spectra),pFind 2.8 (.txt),MSFragger (.pepXML),Sequest (.txt) and PTMiner (.txt)
search_result_format = MSFragger (.pepXML)
#0 = close search, 1 = open search
open_search = 1

# # peptide_tol and peptide_tol_type
# peptide_tol_type MUST be 0 = Da
precursor_matching_tolerance = 300
precursor_matching_tolerance_type = 0

# the four parameters refer to spectrum precision
# 0 = Da, 1 = PPM
precursor_tol = 10
precursor_tol_type = 1
fragment_tol = 10
fragment_tol_type = 1

# fixed modifications in peptide identification process
fixed_mod_number = 0

# variable modifications in peptide identification process
var_mod_number = 0

# dataset
# now only support mgf format
dataset_format = mgf
{mgf_files_string}

# peptide identification results
pep_ident_format = MSFragger (.pepXML)
{pepXML_files_string}

[Fdr]
# do fdr control or not. 0=No, 1=Yes
is_fdr_control = 1

# decoy tag
decoy_tag = {decoytag}
# fdr threshold
fdr_threshold = 0.01

# fdr method, global=1,separate=2,transferred=3
fdr_method = 3

[Localization]
# do localization or not. 0=No, 1=Yes
is_localized = 1

# the minimum number of modifications
min_mod_number = 5

# use the prior probability or not.  0=No, 1=Yes
use_prior = 1

# the method to filter localization results
# 0 = probability, 1 = flr
filter_method = 1

# filter threshold
filter_threshold = 0.01

[Annotation]
# do annotation or not. 0=No, 1=Yes
is_annotated = 1
protein_database = z:\data\{os.path.basename(fasta_output)}

[Output]
# output path
output = Z:\data\
"""

with open("pram.parameters", "w") as out:
    out.write(parm_file)
