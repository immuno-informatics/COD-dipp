"""
    Checks if Denovo intronic peptides are shared between samples
    since authors DOI: 10.1126/scitranslmed.aau5516 claimed that
    TSAs (Tumor Specific Antigens) comming from non coding regiosn are shared
    between patients
"""

import pandas as pd
import os
import random

NETMHCPAN = "/net/archive/groups/plgg_iccvs/tools/netmhccons/netMHCpan-4.0/netMHCpan"
NETMHCPAN_PARSER = "$HOME/git_repositories/covid19/src/NetMHCPan/netmhcpan_parser.py"
POPULATION_COVERAGE_SCRIPT = "/net/archive/groups/plgg_iccvs/tools/population_coverage/calculate_population_coverage.py"
IMMUNOPRED = "/net/archive/groups/plgg_iccvs/tools/immunogenicity/predict_immunogenicity.py"

ALPHABET = "ARNDCEQGHILKMFPSTWYV"


def random_peptide_sequence():
    length = random.choice([8, 9, 10, 11, 12])
    return ''.join(random.choice(ALPHABET) for _ in range(length))


def get_hla_binding(peptide_list, supertypes=None):
    pep_file = "/tmp/pepfile.txt"
    netmhcpan_output = "/tmp/pepfile_netmhcpan"
    netmhcpan_parsed = "/tmp/pepfile_netmhcpan_parsed"

    if not supertypes:
        supertypes = "HLA-A02:01,HLA-A01:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01,HLA-B15:01"

    with open(pep_file, "w") as fh:
        fh.write("\n".join(peptide_list))

    os.system(
        f"{NETMHCPAN} -f {pep_file} -a {supertypes} -inptype 1 -l 8,9,10,11,12 > {netmhcpan_output}"
    )
    os.system(
        f"python {NETMHCPAN_PARSER} {netmhcpan_output} {netmhcpan_parsed}")

    bp_binding = pd.read_csv(netmhcpan_parsed, sep="\t")
    os.system(f"rm -rf {pep_file} {netmhcpan_output} {netmhcpan_parsed}")
    return bp_binding


def main():
    # yapf: disable
    output_root = "/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all/denovo_annotation/overlap_new"
    output_bindings = os.path.join(output_root, "random_peptides_bindings.tsv")
    output_bindings_best = os.path.join(output_root, "random_peptides_best_bindings.tsv")
    # yapf: enable
    random_peptides = [random_peptide_sequence() for i in range(239)]
    bp_binding = get_hla_binding(random_peptides)
    bp_binding_best = bp_binding[~bp_binding["BindLevel"].isnull(
    )].sort_values("%Rank").drop_duplicates('Peptide')

    # write binding affinities to tsv files
    bp_binding.to_csv(output_bindings, sep="\t", header=True, index=False)
    bp_binding_best.to_csv(output_bindings_best,
                           sep="\t",
                           header=True,
                           index=False)


if __name__ == "__main__":
    main()
