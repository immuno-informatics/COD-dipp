"""
    This script is for removing PTM assigned by the open search that are not
    true using PTMProphet post-processing.
"""
import re
import sys

from pyteomics import mass

import numpy as np
from Bio import SeqIO


def verify_unimod(site, position):
    if not re.match("^[A-Z]$", site):
        if site not in ['N-term', 'C-term']:
            raise ValueError('{}: Not A unimod site'.format(site))

    positions = [
        "Anywhere", "Any N-term", "Any C-term", "Protein C-term",
        "Protein N-term"
    ]
    if position not in positions:
        raise ValueError('{}: Not A unimod position'.format(position))


def _unimod_parser_map():
    '''
        create UNIMOD dictionairy {mono_mass -> full name (mono_mass)}
    '''
    unimod_dict = {}
    # connect UNIMOD database
    unimod_db = mass.Unimod(source='http://www.unimod.org/xml/unimod.xml')
    step = 0.0001
    for mod in unimod_db.mods:
        for site in mod['specificity']:
            err = 0.01 # dalton
            for mono_mass in np.arange(round(mod['mono_mass'] - err, 4),
                                       round(mod['mono_mass'] + err, 4) + step,
                                       step):
                key = str(round(mono_mass, 4))
                value = mod['full_name'] + '(' + str(
                    round(mod['mono_mass'], 3)) + ')'
                try:
                    unimod_dict[key].append(value)
                except KeyError:
                    unimod_dict[key] = []
                    unimod_dict[key].append(value)
    return unimod_dict


def _unimod_parser_():
    """
        return:
            (dict) {
                "PTM name/title": {
                    "site": ["amino acid" / "N-term" / C-term"],
                    "Position": ["Anywhere", "Any N-term", "Any C-term",
                                 "Protein N-term", "Protein C-term"]
                    "Mono_mass": (float)
                }
            }
    """
    desc_dict = {}
    # connect UNIMOD database
    unimod_db = mass.Unimod(source='http://www.unimod.org/xml/unimod.xml')
    for mod in unimod_db.mods:
        for name in ['full_name', 'title']:
            desc_dict[mod[name]] = {}
            desc_dict[mod[name]]['mono_mass'] = mod['mono_mass']
            site = [x['site'] for x in mod['specificity']]
            pos = [x['position'] for x in mod['specificity']]
            desc_dict[mod[name]]['site'] = site
            desc_dict[mod[name]]['position'] = pos
    return desc_dict


def prot_term_checker(peptide, mapped_proteins, db_dict):
    """
        checks if the peptide falls on the Protein N/C-term
    """
    conv_dict = {}
    conv_dict['N-term'] = 0
    conv_dict['C-term'] = len(peptide) - 1

    is_prot_nterm = []
    is_prot_cterm = []

    for protein in mapped_proteins:
        try:
            seq = str(db_dict[protein].seq)
        except KeyError:
            return None, None
        try:
            pep_start = seq.index(peptide) + 1
        except ValueError:
            peptide = peptide.replace('I', 'L')
            seq = seq.replace('I', 'L')
            pep_start = seq.index(peptide) + 1
        except ValueError:
            return None, None
        pep_end = pep_start + conv_dict['C-term']
        protein_end = len(seq)
        is_prot_cterm.append(pep_end == protein_end)
        is_prot_nterm.append(pep_start == 1)
    is_prot_nterm = any(is_prot_nterm)
    is_prot_cterm = any(is_prot_cterm)
    return is_prot_nterm, is_prot_cterm


def get_modified_aa(modified_pep):
    """
        Given a PTMProphet format modified peptide, return a (dict)
        {
        "amino acid letter": [position of the modification starting from 1],
         ...
        }
    """
    aa_dict_ = {}
    aa_index = 0
    for index, aa in enumerate(modified_pep):
        if aa.islower():
            aa = aa.upper()
            try:
                aa_dict_[aa].append(aa_index)
            except KeyError:
                aa_dict_[aa] = [index + 1]
    return aa_dict_


def verify_ptm(peptide, unimod, mod_name, aa_dict, is_prot_nterm,
               is_prot_cterm):
    """
        Given a modification name, a dict of a peptide modified positions
        returns :
            - new moodification string format:
             "position modification_name(average mass)"

    """
    new_mod = None

    original_mod = mod_name
    mod_name = mod_name.rsplit('(', 1)[0]
    mod_mass = original_mod.rsplit('(', 1)[1][:-1]
    try:
        mod_sites = unimod[mod_name]['site']
        mod_pos = unimod[mod_name]['position']
    except KeyError:
        return new_mod
    for site, position in zip(mod_sites, mod_pos):
        verify_unimod(site, position)
        original_site = site
        if "N-term" in site:
            site = 0
        elif "C-term" in site:
            site = len(peptide) - 1
        try:
            for pep_pos_loc in aa_dict[site]:
                if 'N-term' in position:
                    if pep_pos_loc != 1:
                        pep_pos_loc = -1
                    if 'Prot' in position and not is_prot_nterm:
                        pep_pos_loc = -1
                elif 'C-term' in position:
                    if pep_pos_loc != len(peptide):
                        pep_pos_loc = -1
                    if 'Prot' in position and not is_prot_cterm:
                        pep_pos_loc = -1
                if pep_pos_loc != -1:
                    new_mod = mod_mass + '@' + str(pep_pos_loc) 
        except KeyError:
            pass
    return new_mod


def unimod_annotate_tsv(input_file, database_file):    
    DB = SeqIO.index(database_file, 'fasta')
    unimod = _unimod_parser_()
    unimod_map = _unimod_parser_map()
    list_of_mods = {}

    # pattern = re.compile('.+\(([0-9]+\.[0-9]+)\)')
    fh = open(input_file)

    cols = {}
    header = fh.readline().strip('\n').split('\t')
    for index, col in enumerate(header):
        cols[col] = index
    modifications = []
    for line in fh:
        line = line.strip('\n').split('\t')
        modified_pep = line[cols['best_locs']]
        peptide = line[cols['peptide']]
        if modified_pep != '':
            try:
                massdiff = f"{round(float(line[cols['massdiff']]), 4):.4f}"
                obs_mods = list(set(unimod_map[massdiff]))
            except KeyError:
                modifications.append([])
                continue
        else:
            modifications.append([])
            continue
        # get the modified amino acids
        aa_dict = get_modified_aa(modified_pep)
        # get a listed of the mapped proteins
        list_mapped_prot = [line[cols['protein']].split(' ')[0]]
        # check if it's a Protein N/C-term peptide
        is_prot_nterm, is_prot_cterm = prot_term_checker(
            peptide, list_mapped_prot, DB)
        # in case something went wrong
        if is_prot_nterm is None or is_prot_cterm is None:
            modifications.append([])
            continue
        new_mod_line = []
        # for each suggested unimod modification
        for mod in obs_mods:
            new_mod = verify_ptm(peptide, unimod, mod, aa_dict,
                                                is_prot_nterm, is_prot_cterm)
            # if the PTM sis verified
            if new_mod:
                new_mod_line.append(new_mod)

        # in case no PTM got verified
        if not new_mod_line:
            modifications.append([])
            continue
        # replace the value of "Observed Modifications"
        modifications.append(new_mod_line)
    fh.close()

    return modifications