 = ''
mztab_template += ''

def create_mztab_output():
    mztab_template = []
    mztab_template.append(['MTD', 'mzTab-version', '1.0 rc5'])
    mztab_template.append(['MTD', 'mzTab-mode', 'Summary'])
    mztab_template.append(['MTD', 'mzTab-type', 'Identification'])
    mztab_template.append(['MTD', 'description', 'mzTab output for Scavager'])
    mztab_template.append(['MTD', 'ms_run[1]-location', mzml_path])
    mztab_template.append(['MTD', 'protein_search_engine_score[1]', '[MS, MS:XXX, Scavager:Protein score]'])
    mztab_template.append(['MTD', 'psm_search_engine_score[1]', '[MS, MS:XXX, Scavager:PSM score]'])
    mztab_template.append(['MTD', 'fixed_mod[1]', '[UNIMOD, UNIMOD:XXX, Carbamido]'])
    mztab_template.append(['MTD', 'variable_mod[1]', '[UNIMOD, UNIMOD:XXX, Oxid]'])