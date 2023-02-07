import os
import csv
from itertools import groupby
from functools import reduce


'''
-------------------------------------------------------------------------------------------------------
                        DEFENSE FINDER (https://github.com/mdmparis/defense-finder)
-------------------------------------------------------------------------------------------------------
'''


def get_best_solution(tmp_dir):
    results = os.listdir(tmp_dir)
    acc = []
    for family_dir in results:
        family_path = os.path.join(tmp_dir, family_dir)
        acc = acc + parse_best_solution(family_path)
    return format_best_solution(acc)


def parse_best_solution(dir):
    # attention of macsyfinder number of header lines output
    tsv_file = open(os.path.join(dir, 'best_solution.tsv'))
    tsv = csv.reader(tsv_file, delimiter='\t')
    data = []
    for row in tsv:
        data.append(row)
    tsv_file.close()
    if "No Systems found" in ' '.join(data[3]):
        return []
    else:
        data = data[4:]
        header = data.pop(0)
        out = []
        for lst in data:
            if not lst:
                continue
            line_as_dict = {}
            for idx, val in enumerate(header):
                line_as_dict[val] = lst[idx]
            out.append(line_as_dict)
        return out


def get_best_solution_keys():
    return [
            'replicon', 'hit_id', 'gene_name',
            'hit_pos', 'model_fqn', 'sys_id', 'sys_loci', 'locus_num',
            'sys_wholeness', 'sys_score', 'sys_occ', 'hit_gene_ref',
            'hit_status', 'hit_seq_len', 'hit_i_eval', 'hit_score',
            'hit_profile_cov', 'hit_seq_cov', 'hit_begin_match',
            'hit_end_match', 'counterpart', 'used_in'
            ]


def format_best_solution(p):
    out = []
    for lst in p:
        gene_ref = lst['model_fqn']
        gene_ref_elements = gene_ref.split('/')
        type = gene_ref_elements[2]
        subtype = gene_ref_elements[3]
        native_keys = list(filter(lambda i: i not in ['type', 'subtype'], get_best_solution_keys()))
        new_line = {key: lst[key] for key in native_keys}
        new_line['type'] = type
        new_line['subtype'] = subtype
        out.append(new_line)
    return out


def export_defense_finder_genes(defense_finder_genes, outdir, orgName, dfmodelV):
    defense_finder_genes_list = defense_finder_genes_to_list(defense_finder_genes)
    write_defense_finder_genes(defense_finder_genes_list, outdir, orgName, dfmodelV)


def defense_finder_genes_to_list(defense_finder_genes):
    header = get_best_solution_keys()
    out = [header]
    for g in defense_finder_genes:
        lst = []
        for key in header:
            lst.append(g[key])
        out.append(lst)
    return out


def write_defense_finder_genes(defense_finder_genes_list, outdir, orgName, dfmodelV):
    filepath = os.path.join(outdir, orgName+"_DefenseFinder_"+dfmodelV+"_genes.tsv")
    with open(filepath, 'w') as defense_finder_genes_file:
        write = csv.writer(defense_finder_genes_file, delimiter='\t')
        write.writerows(defense_finder_genes_list)


def get_hit_sort_attr(hit):
    #
    return hit['hit_id']


def remove_duplicates(hmmer_hits):
    temp_list = []
    for i in hmmer_hits:
        if i not in temp_list:
            temp_list.append(i)
    return temp_list


def export_defense_finder_hmmer_hits(tmp_dir, outdir, orgName, dfmodelV):
    paths = get_hmmer_paths(tmp_dir)
    hmmer_hits = []
    for path in paths:
        d = parse_hmmer_results_file(path)
        hmmer_hits = hmmer_hits + remove_duplicates(d)
    sorted_hmmer_hits = sorted(hmmer_hits, key=get_hit_sort_attr)
    hmmer_hits_list = hmmer_to_list(sorted_hmmer_hits)
    write_defense_finder_hmmer(hmmer_hits_list, outdir, orgName, dfmodelV)


def write_defense_finder_hmmer(hmmer_hits_list, outdir, orgName, dfmodelV):
    filepath = os.path.join(outdir, orgName+"_DefenseFinder_"+dfmodelV+"_hmmer.tsv")
    with open(filepath, 'w') as defense_finder_hmmer_file:
        write = csv.writer(defense_finder_hmmer_file, delimiter='\t')
        write.writerows(hmmer_hits_list)
        defense_finder_hmmer_file.close()


def get_hmmer_keys():
    #
    return ['hit_id', 'replicon', 'position_hit', 'hit_sequence_length', 'gene_name']


def parse_hmmer_results_file(path):
    tsv_file = open(path)
    tsv = csv.reader(tsv_file, delimiter='\t')
    data = []
    for row in tsv:
        if not row[0].startswith('#'):
            data.append(row)
    tsv_file.close()
    out = []
    for lst in data:
        if not lst:
            continue
        line_as_dict = {}
        for idx, val in enumerate(get_hmmer_keys()):
            line_as_dict[val] = lst[idx]
        out.append(line_as_dict)
    return out


def get_hmmer_paths(results_dir):
    family_dirs = os.listdir(results_dir)
    files = []
    for family_dir in family_dirs:
        hmmer_results_dir = os.path.join(results_dir, family_dir, 'hmmer_results')
        with os.scandir(hmmer_results_dir) as it:
            for entry in it:
                if entry.name.endswith('extract') and entry.is_file():
                    files.append(entry)
    return list(map(lambda i: i.path, files))


def hmmer_to_list(hmmer_hits):
    header = get_hmmer_keys()
    out = [header]
    for s in hmmer_hits:
        lst = []
        for key in header:
            lst.append(s[key])
        out.append(lst)
    return out


def export_defense_finder_systems(defense_finder_genes, outdir, orgName, dfmodelV):
    systems = build_defense_finder_systems(defense_finder_genes)
    systems_list = systems_to_list(systems)
    write_defense_finder_systems(systems_list, outdir, orgName, dfmodelV)


def systems_to_list(systems):
    header = get_system_keys()
    out = [header]
    for s in systems:
        lst = []
        for key in header:
            lst.append(s[key])
        out.append(lst)
    return out


def write_defense_finder_systems(systems_list, outdir, orgName, dfmodelV):
    filepath = os.path.join(outdir, orgName+"_DefenseFinder_"+dfmodelV+"_systems.tsv")
    with open(filepath, 'w') as defense_finder_systems_file:
        write = csv.writer(defense_finder_systems_file, delimiter='\t')
        write.writerows(systems_list)


def get_system_keys():
    #
    return ['sys_id', 'type', 'subtype', 'sys_beg', 'sys_end', 'protein_in_syst', 'genes_count', 'name_of_profiles_in_sys']


def projection(val):
    #
    return val['sys_id']


def build_defense_finder_systems(defense_finder_genes):
    system_goups = [list(it) for k, it in groupby(defense_finder_genes, projection)]
    out = []
    for system_group in system_goups:
        item = {}
        first_item = system_group[0]
        last_item = system_group[-1]
        item['sys_id'] = first_item['sys_id']
        item['sys_beg'] = first_item['hit_id']
        item['sys_end'] = last_item['hit_id']
        item['type'] = first_item['type']
        item['subtype'] = first_item['subtype']
        item['protein_in_syst'] = reduce(lambda acc, s: acc + ',' + s['hit_id'], system_group, '')[1:]
        item['genes_count'] = len(system_group)
        item['name_of_profiles_in_sys'] = reduce(lambda acc, s: acc + ',' + s['gene_name'], system_group, '')[1:]
        out.append(item)
    return out
