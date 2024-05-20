import numpy as np
import pandas as pd
import primer3
from primer3.bindings import calcHeterodimer
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from pyensembl import EnsemblRelease
import requests
import urllib.request

from Bio import Entrez
from Bio import SeqIO

def ensemble_info(gene_name):
    """
    Get information from ensemble database.
    Release 109 uses human reference genome GRCh38
    """

    data = EnsemblRelease(109)

    # get gene information
    gene = data.genes_by_name(gene_name)[0]
    strand = gene.strand
    chromosome_name = gene.contig

    gene_dict = gene.to_dict()

    start_gene = gene.start - 1500
    end_gene = gene.end + 1500

    # get sequence
    url = (
        f'https://rest.ensembl.org/sequence/region/human/'
        f'{chromosome_name}:{start_gene}-{end_gene}'
        '?content-type=application/json;version=109'
    )
    data = requests.get(url)
    gene_seq = data.json()["seq"]

    coord_list = list([i for i in range(start_gene, end_gene + 1)])


    # prepare sequence for reverse strand gene location
    if strand == "-":
        compl_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
        gene_seq = [compl_dict[l] for l in gene_seq][::-1]
        coord_list = coord_list[::-1]

    ensemble_gene_seq = "".join(gene_seq)

    del gene_dict["gene_name"]
    del gene_dict["gene_id"]
    del gene_dict["genome"]

    return ensemble_gene_seq, gene_dict, strand

def ncbi_info(ncbi_id, strand, search_sequence):
    Entrez.email = "Your.Name.Here@example.org"
    handle = Entrez.efetch(db="gene", id=ncbi_id, rettype="gb", retmode="text")
    information = handle.read()
    gene_nt_id = 'NC_' + information.split('Annotation: ')[1].split('NC_')[1].split(' ')[0]
    start_gene_nt, end_gene_nt = information.split('Annotation: ')[1].split('(')[1].split(')')[0].split(',')[0].split('..')

    start_gene_nt = int(start_gene_nt) - 1500
    end_gene_nt = int(end_gene_nt) + 1500

    if strand == '+':
        strand_nt = 1
    else:
        strand_nt = 2

    Entrez.email = "Your.Name.Here@example.org"
    handle = Entrez.efetch(db="nucleotide", id=gene_nt_id, rettype="gb", retmode="text", 
                        seq_start=start_gene_nt, seq_stop=end_gene_nt, strand=strand_nt)
    seq_record = [seq_record for seq_record in SeqIO.parse(handle, "gb")][0]
    refseq_sequence = str(seq_record.seq)

    if strand == '+':
        amplicon_start = len(refseq_sequence.split(search_sequence)[0]) + int(start_gene_nt)
        amplicon_end = amplicon_start + len(search_sequence)
    else:
        amplicon_end = end_gene_nt - len(refseq_sequence.split(search_sequence)[0])
        amplicon_start = amplicon_end - len(search_sequence)

    return amplicon_start, amplicon_end, gene_nt_id

def guide_info(guide_seq, strand, ensemble_gene_seq):
    """
    Search guide position
    """

    compl_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    if ';' not in guide_seq:
        right_guide = guide_seq[-20:]
        left_guide = ''.join([compl_dict[i] for i in guide_seq[:20]][::-1])
        guide_seq = ';'.join([left_guide, right_guide])

        right_guide_save = right_guide
        left_guide_save = left_guide

    if ';' not in guide_seq:
        if strand == '-':
            guide = ''.join([compl_dict[i] for i in guide_seq][::-1])
        else:
            guide = guide_seq
    else:
        left_guide, right_guide = guide_seq.split(';')
        left_guide = left_guide.strip()
        right_guide = right_guide.strip()

        right_guide_save = right_guide
        left_guide_save = left_guide
        
        if left_guide not in ensemble_gene_seq:
            left_guide = ''.join([compl_dict[i] for i in left_guide][::-1])
        if right_guide not in ensemble_gene_seq:
            right_guide = ''.join([compl_dict[i] for i in right_guide][::-1])
        
        left_guide_start = len(ensemble_gene_seq.split(left_guide)[0])
        right_guide_start = len(ensemble_gene_seq.split(right_guide)[0])
        
        if left_guide_start>right_guide_start:
            guide = (right_guide 
                    + ensemble_gene_seq.split(right_guide)[1].split(left_guide)[0] 
                    + left_guide)
        else:
            guide = (left_guide 
                    + ensemble_gene_seq.split(left_guide)[1].split(right_guide)[0] 
                    + right_guide)  

    return guide

def blast_primers(search_sequence, 
                  primer5_start, 
                  primer5_end, 
                  primer3_start, 
                  primer3_end, 
                  product_size_min, 
                  product_size_max):
    # url = f'https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?LINK_LOC=bookmark&INPUT_SEQUENCE={search_sequence}&PRIMER5_START={primer5_start}&PRIMER5_END={primer5_end}&PRIMER3_START={primer3_start}&PRIMER3_END={primer3_end}&OVERLAP_5END=7&OVERLAP_3END=4&PRIMER_PRODUCT_MIN={product_size_min}&PRIMER_PRODUCT_MAX={product_size_max}&PRIMER_NUM_RETURN=1000&PRIMER_MIN_TM=57.0&PRIMER_OPT_TM=60.0&PRIMER_MAX_TM=63.0&PRIMER_MAX_DIFF_TM=2&PRIMER_ON_SPLICE_SITE=0&SEARCHMODE=0&SPLICE_SITE_OVERLAP_5END=7&SPLICE_SITE_OVERLAP_3END=4&SPLICE_SITE_OVERLAP_3END_MAX=8&SPAN_INTRON=off&MIN_INTRON_SIZE=1000&MAX_INTRON_SIZE=1000000&SEARCH_SPECIFIC_PRIMER=on&EXCLUDE_ENV=on&EXCLUDE_XM=on&TH_OLOGO_ALIGNMENT=off&TH_TEMPLATE_ALIGNMENT=off&ORGANISM=Homo%20sapiens&PRIMER_SPECIFICITY_DATABASE=refseq_representative_genomes&TOTAL_PRIMER_SPECIFICITY_MISMATCH=1&PRIMER_3END_SPECIFICITY_MISMATCH=1&MISMATCH_REGION_LENGTH=5&TOTAL_MISMATCH_IGNORE=6&MAX_TARGET_SIZE=4000&ALLOW_TRANSCRIPT_VARIANTS=off&HITSIZE=50000&EVALUE=30000&WORD_SIZE=7&MAX_CANDIDATE_PRIMER=500&PRIMER_MIN_SIZE=18&PRIMER_OPT_SIZE=20&PRIMER_MAX_SIZE=25&PRIMER_MIN_GC=20.0&PRIMER_MAX_GC=80.0&GC_CLAMP=0&NUM_TARGETS_WITH_PRIMERS=1000&NUM_TARGETS=20&MAX_TARGET_PER_TEMPLATE=100&POLYX=3&SELF_ANY=8.00&SELF_END=3.00&PRIMER_MAX_END_STABILITY=9&PRIMER_MAX_END_GC=5&PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00&PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00&PRIMER_MAX_SELF_ANY_TH=45.0&PRIMER_MAX_SELF_END_TH=35.0&PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0&PRIMER_PAIR_MAX_COMPL_END_TH=35.0&PRIMER_MAX_HAIRPIN_TH=24.0&PRIMER_MAX_TEMPLATE_MISPRIMING=12.00&PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00&PRIMER_PAIR_MAX_COMPL_ANY=8.00&PRIMER_PAIR_MAX_COMPL_END=3.00&PRIMER_MISPRIMING_LIBRARY=repeat/repeat_9606&NO_SNP=on&LOW_COMPLEXITY_FILTER=on&MONO_CATIONS=50.0&DIVA_CATIONS=1.5&CON_ANEAL_OLIGO=300.0&CON_DNTPS=0.6&SALT_FORMULAR=1&TM_METHOD=1&PRIMER_INTERNAL_OLIGO_MIN_SIZE=18&PRIMER_INTERNAL_OLIGO_OPT_SIZE=20&PRIMER_INTERNAL_OLIGO_MAX_SIZE=27&PRIMER_INTERNAL_OLIGO_MIN_TM=57.0&PRIMER_INTERNAL_OLIGO_OPT_TM=60.0&PRIMER_INTERNAL_OLIGO_MAX_TM=63.0&PRIMER_INTERNAL_OLIGO_MAX_GC=80.0&PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT=50&PRIMER_INTERNAL_OLIGO_MIN_GC=20.0&PICK_HYB_PROBE=off&NEWWIN=on&NEWWIN=on&SHOW_SVIEWER=true'
    url = (
        f'https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?LINK_LOC=bookmark'
        f'&INPUT_SEQUENCE={search_sequence}'
        f'&PRIMER5_START={primer5_start}'
        f'&PRIMER5_END={primer5_end}'
        f'&PRIMER3_START={primer3_start}'
        f'&PRIMER3_END={primer3_end}'
        f'&OVERLAP_5END=7'
        f'&OVERLAP_3END=4'
        f'&PRIMER_PRODUCT_MIN={product_size_min}'
        f'&PRIMER_PRODUCT_MAX={product_size_max}'
        f'&PRIMER_NUM_RETURN=1000'
        f'&PRIMER_MIN_TM=55.0'
        f'&PRIMER_OPT_TM=60.0'
        f'&PRIMER_MAX_TM=65.0'
        f'&PRIMER_MAX_DIFF_TM=3'
        f'&PRIMER_ON_SPLICE_SITE=0'
        f'&SEARCHMODE=0'
        f'&SPLICE_SITE_OVERLAP_5END=7'
        f'&SPLICE_SITE_OVERLAP_3END=4'
        f'&SPLICE_SITE_OVERLAP_3END_MAX=8'
        f'&SPAN_INTRON=off'
        f'&MIN_INTRON_SIZE=1000'
        f'&MAX_INTRON_SIZE=1000000'
        f'&SEARCH_SPECIFIC_PRIMER=on'
        f'&EXCLUDE_ENV=on'
        f'&EXCLUDE_XM=on'
        f'&TH_OLOGO_ALIGNMENT=off'
        f'&TH_TEMPLATE_ALIGNMENT=off'
        f'&ORGANISM=Homo%20sapiens'
        f'&PRIMER_SPECIFICITY_DATABASE=refseq_representative_genomes'
        f'&TOTAL_PRIMER_SPECIFICITY_MISMATCH=1'
        f'&PRIMER_3END_SPECIFICITY_MISMATCH=1'
        f'&MISMATCH_REGION_LENGTH=5'
        f'&TOTAL_MISMATCH_IGNORE=6'
        f'&MAX_TARGET_SIZE=4000'
        f'&ALLOW_TRANSCRIPT_VARIANTS=off'
        f'&HITSIZE=50000'
        f'&EVALUE=30000'
        f'&WORD_SIZE=7'
        f'&MAX_CANDIDATE_PRIMER=500'
        f'&PRIMER_MIN_SIZE=16'
        f'&PRIMER_OPT_SIZE=20'
        f'&PRIMER_MAX_SIZE=25'
        f'&PRIMER_MIN_GC=20.0'
        f'&PRIMER_MAX_GC=80.0'
        f'&GC_CLAMP=0'
        f'&NUM_TARGETS_WITH_PRIMERS=1000'
        f'&NUM_TARGETS=20'
        f'&MAX_TARGET_PER_TEMPLATE=100'
        f'&POLYX=3'
        f'&SELF_ANY=8.00'
        f'&SELF_END=3.00'
        f'&PRIMER_MAX_END_STABILITY=9'
        f'&PRIMER_MAX_END_GC=5'
        f'&PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00'
        f'&PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00'
        f'&PRIMER_MAX_SELF_ANY_TH=45.0'
        f'&PRIMER_MAX_SELF_END_TH=35.0'
        f'&PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0'
        f'&PRIMER_PAIR_MAX_COMPL_END_TH=35.0'
        f'&PRIMER_MAX_HAIRPIN_TH=24.0'
        f'&PRIMER_MAX_TEMPLATE_MISPRIMING=12.00'
        f'&PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00'
        f'&PRIMER_PAIR_MAX_COMPL_ANY=8.00'
        f'&PRIMER_PAIR_MAX_COMPL_END=3.00'
        f'&PRIMER_MISPRIMING_LIBRARY=repeat/repeat_9606'
        f'&NO_SNP=on'
        f'&LOW_COMPLEXITY_FILTER=on'
        f'&MONO_CATIONS=50.0'
        f'&DIVA_CATIONS=1.5'
        f'&CON_ANEAL_OLIGO=300.0'
        f'&CON_DNTPS=0.6'
        f'&SALT_FORMULAR=1'
        f'&TM_METHOD=1'
        f'&PRIMER_INTERNAL_OLIGO_MIN_SIZE=18'
        f'&PRIMER_INTERNAL_OLIGO_OPT_SIZE=20'
        f'&PRIMER_INTERNAL_OLIGO_MAX_SIZE=27'
        f'&PRIMER_INTERNAL_OLIGO_MIN_TM=57.0'
        f'&PRIMER_INTERNAL_OLIGO_OPT_TM=60.0'
        f'&PRIMER_INTERNAL_OLIGO_MAX_TM=63.0'
        f'&PRIMER_INTERNAL_OLIGO_MAX_GC=80.0'
        f'&PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT=50'
        f'&PRIMER_INTERNAL_OLIGO_MIN_GC=20.0'
        f'&PICK_HYB_PROBE=off'
        f'&NEWWIN=on'
        f'&NEWWIN=on'
        f'&SHOW_SVIEWER=true'
    )
 
    return url

def blast_results(job_id, 
                  gene_nt_id, 
                  amplicon_start, 
                  amplicon_end, 
                  search_sequence, 
                  gene_name, 
                  guide_name, 
                  return_all):
    r = requests.get('https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?job_key=' + job_id)
    text = r.text

    all_results_f = []
    for step, primer_text in enumerate(text.split('<div class="prPairInfo">\n<a name = ')[1:]):
        all_in_frame = 0
        all_results = []
        forw_primer = primer_text.split('\n')[3].replace('</td></tr>', '').split('</th><td>')[1].split('</td><td>')[0]
        rev_primer = primer_text.split('\n')[4].replace('</td></tr>', '').split('</th><td>')[1].split('</td><td>')[0]
        for template in primer_text.split('NC_')[1:]:
            template = template.split('>N')[0]
            for template_i in template.split('product length ')[1:]:
                in_frame = 0
                for string_i in template_i.split('\nTemplate        ')[1:]:
                    start_string = int(string_i.split('\n')[0].split('  ')[0])
                    end_string = int(string_i.split('\n')[0].split('  ')[2])
                    if end_string<start_string:
                        new_positions = [start_string, end_string]
                        start_string = new_positions[1]
                        end_string = new_positions[0]
                    if (('NC_' + template.split('</a>')[0] == gene_nt_id)  
                        & (start_string>=amplicon_start) & (end_string<=amplicon_end)):
                        in_frame += 1
                if in_frame<2:
                    res_1 = [np.sum([j=='.' for j in i.split('\n\n')[0].split('  ')[1][-1:]])/1 for i in template_i.split('\nTemplate        ')[1:]]
                    res_2 = [np.sum([j=='.' for j in i.split('\n\n')[0].split('  ')[1][-2:]])/2 for i in template_i.split('\nTemplate        ')[1:]]
                    res_3 = [np.sum([j=='.' for j in i.split('\n\n')[0].split('  ')[1][-3:]])/3 for i in template_i.split('\nTemplate        ')[1:]]
                    res_4 = [np.sum([j=='.' for j in i.split('\n\n')[0].split('  ')[1][-4:]])/4 for i in template_i.split('\nTemplate        ')[1:]]
                    res_5 = [np.sum([j=='.' for j in i.split('\n\n')[0].split('  ')[1][-5:]])/5 for i in template_i.split('\nTemplate        ')[1:]]
                    all_results.append([forw_primer, rev_primer, 
                                        res_1[0], res_1[1], 
                                        res_2[0], res_2[1], 
                                        res_3[0], res_3[1], 
                                        res_4[0], res_4[1], 
                                        res_5[0], res_5[1]])
                elif (in_frame==2) & (all_in_frame>0):
                    res_1 = [np.sum([j=='.' for j in i.split('\n\n')[0].split('  ')[1][-1:]])/1 for i in template_i.split('\nTemplate        ')[1:]]
                    res_2 = [np.sum([j=='.' for j in i.split('\n\n')[0].split('  ')[1][-2:]])/2 for i in template_i.split('\nTemplate        ')[1:]]
                    res_3 = [np.sum([j=='.' for j in i.split('\n\n')[0].split('  ')[1][-3:]])/3 for i in template_i.split('\nTemplate        ')[1:]]
                    res_4 = [np.sum([j=='.' for j in i.split('\n\n')[0].split('  ')[1][-4:]])/4 for i in template_i.split('\nTemplate        ')[1:]]
                    res_5 = [np.sum([j=='.' for j in i.split('\n\n')[0].split('  ')[1][-5:]])/5 for i in template_i.split('\nTemplate        ')[1:]]
                    all_results.append([forw_primer, rev_primer, 
                                        res_1[0], res_1[1], 
                                        res_2[0], res_2[1], 
                                        res_3[0], res_3[1], 
                                        res_4[0], res_4[1], 
                                        res_5[0], res_5[1]])
                else:
                    all_in_frame += 1

        if len(all_results)==0:
            all_results.append([forw_primer, rev_primer, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0])
        all_results_f.append([step+1, forw_primer ,rev_primer, 
                            pd.DataFrame(all_results)[[2, 3]].min(axis=1).max(), 
                            pd.DataFrame(all_results)[[4, 5]].min(axis=1).max(), 
                            pd.DataFrame(all_results)[[6, 7]].min(axis=1).max(),
                            pd.DataFrame(all_results)[[8, 9]].min(axis=1).max(),
                            pd.DataFrame(all_results)[[10, 11]].min(axis=1).max(),
                            np.sum(pd.DataFrame(all_results)[[2, 3]].min(axis=1)==1),
                            np.sum(pd.DataFrame(all_results)[[4, 5]].min(axis=1)==1),
                            np.sum(pd.DataFrame(all_results)[[6, 7]].min(axis=1)==1),
                            np.sum(pd.DataFrame(all_results)[[8, 9]].min(axis=1)==1),
                            np.sum(pd.DataFrame(all_results)[[10, 11]].min(axis=1)==1)])
        
    all_results_f = pd.DataFrame(all_results_f)
    if return_all:
        all_results_f = all_results_f[all_results_f[[3, 4, 5, 6, 7]].min(axis=1)<=1]
    else:
        all_results_f = all_results_f[all_results_f[[3, 4, 5, 6, 7]].min(axis=1)<1]

    all_results_f.columns = ['primer_pair', 'primer_left', 'primer_right', 
                            'score_1nt', 'score_2nt',
                            'score_3nt', 'score_4nt', 'score_5nt',
                            'bad_count_1nt', 'bad_count_2nt',
                            'bad_count_3nt', 'bad_count_4nt', 'bad_count_5nt']
    
    all_primers = []
    for step, primer_text in enumerate(text.split('<div class="prPairInfo">\n<a name = ')[1:]):
        primer_cols = primer_text.split('\n')[2].replace('</th></tr>', '').split('</th><th>')[1:]
        forw_primer = primer_text.split('\n')[3].replace('</td></tr>', '').split('</th><td>')[1].split('</td><td>')
        rev_primer = primer_text.split('\n')[4].replace('</td></tr>', '').split('</th><td>')[1].split('</td><td>')
        
        forw_primer = pd.DataFrame([forw_primer], columns=primer_cols)
        rev_primer = pd.DataFrame([rev_primer], columns=primer_cols)

        forw_primer.columns = [c + '_L' for c in forw_primer.columns]
        rev_primer.columns = [c + '_R' for c in rev_primer.columns]

        primer = pd.concat([forw_primer, rev_primer], axis=1)

        primer['pair_name'] = 'primer_' + guide_name + '_' + str(step+1)
        primer['left_name'] = primer['pair_name'].apply(lambda p: p + '_L')
        primer['right_name'] = primer['pair_name'].apply(lambda p: p + '_R')
        all_primers.append(primer)

    all_primers = pd.concat(all_primers)

    new_cut_site = len(search_sequence)//2
    all_primers['cut_site'] = new_cut_site
    all_primers['cut_site_dist_L'] = all_primers['Start_L'].apply(lambda p: new_cut_site - int(p))
    all_primers['cut_site_dist_R'] = all_primers['Start_R'].apply(lambda p: int(p) - new_cut_site)
    all_primers['amplicon_size'] = all_primers['cut_site_dist_L'] + all_primers['cut_site_dist_R']

    all_primers['COMPL_ANY'] = all_primers['pair_name'].apply(lambda p: calcHeterodimer(all_primers[all_primers['pair_name'] == p]["Sequence (5'->3')_L"][0],
                                                        all_primers[all_primers['pair_name'] == p]["Sequence (5'->3')_R"][0]).tm)
    
    new_columns = ['pair_name', 'left_name', 'right_name',
                "Sequence (5'->3')_L", "Sequence (5'->3')_R",
                'Start_L', 'Start_R', 'Stop_L', 'Stop_R', 'Length_L', 'Length_R',
                'cut_site_dist_L', 'cut_site_dist_R', 'amplicon_size',
                'Tm_L', 'Tm_R', 'GC%_L', 'GC%_R',
                'Self complementarity_L', 'Self complementarity_R',
                "Self 3' complementarity_L", "Self 3' complementarity_R",
                'COMPL_ANY']
    
    all_primers = all_primers[new_columns]
    all_primers = pd.merge(all_primers, all_results_f, left_on=["Sequence (5'->3')_L", "Sequence (5'->3')_R"],
                                            right_on=['primer_left', 'primer_right'])
    all_primers.drop(columns=['primer_left', 'primer_right'], inplace=True)

    all_primers.sort_values(by=['score_1nt', 'score_2nt', 'score_3nt', 'score_4nt', 'score_5nt'], inplace=True)

    all_primers.to_csv('src/static/outputs/' + gene_name + '/' + gene_name + ' ' + guide_name +'_primers.csv', index=None)

    return all_primers


