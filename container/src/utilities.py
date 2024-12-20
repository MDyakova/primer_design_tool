"""Functions for main.py"""

import os
import numpy as np
import pandas as pd
from primer3.bindings import calcHeterodimer
from pyensembl import EnsemblRelease
import requests
from Bio import Entrez, SeqIO

def ensemble_info(gene_name):
    """
    Get information from ensemble database:
    1. Gene sequence.
    2. Gene position on chromosome.
    3. Strand.
    4. Additional information.
    Release 109 uses human reference genome GRCh38.
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
        f"https://rest.ensembl.org/sequence/region/human/"
        f"{chromosome_name}:{start_gene}-{end_gene}"
        "?content-type=application/json;version=109"
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
    """
    Find seqrch sequence coordinates and ncbi id.
    This data is used to looking for off-targets of primers
    in function blast_results.    
    """

    Entrez.email = "Your.Name.Here@example.org"
    handle = Entrez.efetch(db="gene", id=ncbi_id, rettype="gb", retmode="text")
    information = handle.read()
    gene_nt_id = (
        "NC_" + information.split("Annotation: ")[1].split("NC_")[1].split(" ")[0]
    )
    start_gene_nt, end_gene_nt = (
        information.split("Annotation: ")[1]
        .split("(")[1]
        .split(")")[0]
        .split(",")[0]
        .split("..")
    )

    start_gene_nt = int(start_gene_nt) - 1500
    end_gene_nt = int(end_gene_nt) + 1500

    if strand == "+":
        strand_nt = 1
    else:
        strand_nt = 2

    Entrez.email = "Your.Name.Here@example.org"
    handle = Entrez.efetch(
        db="nucleotide",
        id=gene_nt_id,
        rettype="gb",
        retmode="text",
        seq_start=start_gene_nt,
        seq_stop=end_gene_nt,
        strand=strand_nt,
    )
    seq_record = [seq_record for seq_record in SeqIO.parse(handle, "gb")][0]
    refseq_sequence = str(seq_record.seq)

    if strand == "+":
        amplicon_start = len(refseq_sequence.split(search_sequence)[0]) + int(
            start_gene_nt
        )
        amplicon_end = amplicon_start + len(search_sequence)
    else:
        amplicon_end = end_gene_nt - len(refseq_sequence.split(search_sequence)[0])
        amplicon_start = amplicon_end - len(search_sequence)

    return amplicon_start, amplicon_end, gene_nt_id


def guide_info(guide_seq, ensemble_gene_seq):
    """
    Search guide position in ensemble sequence.
    If we have two guides (TGEE) it should be separated by ;.
    If we have one guide (Cas9) it is single sequence.
    Output is sequence between two guides (included guides sequences) 
    on one strand.
    """

    compl_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    if ";" not in guide_seq:
        if len(guide_seq) > 23:
            right_guide = guide_seq[-20:]
            left_guide = "".join([compl_dict[i] for i in guide_seq[:20]][::-1])
            guide_seq = ";".join([left_guide, right_guide])
    if ";" not in guide_seq:
        if guide_seq not in ensemble_gene_seq:
            guide = "".join([compl_dict[i] for i in guide_seq][::-1])
        else:
            guide = guide_seq
    else:
        left_guide, right_guide = guide_seq.split(";")
        left_guide = left_guide.strip()
        right_guide = right_guide.strip()

        if left_guide not in ensemble_gene_seq:
            left_guide = "".join([compl_dict[i] for i in left_guide][::-1])
        if right_guide not in ensemble_gene_seq:
            right_guide = "".join([compl_dict[i] for i in right_guide][::-1])

        left_guide_start = len(ensemble_gene_seq.split(left_guide)[0])
        right_guide_start = len(ensemble_gene_seq.split(right_guide)[0])

        if left_guide_start > right_guide_start:
            guide = (
                right_guide
                + ensemble_gene_seq.split(right_guide)[1].split(left_guide)[0]
                + left_guide
            )
        else:
            guide = (
                left_guide
                + ensemble_gene_seq.split(left_guide)[1].split(right_guide)[0]
                + right_guide
            )

    return guide


def blast_primers(
    search_sequence,
    primer5_start,
    primer5_end,
    primer3_start,
    primer3_end,
    product_size_min,
    product_size_max,
):
    """
    Make a lint to Blast Primers site with nessesary parameters.
    1. search_sequence - sequence where we find primers.
    2. primer5_start, primer5_end, primer3_start, primer3_end - possible coordinates for primers in search sequence.
    """

    url = (
        f"https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?LINK_LOC=bookmark"
        f"&INPUT_SEQUENCE={search_sequence}"
        f"&PRIMER5_START={primer5_start}"
        f"&PRIMER5_END={primer5_end}"
        f"&PRIMER3_START={primer3_start}"
        f"&PRIMER3_END={primer3_end}"
        f"&OVERLAP_5END=7"
        f"&OVERLAP_3END=4"
        f"&PRIMER_PRODUCT_MIN={product_size_min}"
        f"&PRIMER_PRODUCT_MAX={product_size_max}"
        f"&PRIMER_NUM_RETURN=1000"
        f"&PRIMER_MIN_TM=55.0"
        f"&PRIMER_OPT_TM=60.0"
        f"&PRIMER_MAX_TM=65.0"
        f"&PRIMER_MAX_DIFF_TM=3"
        f"&PRIMER_ON_SPLICE_SITE=0"
        f"&SEARCHMODE=0"
        f"&SPLICE_SITE_OVERLAP_5END=7"
        f"&SPLICE_SITE_OVERLAP_3END=4"
        f"&SPLICE_SITE_OVERLAP_3END_MAX=8"
        f"&SPAN_INTRON=off"
        f"&MIN_INTRON_SIZE=1000"
        f"&MAX_INTRON_SIZE=1000000"
        f"&SEARCH_SPECIFIC_PRIMER=on"
        f"&EXCLUDE_ENV=on"
        f"&EXCLUDE_XM=on"
        f"&TH_OLOGO_ALIGNMENT=off"
        f"&TH_TEMPLATE_ALIGNMENT=off"
        f"&ORGANISM=Homo%20sapiens"
        f"&PRIMER_SPECIFICITY_DATABASE=refseq_representative_genomes"
        f"&TOTAL_PRIMER_SPECIFICITY_MISMATCH=1"
        f"&PRIMER_3END_SPECIFICITY_MISMATCH=1"
        f"&MISMATCH_REGION_LENGTH=5"
        f"&TOTAL_MISMATCH_IGNORE=6"
        f"&MAX_TARGET_SIZE=4000"
        f"&ALLOW_TRANSCRIPT_VARIANTS=off"
        f"&HITSIZE=50000"
        f"&EVALUE=30000"
        f"&WORD_SIZE=7"
        f"&MAX_CANDIDATE_PRIMER=500"
        f"&PRIMER_MIN_SIZE=16"
        f"&PRIMER_OPT_SIZE=20"
        f"&PRIMER_MAX_SIZE=25"
        f"&PRIMER_MIN_GC=20.0"
        f"&PRIMER_MAX_GC=80.0"
        f"&GC_CLAMP=0"
        f"&NUM_TARGETS_WITH_PRIMERS=1000"
        f"&NUM_TARGETS=20"
        f"&MAX_TARGET_PER_TEMPLATE=100"
        f"&POLYX=3"
        f"&SELF_ANY=8.00"
        f"&SELF_END=3.00"
        f"&PRIMER_MAX_END_STABILITY=9"
        f"&PRIMER_MAX_END_GC=5"
        f"&PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00"
        f"&PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00"
        f"&PRIMER_MAX_SELF_ANY_TH=45.0"
        f"&PRIMER_MAX_SELF_END_TH=35.0"
        f"&PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0"
        f"&PRIMER_PAIR_MAX_COMPL_END_TH=35.0"
        f"&PRIMER_MAX_HAIRPIN_TH=24.0"
        f"&PRIMER_MAX_TEMPLATE_MISPRIMING=12.00"
        f"&PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00"
        f"&PRIMER_PAIR_MAX_COMPL_ANY=8.00"
        f"&PRIMER_PAIR_MAX_COMPL_END=3.00"
        f"&PRIMER_MISPRIMING_LIBRARY=repeat/repeat_9606"
        f"&NO_SNP=on"
        f"&LOW_COMPLEXITY_FILTER=on"
        f"&MONO_CATIONS=50.0"
        f"&DIVA_CATIONS=1.5"
        f"&CON_ANEAL_OLIGO=300.0"
        f"&CON_DNTPS=0.6"
        f"&SALT_FORMULAR=1"
        f"&TM_METHOD=1"
        f"&PRIMER_INTERNAL_OLIGO_MIN_SIZE=18"
        f"&PRIMER_INTERNAL_OLIGO_OPT_SIZE=20"
        f"&PRIMER_INTERNAL_OLIGO_MAX_SIZE=27"
        f"&PRIMER_INTERNAL_OLIGO_MIN_TM=57.0"
        f"&PRIMER_INTERNAL_OLIGO_OPT_TM=60.0"
        f"&PRIMER_INTERNAL_OLIGO_MAX_TM=63.0"
        f"&PRIMER_INTERNAL_OLIGO_MAX_GC=80.0"
        f"&PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT=50"
        f"&PRIMER_INTERNAL_OLIGO_MIN_GC=20.0"
        f"&PICK_HYB_PROBE=off"
        f"&NEWWIN=on"
        f"&NEWWIN=on"
        f"&SHOW_SVIEWER=true"
    )

    return url


def blast_results(
    job_id,
    gene_nt_id,
    amplicon_start,
    amplicon_end,
    search_sequence,
    gene_name,
    guide_name,
    return_all,
):
    """
    1. Get results from Blast primers site and find off-targets.
    2. Parsing all results found. 
    3. If off-target is inside the search sequence it to be on-target. 
    4. For each pair of primers we count number of off-targets:
        a. score_1nt - proportion of complimentary nucleotides for last one nucleotide. 
           bad_count_1nt - number of off-targets which have complimentary last one nucleotide.
           If score_1nt = 0 it means that all off-targets for this pair of primers 
           have non-complimentary last nucleotides to off-target sequence.
           If score_1nt = 1 it means that at least one off-target for this pair of primers 
           has complimentary last nucleotide to off-target sequence.
        b. score_2nt - proportion of complimentary nucleotides for last two nucleotides. 
           bad_count_2nt - number of off-targets which have complimentary last two nucleotides.
           If score_2nt = 0 it means that all off-targets for this pair of primers 
           have non-complimentary last two nucleotides to off-target sequence.
           If score_2nt = 1 it means that at least one off-target for this pair of primers 
           has complimentary last two nucleotides to off-target sequence. 
           If score_2nt = 0.5 it means that at least one off-target for this pair of primers 
           has complimentary one of two last nucleotides to off-target sequence. 
        c. score_3nt, score_4nt, score_5nt - similar to score_2nt for last three, 
           four and five nucleotides respectively.
    5. Calculate termodynamics parameters for pair of primers.
    6. Calculate distances to cut site or middle of sequence (if we don't have guide sequence).
    7. Calculate amplicon size for each pair of primers.
    """

    request_results = requests.get(
        "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?job_key="
        + job_id
    )
    text = request_results.text

    all_results_f = []
    for step, primer_text in enumerate(
        text.split('<div class="prPairInfo">\n<a name = ')[1:]
    ):
        all_in_frame = 0
        all_results = []
        forw_primer = (
            primer_text.split("\n")[3]
            .replace("</td></tr>", "")
            .split("</th><td>")[1]
            .split("</td><td>")[0]
        )
        rev_primer = (
            primer_text.split("\n")[4]
            .replace("</td></tr>", "")
            .split("</th><td>")[1]
            .split("</td><td>")[0]
        )
        for template in primer_text.split("NC_")[1:]:
            template = template.split(">N")[0]
            for template_i in template.split("product length ")[1:]:
                in_frame = 0
                for string_i in template_i.split("\nTemplate        ")[1:]:
                    start_string = int(string_i.split("\n")[0].split("  ")[0])
                    end_string = int(string_i.split("\n")[0].split("  ")[2])
                    if end_string < start_string:
                        new_positions = [start_string, end_string]
                        start_string = new_positions[1]
                        end_string = new_positions[0]
                    if (
                        ("NC_" + template.split("</a>")[0] == gene_nt_id)
                        & (start_string >= amplicon_start)
                        & (end_string <= amplicon_end)
                    ):
                        in_frame += 1
                if in_frame < 2:
                    res_1 = [
                        np.sum(
                            [j == "." for j in i.split("\n\n")[0].split("  ")[1][-1:]]
                        )
                        / 1
                        for i in template_i.split("\nTemplate        ")[1:]
                    ]
                    res_2 = [
                        np.sum(
                            [j == "." for j in i.split("\n\n")[0].split("  ")[1][-2:]]
                        )
                        / 2
                        for i in template_i.split("\nTemplate        ")[1:]
                    ]
                    res_3 = [
                        np.sum(
                            [j == "." for j in i.split("\n\n")[0].split("  ")[1][-3:]]
                        )
                        / 3
                        for i in template_i.split("\nTemplate        ")[1:]
                    ]
                    res_4 = [
                        np.sum(
                            [j == "." for j in i.split("\n\n")[0].split("  ")[1][-4:]]
                        )
                        / 4
                        for i in template_i.split("\nTemplate        ")[1:]
                    ]
                    res_5 = [
                        np.sum(
                            [j == "." for j in i.split("\n\n")[0].split("  ")[1][-5:]]
                        )
                        / 5
                        for i in template_i.split("\nTemplate        ")[1:]
                    ]
                    all_results.append(
                        [
                            forw_primer,
                            rev_primer,
                            res_1[0],
                            res_1[1],
                            res_2[0],
                            res_2[1],
                            res_3[0],
                            res_3[1],
                            res_4[0],
                            res_4[1],
                            res_5[0],
                            res_5[1],
                        ]
                    )
                elif (in_frame == 2) & (all_in_frame > 0):
                    res_1 = [
                        np.sum(
                            [j == "." for j in i.split("\n\n")[0].split("  ")[1][-1:]]
                        )
                        / 1
                        for i in template_i.split("\nTemplate        ")[1:]
                    ]
                    res_2 = [
                        np.sum(
                            [j == "." for j in i.split("\n\n")[0].split("  ")[1][-2:]]
                        )
                        / 2
                        for i in template_i.split("\nTemplate        ")[1:]
                    ]
                    res_3 = [
                        np.sum(
                            [j == "." for j in i.split("\n\n")[0].split("  ")[1][-3:]]
                        )
                        / 3
                        for i in template_i.split("\nTemplate        ")[1:]
                    ]
                    res_4 = [
                        np.sum(
                            [j == "." for j in i.split("\n\n")[0].split("  ")[1][-4:]]
                        )
                        / 4
                        for i in template_i.split("\nTemplate        ")[1:]
                    ]
                    res_5 = [
                        np.sum(
                            [j == "." for j in i.split("\n\n")[0].split("  ")[1][-5:]]
                        )
                        / 5
                        for i in template_i.split("\nTemplate        ")[1:]
                    ]
                    all_results.append(
                        [
                            forw_primer,
                            rev_primer,
                            res_1[0],
                            res_1[1],
                            res_2[0],
                            res_2[1],
                            res_3[0],
                            res_3[1],
                            res_4[0],
                            res_4[1],
                            res_5[0],
                            res_5[1],
                        ]
                    )
                else:
                    all_in_frame += 1

        if len(all_results) == 0:
            all_results.append([forw_primer, rev_primer, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        all_results_f.append(
            [
                step + 1,
                forw_primer,
                rev_primer,
                pd.DataFrame(all_results)[[2, 3]].min(axis=1).max(),
                pd.DataFrame(all_results)[[4, 5]].min(axis=1).max(),
                pd.DataFrame(all_results)[[6, 7]].min(axis=1).max(),
                pd.DataFrame(all_results)[[8, 9]].min(axis=1).max(),
                pd.DataFrame(all_results)[[10, 11]].min(axis=1).max(),
                np.sum(pd.DataFrame(all_results)[[2, 3]].min(axis=1) == 1),
                np.sum(pd.DataFrame(all_results)[[4, 5]].min(axis=1) == 1),
                np.sum(pd.DataFrame(all_results)[[6, 7]].min(axis=1) == 1),
                np.sum(pd.DataFrame(all_results)[[8, 9]].min(axis=1) == 1),
                np.sum(pd.DataFrame(all_results)[[10, 11]].min(axis=1) == 1),
            ]
        )

    all_results_f = pd.DataFrame(all_results_f)
    if return_all:
        all_results_f = all_results_f[all_results_f[[3, 4, 5, 6, 7]].min(axis=1) <= 1]
    else:
        all_results_f = all_results_f[all_results_f[[3, 4, 5, 6, 7]].min(axis=1) < 1]

    all_results_f.columns = [
        "primer_pair",
        "primer_left",
        "primer_right",
        "score_1nt",
        "score_2nt",
        "score_3nt",
        "score_4nt",
        "score_5nt",
        "bad_count_1nt",
        "bad_count_2nt",
        "bad_count_3nt",
        "bad_count_4nt",
        "bad_count_5nt",
    ]

    all_primers = []
    for step, primer_text in enumerate(
        text.split('<div class="prPairInfo">\n<a name = ')[1:]
    ):
        primer_cols = (
            primer_text.split("\n")[2].replace("</th></tr>", "").split("</th><th>")[1:]
        )
        forw_primer = (
            primer_text.split("\n")[3]
            .replace("</td></tr>", "")
            .split("</th><td>")[1]
            .split("</td><td>")
        )
        rev_primer = (
            primer_text.split("\n")[4]
            .replace("</td></tr>", "")
            .split("</th><td>")[1]
            .split("</td><td>")
        )

        forw_primer = pd.DataFrame([forw_primer], columns=primer_cols)
        rev_primer = pd.DataFrame([rev_primer], columns=primer_cols)

        forw_primer.columns = [c + "_L" for c in forw_primer.columns]
        rev_primer.columns = [c + "_R" for c in rev_primer.columns]

        primer = pd.concat([forw_primer, rev_primer], axis=1)

        primer["pair_name"] = "primer_" + guide_name + "_" + str(step + 1)
        primer["left_name"] = primer["pair_name"].apply(lambda p: p + "_L")
        primer["right_name"] = primer["pair_name"].apply(lambda p: p + "_R")
        all_primers.append(primer)

    all_primers = pd.concat(all_primers)

    new_cut_site = len(search_sequence) // 2
    all_primers["cut_site"] = new_cut_site
    all_primers["cut_site_dist_L"] = all_primers["Start_L"].apply(
        lambda p: new_cut_site - int(p)
    )
    all_primers["cut_site_dist_R"] = all_primers["Start_R"].apply(
        lambda p: int(p) - new_cut_site
    )
    all_primers["amplicon_size"] = (
        all_primers["cut_site_dist_L"] + all_primers["cut_site_dist_R"]
    )

    all_primers["COMPL_ANY"] = all_primers["pair_name"].apply(
        lambda p: calcHeterodimer(
            all_primers[all_primers["pair_name"] == p]["Sequence (5'->3')_L"][0],
            all_primers[all_primers["pair_name"] == p]["Sequence (5'->3')_R"][0],
        ).tm
    )

    new_columns = [
        "pair_name",
        "left_name",
        "right_name",
        "Sequence (5'->3')_L",
        "Sequence (5'->3')_R",
        "Start_L",
        "Start_R",
        "Stop_L",
        "Stop_R",
        "Length_L",
        "Length_R",
        "cut_site_dist_L",
        "cut_site_dist_R",
        "amplicon_size",
        "Tm_L",
        "Tm_R",
        "GC%_L",
        "GC%_R",
        "Self complementarity_L",
        "Self complementarity_R",
        "Self 3' complementarity_L",
        "Self 3' complementarity_R",
        "COMPL_ANY",
    ]

    all_primers = all_primers[new_columns]
    all_primers = pd.merge(
        all_primers,
        all_results_f,
        left_on=["Sequence (5'->3')_L", "Sequence (5'->3')_R"],
        right_on=["primer_left", "primer_right"],
    )
    all_primers.drop(columns=["primer_left", "primer_right"], inplace=True)

    all_primers.sort_values(
        by=[
            "bad_count_1nt",
            "bad_count_2nt",
            "bad_count_3nt",
            "bad_count_4nt",
            "bad_count_5nt",
        ],
        inplace=True,
    )

    files_path = os.path.join(
        "src",
        "static",
        "outputs",
        gene_name,
        gene_name + " " + guide_name + "_primers.csv",
    )
    all_primers.to_csv(files_path, index=None)

    return all_primers

"""
Next templates for gene bank file.
This file allow to automatically make SnapGene file with all features and primers
and open it in SnapGene Viewer without pair SnapGene program. 
"""

TITLE = """LOCUS       {gene_name}        {len_full_sequence} bp DNA     linear   UNA {date_today}
DEFINITION  {gene_name}.
ACCESSION   .
VERSION     .
KEYWORDS    .
SOURCE      synthetic DNA construct
  ORGANISM  synthetic DNA construct
REFERENCE   1  (bases 1 to {len_full_sequence})
  AUTHORS   .
  TITLE     Direct Submission
  JOURNAL   For SnapGene Viewer
            https://www.snapgene.com
FEATURES             Location/Qualifiers"""

FEATURE_SOURCE = '''     source          1..{len_full_sequence}
                     /mol_type="other DNA"
                     /note="color: #ffffff"
                     /organism="synthetic DNA construct"'''

ORIGIN = """ORIGIN
{origin_seq}
//"""


def misc_feature_template(start, end, label, color, direction):
    """Make a template for snapgene features"""

    if direction == "None":
        misc_f = f'''     misc_feature    {start}..{end}
                        /label={label}
                        /note="color: {color}"'''
    else:
        misc_f = f'''     misc_feature    {start}..{end}
                        /label={label}
                        /note="color: {color}; direction: {direction}"'''
    return misc_f


def primer_template(name, seq, date_today, start, end):
    """Make a template for snapgene primers"""
    primer_f = f'''     primer_bind     complement({start}..{end})
                     /label={name}
                     /note="color: black; sequence: 
                     {seq}; added: 
                     {date_today}"'''
    return primer_f


def gene_bank_file(
    gene_name,
    full_sequence,
    date_today,
    elements_list,
    files_name,
    oligos=None,
    title=TITLE,
    feature_sourse=FEATURE_SOURCE,
    origin=ORIGIN,
):
    """Make gene bank file for snapgene tool"""

    if oligos is None:
        oligos = []

    title = title.format(
        gene_name=gene_name,
        full_sequence=full_sequence,
        date_today=date_today,
        len_full_sequence=len(full_sequence),
    )

    feature_sourse = feature_sourse.format(len_full_sequence=len(full_sequence))

    all_misc_feature = ""
    for feature in elements_list:
        start = feature[1]
        end = feature[2]
        name = feature[0]
        color = feature[4]
        direction = feature[3]
        if direction == "+":
            direction = "RIGHT"
        elif direction == "-":
            direction = "LEFT"
        else:
            direction = "None"
        misc_feature = misc_feature_template(start, end, name, color, direction)
        all_misc_feature += misc_feature + "\n"

    all_primers = ""
    for oligo in oligos:
        start = oligo[2]
        end = oligo[3]
        name = oligo[0]
        seq = oligo[1]

        primer_feature = primer_template(name, seq, date_today, start, end)
        all_primers += primer_feature + "\n"

    
    """
    This part make a sequence in Gene Bank format.
    """
    origin_seq = ""
    for i in range(len(full_sequence)):
        if i % 60 == 0:
            origin_seq += "\n" + " " * (9 - len(str(i + 1))) + str(i + 1) + " "
        elif i % 10 == 0:
            origin_seq += " "
        origin_seq += full_sequence[i].lower()
    origin_seq = origin_seq[1:]

    origin = origin.format(origin_seq=origin_seq)

    gbk_file_name = os.path.join(
        "src", "static", "outputs", gene_name, files_name + ".gbk"
    )

    with open(gbk_file_name, "w", encoding="utf-8") as file:
        file.write(title + "\n")
        file.write(feature_sourse + "\n")
        file.write(all_misc_feature)
        file.write(all_primers)
        file.write(origin + "\n")


def find_elements(sequence, guide_seq, guide_name, selected_primers, is_guide=True):
    """
    Make lists with guide and primer positions for SnapGene file.
    If we don't have guide, we don't have features in SnapGene file.
    """

    if is_guide:
        start_guide = len(sequence.split(guide_seq)[0]) + 1
        end_guide = start_guide + len(guide_seq) - 1

        elements_list = [[guide_name, start_guide, end_guide, "None", "#0000EE"]]
    else:
        elements_list = []

    oligos = []
    for name, seq in zip(
        selected_primers["left_name"], selected_primers["Sequence (5'->3')_L"]
    ):
        oligos.append([name, seq, 1, len(seq)])
    for name, seq in zip(
        selected_primers["right_name"], selected_primers["Sequence (5'->3')_R"]
    ):
        oligos.append([name, seq, 1, len(seq)])

    return elements_list, oligos


def primers_pivot_table(
    selected_primers,
    gene_name,
    guide_name,
    guide_full_seq,
    min_dist,
    max_dist,
    min_size,
    max_size,
    insert_seq,
):
    """
    Make a pivot table with amplicon sizes for all possible primers combinations.
    This data used if we want to check all selected primers for one guide or search sequence between each other. 
    In table are pairs which satisfy the conditions of amplicon range and dictances range. 
    """

    all_primers_dist = []

    for name_l, dist_l in zip(
        selected_primers["left_name"], selected_primers["cut_site_dist_L"]
    ):
        for name_r, dist_r in zip(
            selected_primers["right_name"], selected_primers["cut_site_dist_R"]
        ):
            if insert_seq == "":
                if (
                    (dist_l >= min_dist)
                    & (dist_l <= max_dist)
                    & (dist_r >= min_dist)
                    & (dist_r <= max_dist)
                    & ((dist_l + dist_r) >= min_size)
                    & ((dist_l + dist_r) <= max_size)
                ):
                    all_primers_dist.append(
                        [guide_name, name_l, name_r, dist_l + dist_r]
                    )
            elif guide_full_seq != "":
                if (
                    (dist_l >= min_dist)
                    & (dist_l <= max_dist)
                    & (dist_r >= min_dist)
                    & (dist_r <= max_dist)
                    & ((dist_l + dist_r) >= min_size)
                    & ((dist_l + dist_r) <= max_size)
                ):
                    all_primers_dist.append(
                        [guide_name, name_l, name_r, dist_l + dist_r]
                    )
            else:
                if ((dist_l + dist_r) >= min_size) & ((dist_l + dist_r) <= max_size):
                    all_primers_dist.append(
                        [guide_name, name_l, name_r, dist_l + dist_r]
                    )

    all_primers_dist = pd.DataFrame(
        all_primers_dist, columns=("guide_name", "primer_L", "primer_R", "dist")
    )
    all_primers_dist = all_primers_dist.pivot_table(
        "dist", index=["primer_L", "primer_R"], columns=["guide_name"], aggfunc="max"
    )

    file_path = os.path.join(
        "src", "static", "outputs", gene_name, guide_name + "_distances.csv"
    )
    all_primers_dist.to_csv(file_path)


def primers_pivot_table_few_guides(
    primers_table, min_dist, max_dist, min_size, max_size, save_name
):
    """
    Make a pivot table with amplicon sizes for all possible primers combinations for all guides in table.
    This data used if we want to check all selected primers for few guides or search sequence between each other. 
    In table are pairs which satisfy the conditions of amplicon range and dictances range. 
    
    """

    all_primers_dist = []
    all_primers_ampl = []

    for guide_name in pd.unique(primers_table["guide_name"]):
        cut_site = np.max(
            primers_table[primers_table["guide_name"] == guide_name]["cut_site_coords"]
        )
        for name_l, coord_l in zip(
            primers_table["left_name"], primers_table["coords_left"]
        ):
            for name_r, coord_r in zip(
                primers_table["right_name"], primers_table["coords_right"]
            ):
                dist_l = int(cut_site - coord_l)
                dist_r = int(coord_r - cut_site)
                if (
                    (dist_l >= min_dist)
                    & (dist_l <= max_dist)
                    & (dist_r >= min_dist)
                    & (dist_r <= max_dist)
                    & ((dist_l + dist_r) >= min_size)
                    & ((dist_l + dist_r) <= max_size)
                ):
                    all_primers_dist.append(
                        [
                            guide_name,
                            name_l,
                            name_r,
                            "/".join([str(dist_l), str(dist_r)]),
                        ]
                    )
                    all_primers_ampl.append(
                        [guide_name, name_l, name_r, dist_l + dist_r]
                    )

    all_primers_dist = pd.DataFrame(
        all_primers_dist, columns=("guide_name", "primer_L", "primer_R", "dist")
    )
    all_primers_dist = all_primers_dist.pivot_table(
        "dist", index=["primer_L", "primer_R"], columns=["guide_name"], aggfunc="max"
    )

    all_primers_ampl = pd.DataFrame(
        all_primers_ampl,
        columns=("guide_name", "primer_L", "primer_R", "amplicon_size"),
    )
    all_primers_ampl = all_primers_ampl.groupby(
        by=["primer_L", "primer_R"], as_index=False
    ).max()

    all_primers_dist = pd.merge(
        all_primers_dist,
        all_primers_ampl[["primer_L", "primer_R", "amplicon_size"]],
        on=["primer_L", "primer_R"],
    )

    file_path = os.path.join("src", "static", "outputs", save_name + "_distances.csv")
    all_primers_dist.to_csv(file_path, index=False)


def primers_coords(ensemble_gene_seq, selected_primers):
    """Find primers coordinates in ensemble gene sequence"""

    compl_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    oligos_left = []
    for name, seq in zip(
        selected_primers["left_name"], selected_primers["Sequence (5'->3')_L"]
    ):
        start_primer = len(ensemble_gene_seq.split(seq)[0]) + 1
        oligos_left.append([name, seq, start_primer])

    oligos_right = []
    for name, seq in zip(
        selected_primers["right_name"], selected_primers["Sequence (5'->3')_R"]
    ):
        seq_rev = "".join([compl_dict[n] for n in seq][::-1])
        start_primer = len(ensemble_gene_seq.split(seq_rev)[0]) + 1
        end_primer = start_primer + len(seq_rev)
        oligos_right.append([name, seq, end_primer])

    oligos_left = pd.DataFrame(
        oligos_left, columns=("left_name", "Sequence (5'->3')_L", "coords_left")
    )
    oligos_right = pd.DataFrame(
        oligos_right, columns=("right_name", "Sequence (5'->3')_R", "coords_right")
    )

    return oligos_left, oligos_right
