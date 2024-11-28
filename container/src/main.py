"""
Code launch FLASK server to find primers for gene editing tasks
"""
import os
from io import BytesIO
import zipfile
import shutil
from datetime import date
import pandas as pd
import json

from flask import Flask, render_template, request, send_file
from flask_wtf import FlaskForm
from wtforms import StringField, TextAreaField

from utilities import (
    ensemble_info,
    ncbi_info,
    guide_info,
    blast_primers,
    blast_results,
    gene_bank_file,
    find_elements,
    primers_pivot_table,
    primers_coords,
    primers_pivot_table_few_guides,
)

app = Flask(__name__)
app.config["SECRET_KEY"] = "mysecretkey"  # fake key to work with flask server

DATE_TODAY = str(date.today())


# dictionary with all useful variables necessary to save throw whole pipeline

with open('config.json', 'r') as f:
    config = json.load(f)

out_dict_outer = config['initial_values']

# This dict used to clean all forms. 
out_dict_start = out_dict_outer.copy()

out_keys_clean = []

compl_dict = config['compl_dict']


class GeneInfo(FlaskForm):
    """
    Information about ensemble and ncbi gene names or id
    """

    text_field = StringField("Ensemble gene name", default="")
    text_field2 = StringField("  NCBI gene id", default="")
    text_field3 = StringField(
        "Guide sequence", default="", render_kw={"style": "width: 400px;"}
    )
    text_field3_2 = StringField(
        "Guide name or task name", default="", render_kw={"style": "width: 100px;"}
    )
    text_field4 = StringField(
        "Minimal distance", default="", render_kw={"style": "width: 50px;"}
    )
    text_field5 = StringField(
        "Maximal distance", default="", render_kw={"style": "width: 50px;"}
    )
    text_field6 = StringField(
        "Product size min", default="", render_kw={"style": "width: 50px;"}
    )
    text_field7 = StringField(
        "Product size max", default="", render_kw={"style": "width: 50px;"}
    )
    text_field8 = TextAreaField(
        "", default="", render_kw={"style": "width: 700px; height: 100px"}
    )


class BlastInfo(FlaskForm):
    """
    Blast results id
    """

    text_field9 = StringField(
        "Blast results id", default="", render_kw={"style": "width: 400px;"}
    )


def index(out_dict):
    """
    Launch main code to prepare input data for golden gate primers
    """

    gene_info_form = GeneInfo()
    blast_info_form = BlastInfo()

    forms = {"gene_info_form": gene_info_form, "blast_info_form": blast_info_form}

    short_links = {"Blast": out_dict["blast_url"]}

    out_dict["gene_dict"] = ""

    # Show all primers
    checkbox_all_cond = False

    if request.method == "POST":
        if "gene_info_form_submit" in request.form:
            # Data from ensemble and ncbi, guide sequence
            gene_name = gene_info_form.text_field.data
            gene_name = gene_name.upper()
            ncbi_id = gene_info_form.text_field2.data
            ncbi_id = ncbi_id.upper()
            guide_seq = gene_info_form.text_field3.data
            guide_seq = guide_seq.upper()
            guide_name = gene_info_form.text_field3_2.data

            # Parameters for primers location around cut site
            min_dist = gene_info_form.text_field4.data
            max_dist = gene_info_form.text_field5.data
            min_size = gene_info_form.text_field6.data
            max_size = gene_info_form.text_field7.data

            out_dict["guide_seq"] = guide_seq
            out_dict["guide_name"] = guide_name

            try:
                out_dict["min_dist"] = int(min_dist)
                out_dict["max_dist"] = int(max_dist)
                out_dict["min_size"] = int(min_size)
                out_dict["max_size"] = int(max_size)
            except Exception:
                text_error = "primer position parameters must be > 0"
                out_dict["gene_dict"] = (
                    "<span class='red-text'>" + "Error: " + str(text_error) + "</span>"
                )
                return render_template(
                    "home.html",
                    out_dict=out_dict,
                    forms=forms,
                    links=short_links,
                    tables=[],
                    checked=False,
                )

            # Insert sequence (if we don't have guide sequence or gene information)
            insert_seq = gene_info_form.text_field8.data
            out_dict["insert_seq"] = insert_seq
            insert_seq = insert_seq.upper()

            # remove all style elements from insert sequence
            if insert_seq != "":
                if ">" in insert_seq:
                    gene_name = insert_seq.split("\n")[0].replace(">", "").strip()
                    out_dict["gene_name"] = gene_name
                    insert_seq = "".join(insert_seq.split("\n")[1:])
                insert_seq = (
                    insert_seq.replace("\n", "")
                    .replace("\r", "")
                    .replace(" ", "")
                    .strip()
                )
                out_dict["insert_seq"] = insert_seq

                if guide_seq == "":
                    # If we have only search sequence
                    primer5_start = ""
                    primer5_end = ""
                    primer3_start = ""
                    primer3_end = ""
                    out_dict["search_sequence"] = insert_seq

                else:
                    """
                    If we have search sequence and guide sequence.
                    Find cut site position and coordinates of left and right possible primer
                    positions in search_sequence.
                    """
                    out_dict["strand"] = "+"

                    try:
                        guide_full_seq = guide_info(
                            out_dict["guide_seq"],
                            out_dict["insert_seq"],
                        )

                        out_dict["guide_full_seq"] = guide_full_seq

                        left_flank = out_dict["insert_seq"].split(
                            out_dict["guide_full_seq"]
                        )[0][-out_dict["max_dist"] :]
                        right_flank = out_dict["insert_seq"].split(
                            out_dict["guide_full_seq"]
                        )[1][: out_dict["max_dist"]]

                        search_sequence = (
                            left_flank + out_dict["guide_full_seq"] + right_flank
                        )
                        out_dict["search_sequence"] = search_sequence
                        cut_size = (
                            len(search_sequence.split(out_dict["guide_full_seq"])[0])
                            + len(out_dict["guide_full_seq"]) // 2
                        )

                        primer5_start = 1
                        primer5_end = cut_size - out_dict["min_dist"]
                        primer3_start = cut_size + out_dict["min_dist"]
                        primer3_end = len(search_sequence)

                    except Exception:
                        text_error = "check guide sequence"
                        out_dict["gene_dict"] = (
                            "<span class='red-text'>"
                            + "Error: "
                            + str(text_error)
                            + "</span>"
                        )
                        return render_template(
                            "home.html",
                            out_dict=out_dict,
                            forms=forms,
                            links=short_links,
                            tables=[],
                            checked=False,
                        )

                """
                If we don't have gene information we can't distinguish on-target and off-taget.
                In this case all blast hits are off-targets.
                """

                out_dict["amplicon_start"] = -1
                out_dict["amplicon_end"] = -1
                out_dict["gene_nt_id"] = ""
                out_dict["guide_name"] = guide_name
                checkbox_all_cond = True

                # Create output directory or clean if directory is exist.
                output_directory = os.path.join("src", "static", "outputs", gene_name)
                try:
                    shutil.rmtree(output_directory)
                except Exception:
                    pass
                os.makedirs(output_directory, exist_ok=True)

            if insert_seq == "":
                #If we have gene information and guide sequence
                if (gene_name == "") | (ncbi_id == "") | (guide_seq == ""):
                    text_error = "If you want to search new primers, enter all the data"
                    out_dict["gene_dict"] = (
                        "<span class='red-text'>"
                        + "Error: "
                        + str(text_error)
                        + "</span>"
                    )
                    return render_template(
                        "home.html",
                        out_dict=out_dict,
                        forms=forms,
                        links=short_links,
                        tables=[],
                        checked=False,
                    )

                out_dict["gene_name"] = gene_name
                out_dict["ncbi_id"] = ncbi_id
                out_dict["guide_seq"] = guide_seq
                out_dict["guide_name"] = guide_name

                # Create output directory
                output_directory = os.path.join("src", "static", "outputs", gene_name)

                try:
                    shutil.rmtree(output_directory)
                except Exception:
                    pass
                os.makedirs(output_directory, exist_ok=True)

                try:
                    # Get gene information from ensemble database
                    ensemble_gene_seq, gene_dict, strand = ensemble_info(gene_name)
                    out_dict["ensemble_gene_seq"] = ensemble_gene_seq
                    out_dict["strand"] = strand
                    out_dict["gene_dict"] = "Gene info: " + str(gene_dict)
                except Exception:
                    text_error = "check gene name"
                    out_dict["gene_dict"] = (
                        "<span class='red-text'>"
                        + "Error: "
                        + str(text_error)
                        + "</span>"
                    )
                    return render_template(
                        "home.html",
                        out_dict=out_dict,
                        forms=forms,
                        links=short_links,
                        tables=[],
                        checked=False,
                    )

                try:
                    """ Get guide full sequence on one gene strand.
                    It included both guides and sequence between them for TGEE
                    and one guide sequence on gene strand for Cas9.
                    """
                    guide_full_seq = guide_info(
                        out_dict["guide_seq"],
                        out_dict["ensemble_gene_seq"],
                    )

                    out_dict["guide_full_seq"] = guide_full_seq

                    # Make sequence where search primers. 
                    left_flank = ensemble_gene_seq.split(out_dict["guide_full_seq"])[0][
                        -out_dict["max_dist"] :
                    ]
                    right_flank = ensemble_gene_seq.split(out_dict["guide_full_seq"])[
                        1
                    ][: out_dict["max_dist"]]

                    search_sequence = (
                        left_flank + out_dict["guide_full_seq"] + right_flank
                    )
                    out_dict["search_sequence"] = search_sequence
                    cut_size = (
                        len(search_sequence.split(out_dict["guide_full_seq"])[0])
                        + len(out_dict["guide_full_seq"]) // 2
                    )

                    # Possible coordinates for left and right primers
                    primer5_start = cut_size - out_dict["max_dist"]
                    primer5_end = cut_size - out_dict["min_dist"]
                    primer3_start = cut_size + out_dict["min_dist"]
                    primer3_end = cut_size + out_dict["max_dist"]

                except Exception:
                    text_error = "check guide sequence"
                    out_dict["gene_dict"] = (
                        "<span class='red-text'>"
                        + "Error: "
                        + str(text_error)
                        + "</span>"
                    )
                    return render_template(
                        "home.html",
                        out_dict=out_dict,
                        forms=forms,
                        links=short_links,
                        tables=[],
                        checked=False,
                    )

                try:
                    # Get ncbi gene information
                    amplicon_start, amplicon_end, gene_nt_id = ncbi_info(
                        out_dict["ncbi_id"],
                        out_dict["strand"],
                        out_dict["search_sequence"],
                    )
                    out_dict["amplicon_start"] = amplicon_start
                    out_dict["amplicon_end"] = amplicon_end
                    out_dict["gene_nt_id"] = gene_nt_id
                except Exception:
                    text_error = "check ncbi id"
                    out_dict["gene_dict"] = (
                        "<span class='red-text'>"
                        + "Error: "
                        + str(text_error)
                        + "</span>"
                    )
                    return render_template(
                        "home.html",
                        out_dict=out_dict,
                        forms=forms,
                        links=short_links,
                        tables=[],
                        checked=False,
                    )

            # Get url for Blast Primers tool with all nessesary parameters.
            blast_url = blast_primers(
                out_dict["search_sequence"],
                primer5_start,
                primer5_end,
                primer3_start,
                primer3_end,
                out_dict["min_size"],
                out_dict["max_size"],
            )
            out_dict["blast_url"] = blast_url

            short_links = {"Blast": out_dict["blast_url"]}

        if "blast_info_form_submit" in request.form:
            blast_id = blast_info_form.text_field9.data
            out_dict["blast_id"] = blast_id

            checkbox_all_value = request.form.get("checkbox_all")
            return_all = checkbox_all_value == "all_results"

            try:
                # Get all results from Blast Primers tool and compute nessesary parameters
                all_primers = blast_results(
                    out_dict["blast_id"],
                    out_dict["gene_nt_id"],
                    out_dict["amplicon_start"],
                    out_dict["amplicon_end"],
                    out_dict["search_sequence"],
                    out_dict["gene_name"],
                    out_dict["guide_name"],
                    return_all=return_all,
                )
            except Exception:
                text_error = "check blast id"
                out_dict["gene_dict"] = (
                    "<span class='red-text'>" + "Error: " + str(text_error) + "</span>"
                )
                return render_template(
                    "home.html",
                    out_dict=out_dict,
                    forms=forms,
                    links=short_links,
                    tables=[],
                    checked=False,
                )

            # Save all found primers to flask table
            tables = [all_primers.to_html(classes="data")]
            out_dict["tables"] = tables
            out_dict["tables_df"] = [all_primers]

            selected_ids = request.form.getlist("checkbox")

            # Add a column for checkboxes to choose good primers
            all_primers["ID_index"] = all_primers.index
            all_primers["ID"] = all_primers["ID_index"].apply(
                lambda x: f'<input type="checkbox" name="checkbox" value="{x}">'
            )

            # Reorder columns to put 'Select' first
            cols = all_primers.columns.tolist()
            cols = cols[-1:] + cols[:-1]
            all_primers = all_primers[cols]

            html_table = all_primers.to_html(classes="data", escape=False, index=False)

            tables = [html_table]
            out_dict["tables"] = tables

            out_dict["all_primers"] = all_primers

        if "checkbox_table_submit" in request.form:
            if len(out_dict["all_primers"]) == 0:
                text_error = "Make a primer table"
                out_dict["gene_dict"] = (
                    "<span class='red-text'>" + "Error: " + str(text_error) + "</span>"
                )
                return render_template(
                    "home.html",
                    out_dict=out_dict,
                    forms=forms,
                    links=short_links,
                    tables=[],
                    checked=False,
                )
            # Get all selected primers ID and filter table
            selected_ids = request.form.getlist("checkbox")

            all_primers_selected = out_dict["all_primers"][
                out_dict["all_primers"]["ID_index"].apply(
                    lambda p: str(p) in selected_ids
                )
            ]
            all_primers_selected.drop(columns=["ID", "ID_index"], inplace=True)

            all_primers_selected["gene_name"] = out_dict["gene_name"]
            all_primers_selected["guide_name"] = out_dict["guide_name"]

            if (out_dict["ensemble_gene_seq"] != "") & (
                out_dict["guide_seq"] != ""
            ):
                oligos_left, oligos_right = primers_coords(
                    out_dict["ensemble_gene_seq"], all_primers_selected
                )
                all_primers_selected = pd.merge(
                    all_primers_selected,
                    oligos_left,
                    on=["left_name", "Sequence (5'->3')_L"],
                )
                all_primers_selected = pd.merge(
                    all_primers_selected,
                    oligos_right,
                    on=["right_name", "Sequence (5'->3')_R"],
                )

                cut_site_coords = (
                    len(
                        out_dict["ensemble_gene_seq"].split(
                            out_dict["guide_full_seq"]
                        )[0]
                    )
                    + len(out_dict["guide_full_seq"]) // 2
                    + 1
                )
                all_primers_selected["cut_site_coords"] = cut_site_coords

            file_path = output_directory = os.path.join(
                "src",
                "static",
                "outputs",
                out_dict["gene_name"],
                out_dict["guide_name"] + "_selected_primers.csv",
            )
            all_primers_selected.to_csv(file_path, index=None)

            file_names = (
                "Primers_"
                + out_dict["gene_name"]
                + "_"
                + out_dict["guide_name"]
                + "_"
                + DATE_TODAY
            )
            # Lists for guide and primers for SnapGene file
            if (out_dict["insert_seq"] != "") & (out_dict["guide_seq"] == ""):
                elements_list, oligos = find_elements(
                    out_dict["search_sequence"],
                    out_dict["guide_full_seq"],
                    out_dict["guide_name"],
                    all_primers_selected,
                    is_guide=False,
                )
            else:
                elements_list, oligos = find_elements(
                    out_dict["search_sequence"],
                    out_dict["guide_full_seq"],
                    out_dict["guide_name"],
                    all_primers_selected,
                    is_guide=True,
                )

            # Make Gene Bank file for SnapGene tool
            gene_bank_file(
                out_dict["gene_name"],
                out_dict["search_sequence"],
                DATE_TODAY,
                elements_list,
                file_names,
                oligos=oligos,
            )

            # Make a table with amplicon sizes all possible combination between selected primers.
            primers_pivot_table(
                all_primers_selected,
                out_dict["gene_name"],
                out_dict["guide_name"],
                out_dict["guide_full_seq"],
                out_dict["min_dist"],
                out_dict["max_dist"],
                out_dict["min_size"],
                out_dict["max_size"],
                out_dict["insert_seq"],
            )

            files_path = os.path.join(
                "src", "static", "outputs", out_dict["gene_name"]
            )
            output_files = os.listdir(files_path)

            # Create a BytesIO object to store the ZIP file
            zip_buffer = BytesIO()

            # Create a ZipFile object
            with zipfile.ZipFile(
                zip_buffer, "a", zipfile.ZIP_DEFLATED, False
            ) as zip_file:
                for file_name in output_files:
                    # Add the FASTA file to the ZIP file with a custom name
                    file_path = os.path.join(
                        "src", "static", "outputs", out_dict["gene_name"], file_name
                    )
                    zip_file.write(file_path, arcname=file_name)

            # Move the buffer's position to the beginning to ensure all the data is read
            zip_buffer.seek(0)

            # Return the ZIP file as an attachment
            return send_file(
                zip_buffer,
                download_name=file_names + ".zip",
                as_attachment=True,
            )
        if "clear_forms_submit" in request.form:
            # Clear all forms
            for key in out_dict.keys():
                if key in out_dict_start.keys():
                    out_dict[key] = out_dict_start[key]
                else:
                    out_dict[key] = ""

        if "file_upload_submit" in request.form:
            # upload fasta file with found primers for few guides
            if "file" not in request.files:
                return "No file part"

            file = request.files["file"]
            filename = file.filename

            save_name, file_type = (
                ".".join(filename.split(".")[:-1]),
                filename.split(".")[-1],
            )
            if file_type == "csv":
                primers_table = pd.read_csv(file)
            elif file_type == "xlsx":
                primers_table = pd.read_excel(file)

            """
            Make a table with amplicon sizes and distances all possible combination 
            between selected primers for few guides. 
            """

            primers_pivot_table_few_guides(
                primers_table,
                out_dict["min_dist"],
                out_dict["max_dist"],
                out_dict["min_size"],
                out_dict["max_size"],
                save_name,
            )

            file_names = save_name + "_distances.csv" + "_" + DATE_TODAY

            # Create a BytesIO object to store the ZIP file
            zip_buffer = BytesIO()

            # Create a ZipFile object
            with zipfile.ZipFile(
                zip_buffer, "a", zipfile.ZIP_DEFLATED, False
            ) as zip_file:
                file_path = os.path.join(
                    "src", "static", "outputs", save_name + "_distances.csv"
                )
                zip_file.write(file_path, arcname=save_name + "_distances.csv")

            # Move the buffer's position to the beginning to ensure all the data is read
            zip_buffer.seek(0)

            # Return the ZIP file as an attachment
            return send_file(
                zip_buffer,
                download_name=file_names + ".zip",
                as_attachment=True,
            )

        # Save selected parameters to input windows
        gene_info_form.text_field.default = out_dict["gene_name"]
        gene_info_form.text_field2.default = out_dict["ncbi_id"]
        gene_info_form.text_field3.default = out_dict["guide_seq"]
        gene_info_form.text_field3_2.default = out_dict["guide_name"]
        gene_info_form.text_field4.default = out_dict["min_dist"]
        gene_info_form.text_field5.default = out_dict["max_dist"]
        gene_info_form.text_field6.default = out_dict["min_size"]
        gene_info_form.text_field7.default = out_dict["max_size"]
        gene_info_form.text_field8.default = out_dict["insert_seq"]

        gene_info_form.process()

        blast_info_form.text_field9.default = out_dict["blast_id"]

        blast_info_form.process()

        forms = {"gene_info_form": gene_info_form, "blast_info_form": blast_info_form}

        return render_template(
            "home.html",
            out_dict=out_dict,
            forms=forms,
            links=short_links,
            tables=out_dict["tables"],
            checked=checkbox_all_cond
        )

    # Save selected parameters to input windows
    gene_info_form.text_field.default = out_dict["gene_name"]
    gene_info_form.text_field2.default = out_dict["ncbi_id"]
    gene_info_form.text_field3.default = out_dict["guide_seq"]
    gene_info_form.text_field3_2.default = out_dict["guide_name"]
    gene_info_form.text_field4.default = out_dict["min_dist"]
    gene_info_form.text_field5.default = out_dict["max_dist"]
    gene_info_form.text_field6.default = out_dict["min_size"]
    gene_info_form.text_field7.default = out_dict["max_size"]
    gene_info_form.text_field8.default = out_dict["insert_seq"]

    gene_info_form.process()

    blast_info_form.text_field9.default = out_dict["blast_id"]

    blast_info_form.process()

    forms = {"gene_info_form": gene_info_form, "blast_info_form": blast_info_form}

    return render_template(
        "home.html",
        out_dict=out_dict,
        forms=forms,
        links=short_links,
        tables=out_dict["tables"],
        checked=checkbox_all_cond
    )


@app.route("/", methods=["GET", "POST"])
def root():
    """
    Load main page of server
    """
    return index(out_dict_outer)


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000)
