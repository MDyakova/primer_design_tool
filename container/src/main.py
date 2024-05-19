"""
Code launch FLASK server to find primers for gene editing tasks
"""
import os
from io import BytesIO
import zipfile
import pandas as pd
import numpy as np
import shutil
from datetime import date

from flask import Flask, render_template, request, send_file
from flask_wtf import FlaskForm
from wtforms import StringField, TextAreaField

from utilities import ensemble_info, ncbi_info, guide_info, blast_primers, blast_results

app = Flask(__name__)
app.config["SECRET_KEY"] = "mysecretkey"  # fake key to work with flask server

date_today = str(date.today())


# dictionary with all useful variables necessary to save throw whole pipeline
out_dict = {
    "gene_name": "",
    "ncbi_id": "",
    "guide_seq": "",
    "guide_name": "",
    "min_dist": 50,
    "max_dist": 150,
    "min_size": 150,
    "max_size": 300,
    "gene_dict":"",
    "ensemble_gene_seq":'',
    "strand":"",
    "guide_full_seq":"",
    "search_sequence":"",
    "blast_url":"https://www.ncbi.nlm.nih.gov/tools/primer-blast/",
    "blast_id":"",
    "gene_nt_id":"",
    "tables":[],
    "tables_df":[]


}

out_dict_start = out_dict.copy()

out_keys_clean = []

compl_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', '-':'-'}


class GeneInfo(FlaskForm):
    """
    Information about ensemble and ncbi gene names or id
    """
    text_field = StringField("Ensemble gene name", default="")
    text_field2 = StringField("  NCBI gene id", default="")
    text_field3 = StringField(
        "Guide sequence", default="", render_kw={"style": "width: 400px;"}
    )
    text_field3_2 = StringField("Guide name", default="", 
                                render_kw={"style": "width: 100px;"})
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

class BlastInfo(FlaskForm):
    """
    Blast results id
    """
    text_field8 = StringField(
        "Blast results id", default="", render_kw={"style": "width: 400px;"}
    )





def index(out_dict):
    """
    Launch main code to prepare input data for golden gate primers
    """

    gene_info_form = GeneInfo()
    blast_info_form = BlastInfo()

    forms = {
        "gene_info_form":gene_info_form,
        "blast_info_form":blast_info_form
    }

    short_links = {'Blast': out_dict['blast_url']}

    out_dict["error_message"] = ''

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
            min_size= gene_info_form.text_field6.data
            max_size = gene_info_form.text_field7.data


            if (gene_name == '') | (ncbi_id == '') | (guide_seq == ''):
                text_error = 'enter all data'
                out_dict["gene_dict"] = ("<span class='red-text'>" 
                                         + 'Error: ' + str(text_error)
                                         + "</span>")
                return render_template("home.html", out_dict=out_dict, forms=forms, links=short_links)
            
            out_dict["gene_name"] = gene_name
            out_dict["ncbi_id"] = ncbi_id
            out_dict["guide_seq"] = guide_seq
            out_dict["guide_name"] = guide_name

            # Create output directory
            output_directory = "src/static/outputs/" + gene_name
            os.makedirs(output_directory, exist_ok=True)
            
            try:
                out_dict["min_dist"] = int(min_dist)
                out_dict["max_dist"] = int(max_dist)
                out_dict["min_size"] = int(min_size)
                out_dict["max_size"] = int(max_size)
            except Exception as e:
                text_error = 'primer position parameters must be > 0'
                out_dict["gene_dict"] = ("<span class='red-text'>" 
                                         + 'Error: ' + str(text_error)
                                         + "</span>")
                return render_template("home.html", out_dict=out_dict, forms=forms, links=short_links)
            
            try:
                ensemble_gene_seq, gene_dict, strand = ensemble_info(gene_name)
                out_dict["ensemble_gene_seq"] = ensemble_gene_seq
                out_dict["strand"] = strand
                out_dict["gene_dict"] = 'Gene info: ' + str(gene_dict)
            except Exception as e:
                text_error = 'check gene name'
                out_dict["gene_dict"] = ("<span class='red-text'>" 
                                         + 'Error: ' + str(text_error)
                                         + "</span>")
                return render_template("home.html", out_dict=out_dict, forms=forms, links=short_links)
            
            try:
                guide_full_seq = guide_info(
                                            out_dict["guide_seq"], out_dict["strand"], 
                                            out_dict["ensemble_gene_seq"]
                )

                out_dict["guide_full_seq"] = guide_full_seq

                left_flank = ensemble_gene_seq.split(out_dict["guide_full_seq"])[0][-out_dict["max_dist"]:]
                right_flank = ensemble_gene_seq.split(out_dict["guide_full_seq"])[1][:out_dict["max_dist"]]

                search_sequence = left_flank + out_dict["guide_full_seq"] + right_flank
                out_dict['search_sequence'] = search_sequence
                cut_size = len(search_sequence)//2

                primer5_start = 1
                primer5_end = cut_size - out_dict["min_dist"]
                primer3_start = cut_size + out_dict["min_dist"]
                primer3_end = len(search_sequence)

            except Exception as e:
                text_error = 'check guide sequence'
                out_dict["gene_dict"] = ("<span class='red-text'>" 
                                         + 'Error: ' + str(text_error)
                                         + "</span>")
                return render_template("home.html", out_dict=out_dict, forms=forms, links=short_links)
            
            try:
                amplicon_start, amplicon_end, gene_nt_id = ncbi_info(out_dict["ncbi_id"], 
                                                            out_dict["strand"], 
                                                            out_dict['search_sequence'])
                out_dict['amplicon_start'] = amplicon_start
                out_dict['amplicon_end'] = amplicon_end
                out_dict['gene_nt_id'] = gene_nt_id
            except:
                text_error = 'check ncbi id'
                out_dict["gene_dict"] = ("<span class='red-text'>" 
                                         + 'Error: ' + str(text_error)
                                         + "</span>")
                return render_template("home.html", out_dict=out_dict, forms=forms, links=short_links)
            
            blast_url = blast_primers(out_dict['search_sequence'] , 
                                        primer5_start, 
                                        primer5_end, 
                                        primer3_start, 
                                        primer3_end, 
                                        out_dict["min_size"], 
                                        out_dict["max_size"])
            out_dict['blast_url'] = blast_url

            short_links = {'Blast': out_dict['blast_url']}

        if "blast_info_form_submit" in request.form:
            blast_id = blast_info_form.text_field8.data
            out_dict['blast_id'] = blast_id

            checkbox_all_value = request.form.get("checkbox_all")

            if checkbox_all_value == 'all_results':
                return_all = True
            else:
                return_all = False

            all_primers = blast_results(out_dict['blast_id'], 
                                        out_dict['gene_nt_id'], 
                                        out_dict['amplicon_start'], 
                                        out_dict['amplicon_end'], 
                                        out_dict['search_sequence'], 
                                        out_dict["gene_name"], 
                                        out_dict["guide_name"],
                                        return_all=return_all)

            tables=[all_primers.to_html(classes='data')]
            out_dict['tables'] = tables
            out_dict['tables_df'] = [all_primers]
           
            selected_ids = request.form.getlist('checkbox')

            # # Example DataFrame
            # all_primers = pd.DataFrame({
            #     'ID': [1, 2, 3],
            #     'Name': ['Primer1', 'Primer2', 'Primer3'],
            #     'Sequence': ['ATCG', 'CGTA', 'GTAC'],
            #     # Add other columns as necessary
            # })

            # Add a column for checkboxes
            all_primers['ID'] = all_primers.index
            all_primers['ID'] = all_primers['ID'].apply(lambda x: f'<input type="checkbox" name="checkbox" value="{x}">')

            # Reorder columns to put 'Select' first
            cols = all_primers.columns.tolist()
            cols = cols[-1:] + cols[:-1]
            all_primers = all_primers[cols]

            html_table = all_primers.to_html(classes='data', escape=False, index=False)

            tables = [html_table]
            out_dict['tables'] = tables

        if "checkbox_table_submit" in request.form:
            selected_ids = request.form.getlist('checkbox')
            # out_dict["gene_dict"] = selected_ids

            all_primers_selected = all_primers[all_primers['ID'].apply(lambda p: p in selected_ids)]
            all_primers_selected.to_csv('src/static/outputs/' + out_dict["gene_name"]  + '/' + out_dict["guide_name"] 
                               + ' ' + guide_name +'__selected_primers.csv', index=None)
            
            output_files = os.listdir('src/static/outputs/' + gene_name + '/')

            file_names = ('Primers_' 
                            + out_dict["gene_name"] 
                            + '_' 
                            + out_dict["guide_name"]
                            + '_'
                            + date_today)
            
            # Create a BytesIO object to store the ZIP file
            zip_buffer = BytesIO()

            # Create a ZipFile object
            with zipfile.ZipFile(
                zip_buffer, "a", zipfile.ZIP_DEFLATED, False
            ) as zip_file:
                
                for file_name in output_files:
                    # Add the FASTA file to the ZIP file with a custom name
                    zip_file.write('src/static/outputs/' + gene_name + '/' + file_name, 
                                   arcname=file_name)

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
        gene_info_form.text_field4.default = out_dict["min_dist"]
        gene_info_form.text_field5.default = out_dict["max_dist"]
        gene_info_form.text_field6.default = out_dict["min_size"]
        gene_info_form.text_field7.default = out_dict["max_size"]

        gene_info_form.process()

        blast_info_form.text_field8.default = out_dict["blast_id"]

        blast_info_form.process()

        forms = {
            "gene_info_form": gene_info_form,
            "blast_info_form":blast_info_form
        }
        

        return render_template("home.html", out_dict=out_dict, forms=forms, links=short_links,
                               tables=out_dict['tables']
                            #    data=data
                               )

    # Save selected parameters to input windows
    gene_info_form.text_field.default = out_dict["gene_name"]
    gene_info_form.text_field2.default = out_dict["ncbi_id"]
    gene_info_form.text_field3.default = out_dict["guide_seq"]
    gene_info_form.text_field4.default = out_dict["min_dist"]
    gene_info_form.text_field5.default = out_dict["max_dist"]
    gene_info_form.text_field6.default = out_dict["min_size"]
    gene_info_form.text_field7.default = out_dict["max_size"]

    gene_info_form.process()

    blast_info_form.text_field8.default = out_dict["blast_id"]

    blast_info_form.process()

    forms = {
        "gene_info_form": gene_info_form,
        "blast_info_form":blast_info_form
    }

    return render_template("home.html", out_dict=out_dict, forms=forms,links=short_links,
                           tables=out_dict['tables']
                        #    data=data
                           )


@app.route("/", methods=["GET", "POST"])
def root():
    """
    Load main page of server
    """
    return index(out_dict)


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000)
