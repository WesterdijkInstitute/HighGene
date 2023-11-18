#!/usr/bin/env python

"""
HighGene is a Python tool that allows for the visualisation 
of normalised gene counts. It was created for gene class predictions in the 
style of JGI, but also allows for custom files to be visualised.

Example:

    $ highgene [-h] -i Input_folder -o Output_folder -m {0,1,2,3,4,5}

"""

"""Import Statements"""
import os
import pandas as pd
from Bio import SeqIO
from scipy import stats
import plotly.express as px
import plotly.figure_factory as ff

"""Authorship Information"""
__author__ = "Tim Verschuren"
__credits__ = ["Tim Verschuren", "Jérôme Collemare"]

__licence__ = "MIT"
__date__ = "06-09-2023"
__version__ = "0.1.0"
__maintainer__ = "Tim Verschuren"
__email__ = "t.verschuren@wi.knaw.nl"
__status__ = "Development"


class FastaReader:
    """
    This class reads the fasta files in the input folder and calculates the 
    number of genes per species

    Attributes:
    input (str): Name of input folder containing fasta files

    Returns:
    gene_count (dict): Dictionary containing the number of genes per species.
    """
    def __init__(self, input_folder: str):
        self.input = input_folder
    
    # Read content of fasta files and calculate number of protein coding genes
    def read_fasta(self) -> dict:
        gene_count = {}
        for file in os.listdir(self.input):
            if file.endswith("fasta"):
                kegg_dict = {}
                for record in SeqIO.parse(f"{self.input}/{file}", "fasta"):
                    species_name = (record.id).split("|")[1]
                    kegg_dict[record.id] = [record.seq]
                gene_count[species_name] = len(kegg_dict)
        return gene_count


class TfData:
    """
    Class for the visualisation of TF domains from IPR classification files. 

    Attributes:
    input (str): Path to folder containing input files
    """  
    # Load in input folder and output file name
    def __init__(self, input_folder: str):
        self.input = input_folder
        self.tf_data = {"species": [], "factors": [], "counts": []}
    
    def tf_reader(self) -> pd.DataFrame:
        """
        Takes the IPR input files in the input folder and 
        calculates the total number of genes per TF domain before
        dividing the counts by the total number of genes per species. 
        """
        tf_pfam_list = []
        tf_name = {}
        gene_count = FastaReader(self.input).read_fasta()

        # Read content of transcription factor Pfam file
        with open(f"{self.input}/TF_PFAM.tsv", "r") as pfam:
            for line in pfam.readlines():
                pfam_split = line.split("\t")
                tf_pfam_list.append(pfam_split[0])
                tf_name[pfam_split[0]] = pfam_split[1]

        # Count number of transcription factors
        for file in os.listdir(self.input):
            if file.__contains__("IPR"):
                with open(f"{self.input}/{file}", "r") as ipr_file:
                    tf_pfam_counts = {}
                    for line in ipr_file.readlines():
                        if line.startswith("#"):
                            pass
                        else:
                            split_line = (line.split("\t"))
                            if split_line[4] in tf_pfam_list:
                                if split_line[4].__contains__("PF"):
                                    if split_line[4] in tf_pfam_counts:
                                        tf_pfam_counts[split_line[4]] += 1
                                    else:
                                        tf_pfam_counts[split_line[4]] = 1

                # Write extracted data to a dictionary
                for key, value in tf_pfam_counts.items():
                    self.tf_data["species"].append(str(file).split("_")[0])
                    self.tf_data["factors"].append(key)
                    self.tf_data["counts"].append(
                        value/gene_count[str(file).split("_")[0]])

        # Convert dictionary to dataframe
        df = pd.DataFrame(data = self.tf_data)
        df = df.pivot(index="species", columns="factors", values="counts")
        return df.rename(columns=tf_name), "Transcription Factor Domains"

    
class CazyData:
    """
    Class for the visualisation of CAZyme classification files. 

    Attributes:
    input (str): Path to folder containing input files
    """  
    def __init__(self, input_folder: str):
        self.input = input_folder
        self.cazy_data = {"species": [], "cazymes": [], "counts": []}
        self.cazy_names = {"AA": "Auxiliary Activities Family ", 
                           "CBM": "Carbohydrate-Binding Module Family ", 
                           "CE": "Carbohydrate Esterase Family ", 
                           "GH": "Glycoside Hydrolase Family ",
                           "GT": "Glycosyl Transferase Family ",
                           "PL": "Polysaccharide Lyase Family "}

    def cazy_reader(self) -> pd.DataFrame:
        """
        Takes the CAZy input files in the input folder and 
        calculates the total number of genes per family before
        dividing the counts by the total number of genes per species. 
        """
        gene_count = FastaReader(self.input).read_fasta()
        cazy_dict = {}
        # Read content of CAZy files
        for file in os.listdir(self.input):
            if file.__contains__("CAZy"):
                cazy_families = {"AA": 0, "CBM": 0, "CE": 0, 
                                 "GH": 0, "GT": 0, "PL": 0}
                with open(f"{self.input}/{file}", "r") as cazy_file:
                    for line in cazy_file.readlines():
                        if line.startswith("#"):
                            pass
                        else:
                            # Determine CAZyme family
                            split_line = line.split("\t")
                            for item in split_line[2].split("+"):
                                if item.__contains__("AA"):
                                    cazy_families["AA"] += 1
                                if item.__contains__("CBM"):
                                    cazy_families["CBM"] += 1
                                if item.__contains__("CE"):
                                    cazy_families["CE"] += 1
                                if item.__contains__("GH"):
                                    cazy_families["GH"] += 1
                                if item.__contains__("GT"):
                                    cazy_families["GT"] += 1
                                if item.__contains__("PL"):
                                    cazy_families["PL"] += 1

                    # Check if CAZyme data is devided over multiple files
                    if str(file).split("_")[0] in cazy_dict:
                        for key, value in cazy_families.items():
                            cazy_dict[str(file).split("_")[0]][key] += value
                    else:
                        cazy_dict[str(file).split("_")[0]] = cazy_families
        
        # Write extracted data to a dictionary
        for key, value in cazy_dict.items():
            species_name = key
            for key, value in cazy_dict[species_name].items():
                self.cazy_data["species"].append(species_name)
                self.cazy_data["cazymes"].append(key)
                self.cazy_data["counts"].append(value/gene_count[species_name])
        
        # Convert dictionary to dataframe
        df = pd.DataFrame(data = self.cazy_data)
        df = df.pivot(index="species", columns="cazymes", values="counts")
        return df.rename(columns=self.cazy_names), "CAZyme Families"


class KeggData:
    """
    Class for the visualisation of KEGG class files. 

    Attributes:
    input (str): Path to folder containing input files
    """  
    def __init__(self, input_folder: str):
        self.input = input_folder
        self.kegg_data = {"species": [], "KEGG": [], "counts": []}

    def kegg_reader(self) -> pd.DataFrame:
        """
        Takes the KEGG input files in the input folder and 
        calculates the total number of genes per class before
        dividing the counts by the total number of genes per species. 
        """
        gene_count = FastaReader(self.input).read_fasta()
        kegg_out = {}
        # Read content of KEGG files
        for file in os.listdir(self.input):
            if file.__contains__("KEGG"):
                with open(f"{self.input}/{file}", "r") as kegg_file:
                    kegg_dict = {}
                    kegg_counts = {}
                    for line in kegg_file.readlines():
                        if line.startswith("#"):
                            pass
                        else:
                            split_line = line.split("\t")
                            if split_line[0] not in kegg_dict:
                                kegg_dict[split_line[0]] = [split_line[7]]
                            else:
                                kegg_dict[split_line[0]].append(split_line[7])
                    # Determine if a gene has multiple different KEGG_classes
                    for key, value in kegg_dict.items():
                        if len(value) > 1:
                            if value.count(value[0]) == len(value):
                                kegg_dict[key] = value[0]
                            else:
                                kegg_dict[key] = "Undetermined"
                        else:
                            kegg_dict[key] = value[0]
                # Count occurance of each KEGG_class
                for key, value in kegg_dict.items():
                    if value in kegg_counts:
                        kegg_counts[value] += 1
                    else:
                        kegg_counts[value] = 1
                kegg_out[str(file).split("_")[0]] = kegg_counts
        
        # Write extracted data to a dictionary
        for key, value in kegg_out.items():
            species_name = key
            for key, value in kegg_out[species_name].items():
                if not key == "\\N":
                    if key == "SORTING AND DEGRADATION":
                        key = "Sorting and Degradation"
                    self.kegg_data["species"].append(species_name)
                    self.kegg_data["KEGG"].append(key)
                    self.kegg_data["counts"].append(
                        value/gene_count[species_name])

        # Convert dictionary to dataframe
        df = pd.DataFrame(data = self.kegg_data)
        return df.pivot(
            index="species", columns="KEGG", values="counts"), "KEGG Classes"


class KogData:
    """
    Class for the visualisation of KOG class files. 

    Attributes:
    input (str): Path to folder containing input files
    """
    def __init__(self, input_folder: str):
        self.input = input_folder
        self.kog_data = {"species": [], "KOG": [], "counts": []}
        self.column_order = [
            "Cell wall/membrane/envelope biogenesis ",
            "Cell motility ",
            "Posttranslational modification, protein turnover, chaperones ",
            "Signal transduction mechanisms ",
            "Intracellular trafficking, secretion, and vesicular transport ",
            "Defense mechanisms ",
            "Extracellular structures ",
            "Nuclear structure ",
            "Cytoskeleton ",
            "RNA processing and modification ",
            "Chromatin structure and dynamics ",
            "Translation, ribosomal structure and biogenesis ",
            "Transcription ",
            "Replication, recombination and repair ",
            "Energy production and conversion ",
            "Cell cycle control, cell division, chromosome partitioning ",
            "Amino acid transport and metabolism ",
            "Nucleotide transport and metabolism ",
            "Carbohydrate transport and metabolism ",
            "Coenzyme transport and metabolism ",
            "Lipid transport and metabolism ",
            "Inorganic ion transport and metabolism ",
            "Secondary metabolites biosynthesis, transport and catabolism ",
            "General function prediction only ",
            "Function unknown ",
            "Undetermined "]
    
    def kog_reader(self) -> pd.DataFrame:
        """
        Takes the KOG input files in the input folder and 
        calculates the total number of genes per class before
        dividing the counts by the total number of genes per species. 
        """
        gene_count = FastaReader(self.input).read_fasta()
        kog_out = {}
        # Read content of KOG files
        for file in os.listdir(self.input):
            if file.__contains__("KOG"):
                with open(f"{self.input}/{file}", "r") as kog_file:
                    kog_counts = {}
                    kog_dict = {}
                    for line in kog_file.readlines():
                        if line.startswith("#"):
                            pass
                        else:
                            split_line = line.split("\t")
                            if split_line[1] not in kog_dict:
                                kog_dict[split_line[1]] = [split_line[4]]
                            else:
                                kog_dict[split_line[1]].append(split_line[4])
                    # Determine if a gene has multiple different KOG_classes
                    for key, value in kog_dict.items():
                        if len(value) > 1:
                            if value.count(value[0]) == len(value):
                                kog_dict[key] = value[0]
                            else:
                                kog_dict[key] = "Undetermined "
                        else:
                            kog_dict[key] = value[0]
                # Count occurance of each KOG_class
                for key, value in kog_dict.items():
                    if value in kog_counts:
                        kog_counts[value] += 1
                    else:
                        kog_counts[value] = 1
                kog_out[str(file).split("_")[0]] = kog_counts

        # Write extracted data to a dictionary
        for key, value in kog_out.items():
            species_name = key
            for key, value in kog_out[species_name].items():
                self.kog_data["species"].append(species_name)
                self.kog_data["KOG"].append(key)
                self.kog_data["counts"].append(value/gene_count[species_name])

        # Convert dictionary to dataframe
        df = pd.DataFrame(data = self.kog_data)
        df = df.pivot(index="species", columns="KOG", values="counts")
        df = df.rename(columns={"RNA processing and modification  ": 
                                "RNA processing and modification "})
        return df.reindex(columns=self.column_order), "KOG Classes"


class PepData:
    """
    Class for the visualisation of peptidase clans from IPR classification 
    files. 

    Attributes:
    input (str): Path to folder containing input files
    """
    def __init__(self, input_folder: str, mode: int):
        self.input = input_folder
        self.pep_data = {"species": [], "peptidase": [], "counts": []}
        self.mode = int(mode)

    def pep_reader(self) -> pd.DataFrame:
        """
        Takes the IPR input files in the input folder and 
        calculates the total number of genes per peptidase clan before
        dividing the counts by the total number of genes per species. 
        """
        pep_pfam_list = []
        pep_name = {}
        gene_count = FastaReader(self.input).read_fasta()

        # Read content of transcription factor Pfam file
        with open(f"{self.input}/Pep_PFAM.tsv", "r") as pfam:
            for line in pfam.readlines():
                pfam_split = line.split("\t")
                pep_pfam_list.append(pfam_split[2].replace("\n", ""))
                pep_name[pfam_split[2].replace("\n", "")] = \
                pfam_split[self.mode].replace("\n", "")
        # Count number of transcription factors
        for file in os.listdir(self.input):
            if file.__contains__("IPR"):
                with open(f"{self.input}/{file}", "r") as ipr_file:
                    pep_pfam_counts = {}
                    for line in ipr_file.readlines():
                        if line.startswith("#"):
                            pass
                        else:
                            split_line = (line.split("\t"))
                            if split_line[4] in pep_pfam_list:
                                if split_line[4].__contains__("PF"):
                                    if split_line[4] in pep_pfam_counts:
                                        pep_pfam_counts[split_line[4]] += 1
                                    else:
                                        pep_pfam_counts[split_line[4]] = 1       
                                 
                # Write extracted data to a dictionary
                pep_clan_counts = {}
                for key, value in pep_pfam_counts.items():
                    clan = pep_name[key]
                    if clan in pep_clan_counts:
                        pep_clan_counts[clan] += value
                    else:
                        pep_clan_counts[clan] = value
                if "N/A" in pep_clan_counts:
                    pep_clan_counts.pop("N/A")
                for key, value in pep_clan_counts.items():
                    self.pep_data["species"].append(str(file).split("_")[0])
                    self.pep_data["peptidase"].append(key)
                    self.pep_data["counts"].append(
                        value/gene_count[str(file).split("_")[0]])

        # Convert dictionary to dataframe
        df = pd.DataFrame(data = self.pep_data)
        return df.pivot(index="species", 
                        columns="peptidase", 
                        values="counts"), "Peptidase Clan"


class CustomData:
    """
    Class for the visualisation of custom gene class files. 

    Attributes:
    input (str): Path to folder containing input files
    mode (str): (y/n) value based on desired analysis type
    """
    def __init__(self, input_folder: str, mode: str):
        self.input = input_folder
        self.custom_data = {"species": [], "var_fac": [], "counts": []}
        self.mode = mode

    def custom_reader(self) -> pd.DataFrame:
        """
        Takes the custom input files in the input folder and 
        calculates the total number of genes per class before
        dividing the counts by the total number of genes per species. 
        """
        gene_count = FastaReader(self.input).read_fasta()
        custom_out = {}
        for file in os.listdir(self.input):
            if not file.__contains__("fasta"):
                # Read content of custom file
                with open(f"{self.input}/{file}", "r") as custom_file:
                    custom_dict = {}
                    custom_counts = {}
                    for line in custom_file.readlines():
                        # Extract axis titles
                        if line.startswith("#"):
                            split_line = line.split("\t")
                            header = split_line[1]
                        else:
                            split_line = line.split("\t")
                            if split_line[0] not in custom_dict:
                                custom_dict[split_line[0]] = [split_line[1]]
                            else:
                                custom_dict[split_line[0]].append(
                                    split_line[1])
                # Determine counts if multiple factors per gene are allowed
                if self.mode == "y":
                    for key, value in custom_dict.items():
                        for item in value:
                            if item in custom_counts:
                                custom_counts[item.replace("\n", "")] += 1
                            else:
                                custom_counts[item.replace("\n", "")] = 1
                # Determine counts if multiple factors per gene are not allowed
                elif self.mode == "n":
                    for key, value in custom_dict.items():
                        if len(value) > 1:
                            if value.count(value[0]) == len(value):
                                custom_dict[key] = value[0].replace("\n", "")
                            else:
                                custom_dict[key] = "Undetermined"
                        else:
                            custom_dict[key] = value[0].replace("\n", "")
                    for key, value in custom_dict.items():
                        if value in custom_counts:
                            custom_counts[value] += 1
                        else:
                            custom_counts[value] = 1
                custom_out[str(file).split("_")[0]] = custom_counts

        # Write extracted data to a dictionary
        for key, value in custom_out.items():
            species_name = key
            for key, value in custom_out[species_name].items():
                self.custom_data["species"].append(species_name)
                self.custom_data["var_fac"].append(key)
                self.custom_data["counts"].append(
                    value/gene_count[species_name])

        # Convert dictionary to dataframe
        df = pd.DataFrame(data = self.custom_data)
        return df.pivot(index="species", 
                        columns="var_fac", 
                        values="counts"), header


class HeatGenerator:
    """
    This class contains functions which allow the user to generate heatmaps 
    from gene count data.
    The data, in the form of a pandas dataframe, will first be normalised by 
    calculating the Z-score.

    Attributes:
    inputs (list): A list containing a pandas dataframe and the title of the 
    y-axis [dataframe, title].
    output (str): The name of the folder where the figures should be saved.
    """
    def __init__(self, inputs: str, output: str):
        self.dataframe = inputs[0].fillna(0).select_dtypes(
            include='number').apply(stats.zscore).transpose()
        self.output = output
        self.header = inputs[1]

    def generate(self) -> px:
        """
        Takes the normalised dataframe and performs heirarchical clustering 
        across the y-axis. Subsequently, a heatmap is generated for this data 
        with the species along the x-axis and the classes along the y-axis. 
        The heatmap will be saved as an HTML file.
        """
        dendro_side = ff.create_dendrogram(self.dataframe, orientation='right')
        dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
        dendro_leaves = list(map(int, dendro_leaves))

        factor_order = self.dataframe.loc\
            [list(self.dataframe.index[dendro_leaves])]
        fig = px.imshow(
                    factor_order, 
                    labels = dict(
                        x="Species", 
                        y=self.header, 
                        color="Z-score"
                        ), 
                    color_continuous_scale=px.colors.diverging.RdYlBu_r
                    )
        
        fig.update_xaxes(side="top", tickangle=-45)
        fig.update_layout(xaxis_title=None)
        fig.write_html(f"{self.output}.html", auto_open = True)
        