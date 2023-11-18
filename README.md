# HighGene

HighGene is a Python package enabling the visualisation of normalised gene counts. Instead of displaying absolute gene counts,  will generate a figure with normalised counts based on the total amount of genes per species. Subsequently the Z-Score is calculated for each class. HighGene contains analyses built specifically for [JGI's Mycocosm](https://mycocosm.jgi.doe.gov/mycocosm/home) formatted Gene Class Classification results (eg., KOG/KEGG Classes), but allows for the visualisation of custom files.

## Installation

HighGene can be installed by first cloning this GitHub repository. Subsequently, go into the location where the repository was saved and use pip to install the package.
```bash
# Clone repository
git clone https://github.com/WesterdijkInstitute/HighGene.git

# Go into repository location
cd .../HighGene/

# Install package with pip
pip install .
```

## Usage

HighGene requires 3 inputs, the input folder, output folder and file name and the analysis mode.

```bash
usage: HighGene [-h] -i INPUT -o OUTPUT -m {0,1,2,3,4,5}

HighGene - Highlighting gene content varation, For more information see: https://github.com/WesterdijkInstitute/HighGene

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input folder
  -o OUTPUT, --output OUTPUT
                        Output folder + File name
  -m {0,1,2,3,4,5}, --mode {0,1,2,3,4,5}
                        Analysis modes:
                            0: Custom File
                            1: Transcription Factors
                            2: Peptidase Clans
                            3: CAZymes
                            4: KEGG
                            5: KOG
```

In addition to a protein fasta file of the proteome of each species, the different analysis modes all require different files to be present in the input folder. They are the following:

- Custom file:
    - A .tsv file with two columns and headers for each species.
      - 1: Protein ID
      - 2: Class prediction

- Transcription Factors:
    - An InterPro classification file (JGI) for each species. 'IPR' must be included in the file name
    - A .tsv file named 'TF_PFAM.tsv' containing the Pfam Accessions for the transcription factors of interest.

- Peptidase Clans:
    - An InterPro classification file (JGI) for each species. 'IPR' must be included in the file name
    - A .tsv file named 'Pep_PFAM.tsv' containing the Pfam Accessions for the peptidase clans of interest.

- CAZymes:
    - A CAZyme prediction results from [dbCAN3](https://bcb.unl.edu/dbCAN2/blast.php) for each species. 'CAZy' must be included in the file name

- KEGG:
    - A KEGG classification file (JGI) for each species. 'KEGG' must be included in the file name

- KOG:
    - A KOG classification file (JGI) for each species. 'KOG' must be included in the file name

## Important formatting information
- Headers of classification files should always start with a #

- It is important that all files start with the same abbreviated species name that is also used in the protein fasta file followed by an underscore (e.g., Species1_Genes.fasta). Additionally, the abbreviated species names may ***NOT*** contain an underscore (e.g., Species1_A_Genes.fasta). If this is the case in the protein fasta files, change the names accordingly. HighGene uses this to determine which fasta file belongs to which gene classification files.

- When downloading protein fasta files from JGI's Mycocosm, the fasta headers will be formatted as: 'jgi|Species1|0000001|...'. HighGene requires this fasta header format. If a different header is present, change it to something that looks like JGI's format (e.g., abc|Species1).

- Examples of all file format types can be found under the 'Examples' folder.

## Author
Tim Verschuren <br/>
[GitHub](https://github.com/TimVerschuren) [LinkedIn](https://www.linkedin.com/in/tim-verschuren-27082919b/)
