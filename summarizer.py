from Bio import SeqIO
import xlsxwriter
import os
import sys

fasta_dict = {"Protein_ID": ["Length", "SignalP", "TMHMM", "KEGG_Class", "KOG_Class", "CAZy", "Peptidases", "Transcription Factors"]}
workbook = xlsxwriter.Workbook(f"{str(sys.argv[2])}/{str(sys.argv[1])}.xlsx")
worksheet = workbook.add_worksheet()
folder_content = os.listdir(str(sys.argv[2]))


def SP_reader(input_file):
    with open(input_file, "r") as file:
        for line in file.readlines():
            if line.startswith("#"):
                pass
            else:
                split_line = line.split("\t")
                fasta_dict[int(((split_line[0]).split("|"))[2])].append(split_line[1])


def TMHMM_reader(input_files):
    for i, x in enumerate(input_files):
        with open(x, "r") as file:
            for line in file.readlines():
                if line.startswith("#"):
                    pass
                else:
                    split_line = line.split("\t")
                    fasta_dict[int(((split_line[0]).split("|"))[2])].append(int(split_line[4].replace("PredHel=", "")))
    for key, value in fasta_dict.items():
        if not len(value) == 3:
            fasta_dict[key].append("N/A")


def KEGG_reader(input_file):
    with open(input_file, "r") as file:
        for line in file.readlines():
            if line.startswith("#"):
                pass
            else:
                split_line = line.split("\t")
                if not type(fasta_dict[int(split_line[0])][-1]) == list: 
                    fasta_dict[int(split_line[0])].append([split_line[7]])
                else:
                    fasta_dict[int(split_line[0])][-1].append(split_line[7])
    for key, value in fasta_dict.items():
        if not len(value) == 4:
            fasta_dict[key].append(["\\N"])
        if len(value[-1]) > 1:
            if ([value[-1][0]]*len(value[-1]) == value[-1]):
                fasta_dict[key] = value[:-1]
                fasta_dict[key].append(value[-1][0])
            else:
                fasta_dict[key] = value[:-1]
                fasta_dict[key].append("Undetermined")
        else:
            fasta_dict[key] = value[:-1]
            fasta_dict[key].append(value[-1][0])
        if value[-1] == ["\\N"]:
            fasta_dict[key] = value[:-1]
            fasta_dict[key].append("N/A")


def KOG_reader(input_file):
    with open(input_file, "r") as file:
        for line in file.readlines():
            if line.startswith("#"):
                pass
            else:
                split_line = line.split("\t")
                if not type(fasta_dict[int(split_line[1])][-1]) == list: 
                    fasta_dict[int(split_line[1])].append([split_line[4]])
                else:
                    fasta_dict[int(split_line[1])][-1].append(split_line[4])
    for key, value in fasta_dict.items():
        if not len(value) == 5:
            fasta_dict[key].append(["\\N"])
        if len(value[-1]) > 1:
            if ([value[-1][0]]*len(value[-1]) == value[-1]):
                fasta_dict[key] = value[:-1]
                fasta_dict[key].append(value[-1][0])
            else:
                fasta_dict[key] = value[:-1]
                fasta_dict[key].append("Undetermined")
        else:
            fasta_dict[key] = value[:-1]
            fasta_dict[key].append(value[-1][0])
        if value[-1] == ["\\N"]:
            fasta_dict[key] = value[:-1]
            fasta_dict[key].append("N/A")


def CAZy_reader(input_files):
    for i, x in enumerate(input_files):
        with open(x, "r") as file:
            for line in file.readlines():
                if line.startswith("#"):
                    pass
                else:
                    split_line = line.split("\t")
                    CAZy_family = ""
                    for family in (split_line[2]).split("+"):
                        CAZy_family += f'{family.split("(")[0]}, '
                    fasta_dict[int(((split_line[0]).split("|"))[2])].append(CAZy_family[:-2])
    for key, value in fasta_dict.items():
        if not len(value) == 6:
            fasta_dict[key].append("N/A")


def MEROPS_reader(input_file, header_file):
    header_dict = {}
    with open(header_file, "r") as headers:
        for line in headers.readlines():
            header_split = line.split(" - ")
            header_dict[header_split[0]] = header_split[1].replace(" \n", "")
    
    with open(input_file, "r") as file:
        for line in file.readlines():
            split_line = line.split("\t")
            if not type(fasta_dict[int(((split_line[0]).split("|"))[2])][-1]) == list: 
                fasta_dict[int(((split_line[0]).split("|"))[2])].append([header_dict[split_line[1]]])
            else:
                fasta_dict[int(((split_line[0]).split("|"))[2])][-1].append(header_dict[split_line[1]])
    for key, value in fasta_dict.items():
        if not len(value) == 7:
            fasta_dict[key].append(["N/A"])
        if len(value[-1]) > 1:
            fasta_dict[key] = value[:-1]
            fasta_dict[key].append(value[-1][0])
        else:
            fasta_dict[key] = value[:-1]
            fasta_dict[key].append(value[-1][0])


def IPR_reader(input_file, tf_pfam):
    tf_pfam_list = []
    with open(tf_pfam, "r") as pfam:
        for line in pfam.readlines():
            pfam_split = line.split("\t")
            tf_pfam_list.append(pfam_split[0])
    with open(input_file, "r") as file:
        for line in file.readlines():
            if line.startswith("#"):
                pass
            else:
                split_line = (line.split("\t"))
                if split_line[4] in tf_pfam_list:
                    if not type(fasta_dict[int(split_line[0])][-1]) == list: 
                        fasta_dict[int(split_line[0])].append([split_line[5]])
                    else:
                        fasta_dict[int(split_line[0])][-1].append(split_line[5])
                else:
                    pass
    for key, value in fasta_dict.items():
        if not len(value) == 8:
            fasta_dict[key].append(["N/A"])
        tfs = ""
        for tf in value[-1]:
            tfs += f"{tf}; "
        fasta_dict[key] = value[:-1]
        fasta_dict[key].append(tfs[:-2])


def BUSCO_reader(BUSCO_file, gff_file):
    with open(BUSCO_file, "r") as BUSCO_data:
        for line in BUSCO_data.readlines():
            split_line = line.split("\t")
            


tmhmm_list = []
cazy_list = []
for item in folder_content:
    if item.__contains__("fasta"):
        fasta_file = f"{str(sys.argv[2])}/{item}"
    if item.__contains__("SignalP"):
        signalp_file = f"{str(sys.argv[2])}/{item}"
    if item.__contains__("TMHMM"):
        tmhmm_list.append(f"{str(sys.argv[2])}/{item}")
    if item.__contains__("KEGG"):
        kegg_file = f"{str(sys.argv[2])}/{item}"
    if item.__contains__("KOG"):
        kog_file = f"{str(sys.argv[2])}/{item}"
    if item.__contains__("CAZy"):
        cazy_list.append(f"{str(sys.argv[2])}/{item}")
    if item.__contains__("MEROPS"):
        merops_file = f"{str(sys.argv[2])}/{item}"
    if item.__contains__("IPR"):
        IPR_file = f"{str(sys.argv[2])}/{item}"

for record in SeqIO.parse(fasta_file, "fasta"):
    fasta_dict[int(((record.id).split("|"))[2])] = [len(record.seq)]

SP_reader(signalp_file)
TMHMM_reader(tmhmm_list)
KEGG_reader(kegg_file)
KOG_reader(kog_file)
CAZy_reader(cazy_list)
MEROPS_reader(merops_file, "MEROPS_headers.txt")
IPR_reader(IPR_file, "TF_PFAM.tsv")

row = 0
for key, value in fasta_dict.items():
    worksheet.write(row, 0, key)
    worksheet.write(row, 1, value[0])
    worksheet.write(row, 2, value[1])
    worksheet.write(row, 3, value[2])
    worksheet.write(row, 4, value[3])
    worksheet.write(row, 5, value[4])
    worksheet.write(row, 6, value[5])
    worksheet.write(row, 7, value[6])
    worksheet.write(row, 8, value[7])
    row += 1
    
workbook.close()
