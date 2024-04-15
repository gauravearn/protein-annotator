#!/usr/bin/python3
#Author Gaurav Sablok
# Universitat Potsdam
# Date 2024-4-3
import pandas as pd
def generatemRNAs(inputgff, fasta = None, file = None):
    """
    docstring: given a miniprot alignment and the fasta
    it will extract the aligned regions and will make a fasta
    of the same. Only the aligned regions are extracted as a 
    part of the substring.
    usage: generatingAlignments("/home/gaurav/Desktop/sample.gff")
    it optionally takes a fasta and a file name to write the corresponding
    patterns
    @inputgff = gff aligned from the miniprot
    @fasta = fasta sequences from which you want to extract the pattern
    @file = file to which you want to write the sequences of the extracted pattern
    generatingAlignments("/home/gaurav/Desktop/final_code_push/multi.gff", 
                        "/home/gaurav/Desktop/final_code_push/multi.fasta", 
                               "/home/gaurav/Desktop/final_code_push/multiout.fasta")
    """
    readfile = [i for i in open(inputgff, "r").readlines() if "#" not in i]
    with open(inputgff + ".clipped.gff", "w") as writegff:
        writegff.write("col0 \t col1 \t col2 \t col3 \t col4 \t col5 \t col6 \t col7 \t col8 \t col9\n")
        for line in readfile:
            writegff.write(line)
        writegff.close()
    with open(inputgff + ".clipped.gff", "r") as readgff: 
        readdataframe = pd.read_csv(readgff, sep = "\t")
        mRNAannotations = []
        mRNAstart_coordinate = []
        mRNAend_coordinate = []
        for i in range(len(readdataframe.iloc[:,2])):
            if readdataframe.iloc[:,2][i] == "mRNA":
                mRNAannotations.append(readdataframe.iloc[:,2][i])
        for i in range(len(readdataframe.iloc[:,3])):
            if readdataframe.iloc[:,2][i] == "mRNA":
               mRNAstart_coordinate.append(readdataframe.iloc[:,3][i])
        for i in range(len(readdataframe.iloc[:,4])):
            if readdataframe.iloc[:,2][i] == "mRNA":
               mRNAend_coordinate.append(readdataframe.iloc[:,4][i])
    read_transcripts = [i.strip() for i in open(fasta, "r").readlines()]
    fasta_transcript_dict = {}
    for i in read_transcripts:
        if i.startswith(">"):
            path = i.strip()
            if i not in fasta_transcript_dict:
                fasta_transcript_dict[i] = ""
                continue
        fasta_transcript_dict[path] += i.strip()
    fasta_sequences = list(fasta_transcript_dict.values())
    fasta_names = list(fasta_transcript_dict.keys())
    extract_pattern = {}
    for i in range(len(fasta_sequences)):
        extract_pattern[fasta_names[i]] = fasta_sequences[mRNAstart_coordinate[i]:mRNAend_coordinate[i]]
    extractkeys = list(extract_pattern.keys())
    extractvalues = list(extract_pattern.values())
    finalextractvalues = [str(''.join(extractvalues[i])) for i in range(len(list(extractvalues)))]
    with open(file, "w") as fastawrite:
        for i in range(len(extractkeys)):
            fastawrite.write(f">{extractkeys[i]}\n{finalextractvalues[i]}\n")
        fastawrite.close()


def generateCDS(inputgff, fasta, outfile):
    """
    # Author: Gaurav Sablok
    # Universitat Potsdam
    # Date: 2024-4-15
    a coding regions stichter for the genome annotations and can 
    annotated and extract all the coding regions present by the alignment
    of the protein to the genome regions. It writes a fasta and if you want
    then you can invoke this a part of the tokenziers but dont forget to add the 
    padding. such as 
    import keras 
    import sklearn 
    ** add padding for the sequence labelling. 
    """
    readfile = [i for i in open(inputgff, "r").readlines() if "#" not in i]
    with open(inputgff + ".coding.gff", "w") as writegff:
        writegff.write("col0 \t col1 \t col2 \t col3 \t col4 \t col5 \t col6 \t col7 \t col8 \t col9\n")
        for line in readfile:
            writegff.write(line)
        writegff.close()
    iterator = [i.strip().split() for i in open(inputgff + ".coding.gff").readlines() if i.strip().split()[2] == "CDS"]
    iteratorids = list(set([i[0] for i in iterator]))
    iteratorgetter = []
    for i in range(len(iterator)):
        for j in range(len(iteratorids)):
            if iteratorids[j] == str(iterator[i][0]):
                iteratorgetter.append([iterator[i][0],iterator[i][2],iterator[i][3], iterator[i][4]])
    read_transcripts = [i.strip() for i in open(fasta, "r").readlines()]
    fasta_transcript_dict = {}
    for i in read_transcripts:
        if i.startswith(">"):
            path = i.strip()
            if i not in fasta_transcript_dict:
                fasta_transcript_dict[i] = ""
                continue
        fasta_transcript_dict[path] += i.strip()
    fasta_sequences = list(fasta_transcript_dict.values())
    fasta_names = [i.replace(">", "") for i in list(fasta_transcript_dict.keys())]
    extractcoding = []
    for i in range(len(iteratorgetter)):
        for j in range(len(fasta_names)):
            for k in range(len(fasta_sequences)):
                if iteratorgetter[i][0] == fasta_names[j]:
                    extractcoding.append([iteratorgetter[i][0], fasta_sequences[k][int(iteratorgetter[i][2]):int(iteratorgetter[i][3])]])
    with open(outfile, "w") as fastawrite:
        for i in range(len(extractcoding)):
            fastawrite.write(f">{extractcoding[i][0]}\n{extractcoding[i][1]}\n")
        fastawrite.close()