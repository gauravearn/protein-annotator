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
        extract_pattern[fasta_names[i]] = fasta_sequences[i][mRNAstart_coordinate[i]-1:mRNAend_coordinate[i]]
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

def plotCDS(inputfile, outputfile):
    """
    Author: Gaurav Sablok
    Universitat Potsdam
    Date: 2024-4-16
    a coding region plotter for the protein annotations 
    given a reference genome with the proteome aligned,
    it will plot the coding regions.
    @inputfile = file aligned to the genome using the protein hints
    @outputfile = file to which the binary layout for the protein hints should be written. 
    """
    readfile = [i for i in open(inputgff, "r").readlines() if "#" not in i]
    with open(inputgff + ".coding.gff", "w") as writegff:
        writegff.write("col0 \t col1 \t col2 \t col3 \t col4 \t col5 \t col6 \t col7 \t col8 \t col9\n")
        for line in readfile:
            writegff.write(line)
        writegff.close()
    iterator = [i.strip().split() for i in open(inputgff + ".coding.gff").readlines() if i.strip().split()[2] == "CDS"]
    codinglength = []
    for i in range(len(readfile)):
        codinglength.append(int(readfile[i].strip().split()[4])-int(readfile[i].strip().split()[3]))
    binrange_500 = []
    binrange_100 = []
    binrange_1000 = []
    binrange_2000 = []
    for i in range(len(codinglength)):
        if codinglength[i] <= 100:
            binrange_100.append(codinglength[i])
        elif codinglength[i] <= 500 and codinglength[i] >= 100:
            binrange_500.append(codinglength[i])
        elif codinglength[i] <= 1000 and codinglength[i] >= 500:
            binrange_1000.append(codinglength[i])
        elif codinglength[i] <= 2000 and codinglength[i] >= 1000:
            binrange_2000.append(codinglength[i])
        else:
            pass
    lengthbins = [len(binrange_100), len(binrange_500), len(binrange_1000), len(binrange_2000)]
    figure = sns.histplot(pd.DataFrame(lengthbins, columns = ["lengthbins"])).get_figure()
    figure.savefig("save.png")
    with open(outputfile, "w") as writebins:
        writebins.write("The length distribution for the coding regions for the less than 100bp are:")
        for i in range(len(binrange_100)):
            writebins.write(f"{binrange_100[i]}\n")
        writebins.write("The length distribution for the coding regions for the less than 500bp are:")
        for i in range(len(binrange_500)):
            writebins.write(f"{binrange_500[i]}\n")
        writebins.write("The length distribution for the coding regions for the less than 1000bp are:")
        for i in range(len(binrange_1000)):
            writebins.write(f"{binrange_1000[i]}\n")
        writebins.write("The length distribution for the coding regions for the less than 2000bp are:")
        for i in range(len(binrange_2000)):
            writebins.write(f"{binrange_2000[i]}\n")
        writebins.close()


def plotmRNAs(inputgff, outfile):
    """
    Author Gaurav Sablok
    Universitat Potsdam
    Date: 2024-4-15
    plotting the mRNA for the genome annotation and making
    the graph layout of the mRNAs using the graphql API. 
    @inputfile = file aligned to the genome using the protein hints
    @outputfile = file to which the binary layout for the protein hints should be written.
    """
    readfile = [i for i in open(inputgff, "r").readlines() if "#" not in i]
    with open(inputgff + ".code.gff", "w") as writegff:
        writegff.write("col0 \t col1 \t col2 \t col3 \t col4 \t col5 \t col6 \t col7 \t col8 \t col9\n")
        for line in readfile:
            writegff.write(line)
        writegff.close()
    read = [i.strip().split() for i in open(inputgff + ".code.gff").readlines() if i.strip().split()[2] == "mRNA"]
    mRNAplotterstart = []
    mRNAplotterend = []
    for i in range(len(read)):
        mRNAplotterstart.append(int(read[i][3]))
        mRNAplotterend.append(int(read[i][4]))
    lengthestimates = []
    for i in range(len(mRNAplotterstart)):
        lengthestimates.append(mRNAplotterend[i]-mRNAplotterstart[i])
    with open(outfile, "w") as writefile:
        writefile.write(f"These are the length estimates of the protein predicted\n")
        for i in range(len(read)):
            writefile.write(f"{mRNAplotterstart[i]}\t{mRNAplotterend[i]}\n")
        writefile.close()
    difference = []
    for i in range(len(mRNAplotterstart)):
        difference.append(mRNAplotterend[i]-mRNAplotterstart[i])
    dataframe = pd.DataFrame(difference, columns = ["mRNAlength"])
    histogram = sns.histplot(dataframe).get_figure()
    histogram.savefig("mRNAlength.png")

def generateintergenic(inputgff, fasta, outfile):
    """
    # Author: Gaurav Sablok
    # Universitat Potsdam
    # Date: 2024-4-16
    a integenic regions stichter for the genome annotations and can 
    annotated and extract all the integenic regions present by the alignment
    of the protein to the genome regions.
    generateintergenic("/home/gaurav/Desktop/final_code_push/multi.gff", 
                             "/home/gaurav/Desktop/final_code_push/multi.fasta", 
                                          "/home/gaurav/Desktop/final_code_push/final.fasta")
    """
    readfile = [i for i in open(inputgff, "r").readlines() if "#" not in i]
    with open(inputgff + ".coding.gff", "w") as writegff:
        writegff.write("col0 \t col1 \t col2 \t col3 \t col4 \t col5 \t col6 \t col7 \t col8 \t col9\n")
        for line in readfile:
            writegff.write(line)
        writegff.close()
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
    readiterator = [i.strip().split() for i in open(inputgff + ".coding.gff").readlines() if i.strip().split()[2] == "CDS"]
    ids = set([readiterator[i][0] for i in range(len(readiterator))])
    intergenic = []
    for i in range(len(readiterator)):
        if readiterator[i][0] in ids:
            intergenic.append([readiterator[i][0], readiterator[i][3], readiterator[i][4]])
        else:
            pass
    extractintergenic = []
    for i in range(len(intergenic)-1):
        for j in range(len(fasta_names)):
            for k in range(len(fasta_sequences)):
                if intergenic[i][0] == fasta_names[j]:
                    extractintergenic.append([intergenic[i][0], fasta_sequences[k][int(intergenic[i][2]):int(intergenic[i+1][1])]])
    finalwrite = []
    for i in range(len(extractintergenic)):
        if extractintergenic[i][1] == "":
            pass
        else:
            finalwrite.append(extractintergenic[i])
    with open(outfile, "w") as fastawrite:
        for i in range(len(finalwrite)):
            fastawrite.write(f">{finalwrite[i][0]}\n{finalwrite[i][1]}\n")
        fastawrite.close()
