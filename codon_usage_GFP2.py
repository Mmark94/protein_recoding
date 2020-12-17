GFP = 'ATGTCTAAAGGTGAAGAATTATTCACTGGTGTTGTCCCAATTTTGGTTGAATTAGATGGTGATGTTAATGGTCACAAATTTTCTGTCTCCGGTGAAGGTGAAGGTGATGCTACTTACGGTAAATTGACCTTAAAATTTATTTGTACTACTGGTAAATTGCCAGTTCCATGGCCAACCTTAGTCACTACTTTCGGTTATGGTGTTCAATGTTTTGCTAGATACCCAGATCATATGAAACAACATGACTTTTTCAAGTCTGCCATGCCAGAAGGTTATGTTCAAGAAAGAACTATTTTTTTCAAAGATGACGGTAACTACAAGACCAGAGCTGAAGTCAAGTTTGAAGGTGATACCTTAGTTAATAGAATCGAATTAAAAGGTATTGATTTTAAAGAAGATGGTAACATTTTAGGTCACAAATTGGAATACAACTATAACTCTCACAATGTTTACATCATGGCTGACAAACAAAAGAATGGTATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGTTCTGTTCAATTAGCTGACCATTATCAACAAAATACTCCAATTGGTGATGGTCCAGTCTTGTTACCAGACAACCATTACTTATCCACTCAATCTGCCTTATCCAAAGATCCAAACGAAAAGAGAGACCACATGGTCTTGTTAGAATTTGTTACTGCTGCTGGTATTACCCATGGTATGGATGAATTGTACAAATAA'

# A function that translates DNA into a protein sequence.
def translate(DNA):
    # A python dictionary data structure to translate.
    dna_to_pro = {'ATG': 'M', 'GCG': 'A', 'TCA': 'S', 'GAA': 'E', 'GGG': 'G', 'GGT': 'G', 'AAA': 'K', 'GAG': 'E', 'AAT': 'N', 'CTA': 'L',
                  'CAT': 'H', 'TCG': 'S', 'TAG': 'STOP', 'GTG': 'V', 'TAT': 'Y', 'CCT': 'P', 'ACT': 'T', 'TCC': 'S', 'CAG': 'Q', 'CCA': 'P',
                  'TAA': 'STOP', 'AGA': 'R', 'ACG': 'T', 'CAA': 'Q', 'TGT': 'C', 'GCT': 'A', 'TTC': 'F', 'AGT': 'S', 'ATA': 'I', 'TTA': 'L',
                  'CCG': 'P', 'ATC': 'I', 'TTT': 'F', 'CGT': 'R', 'TGA': 'STOP', 'GTA': 'V', 'TCT': 'S', 'CAC': 'H', 'GTT': 'V', 'GAT': 'D',
                  'CGA': 'R', 'GGA': 'G', 'GTC': 'V', 'GGC': 'G', 'TGC': 'C', 'CTG': 'L', 'CTC': 'L', 'CGC': 'R', 'CGG': 'R', 'AAC': 'N',
                  'GCC': 'A', 'ATT': 'I', 'AGG': 'R', 'GAC': 'D', 'ACC': 'T', 'AGC': 'S', 'TAC': 'Y', 'ACA': 'T', 'AAG': 'K', 'GCA': 'A',
                  'TTG': 'L', 'CCC': 'P', 'CTT': 'L', 'TGG': 'W'}
    protein = []
    start = 0
    # Step through the DNA sequence and translate.
    while start + 2 < len(DNA):
        codon = DNA[start:start + 3]
        protein.append(dna_to_pro[codon])
        start += 3
    return ''.join(protein)

def mean_list(List):
    sum = 0
    for number in List:
        sum = sum + number
    mean = sum / len(List)
    return mean

# this function calculate the mean of the codon usage in a DNA sequence
def codon_usage_mean(seq: str):
    Dic_codon_usage2 = {'TTT': 0.59, 'TTC': 0.41, 'TTA': 0.28, 'TTG': 0.29, 'CTT': 0.13, 'CTC': 0.06, 'CTA': 0.14,
                        'CTG': 0.11, 'ATT': 0.46, 'ATC': 0.26, 'ATA': 0.27, 'GTT': 0.39, 'GTC': 0.21, 'GTA': 0.21,
                        'GTG': 0.19, 'TCT': 0.26, 'TCC': 0.16, 'TCA': 0.21, 'TCG': 0.1, 'AGT': 0.16, 'AGC': 0.11,
                        'CCT': 0.31, 'CCC': 0.15, 'CCA': 0.41, 'CCG': 0.12, 'ACT': 0.35, 'ACC': 0.22, 'ACA': 0.3,
                        'ACG': 0.13, 'GCT': 0.38, 'GCC': 0.22, 'GCA': 0.29, 'GCG': 0.11, 'TAT': 0.56, 'TAC': 0.44,
                        'CAT': 0.64, 'CAC': 0.36, 'CAA': 0.69, 'CAG': 0.31, 'AAT': 0.59, 'AAC': 0.41, 'AAA': 0.58,
                        'AAG': 0.42, 'GAT': 0.65, 'GAC': 0.35, 'GAA': 0.71, 'GAG': 0.29, 'TGT': 0.63, 'TGC': 0.37,
                        'CGT': 0.15, 'CGC': 0.06, 'CGA': 0.07, 'CGG': 0.04, 'AGA': 0.48, 'AGG': 0.21, 'GGT': 0.47,
                        'GGC': 0.19, 'GGA': 0.22, 'GGG': 0.12, 'ATG': 1.0, 'TGG': 1.0, 'TAA': 0.61, 'TAG': 0.09, 'TGA': 0.3}
    codon_usage_protein=[]
    first_codon="ATG"
    start = seq.find(first_codon)
    if start == -1:
        start=0
    # Step through the DNA sequence.
    while start + 2 < len(seq):
        codon = seq[start:start + 3]
        codon_usage_protein.append(Dic_codon_usage2[codon])
        start += 3
    Mean = mean_list(codon_usage_protein)
    #print("the mean of codon usage is", Mean)
    #print("the codon usage is", codon_usage_protein)
    return Mean

# seq = DNA sequence to change. The function will change the codon1 with the codon2
def change_codon_start_Dict_score(seq: str, Dic: dict, skip=0):
    if type(Dic) != dict:
        print("Error, the input is not a dictionary")
        return 0
    print("Number of codon type to change =", len(Dic))
    if len(Dic) == 0:
        return seq
    Dic_codon_usage2 = {'TTT': 0.59, 'TTC': 0.41, 'TTA': 0.28, 'TTG': 0.29, 'CTT': 0.13, 'CTC': 0.06, 'CTA': 0.14,
                        'CTG': 0.11, 'ATT': 0.46, 'ATC': 0.26, 'ATA': 0.27, 'GTT': 0.39, 'GTC': 0.21, 'GTA': 0.21,
                        'GTG': 0.19, 'TCT': 0.26, 'TCC': 0.16, 'TCA': 0.21, 'TCG': 0.1, 'AGT': 0.16, 'AGC': 0.11,
                        'CCT': 0.31, 'CCC': 0.15, 'CCA': 0.41, 'CCG': 0.12, 'ACT': 0.35, 'ACC': 0.22, 'ACA': 0.3,
                        'ACG': 0.13, 'GCT': 0.38, 'GCC': 0.22, 'GCA': 0.29, 'GCG': 0.11, 'TAT': 0.56, 'TAC': 0.44,
                        'CAT': 0.64, 'CAC': 0.36, 'CAA': 0.69, 'CAG': 0.31, 'AAT': 0.59, 'AAC': 0.41, 'AAA': 0.58,
                        'AAG': 0.42, 'GAT': 0.65, 'GAC': 0.35, 'GAA': 0.71, 'GAG': 0.29, 'TGT': 0.63, 'TGC': 0.37,
                        'CGT': 0.15, 'CGC': 0.06, 'CGA': 0.07, 'CGG': 0.04, 'AGA': 0.48, 'AGG': 0.21, 'GGT': 0.47,
                        'GGC': 0.19, 'GGA': 0.22, 'GGG': 0.12, 'ATG': 1.0, 'TGG': 1.0, 'TAA': 0.61, 'TAG': 0.09,
                        'TGA': 0.3}
    stop_codons=["TAA","TAG","TGA"]
    new_seq = ""
    codon_changed = 0
    score = 0
    score_temp = 0
    first_codon="ATG"
    start = seq.find(first_codon)
    if start == -1:
        start=0
    # Step through the DNA sequence.

    # if skip is given, it will not change the first aa after the ATG
    start = start + skip*3
    new_seq = seq[:start]
    while start + 2 < len(seq):
        codon = seq[start:start + 3]
        # if the codon is the codon1, change it to codon2
        if codon in Dic:
            new_seq = new_seq + Dic[codon]
            #codon_changed and score
            codon_changed = codon_changed + 1
            score_temp = Dic_codon_usage2[codon] - Dic_codon_usage2[Dic[codon]]
            score = score + score_temp
        #if he finds a stop codon he will stop the changes
        elif codon in stop_codons:
            new_seq = new_seq + seq[start:]
            break
        else:
            new_seq = new_seq + codon
        start += 3

    print("codon_changed =", codon_changed)
    print("Score =", round(score, 4))
    usage_mean1 = round(codon_usage_mean(seq),4)
    usage_mean2 = round(codon_usage_mean(new_seq),4)
    print("the previous codon usage mean is ", usage_mean1, "the codon usage of the new sequence is", usage_mean2)
    translated1 = translate(seq)
    translated2 = translate(new_seq)
    if translated1 != translated2:
        print("The two translated sequences are different")
        print("sequence translation1 =", translated1)
        print("sequence translation2 =", translated2)
    return new_seq

# put the codon to change here in DIC. {"codon_to_change":"codon_to_change_with"}
DIC= {"TTA":"CTC", "TTG":"CTC"}
#print("Changed sequence: ", change_codon_start_Dict_score(GFP,DIC))

#if you want to make the changes and keep the first x aa as they are in the original sequence you can put this number as third argument.
#In this example, it will save the first 5 aa
print("Changed sequence: ", change_codon_start_Dict_score(GFP,DIC, 0))