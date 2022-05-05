####################################################################################################################################

# 1 Import necessary modules
import Bio
from Bio import SeqIO
import csv

####################################################################################################################################

# 2 Upload fasta file with alligned sequences as fasta_file; and clinvar mutations as variants
print("Path to fasta file")
fasta_file=str(input()) # Multiple_alignment.fa, first sequence should be Homo sapiens
print("Path to file with variants in format - NP_006209.2 1 M V")
variants=str(input()) # PIK3CA_ClinVar.txt - txt-file with mutations in special format (e.g. NP_006209.2 1 M V)

#####################################################################################################################################

# 3 Create list of all sequences as seqen
seqen=[]
for seq_record in SeqIO.parse(fasta_file,"fasta"):
    seqen.append(str(seq_record.seq))
    
#####################################################################################################################################

# 4 Create list of all variants as clean_clinvar
clinvar=[]
with open(variants) as f:
    for line in f:
        clinvar.append(str(line))       
clean_clinvar=[]
for i in clinvar:
    new=i.strip("\n")
    clean_clinvar.append(new.split(" "))

#######################################################################################################################################

# 5 Code

header=['Input', 'Prediction', 'Conservative', 'Substitutions', 'Uncertain_substitutions']

with open('clinvar_predictions.csv', 'w', encoding='UTF8') as f:
    writer = csv.writer(f)
    # write the header
    writer.writerow(header)

    for x in clean_clinvar: 
        # Positions
        pos=int(x[1])-1 # position of interest in alligned sequences
        pos_bfr=pos-1 
        pos_aft=pos+1

        # Human positions
        human=seqen[0] # human sequence
        human_pos=seqen[0][pos]
        if pos==0:
            hpos_bfr=0 # special case for the fist amino acid
            hpos_aft=seqen[0][pos+1]
        elif pos==lenght:
            hpos_bfr=seqen[0][pos-1]# special case for the last amino acid
            hpos_aft=0
        else:
            hpos_bfr=seqen[0][pos-1]
            hpos_aft=seqen[0][pos+1]

        # Gaps and X-amino acids
        gap="-"
        X="X"

        # Counters
        count=0 #counts how many amino acids identical to human are present in pos (pos_bfr = hpos_bfr and pos_aft = hpos_aft)
        uncount=0 #counts how many amino acids identical to human are present in pos when pos_bfr != hpos_bfr or pos_aft != hpos_aft
        uncer=0 #counts how many amino acids different from human are present in pos when pos_bfr != hpos_bfr or pos_aft != hpos_aft
        X_count=0
        gap_count=0

        # Lists of amino acids
        sub_list=[] # will collect amino acids when pos_bfr = hpos_bfr and pos_aft = hpos_aft
        uncert_list=[] # will collect amino acids when pos when pos_bfr != hpos_bfr or pos_aft != hpos_aft
        # Human lists
        human_list=[human_pos,gap,X] 
        hpos_bfr_list=[hpos_bfr,gap,X]
        hpos_aft_list=[hpos_aft,gap,X]

        # Lenght
        lenght=len(human)-1
        seq_num=len(seqen)-1

        for i in seqen[1:]: # start from 1, because first(0) sequence is human

            # special case for the first amino acid           
            if pos==0:
                if i[pos] in human_list:
                    # conservative amino acid in concervative fragment
                    if i[pos_aft] in hpos_aft_list:
                        count+=1
                    elif i[pos_aft] not in hpos_aft_list:
                        uncount+=1

                elif i[pos] not in human_list:
                    # not conservative amino acid in conservative fragment
                    if i[pos_aft] in hpos_aft_list:
                        sub_list.append(i[pos]) # Add amino acid to the list of substitution

                    # not conservative amino acid in not conservative fragment   
                    elif i[pos_aft] not in hpos_aft_list:
                        uncer+=1 
                        uncert_list.append(i[pos]) # Add amino acid to the list of uncertain prediction

            # special case for the last amino acid           
            if pos==lenght:
                if i[pos] in human_list:
                    # conservative amino acid in concervative fragment
                    if i[pos_bfr] in hpos_bfr_list :
                        count+=1
                    elif i[pos_bfr] not in hpos_bfr_list :
                        uncount+=1

                elif i[pos] not in human_list:
                    # not conservative amino acid in conservative fragment
                    if i[pos_bfr] in hpos_bfr_list:
                        sub_list.append(i[pos]) # Add amino acid to the list of substitution

                    # not conservative amino acid in not conservative fragment   
                    elif i[pos_bfr] not in hpos_bfr_list:
                        uncer+=1 
                        uncert_list.append(i[pos]) # Add amino acid to the list of uncertain prediction

            # all other    
            elif pos in range (1,lenght-1):   
                if i[pos] in human_list:
                    # conservative amino acid in concervative fragment
                    if (i[pos_bfr] in hpos_bfr_list) and (i[pos_aft] in hpos_aft_list):
                        count+=1
                    elif (i[pos_bfr] not in hpos_bfr_list) or (i[pos_aft] not in hpos_aft_list):
                        uncount+=1

                elif i[pos] not in human_list:
                    # not conservative amino acid in conservative fragment
                    if (i[pos_bfr] in hpos_bfr_list and i[pos_aft] in hpos_aft_list):
                        sub_list.append(i[pos]) # Add amino acid to the list of substitution

                    # not conservative amino acid in not conservative fragment   
                    elif (i[pos_bfr] not in hpos_bfr_list) or (i[pos_aft] not in hpos_aft_list):
                        uncer+=1 
                        uncert_list.append(i[pos]) # Add amino acid to the list of uncertain prediction

#######################################################################################################################
# 5 Output with predictions is writing in the csv file "clinvar_predictions.csv" with header

        total=count+uncount # count the number of identical to human amino acids in pos

        # Editing the lists
        sub_list=list(set(sub_list))
        sub_list=(",".join(sub_list))
        sub_list=sub_list.replace(","," ")

        uncert_list=list(set(uncert_list))
        uncert_list=(",".join(uncert_list))
        uncert_list=uncert_list.replace(","," ")

        x=(",".join(x))
        x=x.replace(","," ")

        # Writing the results and predictions into csv file clinvar_predictions.csv
        if seq_num-total<=1 and uncer==0:
            data=(x,'Damaging',total, 0, 0)
            writer.writerow(data)
        elif seq_num-total>1 and x[-1] in sub_list:
            data=(x,'Benign',total, sub_list, uncert_list)
            writer.writerow(data)
        elif uncer >=0 and x[-1] in uncert_list:# x[-1]- substitute amino acid
            data=(x,'Uncertain significance', total, sub_list, uncert_list)
            writer.writerow(data)
        elif seq_num-total>1 and x[-1] not in sub_list:
            data=(x,'Damaging', total, sub_list, uncert_list)
            writer.writerow(data)

                
