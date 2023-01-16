## importing all the important libraries

import os
import numpy as np
import pandas as pd
import re
from Bio import SeqIO
import subprocess
import warnings



### Function: To fill Mutation Column. split the value in Gene column and use that to fill the Mutation Column.
def Mutation_col(final_out_file):
    counter = 0
    for i in final_out_file['Gene']:
        str_i = i.split('_')[1]
        final_out_file['Mutation'][counter] = str_i
        counter += 1
    return final_out_file

######### Start, End, Orientation Column #######

### Function: Extract start, end, orientation from gff using below function. workflow to fill these column is given below
### extact all the gene rows and specific protein rows from the gff file and save it to single dataframe
###  use this dataframe obtained from gff file to fill columns in main output file

def extract_gff(file):
  gene_info = pd.DataFrame(columns=['Gene','Start', 'End', 'Orientation'])
  gff_df = pd.read_csv(file, sep = '\t', engine='python', header=0, skiprows = [1, 2, 3, 4, 5,6])
  gff_df = gff_df.reset_index()
  gff_df.columns = ['ID', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attribute']
  gff_df.dropna(inplace=True)
  gff_df["Strand"] = np.where(gff_df["Strand"] == "+", 'Forward', 'Backward')

  gff_df_gene = gff_df[gff_df.Attribute.str.contains("gene")]
  counter = 0
  for i in gff_df_gene['Attribute']:
    words = i.split(sep = ';')
    x_gene_ls = ['gene=']
    res = [sub for sub in words
          if all((ele in sub) for ele in x_gene_ls)]
    if len(res) == 1:
      gene_name = res[0].split(sep = '=')[1]
      gene_info = gene_info.append({'Gene': gene_name, 'Start': gff_df_gene['Start'].values[counter], 'End': gff_df_gene['End'].values[counter],
                                    'Orientation': gff_df_gene['Strand'].values[counter]}, ignore_index = True)
    counter += 1
  gene_info = gene_info.drop_duplicates()
  final_genefile = pd.DataFrame(columns = ['Gene','Start', 'End', 'Orientation'])
  gene = []
  for i in gene_info['Gene']:
    if i not in gene:
      gene.append(i)
  for i in gene_info['Gene']:
    df1 = gene_info[gene_info['Gene'] == i]
    df1['Max_diff']  = df1['End'] - df1['Start']
    df1 = df1[df1['Max_diff']==df1['Max_diff'].max()]
    df1 = df1.drop(['Max_diff'], axis=1)
    final_genefile = final_genefile.append(df1)
    
  final_genefile = final_genefile.drop_duplicates()
  final_genefile = final_genefile.reset_index(drop=True)

  return final_genefile

#### Function: TO Extract Start, End, orientation for Proteins
def prot_gff(file):
    prot_info = pd.DataFrame(columns=['Gene','Start', 'End', 'Orientation'])
    gff_df = pd.read_csv(file, sep = '\t', engine='python', header=0, skiprows = [1, 2, 3, 4, 5,6])
    prot_dict = {"proteinase":"protease", "reverse transcriptase":"RT", "integrase":"integrase", "gene=env":"gp160", "product=capsid":"capsid"}  ## line specific for HIV GFF
    word_gff = list(prot_dict.keys())
    gff_df = gff_df.reset_index()
    gff_df.columns = ['ID', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attribute']
    gff_df.dropna(inplace=True)
    gff_df["Strand"] = np.where(gff_df["Strand"] == "+", 'Forward', 'Backward')
    
    gff_df['Pattern'] = ''
    counter = 0
    for i in gff_df['Attribute']:
        wor_lsd = []
        for j in word_gff:
            if j in i:
                wor_lsd.append(j)
        gff_df['Pattern'][counter] = wor_lsd    
        counter += 1

    sliced_data = gff_df[gff_df['Pattern'].apply(lambda x: len(x)) == 1]
    prot_data = sliced_data.copy(deep=True)
    prot_data = prot_data.reset_index(drop=True)
    prot_data = sliced_data.copy(deep=True)
    prot_data = prot_data.reset_index(drop=True)
    prot_data['Gene'] = ''
    counter = 0
    for i in prot_data['Pattern']:
        prot_data['Gene'][counter] = i[0]
        counter += 1

    prot_data = prot_data.drop(['Pattern'], axis=1)
    prot_data['Gene'].replace(prot_dict, inplace=True)

    counter = 0
    for i in prot_data['Gene']:
        prot_info = prot_info.append({'Gene': i, 'Start': prot_data['Start'].values[counter], 'End': prot_data['End'].values[counter],
                                            'Orientation': prot_data['Strand'].values[counter]}, ignore_index = True)
        counter += 1

    final_protfile = pd.DataFrame(columns = ['Gene','Start', 'End', 'Orientation'])
    for i in prot_info['Gene']:
        df1 = prot_info[prot_info['Gene'] == i]
        df1['Max_diff']  = df1['End'] - df1['Start']
        df1 = df1[df1['Max_diff']==df1['Max_diff'].max()]
        df1 = df1.drop(['Max_diff'], axis=1)
        final_protfile = final_protfile.append(df1)
    final_protfile = final_protfile.drop_duplicates()
    final_protfile = final_protfile.reset_index(drop=True)
    return final_protfile

#### Function: To fill Start, End, Orientation Column
def St_end_ori(final_out_dat, gf_file):
    full_gene_file = extract_gff(gf_file)
    full_prot_file = prot_gff(gf_file)
    full_gene_file  = pd.concat([full_gene_file, full_prot_file], axis=0)
    full_gene_file = full_gene_file.reset_index(drop  =True)
    counter = 0
    for i in final_out_dat['Gene']:
        str_i = i.split('_')[0]
        for j in full_gene_file['Gene']:
            if j.lower() == str_i.lower():
                df_j = full_gene_file.loc[full_gene_file['Gene'] == j]
                final_out_dat['Start'][counter] = df_j.iloc[0, 1]
                final_out_dat['End'][counter] = df_j.iloc[0, 2]
                final_out_dat['Orientation'][counter] = df_j.iloc[0, 3]
        counter += 1
    return final_out_dat

######### Start, End, Orientation Column #######

#### Function: To fill Type of Mutation Column: check the Mutation column and fill the type of mutation according to the case of the letter
def Type_of_Mut(fin_out):
    counter = 0
    indel_ls = ['insertion', 'deletion']
    for i in fin_out['Mutation']:
        res = [ele for ele in indel_ls if(ele in i)]
        if len(res) > 0:
            fin_out['Type_of_Mutation'][counter] = 'indel'
        else:
            str_mut = [x for x in i]
            if str_mut[-1].islower() == True:
                fin_out['Type_of_Mutation'][counter] = 'Nucloetide Change'
            else:
                fin_out['Type_of_Mutation'][counter] = 'Codon Change'
        counter += 1
    return fin_out


### Function: To fill Codon Position Column: Take the codon position value from the Mutation column based on type of mutation.
def Codon_Pos_col(fin_dat):
    counter = 0
    for i in fin_dat['Type_of_Mutation']:
        if i == 'Nucleotide Change' or i == 'indel':
            fin_dat['Codon_Position'] = '-'
        else:
            mut_val = fin_dat['Mutation'][counter]                                  
            x = re.findall('[0-9]+', mut_val)
            fin_dat['Codon_Position'][counter] = x[0]
        counter += 1
    return fin_dat

### Function: To fill Codon Middle position Column: use the formula to fill the codon middle position column based on type of mutation.
def Codon_Mid_Col(final_dat_file):
    counter = 0 
    for i in final_dat_file['Type_of_Mutation']:
        if i == 'Nucleotide Change' or i == 'indel':
            final_dat_file['Codon_Middle_Position'][counter] = '-'
        else:
            if final_dat_file['Orientation'].values[counter] == 'Forward':
                final_dat_file['Codon_Middle_Position'][counter] = int(final_dat_file['Start'].values[counter]) + (int(final_dat_file['Codon_Position'].values[counter]) * 3) - 2
            else:
                final_dat_file['Codon_Middle_Position'][counter] = int(final_dat_file['End'].values[counter]) + (int(final_dat_file['Codon_Position'].values[counter]) * 3) + 2
        counter += 1
    return final_dat_file


##################### All Functions to fill Original Codon and Codon Change Column  ######################################

## There are multiple functions to fill these columns. Workflow if same as you proceed down the code.

#*# Function: To convert gff to fasta #*#  Using shell script to convert gff to fasta
def gff2fasta(gff_name, fasta_name, output_name):
    subprocess.run(["./gtf2genes_sh.sh", f"{gff_name}", f"{fasta_name}", f"{output_name}"])
    req_fa_file = f"{output_name}_gene.fa"
    return req_fa_file

### Function: To Extract data from the FASTA file: extracting location and sequence using Bio module of python in different lists
def sequence_extract_fasta(fasta_files):
    # Defining empty list for the Fasta id and fasta sequence variables
    fasta_id = []
    fasta_seq = []

    # opening given fasta file using the file path
    with open(fasta_files, 'r') as fasta_file:
        # extracting multiple data in single fasta file using biopython
        for record in SeqIO.parse(fasta_file, 'fasta'):  # (file handle, file format)
        
            # appending extracted fasta data to empty lists variables
            fasta_seq.append(record.seq)
            fasta_id.append(record.id)
    return fasta_id, fasta_seq


### Function: Split location from fasta identifier obtained from fasta files
def split_fasta_loc(fasta_loc, empty_ls):
    for i in fasta_loc:
        word = i.split(':')[1]
        empty_ls.append(word)
    return empty_ls


### Function: Convert Sequence object to string to get nucleotide sequence
def get_nuc_seq(fasta_seq, empty_seq_ls):
    for i in fasta_seq:
        x = str(i)
        empty_seq_ls.append(x)
    return empty_seq_ls

##### Function: To get all the locations from final output file (from start and end column) and did minus 1 because in the fasta file we obtained. The sequence starts from 1 less then the actual
def get_loc_final(final_file, empty_ls_final_loc):
    for i in range(len(final_file)):
        loc_value = f"{int(final_file['Start'].values[i]) - 1}-{int(final_file['End'].values[i])}"             
        empty_ls_final_loc.append(loc_value)
    empty_ls_final_loc = [*set(empty_ls_final_loc)]
    return empty_ls_final_loc


#### Function: To Translate Nucleotide sequnence to amino acid sequence
def translate_frameshifted(sequence):
  temp_df = pd.DataFrame(columns = ['codon', 'amino'])

  gencode = {
      "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
       "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
       "UAU":"Y", "UAC":"Y", "UAA":"!", "UAG":"!",
       "UGU":"C", "UGC":"C", "UGA":"!", "UGG":"W",
       "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
       "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
       "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
       "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

  for i in range(len(sequence) // 3):
    codon_value = sequence[3*i:3*i+3]
    temp_df = temp_df.append({'codon': codon_value, 'amino': gencode[codon_value]}, ignore_index=True)
  return temp_df


#### Function: If type of mutation is Codon Change
# def if_codon_change(fin_output, list_locat_fin, seqloc_dict):
#     for k in list_locat_fin:
#         if k in seqloc_dict:
#             nu_seq = seqloc_dict[k]
#             df_amino = translate_frameshifted(nu_seq)
#             df_amino['Codon_Middle_Position'] = ''
#             start = int(k.split('-')[0])
#             stop = int(k.split('-')[1])
#             counter = 0
#             for l in range(start + 1, stop + 1, 3):
#                 df_amino['Codon_Middle_Position'][counter] = l+1
#                 counter += 1
#             for m in df_amino['Codon_Middle_Position']:
#                 gencode = {
#                 "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
#        "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
#        "UAU":"Y", "UAC":"Y", "UAA":"!", "UAG":"!",
#        "UGU":"C", "UGC":"C", "UGA":"!", "UGG":"W",
#        "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
#        "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
#        "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
#        "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
#        "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
#        "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
#        "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
#        "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
#        "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
#        "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
#        "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
#        "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
#                 if m in fin_output['Codon_Middle_Position'].values:
#                     idvalue = fin_output[fin_output['Codon_Middle_Position']==m].index.values
#                     idvalue_ls = idvalue.tolist()
#                     for n in idvalue_ls:
#                         ama = fin_output['Mutation'].values[n]
#                         aa = [s for s in ama]
#                         Originalcod_df = df_amino.loc[df_amino['Codon_Middle_Position'] == m, 'codon']
#                         Originalcod_df = Originalcod_df.reset_index(drop=True)
#                         Originalcod = Originalcod_df[0]
#                         fin_output['Original_Codon'][n] = Originalcod
#                         codon_list = [u for u,v in gencode.items() if v == aa[-1]]
#                         fin_output['Codon_Change'][n] = codon_list
#     return fin_output

#####################################################################################################

#### Function: To Fill Original Nucleotide and Nucleotide Change Column & ### Function: To fill Nucleotide Position Column
### to fill these columns, we have to use original codon and codon change column values. 
# #Match all the values and find the difference to fill the Original nucleotide and Nucleotide change value and use codon middle position to fill the nucleotide posiiton value.
def OriC_Cod_Col(final_op, user_gff, fasta_file, output_fastaname):
    final_fasta_filename = gff2fasta(user_gff, fasta_file, output_fastaname)
    counter_1 = 0
    for t in final_op['Type_of_Mutation']:
        if t == 'Codon Change':
            fasta_loc_ls, fasta_seq_ls = sequence_extract_fasta(final_fasta_filename)
            all_loc_ls = []
            for i in fasta_loc_ls:
                word = i.split(':')[1]
                all_loc_ls.append(word)
            all_seq_ls = []
            for i in fasta_seq_ls:
                x = str(i)
                all_seq_ls.append(x)
            loc_seq_dict = {all_loc_ls[i]: all_seq_ls[i] for i in range(len(all_loc_ls))}
            list_loc_finaldata = []
            for i in range(len(final_op)):
                loc_value = f"{int(final_op['Start'].values[i]) - 1}-{int(final_op['End'].values[i])}"
                list_loc_finaldata.append(loc_value)
            list_loc_finaldata = [*set(list_loc_finaldata)]

            for y in list_loc_finaldata:
                if y in loc_seq_dict:
                    nu_seq = loc_seq_dict[y]
                    df_amino = translate_frameshifted(nu_seq)      #### translate sequence
                    df_amino['Codon_Middle_Position'] = ''
                    start = int(y.split('-')[0])
                    stop = int(y.split('-')[1])
                    counter = 0                                                ##### added counter here
                    for l in range(start + 1, stop +1, 3):                    #####code changed: add 1 to start to match middle position and changed +3 to +1 in stop
                        df_amino['Codon_Middle_Position'][counter] = l+1
                        counter += 1
                    for m in df_amino['Codon_Middle_Position']:
                        gencode = {
                        "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
       "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
       "UAU":"Y", "UAC":"Y", "UAA":"!", "UAG":"!",
       "UGU":"C", "UGC":"C", "UGA":"!", "UGG":"W",
       "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
       "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
       "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
       "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
                        if m in final_op['Codon_Middle_Position'].values:
                            idvalue = final_op[final_op['Codon_Middle_Position']==m].index.values
                            idvalue_ls = idvalue.tolist()
                            for j in idvalue_ls:
                                ama = final_op['Mutation'].values[j]
                                aa = [s for s in ama]
                                Originalcod_df = df_amino.loc[df_amino['Codon_Middle_Position'] == m, 'codon']
                                Originalcod_df = Originalcod_df.reset_index(drop=True)
                                Originalcod = Originalcod_df[0]
                                final_op['Original_Codon'][j] = Originalcod         ### changed again
                                codon_list = [k for k,v in gencode.items() if v == aa[-1]]
                                final_op['Codon_Change'][j] = codon_list
        else:
            final_op['Original_Codon'][counter_1] = '-'
            final_op['Codon_Change'][counter_1] = '-'
        counter_1 += 1
    return final_op    



#####################################################################################################

#### Function: To Fill Original Nucleotide and Nucleotide Change Column & ### Function: To fill Nucleotide Position Column
def OriNuc_NucChan(output_file):
    counter = 0
    for i in output_file['Type_of_Mutation']:
        if i == 'Nucleotide Change' or i == 'indel':
            mut_string = [r for r in output_file['Mutation'][counter]]
            mut_ori_val = mut_string[0]
            mut_ch_val = mut_string[-1]
            output_file['Original_Nucleotide'][counter] = mut_ori_val
            output_file['Nucleotide_Change'][counter] = mut_ch_val
            xy = re.findall('[0-9]+', output_file['Mutation'].values[counter])
            if output_file['Orientation'].values[counter] == 'Forward':
                output_file['Nucleotide_Position'][counter] = int(output_file['Start'].values[counter]) + int(xy[0]) - 1
            else:
                output_file['Nucleotide_Position'][counter] = int(output_file['End'].values[counter]) + int(xy[0]) + 1
        elif i == 'Codon Change':
            Ori_Cod_val = output_file['Original_Codon'][counter]
            A = Ori_Cod_val[0:2]
            A_rest = Ori_Cod_val[2]
            B = Ori_Cod_val[1:]
            B_rest = Ori_Cod_val[0]
            C = Ori_Cod_val[0]+Ori_Cod_val[2]
            C_rest = Ori_Cod_val[1]
            ls_cod_ch = output_file['Codon_Change'][counter]
            ls_cod_ch = [ls_cod_ch.strip() for ls_cod_ch in eval(str(ls_cod_ch))]
            ori_nuc = []
            nuc_chan = []
            nuc_pos_val = []
            cod_mid_val = output_file['Codon_Middle_Position'][counter]
            for j in ls_cod_ch:
                D = j[0:2]
                D_rest = j[2]
                E = j[1:]
                E_rest = j[0]
                F = j[0] + j[2]
                F_rest = j[1]
                if A == D:
                    if A_rest not in ori_nuc:
                        ori_nuc.append(A_rest)
                    nuc_chan.append(D_rest)
                    nuc_pos_val.append(int(cod_mid_val) + 1)
                elif B == E:
                    if B_rest not in ori_nuc:
                        ori_nuc.append(B_rest)
                    nuc_chan.append(E_rest)
                    nuc_pos_val.append(int(cod_mid_val) - 1)
                elif C == F:
                    if C_rest not in ori_nuc:
                        ori_nuc.append(C_rest)
                    nuc_chan.append(F_rest)
                    nuc_pos_val.append(int(cod_mid_val))
        output_file['Original_Nucleotide'][counter] = ori_nuc
        output_file['Nucleotide_Change'][counter] = nuc_chan
        #### Nucleotide Position ### 
        output_file['Nucleotide_Position'][counter] = nuc_pos_val
        counter += 1
    return output_file


#### Main Function :: for full catalgue. This is the main function which should be called and execute and within these all the other functions are called.
def catalogue_ready(input_file, gff_file, fasta_file, output_fasta_filename):
    ### all the required data ###
    df_user = pd.read_excel(input_file)

    print("Automation Started......")
    final_output = df_user.copy(deep = True) ### final user output file

    ### add all the required columns ###
    final_output['Mutation'] = ''
    final_output['Start'] = ''
    final_output['End'] = ''
    final_output['Orientation'] = ''
    final_output['Type_of_Mutation'] = ''
    final_output['Codon_Position'] = ''
    final_output['Codon_Middle_Position'] = ''
    final_output['Original_Codon'] = ''
    final_output['Codon_Change'] = ''
    final_output['Original_Nucleotide'] = ''
    final_output['Nucleotide_Change'] = ''
    final_output['Nucleotide_Position'] = ''

    ### add values to mutation column ###
    final_output = Mutation_col(final_output)

    ### extract data about every gene from gff file like start, end , orientation

    print("10% Completed......")
    ### Fill Start, End, orientation column ###
    final_output = St_end_ori(final_output, gff_file)                            

    ### Fill Type of Mutation column ####
    final_output = Type_of_Mut(final_output)                                                   
    print("20% Completed......")

    ### Fill Codon Position Column
    final_output = Codon_Pos_col(final_output)
    
    print("30% Completed......")
    ### Fill Codon Middle Position Column
    final_output = Codon_Mid_Col(final_output)                                           ####-----check- replace 'Stop' with 'End' --done

    print("50% Completed......")
    ### Fill Original Codon And Codon Change Column
    final_output = OriC_Cod_Col(final_output, gff_file, fasta_file, output_fasta_filename)

    print("70% Completed......")
    #### Fill Original Nucleotide and Nucleotide Change
    final_output = OriNuc_NucChan(final_output)



    #### Final Output to excel for user
    final_output.to_excel("rawRNA_result.xlsx", index=False, header=True)
    print("95% Completed......Almost done")
    file_exists = os.path.exists("rawRNA_result.xlsx")          ## checking if the file exist
    if file_exists == True:
        final_output = pd.read_excel("rawRNA_result.xlsx")
        final_output = final_output.replace(r"[][()']", "", regex=True)           #### removing all the unwanted things from the final file
        final_output = final_output.replace(r",", " ", regex=True)
        final_output['Original_Nucleotide'] = final_output['Original_Nucleotide'].str.upper()         ### changing the case of the values
        final_output['Nucleotide_Change'] = final_output['Nucleotide_Change'].str.upper()

        ## after this code:  for filling the missing values in Original nucleotide and nucleotide change column
        miss_data = final_output[final_output['Nucleotide_Change'].isna()]
        miss_data = final_output.reset_index(drop=True)
        counter = 0
        for i,j in zip(miss_data.Original_Codon, miss_data.Codon_Change):
            oc = i
            ccv = j
            ccv = ccv.split('  ')
            cod_mid = int(miss_data['Codon_Middle_Position'][counter])
            o_nuc_v = []
            nuc_ch_v = []
            nuc_pos_v = []
            # print(type(i))
            for k in ccv:
                lst = [t for t in range(len(oc)) if oc[t] != k[t]]
                for d in lst:
                    o_nuc_v.append(oc[d])
                    nuc_ch_v.append(k[d])
                    if d == 0:
                        val_nuc = cod_mid - 1
                        nuc_pos_v.append(val_nuc)
                    elif d == 1:
                        val_nuc = cod_mid
                        nuc_pos_v.append(val_nuc)
                    elif d == 2:
                        val_nuc = cod_mid + 1
                        nuc_pos_v.append(val_nuc)
            temp_df = miss_data.iloc[[counter]]
            temp_df = pd.concat([temp_df]*len(o_nuc_v), ignore_index = True)

            ## creating temp dataframe and appending it in main dataframe
            for b in range(len(o_nuc_v)):
                temp_df['Original_Nucleotide'][b]  = o_nuc_v[b]
                temp_df['Nucleotide_Change'][b]  = nuc_ch_v[b]
                temp_df['Nucleotide_Position'][b] = nuc_pos_v[b]
            final_output = pd.concat([final_output, temp_df], axis=0, ignore_index=True)
            counter += 1
    final_output.to_excel("rawresult_RNA.xlsx", index=False, header=True)
    data = pd.read_excel("rawresult_RNA.xlsx")

    ### drop all the duplicates and drop NAN values from the final dataframe
    final_data2 = data.drop_duplicates()
    final_data2 = final_data2.dropna(subset=['Original_Nucleotide'])
    final_data2.to_excel("final_RNA_result.xlsx", index=False, header=True)
    print("100% Completed. Thank you for waiting. You can see the Output in the excel file named final_result.")

if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    input_file = input("Please give your Input file path: ")
    gff_file = input("Please provide path for refernce gff file to map input file: ")
    fasta_file = input("Please provide path for fasta file to map with gff file: ")
    output_fasta_file_name = input("Please provide File name for bedtool's output: ")
    catalogue_ready(input_file, gff_file, fasta_file, output_fasta_file_name)
    
    
    