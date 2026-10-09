# -*- coding: utf-8 -*-
"""
Created on Sat Jul 16 21:57:49 2022

@author: Wenyuan Su
"""

# initial

current_CID = '0'
CID_tag = "DTXCID"  # DTXCID
Name_tag = 'Preferred_name'  # INCHIKEY
MOLECULAR_FORMULA_tag = 'Mol_Formula'  # MOLECULAR_FORMULA
AVERAGE_MASS_tag = 'Mol_Weight' #AVERAGE_MASS
CAS_tag = 'CASRN' #CAS
PREDICTED_SPECTRUM_tag = 'PREDICTED SPECTRUM' #PREDICTED SPECTRUM


data_dict = dict()
with open("DSSTox_v2000_FINAL_Spectrum.txt") as f:
    for line in f:
        '''
        if SMILES_tag in line:
            current_CID = line.split(sep=' : ')[-1].strip('\n')
            data_dict[current_CID] = []
        '''
        if CID_tag in line:
            current_CID = line.split(sep='DTXCID')[-1].strip('\n')
            data_dict[current_CID] = []
        if Name_tag in line:
            newline = f.readline()
            data_dict[current_CID].append(newline.strip('\n'))
        if CAS_tag in line:
            newline = f.readline()
            data_dict[current_CID].append(newline.strip('\n'))
        if MOLECULAR_FORMULA_tag in line:
            newline = f.readline()
            data_dict[current_CID].append(newline.strip('\n'))
        if AVERAGE_MASS_tag in line:
            newline = f.readline()
            data_dict[current_CID].append(newline.strip('\n'))
        if PREDICTED_SPECTRUM_tag in line:
            block_list = []
            newline = 'start'
            while newline != '\n':
                newline = f.readline()
                block_list.append(newline.strip('\n'))  # 读取一段区块，指针随之下移。
            data_dict[current_CID].extend(block_list)  # SPECTRUM

RI_index = dict()
with open("PAC_Predicted RI.txt") as f:
    for line in f:
        if CID_tag in line:
            current_CID = line.split(sep='\t')[0].strip('DTXCID')
            RI_pre = line.split(sep='\t')[-1].strip('\n').strip(' ')
            data_dict[current_CID].append(RI_pre)
            
with open('DSSTox_v2000_FINAL_annotated.msp','w') as f:
    for key,element in data_dict.items():
        #f.write('CASNo: ' + str(key) + '\n') # write key(ID)
        f.write('Name: DTXCID' + str(key) + '\n')  # write key(ID)
        f.write('Synon: ' + str(element[0]) + '\n')  # 
        f.write('Retention_index: SemiStdNP=' + str(element[len(element)-1]) + '/0/1' + '\n')  # 
        f.write('Formula: ' + str(element[2])  + '\n')  # 
        f.write('ExactMass: ' + str(element[3])  + '\n')
        f.write('CAS#: ' + str(element[1])  + '\n')
        f.write('Num Peaks: ' + str(len(element)-6)  + '\n') # 减去Inchikey,Formula,ExactMass和预测二级质谱最后一行（\n），共4行
        for i in range(4,len(element)-1): 
                f.write(element[i] + '\n')
        f.write('\n')
        
