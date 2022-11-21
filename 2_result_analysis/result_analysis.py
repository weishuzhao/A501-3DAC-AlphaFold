# Requirements and dependencies, defined some functions for PyMOL analysis
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
from pymol import cmd, CmdException
import os,time
from rdkit import Chem

# import some structure analysis functions
from structure_analysis import *


# all possible pdb names
model_names = ["relaxed_model_%d_ptm_pred_0"%i for i in range(1,6)] + ["unrelaxed_model_%d_ptm_pred_0"%i for i in range(1,6)]
model_names_short = ["relaxed_model_%d_ptm"%i for i in range(1,6)] + ["unrelaxed_model_%d_ptm"%i for i in range(1,6)]
unrelaxed_model_names = ["unrelaxed_model_%d_ptm_pred_0"%i for i in range(1,6)]

part = 11
total_protein = 47
# start_id = 500*(part-1)+1
start_id = 500*(3-1)+1+353
end_id = start_id + total_protein




# start of time
start_time = time.time()

# create dataframe object for A501 proteome
df_A501 = pd.DataFrame(columns=["ID","name","sequence"])

# if you used the full proteome, you can change the total_protein number
df_A501["ID"] = ["A501_%d"%(i) for i in range(start_id,end_id)]


# Feature in this step: sequence, structure available, sequence length
structure_avail_list = []
sequence_list = []
sequence_len_list = []

for filename in df_A501["ID"]:
    structure_avail = min([os.path.exists("A501_summary/%s/%s.pdb"%(filename,model_name)) for model_name in model_names])
    if structure_avail:
        structure_avail_list.append(True)    # if all 5 models are available, then True
        cmd.load("A501_summary/%s/relaxed_model_1_ptm_pred_0.pdb"%filename,"A501")  # load the structure
        sequence = ''.join(cmd.get_fastastr("all").split('\n')[1:])     # get sequence
        sequence_list.append(sequence)
        sequence_len_list.append(len(sequence))
        cmd.reinitialize()
    else:
        structure_avail_list.append(False)
        sequence_list.append(None)
        sequence_len_list.append(None)


# add the sequence and structure availability to the dataframe
df_A501["sequence"] = sequence_list
df_A501["sequence_length"] = sequence_len_list
df_A501["structure_avail"] = structure_avail_list


# Feature in this step: plddt for all 5 unrelaxed models (relaxed models' plddt is equal to unrelaxed models)
template_list_plddt = []
for index in df_A501["ID"]:
    if os.path.exists("A501_summary/%s/ranking_debug.json"%index):
        f_json = open("A501_summary/%s/ranking_debug.json"%index,"r")
        json_dict = json.load(f_json)
        
        # read plddt from json file
        plddt_1 = json_dict["plddts"]["model_1_ptm_pred_0"]
        plddt_2 = json_dict["plddts"]["model_2_ptm_pred_0"]
        plddt_3 = json_dict["plddts"]["model_3_ptm_pred_0"]
        plddt_4 = json_dict["plddts"]["model_4_ptm_pred_0"]
        plddt_5 = json_dict["plddts"]["model_5_ptm_pred_0"]
        
        template_list_plddt.append([index,plddt_1,plddt_2,plddt_3,plddt_4,plddt_5])

    else:
        template_list_plddt.append([index,None,None,None,None,None])

# create a dataframe for plddt
template_df_plddt = pd.DataFrame(template_list_plddt,columns=["ID","model_1_ptm_avg_plddt","model_2_ptm_avg_plddt","model_3_ptm_avg_plddt","model_4_ptm_avg_plddt","model_5_ptm_avg_plddt"])
df_A501 = df_A501.merge(template_df_plddt,how="left",on="ID")


# Feature in this step: number of disulfide bonds, hydrogen bonds, salt bridges, relative surface area, secondary structure
# Feature in this step: chirality of Ca, phi, psi, omega
# Calculate protein properties, this may take a while
disu_overall_list = []    # disulfide bonds
hbon_overall_list = []    # hydrogen bonds
sbrg_overall_list = []    # salt bridges
surf_overall_list = []    # relative surface area
sasa_overall_list = []    # solvent accessible surface area
sstr_overall_list = []    # secondary structure
chir_overall_list = []    # chirality, central carbon
phiv_overall_list = []    # phi angle
psiv_overall_list = []    # psi angle
omgv_overall_list = []    # omega angle, current residue to next residue
pldt_overall_list = []    # plddt

for i in range(df_A501.shape[0]):
    disu_list = []
    hbon_list = []
    sbrg_list = []
    surf_list = []
    sasa_list = []
    sstr_list = []
    chir_list = []
    phiv_list = []
    psiv_list = []
    omgv_list = []
    pldt_list = []
    if df_A501["structure_avail"][i]:
        for model_name in model_names:
            cmd.load("A501_summary/%s/%s.pdb"%(df_A501["ID"][i],model_name),"%s"%df_A501["ID"][i])   # load pdb file
    
            cmd.distance("disul","name SG","name SG",cutoff=2.1)     # disulfide bond
            cmd.distance("hbond","all","all",mode=2)                 # hydrogen bond
            cmd.distance("saltbridge","(resn ASP and name OD2) or (resn GLU and name OD2)","(resn ARG and (name NH1 or name NH2)) or (resn LYS and name NZ) or (resn HIS and name NE2)",cutoff=4)  # salt bridge
            cmd.dss()           # secondary structure
            cmd.assign_stereo() # assign chirality
            ss_string = ""
            for a in cmd.get_model(df_A501["ID"][i] +" and n. ca").atom:
                ss_string = ss_string+a.ss 
            
            chir_protein_list = []
            phiv_protein_list = []
            psiv_protein_list = []
            omgv_protein_list = []
            pldt_protein_list = []

            chirality_space = {'chir_protein_list': chir_protein_list}
            cmd.iterate('name CA', 'chir_protein_list.append(stereo)', space=chirality_space)
            plddt_space = {'pldt_protein_list': pldt_protein_list}
            cmd.iterate('name CA', 'pldt_protein_list.append(b)', space=plddt_space)

            for res_num in range(1,int(df_A501["sequence_length"][i])+1):
                # first residue has no phi
                if res_num > 1:
                    phi = cmd.get_dihedral("%d/c"%(res_num-1), "%d/n"%res_num , "%d/ca"%res_num   , "%d/c"%res_num     )
                    phiv_protein_list.append("%.3f"%phi)
                else:
                    phiv_protein_list.append("None")

                # last residue has no psi and omega
                if res_num < df_A501["sequence_length"][i]:
                    psi = cmd.get_dihedral("%d/n"%res_num    , "%d/ca"%res_num, "%d/c"%res_num    , "%d/n"%(res_num+1) )
                    psiv_protein_list.append("%.3f"%psi)
                    omg = cmd.get_dihedral("%d/ca"%res_num   , "%d/c"%res_num , "%d/n"%(res_num+1), "%d/ca"%(res_num+1))
                    omgv_protein_list.append("%.3f"%omg)
                else:
                    psiv_protein_list.append("None")
                    omgv_protein_list.append("None")

    
            disu_list.append(len(get_raw_distances("disul")))       # disulfide bond
            hbon_list.append(len(get_raw_distances("hbond")))       # hydrogen bond
            sbrg_list.append(len(get_raw_distances("saltbridge")))  # salt bridge
            surf_list.append(cmd.get_area('all'))                   # protein surface
            cmd.set('dot_solvent',1)
            sasa_list.append(cmd.get_area('all'))                   # SASA
            sstr_list.append(ss_string)                             # secondary structure
            chir_list.append("".join(chir_protein_list))            # chirality
            phiv_list.append(",".join(phiv_protein_list))           # phi angle
            psiv_list.append(",".join(psiv_protein_list))           # psi angle
            omgv_list.append(",".join(omgv_protein_list))           # omega angle
            pldt_list.append(",".join(["%.2f"%plddt for plddt in pldt_protein_list]))           # plddt
    
            cmd.reinitialize()
    else:
        disu_list = ["None" for i in range(10)]       # disulfide bond
        hbon_list = ["None" for i in range(10)]       # hydrogen bond
        sbrg_list = ["None" for i in range(10)]       # salt bridge
        surf_list = ["None" for i in range(10)]       # protein surface
        sasa_list = ["None" for i in range(10)]       # SASA
        sstr_list = ["None" for i in range(10)]       # secondary structure
        chir_list = ["None" for i in range(10)]       # chirality
        phiv_list = ["None" for i in range(10)]       # phi angle
        psiv_list = ["None" for i in range(10)]       # psi angle
        omgv_list = ["None" for i in range(10)]       # omega angle
        pldt_list = ["None" for i in range(10)]       # plddt

    
    # add to the list
    disu_overall_list.append(disu_list)
    hbon_overall_list.append(hbon_list)
    sbrg_overall_list.append(sbrg_list)
    surf_overall_list.append(surf_list)
    sasa_overall_list.append(sasa_list)
    sstr_overall_list.append(sstr_list)
    chir_overall_list.append(chir_list)
    phiv_overall_list.append(phiv_list)
    psiv_overall_list.append(psiv_list)
    omgv_overall_list.append(omgv_list)
    pldt_overall_list.append(pldt_list)

    # os.system('clear')
    end_time = time.time()
    print("Progress: %d/%d\nTime taken: %.2f seconds"%(i+1,df_A501.shape[0],end_time-start_time),end="")
    



# Add calculated data to dataframe
for i in range(10):
    df_A501["disulfide_bond_%s"%model_names_short[i]] = [disu_overall_list[prot_index][i] for prot_index in range(df_A501.shape[0])]
for i in range(10):
    df_A501["hydrogen_bond_%s"%model_names_short[i]] = [hbon_overall_list[prot_index][i] for prot_index in range(df_A501.shape[0])]
for i in range(10):
    df_A501["salt_bridge_%s"%model_names_short[i]] = [sbrg_overall_list[prot_index][i] for prot_index in range(df_A501.shape[0])]
for i in range(10):
    df_A501["surface_%s"%model_names_short[i]] = [surf_overall_list[prot_index][i] for prot_index in range(df_A501.shape[0])]
for i in range(10):
    df_A501["SASA_%s"%model_names_short[i]] = [sasa_overall_list[prot_index][i] for prot_index in range(df_A501.shape[0])]
for i in range(10):
    df_A501["secondary_structure_%s"%model_names_short[i]] = [sstr_overall_list[prot_index][i] for prot_index in range(df_A501.shape[0])]
for i in range(10):
    df_A501["chirality_%s"%model_names_short[i]] = [chir_overall_list[prot_index][i] for prot_index in range(df_A501.shape[0])]
for i in range(10):
    df_A501["phi_angle_%s"%model_names_short[i]] = [phiv_overall_list[prot_index][i] for prot_index in range(df_A501.shape[0])]
for i in range(10):
    df_A501["psi_angle_%s"%model_names_short[i]] = [psiv_overall_list[prot_index][i] for prot_index in range(df_A501.shape[0])]
for i in range(10):
    df_A501["omega_angle_%s"%model_names_short[i]] = [omgv_overall_list[prot_index][i] for prot_index in range(df_A501.shape[0])]
for i in range(10):
    df_A501["plddt_%s"%model_names_short[i]] = [pldt_overall_list[prot_index][i] for prot_index in range(df_A501.shape[0])]



df_A501.to_csv("A501_results_part%d.tsv"%part,index=False,sep="\t")  # save A501 result to csv file
columns = ["ID","name","sequence","sequence_length","structure_avail","model_1_ptm_avg_plddt","model_2_ptm_avg_plddt","model_3_ptm_avg_plddt","model_4_ptm_avg_plddt","model_5_ptm_avg_plddt"]
# df_A501[columns].to_csv("A501_results_short.tsv",index=False,sep="\t")  # save A501 result to csv file


# end of time
end_time = time.time()
print("\nTime taken: %.2f seconds"%(end_time-start_time))
