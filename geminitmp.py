'''
┌─┐  ┌─┐  ┌┬┐  ┬  ┌┐┌  ┬
│ ┬  ├┤   │││  │  │││  │
└─┘  └─┘  ┴ ┴  ┴  ┘└┘  ┴
---------------------------------------------------
-*- coding: utf-8 -*-                              |
title            : geminitmp.py                    |
description      : one shot scripts                |
author           : dooguypapua                     |
lastmodification : 20210628                        |
version          : 0.1                             |
python_version   : 3.8.5                           |
---------------------------------------------------
'''

from geminini import *
from geminiparse import *
from geminiannot import *
from geminicluster import *
# from geminiphage import *
from geminiplot import *
import re
from datetime import datetime
from Bio import SeqIO,AlignIO
from Bio.Seq import Seq
from Bio import Entrez
import glob
from ete3 import Tree,NCBITaxa,TreeStyle,NodeStyle,TextFace
from tqdm import tqdm
import pysam
import random
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from io import StringIO
import warnings
warnings.filterwarnings("ignore")
import geminiset
import subprocess
import fileinput
import scipy.stats
import fastq as fq
from sklearn.linear_model import LinearRegression
from scipy.stats import gaussian_kde
from sqlitedict import SqliteDict
from dna_features_viewer import GraphicFeature, GraphicRecord
# os.system("rm -rf /tmp/geminitmp_*")
temp = tempfile.mkdtemp(prefix="geminitmp_", dir="/tmp")
geminiset.fctName = "geminitmp"
geminiset.title = "geminitmp"
slurmBool, cpu, memMax, memMin = get_sys_info()
geminiset.slurmBool = slurmBool
geminiset.cpu = cpu
geminiset.memMax = memMax
geminiset.memMin = memMin
geminiset.pathTMP = temp


# lstGBK = glob.glob("/mnt/g/db/vibrioDB/*/*_genomic.gbff.gz")
# pbar = tqdm(total=len(lstFiles),ncols=75,leave=False,desc="",file=sys.stdout,bar_format="  {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}")
# for pathGBK in lstGBK:
#     pathFAA = "/mnt/c/Users/dgoudenege/Works/Vibrio_FAA_reformat_contig/"+os.path.basename(pathGBK).replace("_genomic.gbff.gz",".faa")
#     gbk_to_faa(pathIN=pathGBK, pathOUT=pathFAA, syntaxic="prodigal", boolSplit=False)
#     pbar.update(1)
# pbar.close()
# exit()


dicoCount = {}
for file in os.listdir("/mnt/c/Users/dgoudenege/Desktop/PICMI_RNAseq/htseq_count"):
    if "_count.tsv" in file:
        lstLines = read_file("/mnt/c/Users/dgoudenege/Desktop/PICMI_RNAseq/htseq_count/"+file,yaspinBool=False)
        dicoCount[file] = {}
        for line in lstLines:
            dicoCount[file][line.split("\t")[0]] = line.split("\t")[1]

print(dicoCount.keys())

for ref in ["V115"]:
    dicoGFF = make_gff_dict("/mnt/c/Users/dgoudenege/Desktop/PICMI_RNAseq/reference/"+ref+".gff")
    dicoCDS = {}
    for cds in dicoGFF[ref]['CDS']:
        lt = cds['attributes']['locus_tag']
        dicoCDS[lt] = {}
        dicoCDS[lt]['contig'] = cds['seqID']
        dicoCDS[lt]['start'] = str(cds['start'])
        dicoCDS[lt]['end'] = str(cds['end'])
        dicoCDS[lt]['strand'] = str(cds['strand'])
        try: dicoCDS[lt]['gene'] = cds['attributes']['gene']
        except: dicoCDS[lt]['gene'] = ""
        try: dicoCDS[lt]['protein_id'] = cds['attributes']['protein_id']
        except: dicoCDS[lt]['protein_id'] = ""
        try: dicoCDS[lt]['pseudo'] = cds['attributes']['pseudo']
        except: dicoCDS[lt]['pseudo'] = "false"
        try: dicoCDS[lt]['product'] = cds['attributes']['product']
        except: dicoCDS[lt]['product'] = ""

    dicoLT = {}
    for gene in dicoGFF[ref]['gene']:
        lt = gene['attributes']['locus_tag']
        if lt in dicoCDS:
            dicoLT[lt] = lt+";CDS;"+dicoCDS[lt]['contig']+";"+dicoCDS[lt]['start']+";"+dicoCDS[lt]['end']+";"+dicoCDS[lt]['strand']+";"+dicoCDS[lt]['protein_id']+";"+dicoCDS[lt]['gene']+";"+dicoCDS[lt]['pseudo']+";"+dicoCDS[lt]['product']
        else:
            dicoLT[lt] = lt+";"+gene['attributes']['gene_biotype']+";"+gene['seqID']+";"+str(gene['start'])+";"+str(gene['end'])+";"+str(gene['strand'])+";;;;"

    header = "locustag;biotype;contig;start;end;strand;protein_id;gene;pseudo;product"
    for phage in ["P115","P27"]:
        for time in ["t0","t5","15","30","60"]:
            for replicat in ["replicat1","replicat2","replicat3"]:
                header += ";"+phage+"_"+time+"_"+replicat
    print(header)
    for lt in dicoLT:
        line = dicoLT[lt]
        for phage in ["P115","P27"]:
            for time in ["t0","t5","t15","t30","t60"]:
                for replicat in ["replicat1","replicat2","replicat3"]:
                    line += ";"+dicoCount[replicat+"_"+phage+"_"+time+"_"+ref+"_count.tsv"][lt]
        print(line)
    print("\n\n")
exit()


pathDIR = "/mnt/c/Users/dgoudenege/Desktop/vibrioDB_satellite_finder"
dicoPICMI = {}
dicoAccJson = load_json("/mnt/g/db/vibrioDB/assembly_summary.json")
lstFiles = os.listdir(pathDIR)
pbar = tqdm(total=len(lstFiles),ncols=75,leave=False,desc="",file=sys.stdout,bar_format="  {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}")
for tsv in lstFiles:
    lstLines = read_file(pathDIR+"/"+tsv,yaspinBool=False)
    accession = tsv.split("_")[0]+"_"+tsv.split("_")[1]
    if len(lstLines) > 0:
        dicoPICMI[accession] = {}
        for i in range(1,len(lstLines),1):
            if lstLines[i] != "":
                splitLine = lstLines[i].split("\t")
                hit_id = splitLine[2]
                systNum = splitLine[6].replace("UserReplicon_PICMI_","")
                if not systNum in dicoPICMI[accession]:
                    dicoPICMI[accession][systNum] = {}
                if splitLine[12] == "Int_Tyr_PF00589.25": dicoPICMI[accession][systNum]["Int"] = hit_id
                elif splitLine[12] == "alpA_PF05930.15": dicoPICMI[accession][systNum]["AlpA"] = hit_id
                elif splitLine[12] == "Vibrio_1": dicoPICMI[accession][systNum]["M1"] = hit_id
                elif splitLine[12] == "Vibrio_2": dicoPICMI[accession][systNum]["M2"] = hit_id
                elif splitLine[12] == "Vibrio_3": dicoPICMI[accession][systNum]["M3"] = hit_id
                elif splitLine[12] == "Fis_PF02954.22": dicoPICMI[accession][systNum]["Fis"] = hit_id
    pbar.update(1)
pbar.close()

print("Accession;OrgName;SystNum;Fis;M1;M2;M3;AlpA;Int")
for accession in dicoPICMI:
    for systNum in dicoPICMI[accession]:
        if accession in dicoAccJson: orgName = dicoAccJson[accession]["org_name"]
        elif accession.replace("GCA_","GCF_") in dicoAccJson: orgName = dicoAccJson[accession.replace("GCA_","GCF_")]["org_name"]
        elif accession.replace("GCF_","GCA_") in dicoAccJson: orgName = dicoAccJson[accession.replace("GCF_","GCA_")]["org_name"]
        else: orgName = "ND"
        line = accession+";"+orgName+";"+systNum
        for geneName in ["Fis","M1","M2","M3","AlpA","Int"]:
            if geneName in dicoPICMI[accession][systNum]: line += ";"+dicoPICMI[accession][systNum][geneName]
            else: line += ";"
        print(line)
exit()


dicoFASTA = make_fasta_dict("/mnt/c/Users/dgoudenege/Desktop/VICR.0123.00004.src.gen")
dicoLT = {}
for key in dicoFASTA:
    tRNA = False
    dicoExt = {}
    splitKey = key.split(" | ")
    lt = splitKey[0].split(" ")[0]
    geneName = splitKey[0].split(" ")[2]
    if geneName == "NA": geneName = ""
    product = splitKey[1]
    if product == "NA": product = ""    
    if splitKey[2] == "NA": ecnumber = ""
    else: ecnumber = splitKey[2]
    if "Aragorn" in splitKey[3]:
        tRNA = True
        uniprotkb = ""
    elif splitKey[3] != "NA":
        extDBname = splitKey[3].split(":")[1]
        extDBvalue = splitKey[3].split(":")[2]
        dicoExt[extDBname] = extDBvalue
    if splitKey[4] == "NA": cog = ""
    else: cog = splitKey[4].split(":")[1]
    dicoLT[lt] = {'geneName':geneName, 'product':product, 'ecnumber':ecnumber ,'dicoExt':dicoExt, 'cog':cog, 'tRNA':tRNA}
dump_json(dicoLT,"/mnt/c/Users/dgoudenege/Desktop/VICR.0123.00004.json")
exit()



for file in os.listdir("/tmp/gemini_replicon_distribution_0jj56l9x"):
    if ".fasta" in file:
        dicoFASTA = make_fasta_dict("/tmp/gemini_replicon_distribution_0jj56l9x/"+file)
        maxSize = 0
        for key in dicoFASTA:
            maxSize = max(maxSize, len(dicoFASTA[key]))
        print(file+";"+str(maxSize))
exit()

# dicoFASTA = make_fasta_dict("/tmp/temp.fasta")
# for key in dicoFASTA:
#     TMP = open("/tmp/in.fasta",'w')
#     TMP.write(">"+key+"\n"+dicoFASTA[key]+"\n")
#     TMP.close()
#     os.system("/mnt/c/Users/dgoudenege/Tools/prodigal_v2.6.3 -a /tmp/in.faa -g 11 -i /tmp/in.fasta -o /tmp/in.gbk -p meta")
#     print("macsyfinder --db-type ordered_replicon --sequence-db /tmp/in.faa --models CONJScan_plasmids system --models-dir /mnt/c/Users/dgoudenege/Tools/macsyfinder/data/models")
#     break
# exit()

# for file in os.listdir("/mnt/c/Users/dgoudenege/Desktop/Scaffolding/Final_qcov75"):
#     dicoAssembly = make_fasta_dict("/mnt/c/Users/dgoudenege/Desktop/Scaffolding/Final_qcov75/"+file)
#     orgName = file.replace("_scaffolds.fasta","").replace("Vibrio_crassostreae_","").replace("Vibrio_sp._GVpop11_","").replace("Vibrio_sp._","")
#     for contig in dicoAssembly:
#         if not "scaffold" in contig:
#             overlap = check_circular(dicoAssembly[contig])
#             if overlap!="" and len(dicoAssembly[contig])>1000:

#                 # print(orgName+";"+contig+";"+str(len(dicoAssembly[contig]))+";"+overlap)
#                 print(">"+orgName+"_____"+contig+"\n"+dicoAssembly[contig])
# exit()


lstToAlign = ["chr2_364"]
lstLines = read_file("/mnt/c/Users/dgoudenege/Desktop/22_O_6_vs_chaga_crass.rbh", yaspinBool=False)
for toAlign in lstToAlign:
    lstTarget = []
    for line in lstLines:
        splitLine = line.split("\t")
        query = splitLine[0]
        if query == toAlign:
            target = splitLine[1]
            ltTarget = target.split("|")[0]
            lstTarget.append(ltTarget)
for ffnFile in os.listdir("/mnt/c/Users/dgoudenege/Data/myVibrio/MAGE/FFN"):
    if ".ffn" in ffnFile and ("crassostreae" in ffnFile or "chagasii" in ffnFile):
        dicoFFN = make_fasta_dict("/mnt/c/Users/dgoudenege/Data/myVibrio/MAGE/FFN/"+ffnFile)
        for key in dicoFFN:
            if key.split("|")[0] in lstTarget:
                print(">"+ffnFile.replace(".ffn.gz","")+"\n"+dicoFFN[key])
                break

exit()


dicoGBK = list(make_gbk_dict("/mnt/c/Users/dgoudenege/Desktop/22_O_6.gbk").values())[0]
dicoStrand = {}
dicoSize = {}
for contig in dicoGBK:
    for lt in dicoGBK[contig]['dicoLT']:
        dicoStrand[lt] = dicoGBK[contig]['dicoLT'][lt]['strand']
        dicoSize[lt] = len(dicoGBK[contig]['dicoLT'][lt]['protSeq'])
dicoRef = {}
lstLT = list(make_fasta_dict("/mnt/c/Users/dgoudenege/Desktop/22_O_6.ffn").keys())
for lt in lstLT:
    dicoRef[lt] = {'setOrg':set(), 'setDescr':set(), 'refLT':"", 'min_id':100}
lstLines = read_file("/mnt/c/Users/dgoudenege/Desktop/22_O_6_vs_chaga_crass.rbh")
for line in lstLines:
    splitLine = line.split("\t")
    query = splitLine[0]
    target = splitLine[1]
    ltTarget = target.split("|")[0]
    descrTarget = target.split("|")[1]
    orgTarget = target.split("|")[2]
    pident = float(splitLine[2])
    dicoRef[query]['setOrg'].add(orgTarget)
    dicoRef[query]['setDescr'].add(descrTarget.replace("#"," "))
    dicoRef[query]['min_id'] = min(dicoRef[query]['min_id'],pident)
    if "22_O_6" in orgTarget:
        dicoRef[query]['refLT'] = ltTarget

print("22_0_6 complete LT;22_0_6 MAGE LT;length (pb);strand;nb Org;min nident;descriptions")
for lt in lstLT:
    print(lt+";"+dicoRef[lt]['refLT']+";"+str(dicoSize[lt]*3+3)+";"+str(dicoStrand[lt])+";"+str(len(dicoRef[lt]['setOrg']))+";"+str(dicoRef[lt]['min_id'])+";"+str(dicoRef[lt]['setDescr']))
exit()



dicoFASTA = make_fasta_dict("/tmp/infA.fasta")
for key in dicoFASTA:
    OUT = open("/tmp/temp.fasta",'w')
    OUT.write(">"+key+"\n"+dicoFASTA[key])
    OUT.close()
    cmd = "/mnt/c/Users/dgoudenege/Tools/prodigal_v2.6.3 -a /tmp/temp.faa -d /tmp/temp.ffn -g 11 -i /tmp/temp.fasta -p single -f gff -o /tmp/temp.gff > /dev/null 2>&1"
    os.system(cmd)
    make_gbk_from_fasta(pathIN1="/tmp/temp.fasta", pathIN2="/tmp/temp.ffn", pathIN3="/tmp/temp.faa", pathOUT="/tmp/"+key+".gbk", identifier=key, topology="linear", division="BCT", taxID=0, pathIN4="None")
exit()

seqInfA = "ATGGCTAAAGAAGACGTAATCGAGATGCAAGGCACTGTCCTTGATACTCTACCAAACACAATGTTCCGTGTTGAGCTTGAAAACGGTCACGTAGTGACAGCACACATCTCTGGTAAAATGCGTAAGAACTACATCCGTATTCTTACTGGTGATAAAGTAACTGTTGAGATGACTCCATACGACCTTTCTAAAGGCCGCATCGTCTTCCGTGCTCGTTAA"
# lstFiles = os.listdir("/mnt/c/Users/dgoudenege/Data/myVibrio/Chagasii")
lstFiles = os.listdir("/mnt/c/Users/dgoudenege/Data/myVibrio/Crassostreae")
pbar = tqdm(total=len(lstFiles),ncols=75,leave=False,desc="",file=sys.stdout,bar_format="  {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}")
dicoSeq = {}
for file in lstFiles:
    if ".fna" in file:
        # dicoFNA = make_fasta_dict("/mnt/c/Users/dgoudenege/Data/myVibrio/Chagasii/"+file)
        dicoFNA = make_fasta_dict("/mnt/c/Users/dgoudenege/Data/myVibrio/Crassostreae/"+file)
        for key in dicoFNA:
            if seqInfA in dicoFNA[key]:
                splitregion = dicoFNA[key].split(seqInfA)
                dicoSeq[file.replace(".fna.gz","")] = splitregion[0][-20000:]+seqInfA+splitregion[1][:20000]
            elif seqInfA in reverse_complement(dicoFNA[key]):
                splitregion = reverse_complement(dicoFNA[key]).split(seqInfA)
                dicoSeq[file.replace(".fna.gz","")] = splitregion[0][-20000:]+seqInfA+splitregion[1][:20000]
    pbar.update(1)
pbar.close()
for key in dicoSeq:
    print(">"+key+"\n"+dicoSeq[key])
exit()

# lstGene = set(list(make_fasta_dict("/tmp/chaga_crass/22_O_6.ffn").values()))
# print(len(lstGene))
# for file in os.listdir("/tmp/chaga_crass"):
#     if ".ffn" in file:
#         lstGene2 = set(list(make_fasta_dict("/tmp/chaga_crass/"+file).values()))
#         lstGene = lstGene & lstGene2
#         if len(lstGene) < 10 :
#             print(len(lstGene),lstGene,file)
#         else:
#             print(len(lstGene),file)
exit()


dicoCluster = make_rbhcluster_dict(pathIN="/mnt/c/Users/dgoudenege/Desktop/chaga1_familyMixed/rbh", pathIN2="/mnt/c/Users/dgoudenege/Desktop/chaga1_familyMixed", pathJSON="None", idThrClust=25, covThrClust=30, ext=".rbh", ext2=".faa")
# dicoCluster = load_json("/mnt/c/Users/dgoudenege/Works/Phage_analysis/myocluster_subfamily_fig3d/rbh.json")
cpt = 0
for clusterNum in dicoCluster:
    setOrg = set()
    p115 = ""
    for gene in dicoCluster[clusterNum]:
        orgName = gene.split("[")[1].replace("]", "")
        if orgName=="Vibrio_phage_115E34-1": p115 = gene
        setOrg.add(orgName)
    if len(setOrg) == 12:  # to avoid paralogous or split gene
        # print(p115)
        cpt+=1
print(cpt)
exit()

dicoHost = {"Vibrio_phage_511E55-1":"chaga","Vibrio_phage_466E53-1":"chaga","Vibrio_phage_115E34-1":"chaga","Vibrio_phage_120E34-1":"chaga","Vibrio_phage_207E29-1":"crass","Vibrio_phage_455E52-1":"chaga","Vibrio_phage_98E28-6a":"crass","Vibrio_phage_18E29-1":"crass","Vibrio_phage_91E28-1a":"crass","Vibrio_phage_12E28-1":"crass","Vibrio_phage_217E38-1":"chaga","Vibrio_phage_191E37-1":"chaga"}
# RENDERING tree
for file in os.listdir("/mnt/c/Users/dgoudenege/Desktop/chaga1_familyMixed/core_tree"):
    if ".treefile" not in file: continue
    coreNum = file.replace("_iqtree2.treefile","")
    tree = Tree("/mnt/c/Users/dgoudenege/Desktop/chaga1_familyMixed/core_tree/"+file)
    for node in tree.traverse("postorder"):
        # Do some analysis on node
        node.name = node.name.split("|")[0]
    ts = TreeStyle()
    ts.show_branch_support = True
    for n in tree.traverse():
        nstyle = NodeStyle()
        nstyle["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
        nstyle["hz_line_type"] = 0
        if n.name == "Vibrio_phage_115E34-1":
            tree.set_outgroup(n)
        try: host = dicoHost[n.name]
        except KeyError: n.set_style(nstyle) ; continue
        if host == "crass":
            nstyle["fgcolor"] = "#ff5599"
        elif host == "chaga":
            nstyle["fgcolor"] = "#37c871"
        # nstyle["shape"] = "sphere"
        nstyle["size"] = 10
        n.set_style(nstyle)
    ts.show_leaf_name = True
    ts.show_branch_support = True
    ts.title.add_face(TextFace(coreNum, fsize=12), column=0)
    # tree.show(tree_style=ts)
    # tree.render("/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/individual_tree/"+dicoLT[clusterNum]+".svg", tree_style=ts)
    tree.render("/mnt/c/Users/dgoudenege/Desktop/chaga1_familyMixed/core_tree/"+coreNum+".png", units="mm", tree_style=ts)

exit()


dicoTMP = {"Vibrio_chagasii_43_P_273":{"abi":"VCHA43P273_v1_250014","ass":"VCHA43P273_v1_250013"},"Vibrio_chagasii_43_P_282":{"abi":"VCHA43P282_v1_250050","ass":"VCHA43P282_v1_250051"},"Vibrio_crassostreae_31_O_70":{"abi":"REV070_v1_370010","ass":"REV070_v1_370009"},"Vibrio_crassostreae_42_O_251":{"abi":"REV251_v1_320056","ass":"REV251_v1_320057"},"Vibrio_crassostreae_45_O_294":{"abi":"REV294_v1_340010","ass":"REV294_v1_340009"},"Vibrio_crassostreae_46_O_330":{"abi":"REV330_v1_340055","ass":"REV330_v1_340056"},"Vibrio_crassostreae_46_P_336":{"abi":"REV336_v1_460011","ass":"REV336_v1_460010"},"Vibrio_crassostreae_7D8_10":{"abi":"GV2335_v1_320010","ass":"GV2335_v1_320009"},"Vibrio_crassostreae_7F1_18":{"abi":"GV2333_v1_320057","ass":"GV2333_v1_320058"},"Vibrio_crassostreae_GV1681":{"abi":"GV2153_v1_250002","ass":"GV2153_v1_250001"},"Vibrio_sp._DDZZcladef-J2-29":{"abi":"VCR29J2v1_80057","ass":"VCR29J2v1_80058"},"Vibrio_sp._DDZZcladef-J2-31":{"abi":"VCR31J2v1_1300137","ass":"VCR31J2v1_1300136"},"Vibrio_sp._GVpop11-8T8-11":{"abi":"GV8T811_v1_510091","ass":"GV8T811_v1_510092"},"Vibrio_sp._MITpop18-1F_55":{"abi":"VB1855_v1_920006","ass":"VB1855_v1_920007"}}
for file in os.listdir("/tmp/GBK"):
    pathGBK = "/tmp/GBK/"+file
    dicoGBK = list(make_gbk_dict(pathGBK).values())[0]
    org = file.replace(".gbk","")
    for contig in dicoGBK:
        start = 0
        end = 0
        seqabin = ""
        seqabip = ""
        seqassn = ""
        seqassp = ""
        strand = 0
        for lt in dicoGBK[contig]['dicoLT']:
            if lt == dicoTMP[org]["abi"]:
                strand = dicoGBK[contig]['dicoLT'][lt]['strand']
                if dicoGBK[contig]['dicoLT'][lt]['strand'] == 1: start = dicoGBK[contig]['dicoLT'][lt]['start']
                else: end = dicoGBK[contig]['dicoLT'][lt]['end']
                seqabin = dicoGBK[contig]['dicoLT'][lt]['geneSeq']
                seqabip = dicoGBK[contig]['dicoLT'][lt]['protSeq']
            elif lt == dicoTMP[org]["ass"]:
                if dicoGBK[contig]['dicoLT'][lt]['strand'] == 1: end = dicoGBK[contig]['dicoLT'][lt]['end']
                else: start = dicoGBK[contig]['dicoLT'][lt]['start']
                seqassn = dicoGBK[contig]['dicoLT'][lt]['geneSeq']
                seqassp = dicoGBK[contig]['dicoLT'][lt]['protSeq']                
        if start!=0 and end!=0:
            pathFNA = "/tmp/PLOT/"+org+".fna"
            FNA = open(pathFNA,"w")
            pathFAA = "/tmp/PLOT/"+org+".faa"
            FAA = open(pathFAA,"w")            
            pathFFN = "/tmp/PLOT/"+org+".ffn"
            FFN = open(pathFFN,"w")            
            if strand == 1: FNA.write(">"+org+"\n"+dicoGBK[contig]['seq'][start:end]+"\n")
            else: FNA.write(">"+org+"\n"+reverse_complement(dicoGBK[contig]['seq'][start:end])+"\n")
            FAA.write(">abi\n"+seqabip+"\n")
            FFN.write(">abi\n"+seqabin+"\n")
            FAA.write(">ass\n"+seqassp+"\n")
            FFN.write(">ass\n"+seqassn+"\n")
            FNA.close()
            FAA.close()
            FFN.close()
            make_gbk_from_fasta(pathIN1=pathFNA, pathIN2=pathFFN, pathIN3=pathFAA, pathOUT="/tmp/PLOT/"+org+".gbk", identifier=org, topology="linear", division="BCT", taxID=0, pathIN4="None")

exit()



for org in os.listdir("/mnt/c/Users/dgoudenege/Desktop/PACBIO/Ordered"):
    pathFASTA = "/mnt/c/Users/dgoudenege/Desktop/PACBIO/Ordered/"+org
    dicoFASTA = make_fasta_dict(pathFASTA)
    orgName = org.replace(".fasta","")
    for key in dicoFASTA:
        TMP = open("/mnt/c/Users/dgoudenege/Desktop/PACBIO/Ordered_replicon/"+orgName+"_"+key.split(" ")[0]+".fasta",'w')
        TMP.write(">"+key+"\n"+dicoFASTA[key]+"\n")
        TMP.close()
exit()


for clade in ["Clade1","Clade2","Clade3","Clade4","Clade5","Clade6","Clade7","Clade8","Hors_Clade"]:
    pathDirClade = "/mnt/c/Users/dgoudenege/Desktop/PACBIO/"+clade

    for phage in os.listdir(pathDirClade):
        pathFASTA = pathDirClade+"/"+phage
        dicoFASTA = make_fasta_dict(pathFASTA)
        NEW = open("/mnt/c/Users/dgoudenege/Desktop/PACBIO/Ordered/all_plasmid."+phage,'w')
        for key in dicoFASTA:
            NEW.write(">"+key+"\n"+dicoFASTA[key]+"\n")
            NEW.write(">rc_"+key+"\n"+reverse_complement(dicoFASTA[key])+"\n")
        NEW.close()
exit()



for phage in os.listdir("/mnt/c/Users/dgoudenege/Desktop/B7117/filtered_assembly"):
    if ".fasta" in phage:
        pathFASTA = "/mnt/c/Users/dgoudenege/Desktop/B7117/filtered_assembly/"+phage
        dicoGenome = make_fasta_dict(pathFASTA)
        circularRepeat = list(dicoGenome.keys())[0].split("circular: ")[1]
        phageGenome = list(dicoGenome.values())[0]
        pathFFN = "/mnt/c/Users/dgoudenege/Desktop/B7117/filtered_assembly/"+phage.replace(".fasta",".ffn")
        dicoFFN = make_fasta_dict(pathFFN)
        orgName = phage.replace(".fasta","")
        pathOUT = "/mnt/c/Users/dgoudenege/Desktop/B7117/filtered_assembly/"+phage.replace(".fasta",".fna")
        # print(orgName)
        cmdDiamond = "rm -f /tmp/temp.out ; diamond blastx --query "+pathFFN+" --db /mnt/g/db/DMND/largeTerminase.dmnd --max-target-seqs 1 --id 30 --subject-cover 50 --outfmt 6 qseqid | sort -u > /tmp/temp.out"
        os.system(cmdDiamond)
        lstLines = read_file("/tmp/temp.out")
        if len(lstLines) == 1:
            termLT = lstLines[0]
            termSeq = dicoFFN[termLT]
        else:
            print(orgName,"terminase not Found")
            print(cmdDiamond)
            break
        if termSeq in phageGenome:
            splitPhageGenome = phageGenome.split(termSeq)
            newPhageGenome = termSeq+splitPhageGenome[1].split(circularRepeat)[0]+splitPhageGenome[0]
            OUT = open(pathOUT,'w')
            OUT.write(">"+orgName+"\n"+newPhageGenome+"\n")
            OUT.close()
        else:
            phageGenome = reverse_complement(phageGenome)
            if termSeq in phageGenome:
                splitPhageGenome = phageGenome.split(termSeq)
                newPhageGenome = termSeq+splitPhageGenome[1].split(reverse_complement(circularRepeat))[0]+splitPhageGenome[0]
                OUT = open(pathOUT,'w')
                OUT.write(">"+orgName+"\n"+newPhageGenome+"\n")
                OUT.close()
            else:
                print(orgName,"  terminase seq not Found")
                break

exit()


lstAss = os.listdir("/mnt/g/db/phages")
for ass in lstAss:
    pathASS = "/mnt/g/db/phages/"+ass
    FFN = ""
    phanotate = ""
    FAA = ""
    if os.path.isdir(pathASS):
        for file in os.listdir(pathASS):
            if "_gene.ffn.gz" in file: FFN = pathASS+"/"+file
            if "_phanotate.ffn.gz" in file: phanotate = pathASS+"/"+file
            if "_protein.faa.gz" in file: FAA = pathASS+"/"+file
        if FFN == "": print("no FFN for :"+ass)
        if phanotate == "": print("no phanotate for :"+ass)
        if FAA == "": 
            print(FFN)
            # print("no FAA for :"+ass)

exit()

pathDIR = "/mnt/c/Users/dgoudenege/Desktop/Scaffolding/Clade1_qcov75/from_scaffolding"
pathCHR1 = pathDIR+"/all_chr1.fasta"
CHR1 = open(pathCHR1, 'w')
pathCHR2 = pathDIR+"/all_chr2.fasta"
CHR2 = open(pathCHR2, 'w')
for file in os.listdir(pathDIR):
    dicoFASTA = make_fasta_dict(pathDIR+"/"+file)
    cpt_ctg = 1
    orgName = file.replace("_scaffolds.fasta","").replace("_complete.fasta","")
    pathReformatFASTA = pathDIR+"/"+orgName+".fasta"
    pathReformatFASTACHR1 = pathDIR+"/"+orgName+"_chr1.fasta"
    pathReformatFASTACHR2 = pathDIR+"/"+orgName+"_chr2.fasta"
    REFORMAT = open(pathReformatFASTA, 'w')
    REFORMAT1 = open(pathReformatFASTACHR1, 'w')
    REFORMAT2 = open(pathReformatFASTACHR2, 'w')
    for key in dicoFASTA:
        replicon = key.replace(" ","_").split("_")[0]
        if "chr1" in replicon:
            if "scaffold" in file: formatKey = "chr1|scaffold ["+orgName+"]"
            else: formatKey = "chr1|complete ["+orgName+"]"
            REFORMAT1.write(">"+formatKey+"\n"+dicoFASTA[key]+"\n")
            CHR1.write(">"+formatKey+"\n"+dicoFASTA[key]+"\n")
        elif "chr2" in replicon:
            if "scaffold" in file: formatKey = "chr2|scaffold ["+orgName+"]"
            else: formatKey = "chr2|complete ["+orgName+"]"
            REFORMAT2.write(">"+formatKey+"\n"+dicoFASTA[key]+"\n")
            CHR2.write(">"+formatKey+"\n"+dicoFASTA[key]+"\n")
        else:
            formatKey = "contig"+str(cpt_ctg).zfill(2)+" ["+orgName+"]"
            cpt_ctg += 1
        REFORMAT.write(">"+formatKey+"\n"+dicoFASTA[key]+"\n")
    REFORMAT.close()     
    REFORMAT1.close()     
    REFORMAT2.close()
CHR1.close()
CHR2.close() 
exit()


lstFiles = []
for folder in lst_genome:
    if os.path.isfile("/mnt/g/db/phages/"+folder+"/"+folder+"_genomic.fna.gz"):
        lstFiles.append("/mnt/g/db/phages/"+folder+"/"+folder+"_genomic.fna.gz")

for pathInitialFASTA in lstFiles:
    orgName = os.path.basename(pathInitialFASTA).replace("_genomic.fna.gz","")
    pathReformatFASTA = "/tmp/"+orgName+".fasta"

    # REFORMAT = open(pathReformatFASTA, 'w')
    # dicoFNA = make_fasta_dict(pathInitialFASTA)
    # REFORMAT.write(">"+orgName+"\n"+("N"*100).join(list(dicoFNA.values()))+"\n")
    # REFORMAT.close()

    seqs = []
    with gzip.open(pathInitialFASTA,'rt') as fasta:        
        prev_seq = []
        for line in fasta:
            if line.startswith(">"):
                seqs.append("".join(prev_seq))
                prev_seq = []
            else:
                prev_seq.append(line.rstrip())
    seqs.append("".join(prev_seq))
    REFORMAT = open(pathReformatFASTA, 'w')
    REFORMAT.write(">"+orgName+"\n"+("N"*100).join(seqs)[100:]+"\n")
    REFORMAT.close()
exit()


# lstFiles = ["GCA_003146765.1_ASM314676v1.json","GCA_003925005.1_ASM392500v1.json","GCA_003925205.1_ASM392520v1.json","GCA_003925405.1_ASM392540v1.json","GCA_003925605.1_ASM392560v1.json"]:
lstFiles = os.listdir("/mnt/g/VIRIDIC2/viridic_json")
pbar = tqdm(total=len(lstFiles),ncols=75,leave=False,desc="",file=sys.stdout,bar_format="  {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}")
for file in lstFiles:
    pathJSON = "/mnt/g/VIRIDIC2/viridic_json/"+file
    dicoViridic = load_json(pathJSON)
    sqldicoViridic = SqliteDict("/tmp/test_sqlitedict/"+file.replace(".json",".sqlite"), outer_stack=False)
    for key in dicoViridic:
        sqldicoViridic[key] = dicoViridic[key]
    sqldicoViridic.commit()
    sqldicoViridic.close()
    pbar.update()
pbar.close()

bigDico = {}
pbar = tqdm(total=len(lstFiles),ncols=75,leave=False,desc="",file=sys.stdout,bar_format="  {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}")
for file in lstFiles:
    sqldicoViridic = SqliteDict("/tmp/test_sqlitedict/"+file.replace(".json",".sqlite"), outer_stack=False)
    bigDico[file] = sqldicoViridic
    pbar.update()
pbar.close()

# for file in bigDico:
#     sqliteDicoViridic = bigDico[file]
#     for key1 in sqliteDicoViridic.keys():
#         print(file,key1,sqliteDicoViridic[key1])
exit()

dicoLinkName = {}
for viridic in os.listdir("/mnt/g/VIRIDIC/viridic_json"): # /mnt/g/VIRIDIC/viridic_json/xxxxxxx_GCA_004777725.1
    if "_GCA" in viridic:
        gcaAccession = "GCA_"+viridic.split("_GCA_")[1].replace(".json","")
        dicoLinkName[gcaAccession] = { 'viridic': viridic.replace(".json","") , 'genomeDB':"" }

for phage in os.listdir("/mnt/g/db/phages"): # /mnt/g/db/phages/GCA_004777725.1_ASM477772v1/
    pathDIRphage = "/mnt/g/db/phages/"+phage
    if os.path.isdir(pathDIRphage):
        gcaAccession = "GCA_"+phage.split("_")[1]   
        try : dicoLinkName[gcaAccession]['genomeDB'] = phage
        except KeyError: dicoLinkName[gcaAccession] = { 'viridic': "" , 'genomeDB':phage }

# Using a list comprehension to make a list of the keys to be deleted
delete = [key for key in dicoLinkName if dicoLinkName[key]['viridic']=="" or dicoLinkName[key]['genomeDB']==""]
# delete the key/s
for key in delete: del dicoLinkName[key]

dicoViridicToDB = {}
for gca in dicoLinkName:
    dicoViridicToDB[dicoLinkName[gca]['viridic']] = dicoLinkName[gca]['genomeDB']


missing = set()
pbar = tqdm(total=len(dicoLinkName),ncols=75,leave=False,desc="",file=sys.stdout,bar_format="  {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}")
for gca in dicoLinkName:
    # BLASTout
    pathBLASTOUT = "/mnt/g/VIRIDIC/blastn_out/"+dicoLinkName[gca]['viridic']+".out"
    pathNEWBLASTOUT = "/mnt/g/VIRIDIC2/blastn_out/"+dicoLinkName[gca]['genomeDB']+".out"
    NEWBLASTOUT = open(pathNEWBLASTOUT, 'w')
    lstLines = read_file(pathBLASTOUT, yaspinBool = False)
    for line in lstLines:
        splitLine = line.split("\t")
        phage1 = splitLine[0]
        phage2 = splitLine[1]
        following = "\t".join(splitLine[2:])
        try:
            newphage1 = dicoViridicToDB[phage1]
            newphage2 = dicoViridicToDB[phage2]
        except KeyError: continue
        NEWBLASTOUT.write(newphage1+"\t"+newphage2+"\t"+following+"\n")
    NEWBLASTOUT.close()    
    # BLAST json
    pathBLASTJSON = "/mnt/g/VIRIDIC/blastn_out/"+dicoLinkName[gca]['viridic']+".json"
    pathNEWBLASTJSON = "/mnt/g/VIRIDIC2/blastn_out/"+dicoLinkName[gca]['genomeDB']+".json"
    lstBLAST = load_json(pathBLASTJSON)
    lstNEWBLAST = []
    for phage in lstBLAST:
        try: lstNEWBLAST.append(dicoViridicToDB[phage])
        except KeyError: continue
    dump_json(lstNEWBLAST, pathNEWBLASTJSON)
    # SIM json
    pathSIMJSON = "/mnt/g/VIRIDIC/viridic_sim/"+dicoLinkName[gca]['viridic']+".json"
    pathNEWSIMJSON = "/mnt/g/VIRIDIC2/viridic_sim/"+dicoLinkName[gca]['genomeDB']+".json"
    dicoSIM = load_json(pathSIMJSON)
    dicoNEWSIM = {}
    for phage in dicoSIM:
        try: dicoNEWSIM[dicoViridicToDB[phage]] = dicoSIM[phage]
        except KeyError: continue
    dump_json(dicoNEWSIM, pathNEWSIMJSON)
    # VIRIDIC json
    pathVIRIDICJSON = "/mnt/g/VIRIDIC/viridic_json/"+dicoLinkName[gca]['viridic']+".json"
    pathNEWVIRIDICJSON = "/mnt/g/VIRIDIC2/viridic_json/"+dicoLinkName[gca]['genomeDB']+".json"
    dicoVIRIDIC = load_json(pathVIRIDICJSON)
    dicoNEWVIRIDIC = {}
    for phage in dicoVIRIDIC:
        try: dicoNEWVIRIDIC[dicoViridicToDB[phage]] = dicoVIRIDIC[phage]
        except KeyError: continue
    dump_json(dicoNEWVIRIDIC, pathNEWVIRIDICJSON)
    pbar.update(1)
pbar.close()

exit()



lst = ["REV330_v1_10026","REV330_v1_290004","","REV070_v1_190037","REV070_v1_380002","REV070_v1_380051","","REV014_v1_10260","REV014_v1_70187","","REV015_v1_10575","REV015_v1_60069","","REV183_v1_70042","","REV099_v1_150133","REV099_v1_150137","REV099_v1_150140","REV099_v1_150141","REV099_v1_150151","REV099_v1_220034","","REV024_v1_80151","REV024_v1_80160","REV024_v1_80174","REV024_v1_90001","REV045_v1_100035","REV045_v1_110083","REV045_v1_30279","","REV051_v1_250078","REV051_v1_250079","","GV1646_v1_50165","","REV207_v1_60032","","REV267_v1_10888","REV267_v1_60172","REV267_v1_70138","","REV207_v1_60032","","REV207_v1_11110","REV207_v1_60061","REV207_v1_60066","","REV028_v1_140168","REV028_v1_40035","REV028_v1_40057","REV028_v1_90054"]
dico = {"1":{"strain":"3870", "lt":"GV1646_v1_50165"}, "2":{"strain":"3872", "lt":"GV1646_v1_50165"}, "3":{"strain":"3873", "lt":"GV1646_v1_50165"}, "4":{"strain":"3874", "lt":"GV1646_v1_50165"}, "5":{"strain":"3875", "lt":"GV1646_v1_50165"}, "6":{"strain":"3876", "lt":"GV1646_v1_50165"}, "7":{"strain":"3882", "lt":"REV014_v1_10260"}, "8":{"strain":"3883", "lt":"REV014_v1_10260"}, "9":{"strain":"3884", "lt":"REV014_v1_10260"}, "10":{"strain":"3885", "lt":"REV014_v1_10260"}, "11":{"strain":"3886", "lt":"REV014_v1_10260"}, "12":{"strain":"3888", "lt":"REV014_v1_10260"}, "13":{"strain":"3880", "lt":"REV014_v1_70187"}, "14":{"strain":"3881", "lt":"REV014_v1_70187"}, "15":{"strain":"3882", "lt":"REV014_v1_70187"}, "16":{"strain":"3883", "lt":"REV014_v1_70187"}, "17":{"strain":"3884", "lt":"REV014_v1_70187"}, "18":{"strain":"3886", "lt":"REV014_v1_70187"}, "19":{"strain":"3887", "lt":"REV014_v1_70187"}, "20":{"strain":"3830", "lt":"REV015_v1_10575"}, "21":{"strain":"3831", "lt":"REV015_v1_10575"}, "22":{"strain":"3832", "lt":"REV015_v1_10575"}, "23":{"strain":"3833", "lt":"REV015_v1_10575"}, "24":{"strain":"3835", "lt":"REV015_v1_10575"}, "25":{"strain":"3838", "lt":"REV015_v1_10575"}, "26":{"strain":"3833", "lt":"REV015_v1_60069"}, "27":{"strain":"3836", "lt":"REV015_v1_60069"}, "28":{"strain":"3809", "lt":"REV024_v1_10258"}, "29":{"strain":"3810", "lt":"REV024_v1_10258"}, "30":{"strain":"3812", "lt":"REV024_v1_10258"}, "31":{"strain":"3810", "lt":"REV024_v1_80079"}, "32":{"strain":"3811", "lt":"REV024_v1_80079"}, "33":{"strain":"3817", "lt":"REV024_v1_80079"}, "34":{"strain":"3818", "lt":"REV024_v1_80079"}, "35":{"strain":"3813", "lt":"REV024_v1_80151"}, "36":{"strain":"3814", "lt":"REV024_v1_80160"}, "37":{"strain":"3809", "lt":"REV024_v1_80174"}, "38":{"strain":"3812", "lt":"REV024_v1_80174"}, "39":{"strain":"3810", "lt":"REV024_v1_90001"}, "40":{"strain":"3811", "lt":"REV024_v1_90001"}, "41":{"strain":"3816", "lt":"REV024_v1_90001"}, "42":{"strain":"3817", "lt":"REV024_v1_90001"}, "43":{"strain":"3818", "lt":"REV024_v1_90001"}, "44":{"strain":"3906", "lt":"REV028_v1_140168"}, "45":{"strain":"3899", "lt":"REV028_v1_40001"}, "46":{"strain":"3904", "lt":"REV028_v1_40001"}, "47":{"strain":"3902", "lt":"REV028_v1_40035"}, "48":{"strain":"3905", "lt":"REV028_v1_40035"}, "49":{"strain":"3906", "lt":"REV028_v1_40057"}, "50":{"strain":"3907", "lt":"REV028_v1_40057"}, "51":{"strain":"3903", "lt":"REV028_v1_90054"}, "52":{"strain":"3904", "lt":"REV028_v1_90054"}, "53":{"strain":"3860", "lt":"REV045_v1_100035"}, "54":{"strain":"3862", "lt":"REV045_v1_100035"}, "55":{"strain":"3868", "lt":"REV045_v1_100088"}, "56":{"strain":"3868", "lt":"REV045_v1_100141"}, "57":{"strain":"3860", "lt":"REV045_v1_110083"}, "58":{"strain":"3862", "lt":"REV045_v1_110083"}, "59":{"strain":"3866", "lt":"REV045_v1_110083"}, "60":{"strain":"3868", "lt":"REV045_v1_130035"}, "61":{"strain":"3868", "lt":"REV045_v1_190122"}, "62":{"strain":"3868", "lt":"REV045_v1_20109"}, "63":{"strain":"3868", "lt":"REV045_v1_20115"}, "64":{"strain":"3868", "lt":"REV045_v1_20127"}, "65":{"strain":"3868", "lt":"REV045_v1_20159"}, "66":{"strain":"3868", "lt":"REV045_v1_20296"}, "67":{"strain":"3868", "lt":"REV045_v1_210106"}, "68":{"strain":"3868", "lt":"REV045_v1_280005"}, "69":{"strain":"3859", "lt":"REV045_v1_290078"}, "70":{"strain":"3860", "lt":"REV045_v1_290078"}, "71":{"strain":"3861", "lt":"REV045_v1_290078"}, "72":{"strain":"3861", "lt":"REV045_v1_290078"}, "73":{"strain":"3861", "lt":"REV045_v1_290078"}, "74":{"strain":"3862", "lt":"REV045_v1_290078"}, "75":{"strain":"3862", "lt":"REV045_v1_290078"}, "76":{"strain":"3862", "lt":"REV045_v1_290078"}, "77":{"strain":"3862", "lt":"REV045_v1_290078"}, "78":{"strain":"3863", "lt":"REV045_v1_290078"}, "79":{"strain":"3863", "lt":"REV045_v1_290078"}, "80":{"strain":"3864", "lt":"REV045_v1_290078"}, "81":{"strain":"3864", "lt":"REV045_v1_290078"}, "82":{"strain":"3864", "lt":"REV045_v1_290078"}, "83":{"strain":"3864", "lt":"REV045_v1_290078"}, "84":{"strain":"3865", "lt":"REV045_v1_290078"}, "85":{"strain":"3865", "lt":"REV045_v1_290078"}, "86":{"strain":"3865", "lt":"REV045_v1_290078"}, "87":{"strain":"3866", "lt":"REV045_v1_290078"}, "88":{"strain":"3866", "lt":"REV045_v1_290078"}, "89":{"strain":"3866", "lt":"REV045_v1_290078"}, "90":{"strain":"3866", "lt":"REV045_v1_290078"}, "91":{"strain":"3867", "lt":"REV045_v1_290078"}, "92":{"strain":"3868", "lt":"REV045_v1_290078"}, "93":{"strain":"3868", "lt":"REV045_v1_30264"}, "94":{"strain":"3865", "lt":"REV045_v1_30279"}, "95":{"strain":"3868", "lt":"REV045_v1_30355"}, "96":{"strain":"3868", "lt":"REV045_v1_320019"}, "97":{"strain":"3868", "lt":"REV045_v1_340008"}, "98":{"strain":"3868", "lt":"REV045_v1_400038"}, "99":{"strain":"3868", "lt":"REV045_v1_400038"}, "100":{"strain":"3868", "lt":"REV045_v1_400038"}, "101":{"strain":"3868", "lt":"REV045_v1_400038"}, "102":{"strain":"3868", "lt":"REV045_v1_400038"}, "103":{"strain":"3868", "lt":"REV045_v1_400039"}, "104":{"strain":"3868", "lt":"REV045_v1_400039"}, "105":{"strain":"3868", "lt":"REV045_v1_400039"}, "106":{"strain":"3868", "lt":"REV045_v1_400039"}, "107":{"strain":"3868", "lt":"REV045_v1_400039"}, "108":{"strain":"3868", "lt":"REV045_v1_400039"}, "109":{"strain":"3868", "lt":"REV045_v1_400039"}, "110":{"strain":"3868", "lt":"REV045_v1_400039"}, "111":{"strain":"3868", "lt":"REV045_v1_400039"}, "112":{"strain":"3868", "lt":"REV045_v1_40016"}, "113":{"strain":"3868", "lt":"REV045_v1_40143"}, "114":{"strain":"3868", "lt":"REV045_v1_50053"}, "115":{"strain":"3868", "lt":"REV045_v1_530011"}, "116":{"strain":"3868", "lt":"REV045_v1_80080"}, "117":{"strain":"3868", "lt":"REV045_v1_90088"}, "118":{"strain":"3849", "lt":"REV051_v1_250078"}, "119":{"strain":"3850", "lt":"REV051_v1_250078"}, "120":{"strain":"3851", "lt":"REV051_v1_250078"}, "121":{"strain":"3852", "lt":"REV051_v1_250078"}, "122":{"strain":"3853", "lt":"REV051_v1_250078"}, "123":{"strain":"3854", "lt":"REV051_v1_250078"}, "124":{"strain":"3856", "lt":"REV051_v1_250078"}, "125":{"strain":"3857", "lt":"REV051_v1_250078"}, "126":{"strain":"3858", "lt":"REV051_v1_250078"}, "127":{"strain":"3849", "lt":"REV051_v1_250079"}, "128":{"strain":"3850", "lt":"REV051_v1_250079"}, "129":{"strain":"3851", "lt":"REV051_v1_250079"}, "130":{"strain":"3852", "lt":"REV051_v1_250079"}, "131":{"strain":"3853", "lt":"REV051_v1_250079"}, "132":{"strain":"3854", "lt":"REV051_v1_250079"}, "133":{"strain":"3856", "lt":"REV051_v1_250079"}, "134":{"strain":"3857", "lt":"REV051_v1_250079"}, "135":{"strain":"3858", "lt":"REV051_v1_250079"}, "136":{"strain":"3954", "lt":"REV070_v1_190037"}, "137":{"strain":"3957", "lt":"REV070_v1_190037"}, "138":{"strain":"3953", "lt":"REV070_v1_380002"}, "139":{"strain":"3954", "lt":"REV070_v1_380002"}, "140":{"strain":"3957", "lt":"REV070_v1_380002"}, "141":{"strain":"3960", "lt":"REV070_v1_380002"}, "142":{"strain":"3953", "lt":"REV070_v1_380051"}, "143":{"strain":"3954", "lt":"REV070_v1_380051"}, "144":{"strain":"3960", "lt":"REV070_v1_380051"}, "145":{"strain":"3957", "lt":"REV070_v1_60079"}, "146":{"strain":"3846", "lt":"REV099_v1_150133"}, "147":{"strain":"3847", "lt":"REV099_v1_150133"}, "148":{"strain":"3845", "lt":"REV099_v1_150137"}, "149":{"strain":"3839", "lt":"REV099_v1_150140"}, "150":{"strain":"3840", "lt":"REV099_v1_150140"}, "151":{"strain":"3841", "lt":"REV099_v1_150140"}, "152":{"strain":"3843", "lt":"REV099_v1_150140"}, "153":{"strain":"3844", "lt":"REV099_v1_150140"}, "154":{"strain":"3848", "lt":"REV099_v1_150141"}, "155":{"strain":"3846", "lt":"REV099_v1_150151"}, "156":{"strain":"3839", "lt":"REV099_v1_220034"}, "157":{"strain":"3840", "lt":"REV099_v1_220034"}, "158":{"strain":"3841", "lt":"REV099_v1_220034"}, "159":{"strain":"3844", "lt":"REV099_v1_220034"}, "160":{"strain":"3845", "lt":"REV099_v1_220034"}, "161":{"strain":"3846", "lt":"REV099_v1_220034"}, "162":{"strain":"3843", "lt":"REV099_v1_430003"}, "163":{"strain":"3783", "lt":"REV183_v1_130001"}, "164":{"strain":"3782", "lt":"REV183_v1_30160"}, "165":{"strain":"3781", "lt":"REV183_v1_70042"}, "166":{"strain":"3782", "lt":"REV183_v1_70042"}, "167":{"strain":"3784", "lt":"REV183_v1_70042"}, "168":{"strain":"3785", "lt":"REV183_v1_70042"}, "169":{"strain":"3786", "lt":"REV183_v1_70042"}, "170":{"strain":"3788", "lt":"REV183_v1_70042"}, "171":{"strain":"3781", "lt":"REV183_v1_80005"}, "172":{"strain":"3798", "lt":"REV207_v1_11110"}, "173":{"strain":"3789", "lt":"REV207_v1_20923"}, "174":{"strain":"3790", "lt":"REV207_v1_20923"}, "175":{"strain":"3794", "lt":"REV207_v1_20923"}, "176":{"strain":"3795", "lt":"REV207_v1_20923"}, "177":{"strain":"3797", "lt":"REV207_v1_20923"}, "178":{"strain":"3798", "lt":"REV207_v1_20923"}, "179":{"strain":"3889", "lt":"REV207_v1_20923"}, "180":{"strain":"3890", "lt":"REV207_v1_20923"}, "181":{"strain":"3892", "lt":"REV207_v1_20923"}, "182":{"strain":"3893", "lt":"REV207_v1_20923"}, "183":{"strain":"3894", "lt":"REV207_v1_20923"}, "184":{"strain":"3895", "lt":"REV207_v1_20923"}, "185":{"strain":"3896", "lt":"REV207_v1_20923"}, "186":{"strain":"3897", "lt":"REV207_v1_20923"}, "187":{"strain":"3898", "lt":"REV207_v1_20923"}, "188":{"strain":"3910", "lt":"REV207_v1_20923"}, "189":{"strain":"3911", "lt":"REV207_v1_20923"}, "190":{"strain":"3912", "lt":"REV207_v1_20923"}, "191":{"strain":"3913", "lt":"REV207_v1_20923"}, "192":{"strain":"3914", "lt":"REV207_v1_20923"}, "193":{"strain":"3915", "lt":"REV207_v1_20923"}, "194":{"strain":"3916", "lt":"REV207_v1_20923"}, "195":{"strain":"3917", "lt":"REV207_v1_20923"}, "196":{"strain":"3918", "lt":"REV207_v1_20923"}, "197":{"strain":"3889", "lt":"REV207_v1_60032"}, "198":{"strain":"3890", "lt":"REV207_v1_60032"}, "199":{"strain":"3892", "lt":"REV207_v1_60032"}, "200":{"strain":"3893", "lt":"REV207_v1_60032"}, "201":{"strain":"3895", "lt":"REV207_v1_60032"}, "202":{"strain":"3896", "lt":"REV207_v1_60032"}, "203":{"strain":"3897", "lt":"REV207_v1_60032"}, "204":{"strain":"3898", "lt":"REV207_v1_60032"}, "205":{"strain":"3910", "lt":"REV207_v1_60032"}, "206":{"strain":"3911", "lt":"REV207_v1_60032"}, "207":{"strain":"3913", "lt":"REV207_v1_60032"}, "208":{"strain":"3914", "lt":"REV207_v1_60032"}, "209":{"strain":"3915", "lt":"REV207_v1_60032"}, "210":{"strain":"3916", "lt":"REV207_v1_60032"}, "211":{"strain":"3917", "lt":"REV207_v1_60032"}, "212":{"strain":"3918", "lt":"REV207_v1_60032"}, "213":{"strain":"3792", "lt":"REV207_v1_60061"}, "214":{"strain":"3794", "lt":"REV207_v1_60066"}, "215":{"strain":"3798", "lt":"REV207_v1_60066"}, "216":{"strain":"3821", "lt":"REV267_v1_10888"}, "217":{"strain":"3819", "lt":"REV267_v1_60172"}, "218":{"strain":"3821", "lt":"REV267_v1_60172"}, "219":{"strain":"3824", "lt":"REV267_v1_60172"}, "220":{"strain":"3820", "lt":"REV267_v1_70138"}, "221":{"strain":"3822", "lt":"REV267_v1_70138"}, "222":{"strain":"3823", "lt":"REV267_v1_70138"}, "223":{"strain":"3826", "lt":"REV267_v1_70138"}, "224":{"strain":"3799", "lt":"REV330_v1_10026"}, "225":{"strain":"3800", "lt":"REV330_v1_10026"}, "226":{"strain":"3801", "lt":"REV330_v1_10026"}, "227":{"strain":"3802", "lt":"REV330_v1_10026"}, "228":{"strain":"3803", "lt":"REV330_v1_10026"}, "229":{"strain":"3804", "lt":"REV330_v1_10026"}, "230":{"strain":"3805", "lt":"REV330_v1_10026"}, "231":{"strain":"3806", "lt":"REV330_v1_10026"}, "232":{"strain":"3808", "lt":"REV330_v1_10026"}, "233":{"strain":"3799", "lt":"REV330_v1_290004"}, "234":{"strain":"3805", "lt":"REV330_v1_290004"}, "235":{"strain":"3803", "lt":"REV330_v1_400001"}}
lst_strain=["3779","3780","3781","3782","3783","3784","3785","3786","3787","3788","3789","3790","3792","3794","3795","3797","3798","3799","3800","3801","3802","3803","3804","3805","3806","3808","3809","3810","3811","3812","3813","3814","3815","3816","3817","3818","3819","3820","3821","3822","3823","3824","3826","3830","3831","3832","3833","3835","3836","3838","3839","3840","3841","3843","3844","3845","3846","3847","3848","3849","3850","3851","3852","3853","3854","3856","3857","3858","3859","3860","3861","3862","3863","3864","3865","3866","3867","3868","3870","3872","3873","3874","3875","3876","3880","3881","3882","3883","3884","3885","3886","3887","3888","3889","3890","3892","3893","3894","3895","3896","3897","3898","3899","3900","3901","3902","3903","3904","3905","3906","3907","3908","3910","3911","3912","3913","3914","3915","3916","3917","3918","3953","3954","3957","3960"]
for lt in lst:
    if lt == "": print("")
    else:
        line = ""
        for strain in lst_strain:
            found = 0
            for num in dico:
                if strain == dico[num]["strain"] and lt == dico[num]["lt"]:
                    found = 1
                    break
            line+=";"+str(found)
        print(line)
exit()


# # dicoHost={"413E50-1":"chaga","402E50-1":"chaga","409E50-1":"chaga","521E56-1":"chaga","384E50-1":"chaga","405E50-1":"chaga","464E53-1":"chaga","468E53-1":"chaga","177E37-1":"chaga","131E34-1":"chaga","121E34-1":"chaga","456E52-1":"chaga","416E50-1":"chaga","284E43-1":"chaga","237E40-1":"chaga","22O28-1":"chaga","434O48-1":"chaga","282E43-1":"chaga","417E50-1":"chaga","448O51-1":"chaga","217E38-1":"chaga","455E52-1":"chaga","466E53-1":"chaga","191E37-1":"chaga","115E34-1":"chaga","511E55-1":"chaga","120E34-1":"chaga","495E54-1":"chaga","496E54-1":"chaga","277E43-1":"chaga","193E37-1":"chaga","249E41-1":"chaga","382E49-1":"chaga","141O35-1":"chaga","141E35-1":"chaga","230E39-1":"chaga","137E35-1":"chaga","393E50-1":"chaga","489E54-1":"chaga","424E50-1":"chaga","501E54-1":"chaga","381E49-1":"chaga","242E40-1":"chaga","275E43-1":"chaga","142E35-1":"chaga","168E36-1":"chaga","236O40-1":"chaga","199E37-1":"chaga","150E35-1":"chaga","24E30-2":"crass","66E30-1":"crass","24E35-2":"crass","64E30-1":"crass","44E38-2":"crass","36E38-1":"crass","41E34-2":"crass","75E35-1":"crass","44E38-1":"crass","34E29-1":"crass","70E35-5a":"crass","234P7b":"crass","70E35-2":"crass","294E48-1":"crass","70E38-1":"crass","70E37-6":"crass","70E37-1":"crass","70E35-6":"crass","234P1":"crass","11E33-1":"crass","82E32-1":"crass","82E32-3":"crass","38E33-6a":"crass","82E33-2":"crass","82E32-2":"crass","19E33-1":"crass","219E41-1":"crass","219E41-2":"crass","234P8":"crass","431E48-2":"crass","431E46-1":"crass","431E45-1":"crass","207E29-1":"crass","91E28-1a":"crass","98E28-6a":"crass","18E29-1":"crass","12E28-1":"crass","264E42-1":"crass","51E28-1":"crass","51E28-4":"crass","340E47-2":"crass","5P1a":"crass","69E27-1":"crass","23E28-1":"crass","184E37-3a":"crass","184E37-1":"crass","104E43-1":"crass","103E44-1":"crass","14E30-1":"crass","14E30-2":"crass","207E48-1":"crass","144E46-1":"crass","252E42-2":"crass","268E42-1":"crass","6E35-1a":"crass","15E36-1":"crass","172P1":"crass"}
# # dicoLtTOorg = {}
# # for file in os.listdir("/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/FAA_phages"):
# #     if ".faa" in file:
# #         dicoTMP = make_fasta_dict("/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/FAA_phages/"+file)
# #         for key in dicoTMP:
# #             # dicoLtTOorg[key.split("|")[0]] = file.replace("Vibrio_phage_","").replace(".faa","")
# #             print(key.split("|")[0]+";"+dicoHost[file.replace("Vibrio_phage_","").replace(".faa","")])
# # exit()

# # dicoPairwise = {}
# # for file in os.listdir("/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/family_rbh_25_50"):
# #     pathRBH = "/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/family_rbh_25_50/"+file
# #     lstLines = read_file(pathRBH)
# #     for line in lstLines:
# #         splitLine = line.split("\t")
# #         ref = splitLine[0].split("|")[0]
# #         target = splitLine[1]
# #         targetOrg = dicoLtTOorg[target.split("|")[0]]
# #         pident = float(splitLine[2])
# #         if ref not in dicoPairwise: dicoPairwise[ref] = {}
# #         if targetOrg not in dicoPairwise[ref] or pident>dicoPairwise[ref][targetOrg]: dicoPairwise[ref][targetOrg] = pident
# #     dicoMeanPerHost = {}
# #     for ltRef in dicoPairwise:
# #         lst_pident_crass = []
# #         lst_pident_chaga = []
# #         for targetOrg in dicoPairwise[ltRef]:
# #             if targetOrg in ["91E28-1a","18E29-1","12E28-1","98E28-6a","207E29-1"]: lst_pident_crass.append(dicoPairwise[ltRef][targetOrg])
# #             elif targetOrg in ["115E34-1","455E52-1","217E38-1","120E34-1","511E55-1","191E37-1","466E53-1"]: lst_pident_chaga.append(dicoPairwise[ltRef][targetOrg])
# #             # if dicoHost[targetOrg] == "crass": lst_pident_crass.append(dicoPairwise[ltRef][targetOrg])
# #             # else: lst_pident_chaga.append(dicoPairwise[ltRef][targetOrg])
# #         # if len(lst_pident_crass)==4:
# #         print(ltRef+";"+str(len(lst_pident_crass))+";"+str(np.mean(lst_pident_crass))+";"+str(len(lst_pident_chaga))+";"+str(np.mean(lst_pident_chaga)))
# # exit()


# dicoHost={"413E50-1":"chaga","402E50-1":"chaga","409E50-1":"chaga","521E56-1":"chaga","384E50-1":"chaga","405E50-1":"chaga","464E53-1":"chaga","468E53-1":"chaga","177E37-1":"chaga","131E34-1":"chaga","121E34-1":"chaga","456E52-1":"chaga","416E50-1":"chaga","284E43-1":"chaga","237E40-1":"chaga","22O28-1":"chaga","434O48-1":"chaga","282E43-1":"chaga","417E50-1":"chaga","448O51-1":"chaga","217E38-1":"chaga","455E52-1":"chaga","466E53-1":"chaga","191E37-1":"chaga","115E34-1":"chaga","511E55-1":"chaga","120E34-1":"chaga","495E54-1":"chaga","496E54-1":"chaga","277E43-1":"chaga","193E37-1":"chaga","249E41-1":"chaga","382E49-1":"chaga","141O35-1":"chaga","141E35-1":"chaga","230E39-1":"chaga","137E35-1":"chaga","393E50-1":"chaga","489E54-1":"chaga","424E50-1":"chaga","501E54-1":"chaga","381E49-1":"chaga","242E40-1":"chaga","275E43-1":"chaga","142E35-1":"chaga","168E36-1":"chaga","236O40-1":"chaga","199E37-1":"chaga","150E35-1":"chaga","24E30-2":"crass","66E30-1":"crass","24E35-2":"crass","64E30-1":"crass","44E38-2":"crass","36E38-1":"crass","41E34-2":"crass","75E35-1":"crass","44E38-1":"crass","34E29-1":"crass","70E35-5a":"crass","234P7b":"crass","70E35-2":"crass","294E48-1":"crass","70E38-1":"crass","70E37-6":"crass","70E37-1":"crass","70E35-6":"crass","234P1":"crass","11E33-1":"crass","82E32-1":"crass","82E32-3":"crass","38E33-6a":"crass","82E33-2":"crass","82E32-2":"crass","19E33-1":"crass","219E41-1":"crass","219E41-2":"crass","234P8":"crass","431E48-2":"crass","431E46-1":"crass","431E45-1":"crass","207E29-1":"crass","91E28-1a":"crass","98E28-6a":"crass","18E29-1":"crass","12E28-1":"crass","264E42-1":"crass","51E28-1":"crass","51E28-4":"crass","340E47-2":"crass","5P1a":"crass","69E27-1":"crass","23E28-1":"crass","184E37-3a":"crass","184E37-1":"crass","104E43-1":"crass","103E44-1":"crass","14E30-1":"crass","14E30-2":"crass","207E48-1":"crass","144E46-1":"crass","252E42-2":"crass","268E42-1":"crass","6E35-1a":"crass","15E36-1":"crass","172P1":"crass"}
# dicoAnnot={"cluster0001":"large_terminase","cluster0002":"putative_baseplate_wedge_protein","cluster0004":"unkwown","cluster0005":"baseplate_J_protein","cluster0007":"unkwown","cluster0008":"neck_protein","cluster0009":"unkwown","cluster0011":"putative_tail_sheath_protein","cluster0012":"DNA_helicase","cluster0014":"putative_gp55_tail-fiber","cluster0015":"nucleic_acid-binding","cluster0016":"unkwown","cluster0017":"unkwown","cluster0018":"putative_baseplate_assembly_protein","cluster0019":"minor_capsid_protein","cluster0021":"baseplate_protein_J-like_protein","cluster0023":"unkwown","cluster0027":"DNA_helicase","cluster0029":"DNA_primase-helicase","cluster0030":"unkwown","cluster0031":"unkwown","cluster0032":"tail_assembly_chaperone","cluster0035":"unkwown","cluster0036":"head","cluster0037":"unkwown","cluster0039":"putative_structural_protein","cluster0040":"putative_lytic_transglycosylase","cluster0041":"putative_head_protein","cluster0042":"unkwown","cluster0043":"major_capsid_protein","cluster0044":"portal","cluster0045":"unkwown","cluster0048":"DNA_nuclease","cluster0051":"tail-completion_protein","cluster0053":"putative_tail_tube_protein","cluster0055":"unkwown","cluster0057":"head_completion_adaptor","cluster0058":"unkwown","cluster0059":"head-closure_protein","cluster0060":"unkwown","cluster0061":"unkwown","cluster0063":"unkwown","cluster0064":"metallo-dependent_phosphatase-like_protein"}
# dicoLT = {}
# cpt = 0
# typeseq="faa"
# dicoSeq = {}
# for file in os.listdir("/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/"+typeseq.upper()+"_phages"):
#     if typeseq in file:
#         dicoTMP = make_fasta_dict("/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/"+typeseq.upper()+"_phages/"+file)
#         for key in dicoTMP:
#             dicoSeq[key.split("|")[0]] = dicoTMP[key]
# # dicoCluster = load_json("/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/all_phages_rbh_25_50/all_phages_rbh_25_50.json")
# # dicoCluster = load_json("/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/all_phages_rbh_25_80.json")
# dicoCluster = load_json("/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/family_rbh_25_80.json")
# for clusterNum in dicoCluster:
#     setCrass = set()
#     setChaga = set()
#     for prot in dicoCluster[clusterNum]: # "VP495E541_p0026|-1 [Vibrio_phage_495E54-1]",
#         orgName = prot.split("[")[1].replace("Vibrio_phage_","").replace("]","")
#         if orgName=="115E34-1":  dicoLT[clusterNum] = prot.split("|")[0]
#         if orgName in ["91E28-1a","18E29-1","12E28-1","98E28-6a","207E29-1"]: setCrass.add(orgName)
#         elif orgName in ["115E34-1","455E52-1","217E38-1","120E34-1","511E55-1","191E37-1","466E53-1"]: setChaga.add(orgName)
#     if len(setCrass) == 5 and len(setChaga) == 7: #> 1:
#     # if len(setCrass)+len(setChaga)>1:
#         # print(clusterNum)
#         cpt += 1
#         pathFASTA = "/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/individual_tree/"+clusterNum+".fasta"
#         pathALIGN = "/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/individual_tree/"+clusterNum+"_famsa.fasta"
#         if not os.path.isfile(pathFASTA):
#             OUT = open(pathFASTA,'w')
#             for prot in dicoCluster[clusterNum]:
#                 orgName = prot.split("[")[1].replace("Vibrio_phage_","").replace("]","")
#                 lt = prot.split("|")[0]
#                 OUT.write(">"+orgName+"|"+lt+"\n"+dicoSeq[lt]+"\n")
#             OUT.close()
#         pathTREE = "/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/individual_tree/"+clusterNum
#         # famsa
#         if not os.path.isfile(pathALIGN):
#             os.system("famsa "+pathFASTA+" "+pathALIGN+" 2>/dev/null")

#         # # pairwise identity
#         # dicoALIGN = make_fasta_dict(pathALIGN)
#         # lst_pident_crass = []
#         # lst_pident_chaga = []
#         # lst_pident_both = []
#         # lstLT= []
#         # for orgName1 in dicoALIGN:
#         #     lstLT.append(orgName1.split("|")[1])
#         #     for orgName2 in dicoALIGN:
#         #         if orgName1 != orgName2:
#         #             nbMatch = 0
#         #             for i in range(len(dicoALIGN[orgName1])):
#         #                 if dicoALIGN[orgName1][i] == dicoALIGN[orgName2][i]: nbMatch += 1
#         #             pident = (nbMatch*100) / len(dicoALIGN[orgName1])
#         #             if dicoHost[orgName1.split("|")[0]] == dicoHost[orgName2.split("|")[0]]:
#         #                 if dicoHost[orgName1.split("|")[0]] == "chaga": lst_pident_chaga.append(pident)
#         #                 else: lst_pident_crass.append(pident)
#         #             else: lst_pident_both.append(pident)
#         # if len(lst_pident_crass)==0: min_crass = -1
#         # else: min_crass = min(lst_pident_crass)
#         # if len(lst_pident_chaga)==0: min_chaga = -1
#         # else: min_chaga = min(lst_pident_chaga)
#         # if len(lst_pident_both)==0: max_both = -1
#         # else: max_both = max(lst_pident_both)
#         # print(clusterNum+";"+str(len(setCrass))+";"+str(len(setChaga))+";"+str(min_crass)+";"+str(min_chaga)+";"+str(max_both)+";"+",".join(lstLT))


#         # iqtree
#         if not os.path.isfile(pathTREE+".nwk"):
#             if typeseq == "faa": os.system("iqtree2 -s "+pathALIGN+" -T 12 --mem 25GB -m LG -B 1000 --seqtype AA --prefix "+pathTREE+" --keep-ident > /dev/null 2>&1")
#             else: os.system("iqtree2 -s "+pathALIGN+" -T 12 --mem 25GB -m GTR -B 1000 --seqtype DNA --prefix "+pathTREE+" --keep-ident > /dev/null 2>&1")
#             os.system("rm -f "+pathTREE+".splits.nex "+pathTREE+".log "+pathTREE+".mldist "+pathTREE+".bionj "+pathTREE+".ckp.gz "+pathTREE+".contree "+pathTREE+".iqtree")
#             shutil.move(pathTREE+".treefile",pathTREE+".nwk")
#         # RENDERING tree
#         tree = Tree(pathTREE+".nwk")
#         for node in tree.traverse("postorder"):
#             # Do some analysis on node
#             node.name = node.name.split("|")[0]
#         ts = TreeStyle()
#         ts.show_branch_support = True
#         for n in tree.traverse():
#             nstyle = NodeStyle()
#             nstyle["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
#             nstyle["hz_line_type"] = 0
#             # if n.support != 1.0:
#             #     nstyle["vt_line_width"] = n.support/20
#             #     nstyle["hz_line_width"] = n.support/20    
#             #     nstyle["size"] = 0
#             try: host = dicoHost[n.name]
#             except KeyError: n.set_style(nstyle) ; continue
#             if host == "crass":
#                 if n.name in ["91E28-1a","18E29-1","12E28-1","98E28-6a","207E29-1"]: nstyle["fgcolor"] = "#ff5599"
#                 else: nstyle["fgcolor"] = "#ffd5f6"
#             elif host == "chaga":
#                 if n.name in ["115E34-1","455E52-1","217E38-1","120E34-1","511E55-1","191E37-1","466E53-1"]: nstyle["fgcolor"] = "#37c871"
#                 else: nstyle["fgcolor"] = "#d7f4e3"
#             # nstyle["shape"] = "sphere"
#             nstyle["size"] = 10
#             n.set_style(nstyle)
#         ts.show_leaf_name = True
#         ts.show_branch_support = True
#         ts.title.add_face(TextFace(dicoLT[clusterNum], fsize=12), column=0)
#         # tree.show(tree_style=ts)
#         tree.render("/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/individual_tree/"+dicoLT[clusterNum]+".svg", tree_style=ts)
#         # tree.render("/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/individual_tree/"+dicoLT[clusterNum]+".png", units="mm", tree_style=ts)
#         # input("*****")
# # print(cpt)
# exit()
# nstyle["shape"] = "sphere"
# nstyle["size"] = 10
# nstyle["fgcolor"] = "darkred"

# # ***** PLOT SIMILARITY MATRIX ***** #
# dicoSimMA = load_json("/mnt/c/Users/dgoudenege/Taf/Vibrio_chagasii/VIRIDIC/viridic.json")
# df = pd.DataFrame(dicoSimMA)
# # cg = sns.clustermap(df, cmap = 'crest', figsize = (50, 50), tree_kws = {'linewidths': 2.5}, dendrogram_ratio = 0.15, annot_kws = {"size": 35 / np.sqrt(len(df))}, cbar_kws = {'label': 'similarity %'})
# cmap = sns.color_palette("light:#d40000", as_cmap=True)
# cg = sns.clustermap(df, cmap=cmap, figsize=(50, 50), tree_kws={'linewidths': 2.5}, dendrogram_ratio=0.15, annot_kws={"size": 35 / np.sqrt(len(df))}, cbar_kws={'label': 'similarity %'}, linewidths=0.0, rasterized=True)
# # Retrieve ordered ticks label
# newColums = df.columns[cg.dendrogram_col.reordered_ind]
# newIndexs = df.index[cg.dendrogram_row.reordered_ind]
# newData = df.loc[newIndexs, newColums]
# orderedOrg = list(newData.keys())
# # Plot clustered heatmap
# cg.ax_cbar.tick_params(labelsize=40)
# cg.ax_cbar.yaxis.label.set_size(50)
# plt.savefig("/tmp/matrix.png", dpi=300)
# plt.savefig("/tmp/matrix.svg")
# # ***** WRITE SIMILARITY MATRIX ***** #
# pathOUTmatrix = "/tmp/matrix.tsv"
# OUT = open(pathOUTmatrix, 'w')
# header = "Organism"
# for orgName in orderedOrg:
#     header += "\t"+orgName
# OUT.write(header+"\n")
# for orgName1 in orderedOrg:
#     line = orgName1
#     for orgName2 in orderedOrg:
#         line += "\t"+str(dicoSimMA[orgName1][orgName2]).replace(".", ",")
#     OUT.write(line+"\n")
# OUT.close()


# exit()

# dicoANNOT = {"cluster0001":"STRUCT","cluster0002":"STRUCT","cluster0003":"UNK","cluster0004":"UNK","cluster0005":"STRUCT","cluster0006":"UNK","cluster0007":"UNK","cluster0008":"STRUCT","cluster0009":"UNK","cluster0010":"MT","cluster0011":"STRUCT","cluster0012":"DNA","cluster0013":"UNK","cluster0014":"STRUCT","cluster0015":"OTHERS","cluster0016":"UNK","cluster0017":"UNK","cluster0018":"STRUCT","cluster0019":"STRUCT","cluster0020":"UNK","cluster0021":"STRUCT","cluster0022":"UNK","cluster0023":"UNK","cluster0024":"UNK","cluster0025":"DNA","cluster0026":"UNK","cluster0027":"DNA","cluster0028":"UNK","cluster0029":"DNA","cluster0030":"UNK","cluster0031":"UNK","cluster0032":"STRUCT","cluster0033":"UNK","cluster0034":"UNK","cluster0035":"UNK","cluster0036":"STRUCT","cluster0037":"UNK","cluster0038":"UNK","cluster0039":"STRUCT","cluster0040":"OTHERS","cluster0041":"STRUCT","cluster0042":"UNK","cluster0043":"STRUCT","cluster0044":"STRUCT","cluster0045":"UNK","cluster0046":"STRUCT","cluster0047":"UNK","cluster0048":"DNA","cluster0049":"UNK","cluster0050":"UNK","cluster0051":"STRUCT","cluster0052":"DNA","cluster0053":"STRUCT","cluster0054":"UNK","cluster0055":"UNK","cluster0056":"STRUCT","cluster0057":"STRUCT","cluster0058":"UNK","cluster0059":"STRUCT","cluster0060":"UNK","cluster0061":"UNK","cluster0062":"UNK","cluster0063":"UNK","cluster0064":"OTHERS","cluster0065":"UNK","cluster0066":"UNK","cluster0067":"UNK","cluster0068":"UNK","cluster0069":"UNK","cluster0070":"UNK","cluster0071":"UNK","cluster0072":"UNK","cluster0073":"MT","cluster0074":"UNK","cluster0075":"UNK","cluster0076":"UNK","cluster0077":"UNK","cluster0078":"OTHERS","cluster0079":"STRUCT","cluster0080":"UNK","cluster0081":"DNA","cluster0082":"UNK","cluster0083":"STRUCT","cluster0084":"UNK","cluster0085":"UNK","cluster0086":"UNK","cluster0087":"UNK","cluster0088":"UNK","cluster0089":"UNK","cluster0090":"UNK","cluster0091":"UNK","cluster0092":"UNK","cluster0093":"UNK","cluster0094":"STRUCT","cluster0095":"UNK","cluster0096":"OTHERS","cluster0097":"UNK","cluster0098":"UNK","cluster0099":"UNK","cluster0100":"UNK","cluster0101":"STRUCT","cluster0102":"UNK","cluster0103":"OTHERS","cluster0104":"UNK","cluster0105":"UNK","cluster0106":"UNK","cluster0107":"DNA","cluster0108":"UNK","cluster0109":"UNK","cluster0110":"UNK","cluster0111":"UNK","cluster0112":"UNK","cluster0113":"UNK","cluster0114":"UNK","cluster0115":"OTHERS","cluster0116":"UNK","cluster0117":"UNK","cluster0118":"UNK","cluster0119":"UNK","cluster0120":"UNK","cluster0121":"UNK","cluster0122":"UNK","cluster0123":"OTHERS","cluster0124":"UNK","cluster0125":"UNK","cluster0126":"UNK","cluster0127":"UNK","cluster0128":"DNA","cluster0129":"UNK","cluster0130":"UNK","cluster0131":"UNK","cluster0132":"UNK","cluster0133":"UNK","cluster0134":"STRUCT","cluster0135":"UNK","cluster0136":"UNK","cluster0137":"UNK","cluster0138":"STRUCT","cluster0139":"UNK","cluster0140":"DNA","cluster0141":"UNK","cluster0142":"UNK","cluster0143":"MT","cluster0144":"UNK","cluster0145":"UNK","cluster0146":"UNK","cluster0147":"UNK","cluster0148":"UNK","cluster0149":"UNK","cluster0150":"UNK","cluster0151":"DNA","cluster0152":"DNA","cluster0153":"UNK","cluster0154":"DNA","cluster0155":"UNK","cluster0156":"UNK","cluster0157":"UNK","cluster0158":"UNK","cluster0159":"UNK","cluster0160":"UNK","cluster0161":"UNK","cluster0162":"UNK","cluster0163":"UNK","cluster0164":"UNK","cluster0165":"UNK","cluster0166":"UNK","cluster0167":"UNK","cluster0168":"UNK","cluster0169":"UNK","cluster0170":"UNK","cluster0171":"UNK","cluster0172":"UNK","cluster0173":"UNK","cluster0174":"UNK","cluster0175":"UNK","cluster0176":"UNK","cluster0177":"UNK","cluster0178":"STRUCT","cluster0179":"STRUCT","cluster0180":"STRUCT","cluster0181":"UNK","cluster0182":"UNK","cluster0183":"UNK","cluster0184":"UNK","cluster0185":"OTHERS","cluster0186":"UNK","cluster0187":"UNK","cluster0188":"UNK","cluster0189":"UNK","cluster0190":"UNK","cluster0191":"UNK","cluster0192":"UNK","cluster0193":"UNK","cluster0194":"UNK","cluster0195":"UNK","cluster0196":"UNK","cluster0197":"MT","cluster0198":"UNK","cluster0199":"UNK","cluster0200":"UNK","cluster0201":"DNA","cluster0202":"UNK","cluster0203":"UNK","cluster0204":"UNK","cluster0205":"UNK","cluster0206":"UNK","cluster0207":"UNK","cluster0208":"UNK","cluster0209":"UNK","cluster0210":"MT","cluster0211":"UNK","cluster0212":"UNK","cluster0213":"UNK","cluster0214":"UNK","cluster0215":"UNK","cluster0216":"UNK","cluster0217":"UNK","cluster0218":"UNK","cluster0219":"UNK","cluster0220":"UNK","cluster0221":"UNK","cluster0222":"UNK","cluster0223":"DNA","cluster0224":"UNK","cluster0225":"UNK","cluster0226":"UNK","cluster0227":"UNK","cluster0228":"UNK","cluster0229":"UNK","cluster0230":"UNK","cluster0231":"UNK","cluster0232":"UNK","cluster0233":"UNK","cluster0234":"UNK","cluster0235":"UNK","cluster0236":"UNK","cluster0237":"UNK","cluster0238":"UNK","cluster0239":"STRUCT","cluster0240":"UNK","cluster0241":"OTHERS","cluster0242":"UNK","cluster0243":"UNK","cluster0244":"UNK","cluster0245":"UNK","cluster0246":"UNK","cluster0247":"UNK","cluster0248":"OTHERS","cluster0249":"UNK","cluster0250":"UNK","cluster0251":"UNK","cluster0252":"UNK","cluster0253":"UNK"}
# dicoCOUNTANNOT = {}
# dicoColorCluster = {}
# for clusterNum in dicoANNOT:
#     annot = dicoANNOT[clusterNum]
#     try: dicoCOUNTANNOT[annot] += 1
#     except: dicoCOUNTANNOT[annot] = 1
# for annot in dicoCOUNTANNOT:
#     if annot == "DNA":  gradient = linear_gradient("#00aad4", finish_hex="#d5f6ff", n=dicoCOUNTANNOT[annot])[0]
#     elif annot == "MT": gradient = linear_gradient("#ff0000", finish_hex="#ffd5d5", n=dicoCOUNTANNOT[annot])[0]
#     elif annot == "STRUCT": gradient = linear_gradient("#a05a2c", finish_hex="#e9c6af", n=dicoCOUNTANNOT[annot])[0]
#     elif annot == "OTHERS": gradient = linear_gradient("#2ca05a", finish_hex="#d7f4e3", n=dicoCOUNTANNOT[annot])[0]
#     elif annot == "UNK": gradient = linear_gradient("#cccccc", finish_hex="#cccccc", n=dicoCOUNTANNOT[annot])[0]
#     random.shuffle(gradient)
#     cpt = 0
#     for clusterNum in dicoANNOT:
#         if dicoANNOT[clusterNum] == annot:
#             dicoColorCluster[clusterNum] = gradient[cpt]
#             cpt += 1

# dicoProt = {}
# pathDir = "/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/phage_family_fig4d"
# for faa in os.listdir(pathDir):
#     if ".faa" not in faa: continue
#     dicoFAA = make_fasta_dict("/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/phage_family_fig4d/"+faa)
#     for key in dicoFAA: dicoProt[key.split('|')[0]] = dicoFAA[key]

# dicoColor = {}
# dicoBorder = {}
# setCore = set()
# setLTGenus13 = set()
# dicoCluster = load_json("/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/phage_family_fig4d/rbh_25_50.json")
# for clusterNum in dicoCluster:
#     setAll = set()
#     setGenus13 = set()
#     setOthers = set()
#     for prot in dicoCluster[clusterNum]:
#         orgName = prot.split("[")[1].replace("Vibrio_phage_","").replace("]","")
#         setAll.add(orgName)
#         if orgName in ["91E28-1a","18E29-1","12E28-1","98E28-6a"]:
#             setGenus13.add(orgName)
#         else:
#             setOthers.add(orgName)
#     if len(setAll) == 1: genocomp = "single"
#     else:
#         if len(setAll) == 12: genocomp = "core"
#         else: genocomp = "variable"

#     for prot in dicoCluster[clusterNum]:
#         lt = prot.split("|")[0]
#         orgName = prot.split("[")[1].replace("Vibrio_phage_","").replace("]","")
#         dicoColor[lt] = dicoColorCluster[clusterNum]
#         if genocomp == "core": dicoBorder[lt] = 5
#         else: dicoBorder[lt] = 0

# pathDir = "/mnt/c/Users/dgoudenege/Desktop/reviewing_chaga/phage_family_fig4d"
# for gff in os.listdir(pathDir):
#     if ".gff" not in gff: continue
#     orgName = gff.replace(".gff", "")
#     pathPNG = "/tmp/"+orgName+".png"
#     pathSVG = "/tmp/"+orgName+".svg"
#     features = []
#     # ***** PARSE GFF ***** #
#     dicoGFF = make_gff_dict(pathIN=pathDir+"/"+gff, ext=".gff")
#     startRegion = 0
#     endRegion = 0
#     # ***** BROWSE GENES ***** #
#     for geneType in dicoGFF[orgName]:
#         if geneType != 'length':
#             for geneEntry in dicoGFF[orgName][geneType]:
#                 if geneType == "CDS":
#                     color = "#2a7fff"
#                 elif geneType == "tRNA":
#                     color = "#37c8ab"
#                 else:
#                     continue
#                 if 'locus_tag' in geneEntry['attributes']:
#                     lt = geneEntry['attributes']['locus_tag']
#                     if lt in dicoColor:
#                         geneFeature = GraphicFeature(start=geneEntry['start'], end=geneEntry['end'], strand=int(geneEntry['strand']+"1"), color=dicoColor[lt], linewidth=dicoBorder[lt])                  
#                     else:
#                         geneFeature = GraphicFeature(start=geneEntry['start'], end=geneEntry['end'], strand=int(geneEntry['strand']+"1"), color=color, linewidth=dicoBorder[lt])
#                     features.append(geneFeature)
#     # ***** PLOT GENES ***** #
#     seqLen = 55000 # dicoGFF[orgName]['length']
#     record = GraphicRecord(sequence_length=seqLen+int(seqLen/10), features=features, first_index=startRegion-100)
#     ax, _ = record.plot(figure_width=50)
#     # ax.figure.savefig(pathPNG, dpi=300)
#     ax.figure.savefig(pathSVG)
#     plt.close('all')


# os.system("svg_stack.py --direction=v --margin=1 /tmp/Vibrio_phage_115E34-1.svg /tmp/Vibrio_phage_91E28-1a.svg /tmp/Vibrio_phage_18E29-1.svg /tmp/Vibrio_phage_12E28-1.svg /tmp/Vibrio_phage_98E28-6a.svg /tmp/Vibrio_phage_455E52-1.svg /tmp/Vibrio_phage_217E38-1.svg /tmp/Vibrio_phage_120E34-1.svg /tmp/Vibrio_phage_511E55-1.svg /tmp/Vibrio_phage_207E29-1.svg /tmp/Vibrio_phage_191E37-1.svg /tmp/Vibrio_phage_466E53-1.svg > /mnt/c/Users/dgoudenege/Desktop/merge.svg")
# import cairosvg ; cairosvg.svg2png(url="/mnt/c/Users/dgoudenege/Desktop/merge.svg", write_to="/mnt/c/Users/dgoudenege/Desktop/merge.png") 
# os.system("rm -f /mnt/c/Users/dgoudenege/Desktop/Vibrio_phage*.svg")

# exit()



# BLAST to order reads
# AATACGGCATGAACTAA (GGCATGAACTAA)
pathREF = "/mnt/c/Users/dgoudenege/Taf/Nanopore_phages/PICMI_Vibrio_chagasii_34_P_115.fna"
pathREADS = "/mnt/c/Users/dgoudenege/Taf/Nanopore_phages/PICMI_reads.fasta"
os.system("blastn -task blastn -query "+pathREADS+" -subject "+pathREF+" -max_hsps 1 -outfmt \"6 delim=; qseqid sstart send\" -evalue 0.001 1>/tmp/blast.out 2>/dev/null")
IN = open("/tmp/blast.out",'r')
lstLines = IN.read().split("\n")[:-1]
IN.close()
seqPICMI = "AATACGGCATGAACTAAACTAAACAGTCTAACCGTATGAATAAAAAGAGCTTATTTGAATAAAATAGGCTCTTTTTCTATTTAACTCCCACAATAACCCCCACAACTAAAAAACAACTCCCACATCACCACATGAAAAAATATCATGCGATCACACTTTTTATTTCTGCAAACAAAAAACACTTAGGTCAACAGGACACTTAGGACATTCCCTTATCTGGTAAAGGCTGCGATTGTCCTAAACTCTTTTCTACTTAGAACAAACTATCTATTAAACAGGACACATCAAAGAACAGTTGCATTCTTAGAGGCCAAAAAGGGGAATATTTGAAAGTCATAGAGTCATTCTCGGAGTTCACGTCTCGACCTTGTTACAGTTCGATTGCTTTAGTCCACGGGGGTAAGGAAATTGCTTTTTACTGCTAAGGTTCTTTTAATACTTCGACATCACTTGAGGGTATTCTCGTCTTGTTTCCACACGATATAAAAAAGGACGCTAGGAAACGTCCTGTTAGCTATTTACCTTTTTAGGTGGAAATCCTCACGTTTACAATTAAATCAGCTTACATGCCAGTTTAGATGGCTTAACGCCGGCGATAGCTTTCTTCTTCGCGAACTTAAAGATACAACCCTTTCTATCGACTGGTATTTTCATCAAATTGGGGTAAATACACTCCCAAAATTTAGTGATAGATTTTCGAACACCATAAAATTTAGATATCTATACGCCCAACCTGAGTGAAAAAGAACTATGAAAGGCTTGAGTGCATATGCACAGTATATGAACAATCGCCATGTACAACATAATTAATACTGGTAGGATTCTAATCTCGCTACACACATTGAGGAAAATATGAAAAGTAAGATCTTACTATTAGCAGTCCCTTTAGCTACGGTTCAAGGTTGTGCAACTAACCCTACTGAAGTATCAATGGATATAAACATATCCATTGATAAAAGAGTGTTAAACTTTGAGTCATTAAAAAGTTCATTTAGCTACAACTGCCTCCCACTTTCACAGCATGAACCAGGGGAAGTTAAAGCTAGGTTTTCGCCTAAGCTATATAACCCACCTGTTGGCTATATCTCTATTCAAGAACCTTCAACCGAACTTAGACAAGGAATAAAACTCTCACAAACTGAATCACTTGAATCACTCGCGTCTGTGAACATGTTCACAGCTGAAATAAACGGTATCCACGGTTCAAATCAAACGATTTATCGTTATGACAATTCCATCATTAAACAGTTTATTGAAGAATGCCGAATATCTACGGTTAATGAATGGAAAATAGAAATTGAAAAGGTGAAGAAGGATGAAATTCAATCCAAATTAAAAGCAAAACAACACCAAAACAAGGTGACAACTTCAATCAATAACGCCCTGAAAAAACATAATATAAAGGGGGTTAAAAATGGTGGCGAACTAGCAATACAAAATGTTTTGTTTGAAAAACGCCATCTAGACAAGGATTCTTACCGAAGTTTTCTCGAAAACATCAAGACTGCTTATGTGATAGATTTCTCTCAATATGTTGTCAAACAAAGCTTAGGCAATAACAAATATGTACTGAGACACCTAATCGAAGACTACTATACCCCTGCTATTACCCCTCCTTTAGCTATTATCCTAAATACAAACAAACTATTATTTGAAGGTGAAACGCCAAGTAAGAGTGTCGTTCTGTATACTGGAATAGAACTACACCAAACTGTAGCAGGTTCAAATAGACAACTGGTAGGCTTAACCACTTTAGATTAAGTAAAAGTCCGCACTTTGGCGGACTCTCTTTTTATTGCTGGTATTGATAGAATTCATCAAAGCATCGTTGCATATCTTCTAAAGGTATTGTCCTATAGCAGTTAAGCCTTTGCCCTGAGAAATTGCTCTTTCCTCCATCAGTAACCACACCTAATGTTTGCTTTAAATACTTACCTAATGATTTTGACGTTTTACGCTCATACCCATTAACTTTCATACTTTCACAATACGCAATGTAACTATTAAAAAAAACTAGCTTTGGAATAGTTTTGCTATCGCCTAAGTCATAGACTTGAGCATCATCAGATAAAAAACCATCTTCTAAAACTTTCTGAAACCATTCAGCTTCACTTGGTAATGATTCCATAATCTGTGCTTGTAACGCATCAGTCATTGGTGCTTTATGCAAATTATAAGTGCTAATTTCATAACGCTTAAGAAAGTTAAAAAACTGCCCTATCGCCTCTGGTGAATCAATATCATTATGTAGCGCACGAAAATATTCATTATCTTGTTTTTTACAATCCGACACCTCAAGTACAAAAAAGCGTCTTTCTTGCTGATGTGCTGGTACAGCCCAATCATTATTTGTACACATCAAATAACGGTGATAACTATCAACCTCAATCGCGTCCTTTCCTTTGGCCTCAATAGTTATCGTACTATCAGTAATCAAGGTTCGAAACTTGCCAGCATCTTTATTTGAACCACTCCAAAAAGCTTCTTCAATTGTCACGAAGAGCTTATTTGCCAAGTGCGAGTTAAATGACCCCAGCAAGTGTTTACTGTCTTGAATCCTCATTGCGTACTGTCCTAACAACTTTTCCATTAACACAGAAACGGTACTTTTCCCTGTCCCTCGACCTTCCGACTTAAGCACAACGCAAACGCCTGTTTTCTTTTCCGGTCTCTGTATAATTTGAGCCATCCAAGCCAATAGGTAATGAAAGTGATCTTGATTGCCGTTGCAGACAACATGGAGCAAATGCCACATGATACGATCTAGCTTGTTATCACCGTCAACCGCCTCAATTTCAAAGCCTGACCACATATTGAGAAAACGCGACTCTGTTCTGTTTGATGGATCAAATACCACTCCTTGAAAAGTACTGCGGCTTCCTGCTTTCATCCATAATTCAAAGCAGTTGTATTGTTTAACTTGATTTTGGTCGTTTAAATAAGGGAAGTTCAGGTTTGCATATAACTGCCGCTTTTGCGGTACAGCAGAAAATGTATAACCCGTATTGTCTTTATGATCGGTTTCACGTTCACACACCCACACTTTCCCCCCTTGCTGAACAATAGAGTGGGTAAGGTTCAAAGTATCAATTATATGTTGCTTGCTGTCGCCTAGTGGTTTTGCACCTACCGTCTCGCAGAACATGTTCAATTGCTCATAAGGCGTTAAGTCTTTCTCGATCTCGTTCAAAAATCTTTCAATCTTTGGTGATAATATTGCATTCATGTTTTACACTCTTGTTCAGCTGGTGGCTTCCTGTCGCCAAACTGAGCCACCAGCTATAGATTTTTAACAAAATTTAGATCCGCTCTGCTCATCAGTATGAGCGGTTTCTTTTTCGTAGCACTCGTCCAACTTGGAATACACGGCATATCTAACCGCCTTTTCCAACCCACGAACAAAGAAGTCCAACGCTAAACAACGTTCCCATACCATGAGAAATTTTTGTTCATGGCTAGTGATAGGCTTCCCATTAGCAAACATGACTGGAATTATGTTGCCCTCAATAAAACCGAAAGAATATCCACTTATCTCCATGACCTTATGAATGTTGTCTTGGCTTGCTATTGGTCGATTATCTCGATTACATTCAAACATATTGCTTAACTCCTAAAATACGATGAGCACATTCGAGGAAACGAACCTTTGCATTGATACGCTTTAATATACGTTGTGACTTTAATGATGAGGGGTTAGGCATTAGCAGAGCCTCCTATGTTATTATCAATCCAAGATTGAATATCAGCAGGACGCCAGCCTATACAGACACCAGCAGGGGTTTTAATTGGAGTTGGGAATATACCTTTTTTAACCCAGCTCCAGATCGTACGACGATCTCGGTTAAGCATCACAGTAAGGTCTGCATAGGTTAGTAGGCGGGTTGTTGTTGGTTGTAGCTTTGATAGTTGGTTGTTTTGTGACATATCTTCTTTCTCCTATTCAGTGAAGATACGTCATTATTAGTGATTTAATTTTTTATAGCATCGGAAAGAAATTTAAAATTTCCGTAGACTAATCAAAAGGGTTTTTGTATGAAAATACTGTTTTCTTTACCATCTTGTGCTTTGATTTATAGACAGAACCCCAAAACCTAGCCCTTTCTGCCTCATCAGTAATAGCGGCCGCCTCAAGCTCGATTTTTTCCCAGAGCTCAAGCTGTTTTCTAACTCTGTCATTTGTTATACCGTCCATTTGGGCGATGTGTTCAATTGCAGTTTGCTTCCTATCTGGACTCCAATACTTCCAATAGTGGTAAGCCCGAGCCCCCCAAAGCAGTGGTGAGTCTTTATCATTAGTTCGATATTTCATTTCATGCTTATGTATTTCAAAACAAAGTCTTTGCTCGTACCCAGTTGCAACGCCACTTTTAAAACGCCTCCATGCCGCCACTAAATAAGTATGTTTACGGTGATCTTGAGATAATAAAAACTTTGCATCTTCATCGCGATCGTAATCTTCACTTTGCTGCCGTTCAAACCACTTCCTATGCTCGGATATTATTCTTTGCACTTCATGCAGGCTGATGTCTAACTCAGTAACAGTTTGTGTTATATCTTCCTTGTTCAATTCACTCATCTAATCGCCCTCACATTTGCACCTGTTGAGATGTCGTTTAATTTTGATACCCATATATCAAGCGCTGCTTTCTTCTTTTCAATATGCTGCGAACGGTTATAAATCTTCATTACACCGCCCAGCGAATGCCCCAGCAATTGCTCGACCACGTACGGCTCTATTTGCATGTCGTTGAGAAGTGTTGCAAACGTGCGCCTTAAATCATGAAATGTCCACGCCTCTTGATGGCCTAGCTTTTGCCATAACTTGCCACCGAACCCAGATACCGCTTCTGGTCGTTTCAGTTCCCCCAACAGGAATTCAGAATTAGCTTGAGCAATGACCGATTGAAGATAGATTCTCATCTGCTCAGGTATTGGACGTTTTATTTCTGAACCTGTCTTGCTATTCGACTTGGGTACTGTCCAGATCATATCGTCCATATTCCATTCATCAACTTTTGATAGCCTTAACTCTTGAGTGCGACAACCAAACACTACCAGCAGCAAAGCAAGATTTCGGTAATACCAATTACTTTTTATGTCATTAGCCCAATGCCAAACATCCATAAACTCAGGCCATGATAAAACCCGTTCTTTCTTTCCTTGGGGCTCTCCTACATCAGCAATGTTTAGATCTTCTAATGCTACTGAGAGTGCATACTGACGATTACGGCAAAATTTTAATGCTTGTTTAGCGTTCTGTAGGATGTAGCCGCTTGCGACTGGTGCGGCTCTGTGGTGTGTGCCGCTGGTAACATCATCAAACACCTGCACCCAATGGCGTGTTTCACACTGCTCAAGGGGTAAGTGGCCTATAAACGGGTAAACATGCTTTTGAAATTGCGCCCTGTGCTTTTCGTGGTTAGTTCGTTTCTTTGCAGCAAAGTTATCTAGCCAATATTCGAGCGCATCCTTTACCGTTACCGCTTTGAGTTTCTGGTCTCGCTGTAGCTGTCGCTCAACCTTAGGATTCAAACCATCAGCAAGGTATTGCTTACACTCTGCCGCCTGCTTTCTTGCTTTAGCAAGGGTTAGCCCATCAGAACCGTCCTTTGAGTAATTTCCCAGCGTCAATTCTTGCGCTTTGCCTGCGTAGCGATAACGAAACACAAAAGATAGCTTTCCAGTCTTAAAATGGTAAACGCTCAACCCCTCTCTATCGGAAATTTTTAACGGCTTGTTGCTGTCGTGAGTCTTACCCTGCAACGCCTTTAGTTTTGAATCAGTGATCGCCATAGTTACCCGCCCTGCTAATGTCCTAAGTTGTCCTAAGTGATTTTTTAACTTAGGCCAGCTTGAGTCTTTGCTGTATATGGTTTTGTCTTAAGTGGCTAAGTGTCCTAAGTAAATTCTTAGTTTATTTTATTTTCCAAACCAAGCATTTTAACCCCCACATTAACCCCCACACTTCTCCATCTACCTACCAAGATAGCGCGAGACAAAGCAAGACACAACAATACATTAACCCACTGTAAATAAAGGCATTAGGTAAACTAATTAAGATTTATTGAGACGCAAGAAAATAGTAAAAATGAT"
dicoFASTA = make_fasta_dict(pathREADS)
NEWFASTA = open("/tmp/temp.fasta", 'w')
for line in lstLines:
    query = line.split(";")[0]
    sstart = int(line.split(";")[1])
    send = int(line.split(";")[2])
    if sstart>send: NEWFASTA.write(">"+query+"\n"+reverse_complement(dicoFASTA[query])+"\n")
    else: NEWFASTA.write(">"+query+"\n"+dicoFASTA[query]+"\n")
NEWFASTA.close()

os.system("blastn -task blastn -query /tmp/temp.fasta -subject "+pathREF+" -max_hsps 25 -outfmt \"6 delim=; qseqid qstart qend sstart send\" -evalue 0.001 1>/tmp/blast.out 2>/dev/null")
IN = open("/tmp/blast.out",'r')
lstLines = IN.read().split("\n")[:-1]
IN.close()
dicoResults = {}
for line in lstLines:
    query = line.split(";")[0]
    qstart = int(line.split(";")[1])
    qend = int(line.split(";")[2])
    sstart = int(line.split(";")[3])
    send = int(line.split(";")[4])
    if not query in dicoResults: dicoResults[query] = {'PICMIstart':9999, 'PICMIend':0, 'readStart':9999999, 'readEnd':0}
    if qend>dicoResults[query]['readEnd']:
        dicoResults[query]['readEnd'] = qend
        dicoResults[query]['PICMIend'] = send
    if qstart<dicoResults[query]['readStart']:
        dicoResults[query]['readStart'] = qstart
        dicoResults[query]['PICMIstart'] = sstart

dicoColor = {"barcode01":"#006400","barcode02":"#b03060","barcode03":"#ff4500","barcode04":"#ffff00","barcode05":"#deb887","barcode06":"#00ff00","barcode08":"#00ffff","barcode10":"#ff00ff","barcode12":"#6495ed"}

for barcode in dicoColor:
    features = []
    for query in dicoResults:
        if query.split("_____")[1]==barcode:
            color = dicoColor[query.split("_____")[1]]
            completePICMIsize = len(dicoFASTA[query])-(6110-dicoResults[query]['PICMIstart'])-(dicoResults[query]['PICMIend'])
            completePICMIcount = int(round((completePICMIsize/6110),1))
            seqModifReadAlign = seqPICMI[dicoResults[query]['PICMIstart']-1:]+completePICMIcount*seqPICMI+seqPICMI[:dicoResults[query]['PICMIend']]
            # seqModifReadAlign = "-"*(dicoResults[query]['PICMIstart']-1)+seqPICMI[dicoResults[query]['PICMIstart']-1:]+completePICMIcount*seqPICMI+seqPICMI[:dicoResults[query]['PICMIend']]
            # seqModifReadAlign += "-"*(60947-len(seqModifReadAlign))
            # print(">"+query+"\n"+seqModifReadAlign)
            features.append(GraphicFeature(start=dicoResults[query]['PICMIstart'], end=dicoResults[query]['PICMIstart']+len(seqModifReadAlign), strand=1, color="blue", linewidth=0))
    record = GraphicRecord(sequence_length=80000, features=features)
    ax, _ = record.plot(figure_width=50)
    # ax.figure.savefig("/tmp/temp.png", dpi=300)
    ax.figure.savefig("/tmp/"+barcode+".svg")
    plt.close('all')

exit()


print("Strain;nbCtgChr1;nbCtgChr2;nbCtgOthers")
for name in os.listdir("/mnt/c/Users/dgoudenege/Desktop/31_O_70_ragtag"):
    if os.path.isdir("/mnt/c/Users/dgoudenege/Desktop/31_O_70_ragtag/"+name):
        dicoFASTA = make_fasta_dict("/mnt/c/Users/dgoudenege/Desktop/31_O_70_ragtag/"+name+"/ragtag.scaffold.fasta")
        nbCtgChr1 = 0
        nbCtgChr2 = 0
        nbCtgOthers = 0
        for key in dicoFASTA:
            splitSeq = dicoFASTA[key].split("N")
            for e in splitSeq:
                if len(e)!=0:
                    if "chr1" in key: nbCtgChr1+=1
                    elif "chr2" in key: nbCtgChr2+=1
                    else: nbCtgOthers+=1
        print(name+";"+str(nbCtgChr1)+";"+str(nbCtgChr2)+";"+str(nbCtgOthers))
exit()



pathREF = "/mnt/c/Users/dgoudenege/Desktop/Clade1/Vibrio_crassostreae_7D8_10_complete.fasta"
# pathCTG = "/mnt/c/Users/dgoudenege/Desktop/Clade1/FNA_others_clade1/Vibrio_crassostreae_31_O_70.fna"
# pathCTG = "/tmp/workflow/abyss-pe_scaffolds.fasta"
pathCTG = "/tmp/SPADES_abyssCtg_nocovcutoff/contigs.fasta"
pathFASTQ1 = "/mnt/c/Users/dgoudenege/Desktop/Clade1/FASTQ_others_clade1/Vibrio_crassostreae_31_O_70_strain_3960_R1.fastq.gz"
pathFASTQ2 = "/mnt/c/Users/dgoudenege/Desktop/Clade1/FASTQ_others_clade1/Vibrio_crassostreae_31_O_70_strain_3960_R2.fastq.gz"
pathTMPOUT = "/tmp/temp.out"
pathTMPFASTA1 = "/tmp/temp1.fasta"
pathTMPFASTA2 = "/tmp/temp2.fasta"
# First BlastN to orient contigs
cmdBLASTN = "blastn -task blastn -query "+pathCTG+" -subject "+pathREF+" -out "+pathTMPOUT+" -max_hsps 5 -max_target_seqs 1 -qcov_hsp_perc 90 -outfmt \"6 delim=; qseqid sseqid pident length qstart qend qlen sstart send slen evalue\" -evalue 0.001 2>/dev/null"
os.system(cmdBLASTN)
lstLines = read_file(file_path=pathTMPOUT)
dicoCTG = make_fasta_dict(pathCTG, onlyLTasHeader=True)
dicoCTGCHR = {}
for ctg in dicoCTG:
    dicoCTGCHR[ctg] = set()
for line in lstLines:
    splitLine = line.split(";")
    qseqid = splitLine[0]
    sseqid = splitLine[1]
    sstart = int(splitLine[7])
    send = int(splitLine[8])
    dicoCTGCHR[qseqid].add(sseqid)

for ctg in dicoCTGCHR:
    if len(dicoCTG[ctg])>=200:
        if len(dicoCTGCHR[ctg])==0: dicoCTGCHR[ctg] = set(["chr1","chr2"])
        for chrom in dicoCTGCHR[ctg]:
            TMPFASTA = open("/tmp/SPADES_abyssCtg_nocovcutoff/"+chrom+".fasta", 'a')
            TMPFASTA.write(">"+ctg+"\n"+dicoCTG[ctg]+"\n")
        TMPFASTA.close()

exit()   
# Second BlastN to place contigs
os.remove(pathTMPOUT)
cmdBLASTN = "blastn -task blastn -query "+pathTMPFASTA+" -subject "+pathREF+" -out "+pathTMPOUT+" -max_hsps 5 -max_target_seqs 1 -qcov_hsp_perc 90 -outfmt \"6 delim=; qseqid sseqid pident length qstart qend qlen sstart send slen evalue\" -evalue 0.001 2>/dev/null"
os.system(cmdBLASTN)
lstLines = read_file(file_path=pathTMPOUT)
for line in lstLines:
    splitLine = line.split(";")
    # qseqid = splitLine[0]
    # sstart = int(splitLine[7])
    # send = int(splitLine[8])
    print(line)

# dicoREF = make_fasta_dict(pathREF)
# for refKey in dicoREF:
#     refSeq = dicoREF[refKey]








exit()




dicoALIGN = make_fasta_dict("/mnt/c/Users/dgoudenege/Desktop/NANOFLYE/P115_like_muscle.fasta")
size = len(list(dicoALIGN.values())[0])
cpt = 0
for i in range(size):
    setP115 = set()
    setOthers = set()
    for key in dicoALIGN:
        if key=="P115_illumina":
            ref = dicoALIGN[key][i]
            if ref!="-": cpt+=1
        elif key in ["P115","P115blue","deltaHH"]: setP115.add(dicoALIGN[key][i])
        else:
            setOthers.add(dicoALIGN[key][i])
    if len(setP115)==1 and list(setP115)[0]==ref and ((len(setOthers)==1 and list(setOthers)[0]!=ref) or len(setOthers)!=1):
        print(str(cpt)+";"+str(i)+";"+",".join(setOthers))
exit()


IN = open("/tmp/temp.out",'r')
lstLines = IN.read().split("\n")[:-1]
IN.close()
dicoRead = {}
for line in lstLines:
    splitLine = line.split(";")
    query = splitLine[0]
    subject = splitLine[1]
    qstart = int(splitLine[2])
    qend = int(splitLine[3])
    sstart = int(splitLine[4])
    send = int(splitLine[5])
    pident = float(splitLine[6])
    picmiSize = send-sstart
exit()




dicoFAA = {}
for file in os.listdir("/mnt/c/Users/dgoudenege/Desktop/NANOFLYE"):
    if "FLYE" in file and ".faa" in file:
        dico = make_fasta_dict("/mnt/c/Users/dgoudenege/Desktop/NANOFLYE/"+file)
        for key in dico:
            dicoFAA[key] =dico[key]
dicoCluster = load_json("/mnt/c/Users/dgoudenege/Desktop/NANOFLYE/prot_rbh.json")
for clusterNum in dicoCluster:
    setP115 = set()
    setOthers = set()
    for prot in dicoCluster[clusterNum]:
        bc = prot.split("[")[1].replace("FLYE_barcode","").replace("_PHAGE]","")
        if bc in ["10","11","12","7","9"]:
            setP115.add(bc)
        else:
            setOthers.add(bc)
    # print(clusterNum,len(setP115),len(setOthers))
    if len(setOthers)==7 and len(setP115)==0:
        print("\n")
        for prot in dicoCluster[clusterNum]:
            print(">"+clusterNum+"_____"+prot+"\n"+dicoFAA[prot])
        # print("Othersspe",clusterNum)
    if len(setOthers)==0 and len(setP115)==5:
        print("P115spe",clusterNum)


# IN = open("/tmp/temp.out", 'r')
# lstLines = IN.read().split("\n")[:-1]
# IN.close()
# dicoRead = {}
# for line in lstLines:
#     splitLine = line.split(";")
#     query = splitLine[0]
#     subject = splitLine[1]
#     qstart = int(splitLine[2])
#     qend = int(splitLine[3])
#     sstart = int(splitLine[4])
#     send = int(splitLine[5])
#     pident = float(splitLine[6])
#     picmiSize = send-sstart
#     if picmiSize>=1000:
#         if not query in dicoRead: dicoRead[query] = {}
#         dicoRead[query][qstart] = sstart
#         dicoRead[query][qend] = send

# for read in dicoRead:
#     minStart = min(list(dicoRead[read].keys()))
#     maxEnd = max(list(dicoRead[read].keys()))
#     print(str(dicoRead[read][minStart])+"\t"+str(dicoRead[read][maxEnd]))
# exit()



for pathFASTQ in glob.glob("/mnt/c/Users/dgoudenege/Desktop/NANOPORE/RUN/guppy_pass/barcode*/*"):
    reads = fq.read(pathFASTQ) # Is able to handle compressed files.
    bc = pathFASTQ.split("/")[9]
    OUT = open("/tmp/temp.fasta", 'w')
    dicoSeq = {}
    for read in reads:
        nameR = read.getHead().split(" ")[0].replace("@","")
        seqR = read.getSeq()
        if len(seqR)>=1000:
            OUT.write(">"+nameR+"\n"+seqR+"\n")
            dicoSeq[nameR] = seqR
    OUT.close()
    os.system("blastn -task blastn -num_threads 12 -query /tmp/integrase.ffn -subject /tmp/temp.fasta -max_hsps 10 -max_target_seqs 10000000 -outfmt \"6 delim=; qseqid sseqid qstart qend sstart send pident\" -evalue 0.000001 > /tmp/temp.out 2>/dev/null")
    IN = open("/tmp/temp.out", 'r')
    lstLines = IN.read().split("\n")[:-1]
    IN.close()
    dicoRead = {}
    for line in lstLines:
        splitLine = line.split(";")
        query = splitLine[0]
        subject = splitLine[1]
        qstart = int(splitLine[2])
        qend = int(splitLine[3])
        sstart = int(splitLine[4])
        send = int(splitLine[5])
        pident = float(splitLine[6])
        try: dicoRead[subject]+=1
        except: dicoRead[subject]=1
    maxi = 0
    for read in dicoRead:
        maxi = max(dicoRead[read],maxi)
    for i in range(1,maxi+1,1):
        cpt = 0
        for read in dicoRead:
            if dicoRead[read]==i: cpt+=1
        print(bc+";"+str(i)+";"+str(cpt))
        # if sstart>send: print(">"+read+"_____"+bc+"\n"+reverse_complement(dicoSeq[read]))
        # else: print(">"+read+"_____"+bc+"\n"+dicoSeq[read])
exit()


# pathDir = "/mnt/c/Users/dgoudenege/Desktop/NANOPORE/ASSEMBLY/PICMI"
# pathFASTA1 = "/mnt/c/Users/dgoudenege/Desktop/NANOPORE/PICMI_Vibrio_chagasii_34_P_115.fna"
# lstFILE = os.listdir(pathDir)
# pbar = tqdm(total=len(lstFILE),ncols=75,leave=False,desc="",file=sys.stdout,bar_format="  {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}")
# results = ""
# for file2 in lstFILE:
#     pathFASTA2 = pathDir+"/"+file2
#     os.system("cat "+pathFASTA1+" "+pathFASTA2+" > /tmp/temp.fasta")
#     os.system("famsa -t 12 /tmp/temp.fasta /tmp/temp.align > /dev/null 2>&1")
#     dicoALIGN = make_fasta_dict("/tmp/temp.align")
#     seq1 = list(dicoALIGN.values())[0]
#     seq2 = list(dicoALIGN.values())[1]
#     nbMatch = 0
#     for i in range(len(seq1)):
#         if seq1[i] == seq2[i]:
#             nbMatch += 1
#     pident = (nbMatch*100)/len(seq1)
#     results += file2.replace(".fasta","")+";"+str(pident)+"\n"
#     pbar.update(1)
# pbar.close()
# print(results)
# exit()



# pathDir = "/mnt/c/Users/dgoudenege/Desktop/NANOPORE/ASSEMBLY/PHAGES"
# pathRef = "/mnt/c/Users/dgoudenege/Desktop/NANOPORE/Vibrio_phage_115E34-1.gbk"
# dicoSNP = {}
# for file2 in os.listdir(pathDir):
#     orgName2 = os.path.basename(file2).replace(".fasta","")
#     cmd = "/mnt/c/Users/dgoudenege/Tools/snippy/bin/snippy --force --outdir /tmp/snippy --ref "+pathRef+" --ctgs "+pathDir+"/"+file2+" --quiet > /dev/null 2>&1"
#     os.system(cmd)
#     IN = open("/tmp/snippy/snps.tab",'r')
#     lstLines = IN.read().split("\n")[1:-1]
#     IN.close()
#     dicoSNP[orgName2] = len(lstLines)
#     for line in lstLines:
#         splitLine = line.split("\t")
#         varID = splitLine[1]+"_"+splitLine[3]+"_"+splitLine[4]
#         print(orgName2+";"+varID)

# print("\n\n")
# for org2 in dicoSNP:
#     print(org2+";"+str(dicoSNP[org2]))
# exit()

# dicoAlign = {}
# for alignFile in os.listdir("/tmp/genesTree_25_50/core_align"):
#     if "all" in alignFile: continue
#     dicoFASTA = make_fasta_dict("/tmp/genesTree_25_50/core_align/"+alignFile)
#     print(alignFile)
#     for key in dicoFASTA:
#         orgName = key.split("[")[1].replace("]", "")
#         print(orgName)
#         try:
#             dicoAlign[orgName] += dicoFASTA[key]
#         except KeyError:
#             dicoAlign[orgName] = dicoFASTA[key]
#     for k in dicoAlign: print("  ",len(dicoAlign[k]))
# ALL = open("/tmp/genesTree_25_50/core_align/all_core_align.fasta", 'w')
# for orgName in dicoAlign:
#     ALL.write(">"+orgName+"\n"+dicoAlign[orgName]+"\n")
# ALL.close()

# # for file in os.listdir("/tmp/genesTree_25_50/core_align"):
# #     dico = make_fasta_dict("/tmp/genesTree_25_50/core_align/"+file)
# #     setSize = set()
# #     for key in dico: setSize.add(len(dico[key]))
# #     print(file,len(set(dico.keys())),len(setSize))

# exit()


# dico = make_fasta_dict("/mnt/c/Users/dgoudenege/Desktop/NANOPORE/NECAT/P115_like.fasta")
# for key in dico:
#     if key=="NANOPORE_F12": bc = "barcode01"
#     elif key=="NANOPORE_F15": bc = "barcode02"
#     elif key=="NANOPORE_F17": bc = "barcode03"
#     elif key=="NANOPORE_F18": bc = "barcode04"
#     elif key=="NANOPORE_F20": bc = "barcode05"
#     elif key=="NANOPORE_F21": bc = "barcode06"
#     elif key=="NANOPORE_F27": bc = "barcode07"
#     elif key=="NANOPORE_F35": bc = "barcode08"
#     elif key=="NANOPORE_P27": bc = "barcode09"
#     elif key=="NANOPORE_P115": bc = "barcode11"
#     else: continue
#     OUT = open("/mnt/c/Users/dgoudenege/Desktop/NANOPORE/NECAT/"+bc+"/"+key+".fna",'w')
#     OUT.write(">"+key+"\n"+dico[key]+"\n")
#     OUT.close()
# exit()

# 334;6110;42;5773
# 1;6110;5757;11831
# 1;6110;11815;17917
# 1;6110;17901;23980
# 1;6110;23964;29957
# 7;6110;29947;36007
# 1;6110;35991;42024
# 1;6110;42008;48072
# 1;1690;48056;49732
pathPICMI = "/mnt/c/Users/dgoudenege/Desktop/NANOPORE/PICMI_Vibrio_chagasii_34_P_115.fna"




# pathPICMI = "/mnt/c/Users/dgoudenege/Desktop/NANOPORE/PICMI_Vibrio_chagasii_34_P_115.fna"
# lstFiles = glob.glob("/mnt/c/Users/dgoudenege/Desktop/NANOPORE/RUN/guppy_pass/barcode*/*")
# for pathFASTQ in lstFiles:
#     bc = pathFASTQ.split("/")[9]
#     reads = fq.read(pathFASTQ)
#     OUT = open("/tmp/bc_fasta/"+bc+".fasta", 'a')
#     for read in reads:
#         if len(read.getSeq())>=1000:
#             OUT.write(">"+read.getHead().split(" ")[0].replace("@","")+"_____"+bc+"\n"+read.getSeq()+"\n")
#     OUT.close()
# exit()

pathPICMI = "/mnt/c/Users/dgoudenege/Desktop/NANOPORE/PICMI_Vibrio_chagasii_34_P_115.fna"
lstFiles = glob.glob("/mnt/c/Users/dgoudenege/Desktop/NANOPORE/RUN/fastq_pass/barcode*/*")
for pathFASTQ in lstFiles:
    bc = pathFASTQ.split("/")[9]
    reads = fq.read(pathFASTQ)
    OUT = open("/tmp/temp.fasta", 'w')
    for read in reads:
        nameR = read.getHead().split(" ")[0].replace("@","")
        seqR = read.getSeq()
        if len(seqR)>=1000:
            OUT.write(">"+nameR+"\n"+seqR+"\n")
    OUT.close()
    os.system("blastn -task blastn -qcov_hsp_perc 75 -query "+pathPICMI+" -subject /tmp/temp.fasta -max_hsps 25 -max_target_seqs 1000000000 -outfmt \"6 delim=; qseqid sseqid qstart qend sstart send pident seq\" -evalue 0.000001 > /tmp/temp.out 2>/dev/null")
    dicoFASTA = make_fasta_dict("/tmp/temp.fasta")
    IN = open("/tmp/temp.out", 'r')
    lstLines = IN.read().split("\n")[:-1]
    IN.close()
    dicoRead = {}
    for line in lstLines:
        splitLine = line.split(";")
        query = splitLine[0]
        subject = splitLine[1]
        qstart = int(splitLine[2])
        qend = int(splitLine[3])
        sstart = int(splitLine[4])
        send = int(splitLine[5])
        pident = float(splitLine[6])
        if pident>=90:
            # print(">"+subject+"\n"+dicoFASTA[subject])
            try: dicoRead[subject] += qend-qstart
            except: dicoRead[subject] = qend-qstart
    dicoCount = {}
    for subject in dicoRead:
        count = round(dicoRead[subject]/6093)
        if count != 0:
            try: dicoCount[count] += 1
            except: dicoCount[count] = 1
    print("\n"+bc)
    for count in sorted(dicoCount):
        print("    "+str(count)+" = "+str(dicoCount[count]))
# AATACGGCATGAACTAA
# exit()

# dicoAlpA = load_json("/mnt/c/Users/dgoudenege/Works/AlpA/AlpA_fromHMM_e10-6.json")
# lstFiles = glob.glob(pathDIRHMM+"/*/*")
# cpt = 0
# for pathHMM in lstFiles:
#     subDir = os.path.basename(os.path.dirname(pathHMM))
#     TMP = open(pathHMM, 'r')
#     lstLines = TMP.read().split("\n")
#     TMP.close()
#     for line in lstLines:
#         if line != "" and line[0] != "#":
#             splitLine = re.compile(r"\s+").sub(" ", line).strip().split(" ")
#             if subDir == "genbank_243": target = splitLine[0]
#             else: continue
#             fullSeqEvalue = float(splitLine[4])
#             if fullSeqEvalue <= 0.000001:
#                 cpt+=1
#                 if not target in dicoAlpA:
#                     print(target)
# print(cpt)
# exit()



exit()      


# dicoFAA = make_fasta_dict("/mnt/c/Users/dgoudenege/Works/AlpA/all_vibrionaceae.faa")
# dicoALPA = {}
# for key in dicoFAA:
#     seq = dicoFAA[key]
#     if not seq in dicoALPA and "partial" not in key:
#         lt = key.split("|")[0]
#         profile = key.split("|")[1]
#         evalue = float(key.split("|")[2])
#         descr = key.split("|")[3].split("[")[0].replace("TPA: ","").replace("MULTISPECIES: ","").replace("MAG: ","")
#         if "hypothetical protein" in descr: descr = "hypothetical"
#         org = key.split("|")[3].split("[")[1].replace("]]","")
#         if evalue<=0.000001:
#             dicoALPA[seq] = {'picmi':False, 'descr':descr.replace(" ","_")}

# for org in os.listdir("/mnt/c/Users/dgoudenege/Desktop/PICMI/FAA"):
#     dicoTMP = make_fasta_dict("/mnt/c/Users/dgoudenege/Desktop/PICMI/FAA/"+org)
#     for key in dicoTMP:
#         if dicoTMP[key] in dicoALPA:
#             dicoALPA[dicoTMP[key]]['picmi'] = True

# OUT = open("/mnt/c/Users/dgoudenege/Works/AlpA/alpa.fasta",'w')
# cpt = 1
# for seq in dicoALPA:
#     if dicoALPA[seq]['picmi'] is True: OUT.write(">"+str(cpt).zfill(4)+"_PICMI\n"+seq+"\n")
#     else: OUT.write(">"+str(cpt).zfill(4)+"\n"+seq+"\n")
#     cpt += 1
# OUT.close()

# exit()


dicoChaga = {}
dico = make_fasta_dict("/tmp/34_P_115.ffn")
for key in dico:
    try: dicoChaga[dico[key][:-3]].append(key)
    except: dicoChaga[dico[key][:-3]] = [key]
dicoNum = {}
for file in os.listdir("/tmp/pana"):
    if ".aln" in file:
        num = file.split(".")[1]
        dicoFFN = make_fasta_dict("/tmp/pana/"+file)
        seqChaga = ""
        for key in dicoFFN:
            if "34_P_115" in key:
                seqChaga = dicoFFN[key].replace("-","")
                break
        if seqChaga != "":
            for seq in dicoChaga:
                if seqChaga in seq or seq in seqChaga:
                    dicoNum[num] = " ".join(dicoChaga[seq])

lst = ["68","68","68","68","210","210","210","210","250","250","250","250","354","354","354","354","361","361","361","361","364","364","364","364","367","367","367","367","373","373","373","373","375","375","375","375","377","377","377","377","385","385","385","385","399","399","399","399","400","400","400","400","414","414","414","414","439","439","439","439","458","458","458","458","467","467","467","467","473","473","473","473","477","477","477","477","535","535","535","535","606","606","606","606","741","741","741","741","776","776","776","776","836","836","836","836","998","998","998","998","1008","1008","1008","1008","1045","1045","1045","1045","1478","1478","1478","1478","1480","1480","1480","1480","1495","1495","1495","1495","1499","1499","1499","1499","1502","1502","1502","1502","1524","1524","1524","1524","1528","1528","1528","1528","1567","1567","1567","1567","1648","1648","1648","1648","1685","1685","1685","1685","1779","1779","1779","1779","1786","1786","1786","1786","1833","1833","1833","1833","1835","1835","1835","1835","2065","2065","2065","2065","2149","2149","2149","2149","2186","2186","2186","2186","2189","2189","2189","2189","2334","2334","2334","2334","2367","2367","2367","2367","2394","2394","2394","2394","2560","2560","2560","2560","2721","2721","2721","2721","2817","2817","2817","2817","3090","3090","3090","3090","3092","3092","3092","3092","3100","3100","3100","3100","3172","3172","3172","3172","3216","3216","3216","3216","3241","3241","3241","3241","3251","3251","3251","3251","3266","3266","3266","3266","3290","3290","3290","3290","3297","3297","3297","3297","3302","3302","3302","3302","3304","3304","3304","3304","3472","3472","3472","3472","3473","3473","3473","3473","3480","3480","3480","3480","3603","3603","3603","3603","3740","3740","3740","3740","3749","3749","3749","3749","3765","3765","3765","3765","3767","3767","3767","3767","3924","3924","3924","3924","3925","3925","3925","3925","4000","4000","4000","4000","4050","4050","4050","4050","4089","4089","4089","4089","4095","4095","4095","4095","4123","4123","4123","4123","4222","4222","4222","4222","4232","4232","4232","4232","4256","4256","4256","4256","4362","4362","4362","4362","4372","4372","4372","4372","4376","4376","4376","4376","4390","4390","4390","4390","4774","4774","4774","4774","4801","4801","4801","4801","4809","4809","4809","4809","4967","4967","4967","4967","4968","4968","4968","4968","4975","4975","4975","4975","5043","5043","5043","5043","5064","5064","5064","5064","5117","5117","5117","5117","5205","5205","5205","5205","5239","5239","5239","5239","5366","5366","5366","5366","5401","5401","5401","5401","5654","5654","5654","5654","5671","5671","5671","5671","5702","5702","5702","5702","5705","5705","5705","5705","5747","5747","5747","5747","5748","5748","5748","5748","5759","5759","5759","5759","5760","5760","5760","5760","5763","5763","5763","5763","5792","5792","5792","5792","5836","5836","5836","5836","5840","5840","5840","5840","5850","5850","5850","5850","5854","5854","5854","5854","5893","5893","5893","5893","5900","5900","5900","5900","5930","5930","5930","5930","5935","5935","5935","5935","5947","5947","5947","5947","5965","5965","5965","5965","6234","6234","6234","6234","6282","6282","6282","6282","6284","6284","6284","6284","6298","6298","6298","6298","6323","6323","6323","6323","6331","6331","6331","6331","6450","6450","6450","6450","6474","6474","6474","6474","6484","6484","6484","6484","6491","6491","6491","6491","6497","6497","6497","6497","6503","6503","6503","6503","6519","6519","6519","6519","6538","6538","6538","6538","6584","6584","6584","6584","6632","6632","6632","6632","6636","6636","6636","6636","6653","6653","6653","6653","6707","6707","6707","6707","6729","6729","6729","6729","6778","6778","6778","6778","6854","6854","6854","6854","6889","6889","6889","6889","6979","6979","6979","6979","7013","7013","7013","7013","7030","7030","7030","7030","7042","7042","7042","7042","7053","7053","7053","7053","7092","7092","7092","7092","7131","7131","7131","7131","7145","7145","7145","7145","7209","7209","7209","7209","7213","7213","7213","7213","7232","7232","7232","7232","7273","7273","7273","7273","7279","7279","7279","7279","7283","7283","7283","7283","7388","7388","7388","7388","7398","7398","7398","7398","7465","7465","7465","7465","7707","7707","7707","7707","7713","7713","7713","7713","7785","7785","7785","7785","7800","7800","7800","7800","7858","7858","7858","7858","7870","7870","7870","7870","7945","7945","7945","7945","8051","8051","8051","8051","8255","8255","8255","8255","8335","8335","8335","8335","8349","8349","8349","8349","8400","8400","8400","8400","8433","8433","8433","8433","8440","8440","8440","8440","8460","8460","8460","8460","8473","8473","8473","8473","8491","8491","8491","8491","8495","8495","8495","8495","8642","8642","8642","8642","8654","8654","8654","8654","8708","8708","8708","8708","8734","8734","8734","8734","8770","8770","8770","8770","8784","8784","8784","8784","8822","8822","8822","8822","8826","8826","8826","8826","8873","8873","8873","8873","8982","8982","8982","8982","9104","9104","9104","9104","9143","9143","9143","9143","9150","9150","9150","9150","9318","9318","9318","9318","9663","9663","9663","9663","9664","9664","9664","9664","9697","9697","9697","9697","9823","9823","9823","9823","9876","9876","9876","9876","9959","9959","9959","9959","9978","9978","9978","9978","10030","10030","10030","10030","10032","10032","10032","10032","10199","10199","10199","10199","10241","10241","10241","10241","10425","10425","10425","10425","10451","10451","10451","10451","10533","10533","10533","10533","10537","10537","10537","10537","10543","10543","10543","10543","10556","10556","10556","10556","10562","10562","10562","10562","10579","10579","10579","10579","10613","10613","10613","10613","10617","10617","10617","10617","10751","10751","10751","10751","10890","10890","10890","10890","10904","10904","10904","10904","10913","10913","10913","10913","10939","10939","10939","10939","11017","11017","11017","11017","11027","11027","11027","11027","11036","11036","11036","11036","11084","11084","11084","11084","11167","11167","11167","11167","11234","11234","11234","11234","11286","11286","11286","11286","11422","11422","11422","11422","11456","11456","11456","11456","11476","11476","11476","11476","11680","11680","11680","11680","11702","11702","11702","11702","11704","11704","11704","11704","11801","11801","11801","11801","11831","11831","11831","11831","11834","11834","11834","11834","11908","11908","11908","11908","11948","11948","11948","11948","12105","12105","12105","12105","12106","12106","12106","12106","12110","12110","12110","12110","12131","12131","12131","12131","12144","12144","12144","12144","12145","12145","12145","12145","12302","12302","12302","12302","12342","12342","12342","12342","12484","12484","12484","12484","12565","12565","12565","12565","12583","12583","12583","12583","12597","12597","12597","12597","12604","12604","12604","12604","12765","12765","12765","12765","12891","12891","12891","12891","12892","12892","12892","12892","12893","12893","12893","12893","12898","12898","12898","12898","12917","12917","12917","12917","12930","12930","12930","12930","12948","12948","12948","12948","12951","12951","12951","12951","12952","12952","12952","12952","12965","12965","12965","12965","12994","12994","12994","12994","13026","13026","13026","13026","13096","13096","13096","13096","13132","13132","13132","13132","13134","13134","13134","13134","13139","13139","13139","13139","13411","13411","13411","13411","13425","13425","13425","13425","13438","13438","13438","13438","13453","13453","13453","13453","13455","13455","13455","13455","13470","13470","13470","13470","13477","13477","13477","13477","13480","13480","13480","13480","13549","13549","13549","13549","13609","13609","13609","13609","13633","13633","13633","13633","13722","13722","13722","13722","13725","13725","13725","13725","13756","13756","13756","13756","13766","13766","13766","13766","13811","13811","13811","13811","13848","13848","13848","13848","13976","13976","13976","13976","14016","14016","14016","14016","14059","14059","14059","14059","14091","14091","14091","14091","14132","14132","14132","14132","14148","14148","14148","14148","14152","14152","14152","14152","14156","14156","14156","14156","14160","14160","14160","14160","14185","14185","14185","14185","14241","14241","14241","14241","14307","14307","14307","14307","14389","14389","14389","14389","14398","14398","14398","14398","14429","14429","14429","14429","14651","14651","14651","14651","14660","14660","14660","14660","14666","14666","14666","14666","14798","14798","14798","14798","14924","14924","14924","14924","14974","14974","14974","14974","15031","15031","15031","15031","15045","15045","15045","15045","15075","15075","15075","15075","15093","15093","15093","15093","15180","15180","15180","15180","15261","15261","15261","15261","15277","15277","15277","15277","15286","15286","15286","15286","15301","15301","15301","15301","15322","15322","15322","15322","15392","15392","15392","15392","15507","15507","15507","15507","15554","15554","15554","15554","15573","15573","15573","15573","15645","15645","15645","15645","15674","15674","15674","15674","15690","15690","15690","15690","15721","15721","15721","15721","15806","15806","15806","15806","15821","15821","15821","15821","15836","15836","15836","15836","15868","15868","15868","15868","15872","15872","15872","15872","15984","15984","15984","15984","16019","16019","16019","16019","16023","16023","16023","16023","16033","16033","16033","16033","16072","16072","16072","16072","16083","16083","16083","16083","16099","16099","16099","16099","16110","16110","16110","16110","16196","16196","16196","16196","16201","16201","16201","16201","16219","16219","16219","16219","16289","16289","16289","16289","16346","16346","16346","16346","16405","16405","16405","16405","16492","16492","16492","16492","16493","16493","16493","16493","16538","16538","16538","16538","16550","16550","16550","16550","16637","16637","16637","16637","16887","16887","16887","16887","16889","16889","16889","16889","16950","16950","16950","16950","17027","17027","17027","17027","17072","17072","17072","17072","17121","17121","17121","17121","17259","17259","17259","17259","17279","17279","17279","17279","17365","17365","17365","17365","17391","17391","17391","17391","17448","17448","17448","17448","17498","17498","17498","17498","17676","17676","17676","17676","17692","17692","17692","17692","17717","17717","17717","17717","17812","17812","17812","17812","17983","17983","17983","17983","17985","17985","17985","17985","17989","17989","17989","17989","17995","17995","17995","17995","18010","18010","18010","18010","18012","18012","18012","18012","18041","18041","18041","18041","18169","18169","18169","18169","18188","18188","18188","18188","18189","18189","18189","18189","18230","18230","18230","18230","18371","18371","18371","18371","18489","18489","18489","18489","18490","18490","18490","18490","18510","18510","18510","18510","18597","18597","18597","18597","18611","18611","18611","18611","18614","18614","18614","18614","18616","18616","18616","18616","18622","18622","18622","18622","18627","18627","18627","18627","18646","18646","18646","18646","18703","18703","18703","18703","18831","18831","18831","18831","18891","18891","18891","18891","19064","19064","19064","19064","19067","19067","19067","19067","19147","19147","19147","19147","19150","19150","19150","19150","19160","19160","19160","19160","19199","19199","19199","19199","19218","19218","19218","19218","19253","19253","19253","19253","19343","19343","19343","19343","19475","19475","19475","19475","19530","19530","19530","19530","19573","19573","19573","19573","19579","19579","19579","19579","19626","19626","19626","19626","19646","19646","19646","19646","19750","19750","19750","19750","19765","19765","19765","19765","19768","19768","19768","19768","19892","19892","19892","19892","19898","19898","19898","19898","19899","19899","19899","19899","19902","19902","19902","19902","20037","20037","20037","20037","20040","20040","20040","20040","20044","20044","20044","20044","20047","20047","20047","20047","20111","20111","20111","20111","20206","20206","20206","20206","20264","20264","20264","20264","20336","20336","20336","20336","20449","20449","20449","20449","20688","20688","20688","20688","21157","21157","21157","21157","21232","21232","21232","21232","21235","21235","21235","21235","21361","21361","21361","21361","21367","21367","21367","21367","21406","21406","21406","21406","21424","21424","21424","21424","21435","21435","21435","21435","21444","21444","21444","21444","21456","21456","21456","21456","21498","21498","21498","21498","21538","21538","21538","21538","21567","21567","21567","21567","21579","21579","21579","21579","21589","21589","21589","21589","21598","21598","21598","21598","21649","21649","21649","21649","21659","21659","21659","21659","21693","21693","21693","21693","21736","21736","21736","21736","22000","22000","22000","22000","22008","22008","22008","22008","22074","22074","22074","22074","22163","22163","22163","22163","22320","22320","22320","22320","22354","22354","22354","22354","22420","22420","22420","22420","22426","22426","22426","22426","22475","22475","22475","22475","22541","22541","22541","22541","22590","22590","22590","22590","22602","22602","22602","22602","22765","22765","22765","22765","22782","22782","22782","22782","22825","22825","22825","22825","22846","22846","22846","22846","22904","22904","22904","22904","22945","22945","22945","22945","22954","22954","22954","22954","22976","22976","22976","22976","23052","23052","23052","23052","23057","23057","23057","23057","23058","23058","23058","23058","23125","23125","23125","23125","23210","23210","23210","23210","23238","23238","23238","23238","23247","23247","23247","23247","23298","23298","23298","23298","23339","23339","23339","23339","23393","23393","23393","23393","23467","23467","23467","23467","23657","23657","23657","23657","23735","23735","23735","23735","23754","23754","23754","23754","23767","23767","23767","23767","23832","23832","23832","23832","23843","23843","23843","23843","23852","23852","23852","23852","23915","23915","23915","23915","23923","23923","23923","23923","24380","24380","24380","24380","24488","24488","24488","24488","24589","24589","24589","24589","24617","24617","24617","24617","24678","24678","24678","24678","24704","24704","24704","24704","24802","24802","24802","24802","24814","24814","24814","24814","24815","24815","24815","24815","24834","24834","24834","24834","24852","24852","24852","24852","24855","24855","24855","24855","24857","24857","24857","24857","24894","24894","24894","24894","24921","24921","24921","24921","24938","24938","24938","24938","24957","24957","24957","24957","24962","24962","24962","24962","25172","25172","25172","25172","25184","25184","25184","25184","25484","25484","25484","25484","25515","25515","25515","25515","25567","25567","25567","25567","25581","25581","25581","25581","25620","25620","25620","25620","25633","25633","25633","25633","25796","25796","25796","25796","25891","25891","25891","25891","25907","25907","25907","25907","25995","25995","25995","25995","25997","25997","25997","25997","26083","26083","26083","26083","26090","26090","26090","26090","26358","26358","26358","26358","26359","26359","26359","26359","26363","26363","26363","26363","26374","26374","26374","26374","26539","26539","26539","26539","26544","26544","26544","26544","26556","26556","26556","26556","26565","26565","26565","26565","26577","26577","26577","26577","26595","26595","26595","26595","26610","26610","26610","26610","26612","26612","26612","26612","26628","26628","26628","26628","26635","26635","26635","26635","26655","26655","26655","26655","26689","26689","26689","26689","26690","26690","26690","26690","26700","26700","26700","26700","26759","26759","26759","26759","26804","26804","26804","26804","26809","26809","26809","26809","26871","26871","26871","26871","26873","26873","26873","26873","26879","26879","26879","26879","26980","26980","26980","26980","27036","27036","27036","27036","27044","27044","27044","27044","27129","27129","27129","27129","27375","27375","27375","27375","27392","27392","27392","27392","27393","27393","27393","27393","27884","27884","27884","27884","27889","27889","27889","27889","27892","27892","27892","27892","27897","27897","27897","27897","27899","27899","27899","27899","27905","27905","27905","27905","27915","27915","27915","27915","27917","27917","27917","27917","27918","27918","27918","27918","27919","27919","27919","27919","27932","27932","27932","27932","28047","28047","28047","28047","28113","28113","28113","28113","28122","28122","28122","28122","28178","28178","28178","28178","28217","28217","28217","28217","28343","28343","28343","28343","28380","28380","28380","28380","28442","28442","28442","28442","28475","28475","28475","28475","28488","28488","28488","28488","28534","28534","28534","28534","28607","28607","28607","28607","28610","28610","28610","28610","28611","28611","28611","28611","28613","28613","28613","28613","28629","28629","28629","28629","28632","28632","28632","28632","28703","28703","28703","28703","28705","28705","28705","28705","28716","28716","28716","28716","28772","28772","28772","28772","28789","28789","28789","28789","28814","28814","28814","28814","28847","28847","28847","28847","28899","28899","28899","28899","28903","28903","28903","28903","28950","28950","28950","28950","28952","28952","28952","28952","29010","29010","29010","29010","29066","29066","29066","29066","29194","29194","29194","29194","29229","29229","29229","29229","29361","29361","29361","29361","29780","29780","29780","29780","29858","29858","29858","29858","29872","29872","29872","29872","29906","29906","29906","29906","29928","29928","29928","29928","30068","30068","30068","30068","30100","30100","30100","30100","30208","30208","30208","30208","30247","30247","30247","30247","30309","30309","30309","30309","30341","30341","30341","30341","30511","30511","30511","30511","30525","30525","30525","30525","30659","30659","30659","30659","30668","30668","30668","30668","30760","30760","30760","30760","30822","30822","30822","30822","30865","30865","30865","30865","30932","30932","30932","30932","30933","30933","30933","30933","30991","30991","30991","30991","31013","31013","31013","31013","31025","31025","31025","31025","31028","31028","31028","31028","31076","31076","31076","31076","31419","31419","31419","31419","31477","31477","31477","31477","31490","31490","31490","31490","31508","31508","31508","31508","31561","31561","31561","31561","31616","31616","31616","31616","31674","31674","31674","31674","31757","31757","31757","31757","31780","31780","31780","31780","31792","31792","31792","31792","31802","31802","31802","31802","31827","31827","31827","31827","31886","31886","31886","31886","32001","32001","32001","32001","32042","32042","32042","32042","32080","32080","32080","32080","32087","32087","32087","32087","32104","32104","32104","32104","32119","32119","32119","32119","32121","32121","32121","32121","32281","32281","32281","32281","32324","32324","32324","32324","32439","32439","32439","32439","32508","32508","32508","32508","32510","32510","32510","32510","32539","32539","32539","32539","32545","32545","32545","32545","32593","32593","32593","32593","32610","32610","32610","32610","32704","32704","32704","32704","32834","32834","32834","32834","32844","32844","32844","32844","33053","33053","33053","33053","33149","33149","33149","33149","33320","33320","33320","33320","33398","33398","33398","33398","33554","33554","33554","33554","33557","33557","33557","33557","33567","33567","33567","33567","33585","33585","33585","33585","33603","33603","33603","33603","33616","33616","33616","33616","33677","33677","33677","33677","33680","33680","33680","33680","33782","33782","33782","33782","33803","33803","33803","33803","33822","33822","33822","33822","33842","33842","33842","33842","33860","33860","33860","33860","34002","34002","34002","34002","34026","34026","34026","34026","34151","34151","34151","34151","34168","34168","34168","34168","34365","34365","34365","34365","34387","34387","34387","34387","34521","34521","34521","34521","34562","34562","34562","34562","34587","34587","34587","34587","34603","34603","34603","34603","34604","34604","34604","34604","34634","34634","34634","34634","34679","34679","34679","34679","34798","34798","34798","34798","34814","34814","34814","34814","34915","34915","34915","34915","34924","34924","34924","34924","34928","34928","34928","34928","34961","34961","34961","34961","34980","34980","34980","34980","34986","34986","34986","34986","35019","35019","35019","35019","35027","35027","35027","35027","35038","35038","35038","35038","35063","35063","35063","35063","35082","35082","35082","35082","35086","35086","35086","35086","35102","35102","35102","35102","35105","35105","35105","35105","35115","35115","35115","35115","35183","35183","35183","35183","35226","35226","35226","35226","35405","35405","35405","35405","35727","35727","35727","35727","35756","35756","35756","35756","35780","35780","35780","35780","35794","35794","35794","35794","35885","35885","35885","35885","35917","35917","35917","35917","35923","35923","35923","35923","36006","36006","36006","36006","36017","36017","36017","36017","36106","36106","36106","36106","36120","36120","36120","36120","36128","36128","36128","36128","36156","36156","36156","36156","36230","36230","36230","36230","36404","36404","36404","36404","36416","36416","36416","36416","36426","36426","36426","36426","36427","36427","36427","36427","36444","36444","36444","36444","36453","36453","36453","36453","36557","36557","36557","36557","36744","36744","36744","36744","36790","36790","36790","36790","36800","36800","36800","36800","36971","36971","36971","36971","36986","36986","36986","36986","37130","37130","37130","37130","37238","37238","37238","37238","37302","37302","37302","37302","37335","37335","37335","37335","37347","37347","37347","37347","37433","37433","37433","37433","37435","37435","37435","37435","37631","37631","37631","37631","37688","37688","37688","37688","37716","37716","37716","37716","37724","37724","37724","37724","37795","37795","37795","37795","37938","37938","37938","37938","38002","38002","38002","38002","38039","38039","38039","38039","38102","38102","38102","38102","38108","38108","38108","38108","38165","38165","38165","38165","38244","38244","38244","38244","38336","38336","38336","38336","38390","38390","38390","38390","38391","38391","38391","38391","38468","38468","38468","38468","38501","38501","38501","38501","38516","38516","38516","38516","38524","38524","38524","38524","38532","38532","38532","38532","38537","38537","38537","38537","38580","38580","38580","38580","38581","38581","38581","38581","38600","38600","38600","38600","38615","38615","38615","38615","38617","38617","38617","38617","38628","38628","38628","38628","38630","38630","38630","38630","38700","38700","38700","38700","38823","38823","38823","38823","38881","38881","38881","38881","38882","38882","38882","38882","38885","38885","38885","38885","38931","38931","38931","38931","39032","39032","39032","39032","39033","39033","39033","39033","39103","39103","39103","39103","39172","39172","39172","39172","39300","39300","39300","39300","39368","39368","39368","39368","39425","39425","39425","39425","39670","39670","39670","39670","39684","39684","39684","39684","39688","39688","39688","39688","39702","39702","39702","39702","39706","39706","39706","39706"]
for num in lst:
    try: print(num+";"+dicoNum[num])
    except: print(num)
exit()


# dico = {}
# for prot1 in ["Fis","Integrase","AlpA","Primase","Core"]:
#     pathTREE1 = "/mnt/c/Users/dgoudenege/Desktop/PICMI_tanglegram/by_gene/"+prot1+"_famsa_iqtree2.nwk"
#     dico[prot1] = {}
#     for prot2 in ["Fis","Integrase","AlpA","Primase","Core"]:
#         pathTREE2 = "/mnt/c/Users/dgoudenege/Desktop/PICMI_tanglegram/by_gene/"+prot1+"_famsa_iqtree2.nwk"
#         R = open("/tmp/tanglegram.R", 'w')
#         R.write("library(ape)\n")
#         R.write("library(phytools)\n")
#         R.write("library(dendextend)\n")
#         R.write("library(viridis)\n")
#         R.write("library(dplyr)\n")
#         R.write("library(phylogram)\n")
#         R.write("tree1 <- read.tree(file = \""+pathTREE1+"\")\n")
#         R.write("tree2 <- read.tree(file = \""+pathTREE2+"\")\n")
#         R.write("tree1 <- compute.brlen(tree1)\n")
#         R.write("tree2 <- compute.brlen(tree2)\n")
#         # R.write("tree1 <- midpoint.root(tree1)\n")
#         # R.write("tree2 <- midpoint.root(tree2)\n")
#         R.write("dend1<- as.dendrogram(tree1)\n")
#         R.write("dend2<- as.dendrogram(tree2)\n")
#         R.write("dend12_corrected_1 <- untangle_random_search(dend1, dend2)\n")
#         R.write("dend12_corrected_2 <- untangle_step_rotate_2side(dend12_corrected_1[[1]], dend12_corrected_1[[2]])\n")
#         R.write("entanglement(dend12_corrected_2[[1]], dend12_corrected_2[[2]], L = 2)\n")
#         R.close()
#         os.system("Rscript /tmp/tanglegram.R > /tmp/temp.out 2>/dev/null")
#         OUT = open("/tmp/temp.out",'r')
#         lstLines = OUT.read().split("\n")[:-1]
#         OUT.close()
#         for line in lstLines:
#             dico[prot1][prot2] = float(line.split(" ")[1])


# # Make matrix
# df = pd.DataFrame(dico)
# figsize = len(df)/ np.sqrt(len(df))
# cmap = sns.color_palette("light:#d40000", as_cmap=True)
# cg = sns.clustermap(df, cmap="flare", figsize=(figsize,figsize), tree_kws={'linewidths':1.5}, method="average", metric="correlation", dendrogram_ratio=0.05, xticklabels=True, yticklabels=True, linewidths=0.0, rasterized=True)
# cg.ax_cbar.remove()
# plt.tick_params(axis='both', which='major', labelsize=1.5, width=0.2)
# fig = plt.gcf()  # or by other means, like plt.subplots
# figsize = fig.get_size_inches()
# fig.set_size_inches(figsize * 1.5)  # scale current size by 1.5
# plt.savefig("/mnt/c/Users/dgoudenege/Desktop/entanglement_matrix.png", dpi=300)
# plt.savefig("/mnt/c/Users/dgoudenege/Desktop/entanglement_matrix.svg")
# # Retrieve ordered ticks label
# newColums = df.columns[cg.dendrogram_col.reordered_ind]
# newIndexs = df.index[cg.dendrogram_row.reordered_ind]
# newData = df.loc[newIndexs,newColums]
# orderedOrg = list(newData.keys())
# OUT = open("/mnt/c/Users/dgoudenege/Desktop/entanglement_matrix.tsv", 'w')
# header = "Protein"
# for org in orderedOrg: header += "\t"+org
# OUT.write(header+"\n")
# for org1 in orderedOrg:
#     line = org1
#     for org2 in orderedOrg: line += "\t"+str(dico[org1][org2]).replace(".",",")
#     OUT.write(line+"\n")
# OUT.close()

# exit()




