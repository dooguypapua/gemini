'''
┌─┐  ┌─┐  ┌┬┐  ┬  ┌┐┌  ┬
│ ┬  ├┤   │││  │  │││  │
└─┘  └─┘  ┴ ┴  ┴  ┘└┘  ┴ phage
---------------------------------------------------
-*- coding: utf-8 -*-                              |
title            : geminiphage.py                  |
description      : gemini phage functions          |
author           : dooguypapua                     |
lastmodification : 20210629                        |
version          : 0.1                             |
python_version   : 3.8.5                           |
---------------------------------------------------
'''

import os
import sys
import re
import shutil
import filecmp
import gzip
import geminiset
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap
from ete3 import NCBITaxa
from typing import Tuple
from tqdm import tqdm
from yaspin import yaspin
from yaspin.spinners import Spinners
from Bio.Seq import Seq
from dna_features_viewer import GraphicFeature, GraphicRecord
from geminini import path_converter, get_input_files, printcolor, get_gemini_path, get_sys_info, fct_checker, cat_lstfiles, read_file
from geminini import dump_json, load_json, reverse_complement, launch_threads, title, exit_gemini, longest_common_substring
from geminiparse import unwrap_fasta, make_blast_dict, make_fasta_dict, make_hmmscan_dict, make_gff_dict
from geminiparse import make_trnascanse_dict, make_interpro_dict, make_eggnog_dict, make_gbk_from_fasta
from geminiparse import blastdict_to_annotsort, hmmscandict_to_dictlt, gbk_to_faa, make_gbk_dict, gbk_to_gff
# from geminiparse import get_refseqacc_from_dict, make_uniprotmapping_dict, make_pvogs_desc_dict
from geminiannot import phanotate, transeq, trnascan_se, pvogs, diamond_p, interproscan, eggnog  # , recombinase
from geminicluster import mmseqs_rbh, make_rbhcluster_dict
from geminidl import get_etag, get_ftp_file_size, get_distant_file_tqdm


@fct_checker
def viridic(pathIN: str, pathOUT: str, ext: str = ".fna") -> Tuple[str, str, str]:
    '''
     ------------------------------------------------------------
    |                          VIRIDIC                           |
    |------------------------------------------------------------|
    |              Phages clustering using VIRIDIC               |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input files or folder (required)      |
    |    pathOUT : path of output files (required)               |
    |    ext     : extension of input files (default=.fna)       |
     ------------------------------------------------------------
    |TOOLS: viridic                                              |
     ------------------------------------------------------------
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "viridic", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: viridic]\nAny input files found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if len(lstFiles) == 1:
        printcolor("[ERROR: viridic]\nViridic required a minimum of two genomes\n", 1, "212;64;89", "None", True)
        exit_gemini()
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'viridic' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['viridic'])
    slurmBool, cpu, memMax, memMin = get_sys_info()
    if pathOUT == "":
        printcolor("[ERROR: viridic]\nMissing '-o'pathOUT\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathOUT = path_converter(pathOUT)
    os.makedirs(pathOUT, exist_ok=True)
    pathLOG = pathOUT+"/viridic.log"
    # Concatenate fasta
    printcolor("♊ Cat files"+"\n")
    pathTMPcat = geminiset.pathTMP+"/concat_viridic.fna"
    strLstError = ""
    for pathFNA in lstFiles:
        file = os.path.basename(pathFNA)
        orgName = file.replace(ext, "").replace("."+ext, "")
        dicoFNA = make_fasta_dict(pathFNA)
        seq = ""
        cptSeq = 0
        for key in dicoFNA:
            seq += dicoFNA[key]
            cptSeq += 1
        if cptSeq != 1:
            strLstError += "[ERROR: viridic]\nNot a single contig genome for '"+file+"'\n"
        TMPCAT = open(pathTMPcat, 'a')
        TMPCAT.write(">"+orgName+"\n"+seq+"\n")
        TMPCAT.close()
    if strLstError != "":
        printcolor(strLstError, 1, "212;64;89", "None", True)
        exit_gemini()
    # Launch viridic
    spinner = yaspin(Spinners.aesthetic, text="♊ VIRIDIC", side="right")
    spinner.start()
    title("VIRIDIC", None)
    cmdViridic = dicoGeminiPath['TOOLS']['viridic']+" in = "+pathTMPcat+" projdir = "+pathOUT+" ncor = "+str(cpu)+" > "+pathLOG+" 2>&1"
    os.system(cmdViridic)
    spinner.stop()
    printcolor("♊ VIRIDIC"+"\n")


@fct_checker
def phage_annotation(pathIN: str, pathOUT: str, boolEMBL: bool = False, enaProject: str = "None", pathTAXO: str = "None", idEvalue: float = "0.01", idThr: int = 30, covThr: int = 50, idThrClust: int = 80, covThrClust: int = 80, ext: str = ".fna") -> Tuple[str, str, bool, str, str, float, int, int, int, int, str]:
    '''
     ------------------------------------------------------------
    |                      PHAGE ANNOTATION                      |
    |------------------------------------------------------------|
    |         Phage syntaxic and functionnal annotation          |
    |         (only for finished genome into one contig)         |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN     : path of input files or folder (required)   |
    |    pathOUT    : path of output files (required)            |
    |    boolEMBL   : output .embl files (default=False)         |
    |    enaProject : ENA project to embl files (default=None)   |
    |    pathTAXO   : path of input taxonomy file (default=None) |
    |    idEvalue   : Evalue threshold (default=0.01)            |
    |    idThr      : %identity threshold (default=30)           |
    |    covThr     : %coverage threshold (default=50)           |
    |    idThrClust : %identity clustering threshold (default=80)|
    |    covThrClust: %coverage clustering threshold (default=80)|
    |    ext        : extension of input files (default=.fna)    |
     ------------------------------------------------------------
     all DMND databases must be formatted as : "identifier|product|organism"
    '''
    pathOUT = path_converter(pathOUT)
    lstFiles, maxpathSize = get_input_files(pathIN, "phage_syntaxic_annotation", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: phage_annotation]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if boolEMBL is True:
        pathTAXO = path_converter(pathTAXO)
        if not os.path.isfile(pathTAXO):
            printcolor("[ERROR: phage_annotation]\nAny taxonomy input file found\n", 1, "212;64;89", "None", True)
            exit_gemini()
        if "PRJEB" not in enaProject:
            printcolor("[ERROR: phage_annotation]\nMissing ENA project for .embl files\n", 1, "212;64;89", "None", True)
            exit_gemini()
    os.makedirs(pathOUT, exist_ok=True)
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    # ***** SYNTAXIC Annotation & TRANSLATION ***** #
    phanotate(pathIN=pathIN, pathOUT=pathOUT, ext=".fna")
    transeq(pathIN=pathOUT, pathOUT=pathOUT, boolOrgName=False, ext=".ffn")
    # ***** tRNA Annotation ***** #
    trnascan_se(pathIN=pathIN, pathOUT=pathOUT, model="-B", ext=".fna")
    dicoTRNA = make_trnascanse_dict(pathIN=pathOUT, pathJSON=pathOUT+"/trnascan_se.json", ext=".trnascanse")
    # ***** Recombinase Annotation ***** #
    # recombinase(pathIN=pathOUT, pathOUT=pathOUT, ext=".faa")
    # dicoREC = make_hmmscan_dict(pathIN=pathOUT, pathJSON=pathOUT+"/recombinase.json", idEvalue=idEvalue, ext=".recomb")
    # dicoREC_LT = hmmscandict_to_dictlt(dicoREC)
    # ***** pVOGS profiles Annotation ***** #
    pvogs(pathIN=pathOUT, pathOUT=pathOUT, ext=".faa")
    dicoPVOGS = make_hmmscan_dict(pathIN=pathOUT, pathJSON=pathOUT+"/pvogs.json", idEvalue=idEvalue, ext=".pvogs")
    # dicoPVOGSdescr = make_pvogs_desc_dict(pathIN=pathOUT+"/pvogs.json", pathJSON=pathOUT+"/pvogs_descr.json")
    dicoPVOGS_LT = hmmscandict_to_dictlt(dicoPVOGS)
    # ***** Nahant collection Annotation ***** #
    diamond_p(pathIN=pathOUT, pathDB=dicoGeminiPath['DATABASES']['nahant_dmnd'],  pathOUT=pathOUT, ext=".faa")
    dicoNAHANT = make_blast_dict(pathIN=pathOUT, dbName="nahant", pathJSON=pathOUT+"/nahant.json", idThr=idThr, minLRthr=covThr, maxLRthr=covThr, ext=".tsv")
    dicoNAHANT_ANNOTSORT = blastdict_to_annotsort(dicoNAHANT)
    # ***** Refseq collection Annotation ***** #
    diamond_p(pathIN=pathOUT, pathDB=dicoGeminiPath['DATABASES']['refseq212bct_dmnd'],  pathOUT=pathOUT, ext=".faa")
    dicoREFSEQ = make_blast_dict(pathIN=pathOUT, dbName="refseq212bct", pathJSON=pathOUT+"/refseq.json", idThr=idThr, minLRthr=covThr, maxLRthr=covThr, ext=".tsv")
    dicoREFSEQ_ANNOTSORT = blastdict_to_annotsort(dicoREFSEQ)
    # ***** phageDB collection Annotation ***** #
    diamond_p(pathIN=pathOUT, pathDB=dicoGeminiPath['DATABASES']['phagedb_dmnd'],  pathOUT=pathOUT, ext=".faa")
    dicoPHAGEDB = make_blast_dict(pathIN=pathOUT, dbName="phagedb", pathJSON=pathOUT+"/phagedb.json", idThr=idThr, minLRthr=covThr, maxLRthr=covThr, ext=".tsv")
    dicoPHAGEDB_ANNOTSORT = blastdict_to_annotsort(dicoPHAGEDB)
    # ***** InterProScan Annotation ***** # (JSON output)
    interproscan(pathIN=pathOUT, pathOUT=pathOUT, ext=".faa")
    dicoIPRSCAN = make_interpro_dict(pathIN=pathOUT+"/interproscan.tsv", idEvalue=idEvalue, pathJSON=pathOUT+"/interproscan.json", ext=".tsv")["interproscan"]
    # ***** EggNOG Annotation ***** #
    eggnog(pathIN=pathOUT, pathOUT=pathOUT, idThr=idThr, covThr=covThr, ext=".faa")
    dicoEGGNOG = make_eggnog_dict(pathIN=pathOUT, pathJSON=pathOUT+"/eggnog.json", ext=".annotations")["eggnog"]
    # ***** RBH clustering ***** #
    mmseqs_rbh(pathIN=pathOUT, pathOUT=pathOUT, idThrClust=idThrClust, covThrClust=covThrClust, ext=".faa")
    make_rbhcluster_dict(pathIN=pathOUT, pathIN2=pathOUT, pathJSON=pathOUT+"/rbh_cluster.json", idThrClust=idThrClust, covThrClust=covThrClust, ext=".rbh", ext2=".faa")
    # ***** GENERATE EMBL ***** #
    if boolEMBL is True:
        # ***** TAXONOMY *****#
        dicoTaxo = {}
        TSV = open(pathTAXO, 'r')
        lstLines = TSV.read().split("\n")[:-1]
        TSV.close()
        for line in lstLines:
            dicoTaxo[line.split("\t")[0]] = line.split("\t")[1]
        # # ***** UNIPROT mapping *****#
        # # Get all refseq accession
        # setRefseqAcc = get_refseqacc_from_dict([dicoNAHANT, dicoREFSEQ, dicoPHAGEDB])
        # # Make uniprot_mapping dictionnary
        # dicoUniprotMapping = make_uniprotmapping_dict(pathIN=dicoGeminiPath['uniprot_mapping'], setRefseqAcc=setRefseqAcc, pathJSON=pathOUT+"/uniprot_mapping.json")
        # ***** Make EMBL file *****#
        for pathFNA in lstFiles:
            orgName = os.path.basename(pathFNA).replace(ext, "").replace("."+ext, "")
            # Read FASTAs
            dicoFNA = make_fasta_dict(pathFNA)
            if len(dicoFNA) > 1:
                exit_gemini("Unable to make embl for more than one contigs FASTA")
            genomeSeq = list(dicoFNA.values())[0]
            dicoFAA = make_fasta_dict(pathOUT+"/"+orgName+".faa")
            dicoFFN = make_fasta_dict(pathOUT+"/"+orgName+".ffn")
            LTprefix = list(dicoFAA.keys())[0].split("_p")[0]+"_p"
            # EMBL description header
            description = orgName.replace("Vibrio_phage_", "Vibrio phage ")+" complete genome."
            wrapper = textwrap.TextWrapper(width=70)
            if len(description) <= 70:
                descr = "DE   "+description+"\n"
            else:
                descr_list = wrapper.wrap(text=description)
                descr = "DE   "+descr_list[0]+"\n"
                for i in range(1, len(descr_list), 1):
                    descr += "DE   "+descr_list[i]+"\n"
            # EMBL taxonomy header
            taxonomy = dicoTaxo[orgName]
            wrapper = textwrap.TextWrapper(width=70)
            if len(taxonomy) <= 70:
                taxo = "OC   "+taxonomy+"\n"
            else:
                taxo_list = wrapper.wrap(text=taxonomy)
                taxo = "OC   "+taxo_list[0]+"\n"
                for i in range(1, len(taxo_list), 1):
                    taxo += "OC   "+taxo_list[i]+"\n"
            # EMBL complete header
            embl_header = "ID   "+LTprefix.replace("_p", "")+"; SV 1; linear; genomic DNA; STD; PHG; "+str(len(genomeSeq))+" BP.\n" +\
                          "XX\n" +\
                          "AC   "+LTprefix.replace("_p", "")+"\n" +\
                          "XX\n" +\
                          "PR   Project:"+enaProject+";\n" +\
                          "XX\n" +\
                          descr +\
                          "XX\n" +\
                          "KW   complete genome.\n" +\
                          "XX\n" +\
                          "OS   "+orgName.replace("_", " ")+"\n" +\
                          taxo +\
                          "XX\n" +\
                          "FH   Key             Location/Qualifiers\n" +\
                          "FH\n"
            # Construct features dictionnary (key = start)
            dicoFeatures = {}
            # TRNA FEATURES
            for trna in dicoTRNA[orgName]:
                if dicoTRNA[orgName][trna]['pseudo'] is False:
                    if dicoTRNA[orgName][trna]['strand'] == 1:
                        strPos = str(dicoTRNA[orgName][trna]['start'])+".."+str(dicoTRNA[orgName][trna]['end'])
                    else:
                        strPos = "complement("+str(dicoTRNA[orgName][trna]['start'])+".."+str(dicoTRNA[orgName][trna]['end'])+")"
                    feature = "FT   gene            "+strPos+"\n" +\
                              "FT                   /locus_tag=\""+LTprefix+"tRNA"+trna.zfill(2)+"\"\n" +\
                              "FT   tRNA            "+strPos+"\n" +\
                              "FT                   /locus_tag=\""+LTprefix+"tRNA"+trna.zfill(2)+"\"\n" +\
                              "FT                   /product=\"tRNA-"+dicoTRNA[orgName][trna]['type'][:3]+"\"\n" +\
                              "FT                   /inference=\"profile:tRNAscan:2.0.9\"\n"
                    dicoFeatures[dicoTRNA[orgName][trna]['start']] = feature
            # CDS FEATURES
            for key in dicoFFN:
                lt = key.split(" ")[0].split("|")[0]
                # GENE POSITIONS
                if dicoFFN[key] in genomeSeq:
                    start = genomeSeq.find(dicoFFN[key])+1
                    end = start+len(dicoFFN[key])-1
                    strPos = str(start)+".."+str(end)
                else:
                    start = genomeSeq.find(reverse_complement(str(Seq(dicoFFN[key]))))+1
                    end = start+len(dicoFFN[key])-1
                    strPos = "complement("+str(start)+".."+str(end)+")"
                # GENE PRODUCT based on blast results (order by confidence)
                validateProduct = set()
                validateUnprecisedProduct = set()
                putativeProduct = set()
                putativeUnprecisedProduct = set()
                foundBlast = False
                validateHitAcc = set()
                for dicoBlastSort in [dicoREFSEQ_ANNOTSORT, dicoNAHANT_ANNOTSORT, dicoPHAGEDB_ANNOTSORT]:
                    if lt in dicoBlastSort[orgName]:
                        foundBlast = True
                        for order in dicoBlastSort[orgName][lt]:
                            for hit in dicoBlastSort[orgName][lt][order]:
                                hitProduct = hit.split("_____")[0]
                                hitAcc = hit.split("_____")[1]
                                pident = float(hit.split("_____")[2])
                                if "hypothetical" not in hitProduct:
                                    if pident >= 75.0:
                                        if "putative" not in hitProduct and "coil containing protein" not in hitProduct and "domain-containing protein" not in hitProduct:
                                            validateProduct.add(hitProduct)
                                            validateHitAcc.add(hitAcc)
                                        else:
                                            validateUnprecisedProduct.add(hitProduct)
                                    elif pident >= 30.0:
                                        if "putative" not in hitProduct and "coil containing protein" not in hitProduct and "domain-containing protein" not in hitProduct:
                                            putativeProduct.add(hitProduct)
                                        else:
                                            putativeUnprecisedProduct.add(hitProduct.replace("putative ", ""))
                if len(validateProduct) > 0:
                    product = list(validateProduct)[0]
                elif len(validateUnprecisedProduct) > 0:
                    product = list(validateUnprecisedProduct)[0]
                elif len(putativeProduct) > 0:
                    product = "putative "+list(putativeProduct)[0]
                elif len(putativeUnprecisedProduct) > 0:
                    product = "putative "+list(putativeUnprecisedProduct)[0]
                elif foundBlast is True:
                    product = "conserved hypothetical protein"
                else:
                    product = "hypothetical protein"
                # GENE FEATURES
                feature = "FT   gene            "+strPos+"\n" +\
                          "FT                   /locus_tag=\""+lt+"\"\n" +\
                          "FT   CDS             "+strPos+"\n" +\
                          "FT                   /locus_tag=\""+lt+"\"\n" +\
                          "FT                   /inference=\"ab initio prediction:PHANOTATE:1.5.0\"\n" +\
                          "FT                   /transl_table=11\n"
                # Format product (split if line>75kr == product>43)
                wrapper = textwrap.TextWrapper(width=43)
                if len(product) <= 43:
                    feature += "FT                   /product=\""+product+"\"\n"
                else:
                    product_list = wrapper.wrap(text=product)
                    feature += "FT                   /product=\""+product_list[0]+"\n"
                    for i in range(1, len(product_list), 1):
                        if i == len(product_list)-1:
                            feature += "FT                   "+product_list[i]+"\"\n"
                        else:
                            feature += "FT                   "+product_list[i]+"\n"
                # DB_XREF based on domain and pvogs RESULTS
                if lt in dicoIPRSCAN:
                    for motif in dicoIPRSCAN[lt]:
                        if "IPR" in motif:
                            feature += "FT                   /db_xref=\"InterPro:"+motif+"\"\n"
                if lt in dicoEGGNOG and dicoEGGNOG[lt]["finalOg"] != "":
                    feature += "FT                   /db_xref=\"EggNOG:"+dicoEGGNOG[lt]["finalOg"]+"\"\n"
                if lt in dicoPVOGS_LT[orgName]:
                    for pvog in dicoPVOGS_LT[orgName][lt]:
                        feature += "FT                   /db_xref=\"pVOGs:"+pvog+"\"\n"
                dicoFeatures[start] = feature
            # EMBL genome sequence
            embl_sequence = "SQ   Sequence "+str(len(genomeSeq))+" BP; " + \
                            str(genomeSeq.count("A"))+" A; "+str(genomeSeq.count("C"))+" C; " + \
                            str(genomeSeq.count("G"))+" G; "+str(genomeSeq.count("T"))+" T; " + \
                            str(genomeSeq.count("N"))+" other;\n"
            wrapper = textwrap.TextWrapper(width=60)
            seq_list = wrapper.wrap(text=genomeSeq)
            cpt = 0
            for i in range(len(seq_list)):
                seqLine = "     "
                for j in range(0, len(seq_list[i]), 10):
                    seqLine += seq_list[i][j:j+10].lower()+" "
                finalPos = cpt+len(seq_list[i])
                seqLine += " "*(80-(len(seqLine)+len(str(finalPos))))+str(finalPos)
                embl_sequence += seqLine+"\n"
                cpt += len(seq_list[i])
            embl_sequence += "//"
            # ***** WRITE .embl file ***** #
            pathEMBLOUT = pathOUT+"/"+orgName+".embl"
            EMBLOUT = open(pathEMBLOUT, 'w')
            EMBLOUT.write(embl_header[:-1]+"\n")
            for start in sorted(dicoFeatures.keys()):
                EMBLOUT.write(dicoFeatures[start][:-1]+"\n")
            EMBLOUT.write("XX\n")
            EMBLOUT.write(embl_sequence+"\n")
            EMBLOUT.close()
    # # ***** CREATE GFF3 files ***** #
    # printcolor("♊ Create GFF3"+"\n")
    # pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    # for pathFNA in lstFiles:
    #     orgName = os.path.basename(pathFNA).replace(ext, "").replace("."+ext, "")
    #     pathGFF = pathOUT+"/"+orgName+".gff"
    #     pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
    #     if not os.path.isfile(pathGFF) or os.path.getsize(pathGFF) == 0:
    #         GFF = open(pathGFF, 'w')
    #         # Retrieve genome sequence
    #         dicoFNA = make_fasta_dict(pathFNA)
    #         genomeSeq = list(dicoFNA.values())[0]
    #         doubleGenomeSeq = list(dicoFNA.values())[0]+list(dicoFNA.values())[0]
    #         # Construct locusTag key pVOGS dict
    #         dicoPVOGS_LT = {}
    #         for pvog in dicoPVOGS[orgName]:
    #             for lt in dicoPVOGS[orgName][pvog].keys():
    #                 if lt not in dicoPVOGS_LT:
    #                     dicoPVOGS_LT[lt] = []
    #                 dicoPVOGS_LT[lt].append(pvog)
    #         # Construct locusTag key recombinase dict
    #         dicoREC_LT = {}
    #         for rec in dicoREC[orgName]:
    #             for lt in dicoREC[orgName][rec].keys():
    #                 if lt not in dicoREC_LT:
    #                     dicoREC_LT[lt] = ""
    #                 dicoREC_LT[lt] += rec+","
    #         for lt in dicoREC_LT:
    #             if dicoREC_LT[lt] != "":
    #                 dicoREC_LT[lt] = dicoREC_LT[lt][: -1]
    #         # GFF3 header
    #         GFF.write("##gff-version\t3.1.26\n")
    #         GFF.write("##sequence-region\t"+orgName+"\t"+str(len(genomeSeq))+"\n")
    #         GFF.write("##species\tVibrio phage\n")
    #         # Add genes
    #         dicoFFN = make_fasta_dict(pathOUT+"/"+orgName+".ffn")
    #         cptGene = 0
    #         for header in dicoFFN:
    #             lt = header.split("|")[0]
    #             frame = header.split("|")[1]
    #             if frame == "1":
    #                 seqNucl = str(Seq(dicoFFN[header]))
    #                 strand = "+"
    #             else:
    #                 seqNucl = reverse_complement(str(Seq(dicoFFN[header])))
    #                 strand = "-"
    #             start = len(doubleGenomeSeq.split(seqNucl)[0])
    #             end = start+len(seqNucl)
    #             # Write gene line
    #             GFF.write(orgName+"\tgemini\tgene\t"+str(start)+"\t"+str(end)+"\t.\t"+strand+"\t0\tID = gene"+str(cptGene)+";locus_tag = "+lt+"\n")
    #             # Write tRNA and CDS line
    #             if header in dicoTRNA[orgName] and dicoTRNA[orgName][header]['pseudo'] is False:
    #                 GFF.write(orgName+"\tgemini\ttRNA\t"+str(start)+"\t"+str(end)+"\t.\t"+strand+"\t0\tID = rna"+str(cptGene)+";Parent = gene"+str(cptGene)+";locus_tag = "+lt+";product = tRNA-"+dicoTRNA[orgName][header]['type']+"\n")
    #             else:
    #                 pvogsRef = pvogsDescr = recombRes = nahantRes = iprscanRes = eggnogRes = ""
    #                 GFFcdsLine = orgName+"\tgemini\tCDS\t"+str(start)+"\t"+str(end)+"\t.\t"+strand+"\t0\tID = rna"+str(cptGene)+";Parent = gene"+str(cptGene)+";locus_tag = "+lt
    #                 # pVOGS
    #                 if lt in dicoPVOGS_LT:
    #                     for pvogsID in dicoPVOGS_LT[lt]:
    #                         pvogsRef += pvogsID+","
    #                         pvogsDescr += ",".join(dicoPVOGSdescr[pvogsID])+","
    #                     pvogsRef = pvogsRef[: -1]
    #                     pvogsDescr = pvogsDescr[: -1]
    #                 # recombinase
    #                 if lt in dicoREC_LT:
    #                     recombRes = dicoREC_LT[lt]
    #                 # Nahant
    #                 if lt in dicoNAHANT[orgName]:
    #                     for key in dicoNAHANT[orgName][lt]:
    #                         if key != "length":
    #                             nahantRes += key.replace(",", " ").replace(";", " ")+","
    #                     nahantRes = nahantRes[: -1]
    #                 # InterProScan
    #                 if lt in dicoIPRSCAN["interproscan"]:
    #                     for key in dicoIPRSCAN["interproscan"][lt].keys():
    #                         iprscanRes += key+":"+dicoIPRSCAN["interproscan"][lt][key]['descr']+","
    #                     iprscanRes = iprscanRes[: -1]
    #                 # EggNOG
    #                 if lt in dicoEGGNOG['eggnog'] and dicoEGGNOG['eggnog'][lt]['descr'] != "-":
    #                     eggnogRes = dicoEGGNOG['eggnog'][lt]['descr']
    #                     if dicoEGGNOG['eggnog'][lt]['prefName'] != "-":
    #                         eggnogRes += "("+dicoEGGNOG['eggnog'][lt]['prefName']+")"
    #                 # Write line
    #                 GFFcdsLine += ";pvogsRef = "+pvogsRef+";pvogsDescr = "+pvogsDescr+";recomb = "+recombRes+";nahant = "+nahantRes+";iprscan = "+iprscanRes+";eggnog = "+eggnogRes+"\n"
    #                 GFF.write(GFFcdsLine.replace("  ", " "))
    #             cptGene += 1
    #         GFF.close()
    #     pbar.update(1)
    #     title("Create GFF3", pbar)
    # pbar.close()
    # # ***** CREATE TABLE from GFF ***** #
    # gff_to_table(pathIN=pathOUT, pathOUT=pathOUT, format=".xlsx", maxWidth=50, ext=".gff")


@fct_checker
def phageDB(pathIN: str, pathOUT: str, checkvHQ: float = 75.0) -> Tuple[str, str, float]:
    '''
     ------------------------------------------------------------
    |                   GENBANK PHAGE DATABASE                   |
    |------------------------------------------------------------|
    |                Make GenBank phages database                |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input GV assembly folder (required)   |
    |    pathOUT : path of output files (required)               |
    |    checkv  : checkV High/MediumQuality % (default=75)      |
     ------------------------------------------------------------
    |TOOLS: seqret, python, hmmscan                              |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    if not os.path.isdir(pathIN):
        printcolor("[ERROR: phageDB]\nAny input GV assembly folder found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    # Create folder tree
    os.makedirs(pathOUT, exist_ok=True)
    os.makedirs(pathOUT+"/FNA", exist_ok=True)
    os.makedirs(pathOUT+"/FNA_filtered", exist_ok=True)
    os.makedirs(pathOUT+"/GBK", exist_ok=True)
    os.makedirs(pathOUT+"/FFN", exist_ok=True)
    os.makedirs(pathOUT+"/FAA", exist_ok=True)
    os.makedirs(pathOUT+"/CHECKV", exist_ok=True)
    os.makedirs(pathOUT+"/ANNOT", exist_ok=True)
    # Paths
    pathDBJSON = pathOUT+"/phageDB.json"
    pathLOG = pathOUT+"/log.out"
    pathGenomeReportVirusesTXT = pathOUT+"/genbank_genome_reports_viruses.txt"
    pathGenomeReportVirusesETAG = pathOUT+"/genbank_genome_reports_viruses.etag"
    pathGBKsummaryTXT = pathOUT+"/assembly_summary_genbank.txt"
    pathGBKsummaryETAG = pathOUT+"/assembly_summary_genbank.etag"
    pathENAsummaryTSV = pathOUT+"/assembly_summary_ena.tsv"
    pathCONCATFNA = pathOUT+"/all_filtered_phages.fna"
    pathCONCATFAA = pathOUT+"/all_filtered_phages.faa"
    pathMAPPINGPROT = pathOUT+"/all_filtered_phages_mapping.csv"
    pathExcludedGCA = pathOUT+"/excluded_gca.txt"
    pathTMP = geminiset.pathTMP
    # URLs
    urlGBKreportVirus = "https: //ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/viruses.txt"
    urlGBKFTP = "ftp.ncbi.nlm.nih.gov"
    urlGBKFTPreport = "genomes/GENOME_REPORTS/viruses.txt"
    urlGBKFTPsummary = "genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
    urlGBKHTTPsummary = "https: //ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
    # Gemini variables
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'seqret' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['seqret'])
    if 'python' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['python'])
    if 'hmmscan' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['hmmscan'])
    slurmBool, cpu, memMax, memMin = get_sys_info()
    lstBadKr = ["#", "'", "(", ")", ",", "/", ": ", " "]
    # Summary JSON
    if not os.path.isfile(pathDBJSON) or os.path.getsize(pathDBJSON) == 0:
        dicoDB = {}
    else:
        dicoDB = load_json(pathDBJSON)
    # Excluded GCA
    lstExcludedGCA = []
    if os.path.isfile(pathExcludedGCA):
        IN = open(pathExcludedGCA, 'r')
        lstLines = IN.read().split("\n")
        IN.close()
        for line in lstLines:
            if line != "":
                lstExcludedGCA.append(line)
    # ***** RETRIEVE ALL VIRUS WITH BACTERIA HOST ***** #
    # Get current local/distant etag (or None if empty)
    previousEtag = get_etag(pathGenomeReportVirusesETAG, None)
    newEtag = get_etag(None, urlGBKreportVirus)
    # Download if new distant etag
    if previousEtag == newEtag:
        printcolor("♊ Genbank genome reports   [up-to-date ]"+"\n")
    else:
        printcolor("♊ Genbank genome reports   [downloading]"+"\n")
        # Get filesize for progress bar
        filesize = get_ftp_file_size(urlGBKFTP, urlGBKFTPreport)
        # Download file
        get_distant_file_tqdm(urlGBKreportVirus, pathGenomeReportVirusesTXT, filesize, 1024)
    # Write new etag
    ETAG = open(pathGenomeReportVirusesETAG, 'w')
    ETAG.write(newEtag)
    ETAG.close()
    # Retrieve species list (containing phage or virus)
    REPORT = open(pathGenomeReportVirusesTXT, 'r')
    lstLines = REPORT.read().split("\n")[1: -1]
    REPORT.close()
    setReportPhageSpecie = set()
    for line in lstLines:
        splitLine = line.split("\t")
        orgName = splitLine[0]
        host = splitLine[8]
        if host == "bacteria":
            if "phage" in orgName:
                specie = orgName.split("phage")[0]+"phage"
            elif "virus" in orgName:
                specie = orgName.split("virus")[0]+"virus"
            elif "Phage" in orgName:
                specie = orgName.split("Phage")[0]+"Phage"
            elif "Virus" in orgName:
                specie = orgName.split("Pirus")[0]+"Virus"
            else:
                specie = orgName
            setReportPhageSpecie.add(specie)
    # ***** GENBANK ASSEMBLY SUMMARY ***** #
    # Get current local/distant etag (or None if empty)
    previousEtag = get_etag(pathGBKsummaryETAG, None)
    newEtag = get_etag(None, urlGBKHTTPsummary)
    gbkUpdate = False
    # Download if new distant etag
    if previousEtag == newEtag:
        printcolor("♊ Genbank assembly summary [up-to-date ]"+"\n")
    else:
        printcolor("♊ Genbank assembly summary [downloading]"+"\n")
        # Get filesize for progress bar
        filesize = get_ftp_file_size(urlGBKFTP, urlGBKFTPsummary)
        # Download file
        get_distant_file_tqdm(urlGBKHTTPsummary, pathTMP+"/assembly_summary_genbank.txt", filesize, 1024)
        # Write new etag
        ETAG = open(pathGBKsummaryETAG, 'w')
        ETAG.write(newEtag)
        ETAG.close()
        gbkUpdate = True
    # Filtering Genbank assembly
    printcolor("♊ Filter Genbank"+"\n")
    if gbkUpdate is True and not os.path.isfile(pathGBKsummaryTXT):
        TMPSUMMARY = open(pathTMP+"/assembly_summary_genbank.txt", 'r')
        lstLines = TMPSUMMARY.read().split("\n")
        TMPSUMMARY.close()
        SUMMARY = open(pathGBKsummaryTXT, 'w')
        pbar = tqdm(total=len(lstLines), ncols=50, leave=False, desc="", file=sys.stdout, bar_format="   Filtering {percentage: 3.0f}%|{bar}|")
        # Init NCBITaxa
        ncbi = NCBITaxa()
        for line in lstLines:
            if line == "" or line[0] == "#":
                SUMMARY.write(line+"\n")
            else:
                isPhage = False
                splitLine = line.split("\t")
                taxID = int(splitLine[5])
                specieID = int(splitLine[6])
                # Search by organism name
                if "phage" in line.lower():
                    try:
                        lineage = ncbi.get_lineage(taxID)
                    except ValueError:
                        try:
                            lineage = ncbi.get_lineage(specieID)
                        except ValueError:
                            if " phage " in line:
                                lineage = [10239]
                    if 10239 in lineage:
                        isPhage = True
                # Search using genome reports
                else:
                    for specie in setReportPhageSpecie:
                        if specie.lower() in line.lower():
                            isPhage = True
                            break
                if isPhage is True:
                    SUMMARY.write(line+"\n")
            pbar.update(1)
            title("Filter Genbank", pbar)
        pbar.close()
        SUMMARY.close()
    # Parse Genbank assembly
    printcolor("♊ Genbank summary"+"\n")
    # if gbkUpdate is True or not os.path.isfile(pathDBJSON):
    SUMMARY = open(pathGBKsummaryTXT, 'r')
    lstLines = SUMMARY.read().split("\n")
    SUMMARY.close()
    pbar = tqdm(total=len(lstLines), ncols=50, leave=False, desc="", file=sys.stdout, bar_format="   Parsing {percentage: 3.0f}%|{bar}|")
    for line in lstLines:
        if line != "" and line[0] != "#":
            splitLine = line.split("\t")
            gca = splitLine[0].split(".")[0]
            version = int(splitLine[0].split(".")[1])
            if not gca+"."+str(version) in lstExcludedGCA:
                # Format orgName
                orgName = splitLine[7]
                infName = splitLine[8].replace("strain = ", "")
                isoName = splitLine[9]
                for badKr in lstBadKr:
                    orgName = orgName.replace(badKr, "_")
                # Add to dico
                if orgName not in dicoDB:
                    dicoDB[orgName] = {}
                if gca not in dicoDB[orgName]:
                    dicoDB[orgName][gca] = {}
                if 'genbank' not in dicoDB[orgName][gca] or dicoDB[orgName][gca]['genbank']['version'] < version:
                    dicoDB[orgName][gca]['genbank'] = {'infName': infName, 'isoName': isoName, 'url': splitLine[19], 'version': version}
        pbar.update(1)
        title("Genbank summary", pbar)
    pbar.close()
    # ***** ENA ASSEMBLY SUMMARY ***** #
    printcolor("♊ ENA summary"+"\n")
    # Download
    spinner = yaspin(Spinners.aesthetic, text="♊ ENA assembly summary     ", side="right")
    spinner.start()
    title("ENA", None)
    if not os.path.isfile(pathENAsummaryTSV) or os.path.getsize(pathENAsummaryTSV) == 0:
        cmdENA = "curl -s -X POST -H \"Content-Type: application/x-www-form-urlencoded\" -d 'result = assembly&query = tax_tree(10239)&fields = accession%2Cversion%2Cscientific_name%2Cstrain%2Cstudy_name%2Ctax_id&format = tsv' \"https: //www.ebi.ac.uk/ena/portal/api/search\" | grep \"phage\" > "+pathENAsummaryTSV
        os.system(cmdENA)
        spinner.stop()
        printcolor("♊ ENA assembly summary     [downloading]"+"\n")
    else:
        # Compare File
        cmdENA = "curl -s -X POST -H \"Content-Type: application/x-www-form-urlencoded\" -d 'result = assembly&query = tax_tree(10239)&fields = accession%2Cversion%2Cscientific_name%2Cstrain%2Cstudy_name%2Ctax_id&format = tsv' \"https: //www.ebi.ac.uk/ena/portal/api/search\" | grep \"phage\" > "+pathTMP+"/assembly_summary_ena.tsv"
        os.system(cmdENA)
        sameFile = filecmp.cmp(pathENAsummaryTSV, pathTMP+"/assembly_summary_ena.tsv", shallow=False)
        spinner.stop()
        if sameFile is True:
            printcolor("♊ ENA assembly summary     [up-to-date ]"+"\n")
        else:
            printcolor("♊ ENA assembly summary     [downloading]"+"\n")
            shutil.copyfile(pathTMP+"/assembly_summary_ena.tsv", pathENAsummaryTSV)
    # Parse
    ENA = open(pathOUT+"/assembly_summary_ena.tsv", 'r')
    lstLines = ENA.read().split("\n")[: -1]
    ENA.close()
    # Parse
    pbar = tqdm(total=len(lstLines), ncols=50, leave=False, desc="", file=sys.stdout, bar_format="   Parsing {percentage: 3.0f}%|{bar}|")
    for line in lstLines:
        splitLine = line.split("\t")
        gca = splitLine[0]
        version = int(splitLine[1])
        if not gca+"."+str(version) in lstExcludedGCA:
            orgName = splitLine[2]
            for badKr in lstBadKr:
                orgName = orgName.replace(badKr, "_")
            infName = splitLine[3]
            isoName = splitLine[4]
            taxID = int(splitLine[5])
            url = "https: //www.ebi.ac.uk/ena/browser/api/fasta/"+gca+"."+str(version)+"?download = true&gzip = true"
            if orgName not in dicoDB:
                dicoDB[orgName] = {}
            if gca not in dicoDB[orgName]:
                dicoDB[orgName][gca] = {}
            if 'ena' not in dicoDB[orgName][gca] or dicoDB[orgName][gca]['ena']['version'] < version:
                dicoDB[orgName][gca]['ena'] = {'infName': infName, 'isoName': isoName, 'url': url, 'version': version}
        pbar.update(1)
        title("ENA summary", pbar)
    pbar.close()
    # ***** GV ASSEMBLIES ***** #
    printcolor("♊ GV team genomes"+"\n")
    for gvFNA in os.listdir(pathIN):
        orgName = os.path.basename(gvFNA).replace(".fna", "").replace("-", ".")
        if orgName not in dicoDB:
            dicoDB[orgName] = {}
        dicoDB[orgName]['gv'] = {'infName': "", 'isoName': "", 'url': pathIN+"/"+gvFNA, 'version': None}
        pathFNA = pathOUT+"/FNA/"+orgName+".fna"
        if not os.path.isfile(pathFNA) or os.path.getsize(pathFNA) == 0:
            shutil.copyfile(pathIN+"/"+gvFNA, pathFNA)
    # ***** CHOOSE SOURCE per organism (for not GV genomes)***** #
    maxOrgNameSize = 0
    for orgName in dicoDB:
        maxOrgNameSize = max(maxOrgNameSize, len(orgName))
        if "gv" not in dicoDB[orgName]:
            for gca in dicoDB[orgName]:
                if "genbank" not in dicoDB[orgName][gca]:
                    choice = "ena"
                elif "ena" not in dicoDB[orgName][gca]:
                    choice = "genbank"
                elif dicoDB[orgName][gca]["ena"]['version'] > dicoDB[orgName][gca]["genbank"]['version']:
                    choice = "ena"
                else:
                    choice = "genbank"
                # New genomes
                if "source_choice" not in dicoDB[orgName][gca]:
                    dicoDB[orgName][gca]["source_choice"] = choice
                # Update genomes
                if dicoDB[orgName][gca]["source_choice"] != choice:
                    dicoDB[orgName][gca]["source_choice"] = choice
                    # delete FNA, FAA and FFN
                    os.system("rm -f "+pathOUT+"/F*/"+orgName+"_"+gca+".*.f*")
    # Download FNA files
    nbGCA = 0
    for org in dicoDB:
        nbGCA += len(dicoDB[org])
    printcolor("♊ Retrieve FNA files"+"\n")
    pbar = tqdm(total=nbGCA, ncols=50+maxOrgNameSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    cptNewFNA = 0
    cptNewGBK = 0
    lstFailedFNA = []
    lstFailedGBK = []
    for orgName in dicoDB:
        if "gv" not in dicoDB[orgName]:
            for gca in dicoDB[orgName]:
                srcChoice = dicoDB[orgName][gca]["source_choice"]
                pbar.set_description_str(orgName+" ".rjust(maxOrgNameSize-len(orgName)))
                gcaVersion = dicoDB[orgName][gca][srcChoice]['version']
                # check and remove previous version
                if gcaVersion != 1:
                    for i in range(1, gcaVersion, 1):
                        pathFNAold = pathOUT+"/FNA/"+orgName+"_"+gca+"."+str(i)+".fna"
                        if os.path.isfile(pathFNAold):
                            os.remove(pathFNAold)
                # Download FNA if required
                pathFNA = pathOUT+"/FNA/"+orgName+"_"+gca+"."+str(gcaVersion)+".fna"
                if not os.path.isfile(pathFNA) or os.path.getsize(pathFNA) == 0:
                    if srcChoice == "genbank":
                        pathTMPgzfna = pathTMP+"/"+orgName+"_"+gca+"."+str(gcaVersion)+".fna.gz"
                        cmdWGET = "wget -O - -q "+dicoDB[orgName][gca][srcChoice]['url']+"/"+os.path.basename(dicoDB[orgName][gca][srcChoice]['url'])+"_genomic.fna.gz > "+pathTMPgzfna
                        os.system(cmdWGET)
                        if os.path.isfile(pathTMPgzfna) and os.path.getsize(pathTMPgzfna) != 0:
                            os.system("gzip -c -d "+pathTMPgzfna+" > "+pathFNA)
                    else:
                        cmdWGET = "wget -O "+pathFNA+" -q "+dicoDB[orgName][gca][srcChoice]['url']
                        os.system(cmdWGET)
                    # Unwrap
                    if not os.path.isfile(pathFNA) or os.path.getsize(pathFNA) == 0:
                        lstFailedFNA.append(os.path.basename(pathFNA))
                    else:
                        unwrap_fasta(pathIN=pathFNA, ext=".fna")
                        cptNewFNA += 1
                # Download GBFF if required
                pathGBK = pathOUT+"/GBK/"+orgName+"_"+gca+"."+str(gcaVersion)+".gbk.gz"
                if not os.path.isfile(pathGBK) or os.path.getsize(pathGBK) == 0:
                    if srcChoice == "genbank":
                        cmdWGET = "wget -O "+pathGBK+" -q "+dicoDB[orgName][gca][srcChoice]['url']+"/"+os.path.basename(dicoDB[orgName][gca][srcChoice]['url'])+"_genomic.gbff.gz"
                        os.system(cmdWGET)
                    else:
                        cmdWGET = "wget -O "+pathTMP+"/"+gca+"."+str(gcaVersion)+".embl -q "+dicoDB[orgName][gca][srcChoice]['url'].replace("/fasta/", "/embl/").replace("?download = true&gzip = true", "")
                        os.system(cmdWGET)
                        # Convert embl to gbk
                        cmdSEQRET = dicoGeminiPath['TOOLS']['seqret']+" -sequence "+pathTMP+"/"+gca+"."+str(gcaVersion)+".embl -sformat1 embl -outseq "+pathGBK.replace(".gz", "")+" -osformat2 genbank >> "+pathLOG+" 2>&1"
                        os.system(cmdSEQRET)
                        # Compress gbk
                        if os.path.isfile(pathGBK.replace(".gz", "")) and os.path.getsize(pathGBK.replace(".gz", "")):
                            os.system("gzip "+pathGBK.replace(".gz", ""))
                    if not os.path.isfile(pathGBK) or os.path.getsize(pathGBK) == 0:
                        lstFailedGBK.append(os.path.basename(pathGBK))
                    else:
                        cptNewGBK += 1
        pbar.update(1)
        title("Retrieve FNA", pbar)
    pbar.close()
    # ***** FILTER LOW QUALITY GENOME ***** #
    printcolor("♊ Filter low quality genome"+"\n")
    lstFNA = os.listdir(pathOUT+"/FNA")
    lstFiltered = []
    pbar = tqdm(total=len(lstFNA), ncols=75, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for fna in lstFNA:
        if fna not in lstFailedFNA:
            if "GCA" in fna:
                orgName = fna.split("_GCA")[0]
                gcaAccession = "GCA_"+fna.replace(".fna", "").split("_")[-1].split(".")[0]
                pbar.set_description_str("GCA_"+fna.replace(".fna", "").split("_")[-1])
            else:
                orgName = fna.replace(".fna", "")
                gcaAccession = "gv"
                pbar.set_description_str(fna.replace(".fna", "")[-15:])
            pathFNA = pathOUT+"/FNA/"+fna
            pathCHECKVOUT = pathOUT+"/CHECKV/"+fna.replace(".fna", "")
            # Launch if required
            if not os.path.isfile(pathCHECKVOUT+"/quality_summary.tsv") or os.path.getsize(pathCHECKVOUT+"/quality_summary.tsv") == 0:
                os.system("echo \""+dicoGeminiPath['TOOLS']['python']+" "+dicoGeminiPath['TOOLS']['checkv']+" end_to_end "+pathFNA+" "+pathCHECKVOUT+" -t "+str(cpu)+" -d "+dicoGeminiPath['DATABASES']['checkv_db']+"\" >> "+pathLOG)
                cmdCHECKV = dicoGeminiPath['TOOLS']['python']+" "+dicoGeminiPath['TOOLS']['checkv']+" end_to_end "+pathFNA+" "+pathCHECKVOUT+" -t "+str(cpu)+" -d "+dicoGeminiPath['DATABASES']['checkv_db']+" >> "+pathLOG+" 2>&1"
                os.system(cmdCHECKV)
                # Delete tmp folder
                shutil.rmtree(pathCHECKVOUT+"/tmp")
            # Parse quality_summary file
            if "check" not in dicoDB[orgName][gcaAccession]:
                TSV = open(pathCHECKVOUT+"/quality_summary.tsv", 'r')
                lstLines = TSV.read().split("\n")[1:-1]
                TSV.close()
                # % of high-quality contig
                dicoCountQuality = {}
                for line in lstLines:
                    splitLine = line.split("\t")
                    checkv_quality = splitLine[7]  # Complete/High-quality/Medium-quality/Low-quality/Not-determined
                    try:
                        dicoCountQuality[checkv_quality] += 1
                    except KeyError:
                        dicoCountQuality[checkv_quality] = 1
                # Add to dicoDB
                dicoDB[orgName][gcaAccession]['check'] = {}
                for quality in dicoCountQuality:
                    dicoDB[orgName][gcaAccession]['check'][quality] = (dicoCountQuality[quality]*100)/len(lstLines)
            # Filtering
            percentQuality = 0
            if "Complete" in dicoDB[orgName][gcaAccession]['check']:
                percentQuality += dicoDB[orgName][gcaAccession]['check']['Complete']
            if "High-quality" in dicoDB[orgName][gcaAccession]['check']:
                percentQuality += dicoDB[orgName][gcaAccession]['check']['High-quality']
            if "Medium-quality" in dicoDB[orgName][gcaAccession]['check']:
                percentQuality += dicoDB[orgName][gcaAccession]['check']['Medium-quality']
            if percentQuality < 75.0:
                lstFiltered.append(gcaAccession)
            elif not os.path.isfile(pathOUT+"/FNA_filtered/"+fna):
                shutil.copyfile(pathFNA, pathOUT+"/FNA_filtered/"+fna)
        pbar.update(1)
        title("Filter lowqual", pbar)
    pbar.close()
    # Dump summary JSON
    dump_json(dicoDB, pathDBJSON)
    # ***** STATISTICS ***** #
    printcolor("♊ phageDB statistics"+"\n")
    printcolor("⏩ Found "+str(len(os.listdir(pathOUT+"/FNA")))+" FNA"+"\n")
    printcolor("⏩ Found "+str(len(os.listdir(pathOUT+"/FNA_filtered")))+" FNA filtered"+"\n")
    printcolor("⏩ Found "+str(len(os.listdir(pathOUT+"/GBK")))+" GBK"+"\n")
    if cptNewFNA == 0:
        printcolor("⏩ Any new FNA downloaded"+"\n")
    else:
        printcolor("⏩ "+str(cptNewFNA)+" new FNA downloaded"+"\n")
    if cptNewGBK == 0:
        printcolor("⏩ Any new GBK downloaded"+"\n")
    else:
        printcolor("⏩ "+str(cptNewGBK)+" new GBK downloaded"+"\n")
    if len(lstFiltered) > 0:
        printcolor("⛔ "+str(len(lstFiltered))+" checkV excluded genomes"+"\n")
    if len(lstFailedFNA) > 0:
        printcolor("⛔ "+str(len(lstFailedFNA))+" failed FNA download ("+", ".join(lstFailedFNA)+")"+"\n")
    if len(lstFailedGBK) > 0:
        printcolor("⛔ "+str(len(lstFailedGBK))+" failed GBK download ("+", ".join(lstFailedGBK)+")"+"\n")
    # ***** CONCATENATE ALL FNA (only filtered) ***** #
    printcolor("♊ Concatenate filtered FNA"+"\n")
    if cptNewFNA != 0 or not os.path.isfile(pathCONCATFNA) or os.path.getsize(pathCONCATFNA) == 0:
        lstFNA = os.listdir(pathOUT+"/FNA_filtered")
        pbar = tqdm(total=len(lstFNA), ncols=75, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
        CONCAT = open(pathCONCATFNA, 'w')
        for fna in lstFNA:
            if "GCA" in fna:
                pbar.set_description_str("GCA_"+fna.replace(".fna", "").split("_")[-1])
            else:
                pbar.set_description_str(fna.replace(".fna", "")[-15:])
            dicoFNA = make_fasta_dict(pathOUT+"/FNA_filtered/"+fna)
            seq = ""
            for key in dicoFNA:
                seq += dicoFNA[key]+"N"*100
            CONCAT.write(">"+fna.replace(".fna", "")+"\n"+seq[:-100]+"\n")
            pbar.update(1)
            title("Cat FNA", pbar)
        pbar.close()
        CONCAT.close()
    # ***** SYNTAXIC Annotation & TRANSLATION ***** #
    phanotate(pathIN=pathOUT+"/FNA_filtered", pathOUT=pathOUT+"/FFN", fromPhageDb=True, ext=".fna")
    transeq(pathIN=pathOUT+"/FFN", pathOUT=pathOUT+"/FAA", fromPhageDb=True, ext=".ffn")
    # ***** CONCATENATE ALL FAA ***** #
    printcolor("♊ Concatenate all FAA"+"\n")
    if cptNewFNA != 0 or not os.path.isfile(pathCONCATFAA) or os.path.getsize(pathCONCATFAA) == 0 or not os.path.isfile(pathMAPPINGPROT) or os.path.getsize(pathMAPPINGPROT) == 0:
        lstFAA = os.listdir(pathOUT+"/FAA")
        pbar = tqdm(total=len(lstFAA), ncols=75, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
        CONCAT = open(pathCONCATFAA, 'w')
        MAPPING = open(pathMAPPINGPROT, 'w')
        MAPPING.write("protein_id, contig_id, keywords\n")
        for faa in lstFAA:
            if "GCA" in faa:
                pbar.set_description_str("GCA_"+faa.replace(".faa", "").split("_")[-1])
            else:
                pbar.set_description_str(faa.replace(".faa", "")[-15:])
            pbar.set_description_str("GCA_"+faa.replace(".faa", "").split("_")[-1])
            dicoFAA = make_fasta_dict(pathOUT+"/FAA/"+faa)
            pathFFN = pathOUT+"/FFN/"+faa.replace(".faa", ".ffn")
            if os.path.isfile(pathFFN):
                for key in dicoFAA:
                    CONCAT.write(">"+key+"\n"+dicoFAA[key]+"\n")
                    MAPPING.write(key.split(" [")[0]+", "+key.split(" [")[1].replace("]", "")+", "+key.split("|")[0]+"\n")
            pbar.update(1)
            title("Cat FAA", pbar)
        pbar.close()
        CONCAT.close()
        MAPPING.close()
    # ***** SEARCH LARGE TERMINASE SUBUNIT ***** #
    spinner = yaspin(Spinners.aesthetic, text="♊ Search large terminase", side="right")
    spinner.start()
    title("Terminase", None)
    lstFAA = os.listdir(pathOUT+"/FAA")
    dicoThread = {}
    for faa in lstFAA:
        pathTBLOUT = pathOUT+"/ANNOT/"+faa.replace(".faa", ".tblout")
        if not os.path.isfile(pathTBLOUT) or os.path.getsize(pathTBLOUT) == 0:
            cmdHmmscan = dicoGeminiPath['TOOLS']['hmmscan']+" --tblout "+pathTBLOUT+" "+dicoGeminiPath['DATABASES']['terminase_hmm']+" "+pathOUT+"/FAA/"+faa
            dicoThread[faa.replace(".faa", "")] = {"cmd": cmdHmmscan, "returnstatut": None, "returnlines": []}
    if len(dicoThread) > 0:
        launch_threads(dicoThread, "hmmscan_terminase", cpu, pathTMP)
    spinner.stop()
    printcolor("♊ Search large terminase"+"\n")


@fct_checker
def phageDBsearch(pathIN: str, pathOUT: str, idThr: int = 20, covThr: int = 50, dmndDB: str = "genbank") -> Tuple[str, str, int, int, str]:
    '''
     ------------------------------------------------------------
    |                  SEARCH IN PHAGE DATABASE                  |
    |------------------------------------------------------------|
    |    Search in genbank source or phanotate phage database    |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of query protein file (required)         |
    |    pathOUT : path of search output dir (required)          |
    |    idThr   : %identity threshold (default=20)              |
    |    covThr  : %coverage threshold (default=50)              |
    |    dmndDB  : database genbank|phanotate (default=phanotate)|
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    if not os.path.isfile(pathIN):
        printcolor("[ERROR: phageDBsearch]\nAny input FAA file\n", 1, "212;64;89", "None", True)
        exit_gemini()
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if dmndDB == "genbank":
        pathDB = dicoGeminiPath['DATABASES']['phagedb_dmnd']
    elif dmndDB == "phanotate":
        pathDB = dicoGeminiPath['DATABASES']['phagedb_dmnd_phanotate']
    else:
        printcolor("[ERROR: phageDBsearch]\nInvalid dmnd database, must be \"genbank\" or \"phanotate\"\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathTMP = geminiset.pathTMP
    diamond_p(pathIN=pathIN, pathDB=pathDB, pathOUT=pathTMP+"/diamond.out", boolSeq=True, ext="")
    lstLines = read_file(pathTMP+"/diamond.out", excludeFirstKr=None, yaspinBool=False)
    if len(lstLines) <= 3:
        printcolor("⛔ Any protein found"+"\n")
    else:
        printcolor("⏩ Found "+str(len(lstLines)-3)+" protein(s)"+"\n")
        os.makedirs(pathOUT, exist_ok=True)
        shutil.move(pathTMP+"/diamond.out", pathOUT+"/diamond.out")
        dicoPerQuery = {}
        splitHeader = lstLines[2].replace("# Fields: ", "").split(", ")
        dicoHeader = {}
        for i in range(len(splitHeader)):
            dicoHeader[splitHeader[i]] = i
        for line in lstLines[3:]:
            splitLine = line.split("\t")
            query = splitLine[dicoHeader["Query ID"]]
            subject = splitLine[dicoHeader["Subject title"]]
            pident = float(splitLine[dicoHeader["Percentage of identical matches"]])
            # alen = int(splitLine[dicoHeader["Alignment length"]])
            qlen = int(splitLine[dicoHeader["Query length"]])
            qstart = int(splitLine[dicoHeader["Start of alignment in query"]])
            qend = int(splitLine[dicoHeader["End of alignment in query"]])
            # slen = int(splitLine[dicoHeader["Subject length"]])
            sstart = int(splitLine[dicoHeader["Start of alignment in subject"]])
            send = int(splitLine[dicoHeader["End of alignment in subject"]])
            # evalue = float(splitLine[dicoHeader["Expected value"]])
            sseq = splitLine[dicoHeader["Subject sequence"]]
            qcov = round(((max(qstart, qend)-min(qstart, qend))*100/qlen), 1)
            scov = round(((max(sstart, send)-min(sstart, send))*100/qlen), 1)
            toWrite = ">"+subject+"|query="+query+"|pident="+str(pident)+"|qcov="+str(qcov)+"|scov="+str(scov)+"\n"+sseq
            try:
                dicoPerQuery[query].add(toWrite)
            except KeyError:
                dicoPerQuery[query] = set([toWrite])
        # Write per query results FASTA
        for query in dicoPerQuery:
            OUT = open(pathOUT+"/"+query+"_phageDBsearch.faa", 'w')
            for toWrite in dicoPerQuery[query]:
                OUT.write(toWrite)
            OUT.close()


@fct_checker
def myVIRIDIC(pathIN: str, pathOUT: str, ref: str = "None", thfam: float = 50.0, thgen: float = 70.0, thsp: float = 95.0, boolFromDB: bool = False, ext: str = ".fna") -> Tuple[str, str, str, float, float, float, bool, str]:
    '''
     ------------------------------------------------------------
    |                         myVIRIDIC                          |
    |------------------------------------------------------------|
    |                Modified version of VIRIDIC                 |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input FNA folder (required)           |
    |    pathOUT : path of output folder (required)              |
    |    ref     : reference organism name (default=None)        |
    |    thfam   : threshold for family clustering (default=50)  |
    |    thgen   : threshold for genus clustering (default=70)   |
    |    thsp    : threshold for species clustering (default=95) |
    |    ext     : extension of input files (default=.fna)       |
     ------------------------------------------------------------
    |TOOLS: blastn, rscript                                      |
     ------------------------------------------------------------
    '''
    setAllOrg = set()
    if boolFromDB:
        printcolor("♊ Read input from DB"+"\n")
        lstFiles = []
        maxpathSize = 15
        lst_genome = os.listdir(pathIN)
        pbar = tqdm(total=len(lst_genome), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
        for folder in lst_genome:
            if os.path.isfile(pathIN+"/"+folder+"/"+folder+"_genomic.fna.gz"):
                lstFiles.append(pathIN+"/"+folder+"/"+folder+"_genomic.fna.gz")
                setAllOrg.add(folder)
            pbar.update(1)
            title("Read db", pbar)
        pbar.close()
    else:
        lstFiles, maxpathSize = get_input_files(pathIN, "myVIRIDIC", [ext])
        for file in lstFiles:
            orgName = os.path.basename(file).replace("_genomic.fna.gz", "").replace(ext, "").replace("."+ext, "")
            setAllOrg.add(orgName)
    if len(lstFiles) == 0:
        printcolor("[ERROR: myVIRIDIC]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if ref != "None" and ref not in setAllOrg:
        printcolor("[ERROR: myVIRIDIC]\nReference \""+ref+"\" not found in input files\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathOUT = path_converter(pathOUT)
    os.makedirs(pathOUT, exist_ok=True)
    if "." not in ext:
        ext = "."+ext
    # Gemini variables
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'blastn' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['blastn'])
    if 'rscript' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['rscript'])
    slurmBool, cpu, memMax, memMin = get_sys_info()
    pathVIRIDICgeminiR = os.path.dirname(os.path.abspath(__file__))+"/utils/geminiVIRIDIC.R"
    # ***** OUTPUT *****#
    pathTMP = geminiset.pathTMP
    pathDIRBLASTN = pathOUT+"/viridic_blastn"
    os.makedirs(pathDIRBLASTN, exist_ok=True)
    pathDIRCSV = pathOUT+"/viridic_csv"
    os.makedirs(pathDIRCSV, exist_ok=True)
    pathDIRJSON = pathOUT+"/viridic_json"
    os.makedirs(pathDIRJSON, exist_ok=True)
    pathTMPFASTA = pathTMP+"/FASTA"
    os.makedirs(pathTMPFASTA, exist_ok=True)
    pathTMPMISSING = pathTMP+"/MISSING"
    os.makedirs(pathTMPMISSING, exist_ok=True)
    printcolor("⏩ Found "+str(len(setAllOrg))+" phages"+"\n")
    # ***** CHECK MISSING DISTANCES ***** #
    printcolor("♊ Check JSON"+"\n")
    pbar = tqdm(total=len(setAllOrg), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
    cptMissing = 0
    for orgName in setAllOrg:
        if ref == "None" or orgName == ref:
            pathVIRIDICJSON = pathDIRJSON+"/"+orgName+".json"
            pathMISSINGJSON = pathTMPMISSING+"/"+orgName+".json"
            dicoMissing = {}
            if not os.path.isfile(pathVIRIDICJSON):
                dump_json({orgName: 100.0}, pathVIRIDICJSON)
                dicoMissing = list(setAllOrg - set([orgName]))
            else:
                dicoSimMA = load_json(pathVIRIDICJSON)
                dicoMissing = list(setAllOrg - set(dicoSimMA.keys()) - set([orgName]))
            cptMissing += len(dicoMissing)
            dump_json(dicoMissing, pathMISSINGJSON)
        pbar.update(1)
        title("Check JSON", pbar)
    pbar.close()
    if ref == "None":
        printcolor("⏩ Missing: "+str(cptMissing)+"/"+str((len(setAllOrg)-1)*len(setAllOrg))+" distances"+"\n")
    else:
        printcolor("⏩ Missing: "+str(cptMissing)+"/"+str(len(setAllOrg)-1)+" distances"+"\n")
    # ***** REFORMAT INPUT FASTA ***** # (only organism name in header and merge contigs with 100N)
    printcolor("♊ Reformat FASTA"+"\n")
    pbar = tqdm(total=len(setAllOrg), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
    for orgName in setAllOrg:
        if boolFromDB:
            pathInitialFASTA = pathIN+"/"+orgName+"/"+orgName+"_genomic.fna.gz"
        else:
            pathInitialFASTA = pathIN+"/"+orgName+ext
        pathReformatFASTA = pathTMPFASTA+"/"+orgName+".fasta"
        seqs = []
        if ".gz" in pathInitialFASTA:
            with gzip.open(pathInitialFASTA, 'rt') as fasta:
                prev_seq = []
                for line in fasta:
                    if line.startswith(">"):
                        seqs.append("".join(prev_seq))
                        prev_seq = []
                    else:
                        prev_seq.append(line.rstrip())
        else:
            with open(pathInitialFASTA, 'r') as fasta:
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
        pbar.update(1)
        title("Reformat", pbar)
    pbar.close()
    # ***** LAUNCH MISSING BLASTN ***** #
    printcolor("♊ Launch blastN"+"\n")
    pbar = tqdm(total=len(setAllOrg), ncols=75+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for orgName1 in setAllOrg:
        pathReformatFASTA1 = pathTMPFASTA+"/"+orgName1+".fasta"
        pathORGBLAST = pathDIRBLASTN+"/"+orgName1
        os.makedirs(pathORGBLAST, exist_ok=True)
        if boolFromDB:
            pbar.set_description_str(orgName1[:15])
        else:
            pbar.set_description_str(orgName1+" ".rjust(maxpathSize-len(orgName1)))
        # Construct blastN threads for current organism
        dicoThread = {}
        for orgName2 in setAllOrg:
            pathReformatFASTA2 = pathTMPFASTA+"/"+orgName2+".fasta"
            pathBLASTOUT = pathORGBLAST+"/"+orgName2+".out"
            cmdBLASTN = dicoGeminiPath['TOOLS']['blastn']+" -query "+pathReformatFASTA1+" -subject "+pathReformatFASTA2+" -out "+pathBLASTOUT + \
                " -outfmt \"6 qseqid sseqid evalue bitscore qlen slen qstart qend sstart send qseq sseq nident gaps\"" + \
                " -evalue 1 -max_target_seqs 50000 -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2"
            if orgName1 != orgName2 and (ref == "None" or orgName1 == ref or orgName2 == ref) and not os.path.isfile(pathBLASTOUT):
                dicoThread[orgName2] = {"cmd": cmdBLASTN, "returnstatut": None, "returnlines": []}
        # Launch blastN threads
        if len(dicoThread) > 0:
            launch_threads(dicoThread, "blastn", cpu, pathTMP)
        pbar.update(1)
        title("blastN", pbar)
    pbar.close()
    # ***** LAUNCH VIRIDIC ***** #
    if ref == "None":
        printcolor("♊ Launch VIRIDIC"+"\n")
        pbar = tqdm(total=len(setAllOrg), ncols=75+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    else:
        spinner = yaspin(Spinners.aesthetic, text="♊ Launch VIRIDIC", side="right")
        spinner.start()
    for orgName1 in setAllOrg:
        if (ref == "None" or orgName1 == ref):
            pathORGBLAST = pathDIRBLASTN+"/"+orgName1
            pathVIRIDICJSON1 = pathDIRJSON+"/"+orgName1+".json"
            dicoSimMA1 = load_json(pathVIRIDICJSON1)
            if ref == "None":
                if boolFromDB:
                    pbar.set_description_str(orgName1[:15])
                else:
                    pbar.set_description_str(orgName1+" ".rjust(maxpathSize-len(orgName1)))
            # Concatenate missing blast results
            pathMERGEBLAST = pathTMP+"/"+orgName1+".out"
            lstBlastOut = []
            pathMISSINGJSON1 = pathTMPMISSING+"/"+orgName1+".json"
            dicoMissing1 = load_json(pathMISSINGJSON1)
            for orgName2 in dicoMissing1:
                lstBlastOut.append(pathORGBLAST+"/"+orgName2+".out")
            cat_lstfiles(lstBlastOut, pathMERGEBLAST)
            # Construct VIRIDIC threads for current organism
            # Blast output could be empty
            pathTMPCSV = pathTMP+"/"+orgName1+".csv"
            try:
                os.remove(pathTMPCSV)
            except FileNotFoundError:
                pass
            cmdVIRIDIC = dicoGeminiPath['TOOLS']['rscript']+" "+pathVIRIDICgeminiR+" blastres="+pathMERGEBLAST+" out="+pathTMPCSV+" > /dev/null 2>&1"
            if os.path.getsize(pathMERGEBLAST) != 0:
                os.system(cmdVIRIDIC)
            # Parse VIRIDIC output
            # ATTENTION : if not use reciprocal blast, viridic must be sum of both comparison
            if os.path.isfile(pathTMPCSV):
                lstLines = read_file(pathTMPCSV, yaspinBool=False)[1:-1]
                for line in lstLines:
                    orgName2 = line.replace("\"", "").split("\t")[1]
                    intergSim = float(line.replace("\"", "").split("\t")[9])
                    dicoSimMA1[orgName2] = min(intergSim*2, 100.0)
                    dicoMissing1.remove(orgName2)
                    pathVIRIDICJSON2 = pathDIRJSON+"/"+orgName2+".json"
                    if os.path.isfile(pathVIRIDICJSON2):
                        dicoSimMA2 = load_json(pathVIRIDICJSON2)
                        dicoSimMA2[orgName1] = min(intergSim*2, 100.0)
                        dump_json(dicoSimMA2, pathVIRIDICJSON2)
                        pathMISSINGJSON2 = pathTMPMISSING+"/"+orgName2+".json"
                        dicoMissing2 = load_json(pathMISSINGJSON2)
                        try:
                            dicoMissing2.remove(orgName1)
                            dump_json(dicoMissing2, pathMISSINGJSON2)
                        except ValueError:
                            pass
            for orgName2 in dicoMissing1:
                dicoSimMA1[orgName2] = 0.0
                pathVIRIDICJSON2 = pathDIRJSON+"/"+orgName2+".json"
                if os.path.isfile(pathVIRIDICJSON2):
                    dicoSimMA2 = load_json(pathVIRIDICJSON2)
                    dicoSimMA2[orgName1] = 0.0
                    dump_json(dicoSimMA2, pathVIRIDICJSON2)
                    pathMISSINGJSON2 = pathTMPMISSING+"/"+orgName2+".json"
                    dicoMissing2 = load_json(pathMISSINGJSON2)
                    try:
                        dicoMissing2.remove(orgName1)
                        dump_json(dicoMissing2, pathMISSINGJSON2)
                    except ValueError:
                        pass
            dicoMissing1 = []
            dump_json(dicoMissing1, pathMISSINGJSON1)
            dump_json(dicoSimMA1, pathVIRIDICJSON1)
            if ref == "None":
                pbar.update(1)
                title("viridic", pbar)
    if ref == "None":
        pbar.close()
    else:
        spinner.stop()
        printcolor("♊ Launch VIRIDIC"+"\n")
    # ***** FAMILY, GENUS and SPECIES assignment ***** #
    if ref == "None":
        printcolor("♊ Genus/Specie assignment"+"\n")
        pbar = tqdm(total=len(setAllOrg), ncols=75+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    else:
        spinner = yaspin(Spinners.aesthetic, text="♊ Genus/Specie assignment", side="right")
        spinner.start()
    dicoTaxo = {'family': {}, 'genus': {}, 'specie': {}}
    for orgName1 in setAllOrg:
        if (ref == "None" or orgName1 == ref):
            if boolFromDB:
                pbar.set_description_str(orgName1[:15])
            else:
                pbar.set_description_str(orgName1+" ".rjust(maxpathSize-len(orgName1)))
            pathVIRIDICJSON = pathDIRJSON+"/"+orgName1+".json"
            dicoSimMA = load_json(pathVIRIDICJSON)
            for orgName2 in setAllOrg:
                for orderTuple in [("family", thfam), ("genus", thgen), ("specie", thsp)]:
                    if dicoSimMA[orgName2] >= orderTuple[1]:
                        findOrder = False
                        for orderNum in dicoTaxo[orderTuple[0]]:
                            if orgName1 in dicoTaxo[orderTuple[0]][orderNum] or orgName2 in dicoTaxo[orderTuple[0]][orderNum]:
                                dicoTaxo[orderTuple[0]][orderNum].update({orgName1, orgName2})
                                findOrder = True
                                break
                        if findOrder is False:
                            dicoTaxo[orderTuple[0]][len(dicoTaxo[orderTuple[0]])+1] = set({orgName1, orgName2})
            if ref == "None":
                pbar.update(1)
                title("assign", pbar)
    if ref == "None":
        pbar.close()
    else:
        spinner.stop()
        printcolor("♊ Genus/Specie assignment"+"\n")
    for order in ["family", "genus", "specie"]:
        if ref == "None":
            printcolor("⏩ Found "+str(len(dicoTaxo[order]))+" "+order+"\n")
        else:
            printcolor("⏩ Found "+str(len(list(dicoTaxo[order].values())[0]))+" reference "+order+" members\n")
    # ***** MAKE CLUSTERS TABLE ***** #
    printcolor("♊ Make clusters table"+"\n")
    pathOUTclusters = pathOUT+"/clusters.tsv"
    OUT = open(pathOUTclusters, 'w')
    OUT.write("Organism\tFamily\tGenus\tSpecie\n")
    dicoAssign = {}
    for orgName in setAllOrg:
        dicoAssign[orgName] = {}
        line = orgName
        for order in ["family", "genus", "specie"]:
            # Genus
            for orderNum in dicoTaxo[order]:
                if orgName in dicoTaxo[order][orderNum]:
                    line += "\t"+str(orderNum)
                    dicoAssign[orgName][order] = orderNum
                    break
        if line != orgName:
            OUT.write(line+"\n")
    OUT.close()
    # ***** PLOT SIMILARITY MATRIX ***** # (only if no reference)
    if ref == "None":
        spinner = yaspin(Spinners.aesthetic, text="♊ Plot similarity matrix", side="right")
        spinner.start()
        title("Plotting", None)
        # Convert individual VIRIDIC sqlitedict to one dataframe
        lstdicoViridic = []
        for orgName in setAllOrg:
            pathVIRIDICJSON = pathDIRJSON+"/"+orgName+".json"
            dicoSimMA = load_json(pathVIRIDICJSON)
            lstdicoViridic.append(dicoSimMA)
        df = pd.DataFrame.from_dict(lstdicoViridic)
        cmap = sns.color_palette("light:#d40000", as_cmap=True)
        cg = sns.clustermap(df, cmap=cmap, figsize=(50, 50), tree_kws={'linewidths': 2.5}, dendrogram_ratio=0.15, annot_kws={"size": 35 / np.sqrt(len(df))}, cbar_kws={'label': 'similarity %'}, linewidths=0.0, rasterized=True)
        # Retrieve ordered ticks label
        newColums = df.columns[cg.dendrogram_col.reordered_ind]
        newIndexs = df.index[cg.dendrogram_row.reordered_ind]
        newData = df.loc[newIndexs, newColums]
        orderedOrg = list(newData.keys())
        # Plot clustered heatmap
        cg.ax_cbar.tick_params(labelsize=40)
        cg.ax_cbar.yaxis.label.set_size(50)
        plt.savefig(pathOUT+"/matrix.png", dpi=300)
        plt.savefig(pathOUT+"/matrix.svg")
        spinner.stop()
        printcolor("♊ Plot similarity matrix"+"\n")
        # ***** WRITE SIMILARITY MATRIX ***** #
        printcolor("♊ Write similarity matrix"+"\n")
        pathOUTmatrix = pathOUT+"/matrix.tsv"
        OUT = open(pathOUTmatrix, 'w')
        header = "Organism\tGenus\tSpecie"
        for orgName in orderedOrg:
            header += "\t"+orgName
        OUT.write(header+"\n")
        for orgName1 in orderedOrg:
            pathVIRIDICJSON = pathDIRJSON+"/"+orgName1+".json"
            dicoSimMA = load_json(pathVIRIDICJSON)
            line = orgName1+"\t"+str(dicoAssign[orgName1]['genus'])+"\t"+str(dicoAssign[orgName1]['specie'])
            for orgName2 in orderedOrg:
                line += "\t"+str(dicoSimMA[orgName2]).replace(".", ",")
            OUT.write(line+"\n")
        OUT.close()


@fct_checker
def PhiSpy(pathIN: str, pathOUT: str, nbAdjacent: int = 3, minCtgLen: int = 5000, boolPvogs: bool = False, ext: str = ".gbk") -> Tuple[str, str, int, int, bool, str]:
    '''
     ------------------------------------------------------------
    |                           PhiSpy                           |
    |------------------------------------------------------------|
    |                PhiSpy prophages prediction                 |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN     : path of input GBK folder (required)        |
    |    pathOUT    : path of output folder (required)           |
    |    nbAdjacent : nb of consecutive genes (default=3)        |
    |    minCtgLen  : minimum contig size (default=5000)         |
    |    boolPvogs  : use PVOGS hmm database (default=False)     |
    |    ext        : extension of input files (default=.gbk)    |
     ------------------------------------------------------------
    |TOOLS: phispy                                               |
     ------------------------------------------------------------
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "filter_phage_assembly", [ext, ext+".gz"])
    if len(lstFiles) == 0:
        printcolor("[ERROR: PhiSpy]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathOUT = path_converter(pathOUT)
    os.makedirs(pathOUT, exist_ok=True)
    pathTMP = geminiset.pathTMP
    pathLOG = pathOUT+"/PhiSpy.log"
    if os.path.isfile(pathLOG):
        os.remove(pathLOG)
    if "." not in ext:
        ext = "."+ext
    # Gemini variables
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'phispy' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['phispy'])
    slurmBool, cpu, memMax, memMin = get_sys_info()
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathGBK in lstFiles:
        orgName = os.path.basename(pathGBK).replace(ext, "").replace(".gz", "")
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        pathPHISPYGBK = pathOUT+"/"+orgName+"_phage.gbk"
        pathPHISPYfasta = pathOUT+"/"+orgName+"_phage.fasta"
        if not os.path.isfile(pathPHISPYGBK) or os.path.getsize(pathPHISPYGBK) == 0 or not os.path.isfile(pathPHISPYfasta) or os.path.getsize(pathPHISPYfasta) == 0:
            if os.path.basename(pathGBK)[-3:] == ".gz":
                os.system("gzip -c -d "+pathGBK+" > "+pathTMP+"/"+orgName+".gbk")
                pathGBK = pathTMP+"/"+orgName+".gbk"
            pathTMPLOG = pathTMP+"/"+orgName+".log"
            if boolPvogs is True:
                cmdPHISPY = dicoGeminiPath['TOOLS']['phispy']+" "+pathGBK+" -o "+pathOUT+" --file_prefix "+orgName+" --output_choice 4 " + \
                    " --phmms "+dicoGeminiPath['DATABASES']['pvogs_hmm']+" --number "+str(nbAdjacent)+" --threads "+str(cpu)+" --log "+pathTMPLOG+" >> "+pathLOG+" 2>&1"
            else:
                cmdPHISPY = dicoGeminiPath['TOOLS']['phispy']+" "+pathGBK+" -o "+pathOUT+" --file_prefix "+orgName+" --output_choice 4 " + \
                    " --number "+str(nbAdjacent)+" --threads "+str(cpu)+" --log "+pathTMPLOG+" >> "+pathLOG+" 2>&1"
            os.system(cmdPHISPY)
            os.system("cat "+pathTMPLOG+" >> "+pathLOG)
            os.system("rm -f "+pathOUT+"/"+orgName+"_bacteria*")
        pbar.update(1)
        title("PhiSpy", pbar)
    pbar.close()


@fct_checker
def picmi_finder_gbk(pathIN: str, pathOUT: str, prefix: str, maxLen: int = 50000) -> Tuple[str, str, str, int]:
    '''
     ------------------------------------------------------------
    |                        PICMI FINDER                        |
    |------------------------------------------------------------|
    |                Search PICMI from gbk files                 |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN : path of input GBK file (required)              |
    |    pathOUT: path of output folder (required)               |
    |    prefix : prefix of output files (required)              |
    |    maxLen : maximum PICMI length (default=50000)           |
     ------------------------------------------------------------
    |TOOLS: diamond, hmmsearch                                   |
     ------------------------------------------------------------
    '''
    # Global
    pathTMP = geminiset.pathTMP
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'diamond' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['diamond'])
    if 'hmmsearch' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['hmmsearch'])
    slurmBool, cpu, memMax, memMin = get_sys_info()
    # Paths
    pathGBK = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    os.makedirs(pathOUT, exist_ok=True)
    os.makedirs(pathOUT+"/JSON", exist_ok=True)
    os.makedirs(pathOUT+"/FNA", exist_ok=True)
    os.makedirs(pathOUT+"/FFN", exist_ok=True)
    os.makedirs(pathOUT+"/FAA", exist_ok=True)
    os.makedirs(pathOUT+"/GBK", exist_ok=True)
    os.makedirs(pathOUT+"/GFF", exist_ok=True)
    os.makedirs(pathOUT+"/PLOTS", exist_ok=True)
    pathJSON = pathOUT+"/JSON/"+prefix+".json"
    pathFisDMND = dicoGeminiPath['DATABASES']['picmi-db']+"/fis.dmnd"
    pathHMM = dicoGeminiPath['DATABASES']['picmi-db']+"/fis-hel-alpa-int.hmm"
    # Temp
    pathFAA = pathTMP+"/"+prefix+".faa"
    pathDIAMOND = pathTMP+"/"+prefix+".tsv"
    pathHMMSEARCH = pathTMP+"/"+prefix+".hmmsearch"
    if not os.path.isfile(pathJSON):
        dicoPICMI = {}
        # Convert to FAA
        gbk_to_faa(pathIN=pathGBK, pathOUT=pathFAA)
        # Diamond blastP to find Fis
        cmd = dicoGeminiPath['TOOLS']['diamond']+" blastp -d "+pathFisDMND+" -q "+pathFAA+" -o "+pathDIAMOND+" --threads "+str(cpu)+" --max-target-seqs 1 --id 30 --outfmt 6 qtitle > /dev/null 2>&1"
        os.system(cmd)
        # If Fis found
        if os.path.isfile(pathDIAMOND):
            # Retrieve Fis sequence
            TSV = open(pathDIAMOND, 'r')
            lstLines = TSV.read().split("\n")[:-1]
            TSV.close()
            lstFisLT = []
            for line in lstLines:
                lstFisLT.append(line.split("|")[0])
            # Launch HMMSEARCH
            cmdhmmsearch = dicoGeminiPath['TOOLS']['hmmsearch']+" --cpu "+str(cpu)+" --tblout "+pathHMMSEARCH+" "+pathHMM+" "+pathFAA+" > /dev/null 2>&1"
            os.system(cmdhmmsearch)
            # If hmmsearch results
            dicoHMM = {'hel': set(), 'alpa': set(), 'int': set()}
            if os.path.isfile(pathHMMSEARCH):
                # Retrieve PICMI core genes
                TMP = open(pathHMMSEARCH, 'r')
                lstLines = TMP.read().split("\n")
                TMP.close()
                for line in lstLines:
                    if line != "" and line[0] != "#":
                        splitLine = re.compile(r"\s+").sub(" ", line).strip().split(" ")
                        lttarget = splitLine[0].split("|")[0]
                        query = splitLine[2]
                        fullSeqEvalue = float(splitLine[4])
                        if fullSeqEvalue <= 0.001:
                            if query == "DUF3987":
                                dicoHMM['hel'].add(lttarget)
                            if query == "DUF5906":
                                dicoHMM['hel'].add(lttarget)
                            if query == "Phage_AlpA":
                                dicoHMM['alpa'].add(lttarget)
                            if query == "HTH_17":
                                dicoHMM['alpa'].add(lttarget)
                            if query == "Vibrio_4_HMM_Profile":
                                dicoHMM['alpa'].add(lttarget)
                            if query == "Phage_integrase":
                                dicoHMM['int'].add(lttarget)
                            if query == "Vibrio_6_HMM_Profile":
                                dicoHMM['int'].add(lttarget)
                # Read GBK
                dicoGBK = list(make_gbk_dict(pathGBK).values())[0]
                for contig in dicoGBK:
                    # Order contig genes and tag picmi core genes
                    dicoProt = {"fis": {}, "hel": {}, "alpa": {}, "int": {}, "all": {}}
                    dicoNumToPos = {}
                    cptLT = 1
                    lstIntNum = []
                    for lt in dicoGBK[contig]['dicoLT']:
                        dicoEntry = dicoGBK[contig]['dicoLT'][lt]
                        dicoProt["all"][lt] = cptLT
                        if lt in lstFisLT:
                            dicoProt["fis"][lt] = cptLT
                        for prot in ["hel", "alpa", "int"]:
                            if lt in dicoHMM[prot] or dicoEntry['protein_id'] in dicoHMM[prot]:
                                dicoProt[prot][lt] = cptLT
                                if prot == "int":
                                    lstIntNum.append(cptLT)
                        dicoNumToPos[cptLT] = {'start': dicoEntry['start'], 'end': dicoEntry['end'], 'strand': dicoEntry['strand'], 'protSeq': dicoEntry['protSeq'], 'geneSeq': dicoEntry['geneSeq']}
                        cptLT += 1
                    cptRegion = 1
                    # Analyse picmi core genes synteny
                    if len(dicoProt["fis"]) == 0:
                        continue
                    for fisLT in dicoProt["fis"]:
                        fisNum = dicoProt["fis"][fisLT]
                        for intLT in dicoProt["int"]:
                            intNum = dicoProt["int"][intLT]
                            if (dicoNumToPos[fisNum]['strand'] > 0 and fisNum < intNum) or (dicoNumToPos[fisNum]['strand'] < 0 and fisNum > intNum):
                                for helLT in dicoProt["hel"]:
                                    helNum = dicoProt["hel"][helLT]
                                    for alpaLT in dicoProt["alpa"]:
                                        alpaNum = dicoProt["alpa"][alpaLT]
                                        if (dicoNumToPos[fisNum]['strand'] > 0 and fisNum < helNum and helNum < alpaNum and alpaNum < intNum) or (dicoNumToPos[fisNum]['strand'] < 0 and fisNum > helNum and helNum > alpaNum and alpaNum > intNum):
                                            pathPICMIFNA = pathOUT+"/FNA/"+prefix+"_____"+str(cptRegion)+".fna"
                                            pathPICMIFFN = pathOUT+"/FFN/"+prefix+"_____"+str(cptRegion)+".ffn"
                                            pathPICMIFAA = pathOUT+"/FAA/"+prefix+"_____"+str(cptRegion)+".faa"
                                            pathPICMIGBK = pathOUT+"/GBK/"+prefix+"_____"+str(cptRegion)+".gbk"
                                            pathPICMIGFF = pathOUT+"/GFF/"+prefix+"_____"+str(cptRegion)+".gff"
                                            pathPICMIPNG = pathOUT+"/PLOTS/"+prefix+"_____"+str(cptRegion)+".png"
                                            pathPICMISVG = pathOUT+"/PLOTS/"+prefix+"_____"+str(cptRegion)+".svg"
                                            # PICMI region
                                            if dicoNumToPos[fisNum]['strand'] > 0:
                                                PICMISeq = dicoGBK[contig]['seq'][dicoNumToPos[fisNum]['start']:dicoNumToPos[intNum+1]['end']]
                                            else:
                                                PICMISeq = dicoGBK[contig]['seq'][dicoNumToPos[intNum-1]['start']:dicoNumToPos[fisNum]['end']]
                                            # Filter by size:
                                            if len(PICMISeq) <= maxLen:
                                                # Write FNA
                                                OUT = open(pathPICMIFNA, 'w')
                                                if dicoNumToPos[fisNum]['strand'] > 0:
                                                    OUT.write(">PICMI ["+prefix+"] \n"+PICMISeq+"\n")
                                                else:
                                                    OUT.write(">PICMI r ["+prefix+"] \n"+reverse_complement(PICMISeq)+"\n")
                                                OUT.close()
                                                # Write FFN & FAA
                                                PICMIFFN = open(pathPICMIFFN, 'w')
                                                PICMIFAA = open(pathPICMIFAA, 'w')
                                                satelliteOthersCpt = 1
                                                if dicoNumToPos[fisNum]['strand'] > 0:
                                                    rangeLoop = range(fisNum, intNum+2, 1)
                                                else:
                                                    rangeLoop = range(fisNum, intNum-2, -1)
                                                for num in rangeLoop:
                                                    if num == fisNum:
                                                        PICMIFFN.write(">Fis ["+prefix+"]\n"+dicoNumToPos[num]['geneSeq']+"\n")
                                                        PICMIFAA.write(">Fis ["+prefix+"]\n"+dicoNumToPos[num]['protSeq']+"\n")
                                                    elif num == helNum:
                                                        PICMIFFN.write(">Primase ["+prefix+"]\n"+dicoNumToPos[num]['geneSeq']+"\n")
                                                        PICMIFAA.write(">Primase ["+prefix+"]\n"+dicoNumToPos[num]['protSeq']+"\n")
                                                    elif num == alpaNum:
                                                        PICMIFFN.write(">AlpA ["+prefix+"]\n"+dicoNumToPos[num]['geneSeq']+"\n")
                                                        PICMIFAA.write(">AlpA ["+prefix+"]\n"+dicoNumToPos[num]['protSeq']+"\n")
                                                    elif num == intNum:
                                                        PICMIFFN.write(">Integrase ["+prefix+"]\n"+dicoNumToPos[num]['geneSeq']+"\n")
                                                        PICMIFAA.write(">Integrase ["+prefix+"]\n"+dicoNumToPos[num]['protSeq']+"\n")
                                                    elif (dicoNumToPos[fisNum]['strand'] > 0 and num == intNum+1) or (dicoNumToPos[fisNum]['strand'] < 0 and num == intNum-1):
                                                        if dicoNumToPos[num]['protSeq'] is not None:
                                                            PICMIFFN.write(">Flanking ["+prefix+"]\n"+dicoNumToPos[num]['geneSeq']+"\n")
                                                            PICMIFAA.write(">Flanking ["+prefix+"]\n"+dicoNumToPos[num]['protSeq']+"\n")
                                                    else:
                                                        if dicoNumToPos[num]['protSeq'] is not None:
                                                            PICMIFFN.write(">Satellite"+str(satelliteOthersCpt).zfill(3)+" ["+prefix+"]\n"+dicoNumToPos[num]['geneSeq']+"\n")
                                                            PICMIFAA.write(">Satellite"+str(satelliteOthersCpt).zfill(3)+" ["+prefix+"]\n"+dicoNumToPos[num]['protSeq']+"\n")
                                                        satelliteOthersCpt += 1
                                                PICMIFFN.close()
                                                PICMIFAA.close()
                                                # Make GBK
                                                make_gbk_from_fasta(pathIN1=pathPICMIFNA, pathIN2=pathPICMIFFN, pathIN3=pathPICMIFAA, pathOUT=pathPICMIGBK,
                                                                    identifier="picmi", topology="linear", division="BCT", taxID=2, boolProgress=False)
                                                # Make GFF
                                                gbk_to_gff(pathIN=pathPICMIGBK, pathOUT=pathPICMIGFF)
                                                # Add repeat positions
                                                if dicoNumToPos[fisNum]['strand'] > 0:
                                                    fisSeq = dicoNumToPos[fisNum]['geneSeq']
                                                    intFlankingInterSeq = dicoGBK[contig]['seq'][dicoNumToPos[intNum]['end']:dicoNumToPos[intNum+1]['start']]
                                                else:
                                                    fisSeq = reverse_complement(dicoNumToPos[fisNum]['geneSeq'])
                                                    intFlankingInterSeq = dicoGBK[contig]['seq'][dicoNumToPos[intNum-1]['end']:dicoNumToPos[intNum]['start']]
                                                repeatSize1, i1, j1 = longest_common_substring(fisSeq, intFlankingInterSeq)
                                                repeatSeq1 = fisSeq[i1: i1 + repeatSize1]
                                                repeatSize2, i2, j2 = longest_common_substring(intFlankingInterSeq, fisSeq)
                                                repeatSeq2 = intFlankingInterSeq[i2: i2 + repeatSize2]
                                                if repeatSize1 > repeatSize2:
                                                    repeatSeq = repeatSeq1
                                                    repeatSize = repeatSize1
                                                else:
                                                    repeatSeq = repeatSeq2
                                                    repeatSize = repeatSize2
                                                if dicoNumToPos[fisNum]['strand'] > 0:
                                                    splitRepeat = PICMISeq.split(repeatSeq)
                                                else:
                                                    splitRepeat = reverse_complement(PICMISeq).split(reverse_complement(repeatSeq))
                                                firstRepeatStart = len(splitRepeat[0])+1
                                                firstRepeatEnd = firstRepeatStart+len(repeatSeq)-1
                                                secondRepeatStart = len(splitRepeat[1])+len(splitRepeat[0])+len(repeatSeq)+1
                                                secondRepeatEnd = secondRepeatStart+len(repeatSeq)-1
                                                GFF = open(pathPICMIGFF, "a")
                                                GFF.write("picmi\tGV\tCDS\t"+str(firstRepeatStart)+"\t"+str(firstRepeatEnd)+"\t.\t+\t0\tlocus_tag = repeat1;product = hypothetical protein\n")
                                                GFF.write("picmi\tGV\tCDS\t"+str(secondRepeatStart)+"\t"+str(secondRepeatEnd)+"\t.\t+\t0\tlocus_tag = repeat2;product = hypothetical protein\n")
                                                GFF.close()
                                                # Write to JSON
                                                if dicoNumToPos[fisNum]['strand'] > 0:
                                                    picmiLen = len(dicoGBK[contig]['seq'][dicoNumToPos[fisNum]['end']:dicoNumToPos[intNum]['end']])
                                                else:
                                                    picmiLen = len(dicoGBK[contig]['seq'][dicoNumToPos[intNum]['start']:dicoNumToPos[fisNum]['start']])
                                                dicoPICMI[cptRegion] = {'filename': os.path.basename(pathGBK), 'picmiLen': picmiLen,
                                                                        'repeatLen': repeatSize, 'repeatSeq': repeatSeq,
                                                                        'nbORFs': satelliteOthersCpt+3-1}
                                                # Make plot
                                                features = []
                                                dicoGFF = list(make_gff_dict(pathIN=pathPICMIGFF, ext=".gff").values())[0]
                                                startRegion = 0
                                                for geneType in dicoGFF:
                                                    if geneType != 'length':
                                                        for geneEntry in dicoGFF[geneType]:
                                                            if geneType == "CDS":
                                                                color = "#2a7fff"
                                                            elif geneType == "tRNA":
                                                                color = "#37c8ab"
                                                            else:
                                                                continue
                                                            if geneEntry['attributes']['locus_tag'] == "Fis":
                                                                color = "#ff7f2a"
                                                            if geneEntry['attributes']['locus_tag'] == "Primase":
                                                                color = "#ff5555"
                                                            if geneEntry['attributes']['locus_tag'] == "AlpA":
                                                                color = "#ffdd55"
                                                            if geneEntry['attributes']['locus_tag'] == "Integrase":
                                                                color = "#5fd35f"
                                                            if "repeat" in geneEntry['attributes']['locus_tag']:
                                                                color = "#000000"
                                                            geneFeature = GraphicFeature(start=geneEntry['start'], end=geneEntry['end'], strand=int(geneEntry['strand']+"1"), color=color, linewidth=0)
                                                            features.append(geneFeature)
                                                record = GraphicRecord(sequence_length=maxLen, features=features, first_index=startRegion-100)
                                                ax, _ = record.plot(figure_width=50)
                                                ax.figure.savefig(pathPICMIPNG, dpi=300)
                                                ax.figure.savefig(pathPICMISVG)
                                                plt.close('all')
                                                # IF many region for a same organism
                                                cptRegion += 1
                                                dump_json(dicoPICMI, pathJSON)


@fct_checker
def picmi_finder_databankseq(pathIN: str, pathOUT: str, maxLen: int = 50000) -> Tuple[str, str, int]:
    '''
     ------------------------------------------------------------
    |                        PICMI FINDER                        |
    |------------------------------------------------------------|
    |                Search PICMI from gbk files                 |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN : path of input databank SEQ file (required)     |
    |    pathOUT: path of output folder (required)               |
    |    maxLen : maximum PICMI length (default=50000)           |
     ------------------------------------------------------------
    |TOOLS: diamond, hmmsearch                                   |
     ------------------------------------------------------------
    '''
    # Global
    pathTMP = geminiset.pathTMP
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'diamond' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['diamond'])
    if 'hmmsearch' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['hmmsearch'])
    slurmBool, cpu, memMax, memMin = get_sys_info()
    # Paths
    pathSEQ = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    os.makedirs(pathOUT, exist_ok=True)
    os.makedirs(pathOUT+"/JSON", exist_ok=True)
    os.makedirs(pathOUT+"/FNA", exist_ok=True)
    os.makedirs(pathOUT+"/FFN", exist_ok=True)
    os.makedirs(pathOUT+"/FAA", exist_ok=True)
    os.makedirs(pathOUT+"/GBK", exist_ok=True)
    os.makedirs(pathOUT+"/GFF", exist_ok=True)
    os.makedirs(pathOUT+"/PLOTS", exist_ok=True)
    pathFisDMND = dicoGeminiPath['DATABASES']['picmi-db']+"/fis.dmnd"
    pathHMM = dicoGeminiPath['DATABASES']['picmi-db']+"/fis-hel-alpa-int.hmm"
    pathLSTdone = pathOUT+"/lst_done.json"
    if os.path.isfile(pathLSTdone):
        lstDone = load_json(pathLSTdone)
    else:
        lstDone = []

    # ***** READ DATABANK SEQ FILE ***** #
    spinner = yaspin(Spinners.aesthetic, text="♊ Read SEQ file", side="right")
    spinner.start()
    title("Read SEQ", None)
    dicoGBK = list(make_gbk_dict(pathSEQ).values())[0]
    spinner.stop()
    printcolor("♊ Read SEQ file"+"\n")
    pbar = tqdm(total=len(dicoGBK), ncols=75, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for record in dicoGBK:
        pbar.set_description_str(record)
        orgName = dicoGBK[record]['dicoSource']['orgName']+"_"+record
        if orgName not in lstDone:
            dicoPICMI = {}
            # Extract protein sequences
            pathFAA = pathTMP+"/"+record+".faa"
            FAA = open(pathFAA, 'w')
            cpt_prot = 1
            for lt in dicoGBK[record]['dicoLT']:
                try:
                    FAA.write(">"+lt+"\n"+dicoGBK[record]['dicoLT'][lt]['protSeq']+"\n")
                    cpt_prot += 1
                except TypeError:
                    pass
            FAA.close()
            pathDIAMOND = pathTMP+"/"+record+".tsv"
            # Diamond blastP to find Fis
            cmd = dicoGeminiPath['TOOLS']['diamond']+" blastp -d "+pathFisDMND+" -q "+pathFAA+" -o "+pathDIAMOND+" --threads "+str(cpu)+" --max-target-seqs 1 --id 30 --outfmt 6 qtitle > /dev/null 2>&1"
            os.system(cmd)
            # If Fis found
            if os.path.isfile(pathDIAMOND):
                # Retrieve Fis sequence
                TSV = open(pathDIAMOND, 'r')
                lstLines = TSV.read().split("\n")[:-1]
                TSV.close()
                lstFisLT = []
                for line in lstLines:
                    lstFisLT.append(line.split("|")[0])
                # Launch HMMSEARCH
                pathHMMSEARCH = pathTMP+"/"+record+".hmmsearch"
                cmdhmmsearch = dicoGeminiPath['TOOLS']['hmmsearch']+" --cpu "+str(cpu)+" --tblout "+pathHMMSEARCH+" "+pathHMM+" "+pathFAA+" > /dev/null 2>&1"
                os.system(cmdhmmsearch)
                # If hmmsearch results
                dicoHMM = {'hel': set(), 'alpa': set(), 'int': set()}
                if os.path.isfile(pathHMMSEARCH):
                    # Retrieve PICMI core genes
                    TMP = open(pathHMMSEARCH, 'r')
                    lstLines = TMP.read().split("\n")
                    TMP.close()
                    for line in lstLines:
                        if line != "" and line[0] != "#":
                            splitLine = re.compile(r"\s+").sub(" ", line).strip().split(" ")
                            lttarget = splitLine[0].split("|")[0]
                            query = splitLine[2]
                            fullSeqEvalue = float(splitLine[4])
                            if fullSeqEvalue <= 0.001:
                                if query == "DUF3987":
                                    dicoHMM['hel'].add(lttarget)
                                if query == "DUF5906":
                                    dicoHMM['hel'].add(lttarget)
                                if query == "Phage_AlpA":
                                    dicoHMM['alpa'].add(lttarget)
                                if query == "Vibrio_4_HMM_Profile":
                                    dicoHMM['alpa'].add(lttarget)
                                if query == "HTH_17":
                                    dicoHMM['alpa'].add(lttarget)
                                if query == "Phage_integrase":
                                    dicoHMM['int'].add(lttarget)
                                if query == "Vibrio_6_HMM_Profile":
                                    dicoHMM['int'].add(lttarget)
                        # Order contig genes and tag picmi core genes
                        dicoProt = {"fis": {}, "hel": {}, "alpa": {}, "int": {}, "all": {}}
                        dicoNumToPos = {}
                        cptLT = 1
                        lstIntNum = []
                        for lt in dicoGBK[record]['dicoLT']:
                            dicoEntry = dicoGBK[record]['dicoLT'][lt]
                            dicoProt["all"][lt] = cptLT
                            if lt in lstFisLT:
                                dicoProt["fis"][lt] = cptLT
                            for prot in ["hel", "alpa", "int"]:
                                if lt in dicoHMM[prot] or dicoEntry['protein_id'] in dicoHMM[prot]:
                                    dicoProt[prot][lt] = cptLT
                                    if prot == "int":
                                        lstIntNum.append(cptLT)
                            dicoNumToPos[cptLT] = {'start': dicoEntry['start'], 'end': dicoEntry['end'], 'strand': dicoEntry['strand'], 'protSeq': dicoEntry['protSeq'], 'geneSeq': dicoEntry['geneSeq']}
                            cptLT += 1
                        cptRegion = 1
                        # Analyse picmi core genes synteny
                        if len(dicoProt["fis"]) == 0:
                            continue
                        for fisLT in dicoProt["fis"]:
                            fisNum = dicoProt["fis"][fisLT]
                            for intLT in dicoProt["int"]:
                                intNum = dicoProt["int"][intLT]
                                if (dicoNumToPos[fisNum]['strand'] > 0 and fisNum < intNum) or (dicoNumToPos[fisNum]['strand'] < 0 and fisNum > intNum):
                                    for helLT in dicoProt["hel"]:
                                        helNum = dicoProt["hel"][helLT]
                                        for alpaLT in dicoProt["alpa"]:
                                            alpaNum = dicoProt["alpa"][alpaLT]
                                            if (dicoNumToPos[fisNum]['strand'] > 0 and fisNum < helNum and helNum < alpaNum and alpaNum < intNum) or (dicoNumToPos[fisNum]['strand'] < 0 and fisNum > helNum and helNum > alpaNum and alpaNum > intNum):
                                                pathPICMIFNA = pathOUT+"/FNA/"+orgName+"_____"+str(cptRegion)+".fna"
                                                pathPICMIFFN = pathOUT+"/FFN/"+orgName+"_____"+str(cptRegion)+".ffn"
                                                pathPICMIFAA = pathOUT+"/FAA/"+orgName+"_____"+str(cptRegion)+".faa"
                                                pathPICMIGBK = pathOUT+"/GBK/"+orgName+"_____"+str(cptRegion)+".gbk"
                                                pathPICMIGFF = pathOUT+"/GFF/"+orgName+"_____"+str(cptRegion)+".gff"
                                                pathPICMIPNG = pathOUT+"/PLOTS/"+orgName+"_____"+str(cptRegion)+".png"
                                                pathPICMISVG = pathOUT+"/PLOTS/"+orgName+"_____"+str(cptRegion)+".svg"
                                                pathPICMIJSON = pathOUT+"/JSON/"+orgName+"_____"+str(cptRegion)+".json"
                                                # PICMI region
                                                if dicoNumToPos[fisNum]['strand'] > 0:
                                                    PICMISeq = dicoGBK[record]['seq'][dicoNumToPos[fisNum]['start']:dicoNumToPos[intNum+1]['end']]
                                                else:
                                                    PICMISeq = dicoGBK[record]['seq'][dicoNumToPos[intNum-1]['start']:dicoNumToPos[fisNum]['end']]
                                                # Filter by size:
                                                if len(PICMISeq) <= maxLen:
                                                    # Write FNA
                                                    OUT = open(pathPICMIFNA, 'w')
                                                    if dicoNumToPos[fisNum]['strand'] > 0:
                                                        OUT.write(">PICMI ["+orgName+"] \n"+PICMISeq+"\n")
                                                    else:
                                                        OUT.write(">PICMI r ["+orgName+"] \n"+reverse_complement(PICMISeq)+"\n")
                                                    OUT.close()
                                                    # Write FFN & FAA
                                                    PICMIFFN = open(pathPICMIFFN, 'w')
                                                    PICMIFAA = open(pathPICMIFAA, 'w')
                                                    satelliteOthersCpt = 1
                                                    if dicoNumToPos[fisNum]['strand'] > 0:
                                                        rangeLoop = range(fisNum, intNum+2, 1)
                                                    else:
                                                        rangeLoop = range(fisNum, intNum-2, -1)
                                                    for num in rangeLoop:
                                                        if num == fisNum:
                                                            PICMIFFN.write(">Fis ["+orgName+"]\n"+dicoNumToPos[num]['geneSeq']+"\n")
                                                            PICMIFAA.write(">Fis ["+orgName+"]\n"+dicoNumToPos[num]['protSeq']+"\n")
                                                        elif num == helNum:
                                                            PICMIFFN.write(">Primase ["+orgName+"]\n"+dicoNumToPos[num]['geneSeq']+"\n")
                                                            PICMIFAA.write(">Primase ["+orgName+"]\n"+dicoNumToPos[num]['protSeq']+"\n")
                                                        elif num == alpaNum:
                                                            PICMIFFN.write(">AlpA ["+orgName+"]\n"+dicoNumToPos[num]['geneSeq']+"\n")
                                                            PICMIFAA.write(">AlpA ["+orgName+"]\n"+dicoNumToPos[num]['protSeq']+"\n")
                                                        elif num == intNum:
                                                            PICMIFFN.write(">Integrase ["+orgName+"]\n"+dicoNumToPos[num]['geneSeq']+"\n")
                                                            PICMIFAA.write(">Integrase ["+orgName+"]\n"+dicoNumToPos[num]['protSeq']+"\n")
                                                        elif num == intNum+1:
                                                            if dicoNumToPos[num]['protSeq'] is not None:
                                                                PICMIFFN.write(">Flanking ["+orgName+"]\n"+dicoNumToPos[num]['geneSeq']+"\n")
                                                                PICMIFAA.write(">Flanking ["+orgName+"]\n"+dicoNumToPos[num]['protSeq']+"\n")
                                                        else:
                                                            if dicoNumToPos[num]['protSeq'] is not None:
                                                                PICMIFFN.write(">Satellite"+str(satelliteOthersCpt).zfill(3)+" ["+orgName+"]\n"+dicoNumToPos[num]['geneSeq']+"\n")
                                                                PICMIFAA.write(">Satellite"+str(satelliteOthersCpt).zfill(3)+" ["+orgName+"]\n"+dicoNumToPos[num]['protSeq']+"\n")
                                                            satelliteOthersCpt += 1
                                                    PICMIFFN.close()
                                                    PICMIFAA.close()
                                                    # Make GBK
                                                    make_gbk_from_fasta(pathIN1=pathPICMIFNA, pathIN2=pathPICMIFFN, pathIN3=pathPICMIFAA, pathOUT=pathPICMIGBK,
                                                                        identifier="picmi", topology="linear", division="BCT", taxID=2, boolProgress=False)
                                                    # Make GFF
                                                    gbk_to_gff(pathIN=pathPICMIGBK, pathOUT=pathPICMIGFF)
                                                    # Add repeat positions
                                                    if dicoNumToPos[fisNum]['strand'] > 0:
                                                        fisSeq = dicoNumToPos[fisNum]['geneSeq']
                                                        intFlankingInterSeq = dicoGBK[record]['seq'][dicoNumToPos[intNum]['end']:dicoNumToPos[intNum+1]['start']]
                                                    else:
                                                        fisSeq = reverse_complement(dicoNumToPos[fisNum]['geneSeq'])
                                                        intFlankingInterSeq = dicoGBK[record]['seq'][dicoNumToPos[intNum-1]['end']:dicoNumToPos[intNum]['start']]
                                                    repeatSize1, i1, j1 = longest_common_substring(fisSeq, intFlankingInterSeq)
                                                    repeatSeq1 = fisSeq[i1: i1 + repeatSize1]
                                                    repeatSize2, i2, j2 = longest_common_substring(intFlankingInterSeq, fisSeq)
                                                    repeatSeq2 = intFlankingInterSeq[i2: i2 + repeatSize2]
                                                    if repeatSize1 > repeatSize2:
                                                        repeatSeq = repeatSeq1
                                                        repeatSize = repeatSize1
                                                    else:
                                                        repeatSeq = repeatSeq2
                                                        repeatSize = repeatSize2
                                                    if dicoNumToPos[fisNum]['strand'] > 0:
                                                        splitRepeat = PICMISeq.split(repeatSeq)
                                                    else:
                                                        splitRepeat = reverse_complement(PICMISeq).split(reverse_complement(repeatSeq))
                                                    firstRepeatStart = len(splitRepeat[0])+1
                                                    firstRepeatEnd = firstRepeatStart+len(repeatSeq)-1
                                                    secondRepeatStart = len(splitRepeat[1])+len(splitRepeat[0])+len(repeatSeq)+1
                                                    secondRepeatEnd = secondRepeatStart+len(repeatSeq)-1
                                                    GFF = open(pathPICMIGFF, "a")
                                                    GFF.write("picmi\tGV\tCDS\t"+str(firstRepeatStart)+"\t"+str(firstRepeatEnd)+"\t.\t+\t0\tlocus_tag = repeat1;product = hypothetical protein\n")
                                                    GFF.write("picmi\tGV\tCDS\t"+str(secondRepeatStart)+"\t"+str(secondRepeatEnd)+"\t.\t+\t0\tlocus_tag = repeat2;product = hypothetical protein\n")
                                                    GFF.close()
                                                    # Write to JSON
                                                    if dicoNumToPos[fisNum]['strand'] > 0:
                                                        picmiLen = len(dicoGBK[record]['seq'][dicoNumToPos[fisNum]['end']:dicoNumToPos[intNum]['end']])
                                                    else:
                                                        picmiLen = len(dicoGBK[record]['seq'][dicoNumToPos[intNum]['start']:dicoNumToPos[fisNum]['start']])
                                                    dicoPICMI[cptRegion] = {'orgName': dicoGBK[record]['dicoSource']['orgName'],
                                                                            'picmiLen': picmiLen, 'repeatLen': repeatSize,
                                                                            'repeatSeq': repeatSeq, 'nbORFs': satelliteOthersCpt+3-1}
                                                    # Make plot
                                                    features = []
                                                    dicoGFF = list(make_gff_dict(pathIN=pathPICMIGFF, ext=".gff").values())[0]
                                                    startRegion = 0
                                                    for geneType in dicoGFF:
                                                        if geneType != 'length':
                                                            for geneEntry in dicoGFF[geneType]:
                                                                if geneType == "CDS":
                                                                    color = "#2a7fff"
                                                                elif geneType == "tRNA":
                                                                    color = "#37c8ab"
                                                                else:
                                                                    continue
                                                                if geneEntry['attributes']['locus_tag'] == "Fis":
                                                                    color = "#ff7f2a"
                                                                if geneEntry['attributes']['locus_tag'] == "Primase":
                                                                    color = "#ff5555"
                                                                if geneEntry['attributes']['locus_tag'] == "AlpA":
                                                                    color = "#ffdd55"
                                                                if geneEntry['attributes']['locus_tag'] == "Integrase":
                                                                    color = "#5fd35f"
                                                                if "repeat" in geneEntry['attributes']['locus_tag']:
                                                                    color = "#000000"
                                                                geneFeature = GraphicFeature(start=geneEntry['start'], end=geneEntry['end'], strand=int(geneEntry['strand']+"1"), color=color, linewidth=0)
                                                                features.append(geneFeature)
                                                    graphRecord = GraphicRecord(sequence_length=maxLen, features=features, first_index=startRegion-100)
                                                    ax, _ = graphRecord.plot(figure_width=50)
                                                    ax.figure.savefig(pathPICMIPNG, dpi=300)
                                                    ax.figure.savefig(pathPICMISVG)
                                                    plt.close('all')
                                                    # Dump JSON
                                                    dump_json(dicoPICMI, pathPICMIJSON)
                                                    # IF many region for a same organism
                                                    cptRegion += 1
            lstDone.append(orgName)
        pbar.update(1)
        title("picmi_finder", pbar)
    pbar.close()
    # save record done
    dump_json(lstDone, pathOUT+"/lst_done.json")
