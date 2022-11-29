'''
┌─┐  ┌─┐  ┌┬┐  ┬  ┌┐┌  ┬
│ ┬  ├┤   │││  │  │││  │
└─┘  └─┘  ┴ ┴  ┴  ┘└┘  ┴ assembly
---------------------------------------------------
-*- coding: utf-8 -*-                              |
title            : geminiassembly.py               |
description      : gemini assembly functions       |
author           : dooguypapua                     |
lastmodification : 20220805                        |
version          : 0.1                             |
python_version   : 3.8.5                           |
---------------------------------------------------
'''

import os
import sys
import shutil
import geminiset
from tabulate import tabulate
from typing import Tuple
from tqdm import tqdm
from geminini import path_converter, get_input_files, printcolor, get_gemini_path, get_sys_info, fct_checker
from geminini import title, exit_gemini
from geminiparse import unwrap_fasta, make_blast_dict, make_fasta_dict, check_circular


@fct_checker
def phage_assembly(pathIN: str, pathOUT: str) -> Tuple[str, str]:
    '''
     ------------------------------------------------------------
    |                       PHAGE ASSEMBLY                       |
     ------------------------------------------------------------
    |    Trim FASTQ & de novo assembly using spades (-careful)   |
     ------------------------------------------------------------
    |PARAMETERS                                                  |
    |    pathIN  : path of input files or folder (required)      |
    |    pathOUT : path of output files (required)               |
     ------------------------------------------------------------
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "phage_assembly", [".fastq.gz", ".fastq", ".fq.gz", ".fq"])
    if len(lstFiles) == 0:
        printcolor("[ERROR: phage_assembly]\nAny input files found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    dicoGeminiPath = get_gemini_path()
    slurmBool, cpu, memMax, memMin = get_sys_info()
    if pathOUT == "":
        printcolor("[ERROR: phage_assembly]\nMissing '-o'pathOUT\n", 1, "212;64;89", "None", True)
        exit_gemini()
    pathOUT = path_converter(pathOUT)
    os.makedirs(pathOUT, exist_ok=True)
    trimmomatic_param = "LEADING: 3 TRAILING: 3 SLIDINGWINDOW: 4:15 MINLEN: 36"
    spades_param = "--careful --cov-cutoff auto -k 21, 33, 55, 77 -m 10"
    pbar = tqdm(total=int(len(lstFiles)/2), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathFASTQ in lstFiles:
        file = os.path.basename(pathFASTQ)
        if "_R1" in file:
            sampleName = os.path.basename(pathFASTQ).split("_R1")[0]
            pathFASTA = pathOUT+"/"+sampleName+".fasta"
            if not os.path.isfile(pathFASTA):
                # Construct requierd path
                pathFASTQr1 = pathFASTQ
                pathFASTQr2 = pathFASTQ.replace("_R1", "_R2")
                if not os.path.isfile(pathFASTQr2):
                    printcolor("[ERROR: phage_assembly]\nUnable to find R2 for \""+pathFASTQr1+"\"\n", 1, "212;64;89", "None", True)
                    exit_gemini()
                pathTRIMr1 = geminiset.pathTMP+"/"+os.path.basename(pathFASTQr1).replace("_R1", "_trimmed_R1")
                pathTRIMr2 = geminiset.pathTMP+"/"+os.path.basename(pathFASTQr2).replace("_R2", "_trimmed_R2")
                pathUNTRIMr1 = geminiset.pathTMP+"/"+os.path.basename(pathFASTQr1).replace("_R1", "_untrimmed_R1")
                pathUNTRIMr2 = geminiset.pathTMP+"/"+os.path.basename(pathFASTQr2).replace("_R2", "_untrimmed_R2")
                pathUNTRIM = geminiset.pathTMP+"/"+os.path.basename(pathFASTQr1).replace("_R1", "_untrimmed")
                pathSPADES = geminiset.pathTMP+"/spades_"+sampleName
                pathSPADESlog = geminiset.pathTMP+"/spades_"+sampleName+".log"
                # Launch command
                cmdTrim = "java -jar "+dicoGeminiPath['trimmomatic']+" PE -threads "+str(cpu)+" -quiet "+pathFASTQr1+" "+pathFASTQr2+" "+pathTRIMr1+" "+pathUNTRIMr1+" "+pathTRIMr2+" "+pathUNTRIMr2+" "+trimmomatic_param
                pbar.set_description_str(sampleName+": trimmomatic"+" ".rjust(maxpathSize-len(sampleName)))
                os.system(cmdTrim)
                cmdCat = "cat "+pathUNTRIMr1+" "+pathUNTRIMr2+" > "+pathUNTRIM
                os.system(cmdCat)
                cmdSpades = dicoGeminiPath['spades']+" -1 "+pathTRIMr1+" -2 "+pathTRIMr2+" -s "+pathUNTRIM+" -o "+pathSPADES+" --threads "+str(cpu)+" "+spades_param+" > "+pathSPADESlog
                pbar.set_description_str(sampleName+": spades"+" ".rjust(maxpathSize-len(sampleName)+5))
                os.system(cmdSpades)
                # Copy output contigs FASTA file
                shutil.copyfile(pathSPADES+"/contigs.fasta", pathFASTA)
            # Unwrap FASTA
            unwrap_fasta(pathFASTA)
            pbar.update(1)
            title("Phage ass", pbar)
    pbar.close()


@fct_checker
def filter_phage_assembly(pathIN: str, pathOUT: str, minLen: int = 100, minCov: int = 10, ext: str = ".fasta") -> Tuple[str, str, int, int, str]:
    '''
     ------------------------------------------------------------
    |                     FILTER PHAGE ASSEMBLY                  |
     ------------------------------------------------------------
    |               Filter phage assembly contigs                |
     ------------------------------------------------------------
    |PARAMETERS                                                  |
    |    pathIN : path of input files or folder (required)       |
    |    pathOUT: path of output files (required)                |
    |    minLen : minimum contig length (default=100)            |
    |    minCov : minimum contig coverage (default=10)           |
    |    ext    : extension of input files (default=.fasta)      |
     ------------------------------------------------------------
    '''
    pathOUT = path_converter(pathOUT)
    lstFiles, maxpathSize = get_input_files(pathIN, "filter_phage_assembly", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: filter_phage_assembly]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    os.makedirs(pathOUT, exist_ok=True)
    tabulateTable = []
    dicoGeminiPath = get_gemini_path()
    printcolor("♊ Filter assembly"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathAssembly in lstFiles:
        file = os.path.basename(pathAssembly)
        orgName = file.replace(ext, "").replace("."+ext, "")
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        # Blast contigs against contamination sequence
        cmdBlastn = dicoGeminiPath['blastn']+" -task blastn -query "+pathAssembly+" -db "+dicoGeminiPath['contaPhageDB']+" -out "+geminiset.pathTMP+"/blast.xml -outfmt 5 -perc_identity 80 -qcov_hsp_perc 30 -max_hsps 10 -max_target_seqs 5"
        os.system(cmdBlastn)
        dicoBlast = make_blast_dict(pathIN=geminiset.pathTMP+"/blast.xml", outfmt=".xml", idThr=30, minLRthr=0, maxLRthr=0, ext=".xml")
        lstContamContig = list(dicoBlast.keys())
        # Read contigs FASTA
        dicoAssembly = make_fasta_dict(pathAssembly)
        dicoFiltered = {}
        totalSize = 0
        for contig in dicoAssembly:
            contigCov = float(contig.split("_cov_")[1])
            if contigCov >= minCov and len(dicoAssembly[contig]) >= minLen and contig not in lstContamContig:
                oneNuclContig = False
                for nucl in ["A", "T", "G", "C"]:
                    if dicoAssembly[contig].count(nucl) == len(dicoAssembly[contig]):
                        oneNuclContig = True
                        break
                if not oneNuclContig:
                    overlap = check_circular(dicoAssembly[contig])
                    if overlap == "":
                        dicoFiltered["CONTIG_"+str(len(dicoFiltered)+1)+"|cov: "+str(contigCov)] = dicoAssembly[contig]
                    else:
                        dicoFiltered["CONTIG_"+str(len(dicoFiltered)+1)+"|cov: "+str(contigCov)+"|circular: "+overlap] = dicoAssembly[contig]
                    totalSize += len(dicoAssembly[contig])
        # Write new filtered FASTA
        pathOUTassembly = pathOUT+"/"+file
        OUT = open(pathOUTassembly, 'w')
        for key in dicoFiltered:
            OUT.write(">"+key+"\n"+dicoFiltered[key]+"\n")
        OUT.close()
        # Results table
        tabulateTable.append([orgName, str(len(dicoAssembly)), str(len(dicoFiltered)), overlap != "", str(totalSize)])
        pbar.update(1)
        title("Filter assembly", pbar)
    pbar.close()
    printcolor("♊ Results\n")
    printcolor(tabulate(tabulateTable, headers=["orgName", "Initial contigs", "Final contigs", "Circular", "Final size"], tablefmt="pretty")+"\n")


@fct_checker
def contigs_extension(pathIN: str, pathOUT: str, minLen: int = 100, minCov: int = 10, ext: str = ".fasta") -> Tuple[str, str, int, int, str]:
    '''
     ------------------------------------------------------------
    |                     CONTIGS EXTENSION                      |
     ------------------------------------------------------------
    |     Extend contigs using Illumina reads and references     |
     ------------------------------------------------------------
    |PARAMETERS                                                  |
    |    pathIN : path of input files or folder (required)       |
    |    pathOUT: path of output files (required)                |
    |    ext    : extension of input files (default=.fasta)      |
     ------------------------------------------------------------
    '''
    pathOUT = path_converter(pathOUT)
    lstFiles, maxpathSize = get_input_files(pathIN, "filter_phage_assembly", [ext])

    start = "TACCAATTAAAGTAAAAATATGATCTAATTGAGGCCGTCGAACTCTTACA"
    end = "GGTTTCTGTTGAAGCCATCCCGTGTCAGCGAGGTGGCATTATAGAGATTT"
    refSeq = ""
    fastq = "/mnt/c/Users/dgoudenege/Desktop/Clade1/FASTQ_others_clade1/Vibrio_crassostreae_31_O_70_strain_3960_R1.fastq.gz"
    concatSeq = start
    while end not in concatSeq:
        print(str(len(concatSeq))+"/"+str(len(refSeq)+100))
        cmdZGREP = "zgrep -Eo \""+concatSeq[-50:]+".{50}\" "+fastq+" | cut -d \":\" -f 2 | sort | uniq -c > /tmp/zgrep.txt"
        os.system(cmdZGREP)
        OUT = open("/tmp/zgrep.txt", 'r')
        lstLines = OUT.read().split("\n")[:-1]
        OUT.close()
        lstGrep = []
        setCount = set()
        for line in lstLines:
            splitLine = line.strip().split(" ")
            count = int(splitLine[0])
            seq = splitLine[1]
            if count >= 10:
                setCount.add(count)
                lstGrep.append((count, seq))
        # If only one concat directly
        if len(lstGrep) == 1:
            rightSeq = lstGrep[0][1][-50:]
            concatSeq += rightSeq
        else:
            # Analyse sequences
            lstInRef = []
            for grepRes in lstGrep:
                findInRef = False
                findInRefSamePos = False
                if grepRes[1] in refSeq:
                    findInRef = True
                    if refSeq.find(grepRes[1], len(concatSeq)-100)+100 == len(concatSeq):
                        findInRefSamePos = True
                if findInRef and findInRefSamePos:
                    lstInRef.append(grepRes)
            # If found in sequence concat
            if len(lstInRef) == 1:
                rightSeq = lstInRef[0][1][-50:]
                concatSeq += rightSeq
            # Else display results
            else:
                dicoChoose = {}
                cpt = 1
                for count in sorted(setCount)[::-1]:
                    for grepRes in lstGrep:
                        if grepRes[0] == count:
                            findInRef = False
                            findInRefSamePos = False
                            if grepRes[1] in refSeq:
                                findInRef = True
                                if refSeq.find(grepRes[1], len(concatSeq)-100)+100 == len(concatSeq):
                                    findInRefSamePos = True
                            dicoChoose[cpt] = grepRes
                            print(str(grepRes[0]).rjust(3)+": "+grepRes[1]+" ["+str(cpt).zfill(2)+"]"+" "+str(findInRef)+" "+str(findInRefSamePos))
                            cpt += 1
                choose = input("Selected seq (or q): ")
                if choose == "q":
                    print("\n"+concatSeq)
                    exit()
                else:
                    rightSeq = dicoChoose[int(choose)][1][-50:]
                    concatSeq += rightSeq
    print(concatSeq)
    exit()
