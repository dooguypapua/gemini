'''
┌─┐  ┌─┐  ┌┬┐  ┬  ┌┐┌  ┬
│ ┬  ├┤   │││  │  │││  │
└─┘  └─┘  ┴ ┴  ┴  ┘└┘  ┴ parse
---------------------------------------------------
-*- coding: utf-8 -*-                              |
title            : py                              |
description      : gemini parser functions         |
author           : dooguypapua                     |
lastmodification : 20210628                        |
version          : 0.1                             |
python_version   : 3.8.5                           |
---------------------------------------------------
'''

import os
import re
import sys
import shutil
import gzip
import pandas as pd
import geminiset
from typing import Tuple
from tqdm import tqdm
from datetime import date
from collections import OrderedDict
import xlsxwriter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from ete3 import NCBITaxa
from geminini import path_converter, printcolor, load_json, dump_json, fct_checker, get_input_files
from geminini import get_gemini_path, reverse_complement, title, exit_gemini, read_file, cat_lstfiles


'''
-------------------------------------------------------------------------------------------------------
                                        MANIPULATION FUNCTIONS
-------------------------------------------------------------------------------------------------------
'''


def get_refseqacc_from_dict(lstDico: list) -> Tuple[list]:
    '''
     ------------------------------------------------------------
    |            GET REFSEQ ACCESSION FROM BLAST DICT            |
    |------------------------------------------------------------|
    |   Get refseq accession from a list of blast dictionnary    |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    lstDico     : list of blast dictionnary (required)      |
    |RETURN                                                      |
    |    setRefseqAcc: set of refseq accession                   |
     ------------------------------------------------------------
    '''
    setRefseqAcc = set()
    for dicoBLAST in lstDico:
        for orgName in dicoBLAST:
            for lt in dicoBLAST[orgName]:
                for hit in dicoBLAST[orgName][lt]:
                    refseqAcc = hit.split("|")[0]
                    setRefseqAcc.add(refseqAcc)
    return setRefseqAcc


def make_uniprotmapping_dict(pathIN: str, setRefseqAcc: set = "None", pathJSON: str = "None") -> Tuple[str, set, str]:
    '''
     ------------------------------------------------------------
    |                   PARSE UNIPROT MAPPING                    |
    |------------------------------------------------------------|
    |       Parse Uniprot mapping and create a dictionnary       |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of TAB uniprot file (required)           |
    |    pathJSON: path of output JSON ("")                      |
    |RETURN                                                      |
    |    dicoUNIPROT: {refseqAcc: 'geneName','lstECnumber'}      |
     ------------------------------------------------------------
    '''
    printcolor("♊ Uniprot mapping"+"\n")
    if pathJSON != "None" and os.path.isfile(pathJSON):
        dicoUniprotMapping = load_json(pathJSON)
    else:
        dicoUniprotMapping = {}
        lstLines = read_file(pathIN)
        pbar = tqdm(total=len(lstLines), ncols=75, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
        for line in lstLines:
            splitLine = line.split("\t")
            geneName = splitLine[2].split(" ")[0]
            splitECnumber = splitLine[3].split("; ")
            splitRefseqAcc = splitLine[4].split(";")
            for refseqAcc in splitRefseqAcc:
                if setRefseqAcc != "None" and refseqAcc in setRefseqAcc:
                    dicoUniprotMapping[refseqAcc.split(" ")[0]] = {'geneName': geneName, 'lstECnumber': splitECnumber}
            pbar.update(1)
        pbar.close()
        if pathJSON != "None":
            dump_json(dicoUniprotMapping, pathJSON)
    return dicoUniprotMapping


# ***** Reformat hmmscanDict to key==LT ***** #
def hmmscandict_to_dictlt(dicoHMMSCAN):
    dicoLT = {}
    for orgName in dicoHMMSCAN:
        dicoLT[orgName] = {}
        for entry in dicoHMMSCAN[orgName]:
            for lt in dicoHMMSCAN[orgName][entry]:
                if lt not in dicoLT[orgName]:
                    dicoLT[orgName][lt] = set()
                dicoLT[orgName][lt].add(entry)
        for lt in dicoLT[orgName]:
            dicoLT[orgName][lt] = list(dicoLT[orgName][lt])
    return dicoLT


'''
-------------------------------------------------------------------------------------------------------
                                             FASTA FUNCTIONS
-------------------------------------------------------------------------------------------------------
'''


@fct_checker
def make_fasta_dict(pathIN: str, onlyLTasHeader: bool = False, pathJSON: str = "None") -> Tuple[str, bool, str]:
    '''
     ------------------------------------------------------------
    |                   MAKE FASTA DICTIONNARY                   |
     ------------------------------------------------------------
    |          Parse FASTA file and create a dictionnary         |
     ------------------------------------------------------------
    |PARAMETERS                                                  |
    |    pathIN         : path of input FASTA file (required)    |
    |    onlyLTasHeader : only keep LT as header (False)         |
    |    pathJSON       : path of output JSON ("")               |
    |RETURN                                                      |
    |    dicoFASTA: {id: sequence}                             |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathJSON = path_converter(pathJSON)
    if not os.path.isfile(pathIN):
        printcolor("[ERROR: make_fasta_dict]\nInput file \""+pathIN+"\" not found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if pathJSON != "None" and os.path.isfile(pathJSON):
        dicoFASTA = load_json(pathJSON)
    else:
        dicoFASTA = {}
        lstRecordObj = []
        if os.path.basename(pathIN)[-3:] == ".gz":
            with gzip.open(pathIN, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    lstRecordObj.append(record)
        else:
            for record in SeqIO.parse(open(pathIN, "r"), "fasta"):
                lstRecordObj.append(record)
        for record in lstRecordObj:
            if onlyLTasHeader is True:
                if "|" in record.description:
                    dicoFASTA[record.description.split("|")[0]] = str(record.seq).replace("*", "")
                elif "[" in record.description:
                    dicoFASTA[record.description.split(" [")[0]] = str(record.seq).replace("*", "")
                else:
                    dicoFASTA[record.description] = str(record.seq).replace("*", "")
            else:
                dicoFASTA[record.description] = str(record.seq).replace("*", "")
        if pathJSON != "None":
            dump_json(dicoFASTA, pathJSON)
    return dicoFASTA


@fct_checker
def search_in_fasta(pathIN: str, searchTerm: str, pathOUT: str = "", boolDoublon: bool = False, ext: str = ".faa") -> Tuple[str, str, str, bool, str]:
    '''
     ------------------------------------------------------------
    |                   GET SEQUENCE WITH TERM                   |
     ------------------------------------------------------------
    |   Search term in FASTA folder or file and return sequence  |
     ------------------------------------------------------------
    | PARAMETERS                                                 |
    |    pathIN      : path of input FASTA file (required)       |
    |    searchTerm  : term to search (delim comma) (required)   |
    |    pathOUT     : path of output FASTA file (default="")    |
    |    boolDoublon : remove sequence doublon (default=False)   |
    |    ext         : extension of input files (default=.faa)   |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    if pathOUT != "":
        pathOUT = path_converter(pathOUT)
    lstFiles, maxpathSize = get_input_files(pathIN, "search_in_fasta", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: search_in_fasta]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    found = False
    splitSearchTerm = searchTerm.split(",")
    dicoSearch = {}
    printcolor("♊ Search in FASTA"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    setSeq = set()
    for pathFASTA in lstFiles:
        filename = os.path.basename(pathFASTA)
        orgName = filename.replace(ext, "")
        dicoFASTA = make_fasta_dict(pathFASTA)
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        for key in dicoFASTA:
            for term in splitSearchTerm:
                if term.lower() in key.lower():
                    if boolDoublon is False or (boolDoublon is True and dicoFASTA[key] not in setSeq):
                        try:
                            dicoSearch[term].append(">"+key+"\n"+dicoFASTA[key])
                        except KeyError:
                            dicoSearch[term] = [">"+key+"\n"+dicoFASTA[key]]
                    else:
                        setSeq.add(dicoFASTA[key])
                    found = True
        pbar.update()
    pbar.close()
    if found is True:
        if pathOUT != "":
            OUT = open(pathOUT, 'w')
        for term in sorted(dicoSearch.keys()):
            for seqFind in dicoSearch[term]:
                if pathOUT == "":
                    printcolor(seqFind+"\n")
                else:
                    OUT.write(seqFind+"\n")
        if pathOUT != "":
            OUT.close()
    else:
        printcolor("Any sequence found with term(s) \""+searchTerm+"\""+"\n")


@fct_checker
def unwrap_fasta(pathIN: str, ext: str = ".fna") -> Tuple[str, str]:
    '''
     ------------------------------------------------------------
    |                   UNWRAP FASTA SEQUENCES                   |
    |------------------------------------------------------------|
    |           Unwrap FASTA by removing sequence '\n'           |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN: path of input FASTA file or folder (required)   |
    |    ext   : extension of input files (default=.fna)         |
     ------------------------------------------------------------
    '''
    lstFiles, maxpathSize = get_input_files(pathIN, "unwrap_fasta", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: unwrap_fasta]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    for pathFASTA in lstFiles:
        dicoFASTA = make_fasta_dict(pathFASTA)
        OUT = open(pathFASTA, 'w')
        for key in dicoFASTA:
            OUT.write(">"+key+"\n"+dicoFASTA[key]+"\n")
        OUT.close()


@fct_checker
def split_fasta(pathIN: str, pathOUT: str, term: str = "N", ext: str = ".fna") -> Tuple[str, str, str, str]:
    '''
     ------------------------------------------------------------
    |                      SPLIT SEQUENCES                       |
    |------------------------------------------------------------|
    |                    Split FASTA sequence                    |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN : path of input FASTA file or folder (required)  |
    |    pathOUT: path of output FASTA file or folder (required) |
    |    term   : split term or character (default=N)            |
    |    ext    : extension of input files (default=.fna)        |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    if os.path.isfile(pathIN):
        ext = pathIN.split(".")[-1]
    lstFiles, maxpathSize = get_input_files(pathIN, "split_fasta", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: split_fasta]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if os.path.isdir(pathIN):
        os.makedirs(pathOUT, exist_ok=True)
    for pathFASTA in lstFiles:
        # Load input FASTA
        dicoFASTA = make_fasta_dict(pathFASTA)
        if os.path.isdir(pathIN):
            pathOUTFASTA = pathOUT+"/"+os.path.basename(pathFASTA)
        else:
            pathOUTFASTA = pathOUT
        # Find max split
        maxSplit = 0
        for key in dicoFASTA:
            if term in dicoFASTA[key]:
                maxSplit = max(maxSplit, len(list(filter(None, dicoFASTA[key].split(term)))))
        # Split & write
        OUT = open(pathOUTFASTA, 'w')
        for key in dicoFASTA:
            if term not in dicoFASTA[key]:
                OUT.write(">"+key+"\n"+dicoFASTA[key]+"\n")
            else:
                cpt = 1
                splitSeq = dicoFASTA[key].split(term)
                for subSeq in list(filter(None, splitSeq)):
                    OUT.write(">"+key+" (split"+str(cpt).zfill(len(str(maxSplit)))+")\n"+subSeq+"\n")
                    cpt += 1
        OUT.close()


@fct_checker
def check_circular(seqIN: str, minLen: int = 3) -> Tuple[str, int]:
    '''
     ------------------------------------------------------------
    |                  CHECK CIRCULAR SEQUENCES                  |
    |------------------------------------------------------------|
    |              Check if a sequence is circular               |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN: path of input FASTA file or folder (required)   |
    |    minLen: min repetition lentgh (default=3)               |
    |RETURN                                                      |
    |    Overlap begin-end sequence                              |
     ------------------------------------------------------------
    '''
    overlap = ""
    for i in range(int(len(seqIN)/2)):
        overlap += seqIN[i]
        # Stop if just find one time
        if seqIN.count(overlap) == 1:
            break
    search = re.search("^"+overlap[: -1]+".+"+overlap[: -1]+"$", seqIN)
    if search and len(overlap)-1 >= minLen:
        return overlap[: -1]
    else:
        return ""


@fct_checker
def fasta_stats(pathIN: str, pathOUT: str = "None", boolSort: bool = True, ext: str = ".fna") -> Tuple[str, str, bool, str]:
    '''
     ------------------------------------------------------------
    |                   MAKE FASTA STATS SIZE                    |
    |------------------------------------------------------------|
    |          Create a table with various FASTA stats           |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input files or folder (required)      |
    |    pathOUT : path of output CSV (default=None [stdout])    |
    |    boolSort: sort output by filename (default=True)        |
    |    ext     : extension of input files (default=.fna)       |
     ------------------------------------------------------------
    '''
    dicoSize = {}
    dicoNbCtg = {}
    if pathOUT != "None":
        pathOUT = path_converter(pathOUT)
    lstFiles, maxpathSize = get_input_files(pathIN, "fasta_genome_size", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: fasta_genome_size]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    for pathFNA in lstFiles:
        file = os.path.basename(pathFNA)
        org = file.replace(ext, "")
        dicoSize[org] = 0
        dicoNbCtg[org] = 0
        for record in SeqIO.parse(open(pathFNA, "r"), "fasta"):
            dicoSize[org] += len(record.seq)
            dicoNbCtg[org] += 1
    if boolSort is True:
        dicoSize = dict(sorted(dicoSize.items(), key=lambda item: item[0]))
        dicoNbCtg = dict(sorted(dicoNbCtg.items(), key=lambda item: item[0]))
    df = pd.DataFrame({'Org': list(dicoNbCtg.keys()), 'NbCtg': list(dicoNbCtg.values()), 'Size': list(dicoSize.values())})
    if pathOUT == "None":
        printcolor(df.to_string(index=False)+"\n")
    else:
        df.to_csv(pathOUT, sep='\t', encoding='utf-8', index=False)


'''
-------------------------------------------------------------------------------------------------------
                                            GENBANK FUNCTIONS
-------------------------------------------------------------------------------------------------------
'''


@fct_checker
def make_gbk_dict(pathIN: str, pathJSON: str = "None", boolSort: bool = True, boolPseudo: bool = False) -> Tuple[str, str, bool, bool]:
    '''
     ------------------------------------------------------------
    |                  MAKE GENBANK DICTIONNARY                  |
    |------------------------------------------------------------|
    |          Parse GBK files and create a dictionnary          |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input files or folder (required)      |
    |    pathJSON: path of output JSON ("")                      |
    |    boolSort: sort output by filename (default=True)        |
    |RETURN                                                      |
    |    dicoGBK: {org: dicoFeatures}                          |
    |    dicoFeatures: {locustag: dicoFeature}                  |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathJSON = path_converter(pathJSON)
    dicoGBK = {}
    lstExcludeGBKtype = ["assembly_gap", "misc_feature", "regulatory", "repeat_region", "gap", "misc_binding", "source"]
    lstFiles, maxpathSize = get_input_files(pathIN, "make_gbk_dict", [".gbff", ".gb", ".gbk", ".gbff.gz", ".gb.gz", ".gbk.gz", ".seq"])
    if len(lstFiles) == 0:
        printcolor("[ERROR: make_gbk_dict]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if pathJSON != "None" and os.path.isfile(pathJSON):
        dicoGBK = load_json(pathJSON)
    else:
        for pathGBK in lstFiles:
            file = os.path.basename(pathGBK)
            orgName = file.replace(".gbff", "").replace(".gbk", "").replace(".gb", "").replace(".gz", "").replace("_genomic", "")
            dicoGBK[orgName] = {}
            lstRecordObj = []
            # Retrieve list of record objects
            if ".gz" in file:
                with gzip.open(pathGBK, "rt") as handle:
                    try:
                        for record in SeqIO.parse(handle, "gb"):
                            lstRecordObj.append(record)
                    except ValueError:
                        pass
            else:
                try:
                    for record in SeqIO.parse(pathGBK, "gb"):
                        lstRecordObj.append(record)
                except ValueError:
                    pass
            # Browse records
            cptLT = 1
            for record in lstRecordObj:
                dicoGBK[orgName][record.id] = {'seq': str(record.seq), 'dicoSource': {}, 'dicoLT': OrderedDict(), 'description': record.description, 'annotations': record.annotations}
                for feature in record.features:
                    if feature.type == "source":
                        for field in ["organism", "strain", "serotype"]:
                            if field in feature.qualifiers:
                                dicoGBK[orgName][record.id]["dicoSource"][field] = feature.qualifiers[field][0]
                            else:
                                dicoGBK[orgName][record.id]["dicoSource"][field] = None
                        # Format orgName
                        try:
                            formatorgName = dicoGBK[orgName][record.id]['dicoSource']['organism'].replace(" ", "_")
                        except AttributeError:
                            formatorgName = orgName
                        strain = dicoGBK[orgName][record.id]['dicoSource']['strain']
                        serotype = dicoGBK[orgName][record.id]['dicoSource']['serotype']
                        if strain is not None and strain not in formatorgName:
                            formatorgName += "_"+strain
                        if serotype is not None and serotype not in formatorgName:
                            formatorgName += "_"+serotype
                        formatorgName = formatorgName.replace(" ", "_").replace("/", "_").replace(" = ", "_").replace("(", "_").replace(")", "").replace(":", "_").replace("[", "_").replace("]", "").replace(";", "").replace("\'", "")
                        dicoGBK[orgName][record.id]["dicoSource"]["orgName"] = formatorgName
                        # DBXREF
                        if "db_xref" in feature.qualifiers:
                            for db_xref in feature.qualifiers["db_xref"]:
                                name = db_xref.split(":")[0]
                                value = db_xref.split(":")[1]
                                if name == "taxon":
                                    dicoGBK[orgName][record.id]["dicoSource"]['taxon'] = int(value)
                        if "host" in feature.qualifiers:
                            dicoGBK[orgName][record.id]["host"] = feature.qualifiers["host"][0]
                        if "lab_host" in feature.qualifiers:
                            dicoGBK[orgName][record.id]["lab_host"] = feature.qualifiers["lab_host"][0]                            
                    elif feature.type not in lstExcludeGBKtype and (boolPseudo is True or 'pseudo' not in feature.qualifiers):
                        try:
                            locustag = feature.qualifiers["locus_tag"][0]
                        except KeyError:
                            locustag = "gene_"+str(cptLT).zfill(4)
                        dicoGBK[orgName][record.id]['dicoLT'][locustag] = {}
                        if 'pseudo' in feature.qualifiers:
                            dicoGBK[orgName][record.id]['dicoLT'][locustag]['type'] = feature.type+"_pseudo"
                        else:
                            dicoGBK[orgName][record.id]['dicoLT'][locustag]['type'] = feature.type
                        try:
                            dicoGBK[orgName][record.id]['dicoLT'][locustag]['product'] = feature.qualifiers["product"][0].replace(" ", "#").replace(";", "#").replace("|", "#").replace(", ", "#")
                        except KeyError:
                            dicoGBK[orgName][record.id]['dicoLT'][locustag]['product'] = None
                        try:
                            dicoGBK[orgName][record.id]['dicoLT'][locustag]['protSeq'] = str(feature.qualifiers["translation"][0])
                        except KeyError:
                            dicoGBK[orgName][record.id]['dicoLT'][locustag]['protSeq'] = None
                        dicoGBK[orgName][record.id]['dicoLT'][locustag]['start'] = feature.location.nofuzzy_start
                        dicoGBK[orgName][record.id]['dicoLT'][locustag]['end'] = feature.location.nofuzzy_end
                        dicoGBK[orgName][record.id]['dicoLT'][locustag]['strand'] = feature.strand
                        try:
                            dicoGBK[orgName][record.id]['dicoLT'][locustag]['protein_id'] = feature.qualifiers["protein_id"][0]
                        except KeyError:
                            dicoGBK[orgName][record.id]['dicoLT'][locustag]['protein_id'] = None
                        if feature.strand == 1:
                            dicoGBK[orgName][record.id]['dicoLT'][locustag]['geneSeq'] = str(record.seq[feature.location.nofuzzy_start: feature.location.nofuzzy_end])
                        else:
                            dicoGBK[orgName][record.id]['dicoLT'][locustag]['geneSeq'] = str(record.seq[feature.location.nofuzzy_start: feature.location.nofuzzy_end].reverse_complement())
                        cptLT += 1
        if boolSort is True:
            if pathJSON != "None":
                dump_json(dict(sorted(dicoGBK.items(), key=lambda item: item[0])), pathJSON)
        else:
            if pathJSON != "None":
                dump_json(dicoGBK, pathJSON)
    return dict(sorted(dicoGBK.items(), key=lambda item: item[0]))
    return dicoGBK


@fct_checker
def make_gff_dict(pathIN: str, pathJSON: str = "None", ext: str = ".gff") -> Tuple[str, str, str]:
    '''
     ------------------------------------------------------------
    |                   MAKE GFF3 DICTIONNARY                    |
    |------------------------------------------------------------|
    |          Parse GFF3 files and create a dictionnary         |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input files or folder (required)      |
    |    pathJSON: path of output JSON ("")                      |
    |    boolSort: sort output by filename (default=True)        |
    |RETURN                                                      |
    |    dicoGFF: {org: dicoFeatures}                          |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathJSON = path_converter(pathJSON)
    dicoGFF = {}
    lstFiles, maxpathSize = get_input_files(pathIN, "make_gff_dict", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: make_gff_dict]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if pathJSON != "None" and os.path.isfile(pathJSON):
        dicoGFF = load_json(pathJSON)
    else:
        for pathGFF in lstFiles:
            file = os.path.basename(pathGFF)
            orgName = file.replace(ext, "")
            dicoGFF[orgName] = {}
            GFF = open(pathGFF, 'r')
            lstLines = GFF.read().split("\n")
            GFF.close()
            for line in lstLines:
                if line != "":
                    splitLine = line.split("\t")
                    # retrieve header
                    if line[0] == "#":
                        # Get length
                        if "sequence-region" in splitLine[0]:
                            try:
                                dicoGFF[orgName]['length'] = int(splitLine[0].split(" ")[3])
                            except IndexError:
                                dicoGFF[orgName]['length'] = int(splitLine[2])
                        elif "Sequence Data" in splitLine[0]:
                            splitHeader = splitLine[0].split(";")
                            for header in splitHeader:
                                if "seqlen" in header:
                                    dicoGFF[orgName]['length'] = int(header.split(" = ")[1])
                    else:
                        seqType = splitLine[2]
                        attributes = splitLine[8]
                        # retrieve attributes
                        dicoAttributes = {}
                        for field in attributes.split(";"):
                            if "=" in field:
                                dicoAttributes[field.split("=")[0]] = field.split("=")[1].replace("\"", "")
                        # create entry subdictionnary
                        dicoEntry = {
                                    'seqID': splitLine[0], 'source': splitLine[1],
                                    'start': int(splitLine[3]), 'end': int(splitLine[4]),
                                    'score': splitLine[5], 'strand': splitLine[6],
                                    'phase': splitLine[7],
                                    'attributes': dicoAttributes
                                   }
                        # Add to dictionnary
                        try:
                            dicoGFF[orgName][seqType].append(dicoEntry)
                        except KeyError:
                            dicoGFF[orgName][seqType] = [dicoEntry]
        if pathJSON != "None":
            dump_json(dicoGFF, pathJSON)
    return dicoGFF


@fct_checker
def gbk_to_faa(pathIN: str, pathOUT: str, syntaxic: str = "prodigal", boolSplit: bool = False) -> Tuple[str, str, str, bool]:
    '''
     ------------------------------------------------------------
    |                   CONVERT GENBANK TO FAA                   |
    |------------------------------------------------------------|
    |  Retrieve protein sequence in a GBK file and create a FAA  |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN   : path of input GBK file (required)            |
    |    pathOUT  : path of output FAA file (required)           |
    |    syntaxic : syntaxic tool for missing (default=prodigal) |
    |    boolSplit: create one output per contig (default=False) |
     ------------------------------------------------------------
    |TOOLS: prodigal                                             |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    dicoGBK = make_gbk_dict(pathIN)
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'prodigal' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['prodigal'])
    org = list(dicoGBK.keys())[0]
    if boolSplit is False:
        OUT = open(pathOUT, 'w')
    for contig in dicoGBK[org]:
        if boolSplit is True and os.path.isfile(pathOUT+"/"+contig+".faa"):
            os.makedirs(pathOUT, exist_ok=True)
            os.remove(pathOUT+"/"+contig+".faa")
        # Case of any syntaxic annotation
        if len(dicoGBK[org][contig]['dicoLT']) == 0:
            TMPFASTA = open(geminiset.pathTMP+"/gbk_to_faa_temp.fna", 'w')
            TMPFASTA.write(">"+contig+"\n"+dicoGBK[org][contig]['seq'])
            TMPFASTA.close()
            if syntaxic == "prodigal":
                cmdSYNTAXIC = dicoGeminiPath['TOOLS']['prodigal']+" -d "+geminiset.pathTMP+"/gbk_to_faa_temp.ffn -g 11 -i "+geminiset.pathTMP+"/gbk_to_faa_temp.fna -q > /dev/null 2>&1"
                os.system(cmdSYNTAXIC)
                dicoFFN = make_fasta_dict(geminiset.pathTMP+"/gbk_to_faa_temp.ffn")
                FAA = open(geminiset.pathTMP+"/gbk_to_faa_temp.faa", 'w')
                cpt = 1
                for key in dicoFFN:
                    cdna = Seq(dicoFFN[key])
                    seqProt = str(cdna.translate(table=11, cds=False, to_stop=True))
                    FAA.write(">"+key+"\n"+seqProt+"\n")
                    cpt += 1
                FAA.close()
                dicoFAA = make_fasta_dict(geminiset.pathTMP+"/gbk_to_faa_temp.faa")
                os.remove(geminiset.pathTMP+"/gbk_to_faa_temp.ffn")
                os.remove(geminiset.pathTMP+"/gbk_to_faa_temp.faa")
                for header in dicoFAA:
                    splitHeader = header.split(" # ")
                    dicoGBK[org][contig]['dicoLT'][splitHeader[0]] = {'type': "CDS", 'product': None, 'start': int(splitHeader[1]), 'end': int(splitHeader[2]),
                                                                      'strand': int(splitHeader[3]), 'geneSeq': None, 'protSeq': dicoFAA[header]}
        for lt in dicoGBK[org][contig]['dicoLT']:
            if dicoGBK[org][contig]['dicoLT'][lt]['protSeq'] is not None:
                toWrite = ">"+lt+"|"+contig+"|"+str(dicoGBK[org][contig]['dicoLT'][lt]['product'])+"|"+org+"\n"+dicoGBK[org][contig]['dicoLT'][lt]['protSeq']+"\n"
                if boolSplit is True:
                    OUT = open(pathOUT+"/"+contig+".faa", 'a')
                OUT.write(toWrite)
        if boolSplit is True and os.path.isfile(pathOUT+"/"+contig+".faa"):
            OUT.close()
    if boolSplit is False:
        OUT.close()


@fct_checker
def gbk_to_ffn(pathIN: str, pathOUT: str, syntaxic: str = "prodigal") -> Tuple[str, str, str]:
    '''
     ------------------------------------------------------------
    |                   CONVERT GENBANK TO FFN                   |
    |------------------------------------------------------------|
    |   Retrieve gene sequence in a GBK file and create a FFN    |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input GBK file (required)             |
    |    pathOUT : path of output FFN file (required)            |
    |    syntaxic: syntaxic missing annot tool (default=prodigal)|
     ------------------------------------------------------------
    |TOOLS: prodigal                                             |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    dicoGBK = make_gbk_dict(pathIN)
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'prodigal' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['prodigal'])
    org = list(dicoGBK.keys())[0]
    OUT = open(pathOUT, 'w')
    for contig in dicoGBK[org]:
        # Case of any syntaxic annotation
        if len(dicoGBK[org][contig]['dicoLT']) == 0:
            TMPFASTA = open(geminiset.pathTMP+"/gbk_to_ffn_temp.fna", 'w')
            TMPFASTA.write(">"+contig+"\n"+dicoGBK[org][contig]['seq'])
            TMPFASTA.close()
            if syntaxic == "prodigal":
                cmdSYNTAXIC = dicoGeminiPath['TOOLS']['prodigal']+" -d "+geminiset.pathTMP+"/gbk_to_ffn_temp.ffn -g 11 -i "+geminiset.pathTMP+"/gbk_to_ffn_temp.fna -q > /dev/null 2>&1"
                os.system(cmdSYNTAXIC)
                dicoFFN = make_fasta_dict(geminiset.pathTMP+"/gbk_to_ffn_temp.ffn")
                for header in dicoFFN:
                    splitHeader = header.split(" # ")
                    dicoGBK[org][contig]['dicoLT'][splitHeader[0]] = {'type': "CDS", 'product': None, 'start': int(splitHeader[1]), 'end': int(splitHeader[2]),
                                                                      'strand': int(splitHeader[3]), 'geneSeq': dicoFFN[header], 'protSeq': None}
        for lt in dicoGBK[org][contig]['dicoLT']:
            if dicoGBK[org][contig]['dicoLT'][lt]['type'] in ["CDS", "tRNA", "rRNA"]:
                toWrite = ">"+lt+"|"+str(dicoGBK[org][contig]['dicoLT'][lt]['product'])+"|"+org+"\n"+dicoGBK[org][contig]['dicoLT'][lt]['geneSeq']+"\n"
                OUT.write(toWrite)
    OUT.close()


@fct_checker
def gbk_to_fna(pathIN: str, pathOUT: str) -> Tuple[str, str]:
    '''
     ------------------------------------------------------------
    |                   CONVERT GENBANK TO FNA                   |
    |------------------------------------------------------------|
    |   Retrieve genome sequence in a GBK file and create a FNA  |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN : path of input GBK file (required)              |
    |    pathOUT: path of output FNA file (required)             |
     ------------------------------------------------------------
     '''
    pathIN = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    name = os.path.basename(pathIN).replace("_genomic.gbff.gz", "").replace(".gbk.gz", "")
    if not os.path.isfile(pathOUT):
        dicoGBK = list(make_gbk_dict(pathIN).values())[0]
        # WRITE FNA
        OUT = open(pathOUT, 'w')
        for contig in dicoGBK:
            toWrite = ">"+contig+" ["+name+"]\n"+dicoGBK[contig]['seq']+"\n"
            OUT.write(toWrite)
        OUT.close()


@fct_checker
def gbk_to_annotFAA(pathIN: str, pathOUT: str) -> Tuple[str, str]:
    '''
     ------------------------------------------------------------
    |              CONVERT GENBANK TO ANNOTATED FAA              |
    |------------------------------------------------------------|
    |  Get protein in GBK files and create FAAs with annotation  |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN : path of input GBK file (required)              |
    |    pathOUT: path of output FNA file (required)             |
     ------------------------------------------------------------
     '''
    pathIN = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    lstFiles, maxpathSize = get_input_files(pathIN, "gbk_to_annotFAA", [".gbff", ".gb", ".gbk", ".gbff.gz", ".gb.gz", ".gbk.gz", ".seq"])
    if len(lstFiles) == 0:
        printcolor("[ERROR: gbk_to_annotFAA]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    os.makedirs(pathOUT, exist_ok=True)
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathGBK in lstFiles:
        file = os.path.basename(pathGBK)
        for ext in [".gbff", ".gb", ".gbk", ".gbff.gz", ".gb.gz", ".gbk.gz", ".seq"]:
            orgName = file.replace(ext, "")
        pathFAA = pathOUT+"/"+orgName+".faa"
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        if not os.path.isfile(pathFAA):
            FAA = open(pathFAA, 'w')
            if file[-3:] == ".gz":
                with gzip.open(pathGBK, "rt") as handle:
                    lstRecords = handle.read().replace("://", ":").replace("//\n ", ":").split("//\n")
            else:
                IN = open(pathGBK, 'r')
                lstRecords = IN.read().replace("://", ":").replace("//\n ", ":").split("//\n")
                IN.close()
            for record in lstRecords:
                if "CDS" in record:
                    splitCDS = record.split("FEATURES")[1].split("     CDS             ")
                    for cds in splitCDS[1:]:
                        dataCDS = cds.replace("\n", "").replace("                     ", "")
                        searchProt = re.search("product=\"([^\"]+)\"/protein_id=\"([^\"]+)\".*/translation=\"([^\"]+)\"", dataCDS)
                        if searchProt:
                            FAA.write(">"+searchProt.group(2)+"|"+searchProt.group(1)+"\n"+searchProt.group(3)+"\n")
                        else:
                            searchProt2 = re.search("protein_id=\"([^\"]+)\".*/translation=\"([^\"]+)\"", dataCDS)
                            if searchProt2:
                                FAA.write(">"+searchProt2.group(1)+"|hypothetical protein\n"+searchProt2.group(2)+"\n")
                            else:
                                if "/pseudo" not in dataCDS and "GCA_" in file:
                                    printcolor("[ERROR: gbk_to_annotFAA]\nParsing \""+file+"\"\n", 1, "212;64;89", "None", True)
                                    printcolor(dataCDS+"\n", 1, "212;64;89", "None", True)
                                    exit_gemini()
            FAA.close()
        pbar.update(1)
        title("extract", pbar)
    pbar.close()


@fct_checker
def gbk_to_all(pathIN: str, pathOUT: str, syntaxic: str = "prodigal") -> Tuple[str, str, str]:
    '''
     ------------------------------------------------------------
    |                   CONVERT GENBANK TO ALL                   |
    |------------------------------------------------------------|
    |   Retrieve genome sequence in a GBK file and create all    |
    |------------------------------------------------------------|
    | PARAMETERS                                                 |
    |    pathIN : path of input GBK file (required)              |
    |    pathOUT: path of output folder (required)               |
    |    syntaxic: syntaxic missing annot tool (default=prodigal)|
     ------------------------------------------------------------
    |TOOLS: prodigal                                             |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    if 'prodigal' in dicoGeminiModule:
        os.system("module load "+dicoGeminiModule['prodigal'])
    os.makedirs(pathOUT, exist_ok=True)
    # Construct paths
    dicoGBK = list(make_gbk_dict(pathIN).values())[0]
    filename = os.path.basename(pathIN)
    orgName = list(dicoGBK.values())[0]['dicoSource']['orgName']
    if "GCA" in filename:
        orgName = orgName+"__GCA"+filename.split("_")[1]
    if ".gz" in filename[-3:]:
        pathGBKtmp = geminiset.pathTMP+"/"+orgName+".gbk.gz"
    else:
        pathGBKtmp = geminiset.pathTMP+"/"+orgName+".gbk"
    pathFNA = pathOUT+"/"+orgName+".fna"
    pathFFN = pathOUT+"/"+orgName+".ffn"
    pathFAA = pathOUT+"/"+orgName+".faa"
    if not os.path.isfile(pathFNA) or os.path.getsize(pathFNA) == 0 or not os.path.isfile(pathFFN) or os.path.getsize(pathFFN) == 0 or not os.path.isfile(pathFAA) or os.path.getsize(pathFAA) == 0:
        shutil.copyfile(pathIN, pathGBKtmp)
        FNA = open(pathFNA, 'w')
        FFN = open(pathFFN, 'w')
        FAA = open(pathFAA, 'w')
        pathTMP = geminiset.pathTMP
        # Launch conversion
        for contig in dicoGBK:
            FNA.write(">"+contig+" ["+orgName+"]\n"+dicoGBK[contig]['seq']+"\n")
            # Case of any syntaxic annotation
            if len(dicoGBK[contig]['dicoLT']) == 0:
                TMPFASTA = open(pathTMP+"/"+orgName+"_tmp.fna", 'w')
                TMPFASTA.write(">"+contig+"\n"+dicoGBK[contig]['seq'])
                TMPFASTA.close()
                if syntaxic == "prodigal":
                    cmdSYNTAXIC = dicoGeminiPath['TOOLS']['prodigal']+" -d "+pathTMP+"/"+orgName+"_tmp.ffn -g 11 -i "+pathTMP+"/"+orgName+"_tmp.fna -q > /dev/null 2>&1"
                    os.system(cmdSYNTAXIC)
                    exit()
                    dicoFFN = make_fasta_dict(pathTMP+"/"+orgName+"_tmp.ffn")
                    TMPFAA = open(pathTMP+"/"+orgName+"_tmp.faa", 'w')
                    for key in dicoFFN:
                        seqProt = str(Seq(dicoFFN[key]).translate(table=11, cds=False, to_stop=True))
                        TMPFAA.write(">"+key+"\n"+seqProt+"\n")
                    TMPFAA.close()
                    dicoFAA = make_fasta_dict(pathTMP+"/"+orgName+"_tmp.faa")
                    for header in dicoFFN:
                        splitHeader = header.split(" # ")
                        dicoGBK[contig]['dicoLT'][splitHeader[0]] = {'type': "CDS", 'product': None, 'start': int(splitHeader[1]), 'end': int(splitHeader[2]),
                                                                     'strand': int(splitHeader[3]), 'geneSeq': dicoFFN[header], 'protSeq': dicoFAA[header]}
            # Browse locus_tag
            for lt in dicoGBK[contig]['dicoLT']:
                if dicoGBK[contig]['dicoLT'][lt]['type'] in ["CDS", "tRNA", "rRNA"]:
                    toWrite = ">"+lt+"|"+str(dicoGBK[contig]['dicoLT'][lt]['product'])+" ["+orgName+"]\n"+dicoGBK[contig]['dicoLT'][lt]['geneSeq']+"\n"
                    FFN.write(toWrite)
                    if dicoGBK[contig]['dicoLT'][lt]['type'] == "CDS":
                        toWrite = ">"+lt+"|"+str(dicoGBK[contig]['dicoLT'][lt]['product'])+" ["+orgName+"]\n"+dicoGBK[contig]['dicoLT'][lt]['protSeq']+"\n"
                        FAA.write(toWrite)
        FNA.close()
        FFN.close()
        FAA.close()


@fct_checker
def gbk_to_gff(pathIN: str, pathOUT: str) -> Tuple[str, str]:
    '''
     ------------------------------------------------------------
    |                   CONVERT GENBANK TO GFF                   |
    |------------------------------------------------------------|
    |   Retrieve genome sequence in a GBK file and create a GFF  |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN : path of input GBK file (required)              |
    |    pathOUT: path of output GFF file (required)             |
     ------------------------------------------------------------
     '''
    pathIN = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    dicoGBK = make_gbk_dict(pathIN)
    org = list(dicoGBK.keys())[0]
    OUT = open(pathOUT, 'w')
    OUT.write("##gff-version 3\n")
    for contig in dicoGBK[org]:
        OUT.write("##sequence-region "+contig+" 1 "+str(len(dicoGBK[org][contig]['seq']))+"\n")
        for lt in dicoGBK[org][contig]['dicoLT']:
            if dicoGBK[org][contig]['dicoLT'][lt]['strand'] == 1:
                frame = "+"
            else:
                frame = "-"
            if dicoGBK[org][contig]['dicoLT'][lt]['product'] is None:
                product = "hypothetical protein"
            else:
                product = dicoGBK[org][contig]['dicoLT'][lt]['product']
            line = contig+"\tGV\tCDS\t"+str(dicoGBK[org][contig]['dicoLT'][lt]['start'])+"\t"+str(dicoGBK[org][contig]['dicoLT'][lt]['end']) + \
                "\t.\t"+frame+"\t0\tlocus_tag="+lt+";product="+product+"\n"
            OUT.write(line)
    OUT.close()


@fct_checker
def gff_to_table(pathIN: str, pathOUT: str, format: str = ".xlsx", maxWidth: int = 50, ext: str = ".gff") -> Tuple[str, str, str, int, str]:
    '''
     ------------------------------------------------------------
    |           GFF to annotation table (.tsv or .xlsx)          |
    |------------------------------------------------------------|
    |          Create a annotation table from GFF files          |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input GFF file (required)             |
    |    pathOUT : path of output file (required)                |
    |    format  : output format .xlsx or .tsv (default=.xlsx)   |
    |    maxWidth: maximum column width (default=50)             |
    |    ext     : extension of input files (default=.gff)       |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathOUT = path_converter(pathOUT)
    lstFiles, maxpathSize = get_input_files(pathIN, "gff_to_table", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: gff_to_table]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    # Check output format
    outFormat = format.replace(".", "").lower()
    if outFormat not in ["tsv", "xlsx", "xls"]:
        printcolor("[ERROR: gff_to_table]\nInvalid output format, must be \".tsv\" or \".xlsx\"\n", 1, "212;64;89", "None", True)
        exit_gemini()
    printcolor("♊ GFF to table"+"\n")
    pbar = tqdm(total=len(lstFiles), ncols=50+maxpathSize, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
    for pathGFF in lstFiles:
        file = os.path.basename(pathGFF)
        orgName = file.replace(ext, "").replace("."+ext, "")
        dicoGFF = make_gff_dict(pathIN=pathGFF)
        dicoAnnot = {}
        dicoAttributes = {}
        pbar.set_description_str(orgName+" ".rjust(maxpathSize-len(orgName)))
        if not os.path.isdir(pathOUT):
            pathTABLE = pathOUT
        else:
            pathTABLE = pathOUT+"/"+orgName+"."+outFormat
        # ***** Get annotations features ***** #
        for seqType in dicoGFF[orgName]:
            if seqType[0] != "#" and seqType not in ["gene", "length"]:
                for dicoEntry in dicoGFF[orgName][seqType]:
                    lt = dicoEntry['attributes']['locus_tag']
                    if seqType == "tRNA":
                        dicoAnnot[lt] = {'product': dicoEntry['attributes']['product']}
                    elif seqType == "CDS":
                        for attrName in dicoEntry['attributes']:
                            if attrName not in ["ID", "Parent", "locus_tag"]:
                                if attrName not in dicoAttributes:
                                    dicoAttributes[len(dicoAttributes)] = attrName
                                value = dicoEntry['attributes'][attrName]
                                if attrName == "pvogs" and value != "":
                                    newValue = ""
                                    for pvogID in value.split(','):
                                        newValue += pvogID+", "
                                    value = newValue[: -1]
                                if attrName == "nahant" and value != "":
                                    setProduct = set()
                                    for nahant in value.split(','):
                                        product = " ".join(nahant.split(" [")[0].split(" ")[1:]).replace("MULTISPECIES: ", "")
                                        if "hypothetical protein" in product:
                                            setProduct.add("hypothetical protein")
                                        else:
                                            setProduct.add(product)
                                    newValue = ""
                                    for product in setProduct:
                                        newValue += product+", "
                                    value = newValue[: -1]
                                try:
                                    dicoAnnot[lt][attrName] = value
                                except KeyError:
                                    dicoAnnot[lt] = {'product': "", attrName: value}
        # ***** Write to TSV ***** #
        if outFormat == "tsv":
            OUT = open(pathTABLE, 'w')
            header = "Locus tag\tProduct"
            for i in range(len(dicoAttributes)):
                header += "\t"+dicoAttributes[i]
            OUT.write(header+"\n")
            for lt in sorted(dicoAnnot.keys()):
                line = lt+"\t"+dicoAnnot[lt]['product']
                for i in range(len(dicoAttributes)):
                    try:
                        line += "\t"+str(dicoAnnot[lt][dicoAttributes[i]])
                    except KeyError:
                        line += "\t."
                OUT.write(line+"\n")
            OUT.close()
        # ***** Write to XLSX ***** #
        if outFormat == "xlsx":
            workbook = xlsxwriter.Workbook(pathTABLE)
            worksheet = workbook.add_worksheet()
            headerFormat = workbook.add_format({'bold': True, 'align': 'center', 'valign': 'vcenter', 'border': 1, 'font_size': 12, 'font_color': "#ffffff", 'bg_color': "#333333"})
            rowFormat = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'border': 1, 'font_size': 11, 'text_wrap': True})
            dicoWidth = {0: len("locusTag")}
            dicoHeight = {0: 12}
            # Headers
            worksheet.write(0, 0, "locusTag", headerFormat)
            col = 1
            for i in range(len(dicoAttributes)):
                worksheet.write(0, col, dicoAttributes[i], headerFormat)
                dicoWidth[col] = len(dicoAttributes[i])
                col += 1
            # Cells
            row = 1
            for lt in sorted(dicoAnnot.keys()):
                worksheet.write(row, 0, lt, rowFormat)
                dicoWidth[0] = max(dicoWidth[0], len(lt))
                col = 1
                for i in range(len(dicoAttributes)):
                    strValue = ""
                    try:
                        splitValue = dicoAnnot[lt][dicoAttributes[i]].split(",")
                        try:
                            dicoHeight[row] = max(dicoHeight[row], len(splitValue)*12)
                        except KeyError:
                            dicoHeight[row] = len(splitValue)*12
                        for subValue in splitValue:
                            strValue += subValue+"\n"
                            dicoWidth[col] = max(dicoWidth[col], len(subValue))
                    except KeyError:
                        pass
                    worksheet.write(row, col, strValue[: -1], rowFormat)
                    col += 1
                row += 1
            # Adjust row height and column width
            for row in dicoHeight:
                worksheet.set_row(row, dicoHeight[row]*1.2)
            for col in dicoWidth:
                worksheet.set_column(col, col, min(dicoWidth[col], maxWidth)*1.2)
            # Freeze pane on the top row.
            worksheet.freeze_panes(1, 0)
            workbook.close()
        pbar.update(1)
        title("GFF2table", pbar)
    pbar.close()


@fct_checker
def make_gbk_from_fasta(pathIN1: str, pathIN2: str, pathIN3: str, pathOUT: str, identifier: str, topology: str, division: str, taxID: int = 0, pathIN4: str = None, boolProgress: bool = True) -> Tuple[str, str, str, str, str, str, str, int, str, bool]:
    '''
     ------------------------------------------------------------
    |                    MAKE GBK FROM FASTAS                    |
    |------------------------------------------------------------|
    |      Make a GBK file from FASTA (.ffn + .faa + .fna)       |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN1     : path of input FNA file (required)         |
    |    pathIN2     : path of input FFN file (required)         |
    |    pathIN3     : path of input FAA file (required)         |
    |    pathOUT     : path of output GBK file (required)        |
    |    identifier  : accession (required)                      |
    |    topology    : topology [linear or circular] (required)  |
    |    division    : division [BCT or PHG] (required)          |
    |    taxID       : genbank taxonomy ID (default=0)           |
    |    pathIN4     : path of trnascanse file (default=None)    |
    |    boolProgress: display progress (default=True)           |
     ------------------------------------------------------------
    |   taxID: Caudo=28883, Myo=10662, Podo=10744, Sipho=10699   |
     ------------------------------------------------------------
    '''
    pathIN1 = path_converter(pathIN1)
    pathIN2 = path_converter(pathIN2)
    pathIN3 = path_converter(pathIN3)
    if pathIN4 is not None and pathIN4 != "None":
        pathIN4 = path_converter(pathIN4)
    pathOUT = path_converter(pathOUT)
    if topology not in ["linear", "circular"]:
        printcolor("[ERROR: make_gbk_from_fasta]\tTopology must be 'linear' or 'circular'\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if division not in ["BCT", "PHG"]:
        printcolor("[ERROR: make_gbk_from_fasta]\nDivision must be 'BCT' for bacteria or 'PHG' for phage\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if boolProgress is True:
        printcolor("♊ FASTAs to GBK"+"\n")
    dicoFNA = make_fasta_dict(pathIN1, onlyLTasHeader=True)
    dicoFFN = make_fasta_dict(pathIN2, onlyLTasHeader=True)
    dicoFAA = make_fasta_dict(pathIN3, onlyLTasHeader=True)
    dicoDuplicatedSeq = {}
    dicoTemp = {}
    # Search duplicated gene sequence
    for key in dicoFFN:
        try:
            dicoTemp[dicoFFN[key]].append(key)
        except KeyError:
            dicoTemp[dicoFFN[key]] = [key]
    for key in dicoTemp:
        if len(dicoTemp[key]) > 1:
            dicoDuplicatedSeq[len(dicoDuplicatedSeq)] = dicoTemp[key]
    dicoTemp.clear()
    # Get trna if done
    if pathIN4 is not None and pathIN4 != "None":
        dicoTRNA = list(make_trnascanse_dict(pathIN=pathIN4, pathJSON="None", ext="."+pathIN4.split(".")[-1]).values())[0]
    else:
        dicoTRNA = {}
    lstGBKfiles = []
    if boolProgress is True:
        pbar = tqdm(total=len(dicoFFN), ncols=75, leave=True, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
    for contig in dicoFNA:
        contigName = contig.split(" ")[0].split("|")[0]
        if len(dicoFNA) > 1:
            identifier = identifier+"_"+contigName
        # Genome Sequence
        seqGenome = Seq(dicoFNA[contig].upper())
        # ***** MAIN FEATURES ***** #
        orgName = os.path.basename(pathIN1).replace(".fna", "").replace(".fasta", "")
        if len(dicoFNA) > 1:
            description = orgName.replace("_", " ")+", WGS linear, "+contigName+", whole genome shotgun"
        else:
            description = orgName.replace("_", " ")+", complete genome"
        if taxID == 0:
            taxonomy = []
        else:
            ncbi = NCBITaxa()
            taxonomy = []
            lstSubTaxID = ncbi.get_lineage(taxID)
            for subTaxID in lstSubTaxID:
                name = ncbi.get_taxid_translator([subTaxID])
                if name[subTaxID] != "root":
                    taxonomy.append(name[subTaxID])
        # ***** ANNOTATION dictionnary ***** #
        dicoAnnot = {'molecule_type': 'DNA', 'topology': topology, 'data_file_division': division,
                     'date': date.today().strftime("%d-%b-%Y").upper(),
                     'accessions': [identifier], 'sequence_version': 1,
                     'keywords': [''],
                     'source': orgName.replace("_", " "),
                     'organism': orgName.replace("_", " "),
                     'taxonomy': taxonomy,
                     'references': []}
        # ***** CREATE RECORD ***** #
        record = SeqRecord(seqGenome, id=identifier, name=identifier, description=description, annotations=dicoAnnot, features=None)
        dicoFeatures = {}
        # ***** FEATURES dictionnary ***** #
        for geneLT in dicoFFN:
            if len(dicoFNA) == 1 or contigName in geneLT:
                # Get gene Locations
                seq = dicoFFN[geneLT]
                dicoFindLocation = {}
                resFor = [i.start() for i in re.finditer(seq.upper(), str(seqGenome))]
                resRev = [i.start() for i in re.finditer(reverse_complement(seq).upper(), str(seqGenome))]
                for res in resFor:
                    dicoFindLocation[res] = 1
                for res in resRev:
                    dicoFindLocation[res] = -1
                if len(dicoFindLocation) == 0:
                    printcolor("[ERROR: make_gbk_from_fasta]\nUnable to find gene \""+geneLT+"\"\n", 1, "212;64;89", "None", True)
                    continue
                    # exit_gemini()
                elif len(dicoFindLocation) == 1:
                    start = list(dicoFindLocation.keys())[0]
                    end = start+len(seq)
                    strand = list(dicoFindLocation.values())[0]
                else:
                    for resStart in dicoFindLocation:
                        if resStart >= end-100 and resStart <= end+10000:
                            start = resStart
                            end = resStart+len(seq)
                            strand = dicoFindLocation[start]
                            break
                featureLocation = FeatureLocation(start, end, strand=strand)
                # Get gene Qualifiers (OrderedDict)
                orderGeneDicoQualifiers = OrderedDict([('locus_tag', [geneLT])])
                orderCDSDicoQualifiers = OrderedDict([('locus_tag', [geneLT]),
                                                      ('codon_start', ['1']),
                                                      ('transl_table', ['11']),
                                                      ('protein_id', [geneLT]),
                                                      ('translation', [dicoFAA[geneLT]])])
                geneSeqFeature = SeqFeature(location=featureLocation, type="gene", strand=strand, id=geneLT, qualifiers=orderGeneDicoQualifiers)
                CDSSeqFeature = SeqFeature(location=featureLocation, type="CDS", strand=strand, id=geneLT, qualifiers=orderCDSDicoQualifiers)
                dicoFeatures[start] = (geneSeqFeature, CDSSeqFeature)
                if boolProgress is True:
                    pbar.update(1)
                    title("writing", pbar)
        # ***** tRNA ***** #
        if contig in dicoTRNA:
            for trna in dicoTRNA[contig]:
                trnaLT = "tRNA_"+str(trna).zfill(3)
                if dicoTRNA[contig][trna]['pseudo'] is False:
                    featureLocation = FeatureLocation(dicoTRNA[contig][trna]['start'], dicoTRNA[contig][trna]['end'], strand=dicoTRNA[contig][trna]['strand'])
                    orderGeneDicoQualifiers = OrderedDict([('locus_tag', [trnaLT])])
                    orderTRNADicoQualifiers = OrderedDict([('locus_tag', [trnaLT]),
                                                           ('codon_start', ['1']),
                                                           ('note', ["tRNA "+dicoTRNA[contig][trna]['type']+" anticodon "+dicoTRNA[contig][trna]['codon']+" , score "+str(dicoTRNA[contig][trna]['score'])]),
                                                           ('product', [dicoTRNA[contig][trna]['type']+" tRNA"])])
                    geneSeqFeature = SeqFeature(location=featureLocation, type="gene", strand=dicoTRNA[contig][trna]['strand'], id=trnaLT, qualifiers=orderGeneDicoQualifiers)
                    TRNASeqFeature = SeqFeature(location=featureLocation, type="tRNA", strand=dicoTRNA[contig][trna]['strand'], id=trnaLT, qualifiers=orderTRNADicoQualifiers)
                    dicoFeatures[dicoTRNA[contig][trna]['start']] = (geneSeqFeature, TRNASeqFeature)
        # ***** Add sort features ***** #
        for pos in sorted(dicoFeatures):
            record.features.append(dicoFeatures[pos][0])
            record.features.append(dicoFeatures[pos][1])
        # ***** Save as GenBank file ***** #
        if len(dicoFNA) == 1:
            output_file = open(pathOUT, 'w')
        else:
            output_file = open(geminiset.pathTMP+"/"+contigName+".gbk", 'w')
            lstGBKfiles.append(geminiset.pathTMP+"/"+contigName+".gbk")
        SeqIO.write(record, output_file, 'genbank')
        output_file.close()
    if boolProgress is True:
        pbar.close()
    # Merge Genbank
    if len(dicoFNA) > 1:
        cat_lstfiles(lstGBKfiles, pathOUT)


'''
-------------------------------------------------------------------------------------------------------
                                             BLAST FUNCTIONS
-------------------------------------------------------------------------------------------------------
'''


@fct_checker
def make_blast_dict(pathIN: str, dbName: str, pathJSON: str = "None", idThr: int = 20, minLRthr: int = 50, maxLRthr: int = 50, ext: str = ".xml") -> Tuple[str, str, int, int, int, str]:
    '''
     ------------------------------------------------------------
    |            PARSE BLAST RESULTS (XML or TABULAR)            |
    |------------------------------------------------------------|
    |         Parse blast output and create a dictionnary        |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of XML blast file (required)             |
    |    dbName  : name of blast database (required)             |
    |    pathJSON: path of output JSON ("")                      |
    |    idThr   : %identity threshold (default=20)              |
    |    minLRthr: %minLrap threshold (default=50)               |
    |    maxLRthr: %maxLrap threshold (default=50)               |
    |    ext     : extension of input files (default=.xml)       |
    |RETURN                                                      |
    |    dicoBlast: {query: {subject: {hsp: {......}}}}      |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathJSON = path_converter(pathJSON)
    lstFiles, maxpathSize = get_input_files(pathIN, "make_blast_dict", [ext], filter=dbName)
    if len(lstFiles) == 0:
        printcolor("[ERROR: make_blast_dict]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    dicoBLAST = {}
    if pathJSON != "None" and os.path.isfile(pathJSON):
        dicoBLAST = load_json(pathJSON)
    else:
        minTabularField = ["qseqid", "sseqid|stitle", "pident", "length", "qlen", "qstart", "qend", "slen", "sstart", "send", "evalue"]
        for pathBLAST in lstFiles:
            file = os.path.basename(pathBLAST)
            orgName = file.replace(ext, "").replace("."+ext, "").replace("_diamond", "").replace("_blast", "").replace("_"+dbName, "")
            dicoBLAST[orgName] = {}
            # ***** PARSE XML BLAST ***** #
            if ext == ".xml" or ext == "xml":
                blastRecords = NCBIXML.parse(open(pathBLAST, "r"))
                # Browse blast query results
                for blastRec in blastRecords:
                    # Browse blast subject results
                    query = blastRec.query.split(" ")[0]
                    # replace transeq frame in header
                    # query = query.split(" ")[0]
                    # if query[-2:] == "_1": query = query[: -2]
                    queryLen = blastRec.query_length
                    for align in blastRec.alignments:
                        subject = align.title
                        subjectLen = align.length
                        # Browse HSPs
                        cptHsp = 0
                        for hsp in align.hsps:
                            # Compute percent identity and Lrap
                            percentId = (hsp.identities/hsp.align_length)*100
                            minLrap = (hsp.align_length/min(queryLen, subjectLen, hsp.align_length))*100
                            maxLrap = (hsp.align_length/max(queryLen, subjectLen))*100
                            # Filter on input threshold
                            if percentId >= idThr and minLrap >= minLRthr and maxLrap >= maxLRthr:
                                # Add to dictionnary
                                if query not in dicoBLAST[orgName]:
                                    dicoBLAST[orgName][query] = {'length': queryLen}
                                if subject not in dicoBLAST[orgName][query]:
                                    dicoBLAST[orgName][query][subject] = {'length': subjectLen}
                                # Add HSP
                                dicoBLAST[orgName][query][subject]["HSP"+str(cptHsp)] = {
                                                                                          'alignLen': hsp.align_length, 'evalue': hsp.expect,
                                                                                          'qstart': hsp.query_start, 'qend': hsp.query_end,
                                                                                          'sstart': hsp.sbjct_start, 'send': hsp.sbjct_end,
                                                                                          'match': hsp.identities, 'gap': hsp.gaps, 'mismatch': (hsp.align_length-hsp.gaps-hsp.identities),
                                                                                          'qalignseq': hsp.query, 'salignseq': hsp.sbjct,
                                                                                          'pident': percentId, 'minLrap': minLrap, 'maxLrap': maxLrap
                                                                                        }
                                cptHsp += 1
            # ***** PARSE TABULAR BLAST ***** #
            if ext in [".tsv", "tsv", "tabular", ".nahant"]:
                BLAST = open(pathBLAST, 'r')
                lstLines = BLAST.read().split("\n")
                BLAST.close()
                # retrieve header for tabular fields and check mandatory tabular fields
                searchFields = re.search("--outfmt 6 ([^-]+) -", lstLines[1])
                lstFields = searchFields.group(1).split(" ")
                for field in minTabularField:
                    splitAltField = field.split("|")
                    findField = False
                    for altfield in splitAltField:
                        if altfield in lstFields:
                            findField = True
                            break
                    if findField is False:
                        printcolor("[ERROR: make_blast_dict]\nMissing mandatory tabular field \""+field+"\"\n", 1, "212;64;89", "None", True)
                        exit_gemini()
                # parse results
                for i in range(3, len(lstLines), 1):
                    if lstLines[i] != "":
                        splitLine = lstLines[i].split("\t")
                        for j in range(len(splitLine)):
                            if lstFields[j] == "qseqid":
                                query = splitLine[j].replace("|-1", "").replace("|1", "")
                            elif lstFields[j] == "sseqid" or lstFields[j] == "stitle":
                                subject = splitLine[j]
                            elif lstFields[j] == "pident":
                                pident = float(splitLine[j])
                            elif lstFields[j] == "length":
                                length = int(splitLine[j])
                            elif lstFields[j] == "qlen":
                                qlen = int(splitLine[j])
                            elif lstFields[j] == "qstart":
                                qstart = int(splitLine[j])
                            elif lstFields[j] == "qend":
                                qend = int(splitLine[j])
                            elif lstFields[j] == "slen":
                                slen = int(splitLine[j])
                            elif lstFields[j] == "sstart":
                                sstart = int(splitLine[j])
                            elif lstFields[j] == "send":
                                send = int(splitLine[j])
                            elif lstFields[j] == "evalue":
                                evalue = float(splitLine[j])
                        # Maxrap/Minrap
                        minLrap = (length/min(qlen, slen))*100
                        maxLrap = (length/max(qlen, slen))*100
                        # Filter on input threshold
                        if pident >= idThr and minLrap >= minLRthr and maxLrap >= maxLRthr:
                            # Add to dictionnary
                            if query not in dicoBLAST[orgName]:
                                dicoBLAST[orgName][query] = {'length': qlen}
                            if subject not in dicoBLAST[orgName][query]:
                                dicoBLAST[orgName][query][subject] = {'length': slen}
                            # Add HSP
                            numHSP = len(dicoBLAST[orgName][query][subject])
                            dicoBLAST[orgName][query][subject]["HSP"+str(numHSP)] = {
                                                                                      'alignLen': length, 'evalue': evalue,
                                                                                      'qstart': qstart, 'qend': qend,
                                                                                      'sstart': sstart, 'send': send,
                                                                                      'pident': pident, 'minLrap': minLrap, 'maxLrap': maxLrap
                                                                                    }
        if pathJSON != "None":
            dump_json(dicoBLAST, pathJSON)
    return dicoBLAST


@fct_checker
def show_blast_dict(dicoBlast: dict) -> dict:
    '''
     ------------------------------------------------------------
    |                  DISPLAY BLAST DICTIONNARY                 |
    |------------------------------------------------------------|
    |             Display blast dictionnary results              |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    dicoBlast: blast dictionnary (required)                 |
     ------------------------------------------------------------
    '''
    for query in dicoBlast:
        printcolor("\nQuery: "+query+"\n")
        for subject in dicoBlast[query]:
            if subject != "length":
                printcolor("  "+subject.split("|")[0]+" ["+subject.split("|")[2]+"]"+"\n")
                for hsp in dicoBlast[query][subject]:
                    if hsp != "length":
                        printcolor("    HSP: %id = "+str(dicoBlast[query][subject][hsp]['percentId'])+" / maxLrap = "+str(dicoBlast[query][subject][hsp]['minLrap'])+" / minLrap = "+str(dicoBlast[query][subject][hsp]['maxLrap'])+"\n")


# ***** Reformat hmmscanDict to key==LT (and order by identity evalue) ***** #
def blastdict_to_annotsort(dicoBLAST):
    dicoANNOT = {}
    for orgName in dicoBLAST:
        dicoANNOT[orgName] = {}
        for lt in dicoBLAST[orgName]:
            dicoANNOT[orgName][lt] = {}
            for hit in dicoBLAST[orgName][lt]:
                if hit != "length":
                    if "|" in hit:
                        hitProduct = hit.split("|")[1].replace("MULTISPECIES: ", "")
                    else:
                        hitProduct = " ".join(hit.split(" [")[0].split(" ")[1:]).replace("MULTISPECIES: ", "")
                    for hsp in dicoBLAST[orgName][lt][hit]:
                        if hsp != "length":
                            evalue = dicoBLAST[orgName][lt][hit][hsp]['evalue']
                            pident = dicoBLAST[orgName][lt][hit][hsp]['pident']
                            try:
                                dicoANNOT[orgName][lt][evalue].add(hitProduct+"_____"+hit.split("|")[0]+"_____"+str(pident))
                            except KeyError:
                                dicoANNOT[orgName][lt][evalue] = set([hitProduct+"_____"+hit.split("|")[0]+"_____"+str(pident)])
    # Sort by evalue
    dicoANNOTSORT = {}
    for orgName in dicoANNOT:
        dicoANNOTSORT[orgName] = {}
        for lt in dicoANNOT[orgName]:
            dicoANNOTSORT[orgName][lt] = {}
            for evalue in sorted(dicoANNOT[orgName][lt].keys()):
                dicoANNOTSORT[orgName][lt][len(dicoANNOTSORT[orgName][lt])] = list(dicoANNOT[orgName][lt][evalue])
    return dicoANNOTSORT


'''
-------------------------------------------------------------------------------------------------------
                                             PARSER FUNCTIONS
-------------------------------------------------------------------------------------------------------
'''


@fct_checker
def make_hmmscan_dict(pathIN: str, pathJSON: str = "None", idEvalue: float = 0.01, ext: str = ".tblout") -> Tuple[str, str, float]:
    '''
     ------------------------------------------------------------
    |             PARSE TBLOUT HMMSCAN BLAST RESULTS             |
    |------------------------------------------------------------|
    |     Parse HMMSCAN tblout output and create a dictionnary   |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of TBLOUT hmmscan file (required)        |
    |    pathJSON: path of output JSON ("")                      |
    |    idEvalue: full sequence Evalue threshold (default=0.01) |
    |    ext     : extension of input files (default=.tblout)    |
    |RETURN                                                      |
    |    dicoHMMSCAN: {orgName: {query: {target: {......}}}}  |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathJSON = path_converter(pathJSON)
    lstFiles, maxpathSize = get_input_files(pathIN, "make_hmmscan_dict", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: make_hmmscan_dict]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    dicoHMMSCAN = {}
    if pathJSON != "None" and os.path.isfile(pathJSON):
        dicoHMMSCAN = load_json(pathJSON)
    else:
        for pathHMMSCAN in lstFiles:
            file = os.path.basename(pathHMMSCAN)
            orgName = file.replace(ext, "").replace("."+ext, "").replace("_hmmscan", "")
            dicoHMMSCAN[orgName] = {}
            OUT = open(pathHMMSCAN, 'r')
            lstLines = OUT.read().split("\n")
            OUT.close()
            for line in lstLines:
                if line != "" and line[0] != "#":
                    splitLine = re.compile(r"\s+").sub(" ", line).strip().split(" ")
                    target = splitLine[0].replace("|-1", "").replace("|1", "")
                    query = splitLine[2].replace("|-1", "").replace("|1", "")
                    # replace transeq frame in header
                    if query[-2:] == "_1":
                        query = query[: -2]
                    try:
                        fullSeqEvalue = float(splitLine[4])
                        if fullSeqEvalue <= idEvalue:
                            if query not in dicoHMMSCAN:
                                dicoHMMSCAN[orgName][query] = {}
                            dicoHMMSCAN[orgName][query][target] = {
                                                              'targetAcc': splitLine[1], 'queryAcc': splitLine[3], 'targetDescr': " ".join(splitLine[18:]).replace("\n", ""),
                                                              'fullSeqEvalue': fullSeqEvalue, 'fullSeqScore': float(splitLine[5]), 'fullSeqBias': float(splitLine[6]),
                                                              'bestDomEvalue': float(splitLine[7]), 'bestDomScore': float(splitLine[8]), 'bestDomBias': float(splitLine[9]),
                                                              'domCount_exp': float(splitLine[10]), 'domCount_reg': float(splitLine[11]), 'domCount_clu': float(splitLine[12]),
                                                              'domCount_ov': float(splitLine[13]), 'domCount_env': float(splitLine[14]), 'domCount_dom': float(splitLine[15]),
                                                              'domCount_rep': float(splitLine[16]), 'domCount_inc': float(splitLine[17])
                                                              }
                    except IndexError:
                        printcolor("[ERROR: make_hmmscan_dict]\nInvalid line for \""+file+"\"", 1, "212;64;89", "None", True)
                        exit_gemini()
        if pathJSON != "None":
            dump_json(dicoHMMSCAN, pathJSON)
    return dicoHMMSCAN


@fct_checker
def make_pvogs_desc_dict(pathIN: str = "None", pathJSON: str = "None") -> Tuple[str, str]:
    '''
     ------------------------------------------------------------
    |                  PARSE pVOGs TABLE FILES                   |
    |------------------------------------------------------------|
    |      Parse pVOGs table files and create a dictionnary      |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input pvogs annotation JSON ("")      |
    |    pathJSON: path of output JSON ("")                      |
    |RETURN                                                      |
    |    dicoPVOGS: {pVOGsID: set(descr, descr...)               |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathJSON = path_converter(pathJSON)
    dicoGeminiPath, dicoGeminiModule = get_gemini_path()
    lstHypoSynonym = ["hypothetical", "putative predicted product", "hyphothetical", "hypotheical", "unknown", "uncharacterised"]
    lstNoRealDescr = ["vs.", "77orf", "ab1gp", "bbp", "bcepgomrgp", "bcepny", "hef", "pas", "phg", "pfwmp", "phi eta orf", "phi eta orf", "phi mu50b-like", "phi pv83 orf", "phi pvl orf", "phi pvl-like", "phi slt orf", "phi105", "phi92_gp", "phi-eta", "phipvl", "rb32orf", "rb69orf", "upf", "vhs", "vpf", "vr7orf"]
    lstRegexpNoRealDescr = ["^gp", r"^[0-9]\.*[0-9]*$", r"^putative protein [0-9]\.*[0-9]*$", "^ea*a*[0-9]+", r"^gene\s*[0-9]+", r"^orf\s*.*[0-9]+", r"^p[0-9]\.*[0-9]*$"]
    setResultsPvogs = set()
    if pathIN != "None" and os.path.isfile(pathIN):
        dicoPVOGSresults = load_json(pathIN)
        for orgName in dicoPVOGSresults:
            setResultsPvogs.update(list(dicoPVOGSresults[orgName].keys()))
    else:
        for tableFile in os.listdir(dicoGeminiPath['DATABASES']['pvogs_table']):
            setResultsPvogs.add(tableFile.replace(".txt", ""))
    if pathJSON != "None" and os.path.isfile(pathJSON):
        dicoPVOGS = load_json(pathJSON)
    else:
        dicoPVOGS = {}
        pbar = tqdm(total=len(setResultsPvogs), ncols=75, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt} [{desc}]")
        for pvogsID in setResultsPvogs:
            pbar.set_description_str(pvogsID)
            pathTable = dicoGeminiPath['DATABASES']['pvogs_table']+"/"+pvogsID+".txt"
            TABLE = open(pathTable, 'r')
            lstLines = TABLE.read().split("\n")
            TABLE.close()
            pvogSetDescr = set()
            for line in lstLines:
                if line != "":
                    if line[0] != "#":
                        descr = line.lower().split("\t")[6].replace(";", ",")
                        reformat = False
                        if descr == "" or descr == "-" or descr == "phage protein" or descr == "conserved phage protein" or descr == "phage related protein" or descr == "cyanophage protein" or descr == "conserved cyanophage protein":
                            descr = "hypothetical protein"
                            reformat = True
                        else:
                            # Format redundant hypothetical
                            for hypoSynonym in lstHypoSynonym:
                                if hypoSynonym in descr and ", " not in descr:
                                    if "membrane" in descr:
                                        descr = "membrane protein"
                                    elif "structural" in descr:
                                        descr = "structural protein"
                                    elif "tail" in descr:
                                        descr = "tail protein"
                                    elif "transcriptional regulator" in descr:
                                        descr = "transcriptional regulator protein"
                                    else:
                                        descr = "hypothetical protein"
                                    reformat = True
                                    break
                            # Exclude ORFxx and others any clear description
                            if reformat is False:
                                for noRealDescr in lstNoRealDescr:
                                    if noRealDescr in descr:
                                        descr = "hypothetical protein"
                                    reformat = True
                                    break
                            if reformat is False:
                                for regexpNoRealDescr in lstRegexpNoRealDescr:
                                    searchRegexp = re.search(regexpNoRealDescr, descr)
                                    if searchRegexp:
                                        descr = "hypothetical protein"
                                        reformat = True
                                        break
                            # Format domain == domain containing == domain-containing
                            if reformat is False and "duf" in descr:
                                searchDuf = re.search("(duf[0-9]+)", descr)
                                descr = searchDuf.group(1)+" domain protein"
                        # Add conserved feature
                        if reformat is True and "conserved" in descr:
                            descr = "conserved "+descr
                        # Add description
                        pvogSetDescr.add(descr.replace(", ", " ").replace(";", " "))
            # Reduce descriptions per table
            if len(pvogSetDescr) > 1:
                try:
                    pvogSetDescr.remove("hypothetical protein")
                except KeyError:
                    pass
                pvogSetDescr2 = pvogSetDescr.copy()
                for descr in pvogSetDescr2:
                    if "putative" in descr and descr.replace("putative", "") in pvogSetDescr:
                        pvogSetDescr.remove(descr)
                    for descr2 in pvogSetDescr2:
                        if descr != descr2 and descr in descr2:
                            pvogSetDescr.remove(descr)
                            break
            # Add description list to pVOGs dico
            dicoPVOGS[pvogsID] = list(pvogSetDescr)
            pbar.update(1)
            title("pvogsDescDict", pbar)
        pbar.close()
        # Dump & return
        if pathJSON != "None":
            dump_json(dicoPVOGS, pathJSON)
    return dicoPVOGS


@fct_checker
def make_trnascanse_dict(pathIN: str, pathJSON: str = "None", ext: str = ".trnascanse") -> Tuple[str, str, str]:
    '''
     ------------------------------------------------------------
    |                PARSE tRNAscan-SE TABLE FILES               |
    |------------------------------------------------------------|
    |      Parse tRNAscan-SE output and create a dictionnary     |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input files or folder (required)      |
    |    pathJSON: path of output JSON ("")                      |
    |    ext     : extension of input files (default=.trnascanse)|
    |RETURN                                                      |
    |    dicoTRNA: {orgName: {lt: tRNA, ...} ...}              |
     ------------------------------------------------------------
    '''
    pathJSON = path_converter(pathJSON)
    lstFiles, maxpathSize = get_input_files(pathIN, "make_trnascanse_dict", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: make_trnascanse_dict]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if pathJSON != "None" and os.path.isfile(pathJSON):
        dicoTRNA = load_json(pathJSON)
    else:
        dicoTRNA = {}
        for pathTRNASCANSE in lstFiles:
            file = os.path.basename(pathTRNASCANSE)
            orgName = file.replace(ext, "").replace("."+ext, "")
            dicoTRNA[orgName] = {}
            if os.path.getsize(pathTRNASCANSE) != 0:
                TRNASCANSE = open(pathTRNASCANSE, 'r')
                lstLines = TRNASCANSE.read().split("\n")
                TRNASCANSE.close()
                for i in range(3, len(lstLines)-1, 1):
                    splitLine = lstLines[i].split("\t")
                    contig = splitLine[0].strip()
                    if contig not in dicoTRNA[orgName]:
                        dicoTRNA[orgName][contig] = {}
                    start = int(splitLine[2])
                    end = int(splitLine[3])
                    score = float(splitLine[8])
                    if start > end:
                        strand = -1
                        start = int(splitLine[3])
                        end = int(splitLine[2])
                    else:
                        strand = 1
                    aa = splitLine[4]
                    codon = splitLine[5]
                    pseudoBool = splitLine[9] == "pseudo"
                    dicoTRNA[orgName][contig][len(dicoTRNA[orgName][contig])+1] = {'start': start, 'end': end, 'type': aa, 'codon': codon, 'strand': strand, 'pseudo': pseudoBool, 'score': score}
        if pathJSON != "None":
            dump_json(dicoTRNA, pathJSON)
    return dicoTRNA


@fct_checker
def reformat_phanotate(pathIN: str) -> Tuple[str]:
    '''
     ------------------------------------------------------------
    |                FORMAT Phanotate FASTA header               |
    |------------------------------------------------------------|
    |        Reformat ORF name in Phanotate output FASTA         |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN: path of input Phanotate FASTA file (required)   |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    dicoFasta = make_fasta_dict(pathIN)
    cptOrf = 1
    OUT = open(pathIN, 'w')
    for key in dicoFasta:
        newHeader = key.split("|")[0]+"__ORF"+str(cptOrf)
        OUT.write(">"+newHeader+"\n"+dicoFasta[key]+"\n")
        cptOrf += 1
    OUT.close()


@fct_checker
def reformat_panacota(pathIN: str) -> Tuple[str]:
    '''
     ------------------------------------------------------------
    |           REFORMAT PanACoTA organism name output           |
    |------------------------------------------------------------|
    |       Reformat organism name in PanACoTA output files      |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN: path of PanACoTA results folder (required)      |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    if not os.path.isdir(pathIN):
        printcolor("[ERROR: reformat_panacota]\nPanACoTA results folder not found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    # Read organism/id association file (LSTINFO-list_file)
    printcolor("♊ Read association file"+"\n")
    pathLSTINFO = pathIN+"/LSTINFO-list_file.lst"
    dicoIDtoOrgName = {}
    if not os.path.isfile(pathLSTINFO):
        printcolor("[ERROR: reformat_panacota]\nPanACoTA \""+pathLSTINFO+"\" not found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    LSTINFO = open(pathLSTINFO, 'r')
    lstLines = LSTINFO.read().split("\n")[1: -1]
    LSTINFO.close()
    for line in lstLines:
        panacotaID = line.split("\t")[0]
        orgName = os.path.basename(line.split("\t")[1]).replace(".fna", "").replace(".fa", "").replace(".fasta", "")
        dicoIDtoOrgName[panacotaID] = orgName
    printcolor("⏩ Found "+str(len(dicoIDtoOrgName))+" identifiers"+"\n")
    # Retrieve all files
    lstFiles = []
    for root, directories, file in os.walk(pathIN):
        for file in file:
            lstFiles.append(os.path.join(root, file))
    # Rename all files containing PanACoTA identifiers
    printcolor("♊ Rename Files"+"\n")
    cptFileRename = 0
    pbar = tqdm(total=len(lstFiles), ncols=75, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
    for filepath in lstFiles:
        for key in dicoIDtoOrgName:
            if key in filepath:
                pathDist = filepath.replace(key, dicoIDtoOrgName[key])
                shutil.move(filepath, pathDist)
                cptFileRename += 1
                break
        pbar.update()
    pbar.close()
    printcolor("⏩ "+str(cptFileRename)+" files renamed"+"\n")
    # Replace PanACoTA identifiers inside files
    printcolor("♊ Replace Files"+"\n")
    lstFiles = []
    for root, directories, file in os.walk(pathIN):
        for file in file:
            lstFiles.append(os.path.join(root, file))
    cptFileReplace = 0
    pbar = tqdm(total=len(lstFiles), ncols=75, leave=False, desc="", file=sys.stdout, bar_format="  {percentage: 3.0f}%|{bar}| {n_fmt}/{total_fmt}")
    for filepath in lstFiles:
        if filepath != pathLSTINFO:
            SRC = open(filepath, 'r')
            try:
                data = SRC.read()
            except UnicodeDecodeError:
                continue
            SRC.close()
            boolReplace = False
            for key in dicoIDtoOrgName:
                if key in data:
                    boolReplace = True
                    data = data.replace(key, dicoIDtoOrgName[key])
            if boolReplace is True:
                SRC = open(filepath, 'w')
                SRC.write(data)
                SRC.close()
                cptFileReplace += 1
        pbar.update()
    pbar.close()
    printcolor("⏩ "+str(cptFileReplace)+" files replaced"+"\n")


@fct_checker
def make_mmseqs_cluster_dict(pathIN: str, pathJSON: str = "None") -> Tuple[str, str]:
    '''
     ------------------------------------------------------------
    |               MAKE MMSEQS CLUSTER DICTIONNARY              |
    |------------------------------------------------------------|
    |   Parse MMSEQS cluster tsv file and create a dictionnary   |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input files or folder (required)      |
    |    pathJSON: path of output JSON ("")                      |
    |RETURN                                                      |
    |    dicoCLUSTER: {clusterNum: [lt, ..., lt], ...}         |
     ------------------------------------------------------------
    '''
    dicoCLUSTER = {}
    pathIN = path_converter(pathIN)
    pathJSON = path_converter(pathJSON)
    if not os.path.isfile(pathIN):
        printcolor("[ERROR: make_mmseqs_cluster_dict]\nInput file not found\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if pathJSON != "None" and os.path.isfile(pathJSON):
        dicoCLUSTER = load_json(pathJSON)
    else:
        # Read mmseqs TSV file
        MMSEQS = open(pathIN, 'r')
        lstLines = MMSEQS.read().split("\n")
        MMSEQS.close()
        # Parse and assign cluster
        for line in lstLines:
            if line != "":
                lt1 = line.split("\t")[0].replace("|1", "").replace("|-1", "")
                lt2 = line.split("\t")[1].replace("|1", "").replace("|-1", "")
                findCluster = False
                for clusterNum in dicoCLUSTER.keys():
                    if lt1 in dicoCLUSTER[clusterNum]:
                        dicoCLUSTER[clusterNum].add(lt2)
                        findCluster = True
                        break
                    if lt2 in dicoCLUSTER[clusterNum]:
                        dicoCLUSTER[clusterNum].add(lt1)
                        findCluster = True
                        break
                if findCluster is False:
                    dicoCLUSTER["cluster"+str(len(dicoCLUSTER)+1).zfill(4)] = set([lt1, lt2])
        # Convert set to list
        for clusterNum in dicoCLUSTER.keys():
            dicoCLUSTER[clusterNum] = list(dicoCLUSTER[clusterNum])
        # Dump and return dictionnary
        if pathJSON != "None":
            dump_json(dicoCLUSTER, pathJSON)
    return dicoCLUSTER


@fct_checker
def make_eggnog_dict(pathIN: str, pathJSON: str = "None", ext: str = ".annotations") -> Tuple[str, str, str]:
    '''
     ------------------------------------------------------------
    |                     PARSE EggNOG FILES                     |
    |------------------------------------------------------------|
    |        Parse EggNOG output and create a dictionnary        |
    |------------------------------------------------------------|
    |PARAMETERS                                                  |
    |    pathIN  : path of input files or folder (required)      |
    |    pathJSON: path of output JSON ("")                      |
    |    ext     : ext of input files (default=.annotations)     |
    |RETURN                                                      |
    |    dicoEGGNOG: {}                                         |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathJSON = path_converter(pathJSON)
    lstFiles, maxpathSize = get_input_files(pathIN, "make_eggnog_dict", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: make_eggnog_dict]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if pathJSON != "None" and os.path.isfile(pathJSON):
        dicoEGGNOG = load_json(pathJSON)
    else:
        dicoEGGNOG = {}
        for pathEGGNOG in lstFiles:
            file = os.path.basename(pathEGGNOG)
            name = file.replace(ext, "").replace("."+ext, "")
            dicoEGGNOG[name] = {}
            EGGNOG = open(pathEGGNOG, 'r')
            lstLines = EGGNOG.read().split("\n")
            EGGNOG.close()
            for line in lstLines:
                if line != "" and line[0] != "#":
                    splitLine = line.split("\t")
                    lt = splitLine[0].replace("|-1", "").replace("|1", "")
                    ortholog = splitLine[1]
                    evalue = splitLine[2]
                    ogs = splitLine[4]
                    finalOg = ""
                    for og in ogs.split(","):
                        identifier = og.split("@")[0]
                        taxonomy = og.split("|")[1]
                        if "Bacteria" in ogs or "Viruses" in ogs:
                            if taxonomy in ["Bacteria", "Viruses"]:
                                finalOg = identifier
                                break
                        elif taxonomy in ["Caudovirales", "Myoviridae", "Podoviridae", "Siphoviridae"]:
                            finalOg = identifier
                            break
                    maxAnnotLvl = splitLine[5]
                    COGcat = splitLine[6]
                    description = splitLine[7].replace(", ", "_").replace(";", "_")
                    prefName = splitLine[8].replace(", ", " ").replace(";", " ")
                    dicoEGGNOG[name][lt] = {
                                            'ortholog': ortholog, 'evalue': evalue, 'ogs': ogs, 'finalOg': finalOg, 'COGcat': COGcat,
                                            'maxAnnotLvl': maxAnnotLvl, 'descr': description, 'prefName': prefName
                                           }
        if pathJSON != "None":
            dump_json(dicoEGGNOG, pathJSON)
    return dicoEGGNOG


@fct_checker
def make_interpro_dict(pathIN: str, idEvalue: float = 0.01, pathJSON: str = "None", ext: str = ".tsv") -> Tuple[str, float, str, str]:
    '''
     ------------------------------------------------------------
    |                PARSE INTERPROSCAN TSV FILES                |
    |------------------------------------------------------------|
    |     Parse InterProScan output and create a dictionnary     |
    |------------------------------------------------------------|
    | PARAMETERS                                                 |
    |    pathIN  : path of input files or folder (required)      |
    |    pathJSON: path of output JSON ("")                      |
    |    idEvalue: evalue threshold (default=0.01)               |
    |    ext     : extension of input files (default=.tsv)       |
    |RETURN                                                      |
    |    dicoIPRSCAN: {}                                        |
     ------------------------------------------------------------
    '''
    pathIN = path_converter(pathIN)
    pathJSON = path_converter(pathJSON)
    lstFiles, maxpathSize = get_input_files(pathIN, "make_interpro_dict", [ext])
    if len(lstFiles) == 0:
        printcolor("[ERROR: make_interpro_dict]\nAny input files found, check extension\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if pathJSON != "None" and os.path.isfile(pathJSON):
        dicoIPRSCAN = load_json(pathJSON)
    else:
        dicoIPRSCAN = {}
        for pathIPRSCAN in lstFiles:
            file = os.path.basename(pathIPRSCAN)
            name = file.replace(ext, "").replace("."+ext, "")
            dicoIPRSCAN[name] = {}
            IPRSCAN = open(pathIPRSCAN, 'r')
            lstLines = IPRSCAN.read().split("\n")
            IPRSCAN.close()
            for line in lstLines:
                if line != "" and line[0] != "#":
                    splitLine = line.split("\t")
                    lt = splitLine[0].replace("|-1", "").replace("|1", "")
                    # md5 = splitLine[1]
                    # length = splitLine[2]
                    # analysis = splitLine[3]  # (e.g. Pfam / PRINTS / Gene3D)
                    signatureAcc = splitLine[4]
                    signatureDescr = splitLine[5].replace(", ", " ").replace(";", " ")
                    # start = splitLine[6]
                    # end = splitLine[7]
                    score = splitLine[8]
                    # matchStatus = splitLine[9]
                    # date = splitLine[10]
                    interproAnnot = splitLine[11]  # (e.g. IPR002093)
                    interproDescr = splitLine[12].replace(", ", " ").replace(";", " ")  # (e.g. BRCA2 repeat)
                    if lt not in dicoIPRSCAN[name]:
                        dicoIPRSCAN[name][lt] = {}
                    # Add interpro annotation or signature if not available
                    if interproAnnot == "-":
                        dicoIPRSCAN[name][lt][signatureAcc] = {'descr': signatureDescr, 'evalue': score}
                    else:
                        dicoIPRSCAN[name][lt][interproAnnot] = {'descr': interproDescr, 'evalue': score}
        if pathJSON != "None":
            dump_json(dicoIPRSCAN, pathJSON)
    return dicoIPRSCAN
