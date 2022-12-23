'''
┌─┐  ┌─┐  ┌┬┐  ┬  ┌┐┌  ┬
│ ┬  ├┤   │││  │  │││  │
└─┘  └─┘  ┴ ┴  ┴  ┘└┘  ┴ ni
---------------------------------------------------
-*- coding: utf-8 -*-                              |
title            : geminini.py                     |
description      : Gemini internal functions       |
author           : dooguypapua                     |
lastmodification : 20210626                        |
version          : 0.1                             |
python_version   : 3.8.5                           |
---------------------------------------------------
'''

import os
import sys
import ast
# import codecs
import json
import tempfile
import psutil
import pprint
import subprocess
import threading
import itertools
import pickle
import sqlite3
import blosc2
import gzip
# import time
import glob
import torch
import torch.nn.functional as F
import geminiset
from tqdm import tqdm
from yaspin import yaspin
from yaspin.spinners import Spinners
from functools import wraps, lru_cache
from math import floor
from random import random
from rich.console import Console
from rich.syntax import Syntax
from scipy.special import expit
from operator import itemgetter


'''
-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------
                                          INTERNAL FUNCTIONS
-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------
'''


# ***** Init internal gemini dico ***** #
def init_gemini_dico():
    pathGeminiJson = os.path.dirname(os.path.abspath(__file__))+"/conf/geminison.json"
    global dicoGemini
    dicoGemini = load_json(pathGeminiJson)
    return dicoGemini


# ***** Gemini exit ***** #
def exit_gemini():
    title(text="⛔", pbar=None)
    exit("\n")


# ***** Load JSON file ***** #
def load_json(pathJSON):
    try:
        if pathJSON.endswith(".gz"):
            with gzip.open(pathJSON, "r") as json_file:
                dico = json.load(json_file)
        else:
            with open(pathJSON, 'r') as json_file:
                dico = json.load(json_file)
        return dico
    except (UnicodeDecodeError, AttributeError):
        printcolor("[ERROR: load_json]\nUnable to load json file \""+pathJSON+"\"\n", 1, "212;64;89", "None", True)
        exit_gemini()


# ***** Dump JSON file ***** #
def dump_json(dico, pathJSON, indent=4):
    try:
        if pathJSON.endswith(".gz"):
            with gzip.open(pathJSON, "wt") as outfile:
                json.dump(dico, outfile, indent=indent)
        else:
            with open(pathJSON, 'w') as outfile:
                json.dump(dico, outfile, indent=indent)
    except (UnicodeDecodeError, AttributeError):
        if os.path.isfile(pathJSON):
            os.remove(pathJSON)
        printcolor("[ERROR: dump_json]\nUnable to dump json file \""+pathJSON+"\"\n", 1, "212;64;89", "None", True)
        exit_gemini()


# ***** FileReader ***** #
def read_file(file_path: str, excludeFirstKr: str = "#", yaspinBool: bool = True) -> []:
    if not os.path.isfile(file_path):
        printcolor("File not found \""+file_path+"\"\n", 1, "212;64;89", "None", True)
        exit_gemini()
    if yaspinBool:
        spinner = yaspin(Spinners.simpleDotsScrolling, text="reading", side="right")
        spinner.start()
    lstLines = []
    with open(file_path) as f:
        for line in f:
            if line[0] != excludeFirstKr:
                lstLines.append(line[:-1])
    if yaspinBool:
        spinner.stop()
    return lstLines


# ***** Concatenate files ***** #
def cat_lstfiles(lstFiles, pathOUT):
    os.system("cp "+lstFiles[0]+" "+pathOUT)
    for i in range(1, len(lstFiles), 1):
        os.system("cat "+lstFiles[i]+" >> "+pathOUT)


# ***** Instanciate internal gemini dico ***** #
def build_gemini_dico(lstArgv):
    fctName = lstArgv[1]
    slurmBool, cpu, memMax, memMin = get_sys_info()
    # temporary folder
    if dicoGemini[fctName]['pathTMP'] is True:
        if slurmBool is True:
            tmpDir = "/scratch2/sbr/dynamic/tmp"
        else:
            tmpDir = "/tmp"
        value = tempfile.mkdtemp(prefix="gemini_"+fctName+"_", dir=tmpDir)
        dicoGemini[fctName]['pathTMP'] = value
    else:
        dicoGemini[fctName]['pathTMP'] = ""
    for i in range(2, len(lstArgv)-1, 2):
        tagArg = lstArgv[i]
        if tagArg == "--cpu":
            continue
        valueArg = lstArgv[i+1]
        findTag = False
        # Apply value
        for fctArg in dicoGemini[fctName]['dicoArgs']:
            if dicoGemini[fctName]['dicoArgs'][fctArg]['tag'] == tagArg:
                value = valueArg
                if value == "true" or value == "True":
                    value = True
                elif value == "false" or value == "False":
                    value = False
                elif value == "None" or value == "none":
                    value = None
                else:
                    try:
                        value = ast.literal_eval(value)
                    except (ValueError, SyntaxError):
                        pass
                dicoGemini[fctName]['dicoArgs'][fctArg]['value'] = value
                findTag = True
        # Check tag
        if not findTag:
            infofct(fctName)
            printcolor("[ERROR]\nInvalid arguments \""+tagArg+"\"\n", 1, "212;64;89", "None", True)
            exit_gemini()


# ***** Retrieve tools path ***** #
def get_gemini_path():
    dicoGeminiPath = {}
    slurmBool, cpu, memMax, memMin = get_sys_info()
    if slurmBool is True:
        pathGeminiPathTXT = os.path.dirname(os.path.abspath(__file__))+"/conf/geminipathslurm.txt"
    else:
        pathGeminiPathTXT = os.path.dirname(os.path.abspath(__file__))+"/conf/geminipath.txt"
    IN = open(pathGeminiPathTXT, 'r')
    lstLines = IN.read().split("\n")
    IN.close()
    for line in lstLines:
        if line != "" and line[0] != "#":
            dicoGeminiPath[line.split(":")[0]] = line.split(":")[1]
    return dicoGeminiPath


# ***** Retrieve system info ***** #
def get_sys_info():
    # return slurmBool, cpu, memMax, memMin
    if "isslurm" in os.environ and os.environ["isslurm"] == "true":
        if "geminicpu" in os.environ and "geminimem" in os.environ:
            return True, int(os.environ["geminicpu"]), int(os.environ["geminimem"]), int(int(os.environ["geminimem"])/2)
        else:
            return True, 32, 64, 48
    else:
        memMax = int(round((psutil.virtual_memory().total/(1024*1024*1024)), 0))
        memMin = int(memMax*75/100)
        return False, 12, memMax, memMin


# ***** Convert arguments to parameters ***** #
def fct_params_from_args(fctName):
    lstArgValue = []
    for fctArg in dicoGemini[fctName]['dicoArgs']:
        lstArgValue.append(dicoGemini[fctName]['dicoArgs'][fctArg]['value'])
    return lstArgValue


# ***** Check function parameters ***** #
def fct_checker(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        dicoGemini = init_gemini_dico()
        dicoFunc = dicoGemini[func.__name__]["dicoArgs"]
        lstArgsFunc = list(dicoFunc.keys())
        # Arguments are passed without name
        if len(args) > 0:
            for i in range(len(args)):
                dicoArg = dicoFunc[lstArgsFunc[i]]
                if type(args[i]).__name__ == dicoArg['type']:
                    dicoFunc[lstArgsFunc[i]]["checktype"] = ""
                else:
                    dicoFunc[lstArgsFunc[i]]["checktype"] = (dicoArg['type'], type(args[i]).__name__)
        # Arguments are passed using name
        else:
            for argName1 in kwargs:
                for argName2 in dicoFunc:
                    if argName1 == argName2:
                        if type(kwargs[argName1]).__name__ == dicoFunc[argName2]['type']:
                            dicoFunc[argName2]["checktype"] = ""
                        else:
                            dicoFunc[argName2]["checktype"] = (dicoArg['type'], type(args[i]).__name__)
                        break
        # Check missing required parameters and incompatible type
        lstMissingArgs = []
        lstBadtypeArgs = []
        for argName in dicoFunc:
            # Missing argument
            if "checktype" not in dicoFunc[argName] and dicoFunc[argName]['value'] == "":
                lstMissingArgs.append(argName)
            if "checktype" in dicoFunc[argName] and dicoFunc[argName]["checktype"] != "":
                lstBadtypeArgs.append("Invalid type for \""+argName+"\", found \""+dicoFunc[argName]["checktype"][1]+"\" but \""+dicoFunc[argName]["checktype"][0]+"\" is required")
        # Display errors
        if len(lstMissingArgs) > 0 or len(lstBadtypeArgs) > 0:
            infofct(func.__name__)
            printcolor("[ERROR: "+func.__name__+"]\n", 1, "212;64;89", "None", True)
            if len(lstMissingArgs) > 0:
                printcolor("Missing parameters: "+", ".join(lstMissingArgs)+"\n", 1, "212;64;89", "None", True)
            if len(lstBadtypeArgs) > 0:
                printcolor("\n".join(lstBadtypeArgs)+"\n", 1, "212;64;89", "None", True)
            exit_gemini()
        return func(*args, **kwargs)
    return wrapper


# ***** Convert windows to unix path ***** #
def path_converter(pathToConvert):
    if pathToConvert is None:
        pathConverted = pathToConvert
    elif "\\" in pathToConvert:
        pathConverted = pathToConvert.replace("C:\\", "/mnt/c/").replace("\\", "/")
    else:
        pathConverted = pathToConvert
    if pathConverted[-1] == "/":
        return pathConverted[:-1]
    else:
        return pathConverted


# ***** Get input files from file path or folder path, count and pathsize ***** #
def get_input_files(pathIN, fctName, fileExt=[""], filter=None):
    maxpathSize = 0
    lstAllFiles = []
    setFiles = set()
    pathIN = path_converter(pathIN)
    if os.path.isfile(pathIN):
        lstAllFiles.append(pathIN)
    elif os.path.isdir(pathIN):
        for file in os.listdir(pathIN):
            lstAllFiles.append(pathIN+"/"+file)
    else:
        printcolor("[ERROR: "+fctName+"]\nUnable to find input file \""+pathIN+"\"\n", 1, "212;64;89", "None", True)
        exit_gemini()
    # Check extension
    for file in lstAllFiles:
        for ext in fileExt:
            if file.endswith(ext) and (filter is None or filter in file):
                setFiles.add(file)
                maxpathSize = max(len(os.path.basename(file).replace(ext, "")), maxpathSize)
    return list(sorted(setFiles)), maxpathSize


# ***** Reverse complement ***** #
def reverse_complement(seq):
    tab = str.maketrans("ACTG", "TGAC")
    return seq.translate(tab)[::-1]


# ***** Parse tabulated file with group definition ***** #
def read_group(pathGROUP):
    dicoGroup = {}
    maxGroupStrSize = 0
    GROUP = open(pathGROUP, 'r')
    lstLines = GROUP.read().split("\n")
    GROUP.close()
    for line in lstLines:
        if line != "":
            orgName = line.split("\t")[0]
            groupName = line.split("\t")[1]
            maxGroupStrSize = max(len(groupName), maxGroupStrSize)
            try:
                dicoGroup[groupName].append(orgName)
            except KeyError:
                dicoGroup[groupName] = [orgName]
    return dicoGroup, maxGroupStrSize


# ***** Group integer list by range ***** #
def to_ranges(listIN):
    listIN = sorted(set(listIN))
    for key, group in itertools.groupby(enumerate(listIN),
                                        lambda t: t[1] - t[0]):
        group = list(group)
        yield group[0][1], group[-1][1]


# ***** FIND longest substring between two strings ***** #
def longest_common_substring(x: str, y: str) -> (int, int, int):
    # function to find the longest common substring
    # Memorizing with maximum size of the memory as 1
    @lru_cache(maxsize=1)
    # function to find the longest common prefix
    def longest_common_prefix(i: int, j: int) -> int:
        if 0 <= i < len(x) and 0 <= j < len(y) and x[i] == y[j]:
            return 1 + longest_common_prefix(i + 1, j + 1)
        else:
            return 0
    # diagonally computing the subproblems
    # to decrease memory dependency

    def digonal_computation():
        # upper right triangle of the 2D array
        for k in range(len(x)):
            yield from ((longest_common_prefix(i, j), i, j)
                        for i, j in zip(range(k, -1, -1), range(len(y) - 1, -1, -1)))
        # lower left triangle of the 2D array
        for k in range(len(y)):
            yield from ((longest_common_prefix(i, j), i, j)
                        for i, j in zip(range(k, -1, -1), range(len(x) - 1, -1, -1)))
    # returning the maximum of all the subproblems
    return max(digonal_computation(), key=itemgetter(0), default=(0, 0, 0))


# ***** sqlitedict ***** #
def my_encode(obj):
    return sqlite3.Binary(blosc2.compress2(pickle.dumps(obj, pickle.HIGHEST_PROTOCOL)))


def my_decode(obj):
    return pickle.loads(blosc2.decompress2(bytes(obj)))


'''
-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------
                                         SUBPROCESS FUNCTIONS
-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------
'''


# ***** Simple threading subprocess ***** #
def launch_subprocess_thread(dicoThread, thread_name, pathLOG):
    process = subprocess.Popen(dicoThread[thread_name]["cmd"], stdout=subprocess.PIPE, stderr=open(pathLOG, 'a'), shell=True)
    process.wait()
    for line in process.stdout.readlines():
        dicoThread[thread_name]["returnlines"].append(str(line.decode('utf-8').strip()))
    dicoThread[thread_name]["returnstatut"] = process.returncode


# ***** Multithread ***** #
def launch_threads(dicoThread, job_prefix, maxThread, pathTMP, spinner=None):
    lstThread = []
    lstThreadcopy = []
    PrevNbFinishThread = 0
    for thread_name in dicoThread.keys():
        pathLOG = pathTMP+"/"+job_prefix+"_"+thread_name+".log"
        LOG = open(pathLOG, 'w')
        LOG.write("#CMD: "+dicoThread[thread_name]["cmd"]+"\n\n")
        LOG.close()
        lstThread.append(threading.Thread(target=launch_subprocess_thread, args=(dicoThread, thread_name, pathLOG), name=job_prefix+"_"+thread_name))
    # Tasks variables
    NbTasks = len(lstThread)
    lstThreadcopy = lstThread.copy()
    # Launch threads
    if NbTasks > 0:
        # Continue while some tasks are availables
        while (len(lstThread) != 0 or job_prefix in str(threading.enumerate())):
            nbThreadActive = str(threading.enumerate()).count(job_prefix)
            if nbThreadActive < maxThread:
                try:
                    thread = lstThread.pop()
                    thread.start()
                    nbThreadActive += 1
                except IndexError:
                    pass
            NbFinishThread = NbTasks-len(lstThread)-nbThreadActive
            if PrevNbFinishThread != NbFinishThread:
                PrevNbFinishThread = NbFinishThread
            # check error
            for thread_name in dicoThread.keys():
                if dicoThread[thread_name]["returnstatut"] != 0:
                    break
            # time.sleep(1)
        for thread in lstThreadcopy:
            try:
                thread.join()
            except RuntimeError:
                pass
    # Check dicoThread
    error = False
    for thread_name in dicoThread.keys():
        pathLOG = pathTMP+"/"+job_prefix+"_"+thread_name+".log"
        if dicoThread[thread_name]["returnstatut"] != 0:
            printcolor("\n[ERROR: thread \""+job_prefix+"\":\""+thread_name+"\"]\nCheck log file in \""+pathLOG+"\"\n", 1, "212;64;89", "None", True)
            error = True
        else:
            os.remove(pathLOG)
    if error is True:
        if spinner is not None:
            spinner.stop()
        exit_gemini()


'''
-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------
                                           BALROG FUNCTIONS
-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------
'''


# Convert amino acid letters to integers.
def tokenize_aa_seq(aa_seq):
    table = {"L": 1, "V": 2, "I": 3, "M": 4, "C": 5, "A": 6, "G": 7, "S": 8, "T": 9, "P": 10, "F": 11,
             "Y": 12, "W": 13, "E": 14, "D": 15, "N": 16, "Q": 17, "K": 18, "R": 19, "H": 20, "*": 0, "X": 0}
    tokenized = torch.tensor([table[aa] for aa in aa_seq])
    return tokenized


def get_start_codon(seq, orfcoords, strand):
    # forward strand
    if strand == 1:
        startcoord = orfcoords[0]
        return seq[startcoord-3:startcoord]
    # reverse strand
    else:
        startcoord = orfcoords[1]
        return seq[startcoord:startcoord+3].reverse_complement()


# Find positions of all open reading frames in given reading frame
def find_ORFs(nuc_seq, minimum_length, translation_table):
    if translation_table == 11:
        starts = set(["ATG", "GTG", "TTG"])
        stops = set(["TAA", "TAG", "TGA"])
    elif translation_table == 4:
        starts = set(["ATG", "GTG", "TTG"])
        stops = set(["TAA", "TAG"])
    else:
        printcolor("Translation table ", translation_table, " not implemented. Please open a GitHub issue if this is a problem.")
        sys.exit()
    ORF_startstop = []
    temp_starts = []
    for i in range(0, len(nuc_seq), 3):
        if i == 0 or nuc_seq[i:i+3] in starts:
            temp_starts.append(i)
            continue
        if ((nuc_seq[i:i+3] in stops) or (i+3 == len(nuc_seq))) and len(temp_starts) != 0:
            for start in temp_starts:
                if (i-start >= minimum_length):
                    ORF_startstop.append((start, i))
            temp_starts = []
    return ORF_startstop


def get_ORF_info(seq_list, min_orf_length, translation_table):
    ORF_seq = []
    ORF_coord = []
    ORF_nucseq = []
    for i, seq in enumerate(seq_list[:]):
        # frame 0: starts at 0 / frame 1: starts at 1 / frame 2: starts at 2 / frame r0: ends at 0, MAY NOT START AT THE LAST COORD DUE TO MULTIPLE OF 3 DIFFERENCES / frame r1: ends at 1 / frame r2: ends at 2
        seqstr = str(seq)
        seq_c = seq.complement()
        seqstr_c = str(seq_c)
        frame_0_end = (len(seqstr)-0)-(len(seqstr)-0) % 3+0
        frame_1_end = (len(seqstr)-1)-(len(seqstr)-1) % 3+1
        frame_2_end = (len(seqstr)-2)-(len(seqstr)-2) % 3+2
        frame_0 = find_ORFs(seqstr[0:frame_0_end], min_orf_length, translation_table)
        frame_1 = find_ORFs(seqstr[1:frame_1_end], min_orf_length, translation_table)
        frame_2 = find_ORFs(seqstr[2:frame_2_end], min_orf_length, translation_table)
        frame_r0 = find_ORFs(seqstr_c[0:frame_0_end][::-1], min_orf_length, translation_table)
        frame_r1 = find_ORFs(seqstr_c[1:frame_1_end][::-1], min_orf_length, translation_table)
        frame_r2 = find_ORFs(seqstr_c[2:frame_2_end][::-1], min_orf_length, translation_table)
        # standardize coords
        ORF_0f_standard_nuccoord = [(x[0]+3, x[1]) for x in frame_0]
        ORF_1f_standard_nuccoord = [(x[0]+4, x[1]+1) for x in frame_1]
        ORF_2f_standard_nuccoord = [(x[0]+5, x[1]+2) for x in frame_2]
        ORF_0r_standard_nuccoord = [(frame_0_end-x[1], frame_0_end-x[0]-3) for x in frame_r0]
        ORF_1r_standard_nuccoord = [(frame_1_end-x[1], frame_1_end-x[0]-3) for x in frame_r1]
        ORF_2r_standard_nuccoord = [(frame_2_end-x[1], frame_2_end-x[0]-3) for x in frame_r2]
        # translate once per frame, then slice
        aa_0 = str(seq[0:frame_0_end].translate(table=translation_table, to_stop=False))
        aa_1 = str(seq[1:frame_1_end].translate(table=translation_table, to_stop=False))
        aa_2 = str(seq[2:frame_2_end].translate(table=translation_table, to_stop=False))
        aa_r0 = str(seq_c[0:frame_0_end][::-1].translate(table=translation_table, to_stop=False))
        aa_r1 = str(seq_c[1:frame_1_end][::-1].translate(table=translation_table, to_stop=False))
        aa_r2 = str(seq_c[2:frame_2_end][::-1].translate(table=translation_table, to_stop=False))
        # reversed because model is trained with first amino acid directly upstream of stop codon
        ORF_0f_aa = [aa_0[slice(*tuple(int(idx/3) for idx in x))][::-1] for x in frame_0]
        ORF_1f_aa = [aa_1[slice(*tuple(int(idx/3) for idx in x))][::-1] for x in frame_1]
        ORF_2f_aa = [aa_2[slice(*tuple(int(idx/3) for idx in x))][::-1] for x in frame_2]
        ORF_0r_aa = [aa_r0[slice(*tuple(int(idx/3) for idx in x))][::-1] for x in frame_r0]
        ORF_1r_aa = [aa_r1[slice(*tuple(int(idx/3) for idx in x))][::-1] for x in frame_r1]
        ORF_2r_aa = [aa_r2[slice(*tuple(int(idx/3) for idx in x))][::-1] for x in frame_r2]
        ORF_seq.append([ORF_0f_aa, ORF_1f_aa, ORF_2f_aa, ORF_0r_aa, ORF_1r_aa, ORF_2r_aa])
        ORF_coord.append([ORF_0f_standard_nuccoord, ORF_1f_standard_nuccoord, ORF_2f_standard_nuccoord, ORF_0r_standard_nuccoord, ORF_1r_standard_nuccoord, ORF_2r_standard_nuccoord])
        ORF_nucseq.append([str(seq[0:frame_0_end]), str(seq[1:frame_1_end]), str(seq[2:frame_2_end]), str(seq_c[0:frame_0_end][::-1]), str(seq_c[1:frame_1_end][::-1]), str(seq_c[2:frame_2_end][::-1])])
    return ORF_seq, ORF_nucseq, ORF_coord


def analyze_overlap(coords0, coords1, strand0, strand1, unidirectional_penalty_per_base, convergent_penalty_per_base, divergent_penalty_per_base, max_gene_overlap):
    # TODO account for fully overlapped gene
    overlap = coords0[1] - coords1[0]
    if overlap <= 0:
        return True, 0
    if overlap > max_gene_overlap:
        return False, 0
    # get prime locations
    if strand0 == 1:
        threeprime0 = coords0[1]
        fiveprime0 = coords0[0]
    else:
        threeprime0 = coords0[0]
        fiveprime0 = coords0[1]
    if strand1 == 1:
        threeprime1 = coords1[1]
        fiveprime1 = coords1[0]
    else:
        threeprime1 = coords1[0]
        fiveprime1 = coords1[1]
    # exclude ORFs in same frame sharing same stop codon
    if strand0 == strand1 and threeprime0 == threeprime1:
        return False, 0
    # unidirectional overlap
    if (threeprime0 < fiveprime0) == (threeprime1 < fiveprime1):
        return True, overlap * unidirectional_penalty_per_base
    # convergent overlap
    if (fiveprime0 < threeprime1 <= threeprime0) or (fiveprime1 < threeprime0 <= threeprime1):
        return True, overlap * convergent_penalty_per_base
    # divergent overlap
    if (threeprime0 < fiveprime1 <= fiveprime0) or (threeprime1 < fiveprime0 <= fiveprime1):
        return True, overlap * divergent_penalty_per_base
    # edge case of exactly 1 ORF
    return True, 0


def predict(X, model):
    model.eval()
    with torch.no_grad():
        X_enc = F.one_hot(X, 21).permute(0, 2, 1).float()
        probs = expit(model(X_enc).cpu())
    return probs


def predict_tis(X, model_tis):
    model_tis.eval()
    with torch.no_grad():
        X_enc = F.one_hot(X, 4).permute(0, 2, 1).float()
        probs = expit(model_tis(X_enc).cpu())
    return probs


def tensor_to_seq(tensor):
    table = {0: "X", 1: "L", 2: "V", 3: "I", 4: "M", 5: "C", 6: "A", 7: "G", 8: "S", 9: "T", 10: "P", 11: "F", 12: "Y", 13: "W", 14: "E", 15: "D", 16: "N", 17: "Q", 18: "K", 19: "R", 20: "H"}
    return "".join([table[x] for x in tensor])


def kmerize(seq, k):
    kmerset = set()
    for i in range(len(seq) - k + 1):
        kmer = tuple(seq[i: i + k].tolist())
        kmerset.add(kmer)
    return kmerset


'''
-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------
                                          DISPLAY FUNCTIONS
-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------
'''


# ***** Print color ***** #
def printcolor(text: str, style: str = None, rgbFG: str = None, rgbBG: str = None, colorBool: bool = False):
    # decodeText = codecs.decode(text, "unicode_escape")
    if colorBool is True:
        if rgbBG != "None":
            print("\x1b["+str(style)+";38;2;"+rgbFG+";48;2;"+rgbBG+"m"+text+"\x1b[0m", end='')
        else:
            print("\x1b["+str(style)+";38;2;"+rgbFG+"m"+text+"\x1b[0m", end='')
    else:
        print(text, end='')
    # Terminal title
    if "ERROR" in text or "⛔" in text:
        title(text="⛔", pbar=None)
    elif "♊" in text:
        formatText = ""
        for word in text.replace("♊ ", "").replace("\n", "").split(" "):
            formatText += word.capitalize()
        title(text=formatText, pbar=None)


# ***** gemini title ***** #
def gemini_header():
    title = "\n┌─┐  ┌─┐  ┌┬┐  ┬  ┌┐┌  ┬\n│ ┬  ├┤   │││  │  │││  │\n└─┘  └─┘  ┴ ┴  ┴  ┘└┘  ┴\n"
    syntax = Syntax(title, "ada", theme="material", background_color="default", line_numbers=False)
    console = Console()
    console.print(syntax)


# ***** gemini function header ***** #
def gemini_fct_header(fctName):
    toDisplay = "# "+dicoGemini[fctName]['descr']+" (v"+dicoGemini[fctName]['version']+")\n"
    toDisplay += "def "+fctName+"("
    for arg in dicoGemini[fctName]['dicoArgs']:
        if arg == "pathTMP":
            continue
        if dicoGemini[fctName]['dicoArgs'][arg]['type'] == "str":
            toDisplay += "\n"+" "*(len(fctName)+5)+arg+"=\""+str(dicoGemini[fctName]['dicoArgs'][arg]['value'])+"\","
        else:
            toDisplay += "\n"+" "*(len(fctName)+5)+arg+"="+str(dicoGemini[fctName]['dicoArgs'][arg]['value'])+","
    toDisplay = toDisplay[:-1]+"\n"+" "*(len(fctName)+4)+")\n"
    syntax = Syntax(toDisplay, "python", theme="monokai", background_color="default", line_numbers=False)
    console = Console()
    console.print(syntax)


# ***** Show all gemini functions ***** #
def listfct():
    gemini_header()
    lstFctTag = ["sequence", "parser", "annotation", "cluster", "phage", "phylo", "plot", "download"]
    toDisplay = ""
    for fctTag in lstFctTag:
        toDisplay += "\"\"\" "+fctTag.upper()+" functions \"\"\"\n"
        for fctName in dicoGemini:
            if dicoGemini[fctName]['tag'] == fctTag:
                toDisplay += "# "+dicoGemini[fctName]['descr']+"\n"
                toDisplay += "def "+fctName+"("
                for arg in dicoGemini[fctName]['dicoArgs']:
                    if arg == "pathTMP":
                        continue
                    if dicoGemini[fctName]['dicoArgs'][arg]['value'] != "":
                        if dicoGemini[fctName]['dicoArgs'][arg]['type'] == "str":
                            toDisplay += "'"+dicoGemini[fctName]['dicoArgs'][arg]['tag']+"'"+arg+"=\""+str(dicoGemini[fctName]['dicoArgs'][arg]['value'])+"\", "
                        else:
                            toDisplay += "'"+dicoGemini[fctName]['dicoArgs'][arg]['tag']+"'"+arg+"="+str(dicoGemini[fctName]['dicoArgs'][arg]['value'])+", "
                    else:
                        toDisplay += "'"+dicoGemini[fctName]['dicoArgs'][arg]['tag']+"'"+arg+", "
                toDisplay = toDisplay[:-2]+")\n"
        toDisplay += "\n"
    syntax = Syntax(toDisplay[:-1], "python", theme="monokai", background_color="default", word_wrap=True, line_numbers=False)
    console = Console()
    console.print(syntax)


# ***** Show one gemini function ***** #
def infofct(fctName):
    gemini_header()
    for file in glob.glob(os.path.dirname(os.path.abspath(__file__))+"/*.py"):
        PY = open(file, 'r')
        data = PY.read()
        PY.close()
        if "def "+fctName in data:
            funcDoc = "'''"
            for docLine in data.split("def "+fctName)[1].split("'''")[1].split("\n"):
                if docLine != "" and docLine != "    ":
                    if "|" in docLine:
                        funcDoc += "\n"+docLine.strip()
                    else:
                        funcDoc += "\n "+docLine.strip()
            toDisplay = funcDoc+"\n'''\n"
    toDisplay += "def "+fctName+"("
    for arg in dicoGemini[fctName]['dicoArgs']:
        if arg == "pathTMP":
            continue
        if dicoGemini[fctName]['dicoArgs'][arg]['value'] != "":
            if dicoGemini[fctName]['dicoArgs'][arg]['type'] == "str":
                toDisplay += "'"+dicoGemini[fctName]['dicoArgs'][arg]['tag']+"'"+arg+"=\""+str(dicoGemini[fctName]['dicoArgs'][arg]['value'])+"\", "
            else:
                toDisplay += "'"+dicoGemini[fctName]['dicoArgs'][arg]['tag']+"'"+arg+"="+str(dicoGemini[fctName]['dicoArgs'][arg]['value'])+", "
        else:
            toDisplay += "'"+dicoGemini[fctName]['dicoArgs'][arg]['tag']+"'"+arg+", "
    toDisplay = toDisplay[:-2]+")\n"
    toDisplay += "\n"
    syntax = Syntax(toDisplay[:-1], "python", theme="monokai", background_color="default", word_wrap=True, line_numbers=False)
    console = Console()
    console.print(syntax)


# ***** Terminal title ***** #
def title(text: str, pbar: tqdm = None):
    # Reformat title text to block tab size
    formatText = geminiset.title
    if text != "":
        formatText += " | "
        for word in text.replace("♊ ", "").replace("\n", "").split(" "):
            formatText += word.capitalize()
    if pbar is not None:
        formatText += " ("+str(int(round(((pbar.format_dict['n']*100)/pbar.format_dict['total']), 0)))+"%)"
    if len(formatText) > 56:
        justifyText = "."+formatText[:56]+"."
    elif (58-len(formatText)) % 2 == 0:
        justifyText = ".".ljust(int((58-len(formatText))/2))+formatText+".".rjust(int((58-len(formatText))/2))
    else:
        justifyText = ".".ljust(int(((58-len(formatText))/2)+1))+formatText+".".rjust(int((58-len(formatText))/2))
    print("\33]0;"+justifyText+"\a", end="", flush=True)


# ***** Pretty dictionnary print ***** #
def printDico(dico, indent=4):
    pp = pprint.PrettyPrinter(indent=indent)
    pp.pprint(dico)


# ***** Generate random hex color ***** #
def random_hex_color(uppercase=False):
    hex_parts = '0123456789abcdef'
    color = '#'
    for i in range(6):
        color += hex_parts[floor(random() * 16)]
    if uppercase is True:
        return color.upper()
    else:
        return color


# ***** Color manipulation ***** #
def hex_to_RGB(hex):
    # return list
    return [int(hex[i:i+2], 16) for i in range(1, 6, 2)]


def RGB_to_hex(RGB):
    RGB = [int(x) for x in RGB]
    return "#"+"".join(["0{0:x}".format(v) if v < 16 else "{0:x}".format(v) for v in RGB])


def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
    s = hex_to_RGB(start_hex)
    f = hex_to_RGB(finish_hex)
    HEX_list = [start_hex]
    RBG_list = [s]
    for t in range(1, n):
        curr_vector = [
          int(s[j] + (float(t)/(n-1))*(f[j]-s[j]))
          for j in range(3)
        ]
        RBG_list.append(curr_vector)
        HEX_list.append(RGB_to_hex(curr_vector))
    return HEX_list, RBG_list
