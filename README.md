![gemini logo](conf/gemini_icon.png)

```abap
python gemini.py --help
```

```mupad
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '\n'
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'pathIN, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
# Search term in FASTA folder or file and return sequence
def search_in_fasta('-i'pathIN, '-term'searchTerm, '-o'pathOUT, '-nodoublon'boolDoublon=False, '-ext'ext=".faa")
# Check if a sequence is circular
def check_circular('-i'seqIN, '-minlen'minLen=3)
# Create a table with the size (pb) of genomes
def fasta_genome_size('-i'pathIN, '-o'pathOUT, '-s'boolSort=True, '-ext'ext=".fna")
# Retrieve genome sequence in a GBK file and create a FNA
def gbk_to_fna('-i'pathIN, '-o'pathOUT)
# Retrieve genome sequence in a GBK file and create a FFN
def gbk_to_ffn('-i'pathIN, '-o'pathOUT, '-syntaxic'syntaxic="prodigal")
# Retrieve protein sequence in a GBK file and create a FAA
def gbk_to_faa('-i'pathIN, '-o'pathOUT, '-syntaxic'syntaxic="prodigal")
# Retrieve protein sequence in a GBK file and create a FAA with annotation
def gbk_to_annotFAA('-i'pathIN, '-o'pathOUT)
# Retrieve protein sequence in a GBK file and create all
def gbk_to_all('-i'pathIN, '-o'pathOUT, '-syntaxic'syntaxic="prodigal")
# Retrieve protein sequence in a GBK file and create a GFF
def gbk_to_gff('-i'pathIN, '-o'pathOUT)
# Create a annotation table from a GFF file
def gff_to_table('-i'pathIN, '-o'pathOUT, '-format'format=".xlsx", '-width'maxWidth=50, '-ext'ext=".gff")

""" PARSER functions """
# Parse GBK files and create a dictionnary
def make_fasta_dict('-i'pathIN, '-ltheader'onlyLTasHeader=False, '-j'pathJSON="None")
# Parse GBK files and create a dictionnary
def make_gbk_dict('-i'pathIN, '-j'pathJSON="None", '-sort'boolSort=True, '-pseudo'boolPseudo=False)
# Parse GFF3 files and create a dictionnary
def make_gff_dict('-i'pathIN, '-j'pathJSON="None", '-ext'ext=".fna")
# Make a GBK file from FASTA (.ffn + .faa + .fna)
def make_gbk_from_fasta('-i1'pathIN1, '-i2'pathIN2, '-i3'pathIN3, '-o'pathOUT, '-identifier'identifier, '-topo'topology, '-div'division, '-taxid'taxID=0, '-i4'pathIN4="None", '-progress'boolProgress=True)
# Parse blast output and create a dictionnary
def make_blast_dict('-i'pathIN, '-j'pathJSON="None", '-pid'idThr=20, '-minlr'minLRthr=50, '-maxlr'maxLRthr=50, '-ext'ext=".xml")
# Parse HMMSCAN tblout output and create a dictionnary
def make_hmmscan_dict('-i'pathIN, '-j'pathJSON="None", '-e'idEvalue=0.01, '-ext'ext=".tblout")
# Parse pVOGs table files and create a dictionnary
def make_pvogs_desc_dict('-i'pathIN="None", '-j'pathJSON="None")
# Parse tRNAscan-SE output and create a dictionnary
def make_trnascanse_dict('-i'pathIN, '-j'pathJSON="None", '-ext'ext=".trnascanse")
# Parse EggNOG output and create a dictionnary
def make_eggnog_dict('-i'pathIN, '-j'pathJSON="None", '-ext'ext=".annotations")
# Parse InterProScan output and create a dictionnary
def make_interpro_dict('-i'pathIN, '-e'idEvalue=0.01, '-j'pathJSON="None", '-ext'ext=".tsv")
# Reformat ORF name in Phanotate output FASTA
def reformat_phanotate('-i'pathIN)
# Reformat organism name in PanACoTA output files
def reformat_panacota('-i'pathIN)
# Parse MMSEQS cluster tsv file and create a dictionnary
def make_mmseqs_cluster_dict('-i'pathIN, '-j'pathJSON="None")

""" ANNOTATION functions """
# Phanotate syntaxic annotation from phage genome files
def phanotate('-i'pathIN, '-o'pathOUT, '-len'minLen=0, '-bool'fromPhageDb=False, '-ext'ext=".fna")
# Balrog syntaxic annotation
def balrog('-i'pathIN, '-o'pathOUT, '-topo'topology, '-div'division, '-taxid'taxID=0, '-len'minLen=30, '-mmseqs'boolMmseqs=True, '-ext'ext=".fna")
# Translate gene CDS to protein (FFN>FAA)
def transeq('-i'pathIN, '-o'pathOUT, '-orgheader'boolOrgName=True, '-phagedb'fromPhageDb=False, '-ext'ext=".ffn")
# Launch tRNAscan-SE to predict tRNA genes
def trnascan_se('-i'pathIN, '-o'pathOUT, '-model'model="-B", '-ext'ext=".fna")
# diamond blastP
def diamond_p('-i'pathIN, '-d'pathDB, '-o'pathOUT, '-addseq'boolSeq=False, '-ext'ext=".faa")
# Annotate protein FASTA with InterProScan
def interproscan('-i'pathIN, '-o'pathOUT, '-ext'ext=".faa")
# Annotate protein FASTA with EggNOG
def eggnog('-i'pathIN, '-o'pathOUT, '-pid'idThr=20, '-cov'covThr=50, '-ext'ext=".faa")
# hmmscan against pVOGS profiles database
def pvogs('-i'pathIN, '-o'pathOUT, '-ext'ext=".faa")
# hmmscan against VirFam/PFAM recombinase profiles database
def recombinase('-i'pathIN, '-o'pathOUT, '-ext'ext=".faa")
# Search defense systems using DefenseFinder & PADLOC
def defense_system('-i'pathIN, '-o'pathOUT, '-dfmodel'dfmodelV="1.1.0", '-plmodel'plmodelV="1.4.0", '-ext'ext=".faa")
# Search phage satellites using SatelliteFinder
def satellite_finder('-i'pathIN, '-o'pathOUT, '-model'model="ALL", '-ext'ext=".faa")

""" CLUSTER functions """
# MMSEQS easycluster: sensitive clustering
def mmseqs_easycluster('-i'pathIN, '-o'pathOUT, '-pid'idThr=30, '-maxlr'maxLRthr=80, '-ext'ext=".faa")
# MMSEQS easyrbh: find reciprocal best hit
def mmseqs_rbh('-i'pathIN, '-o'pathOUT, '-ref'reference="None", '-pid'idThrClust=80, '-cov'covThrClust=80, '-nucl'boolNucl=False, '-ext'ext=".faa")
# RBH organism clustering and create a dictionnary
def make_rbhcluster_dict('-i'pathIN, '-i2'pathIN2, '-j'pathJSON, '-pid'idThrClust=80, '-cov'covThrClust=80, '-ext'ext=".rbh", '-ext2'ext2=".faa")
# Construct group core alignment from rbh clusters
def make_group_core_align('-i'pathIN, '-j'pathJSON, '-group'pathGROUP, '-o'pathOUT, '-gene'boolGene, '-prot'boolProt, '-extN'extN=".ffn", '-extP'extP=".faa")
# Compute pan, core, variable genome by organism group
def pan_core_group('-j'pathJSON, '-i'pathDIST, '-group'pathGROUP, '-o'pathOUT="stdout")
# Construct protein core from vibrio faa files
def make_vibrio_core('-i'pathIN, '-i2'pathIN2, '-o'pathOUT, '-ref'ref, '-pid'idThrClust=30, '-cov'covThrClust=80, '-mash'maxMash=0.3, '-cut'cutN=5, '-maxcontig'maxContig=1000, '-l90'maxL90=100,
'-persRatio'persRatio=0.9, '-minPersPart'minPersPart=0.75, '-minsize'minsize=1.0, '-ext'ext=".faa")
# PPanGGOLiN RGP analysis
def ppanggolin('-i'pathIN, '-i2'pathIN2, '-o'pathOUT, '-maxrgp'maxRGP=-1, '-prefix'prefix="None", '-ext'ext=".gbk.gz")

""" PHAGE functions """
# Phages clustering using VIRIDIC
def viridic('-i'pathIN, '-o'pathOUT, '-ext'ext=".fna")
# Search terminase from FAA files
def search_terminase('-i'pathIN, '-d'pathDMND, '-o'pathOUT, '-j'pathJSON="None", '-pid'idThr=20, '-minlr'minLRthr=50, '-maxlr'maxLRthr=50, '-ext'ext=".faa")
# Phage syntaxic and functionnal annotation
def phage_annotation('-i'pathIN, '-o'pathOUT, '-embl'boolEMBL=False, '-project'enaProject="None", '-taxo'pathTAXO="None", '-e'idEvalue=0.01, '-pid'idThr=30, '-cov'covThr=50, '-pid2'idThrClust=80,
'-cov'covThrClust=80, '-ext'ext=".fna")
# Make GenBank phages database
def phageDB('-i'pathIN, '-o'pathOUT, '-checkv'checkvHQ=75.0)
# Modified version of VIRIDIC
def myVIRIDIC('-i'pathIN, '-o'pathOUT, '-thfam'thfam=50.0, '-thgen'thgen=70.0, '-thsp'thsp=95.0, '-ext'ext=".fna")
# PhiSpy prophages prediction
def PhiSpy('-i'pathIN, '-o'pathOUT, '-nb'nbAdjacent=3, '-len'minCtgLen=5000, '-ext'ext=".gbk")
# Search PICMI from gbk files
def picmi_finder_gbk('-i'pathIN, '-o'pathOUT, '-prefix'prefix, '-len'maxLen=50000)
# Search PICMI from databank SEQ files
def picmi_finder_databankseq('-i'pathIN, '-o'pathOUT, '-len'maxLen=50000)

""" PHYLO functions """
# Mash distance matrix
def mash_matrix('-i'pathIN, '-o'pathOUT, '-sketch'sketchSize=10000, '-ext'ext=".fna")
# Create fastANI JSON database
def fastani_db('-i'pathIN, '-i2'pathIN2, '-j'pathJSON, '-fragLen'fragLen=3000, '-ext'ext=".fna")
# Construct tree and search best group topology
def best_gene_tree_topology('-i'pathIN1, '-i2'pathIN2, '-i3'pathIN3, '-o'pathOUT, '-outgrp'outgroup, '-pid'idThrClust=80, '-cov'covThrClust=80, '-ext1'ext1=".ffn", '-ext2'ext2=".faa")
# Search kmers specific to a organism group
def specific_kmers('-i'pathIN, '-i2'pathIN2, '-o'pathOUT, '-len'kmerLen=25, '-ext'ext=".fna")
# Make core proteins tree
def core_prot_tree('-i'pathIN, '-o'pathOUT, '-pid'idThr=30, '-cov'covThr=80, '-ext'ext=".faa")
# Make individual and core genes tree
def genes_tree('-i'pathIN, '-o'pathOUT, '-pid'idThr=30, '-cov'covThr=80, '-ext'ext=".ffn")
# Create similarity matrix for one protein
def protein_similarity_matrix('-i'pathIN, '-o'pathOUT, '-lt'locusTag, '-pid'idThr=30, '-cov'covThr=80, '-ext'ext=".faa")
# Make flexible genes tree form PanACoTA results
def panacota_flexible_tree('-i'pathIN, '-o'pathOUT, '-filter'filterOrg="None")
# Call pairwise SNPs using Snippy
def snippy('-i'pathIN, '-o'pathOUT, '-ext'ext=".gbk")
# wGRR distance matrix
def wgrr_matrix('-i'pathIN, '-o'pathOUT, '-ext'ext=".faa")

""" PLOT functions """
# Create linear gene plot from GFF3 file
def gff_to_linear_geneplot('-i'pathIN, '-o'pathOUT, '-lt'pathLT="None", '-len'length=-1, '-ext'ext=".gff")
# Create linear gene plot from GFF3 file and merge per group
def gff_to_linear_group_geneplot('-i'pathIN, '-cluster'pathCLUSTER, '-group'pathGROUP, '-o'pathOUT, '-ext'ext=".gff")
# Transform SVG plot from dna_features_viewer
def svg_dna_transform('-i'pathIN, '-o'pathOUT)
# Convert xlsx table to heatmap table
def xlsx_to_heatmap('-i'pathIN, '-o'pathOUT, '-cstart'colorStart="FFFFFF", '-cend'colorEnd="FF0000", '-row'headerRow=-1, '-col'headerCol=-1)
# Circular circos plot from genbank file
def circos_plot('-i'pathIN, '-o'pathOUT, '-i2'pathIN2="None", '-pid'pident=30, '-cov'cov=80)
# Circular circos plot for alignment
def circos_align('-i'pathIN, '-o'pathOUT)

""" DOWNLOAD functions """
# Download all bacteria Genbank files
def dl_genbank_bacteria('-section'section, '-tax'taxonomyID, '-o'pathOUT, '-chunk'chunkSize=100)
```

---------------------------------------

#### python3 requirements
- ast<br/>
- biopython<br/>
- cairosvg<br/>
- dna_features_viewer<br/>
- ete3<br/>
- matplotlib<br/>
- numpy<br/>
- openpyxl<br/>
- pandas<br/>
- psutil<br/>
- requests<br/>
- rich<br/>
- scipy<br/>
- seaborn<br/>
- svgpathtools<br/>
- tabulate<br/>
- torch<br/>
- tqdm<br/>
- xlsxwriter<br/>
- yaspin<br/>

```Cucumber
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'Cucumber, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```abap
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'abap, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```ada
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'ada, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```ahk
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'ahk, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```apacheconf
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'apacheconf, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```applescript
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'applescript, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```as
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'as, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```as3
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'as3, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```asy
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'asy, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```bash
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'bash, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```bat
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'bat, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```befunge
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'befunge, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```blitzmax
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'blitzmax, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```boo
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'boo, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```brainfuck
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'brainfuck, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```c
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'c, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```cfm
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'cfm, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```cheetah
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'cheetah, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```cl
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'cl, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```clojure
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'clojure, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```cmake
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'cmake, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```coffeescript
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'coffeescript, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```console
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'console, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```control
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'control, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```cpp
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'cpp, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```csharp
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'csharp, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```css
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'css, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```cython
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'cython, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```d
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'd, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```delphi
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'delphi, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```diff
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'diff, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```dpatch
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'dpatch, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```duel
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'duel, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```dylan
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'dylan, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```erb
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'erb, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```erl
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'erl, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```erlang
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'erlang, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```evoque
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'evoque, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```factor
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'factor, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```felix
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'felix, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```fortran
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'fortran, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```gas
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'gas, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```genshi
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'genshi, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```glsl
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'glsl, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```gnuplot
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'gnuplot, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```go
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'go, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```groff
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'groff, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```haml
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'haml, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```haskell
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'haskell, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```html
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'html, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```hx
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'hx, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```hybris
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'hybris, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```ini
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'ini, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```io
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'io, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```ioke
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'ioke, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```irc
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'irc, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```jade
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'jade, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```java
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'java, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```js
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'js, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```jsp
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'jsp, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```lhs
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'lhs, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```llvm
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'llvm, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```logtalk
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'logtalk, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```lua
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'lua, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```make
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'make, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```mako
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'mako, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```maql
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'maql, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```mason
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'mason, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```markdown
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'markdown, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```modelica
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'modelica, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```modula2
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'modula2, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```moocode
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'moocode, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```mupad
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'mupad, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```mxml
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'mxml, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```myghty
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'myghty, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```nasm
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'nasm, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```newspeak
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'newspeak, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```objdump
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'objdump, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```objectivec
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'objectivec, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```objectivej
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'objectivej, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```ocaml
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'ocaml, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```ooc
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'ooc, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```perl
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'perl, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```php
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'php, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```postscript
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'postscript, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```pot
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'pot, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```pov
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'pov, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```prolog
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'prolog, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```properties
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'properties, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```protobuf
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'protobuf, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```py3tb
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'py3tb, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```pytb
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'pytb, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```python
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'python, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```rb
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'rb, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```rconsole
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'rconsole, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```rebol
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'rebol, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```redcode
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'redcode, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```rhtml
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'rhtml, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```rst
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'rst, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```sass
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'sass, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```scala
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'scala, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```scaml
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'scaml, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```scheme
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'scheme, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```scss
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'scss, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```smalltalk
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'smalltalk, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```smarty
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'smarty, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```sourceslist
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'sourceslist, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```splus
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'splus, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```sql
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'sql, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```sqlite3
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'sqlite3, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```squidconf
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'squidconf, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```ssp
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'ssp, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```tcl
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'tcl, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```tcsh
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'tcsh, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```tex
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'tex, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```text
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'text, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```v
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'v, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```vala
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'vala, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```vbnet
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'vbnet, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```velocity
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'velocity, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```vim
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'vim, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```xml
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'xml, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```xquery
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'xquery, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```xslt
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'xslt, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```
```yaml
""" SEQUENCE functions """
# Unwrap FASTA by removing sequence '
def unwrap_fasta('-i'pathIN, '-ext'ext=".fna")
# Split FASTA sequence
def split_fasta('-i'yaml, '-o'pathOUT, '-term'term="N", '-ext'ext=".fna")
```