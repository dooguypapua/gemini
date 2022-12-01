<p style="text-align: center;">![gemini logo](conf/gemini_icon.png)</p>

```python
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