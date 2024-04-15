_geminicomplementation()
{
    local cur prev module

    if [ $COMP_CWORD -gt 3 ]
        then
        if [ "$(( $COMP_CWORD % 2 ))" -eq 0 ]
            then
            prev=${COMP_WORDS[1]}
        else
            prev=${COMP_WORDS[$COMP_CWORD-1]}        
        fi
        cur=${COMP_WORDS[COMP_CWORD]}
    else
        cur=${COMP_WORDS[COMP_CWORD]}
        prev=${COMP_WORDS[COMP_CWORD-1]}
    fi

    # Read previous used arguments
    readarray -d " " -t tmpAlreadyOptions<<<"${COMP_LINE}"
    alreadyOptions=()
    for arg in "${tmpAlreadyOptions[@]}"
        do 
        if [[ $arg == -[a-z]* ]]
        then
            alreadyOptions+=(${arg})
        fi
    done

    # Switch case
    case ${COMP_CWORD} in
        1)
            COMPREPLY=($(compgen -W "unwrap_fasta split_fasta search_in_fasta check_circular fasta_genome_size gbk_to_fna gbk_to_ffn gbk_to_faa gbk_to_annotFAA gbk_to_all gbk_to_gff gff_to_table make_fasta_dict make_gbk_dict make_gff_dict make_gbk_from_fasta slice_genes_genbank make_blast_dict make_hmmscan_dict make_pvogs_desc_dict make_trnascanse_dict make_eggnog_dict make_interpro_dict reformat_phanotate reformat_panacota make_mmseqs_cluster_dict phanotate balrog transeq trnascan_se diamond_p interproscan eggnog pvogs recombinase defense_system satellite_finder phage_assembly filter_phage_assembly replicon_distribution viridic phage_annotation phageDB phageDBsearch myVIRIDIC PhiSpy picmi_finder_gbk picmi_finder_databankseq mmseqs_easycluster mmseqs_cluster mmseqs_rbh make_rbhcluster_dict make_group_core_align pan_core_group make_vibrio_core ppanggolin gff_to_linear_geneplot rbh_linear_plot gff_to_linear_group_geneplot svg_dna_transform xlsx_to_heatmap circos_plot circos_align circos_rbh_plot mash_matrix fastani_db best_gene_tree_topology specific_kmers core_prot_tree individual_core_tree protein_similarity_matrix panacota_flexible_tree snippy wgrr_matrix dl_genbank_bacteria" -- ${cur}))
            ;;
        2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 | 22 | 24 | 26 | 28 | 30 | 32 | 34 | 36 | 38 | 40)
            case ${prev} in
                unwrap_fasta)
                options="-i -ext"
                ;;
                split_fasta)
                options="-i -o -term -ext"
                ;;
                search_in_fasta)
                options="-i -o -term -nodoublon -ext"
                ;;
                check_circular)
                options="-i -minlen"
                ;;
                fasta_genome_size)
                options="-i -o -s -ext"
                ;;
                gbk_to_fna | gbk_to_annotFAA | gbk_to_gff | phage_assembly | svg_dna_transform)
                options="-i -o"
                ;;
                gbk_to_ffn | gbk_to_faa | gbk_to_all)
                options="-i -o -syntaxic -split"
                ;;
                gff_to_table)
                options="-i -o -format -width -ext"
                ;;
                make_fasta_dict)
                options="-i -ltheader -j"
                ;;
                make_gbk_dict)
                options="-i -j -sort"
                ;;
                make_gbk_from_fasta)
                options="-i1 -i2 -i3 -i4 -o -id -topo -div -taxid -progress"
                ;;
                slice_genes_genbank)
                options="-i -o -lt1 -lt2"
                ;;                
                make_blast_dict)
                options="-i -j -pid -minlr -maxlr -ext"
                ;;
                make_hmmscan_dict)
                options="-i -j -e -ext"
                ;;
                make_pvogs_desc_dict | make_mmseqs_cluster_dict)
                options="-i -j"
                ;;
                make_trnascanse_dict | make_eggnog_dict | make_gff_dict)
                options="-i -j -ext"
                ;;
                make_interpro_dict)
                options="-i -e -j -ext"
                ;;
                reformat_phanotate | reformat_panacota)
                options="-i"
                ;;
                phanotate)
                options="-i -o -len -ext"
                ;;
                balrog)
                options="-i -o -topo -div -taxid -len -mmseqs -ext"
                ;;
                transeq)
                options="-i -o -orgheader -phagedb -ext"
                ;;
                trnascan_se)
                options="-i -o -model -ext"
                ;;
                diamond_p)
                options="-i -d -o -addseq -ext"
                ;;
                interproscan | pvogs | recombinase | viridic | snippy | wgrr_matrix)
                options="-i -o -ext"
                ;;
                eggnog)
                options="-i -o -pid -cov -ext"
                ;;
                defense_system)
                options="-i -o -pid -dfmodel -plmodel -ext"
                ;;
                satellite_finder)
                options="-i -o -model -ext"
                ;;
                filter_phage_assembly)
                options="-i -o -minlen -mincov -ext"
                ;;
                replicon_distribution)
                options="-i -ref -o -pid -ext1 -ext2"
                ;;
                phage_annotation)
                options="-i -o -embl -project -taxo -e -pid -cov -pid2 -ext"
                ;;
                phageDB)
                options="-i -o -checkv"
                ;;
                phageDBsearch)
                options="-i -o -pid -cov -db"
                ;;
                myVIRIDIC)
                options="-i -o -ref -thfam -thgen -thsp -db -ext"
                ;;
                PhiSpy)
                options="-i -o -nb -len -ext"
                ;;
                picmi_finder_gbk)
                options="-i -o -prefix -len"
                ;;
                picmi_finder_databankseq)
                options="-i -o -len"
                ;;
                mmseqs_easycluster)
                options="-i -o -pid -maxlr -ext"
                ;;
                mmseqs_cluster)
                options="-i -o -pid -cov -ext"
                ;;
                mmseqs_rbh)
                options="-i -o -ref -pid -cov -nucl -ext"
                ;;
                make_rbhcluster_dict)
                options="-i -i2 -j -pid -cov -ext -ext2"
                ;;
                make_group_core_align)
                options="-i -j -group -o -gene -prot -extN -extP"
                ;;
                pan_core_group)
                options="-j -i -group -o"
                ;;
                make_vibrio_core)
                options="-i -i2 -o -ref -pid -cov -mash -cut -maxcontig -l90 -persRatio -minPersPart -minsize -ext"
                ;;
                ppanggolin)
                options="-i -i2 -o -maxrgp -prefix -ext"
                ;;
                gff_to_linear_geneplot)
                options="-i -o -lt -len -ext"
                ;;
                rbh_linear_plot)
                options="-i -cluster -o -distinctcolor -i2 -ext"
                ;;                
                gff_to_linear_group_geneplot)
                options="-i -cluster -group -o -ext"
                ;;
                xlsx_to_heatmap)
                options="-i -o -cstart -cend -row -col"
                ;;
                circos_plot)
                options="-i -o -i2 -pid -cov"
                ;;
                circos_align)
                options="-i -o"                
                ;;
                circos_rbh_plot)
                options="-i -j -o -maxcore"                
                ;;
                mash_matrix)
                options="-i -o -sketch -ext"
                ;;
                fastani_db)
                options="-i -i2 -j -fragLen -ext"
                ;;
                best_gene_tree_topology)
                options="-i -i2 -i3 -o -outgrp -pid -cov -ext1 -ext2"
                ;;
                specific_kmers)
                options="-i -i2 -o -len -ext"
                ;;
                core_prot_tree | individual_core_tree)
                options="-i -o -pid -cov -nucl -ext"
                ;;
                protein_similarity_matrix)
                options="-i -o -lt -pid -cov -ext"
                ;;
                panacota_flexible_tree)
                options="-i -o -filter"
                ;;
                dl_genbank_bacteria)
                options="-section -tax -o -chunk"
                ;;
            esac
            # Replace already done arguments
            for value in "${alreadyOptions[@]}"
                do options=${options/"${value}"}
            done
            # Apply words to comgen
            COMPREPLY=($(compgen -W "${options}" -- ${cur}))
            ;;
        3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | 21 | 23 | 25 | 27 | 29 | 31 | 33 | 35 | 37 | 39 | 41)
            case ${prev} in
                # PATHS
                -i | -o | -j | -i1 | -i2 | -i3 | -i4 | -d)
                COMPREPLY=($(compgen -f -- ${cur}))
                ;;
                # EXTENSION FILTERS
                -extN | -extP | -ext | -ext1 | -ext2 | -format)
                COMPREPLY=($(compgen -W ".annotations .faa .faa.gz .fasta .fasta.gz .ffn .ffn.gz .fna .fna.gz .gbk .gbk.gz .gff .rbh .tblout .trnascanse .tsv .xml .xlsx .xls" -- ${cur}))
                ;;
                # BOOLEAN
                -mmseqs | -ltheader | -nodoublon | -db | -sort | -orgheader | -phagedb | -nucl | -gene | -prot | -embl | -addseq | -split | -progress)
                COMPREPLY=($(compgen -W "true false" -- ${cur}))
                ;;                
                # INTEGER
                -minlen | -width | -pid | -minlr | -maxlr | -len | -cov | -mincov | -pid2 | -nb | -taxid | -cut | -maxcontig | -l90 | -maxrgp | -row | -col | -sketch | -fragLen | -tax | -chunk)
                COMPREPLY=($(compgen -W "<int>" -- ${cur}))
                ;;
                # FLOAT
                -e | -checkv | -thfam | -thgen | -thsp | -mash | -persRatio | -minPersPart | -minsize)
                COMPREPLY=($(compgen -W "<float>" -- ${cur}))
                ;;
                # STRING
                -identifier | -term | -ref | -group | -prefix | -lt | -lt1 | -lt2 | -cluster | -cstart | -cend | -outgrp | -filter | -project | -taxo | -dfmodel | -plmodel)
                COMPREPLY=($(compgen -W "<str>" -- ${cur}))
                ;;
                # PREDEFINED CHOICES
                -syntaxic)
                COMPREPLY=($(compgen -W "prodigal" -- ${cur}))
                ;;
                -topo)
                COMPREPLY=($(compgen -W "linear circular" -- ${cur}))
                ;;
                -div)
                COMPREPLY=($(compgen -W "BCT PHG" -- ${cur}))
                ;;
                -model)
                COMPREPLY=($(compgen -W "-E -B -A -O -G" -- ${cur}))
                ;;
                -section)
                COMPREPLY=($(compgen -W "refseq genbank" -- ${cur}))
                ;;
            esac
            ;;
        *)
            COMPREPLY=()
            ;;
    esac
}

# complete -o nosort -o nospace -F _geminicomplementation gemini

