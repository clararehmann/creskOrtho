#!/bin/bash

###  shell script to download genome and peptide data from refseq
###
###  INPUT: space-separated list of species accession names
###     bash download_refseq.sh Danio_rerio Homo_sapiens Anopheles_gambiae
###  OUTPUT: refseq protein and genome data for each species
###      out/transcripts/Danio_rerio/Danio_rerio/annotation/Danio_rerio_pep.faa.gz
###      out/transcripts/Danio_rerio/Danio_rerio/annotation/Danio_rerio_genomic.gff.gz

species=("$@")

find_refseq () {
    
    spID=$(echo $1 | sed "s/_//g" )

    # find most recent assembly accession
    cout=$(curl -s ftp.ncbi.nlm.nih.gov/genomes/refseq/{archaea,bacteria,fungi,invertebrate,mitochondrion,plant,plasmid,plastid,protozoa,vertebrate_mammalian,vertebrate_other,viral}/$1\/latest_assembly_versions/)
    assembly=${cout##*all}

    mkdir -p out/transcripts/$spID\/$spID\/annotation

    # download genomic gff
    wget -O out/transcripts/$spID\/$spID\/annotation/$spID\_genomic.gff.gz ftp.ncbi.nlm.nih.gov/genomes/all$assembly/${cout##*/}_genomic.gff.gz 
    # download peptide gff
    wget -O out/transcripts/$spID\/$spID\/annotation/$spID\_pep.faa.gz ftp.ncbi.nlm.nih.gov/genomes/all$assembly/${cout##*/}_translated_cds.faa.gz
}

for f in ${species[@]}; do find_refseq $f; done
