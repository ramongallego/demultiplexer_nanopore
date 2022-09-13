#!/usr/bin/env bash


################################################################################
# INPUT
################################################################################
# What is the file path to the directory containing all of the libraries/reads?

PARENT_DIR="${MAIN_DIR}"/data_Dani

# Where is the sequencing metadata file? (SEE FORMATTING GUIDELINES IN README!)
SEQUENCING_METADATA="${PARENT_DIR}"/metadata_Run003.csv



################################################################################
# OUTPUT
################################################################################
# This script will generate a directory (folder) containing the output of the script.
# Where do you want this new folder to go?
OUTPUT_DIRECTORY=/mnt/c/Users/RG.5015511/Documents/Projects/demultiplexer_nanopore/test_demult #"${PARENT_DIR%/*}"


################################################################################
# METADATA DETAILS
################################################################################
# Specify columns for raw sequencing files:
COLNAME_FILE1="file"
#COLNAME_FILE2="file2"

# MUST be unique for each row!
COLNAME_SAMPLE_ID="sample_id"


# Your metadata must have a column corresponding to the subfolders containing the raw reads.
# In order to make this flexible across both multiple and single library preps, you must include this even if you only sequenced one library (sorry!).
COLNAME_ID1_NAME="plate_name.p5"
COLNAME_ID1_SEQ="barcode.p5"

#################################################################################
# MULTI:CORE
################################################################################

# In a Linux system works much faster if you run cutadapt with multicore support
# Type here the number of cores you want to dedicate to the task


N_CORES=8


################################################################################
# DEMULTIPLEXING
################################################################################

# Do the reads contain index sequences which identifies their sample of origin?
SECONDARY_INDEX="YES"

# Specify the nucleotide sequences that differentiate multiplexed samples
# (sometimes, confusingly referred to as "tags" or "barcodes")
# these are the secondary index -- the primary index added with the sequencing adapters should not be in the sequence data
# You can grab these from the file specified above (SEQUENCING_METADATA) by specifying the column name of index sequences.
COLNAME_ID2_SEQ="barcode.p7"
COLNAME_ID2_WELL="Well.p7"

# How many nucleotides pad the 5' end of the tag sequence?
# TODO build in flexibility (this number is unused right now)


################################################################################
# PRIMER REMOVAL
################################################################################
# Specify the primers used to generate these amplicons.
# As with the multiplex indexes, Banzai will grab these from the file SEQUENCING_METADATA.
# You must indicate the column names of the forward and reverse primers
COLNAME_PRIMER1="PrimerF"
COLNAME_PRIMER2="PrimerR"
COLNAME_LOCUS="Locus"

################################################################################
# USE HASH
################################################################################
# Should the sequence ID after dereplication be the output of a hash algorithm?

USE_HASH="YES"

################################################################################
# CLUSTER OTUs: USING decona or NGspeciesID
################################################################################

SEARCH_ASVs="decona"


## TODO: Add variables to control the behaviour of decona

MIN_LENGTH="800"

MAX_LENGTH="1900"

CLUSTER_SIM="0.8"

## ADDED on 20220909
## DO you want the final fastqs to also be split by locus (yes) or do you want just
## one fastq per index (Plate_well)
DEMULT_BY_PRIMER="YES"


################################################################################
# GENERAL SETTINGS
################################################################################
# Would you like to save every single intermediate file as we go? YES | NO
# recommendation: NO, unless testing or troubleshooting
HOARD="YES"
