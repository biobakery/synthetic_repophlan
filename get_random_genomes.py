#! /usr/bin/env python

"""
Author: Timothy Tickle
Summary: Randomly select reference genomes (fna) from a tar.gz and make input files for synthetic community creation.
"""

import argparse
from Bio import SeqIO
import csv
import glob
import random
import os
import shutil
import string
import subprocess

c_sTaxonTableFileName = "TaxonTable"
c_sCommunityPercentFileName = "Input"
c_sFileExtension  = "*.fna"
c_sGenusKey = "g__"
c_sGoldStandardFileName = "GoldStandard"
c_sKingdomKey = "k__"
c_sLognormalMu = 1
c_sLognormalSigma = 2.7
c_sSpeciesKey = "s__"
c_sStrainKey = "t__"

argp = argparse.ArgumentParser( prog = "get_random_genomes.py", description = """Given a supset compressed tar of reference genomes, select random genomes and place in an output directory.""" )
argp.add_argument("-n", "--bugCount", dest = "iBugCount", default = 100, type = int, action = "store", help = "The number of genomes to select.")
argp.add_argument("-s", "--sampleCount", dest = "iSampleCount", default = 10, type = int, action = "store", help = "The number of samples which will be generated." )
argp.add_argument("-l", "--lognormal", dest = "fLognormalDistribution", default = False, action = "store_true", help = "Indicates a lognormal distribution is desired." )
argp.add_argument("--t", "--taxonomy", dest = "sTaxonomyFile", default = "taxonomy_reduced.txt", action = "store", help = "The file from which taxonomy will be retrieved.")
argp.add_argument(dest = "sOutputDirectory", default = "Outputs", action = "store", help = "Output directory to place genomes.)")
argp.add_argument(dest = "sInputDirectories", nargs = "*", action = "store", help = "One or many input compressed directories.")
args = argp.parse_args()

# Holds the file and the archives that they come from
# {"genome1":"archive1", "genome2:archive2", "genome3:archive1" }
# So you later know where to uncompress selecte genomes from
dictGenomes = {}

# Holds the file names for the genomes that are selected
lsGenomeFilesSelected = []

# Holds the information for the taxon table to link the Files to the genome ids
# 15 tab delimited columns all which can be blank except
# Column 1: File name base
# Column 3: enum(Finished,Draft,Permenant Draft)
# Column 4: consensus name for genome genus+species+strain
# Column 13: genus
# Column 14: species
# Column 15: strain
lsTaxonTable = []

# Holds the file information for the community percentages
# bug name percent
# Percent will be made proportional
lsPercentageFile = []

# Holds the gold standard truth to the selected community abundances
lsGoldStandard = []

# Holds the taxonomy files
dictTaxonomy = {}

# The percentage for an even community
iCommunity = 100.0 / args.iBugCount

# Get all files and store indicating their archive
for sGenomeArchive in args.sInputDirectories:
  lsGenomes = [ sFile for sFile in glob.glob( os.path.join( sGenomeArchive, c_sFileExtension ) ) ]
  dictGenomes.update( dict( zip( lsGenomes, [ sGenomeArchive ] * len( lsGenomes ) ) ) )

# Select files
lsSelectedGenomes = random.sample( dictGenomes.keys(), args.iBugCount )

# Reduce dict to just the selected files
dictGenomes = dict( ( sSelectedGenomes, dictGenomes[ sSelectedGenomes ] ) for sSelectedGenomes in lsSelectedGenomes )

# Copy files
# Remove spaces from ids
for sGenomeFastaFile in lsSelectedGenomes:
  sOutFileLocation = os.path.join( args.sOutputDirectory, os.path.basename( sGenomeFastaFile ) )
  shutil.copyfile( sGenomeFastaFile, sOutFileLocation )
  subprocess.call( ["sed", "-i", "s; ;_;g", sOutFileLocation] )

# Read in taxoomy file
# Expected a file formated as follows
# bug_name\tinteger\ttaxonomy|pipe|delimited
with open( args.sTaxonomyFile, 'rU' ) as handleTaxonomy:
  csvr = csv.reader( handleTaxonomy, delimiter = "\t" )
  for lsTaxonomy in csvr:
    if len( dictTaxonomy.get( lsTaxonomy[ 0 ], "" ) ) < len( lsTaxonomy[ 2 ] ):
      dictTaxonomy[ lsTaxonomy[ 0 ] ] = lsTaxonomy[ 2 ]

# Make file indicating an even percent of sampling to occur per genome
# Normalized by genome length
# Assumes one genome in each file
for sGenomeFastaFile in lsSelectedGenomes:
  with open( os.path.join( args.sOutputDirectory,  os.path.basename( sGenomeFastaFile ) ), "rU" ) as handleFasta:

    # Genus, species, strain
    sGenusCur = "genus"
    sKingdomCur = "NA"
    sSpeciesCur = "species"
    sStrainCur = "strain"

    # Get bug name from the sequence information
    sDescription = SeqIO.parse( handleFasta, "fasta" ).next().description # Get description
    sDescription = sDescription.split( "|" )[-1] # Look at description and not sequence info
    sName = sDescription.split( "," )[0] # Remove genome completeness information
    sName = sName[1:] # Remove first character which will be underscore (after above sed)
    sGenomeState = sDescription.split( "," )[1] # Parse genome information starts with _complete_genome
    sGenomeState = sGenomeState.split( "_" )[1]
    sCleanedName  = string.replace( sName,"-","_" )  # The dashes are underscores in the taxonomy table so you have to search with a transformed name.
    
    # Potentially get genus, species, and strain
    if sCleanedName in dictTaxonomy:
      lsTaxonomy = dictTaxonomy[ sCleanedName ].split("|")
      
      for sTaxonomyItem in lsTaxonomy:
        if c_sGenusKey == sTaxonomyItem[0: len( c_sGenusKey )]:
          sGenusCur = sTaxonomyItem[len( c_sGenusKey ):]
        elif c_sSpeciesKey == sTaxonomyItem[0: len( c_sSpeciesKey )]:
          sSpeciesCur = sTaxonomyItem[len( c_sSpeciesKey ):]
        elif c_sStrainKey == sTaxonomyItem[0: len( c_sStrainKey )]:
          sStrainCur = sTaxonomyItem[len( c_sStrainKey ):]
        elif c_sKingdomKey == sTaxonomyItem[0: len( c_sKingdomKey )]:
          sKingdomCur =  sTaxonomyItem[len( c_sKingdomKey ):]
    else:
      print "Could not find: " + sCleanedName

    # Update taxon table
    lsTaxonTable.append( [ os.path.splitext( os.path.basename( sGenomeFastaFile ) )[ 0 ],sKingdomCur,sGenomeState,sName,"NA","NA","NA","NA","NA","NA","NA","NA",sGenusCur,sSpeciesCur,sStrainCur ] + ( [ "NA" ] * 132 ))
    lsPercentageFile.append( " ".join([ sGenusCur, sSpeciesCur ]) )
    lsGoldStandard.append( c_sGenusKey + sGenusCur )

# Input File: Taxon File
# Make Taxon Table file
with open( os.path.join( args.sOutputDirectory, c_sTaxonTableFileName ), "w" ) as handleTaxonTable:
  csv.writer( handleTaxonTable, delimiter = "\t" ).writerows( lsTaxonTable )

# Input file: Community Percents
# Make community proportion file
for iSample in xrange( 1, args.iSampleCount + 1 ):
  with open( os.path.join( args.sOutputDirectory, c_sCommunityPercentFileName + "_" + str( iSample ) + ".txt" ), "w" ) as handleCommunityPercent:
    # Input file: Gold standard for later evaluation
    # Make gold standard file
    with open( os.path.join( args.sOutputDirectory, c_sGoldStandardFileName + "_" + str( iSample ) + ".txt" ), "w" ) as handleGoldStandard:
      for iBugIndex, sBug in enumerate( lsPercentageFile ):
        # If needing a lognormal distribution randomly draw a community proportion
        if args.fLognormalDistribution:
          iCommunity = random.lognormvariate( mu = c_sLognormalMu, sigma = c_sLognormalSigma )
        csv.writer( handleCommunityPercent, delimiter = "\t" ).writerow( [ sBug, iCommunity ] )
        csv.writer( handleGoldStandard, delimiter = "\t" ). writerow( [ lsGoldStandard[iBugIndex], iCommunity ] )
