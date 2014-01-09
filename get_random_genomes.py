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
import subprocess
import tarfile

c_sTaxonTableFileName = "TaxonTable.tsv"
c_sCommunityPercentFileName = "Input"
c_sGoldStandardFileName = "GoldStandard.tsv"

argp = argparse.ArgumentParser( prog = "get_random_genomes.py", description = """Given a supset compressed tar of reference genomes, select random genomes and place in an output directory.""" )
argp.add_argument("-n", "--bugCount", dest = "iBugCount", default = 100, type = int, action = "store", help = "The number of genomes to select.")
argp.add_argument("-s", "--sampleCount", dest = "iSampleCount", default = 10, type = int, action = "store", help = "The number of samples which will be generated." )
argp.add_argument(dest = "sOutputDirectory", default = "Outputs", action = "store", help = "Output directory to place genomes.)")
argp.add_argument(dest = "sInputDirectories", nargs = "*", action = "store", help = "One or many input compressed directories.")
args = argp.parse_args()
print args

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

# The percentage for an even community
iCommunity = 100.0 / args.iSampleCount

# For each tar file 
for sGenomeArchive in args.sInputDirectories:
  print "sGenomeArchive"
  print sGenomeArchive
  tarProcessor = tarfile.open( sGenomeArchive )
  lsGenomes = tarProcessor.getnames()
  dictGenomes.update( dict( zip( lsGenomes, [ sGenomeArchive ] * len( lsGenomes ) ) ) )
  print len( dictGenomes )
  tarProcessor.close()

# Select files
lsSelectedGenomes = random.sample( dictGenomes.keys(), args.iBugCount )
print lsSelectedGenomes
# Reduce dict to just the selected files
dictGenomes = dict( ( sSelectedGenomes, dictGenomes[ sSelectedGenomes ] ) for sSelectedGenomes in lsSelectedGenomes )

# Extract the selected files
for sArchive in set( dictGenomes.itervalues() ):
  print sArchive
  lsSelectedGenomes = [ key for key, value in dictGenomes.iteritems() if value == sArchive ]
  print len( lsSelectedGenomes )
  tarProcessor = tarfile.open( sArchive )
  for sSelectedGenome in lsSelectedGenomes:
    tarProcessor.extract( sSelectedGenome, args.sOutputDirectory )
  tarProcessor.close()

# Flatten the file structure
ltplTerminalFiles = [ sFile for sFile in os.walk( args.sOutputDirectory ) if sFile[-1] ]
print ltplTerminalFiles
for tplTerminalFile in ltplTerminalFiles:
  for sFile in tplTerminalFile[-1]:
    sTerminalFile = os.path.join( tplTerminalFile[ 0 ], sFile)
    lsGenomeFilesSelected.append( sFile )
    shutil.move( sTerminalFile, args.sOutputDirectory )
  shutil.rmtree( tplTerminalFile[ 0 ] )

# Remove spaces from ids
for sGenomeFastaFile in lsGenomeFilesSelected:
  subprocess.call( ["sed", "-i", 's/ /_/g', os.path.join( args.sOutputDirectory, sGenomeFastaFile )] )

# Make file indicating an even percent of sampling to occur per genome
# Normalized by genome length
# Assumes one genome in each file
print lsGenomeFilesSelected
for sGenomeFastaFile in lsGenomeFilesSelected:
  with open( os.path.join( args.sOutputDirectory, sGenomeFastaFile ), "rU" ) as handleFasta:
    iReferenceLength = 0
    sName = ""
    for SeqCur in SeqIO.parse( handleFasta, "fasta" ):
      iReferenceLength += len( SeqCur )
      sName = SeqCur.id.split("|")[-1].strip("_")
    sName = sName.split(",")[0]
    print "NAME"
    print sName
    print iReferenceLength
    # Update taxon table
    lsTaxonTable.append( [ os.path.splitext( sGenomeFastaFile )[ 0 ],"","Finished",sName.split(",")[0],"","","","","","","","","genus","species","strain" ] )
    print iReferenceLength
    print iCommunity
    print iCommunity / iReferenceLength
    lsPercentageFile.append( [ sName.split(",")[0], iCommunity / iReferenceLength ]  )
    lsGoldStandard.append( [ sName.split(",")[0], iCommunity ])

# Make Taxon Table file
with open( os.path.join( args.sOutputDirectory, c_sTaxonTableFileName ), "w" ) as handleTaxonTable:
  csv.writer( handleTaxonTable, delimiter = "\t" ).writerows( lsTaxonTable )

# Make community proportion file
for iSample in xrange( args.iSampleCount ):
  with open( os.path.join( args.sOutputDirectory, c_sCommunityPercentFileName + "_" + str( iSample ) ), "w" ) as handleCommunityPercent:
    csv.writer( handleCommunityPercent, delimiter = "\t" ).writerows( lsPercentageFile )

# Make gold standard file
with open( os.path.join( args.sOutputDirectory, c_sGoldStandardFileName ), "w" ) as handleGoldStandard:
  csv.writer( handleGoldStandard, delimiter = "\t" ). writerows( lsGoldStandard )
