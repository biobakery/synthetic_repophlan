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
import os.path
import shutil
import string
import subprocess

c_sTaxonTableFileName = "taxontable"
c_sCommunityPercentFileName = "Input"
c_sFileExtension  = "*.fna"
c_sGenusKey = "g__"
c_sGoldStandardFileName = "GoldStandard"
c_sKingdomKey = "k__"
c_sSpeciesKey = "s__"
c_sStrainKey = "t__"

# Header for the taxon table
lsTaxonTableHeader = ["taxon_oid","Domain","Status","Genome Name","Taxon Object ID","NCBI Taxon ID","RefSeq Project ID","GenBank Project ID","Phylum","Class","Order","Family","Genus","Species","Strain","Sequencing Center","Funding Agency","Add Date","Is Public","Release Date","IMG Release","IMG Product Assignment","IMG Submission ID","Proposal GOLD ID","Proposal Name","Biotic Relationships","Body Site","Body Subsite","Cell Arrangement","Cell Shape","Diseases","Energy Source","Ecosystem","Ecosystem Category","Ecosystem Type","Ecosystem Subtype","Isolation","Specific Ecosystem","Gram Staining","Host Name","Host Gender","Motility","Metabolism","Oxygen Requirement","Phenotype","Relevance","Salinity","Sporulation","Temperature Range","IMG Sample ID","IMG Project ID","Contact Name","Contact Email","Funding Program","Scaffold Count","CRISPR Count","GC Count","GC %","Coding Base Count","Genome Size","Gene Count","CDS Count","CDS %","RNA Count","rRNA Count","5S rRNA Count","16S rRNA Count","18S rRNA Count","23S rRNA Count","28S rRNA Count","tRNA Count","Other RNA Count","Fused Count","Fused %","Fusion Component Count","Fusion component %","Pseudo Genes Count","Pseudo Genes %","Unchar Count","Unchar %","Obsolete Count","Obsolete %","Revised Count","Revised %","w/ Func Pred Count","w/ Func Pred %","w/o Func Pred Sim Count","w/o Func Pred Sim %","w/o Func Pred No Sim Count","w/o Func Pred No Sim %","Signal Peptide Count","Signal Peptide %","Transmembrane Count","Transmembrane %","Orthologs Count","Orthologs %`","Paralogs Count","Paralogs %","Ortholog Group Count","Paralog Group Count","SwissProt Count","SwissProt %","Not SwissProt Count","Not SwissProt %","SEED Count","SEED %","Not SEED Count","Not SEED %","COG Count","COG %","KOG Count","KOG %","Pfam Count","Pfam %","TIGRfam Count","TIGRfam %","COG Cluster Count","KOG Cluster Count","Pfam Cluster Count","TIGRfam Cluster Count","InterPro Count","InterPro %","Enzyme Count","Enzyme %","TC Count","TC %","KEGG Count","KEGG %","Not KEGG Count","Not KEGG %","KO Count","KO %","Not KO Count","Not KO %","MetaCyc Count","MetaCyc %","Not MetaCyc Count","Not MetaCyc %","IMG Term Count","IMG Term %","IMG Pathwawy Count","IMG Pathway %","IMG Parts List Count","IMG Parts List %","Chromosomal Cassette Gene Count","Chromosomal Cassette Gene %","Chromosomal Cassette Count"]

# Note that the taxon table information being constructed here hold the following information
# The taxontable is needed to link reference genome files to the genome ids
# 15 tab delimited columns all which can be blank except
# Column 1: File name base
# Column 3: enum(Finished,Draft,Permenant Draft)
# Column 4: consensus name for genome genus+species+strain
# Column 13: genus
# Column 14: species
# Column 15: strain

# Converts prefixes in the header to the correct genome completion state
dictGenomeCompletion = {"AC":"Finished", "NC":"Finished", "NG":"Draft", "NT":"Draft", "NW":"Draft", "NS":"Draft", "NZ":"Draft" }

def funcReadTaxonomyFiles( lsTaxonomyFiles, dictTaxonomy = {}):
  """ Read in the taxonomy file """
  # Expected a file formated as follows
  # bug_name\tinteger\ttaxonomy|pipe|delimited

  for sTaxonomyFile in lsTaxonomyFiles:
    with open( sTaxonomyFile, 'rU' ) as handleTaxonomy:
      print "Get Random Genomes:: Reading taxonomy file "+sTaxonomyFile
      csvr = csv.reader( handleTaxonomy, delimiter = "\t" )
      dictTaxonomy.update( dict([[sFile, sTaxonomy] for sFile, sTaxonomy in csvr ]))
  return dictTaxonomy

def funcReadInArchives( lsGenomeArchives ):
  dictGenomes = {}
  for sGenomeArchive in lsGenomeArchives:
    print "Get Random Genomes:: Reading in archive " + sGenomeArchive
    lsGenomes = [ sFile for sFile in glob.glob( os.path.join( sGenomeArchive, c_sFileExtension ) ) ]
    print "Get Random Genomes:: The archive had a total of " + str( len( lsGenomes ) ) + " files."
    lsGenomes = [ sFile for sFile in lsGenomes if os.path.splitext( os.path.basename( sFile ) )[0] in dictTaxonomy]
    print "Get Random Genomes:: Of the files in the archive, " + str( len( lsGenomes ) )  + " were in the taxonomy file ( these were kept )."
    dictGenomes[ sGenomeArchive ] = lsGenomes
  return dictGenomes    

def funcSelectGenomes( iBugCount, iForced, dictGenomes, dictTaxonomy ):
  """ Return a set of selected files, forcing a certain number of files per archive if indicated."""
  # The genomes that are selected.
  lsSelectedGenomes = []
  # Holds all genomes not yet selected
  lsRemainingGenomes = []

  if iForced > 0:
    for sGenomeArchive, lsGenomes in dictGenomes.iteritems():
      print sGenomeArchive + " has " + str( len( lsGenomes ) ) + " genomes."   
      # If there is a forced amount of genomes to include, include.
      lsSelectedGenomesNew = random.sample( lsGenomes, iForced )
      lsSelectedGenomes.extend( lsSelectedGenomesNew )
      lsRemainingGenomes.extend( set(lsGenomes) - set( lsSelectedGenomesNew ) )

  # Select the rest of the files from the remaining genomes.
  print "Get Random Genomes:: Selecting "+str(iBugCount)+" random genomes from remaining " + str( len( lsRemainingGenomes ) ) + " genomes."
  lsSelectedGenomes.extend( random.sample( lsRemainingGenomes, max( iBugCount - len( lsSelectedGenomes ), 0 ) ))
  return lsSelectedGenomes

def funcPrepareSelectedGenomes( lsSelectedGenomes, sMoveToDirectory ):
  """ Copy files and remove spaces from ids """

  print "Get Random Genomes:: Moving selected genomes."
  # If the directory does not exist, make it.
  if not os.path.isdir( sMoveToDirectory ):
    os.mkdir( sMoveToDirectory )
  # Move the reference genomes if they do not exist.
  for sGenomeFastaFile in lsSelectedGenomes:
    sOutFileLocation = os.path.join( sMoveToDirectory, os.path.basename( sGenomeFastaFile ) )
    if not os.path.exists( sOutFileLocation ):
      shutil.copyfile( sGenomeFastaFile, sOutFileLocation )
      subprocess.call( ["sed", "-i", "s; ;_;g", sOutFileLocation] )

def funcSampleCommunity( lsSelectedGenomes, sReferenceFileArchive, sKingdomKey=c_sKingdomKey, sGenusKey=c_sGenusKey, sSpeciesKey=c_sSpeciesKey, sStrainKey=c_sStrainKey ):
  """ Make file indicating an even percent of sampling to occur per genome
      Assumes one genome in each file 
      Makes the community percentages, the gold standard, and taxontable information and returns it in lists."""

  lsPercentageFile = []
  lsGoldStandard = []
  lsTaxonTable = []
  print "Get Random Genomes:: Compiling taxonomy information for selected genomes."
  for sGenomeFastaFile in lsSelectedGenomes:
    with open( os.path.join( sReferenceFileArchive,  os.path.basename( sGenomeFastaFile ) ), "rU" ) as handleFasta:

      # Genus, species, strain
      sGenusCur = "genus"
      sKingdomCur = "NA"
      sSpeciesCur = "species"
      sStrainCur = "strain"

      # Get sequence completion from the sequence information
      sGenomeState = SeqIO.parse( handleFasta, "fasta" ).next().id.split("|")
      if len( sGenomeState ) > 1:
        sGenomeState = sGenomeState[3]
      else:
        sGenomeState = sGenomeState[0]
      sGenomeState = sGenomeState.split("_")[0]
      if sGenomeState not in dictGenomeCompletion:
        print "Get Random Genomes:: Did not recognize NCBI sequence prefix ( May be a Genbank entry ). Prefix ="+sGenomeState
        sGenomeState = "Draft"
      else:
        sGenomeState = dictGenomeCompletion[ sGenomeState ]
      sCleanedName  = os.path.splitext( os.path.basename( sGenomeFastaFile ) )[0]  # Using file name without extension
    
      # Potentially get genus, species, and strain
      sTaxonomy = dictTaxonomy[ sCleanedName ]
      lsGoldStandard.append( sTaxonomy )
      for sTaxonomyItem in sTaxonomy.split("|"):
        if c_sGenusKey == sTaxonomyItem[0: len( c_sGenusKey )]:
          sGenusCur = sTaxonomyItem[len( c_sGenusKey ):]
        elif c_sSpeciesKey == sTaxonomyItem[0: len( c_sSpeciesKey )]:
          sSpeciesCur = sTaxonomyItem[len( c_sSpeciesKey ):]
        elif c_sStrainKey == sTaxonomyItem[0: len( c_sStrainKey )]:
          sStrainCur = sTaxonomyItem[len( c_sStrainKey ):]
        elif c_sKingdomKey == sTaxonomyItem[0: len( c_sKingdomKey )]:
          sKingdomCur =  sTaxonomyItem[len( c_sKingdomKey ):]

      # Update taxon table
      lsTaxonTable.append( [ os.path.splitext( os.path.basename( sGenomeFastaFile ) )[ 0 ],sKingdomCur,sGenomeState,"_".join( [ sGenusCur,sSpeciesCur ] ) ] + ( [ "NA" ] * 8 ) + [ sGenusCur,sSpeciesCur,sStrainCur ] + ( [ "NA" ] * 132 ) )
      lsPercentageFile.append( " ".join([ sGenusCur, sSpeciesCur ] ) )

  return ( lsPercentageFile, lsGoldStandard, lsTaxonTable )

def funcAppendTaxonFile( lsTaxonTable, sTaxonomyFile ):
  """ Write the taxon table to a file.
      If one already exists, append new contents only """

  print "Get Random Genomes:: Making TaxonTable"
  lsTaxonomyContents = []
  sTaxonomyFile = os.path.join( args.sOutputDirectory, c_sTaxonTableFileName )
  if os.path.exists( sTaxonomyFile ):
    with open( sTaxonomyFile, "rU" ) as handleTaxonTableReader:
      lsTaxonomyContents = [ lsLine[0] for lsLine in csv.reader( handleTaxonTableReader, delimiter = "\t" )]
  with open( sTaxonomyFile, "a" ) as handleTaxonTable:
    # If the taxonomy file does not exist, write the header
    if not lsTaxonomyContents:
      csv.writer( handleTaxonTable, delimiter = "\t" ).writerow( lsTaxonTableHeader )

    # Get new contents, not contents already there and append to taxonomy file
    lsTaxonomyContent = [ lsTaxonInformation for lsTaxonInformation in lsTaxonTable if lsTaxonInformation[0] not in lsTaxonomyContents ]
    csv.writer( handleTaxonTable, delimiter = "\t" ).writerows( lsTaxonomyContent )

def funcMakeInputFile( lsPercentages, lsGoldStandard, sOutputDirectory, sInputFilePrefix, sGoldStandardFilePrefix, iFileIndex, fIsLognormal=False ):
  """ Make the input files to indicate community composition and gold standard information. """
  print "Get Random Genomes:: Making gold standards and input files."
  # Get current number of input files
  iEvenCommunityPercent = 1.0 / len( lsPercentages )
  with open( sInputFilePrefix + str( iFileIndex ) + ".txt", "w" ) as handleCommunityPercent:
    # Input file: Gold standard for later evaluation
    # Make gold standard file
    with open( sGoldStandardFilePrefix + "_" + str( iFileIndex ) + ".txt", "w" ) as handleGoldStandard:
      lsWriteActual = []
      lsWriteGold = []
      for iBugIndex, sBug in enumerate( lsPercentages ):
        # If needing a lognormal distribution randomly draw a community proportion
        iCommunity = random.lognormvariate( mu=1, sigma=2.7 ) if fIsLognormal else iEvenCommunityPercent
        lsWriteActual.append( [ sBug, iCommunity ] )
        lsWriteGold.append( [ lsGoldStandard[ iBugIndex ], iCommunity ] )
      csv.writer( handleCommunityPercent, delimiter = "\t" ).writerows( lsWriteActual )
      csv.writer( handleGoldStandard, delimiter = "\t" ). writerows( lsWriteGold )

def funcMinimumLengthFilter( dictGenomes, iMinimumLength ):
  """ Removes genomes from the list that do not have a certain length. """
  print "Get Random Genomes:: Filter for minimum contig length."
  for sArchive, lsGenomes in dictGenomes.iteritems():
    print "Get Random Genomes:: " + sArchive + " started with " + str( len( lsGenomes ) ) + " genomes."
    iGenome = 0
    # Go through each genome and make sure that a contig is longer than the given minimum length
    # If no contig in the genome is longer or equal to the minimum length,
    # Remove from possible genomes to be selected
    # If an archive has no genome of the needed length, remove it as well.
    for sGenome in lsGenomes:
      iGenome = iGenome + 1
      fPassMinimumLength = False
      for rcRecord in SeqIO.parse(sGenome, "fasta"):
        if len( rcRecord ) >= iMinimumLength:
          fPassMinimumLength = True
          break
      if not fPassMinimumLength:
        lsGenomes.remove( sGenome )
        if len( lsGenomes ) < 1 :
          del dictGenomes[ sArchive ]
      if iGenome % 100 == 0:
        print iGenome
    if sArchive in dictGenomes:
      print "Get Random Genomes:: " + sArchive + " has " + str( len( lsGenomes ) ) + " genomes after filtering."
      print str( lsGenomes )
    else:
      print "Get Random Genomes:: " + sArchive + " was deleted."
  return dictGenomes

argp = argparse.ArgumentParser( prog = "get_random_genomes.py", description = """Given a directory(ies) of reference genomes, select random genomes and place in an output directory.""" )
argp.add_argument( "-n", "--bugCount", dest = "iBugCount", default = 100, type = int, action = "store", help = "The number of genomes to select.")
argp.add_argument( "-s", "--sampleCount", dest = "iSampleCount", default = 10, type = int, action = "store", help = "The number of samples which will be generated. These samples will have the same comunity members but different abundances." )
argp.add_argument( "-l", "--lognormal", dest = "fLognormalDistribution", default = False, action = "store_true", help = "Indicates a lognormal distribution is desired." )
argp.add_argument( "-f", "--forced", dest = "iForcedAmount", default = 1, type = int, action = "store", help = "The minimal number of files from each archive which is guarenteed to exist.")
argp.add_argument( "-c","--contigLength", dest = "iMinimumLength", default = 400, action = "store", help = "A file must have atleast one contig that is of this length or equal. Default = 400." )
argp.add_argument( dest = "sOutputDirectory", default = "Outputs", action = "store", help = "Output directory to place genomes.)")
argp.add_argument( dest = "sInputDirectories", nargs = "*", action = "store", help = "One or many input repophlan directories with associated taxonomy files The taxonomy files are generated when repophlan makes the archive. Should be tab delimited and 2 entries, file name and then taxonomy.")
args = argp.parse_args()
print args

# Directory for reference genomes
sReferenceGenomesInput = os.path.join( args.sOutputDirectory, "user_genome" )

# Split up input files into taxonomy and reference genome archives
lsTaxonomyFiles = [ sReferenceAndTaxonomyFile.split(":")[1] for sReferenceAndTaxonomyFile in args.sInputDirectories ]
lsReferenceFiles = [ sReferenceAndTaxonomyFile.split(":")[0] for sReferenceAndTaxonomyFile in args.sInputDirectories ]

# Make output directory if not existing
if not os.path.exists( args.sOutputDirectory ):
  os.mkdir( args.sOutputDirectory )

# Get the combined taxonomy information
dictTaxonomy = funcReadTaxonomyFiles( lsTaxonomyFiles )
# Read in the files from the archives
dictGenomes = funcReadInArchives( lsReferenceFiles )

# Exclude genomes without contigs of a certain size.
dictGenomes = funcMinimumLengthFilter( dictGenomes, args.iMinimumLength )

for iSample in xrange( 1, args.iSampleCount+1 ):
  # Select genomes
  lsSelectedGenomes = funcSelectGenomes( iBugCount=args.iBugCount, iForced=args.iForcedAmount, dictGenomes=dictGenomes, dictTaxonomy=dictTaxonomy )
  # Make genomes ready and move associated genomes
  funcPrepareSelectedGenomes( lsSelectedGenomes, sReferenceGenomesInput )
  # Sample from the selected genomes and create content for required input files
  lsPercentages, lsGoldStandard, lsTaxonTable = funcSampleCommunity( lsSelectedGenomes, sReferenceGenomesInput )
  # Write the taxonomy table
  sTaxonomyFile = os.path.join( args.sOutputDirectory, c_sTaxonTableFileName )
  funcAppendTaxonFile( lsTaxonTable, sTaxonomyFile )
  # Make other input files
  sGoldStandardFilePrefix = os.path.join( args.sOutputDirectory, c_sGoldStandardFileName )
  sInputFilePrefix = os.path.join( args.sOutputDirectory, c_sCommunityPercentFileName + "_")
  funcMakeInputFile( lsPercentages=lsPercentages, lsGoldStandard=lsGoldStandard, sOutputDirectory=args.sOutputDirectory, sInputFilePrefix=sInputFilePrefix, sGoldStandardFilePrefix=sGoldStandardFilePrefix, iFileIndex=iSample, fIsLognormal=args.fLognormalDistribution )
