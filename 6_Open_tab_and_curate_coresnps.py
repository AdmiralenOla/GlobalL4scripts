# Script for opening a core.tab file from snippy-core, as well as the mummy fasta alignment and the core SNP alignment, and the restricted sites to be excluded,
# then using this info to evaluate whether to keep or remove a SNP

import Bio
import sys
import csv
from Bio import AlignIO

# USAGE: python <name.py> core.snps.aln core.tab mummies.fasta problematic-regions.bed

def CreateNewBioSeqObject(al, mummies):
	templist = []
	for record in al:
		newSeq = Bio.Seq.Seq("",Bio.Alphabet.SingleLetterAlphabet())
		newRec = Bio.SeqRecord.SeqRecord(seq=newSeq,id=record.id,name=record.name,description=record.description,dbxrefs=record.dbxrefs)
		templist.append(newRec)
	for record in mummies:
		newseq = Bio.Seq.Seq("",Bio.Alphabet.SingleLetterAlphabet())
		newRec = Bio.SeqRecord.SeqRecord(seq=newSeq,id=record.id,name=record.name,description=record.description,dbxrefs=record.dbxrefs)
		templist.append(newRec)
	return Bio.Align.MultipleSeqAlignment(records=templist,alphabet=Bio.Alphabet.SingleLetterAlphabet())

def AddBaseToSNPS(snps,base):
	counter = 0
	for rec in snps:
		rec.seq += base[counter]
		counter += 1
		
def IncludeCoord(coord, restrictedsites):
	includefromrow = 0
	for c in xrange(len(restrictedsites)):
		#if coord > restrictedsites[c][1]:
		#	includefromrow = c
			#continue
		if coord >= restrictedsites[c][0] and coord <= restrictedsites[c][1]:
			# Dont include
			includefromrow = c
			restrictedsites = restrictedsites[includefromrow:]
			return False, restrictedsites
	return True, restrictedsites
		
def FindBaseInMummies(coord, mums):
	mummiebase = []
	for record in mums:
		base = record.seq[coord-1]
		#print("Mummy %s base: %s %s" % (record.id, str(coord), str(base)))
		mummiebase.append(base)
	return mummiebase

snpcoords = []
newtab = {}
with open(sys.argv[2], "rU") as tabfile:
	tablines = tabfile.readlines()
    # Prepare new tab file
	#newtab.append(tablines[0])
	newtab[0] = tablines[0]
	tablines = tablines[1:]
	for line in tablines: 
		coord = line.split("\t")[1]
		snpcoords.append(int(coord))
		newtab[int(coord)] = line

restrictedsites = []
with open(sys.argv[4], "rU") as restrictfile:
	restrictlines = restrictfile.readlines()
	for line in restrictlines:
		start, stop = int(line.split("\t")[1]), int(line.split("\t")[2])
		restrictedsites.append((start, stop))

with open(sys.argv[1],"rU") as alignmentfile:
	al = Bio.AlignIO.read(alignmentfile,"fasta")
	print("Finished reading alignment file")
	totallength = al.get_alignment_length()
	
	with open(sys.argv[3], "rU") as mummyfile, open("core_script6.tab","w") as outtab:
		mums = Bio.AlignIO.read(mummyfile, "fasta")
		snps = CreateNewBioSeqObject(al, mums)
		# Patch in names of mummies in tab header:
		entrysplit = newtab[0].split("\t")
		entry = entrysplit[:-4] + [x.id for x in mums] + entrysplit[-4:]
		outtab.write("\t".join(entry))
		del newtab[0]
		#finaltab[0] += "\t".join([x.id for x in mums])

		for coord in xrange(1,totallength+1):
			# PROGRESS
			sys.stdout.write("\r %s " % (float(coord)/totallength))
			# 1. First see if base is filtered because of problematic region
			# 2. Evaluate if should be included after filtering of N/-
			# 3. If included, find base in mummies
			# 4. Add base to new snps
			
			# 1
			Include, restrictedsites = IncludeCoord(snpcoords[coord-1], restrictedsites)
			if Include:
				
				# 2
				# Is variable and not too many missing
				base = [record.seq[coord-1] for record in al]
				baseset = set(base)
				baseset.discard("N")
				num_N = base.count("N")
				prop_N = float(num_N)/len(base)
				baseset.discard("-")
				num_dash = base.count("-")
				prop_dash = float(num_dash)/len(base)
				prop_uninformative = prop_N + prop_dash
			
				if len(baseset) >= 2 and prop_uninformative < 0.01:
					# 3 Find base in mummies
					MummieBase = FindBaseInMummies(snpcoords[coord-1], mums)
					base += MummieBase
					AddBaseToSNPS(snps,base)
					# Add line to tab file
					entrysplit = newtab[snpcoords[coord-1]].split("\t")
					# Patch in mummy bases
					entry = entrysplit[:-4] + MummieBase + entrysplit[-4:]
					outtab.write("\t".join(entry))
					del newtab[snpcoords[coord-1]]
					#finaltab[snpcoords[coord-1]] = entry
					#tabtemp = ["NC_000962", str(snpcoords[coord-1])] + base
					#newtab.append("\t".join(tabtemp))
		with open("core_snp_alignment_script6.fasta","w") as outfile:
			Bio.AlignIO.write(snps, outfile, "fasta")



