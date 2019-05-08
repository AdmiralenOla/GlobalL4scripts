#!/usr/bin/env python

'''Script for opening a BEAST tree simplified into newick format, along with a csv file with info about all the nodes and leaves

Requirements:
- Newick tree file - Internal nodes must be labeled with name/number (Use Insert_node_number.py - can be used on Nexus file)
- Infodic. Contains, column-wise:
(0) The unique node name
(1) The height of the node, meaning the rightmost part, including tips 
(2) The length of the node, meaning branch distance from this node up to the parent node, defined as 0 for the root
(3) The location of this node (or tip)
(4) The probability of the estimate in (3) - not necessary
(5) Tip name, if any - not necessary
- Note - Under no circumstances should a nodes height over youngest tip + distance to root be greater than the total tree height

The output is a calendar with the columns being dates (years) and the rows being deme switches (including self-switches)
For each year, the number corresponds to what one would see by taking a cross-section of the tree at that particular year
The numbers in the cells represent the number of branches of that transition that exists at that year
Note that since there are no ways to know if a transition between two different demes happens on the left-most end of a branch,
on the right-most end, or somewhere in the middle, a branch is counted as a transition for its entire length. Ie, lets say we have a branch
with location A and parent location B, that starts in year X and ends in year Y. This branch would add +1 to the count of transition
"A_to_B" for all years from X to Y.

Usage:
    python Get_migration_matrix.py treefile.nwk infofile.csv outfile.csv

'''
from ete3 import Tree
import sys
import csv
import math

def FindPlacement(height, length):
    ''' Find out where a branch starts and where it ends, calender-wise'''
    end = int( math.ceil(rootheight) - math.ceil(height))
    start = int( end - math.floor(length) )
    return start, end

def FindRow(parentloc, thisloc, unique_locs):
    ''' Find out which row of the calendar this branch belongs to, as defined by parent and current location '''
    try:
        parentloc_index = unique_locs.index(parentloc)
        this_index = unique_locs.index(thisloc)
    except ValueError:
        sys.exit("Invalid parent-daughter: %s, %s" % (parentloc, thisloc))

    return parentloc + "_to_" + thisloc

def FindLocationOfParent(node):
    ''' Given a node, find the location of its parent node '''
    parent = node.up
    try:
        return infodic[parent.name]["location"]
    except AttributeError:
        sys.exit("Error finding parent of %s" % node.name)

# Read the tree and calculate root height. Format = 1 means flexible with internal node names
tree = Tree(sys.argv[1], format=1)
rootheight = tree.get_farthest_node()[1]
print("Tree height: %s" % rootheight)



infodic = {}

with open(sys.argv[2],"rU") as infofile:
    csvreader = csv.reader(infofile,delimiter=",")
    header = next(csvreader)
    for row in csvreader:
        infodic[row[0]] = {"height":row[1], "length":row[2],"location":row[3], "countryprob":row[4],"isolate":row[5]}

# Get the unique locations
unique_locs = list( sorted( set( [infodic[sub]["location"] for sub in infodic] ) ) )
num_unique_locs = len(unique_locs)

# Initiate the calendar empty, meaning set all migrations to 0 at first.
# Consists of TREEHEIGHT+1 number of columns
# Consists of NUM_UNIQUE_LOCS ** 2 rows  - (all pairs including self-pairs)

# Assume years are contained in the prefix of the input names
isolates = list( sorted( set( [infodic[sub]["isolate"] for sub in infodic] ) ) )
# The last is always going to be NA, thus, the second last will have the most recent isolate year in the prefix
mostcurrentyear = int(isolates[-2][:4])
rootyear = mostcurrentyear - math.ceil(rootheight)

print("Most current year: %s" % mostcurrentyear)
print("Year of the root: %s" % rootyear)

row_keys = [fromloc + "_to_" + toloc for fromloc in unique_locs for toloc in unique_locs]
column_keys = list(range(rootyear,mostcurrentyear+1))

calendar = {r: {c:0 for c in column_keys} for r in row_keys}

# Iterate over all nodes in the tree
for node in tree.iter_descendants("preorder"):
    try:
        # Get height of current node
        nodeheight = float(infodic[node.name]["height"])
    except KeyError:
        sys.exit("Could not find in infodic: %s" % node.name)
    # Get length of current node
    nodelength = float(infodic[node.name]["length"])
    # Get location of current node
    nodelocation = infodic[node.name]["location"]

    # Get the starting time and the ending time of the branch extending to the current node
    nodestart, nodeend = FindPlacement(nodeheight,nodelength)
    # Convert to year-information
    nodestart += rootyear
    nodeend += rootyear

    # Get the location of the parent of the current node
    nodeparentloc = FindLocationOfParent(node)
    # Find out which row of the calendar this parent-daughter pair corresponds to
    row = FindRow(nodeparentloc, nodelocation, unique_locs)
    for z in range(nodestart,nodeend+1):
        try:
            calendar[row][z] += 1
        except KeyError:
            print(z)
            print(node)
            sys.exit(-1)


with open(sys.argv[3],"w") as outfile:
    csvwriter = csv.writer(outfile,delimiter=",")
    csvwriter.writerow([""] + column_keys)
    for transition in row_keys:
        writelist = [transition]
        #csvwriter.writerow(calendar[transition])
        for year in column_keys:
            writelist.append(calendar[transition][year])
        csvwriter.writerow(writelist)



