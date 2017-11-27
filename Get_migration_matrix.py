#!/usr/bin/env python

'''Script for opening a BEAST tree simplified into newick format, along with a csv file with info about all the nodes and leaves
For each year in the total tree height, find out the number of branches in that of different types:
1. Read newick file
2. Read csv info file
3. Create a tree matrix like this:

Africa_to_Africa
Africa_to_EastAsia
Africa_to_Eurasia
Africa_to_NorthAmr
Africa_to_SouthAmr
EastAsia_to_Africa
...
SouthAmr_to_SouthAmr        0   
year                1000    1001    1002

4. Iterate over tree
5. For each node in tree
	- Find location
	- Find location of parent node
	- Find branch length
	- Add branch to tree matrix in correct years

Usage:
    python Get_migration_matrix.py treefile.nwk infofile.csv outfile.csv

'''
from ete3 import Tree
import sys
import csv
import math

def FindPlacement(height, length):
    end = int( math.ceil(rootheight) - math.ceil(height))
    start = int( end - math.ceil(length) )
    return start, end

def FindRow(parentloc, thisloc):
    if parentloc == "Africa":
        base = 0
    elif parentloc == "EastAsia":
        base = 1
    elif parentloc == "Eurasia":
        base = 2
    elif parentloc == "NorthAmr":
        base = 3
    elif parentloc == "SouthAmr":
        base = 4
    else:
        sys.exit("Invalid parent-daughter: %s, %s" % (parentloc, thisloc))

    if thisloc == "Africa":
        add = 0
    elif thisloc == "EastAsia":
        add = 1
    elif thisloc == "Eurasia":
        add = 2
    elif thisloc == "NorthAmr":
        add = 3
    elif thisloc == "SouthAmr":
        add = 4
    else:
        sys.exit("Invalid parent-daughter: %s, %s" % (parentloc, thisloc))

    return 5*base + add

def FindLocationOfParent(node):
    parent = node.up
    try:
        return infodic[parent.name]["location"]
    except AttributeError:
        sys.exit("Error finding parent of %s" % node.name)

tree = Tree(sys.argv[1], format=1)
rootheight = tree.get_farthest_node()[1]
print("Tree height: %s" % rootheight)

infodic = {}

with open(sys.argv[2],"rU") as infofile:
    csvreader = csv.reader(infofile,delimiter=",")
    header = csvreader.next()
    for row in csvreader:
        infodic[row[0]] = {"height":row[1], "length":row[2],"location":row[3]}

calendar = [[0 for year in xrange(int( math.ceil(rootheight)) + 1 ) ] for direction in xrange(25)]

for node in tree.iter_descendants("preorder"):
    try:
        nodeheight = float(infodic[node.name]["height"])
    except KeyError:
        sys.exit("Could not find in infodic: %s" % node.name)
    nodelength = float(infodic[node.name]["length"])
    nodelocation = infodic[node.name]["location"]

    nodestart, nodeend = FindPlacement(nodeheight,nodelength)
    nodeparentloc = FindLocationOfParent(node)
    row = FindRow(nodeparentloc, nodelocation)
    for z in xrange(nodestart,nodeend):
        calendar[row][z] += 1

with open(sys.argv[3],"w") as outfile:
    csvwriter = csv.writer(outfile,delimiter=",")
    continents = ["Africa","EastAsia","Eurasia","NorthAmr","SouthAmr"]
    for i in range(5):
        for j in range(5):
            from_to = continents[i] + "_to_" + continents[j]
            csvwriter.writerow([from_to] + calendar[(i*5+j)])


