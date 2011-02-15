# -*- coding: utf-8 -*-

import csv
import glob
import os
import sys

#
# Illumina Native format of input filenames: s_(\d+)_(d+)_(\d+)_qseq.txt, where the numbers are lane, index, and tile
# index of 1 contains the sequence reads, index 2 contains the tags
#

def parse_tag_file(tag_filename):
    """Parses tag file into dicts. Import keys are 'Index' (the tag sequence), 'SampleID', and 'Lane'"""
    filename = "SampleSheet.csv"
    handle = open(filename)
    lines = handle.readlines()
    reader = csv.reader(lines)
    header = reader.next()
    dicts = [dict(zip(header, row)) for row in reader]
    return dicts

def parse_qseq(lines):
    header = ["MachineID", "run#", "lane#", "tile#", "x-coord", "y-coord", "index", "read#", "sequence", "q-scores", "p/f flag"]
    #handle = open(filename, "rU")
    dicts = [dict(zip(header, line.split('\t'))) for line in lines]
    return dicts

def dist(s, t):
    """Compute the distance between two strings, assuming they are of the same length."""
    d = 0
    try:
        for i in range(len(s)):
            if s[i] != t[i]:
                d += 1
    except IndexError:
        return float('inf')
    return d

def process_tile(seq_data, tag_data, tags):
    pools = {}
    for tag in tags:
        pools[tag] = []
    pools['other'] = []
    for i in range(len(seq_data)):
        seq = seq_data[i]['sequence']
        tag = tag_data[i]['sequence'][:6] # It's length 7 for some reason that I do not know, clip to 6.
        found = False
        for t in tags:
            d = dist(tag, t)
            if d < 2:
                pools[t].append(seq_data[i])
                found = True
                break
        if not found:
            pools['other'].append(seq_data[i])
    return pools

def discover_filenames(data_directory):
    filenames = glob.glob(os.path.join(data_directory, "s_*_qseq.txt"))
    s = set()
    for filename in filenames:
        split = filename.split('/')
        name = split[-1]
        _, lane, index, tile, _ = name.split('_')
        s.add((lane, tile))
    s = list(s)
    s.sort()
    return s

def main(argv):
    tag_filename, data_directory = argv[1], argv[2]
    # Load Tags from tag info file.
    tag_info = parse_tag_file(tag_filename)
    # The tags are in triplicate in the file, so compress to a set.
    tags = list(set([x['Index'] for x in tag_info]))
    tags.sort()
    # Discover data files: lanes and tiles. Not actually filenames, rather lane and tile pairs.
    filenames = discover_filenames(data_directory)
    keys = ["MachineID", "run#", "lane#", "tile#", "x-coord", "y-coord", "index", "read#", "sequence", "q-scores", "p/f flag"]
    # Prep data for processing for each file pair
    for (lane, tile) in filenames:
        seq_handle = open(os.path.join(data_directory, "s_%s_%s_%s_qseq.txt" % (lane, '1', tile)))
        tag_handle = open(os.path.join(data_directory, "s_%s_%s_%s_qseq.txt" % (lane, '2', tile)))
        seq_data = parse_qseq(seq_handle.readlines())
        tag_data = parse_qseq(tag_handle.readlines())
        # Identify tags and sort.
        pools = process_tile(seq_data, tag_data, tags)
        # Append to tag files
        for tag, value in pools.items():
            f = open("s_%s_%s_qseq.txt" % (lane, tag), 'a')
            lines = []
            for d in value:
                lines.append("\t".join([d[k] for k in keys]))
            f.writelines(lines)
            f.close()
        # Report results
        print lane, tile, ': '
        for tag, value in pools.items():
            print tag, len(value)

if __name__ == "__main__":
    main(sys.argv)

