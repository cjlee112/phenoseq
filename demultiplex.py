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

def parse_qseq(file_handle):
    """Generator for access to sequencer reads."""
    header = ["MachineID", "run#", "lane#", "tile#", "x-coord", "y-coord", "index", "read#", "sequence", "q-scores", "p/f flag"]
    for line in file_handle:
        yield line.split('\t')

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

def process_tile(seq_data, tag_data, tags, max_count=100000):
    """Generator for buffered demultiplexing to handle large qseq files, limiting memory usage."""
    def init_pools_dict():
        pools = {}
        for tag in tags:
            pools[tag] = []
        pools['other'] = []
        return pools
    pools = init_pools_dict()
    count = 0
    for x in seq_data:
        seq = x[8]
        y = tag_data.next()
        tag = y[8][:6] # It's length 7 for some reason that I do not know, clip to 6.
        found = False
        for t in tags:
            d = dist(tag, t)
            if d < 2:
                pools[t].append(x)
                found = True
                break
        if not found:
            pools['other'].append(x)
        count += 1
        if count >= max_count:
            count = 0
            yield pools
            pools = init_pools_dict()
    yield pools

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

def write_output(output_dir, lane, tile, pools):
    """Helper function for buffered output to limit memory consumption for large qseq files."""
    for tag, value in pools.items():
        f = open(os.path.join(output_dir, "s_%s_%s_qseq.txt") % (lane, tag), 'a')
        lines = []
        for d in value:
            lines.append("\t".join(d))
        f.writelines(lines)
        f.close()


def main(argv):
    ### Tagfile "SampleSheet.csv" not available on second run so tags must be given manually
    #tag_filename, data_directory = argv[1], argv[2]
    ## Load Tags from tag info file.
    #tag_info = parse_tag_file(tag_filename)
    ## The tags are in triplicate in the file, so compress to a set.
    #tags = list(set([x['Index'] for x in tag_info]))

    data_directory = argv[1]
    #e.g. tags = ATCACG CGATGT TTAGGC TGACCA ACAGTG GCCAAT CAGATC ACTTGA GATCAG TAGCTT GGCTAC CTTGTA
    tags = argv[2:]
    output_dir = data_directory + '-demultiplexed'
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    # Discover data files: lanes and tiles. Not actually filenames, rather lane and tile pairs.
    filenames = discover_filenames(data_directory)
    keys = ["MachineID", "run#", "lane#", "tile#", "x-coord", "y-coord", "index", "read#", "sequence", "q-scores", "p/f flag"]
    # Prep data for processing for each file pair
    for (lane, tile) in filenames:
        seq_handle = open(os.path.join(data_directory, "s_%s_%s_%s_qseq.txt" % (lane, '1', tile)))
        tag_handle = open(os.path.join(data_directory, "s_%s_%s_%s_qseq.txt" % (lane, '2', tile)))
        seq_data = parse_qseq(seq_handle)
        tag_data = parse_qseq(tag_handle)
        # Identify tags and sort. Append to tag files. Pass optional max_counts argument to change number of reads processed per output, default = 100,000 reads per output.
        pools_gen = process_tile(seq_data, tag_data, tags)
        pool_counts = dict()
        for pools in pools_gen:
            for tag, value in pools.items():
                try:
                    pool_counts[tag] += len(value)
                except KeyError:
                    pool_counts[tag] = len(value)
            write_output(output_dir, lane, tile, pools)
        # Report tagging distribution.
        print lane, tile, ': '
        for tag, value in pool_counts.items():
            print tag, value

if __name__ == "__main__":
    main(sys.argv)
