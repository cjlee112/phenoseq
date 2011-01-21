import csv
import glob
import os
import tarfile
import numpy
from Bio import SeqIO

##### 
# Format of input filenames: s_(\d+)_(d+)_(\d+)_qseq.txt, where the numbers are lane, inde
x, and tile
# index of 1 is the sequence data index 2 is the tag
#####

def parse_tag_file():
    """Parses tag file into dicts. Import keys are 'Index' (the tag sequence), 'SampleID',
 and 'Lane'"""
    filename = "SampleSheet.csv"
    handle = open(filename)
    lines = handle.readlines()
    reader = csv.reader(lines)
    header = reader.next()
    dicts = [dict(zip(header, row)) for row in reader]
    return dicts

def parse_qseq(lines):
    header = ["MachineID", "run#", "lane#", "tile#", "x-coord", "y-coord", "index", "read#
", "sequence", "q-scores", "p/f flag"]
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
        tag = tag_data[i]['sequence'][:6] # It's length 7 for some reason that I do not kn
ow, clip to 6.
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
    #m = Multiset(starts)
    #for key, value in m.items():
        #print key, value

class TarInterface(object):
    def __init__(self, filename):
        self.filename=filename
        self.tar = tarfile.open(filename, 'r:gz')
        
    def filenames(self):
        filenames = []
        members = self.tar.getmembers()
        for member in members:
            if member.isfile:
                filenames.append(member.name)
        return filenames

    def extractfile(self, filename):
        return self.tar.extractfile(filename)

def discover_filenames(tar):
    #filenames = glob.glob("s_*_qseq.txt")
    filenames = tar.filenames()
    s = set()
    for filename in filenames:
        _, lane, index, tile, _ = filename.split('_')
        s.add((lane, tile))
    s = list(s)
    s.sort()
    return s

def init():
    # Create library directories.
    pass

def main():
    # Load Tags from tag info file.
    tag_info = parse_tag_file()
    tags = list(set([x['Index'] for x in tag_info])) # the tags are in triplicate in the f
ile, so compress to a set.
    tags.sort()
    # Discover Files
    tar = TarInterface('liao_data.tar.gz')
    filenames = discover_filenames(tar)
    # Prep data for processing for each file pair
    for (lane, tile) in filenames:
        seq_handle = tar.extractfile("s_%s_%s_%s_qseq.txt" % (lane, '1', tile))
        tag_handle = tar.extractfile("s_%s_%s_%s_qseq.txt" % (lane, '2', tile))
        seq_data = parse_qseq(seq_handle.readlines())
        tag_data = parse_qseq(tag_handle.readlines())
        # Identify tags and sort
        pools = process_tile(seq_data, tag_data, tags)
        # Report results
        print lane, tile, ': '
        for key, value in pools.items():
            print key, len(value)
        # Append to tag files
        keys = ["MachineID", "run#", "lane#", "tile#", "x-coord", "y-coord", "index", "rea
d#", "sequence", "q-scores", "p/f flag"]
        for key, value in pools.items():
            f = open("s_%s_%s_qseq.txt" % (lane, key), 'a')
            lines = []
            for d in value:
                lines.append("\t".join([d[k] for k in keys]))
            f.writelines(lines)
            f.close()

def convert_to_fastq():
    filenames = glob.glob("s_*_qseq.txt")
    for filename in filenames:
        _, lane, tag, _ = filename.split('_')
        new_f = open("s_%s_%s.fastq" % (lane, tag), 'w')
        f = open(filename)
        qseq2fastq(new_f, [f])


#if __name__ == "__main__":
    #main()
