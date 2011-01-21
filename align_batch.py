import glob
import subprocess
import sys

"""
usage: python align_batch.py <novoalign_executable> <index_filename>

e.g.   python align_batch.py novoalign/novocraft/novoalign ssuis.nix
"""

if __name__ == '__main__':
    novoalign_executable = sys.argv[1]
    index_filename = sys.argv[2]
    filenames = glob.glob("s_*?_*?_qseq.txt")
    for filename in filenames:
        _, lane, tag, _ = filename.split('_')
        # Skip the file of undecipherable reads.
        if tag == "other":
            continue
        cmd = "./%s -d %s -o SAM -f s_%s_%s_qseq.txt > aligned_s_%s_%s.sam" % (novoalign_executable, index_filename, lane, tag, lane, tag)
        subprocess.call(cmd)
