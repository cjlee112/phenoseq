import os, subprocess, sys

def sam2pileup(samfile, reffile):
    path, filename = os.path.split(samfile)
    filestem = filename.split('.')[0]
    subprocess.check_call(['samtools', 'view', '-bS', '-o', 
                           filestem + '_tmp.bam', samfile])
    subprocess.check_call(['samtools', 'sort', filestem + '_tmp.bam', 
                           filestem])
    pileupFile = open(filestem + '.pileup', 'w')
    subprocess.check_call(['samtools', 'pileup', '-f', reffile, 
                           filestem + '.bam'], stdout=pileupFile)
    pileupFile.close()
    os.remove(filestem + '_tmp.bam')

if __name__ == '__main__':
    reffile = sys.argv[1]
    for samfile in sys.argv[2:]:
        print 'processing', samfile
        sam2pileup(samfile, reffile)
