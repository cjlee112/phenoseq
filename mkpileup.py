import os, subprocess, sys

def sam2bam(samfile):
    path, filename = os.path.split(samfile)
    filestem = filename.split('.')[0]
    subprocess.check_call(['samtools', 'view', '-bS', '-o', 
                           filestem + '_tmp.bam', samfile])
    subprocess.check_call(['samtools', 'sort', filestem + '_tmp.bam', 
                           filestem])
    os.remove(filestem + '_tmp.bam')
    return filestem + '.bam'

def bam2pileup(bamfile, reffile):
    filestem = bamfile.split('.')[0]
    pileupFile = open(filestem + '.pileup', 'w')
    subprocess.check_call(['samtools', 'pileup', '-f', reffile, 
                           filestem + '.bam'], stdout=pileupFile)
    pileupFile.close()

def merge_bamfiles(bamfiles, tagfilter):
    d = {}
    for filename in bamfiles:
        tag = tagfilter(filename)
        d.setdefault(tag, []).append(filename)
    for tag, bamlist in d.items():
        print 'merging:', bamlist
        subprocess.check_call(['samtools', 'merge', tag + '.bam'] 
                              + bamlist)

def find_snps(bamfiles, reffile):
    for filename in bamfiles:
        filestem = filename.split('.')[0]
        bcfFile = open(filestem + '.bcf', 'w')
        subprocess.check_call(['samtools', 'mpileup', '-gf', reffile, 
                               filename], stdout=bcfFile)
        bcfFile.close()
        snpFile = open(filestem + '.vcf', 'w')
        subprocess.check_call(['bcftools', 'view', '-vcG', 
                               filestem + '.bcf'], stdout=snpFile)
        snpFile.close()

if __name__ == '__main__':
    reffile = sys.argv[1]
    for samfile in sys.argv[2:]:
        print 'processing', samfile
        bamfile = sam2bam(samfile)
        bam2pileup(bamfile, reffile)
