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

def get_tag(filename):
    'extract tag string from filename'
    return filename.split('.')[0].split('_')[-1]

def get_merge_name(tag, bamlist):
    'create filename for merging list of bam files'
    l = [s.split('_')[2] for s in bamlist]
    l.sort()
    return '_'.join([tag] + l) + '.bam'

def generate_sublists(bamlist):
    'generate leave-one-out jackknife samples'
    for i in range(len(bamlist)):
        yield bamlist[:i] + bamlist[i + 1:]

def merge_bamfiles(bamfiles, tagfilter=get_tag, nameFunc=None, 
                   sublistFunc=lambda x:(x)):
    d = {}
    for filename in bamfiles:
        tag = tagfilter(filename)
        d.setdefault(tag, []).append(filename)
    for tag, bamlist in d.items():
        for sublist in sublistFunc(bamlist):
            print 'merging:', sublist
            if nameFunc is not None: # get customized name from user function
                outFile = nameFunc(tag, sublist)
            else:
                outFile = tag + '.bam'
            subprocess.check_call(['samtools', 'merge', outFile] 
                                  + sublist)

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
        find_snps([bamfile], reffile)
        #bam2pileup(bamfile, reffile)
