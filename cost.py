import glob
import subprocess
from math import sqrt
try:
    import numpy
    from matplotlib import pyplot
except ImportError:
    pass

def get_nstrain(tag):
    'return #strains in each pool'
    if 'ACAGTG' in tag or 'CTTGTA' in tag:
        return 4
    else:
        return 3

def sum_nstrain(tags):
    return sum([get_nstrain(tag) for tag in tags])

def sum_nhit(besthits, trueHits=('acrB', 'marC', 'acrA')):
    return sum([1 for t in besthits if t[1] in trueHits])

def nstrain_nhit_dict(subsetDict, nbest=20, results=None):
    if results is None:
        results = {}
    for tags,hits in subsetDict.items():
        nstrain = sum_nstrain(tags)
        nhit = sum_nhit(hits[:nbest])
        results.setdefault(nstrain, []).append(nhit)
    return results

def get_mean_stderr(x):
    'get maximum likelihood estimates of mu, sigma from sample x'
    x = numpy.array(x)
    mean = numpy.mean(x)
    diff = x - mean
    var = numpy.mean(diff * diff) / len(x) # get std error
    return mean, sqrt(var)

def nstrain_nhit_means(d):
    'get sorted list of [(nstrain, meanHits, stdError)]'
    l = []
    for nstrain, nhits in d.items():
        l.append((nstrain,) + get_mean_stderr(nhits))
    l.sort()
    return l

# cost for 3 lane expt is 3 * 1400 / 32. + 10 * 200 / 32. per strain

def run_subexperiments(annodb, al, dna):
    # 3 lane subexperiments, requiring each SNP replicated in 2 of 3 lanes
    l = glob.glob('[ACGT]*.vcf')
    d = analyze.generate_subsets(l, annodb, al, dna)

    # 1 lane subexperiments (probably should use a quality cutoff)
    vcf1 = glob.glob('aligned_s_1_*.vcf')
    d1 = analyze.generate_subsets(vcf1, annodb, al, dna,
                                  lambda s:s.split('.')[0].split('_')[-1],
                                  lambda s:(), 0)
    vcf2 = glob.glob('aligned_s_2_*.vcf')
    d2 = analyze.generate_subsets(vcf2, annodb, al, dna,
                                  lambda s:s.split('.')[0].split('_')[-1],
                                  lambda s:(), 0)
    vcf3 = glob.glob('aligned_s_3_*.vcf')
    d3 = analyze.generate_subsets(vcf3, annodb, al, dna,
                                  lambda s:s.split('.')[0].split('_')[-1],
                                  lambda s:(), 0)

    # 2 lane subexperiments, requiring each SNP replicated in both lanes
    vcf_1_2 = glob.glob('[ACGT]*_1_2.vcf')
    d_1_2 = analyze.generate_subsets(vcf_1_2, annodb, al, dna,
                                     lambda s:s.split('_')[0],
                                     lambda s:('aligned_s_1_' + s + '.vcf',
                                               'aligned_s_2_' + s + '.vcf'), 2)
    vcf_1_3 = glob.glob('[ACGT]*_1_3.vcf')
    d_1_3 = analyze.generate_subsets(vcf_1_3, annodb, al, dna,
                                     lambda s:s.split('_')[0],
                                     lambda s:('aligned_s_1_' + s + '.vcf',
                                               'aligned_s_3_' + s + '.vcf'), 2)
    vcf_2_3 = glob.glob('[ACGT]*_2_3.vcf')
    d_2_3 = analyze.generate_subsets(vcf_2_3, annodb, al, dna,
                                     lambda s:s.split('_')[0],
                                     lambda s:('aligned_s_2_' + s + '.vcf',
                                               'aligned_s_3_' + s + '.vcf'), 2)
    return (d1, d2, d3), (d_1_2, d_1_3, d_2_3), d

def plot_cost(data, costPerStrain=3 * 1400 / 32. + 10 * 200 / 32., **kwargs):
    pyplot.errorbar([(t[0] * costPerStrain) for t in data],
                    [t[1] for t in data],
                    [t[2] for t in data], **kwargs)

def merge_names(names):
    tags = [s.split('.')[0].split('_')[-1] for s in names]
    lanes = [s.split('_')[2] for s in names]
    l = [(tags[i] + '_' + lanes[i]) for i in range(len(names))]
    return '_'.join(l) + '.bam'

def generate_tag_pairs(tags):
    for i, tag1 in enumerate(tags):
        for j in range(i + 1, len(tags)):
            yield tag1, tags[j]

def get_tag_lane_file(tag, lane):
    return 'aligned_s_%d_%s.bam' % (lane, tag)
    
def merge_libs(tags, lanes=range(1,4), bamFunc=get_tag_lane_file, 
               outFunc=merge_names):
    for i in lanes:
        for j in lanes:
            for tag1,tag2 in generate_tag_pairs(tags):
                bam1 = bamFunc(tag1, i)
                bam2 = bamFunc(tag2, j)
                outFile = outFunc((bam1, bam2))
                print 'merging', outFile
                subprocess.check_call(('samtools', 'merge', outFile,
                                       bam1, bam2))
