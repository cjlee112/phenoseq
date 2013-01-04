from math import log, exp, sqrt
from optparse import OptionParser
from pathways import load_func_assoc
from pygr import sequtil
import numpy
from matplotlib import pyplot
from scipy.stats import binom,fisher_exact
import analyze
from pathways import parse_args

rc = dict(A='T', C='G', G='C', T='A')

## Helpers

def normalize_dict(d_):
    s = sum(map(float, d_.values()))
    d = dict()
    for k,v in d_.items():
        d[k] = float(v) / s
    return d

## Mutation spectra.

def spontaneous_mutations():
    d = {('AT','GC'): 0.633, ('GC', 'AT'): 0.326, ('AT','TA'): 0., ('GC', 'TA'): 0.0247, ('AT', 'CG'): 0.011, ('GC', 'CG'): 0.0055}
    return d

def empirical_mutations():
    d = {('AT','GC'): 89./4099, ('GC', 'AT'): 3960./4099, ('AT','TA'): 3./4099, ('GC', 'TA'): 3./4099, ('AT', 'CG'): 19./4099, ('GC', 'CG'): 25./4099}
    return d

mutation_types = [('AT','GC'), ('GC', 'AT'), ('AT','TA'), ('GC', 'TA'), ('AT', 'CG'), ('GC', 'CG')]



## Ka/Ks ratio ##

def codon_weighting(mutation_dist, codon_dist=None, mu=1.):
    # Generate all possible codons.
    codons = []
    for a in 'ATCG':
        for b in 'ATCG':
            for c in 'ATCG':
                codons.append(a + b + c)
    amino_acids = ['*', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # If no codon distribution is given, assume it is uniform.
    if not codon_dist:
        codon_dist = {}
        l = 1./ float(len(codons))
        for codon in codons:
            codon_dist[codon] = l
    at_mutations = {}
    gc_mutations = {}
    for k, v in mutation_dist.items():
        if k[0] == 'AT':
            at_mutations[k[1]] = v
        else:
            gc_mutations[k[1]] = v
    normed_mutations = {}
    normed_mutations['AT'] = normalize_dict(at_mutations)
    normed_mutations['GC'] = normalize_dict(gc_mutations)

    # For each codon, mutate each of the three bases according to the mutation_dist, and record the resultant codon and its probability. For each codon the result is a probability distribution over the possible amino acids.
    codon_transitions = {}
    for codon in codons:
        for i in [0,1,2]:
            for k,v in mutation_dist.items():
                for j in [0,1]:
                    if codon[i] == k[0][j]:
                        result_codon = list(codon)
                        result_codon[i] = k[1][j]
                        result_codon = "".join(result_codon)
                        # 1/3 assumes mutation is equally likely at each position in the codon.
                        codon_transitions[(codon, result_codon)] = mu * 1./3 * normed_mutations[k[0]][k[1]]
        codon_transitions[(codon, codon)] = 1.-mu
    return codon_transitions



def sites_count_dict(codon_transitions=codon_weighting(empirical_mutations())):
    """Count sites for each codon, weighted by a mutation distribution."""
    codons = []
    for a in 'ATCG':
        for b in 'ATCG':
            for c in 'ATCG':
                codons.append(a + b + c)
    d = dict()
    for codon in codons:
        a = 0
        s = 0
        for i in range(len(codon)):
            bases = list('ATCG')
            bases.remove(codon[i])
            aa = sequtil.translate_orf(codon)
            for base in bases:
                new_codon = list(codon)
                new_codon[i] = base
                new_codon = "".join(new_codon)
                new_aa = sequtil.translate_orf(new_codon)
                if new_aa == aa:
                    s += 3 * codon_transitions[(codon, new_codon)]
                else:
                    a += 3 * codon_transitions[(codon, new_codon)]
        d[codon] = (a,s)
    return d



def count_codons(seq=None):
    """ Counts codons in dna sequence."""
    codons = {}
    for i in range(0,len(seq), 3):
        codon = seq[i:i+3]
        try:
            codons[codon] += 1
        except KeyError:
            codons[codon] = 1
    return codons



def count_sites(codons, d=sites_count_dict()):
    """Counts total sites in a dict (multiset) of codons."""
    Na = 0
    Ns = 0
    for codon, count in codons.items():
        try:
            (a,s) = d[codon]
            Na += a * count
            Ns += s * count
        except KeyError:
            pass
    return (Na, Ns)



def genes_sites_dict(annodb):
    """Compute site counts for all genes and return a dictionary with the counts."""
    site_counts = {}
    for k in annodb.keys():
        seq = str(annodb[k].sequence)
        codons = count_codons(seq)
        (Na, Ns) = count_sites(codons)
        site_counts[k] = (len(seq), Na, Ns)
    return site_counts



def ka_ks_counts_for_gene_group(genes, snp_counts, site_counts, pseudo=1.0):
    """Computes the contigency table for a group of genes."""
    Na, Ns = (0, 0)
    Ga, Gs = (0, 0)
    for gene in genes:
        try:
            l, a, s = site_counts[gene]
            Na += a
            Ns += s
        except KeyError:
            continue
        try:
            a, s = snp_counts[gene]
            Ga += a
            Gs += s
        except KeyError:
            pass

    a,b,c,d = map(lambda x: int(x + pseudo), [Ga, Gs, Na, Ns])
    return a, b, c, d
    


## Fisher Exact Test    
def confidence_interval(a, b, c, d, z=1.65):
    """Confidence interval fo Fisher Exact test."""
    try:
        s = sqrt(1./a + 1./b + 1./c + 1./d)
    except ZeroDivisionError:
        return (0,0)
    try:
        l = log(a * d / (b * c))
    except ValueError:
        return (0,0)

    return (exp(l - z * s), exp(1 + z*s))



def fisher_test(a, b, c, d):
    oddsratio, pvalue = fisher_exact([[a, b], [c, d]],
                                     alternative='greater')
    lower, upper = confidence_interval(a, b, c, d)
    return pvalue, lower, upper, oddsratio
    # rpy2 version
    #import rpy2.robjects as robjects
    #robjects.r.options(warn=-1) # supress warnings
    #fisher_test = robjects.r['fisher.test']
    #matrix = robjects.r['matrix']
    #v = robjects.IntVector([Ga, Gs, int(Ga_exp+0.5), int(Gs_exp+0.5)]) 
    #fet = fisher_test(matrix(v, nrow=2), alternative="greater")
    #print Ga_exp, Gs_exp, fet[0][0], fet[1][0], fet[1][1] # p-value, confidence interval (95%)



## Binomial Test

def binomial_test(Ga, Gs, Na, Ns, z=1.65):
    """Returns binomial p-value and Wilson confidence interval."""
    if Gs == 0:
        return None
    theta = float(Na) / (Na + Ns)
    p = float(Ga) / (Ga + Gs)
    n = int(Ga + Gs)
    # p-value
    rv = binom(n, theta)
    pvalue = 1. - rv.cdf(int(Ga))
    # confidence interval (Wilson) http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    b = z * sqrt(p * (1.-p) / n + z*z / (4. * n*n))
    lower = (p + 1. /(2*n)*z*z - b) / (1. + 1./n * z * z)
    upper = (p + 1. /(2*n)*z*z + b) / (1. + 1./n * z * z)
    return (pvalue, lower, upper, p, theta)
    


def perform_tests(gene_groups, snp_counts, site_counts, test_func=fisher_test, pseudo=1.0):
    """Performs Fisher's Exact test for all genes individually."""
    tests = []
    for gene_group in gene_groups:
        a, b, c, d = ka_ks_counts_for_gene_group(gene_group, snp_counts, site_counts, pseudo=1.0)
        tests.append(test_func(a, b, c, d))
    return tests
    


def run_tests(site_counts, snp_counts, genes, functional_groups,
              test_func=fisher_test):
    ## Genome-wide test
    genome_test = perform_tests([genes], snp_counts, site_counts, test_func=test_func)
    print "Genome-wide", genome_test[0]
    
    ## Individual tests on gene groups [gene]
    print "Individual gene tests"
    individual_tests = perform_tests([[gene] for gene in genes], snp_counts, site_counts, test_func=test_func)
    test_dict = dict(zip(genes, individual_tests))
    results = [(individual_tests[i][0], genes[i]) for i in range(len(genes))]
    results.sort()
    for p, gene in results:
        print "%s,%s,%s" % (p, gene, ",".join(map(str, test_dict[gene])))

    ## Functionally associate groups test
    # Don't bother with empty groups.
    print "Association group tests"
    items = [(k, v) for (k,v) in functional_groups.items() if len(v[1]) > 0]
    group_tests = perform_tests([v[1] for (k, v) in items], snp_counts, site_counts, test_func=test_func)
    test_dict = dict(zip([k for (k,v) in items], group_tests))
    results = [(group_tests[i][0], items[i][0]) for i in range(len(items))]
    results.sort()
    for p, group_name in results:
        print "%s,%s,%s,%s" % (p, group_name, " ".join(map(str,functional_groups[group_name][1])), ",".join(map(str,test_dict[group_name])))
    

def run_all(gbfile, groupfile, transfile, tagFiles):
    # does this filter out replicates that appear in every tag?
    annodb, al, dna = analyze.read_genbank_annots(gbfile)
    snps = analyze.read_tag_files(tagFiles)
    gsd = analyze.map_snps_chrom1(snps, al, dna)

    # Count nonsyn vs. syn. sites for each gene
    site_counts = genes_sites_dict(annodb)
    # Count nonsyn vs. syn. snps for each gene
    snp_counts = analyze.get_gene_na_ns(gsd)
    
    genes = annodb.keys()
    functional_groups = load_func_assoc(groupfile, transfile)

    #binomial_tests(snp_counts, site_counts)
    #main(site_counts, snp_counts, genes, functional_groups, test_func=binomial_test)
    for name, test_func in [("Fisher test", fisher_test), ("Binomial Test", binomial_test)]:
        print name
        run_tests(site_counts, snp_counts, genes, functional_groups,
                  test_func=test_func)
        print
    ## Just for output for spreadsheet
    #for gene in genes:
        #counts = ka_ks_counts_for_gene_group([gene], snp_counts, site_counts, pseudo=0.0)
        #print "%s,%s" % (gene, ",".join(map(str,counts)))
    #items = [(k, v) for (k,v) in functional_groups.items() if len(v[1]) > 0]    
    #for k, v in items:
        #counts = ka_ks_counts_for_gene_group(v[1], snp_counts, site_counts, pseudo=0.0)
        #print "%s,%s,%s, %s" % (k, v[0], " ".join(v[1]), ",".join(map(str,counts)))
    
def main():
    # entry point for phenoseq_kaks command defined in setup.py
    gbfile, groupfile, transfile, tagFiles = parse_args()
    run_all(gbfile, groupfile, transfile, tagFiles)

if __name__ == '__main__':
    main()
