import analyze
from math import log, exp, sqrt
from pathways import count_snps_per_gene, load_data, load_func_assoc

from scipy.stats import hypergeom

def phenoseq_top_genes(tagFiles):
    (annodb, al, dna, snps, gsd) = load_data(tagFiles)
    results = analyze.analyze_nonsyn(snps, annodb, al, dna)
    return results

def p_value(num_genes, num_genes_int_top_list, num_top_genes, total_genes=4000):
    rv = hypergeom(total_genes, num_top_genes, num_genes)
    p = 1. - rv.cdf(num_genes_int_top_list,)
    return p

if __name__ == '__main__':
    try:
        N = int(sys.argv[1])
    except IndexError:
        N = 50
    tagFiles = ["aligned_s_8_%s.vcf" % x for x in ['ATCACG','CGATGT','TTAGGC','TGACCA', 'ACAGTG', 'GCCAAT', 'CAGATC', 'ACTTGA']]
    top_genes = phenoseq_top_genes(tagFiles)
    pathway_dict = load_func_assoc()
    top_genes_subset = [y for (x,y) in top_genes[:N]]
    #print top_genes_subset
    results = []
    for name, genes_ in pathway_dict.items():
        genes = genes_[1]
        num_genes_int_top_list = len([g for g in genes if g in top_genes_subset])
        
        if num_genes_int_top_list:
            results.append( (len(pathway_dict)*p_value(len(genes), num_genes_int_top_list, len(top_genes_subset)), name, len(genes), genes))
    results.sort()
    for p, name, n, genes in results:
        print ",".join(map(str, [p, name, n, " ".join(genes)]))

    
    
