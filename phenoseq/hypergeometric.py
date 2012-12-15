import analyze
from pathways import load_func_assoc
import sys

from scipy.stats import hypergeom

def phenoseq_top_genes(gbfile, tagFiles):
    annodb, al, dna = analyze.read_genbank_annots(gbfile)
    snps = analyze.read_tag_files(tagFiles)
    results = analyze.analyze_nonsyn(snps, annodb, al, dna)
    return results

def p_value(num_genes, num_genes_int_top_list, num_top_genes, total_genes=4000):
    rv = hypergeom(total_genes, num_top_genes, num_genes)
    p = 1. - rv.cdf(num_genes_int_top_list,)
    return p



def main():
    N = int(sys.argv[1])
    gbfile = sys.argv[2]
    tagFiles = sys.argv[3:]

    top_genes = phenoseq_top_genes(gbfile, tagFiles)
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

if __name__ == '__main__':
    main()
    
    
