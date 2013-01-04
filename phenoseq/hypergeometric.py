from math import isnan
from optparse import OptionParser
import warnings
from scipy.stats import hypergeom
import analyze
from pathways import load_func_assoc

def phenoseq_top_genes(gbfile, tagFiles):
    annodb, al, dna = analyze.read_genbank_annots(gbfile)
    snps = analyze.read_tag_files(tagFiles)
    results = analyze.analyze_nonsyn(snps, annodb, al, dna)
    return results

def p_value(num_genes, num_genes_int_top_list, num_top_genes, total_genes=4000):
    rv = hypergeom(total_genes, num_top_genes, num_genes)
    p = rv.sf(num_genes_int_top_list - 1)
    if isnan(p): # old version of hypergeom.sf() gives NaN, yuck
        p = rv.pmf(range(num_genes_int_top_list, num_genes + 1)).sum()
    return p

########################################################################
# command line interface

def hypergeom_cmd(gbfile, groupfile, transfile, tagFiles, N=30):
    top_genes = phenoseq_top_genes(gbfile, tagFiles)
    pathway_dict = load_func_assoc(groupfile, transfile)
    top_genes_subset = [y for (x,y) in top_genes[:N]]
    #print top_genes_subset
    results = []
    for name, genes_ in pathway_dict.items():
        genes = genes_[1]
        num_genes_int_top_list = len([g for g in genes if g in top_genes_subset])
        
        if num_genes_int_top_list:
            pval = p_value(len(genes), num_genes_int_top_list,
                           len(top_genes_subset))
            if isnan(pval):
                warnings.warn('ignoring invalid NaN pvalues...')
            else:
                results.append( (len(pathway_dict) * pval, name,
                                 len(genes), genes))
    results.sort()
    for p, name, n, genes in results:
        print ",".join(map(str, [p, name, n, " ".join(genes)]))

def parse_opt():
    usage = "usage: %prog [options] vcf_files"
    parser = OptionParser(usage=usage)
    parser.add_option("-g", "--genbank", dest="gbfile",
                    help="Genbank filename of reference genome")
    parser.add_option("-n", dest="N", default=30,
                    help="Number of genes in top list to test for enrichment")
    parser.add_option("-f", "--functionals", dest="groupfile",
                    help="EcoCyc functionally associated groups filename")
    parser.add_option("-t", "--translation", dest="transfile",
                    help="EcoCyc genes filename for gene id translation")
    (options, tagFiles) = parser.parse_args()

    if len(tagFiles) < 1:
        parser.error("vcf files are required")

    gbfile = options.gbfile
    groupfile = options.groupfile
    transfile = options.transfile
    N = int(options.N)

    for x in [gbfile, groupfile, transfile]:
        if not x:
            print "Missing required files. Try --help for more information."
            exit()
    return gbfile, groupfile, transfile, tagFiles, N

def main():
    # entry point for phenoseq_hypergeom command defined in setup.py
    gbfile, groupfile, transfile, tagFiles, N = parse_opt()
    main(gbfile, groupfile, transfile, tagFiles, N)

if __name__ == '__main__':
    main()

