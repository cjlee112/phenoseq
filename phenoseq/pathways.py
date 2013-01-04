from optparse import OptionParser

import analyze

## For pathways file... here in case it's useful later.
#def load_pathways(filename):
    #"""Loads Ecocyc pathway file. See example at: http://bioinformatics.ai.sri.com/ptools/flatfile-samples/pathways.col"""
    ##pathways = []
    #pathway_dict = dict()
    #with open(filename) as handle:
        #for line in handle:
            ## Toss comments.
            #if line.startswith('#'):
                #continue
            ## Toss headers.
            #skip = False
            #for h in ['UNIQUE', 'NAME', 'GENE', 'EG']:
                #if line.startswith(h):
                    #skip = True
                    #break
            #if skip:
                #continue
            ## Parse lines.
            ##print line
            #split = line.strip().split('\t')
            #genes = []
            ## ID, Name, genes by four letter name. Gene IDs end the list but are not needed.
            #genes = list(sorted([x for x in split[2:] if len(x) in [3,4]]))
            #pathway_dict[split[0]] = (split[1], genes)
            ##p = Pathway(split[0], split[1], [x for x in split[2:] if len(x) in [3,4]])
            ##pathways.append(p)
    #return pathway_dict


def load_gene_translations(filename="16.1/data/genes.col"):
    """Loads Ecocyc genes file. Needed to process other files (to get the correct gene names). See example at: http://bioinformatics.ai.sri.com/ptools/flatfile-samples/genes.col"""
    #pathways = []
    d = dict()
    with open(filename) as handle:
        for line in handle:
            # Toss comments.
            if line.startswith('#'):
                continue
            # Toss headers.
            skip = False
            for h in ['UNIQUE']:
                if line.startswith(h):
                    skip = True
                    break
            if skip:
                continue
            # Parse lines.
            split = line.strip().split('\t')
            d[split[1]] = split[2]
    return d    
    

def load_func_assoc(filename="16.1/data/func-associations.col",
                    transfile="16.1/data/genes.col"):
    """Loads Ecocyc functional associations file. See example at: http://bioinformatics.ai.sri.com/ptools/flatfile-samples/func-associations.col"""
    pathway_dict = dict()
    trans = load_gene_translations(transfile)
    with open(filename) as handle:
        for line in handle:
            # Toss comments.
            if line.startswith('#'):
                continue
            # Toss headers.
            skip = False
            for h in ['UNIQUE', 'NAME', 'GENE', 'EG']:
                if line.startswith(h):
                    skip = True
                    break
            if skip:
                continue
            # Parse lines.
            #print line
            split = line.strip().split('\t')
            #genes = []
            # ID, Name, genes by four letter name. Gene IDs end the list but are not needed.
            translated_genes = []
            genes = list(sorted([x for x in split[2:] if x.startswith('b')]))
            for gene in genes:
                try:
                    mapped = trans[gene]
                    translated_genes.append(mapped)
                except KeyError:
                    pass
            
            pathway_dict[split[0]] = (split[1], translated_genes)
            #p = Pathway(split[0], split[1], [x for x in split[2:] if len(x) in [3,4]])
            #pathways.append(p)
    return pathway_dict    


def promoter_snps(promoter_offset=500):
    tagFiles = ["aligned_s_8_%s.vcf" % x for x in ['ATCACG','CGATGT','TTAGGC','TGACCA', 'ACAGTG', 'GCCAAT', 'CAGATC', 'ACTTGA']]
    (annodb, al, dna, snps, gsd) = load_data(tagFiles)
    snp_locs = [snp.pos for snp in snps]
    results = []
    for gene in annodb.keys():
        seq = annodb[gene].sequence
        if seq.orientation > 0:
            start = seq.start
            l = [s for s in snp_locs if (s >= start - promoter_offset) and (s < start)]
        else:
            #stop = seq.stop
            #print seq.start, seq.stop
            stop = -seq.start
            l = [s for s in snp_locs if (s > stop) and (s <= stop + promoter_offset)]
        if len(l) > 0:
            results.append((gene, len(l)))
            #print gene, len(l), seq.orientation
    return results

def parse_args():
    usage = "usage: %prog [options] vcf_files"
    parser = OptionParser(usage=usage)
    parser.add_option("-g", "--genbank", dest="gbfile",
                    help="Genbank filename of reference genome")
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

    for x in [gbfile, groupfile, transfile]:
        if not x:
            print "Missing required files. Try --help for more information."
            exit()
    return gbfile, groupfile, transfile, tagFiles
    
    
def pathways_cmd(gbfile, groupfile, transfile, tagFiles):
    pathway_dict = load_func_assoc(groupfile, transfile)
    for k,v in pathway_dict.items():
        pathway_dict[k] = v[1] # only keep the gene list
        
    annodb, al, dna = analyze.read_genbank_annots(gbfile)
    snps = analyze.read_tag_files(tagFiles)
    results = analyze.analyze_nonsyn_groups(pathway_dict, snps, annodb,
                                            al, dna)
    for k, v in results:
        print "%s,%s,%s" % (k, v, " ".join(pathway_dict[v]))

def main():
    # entry point for phenoseq_pathways command defined in setup.py
    (gbfile, groupfile, transfile, tagFiles) = parse_args()
    pathways_cmd(gbfile, groupfile, transfile, tagFiles)

if __name__ == '__main__':
    main()
