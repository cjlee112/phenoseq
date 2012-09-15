from phenoseq.analyze import *

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
    
def load_func_assoc(filename="16.1/data/func-associations.col"):
    """Loads Ecocyc functional associations file. See example at: http://bioinformatics.ai.sri.com/ptools/flatfile-samples/func-associations.col"""
    pathway_dict = dict()
    trans = load_gene_translations()
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

def load_data(tagFiles, count_syn=True):
    annodb, al, dna = read_genbank_annots("NC_000913.gbk")
    snps = read_tag_files(tagFiles)
    #gsd = GeneSNPDict(tagDict, annodb, al, dna, count_syn=count_syn)
    gsd = map_snps_chrom1(snps, al, dna)
    return (annodb, al, dna, snps, gsd)    

def count_snps_per_gene(snps, al, dna):
    """Counts snps for each gene."""
    rc = dict(A='T', C='G', G='C', T='A')
    snp_counts = {}
    for snp in snps:
        pos = snp.pos
        try:
            geneSNP = al[dna[pos - 1:pos]].keys()[0]
        except IndexError:
            continue
        else:
            Na = 0
            Ns = 0
            # Determine if snp is synonymous or not.
            if geneSNP.orientation < 0:
                geneSNP = -geneSNP
            gene = geneSNP.id
            ipos = geneSNP.start % 3
            codonStart = geneSNP.start - ipos
            codon = str(geneSNP.path.sequence[codonStart:codonStart + 3])
            b = snp.alt # substitution letter
            if geneSNP.sequence.orientation < 0: # must complement
                b = rc[b.upper()]
            codonAlt = codon[:ipos] + b + codon[ipos + 1:]
            snp.aaRef = sequtil.translate_orf(codon)
            snp.aaAlt = sequtil.translate_orf(codonAlt)
            if snp.aaRef == snp.aaAlt:
                Ns += 1
            else:
                Na += 1
        try:
            a, b = snp_counts[gene]
            Na += a
            Ns += b
        except KeyError:
            pass
        snp_counts[gene] = (Na, Ns)
    return snp_counts

def count_snps_per_pathway(snp_counts, pathway_dict):
    pathway_counts = dict()
    for k in pathway_dict.keys():
        count = 0
        genes = pathway_dict[k][1]
        for gene in genes:
            try:
                count += snp_counts[gene][0]
            except KeyError:
                pass
        pathway_counts[k] = count
    return pathway_counts

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
    
    
if __name__ == '__main__':

    # First EXP
    #tagFiles = ['ACAGTG.vcf', 'ACTTGA.vcf', 'ATCACG.vcf', 'CAGATC.vcf', 'CGATGT.vcf', 'CTTGTA.vcf', 'GATCAG.vcf', 'GCCAAT.vcf', 'TGACCA.vcf', 'TTAGGC.vcf']

    # Luisa's EXP
    tagFiles = ["aligned_s_8_%s.vcf" % x for x in ['ATCACG','CGATGT','TTAGGC','TGACCA', 'ACAGTG', 'GCCAAT', 'CAGATC', 'ACTTGA']]

    pathway_dict = load_func_assoc()
        
    (annodb, al, dna, snps, gsd) = load_data(tagFiles)
    snp_counts = count_snps_per_gene(snps, al, dna)
    pathway_counts = count_snps_per_pathway(snp_counts, pathway_dict)
    results = [(float(v) / len(pathway_dict[k][1]),k) for k, v in pathway_counts.items() if len(pathway_dict[k][1])]
    results.sort()
    #for i, (v, k) in list(enumerate(reversed(results)))[:20]:
        #if v > 0:
            #print i, v, k, pathway_dict[k][1]
    
    results = analyze_nonsyn_groups(pathway_dict, snps, annodb, al, dna)
    for k, v in results:
        print "%s,%s,%s" % (k, v, " ".join(pathway_dict[v][1]))

