from phenoseq.analyze import *
import math

import numpy
from matplotlib import pyplot

def normalize_dict(d_):
    s = sum(map(float, d_.values()))
    d = dict()
    for k,v in d_.items():
        d[k] = float(v) / s
    return d

def entropy(d):
    s = 0.
    for k,v in d.items():
        try:
            log = math.log(v)
        except ValueError:
            pass
        else:
            s += v * log
    return -s

def load_data(count_syn=True):
    annodb, al, dna = read_genome_annots("NC_000913.gbk")
    tagFiles = ['ACAGTG.vcf', 'ACTTGA.vcf', 'ATCACG.vcf', 'CAGATC.vcf', 'CGATGT.vcf', 'CTTGTA.vcf', 'GATCAG.vcf', 'GCCAAT.vcf', 'TGACCA.vcf', 'TTAGGC.vcf']
    tagDict = read_tag_files(tagFiles)
    gsd = GeneSNPDict(tagDict, annodb, al, dna, count_syn=count_syn)
    return (annodb, al, dna, tagDict, gsd)

def mutation_counts():
    (annodb, al, dna, tagDict, gsd) = load_data()
    mutation_types = [('AT','GC'), ('GC', 'AT'), ('AT','TA'), ('GC', 'TA'), ('AT', 'CG'), ('GC', 'CG')]
    mutation_map = dict()
    mutations = dict()
    for k in mutation_types:
        a, b = k
        mutation_map[(a[0], b[0])] = k
        mutation_map[(a[1], b[1])] = k
        mutations[k] = []
    prev_counts = dict(zip('ACGT', [0,0,0,0]))
    next_counts = dict(zip('ACGT', [0,0,0,0]))
    # Find and classify mutations
    for v in tagDict.get_snps():
        _, _, j = v
        ori, new = j.ref, j.alt
        pos = j.pos
        prev = str(dna[pos-2:pos-1])
        next = str(dna[pos:pos+1])
        prev_counts[prev] += 1
        next_counts[next] += 1
        mutations[mutation_map[(ori, new)]].append((prev, next))
    for k in mutation_types:
        print k, len(mutations[k])

def NTG_bias():
    """Analyzes bias in context of NTG mutagenesis."""
    (annodb, al, dna, tagDict, gsd) = load_data()
    mutations = list()
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

    ## Find and classify mutations.
    #for k, v in gsd.items():
        #for i, j in v:
            #ori, new = j.ref, j.alt
            #if j.qual < 99:
                #continue
            ## Find surrounding bases.
            #pos = j.pos
            #prev = str(dna[pos-2:pos-1])
            #next = str(dna[pos:pos+1])
            #mutations.append((ori, new, prev, next))

    # Use all SNPs, not just coding region SNPs.
    for v in tagDict.get_snps():
        _, _, j = v
        ori, new = j.ref, j.alt
        pos = j.pos
        prev = str(dna[pos-2:pos-1])
        next = str(dna[pos:pos+1])
        #if j.orientation < 0:
            #prev, next = next, prev
        mutations.append((ori, new, prev, next))
        #prev_counts[prev] += 1
        #next_counts[next] += 1
        #mutations[mutation_map[(ori, new)]].append((prev, next))


    ## Previous Context Counts##
    prev_context_positions = dict()
    dna_str = str(dna[0:-1])
    for i in range(1, len(dna_str)):
        site, prev = dna_str[i],dna_str[i-1]
        try:
            prev_context_positions[(site, prev)] += 1
        except KeyError:
            prev_context_positions[(site, prev)] = 1

    prev_context_counts = dict()

    ## Previous Context
    for context, positions in prev_context_positions.items():
        counts = dict(zip('ACGT', [0,0,0,0]))
        orig, prev = context
        counts[orig] += 32 * positions
        for (orig_, new, prev_, next) in mutations:
            if (orig, prev) == (orig_, prev_):
                counts[orig] -= 1
                counts[new] += 1
        prev_context_counts[context] = counts
    prev_results = dict()
    prev_results_compl = dict()
    total_counts = 0
    for context, dist in prev_context_counts.items():
        orig, prev = context
        sum_ = 0
        for key, value in dist.items():
            sum_ += value
        total_counts += sum_
        for key, value in dist.items():
            #print "P(%s|%s, %s) = %s" % (key, orig, prev, float(value)/sum_)
            prev_results[(orig, key, prev)] = float(value)/sum_
            prev_results_compl[(complement[orig], complement[key], complement[prev])] = float(value)/sum_


    ## Post Context##
    post_context_positions = dict()
    for i in range(0, len(dna_str) - 1):
        site, next = dna_str[i],dna_str[i+1]
        try:
            post_context_positions[(site, next)] += 1
        except KeyError:
            post_context_positions[(site, next)] = 1
    post_context_counts = dict()
    for context, positions in post_context_positions.items():
        counts = dict(zip('ACGT', [0,0,0,0]))
        orig, next = context
        counts[orig] += 32 * positions
        for (orig_, new, prev_, next_) in mutations:
            if (orig, next) == (orig_, next_):
                counts[orig] -= 1
                counts[new] += 1
        post_context_counts[context] = counts
    total_counts = 0
    post_results = dict()
    post_results_compl = dict()
    for context, dist in post_context_counts.items():
        orig, next = context
        sum_ = 0
        for key, value in dist.items():
            sum_ += value
        total_counts += sum_
        for key, value in dist.items():
            post_results[(complement[orig], complement[key], complement[next])] = float(value)/sum_
            #print "P(%s|%s, %s) = %s" % (key, orig, next, float(value)/sum_)
            post_results_compl[(orig, key, next)] = float(value)/sum_

    scale = 1.0e-8

    print "Pre-context"
    for a in 'GATC':
        for b in 'ATCG':
            if a == b:
                continue
            for c in 'ATCG':
                key = (a,b,c)
                print prev_results[key]/scale, post_results[key]/scale,
            print ""

    print "Post-context"
    for a in 'GATC':
        for b in 'ATCG':
            if a == b:
                continue
            for c in 'ATCG':
                key = (a,b,c)
                print post_results_compl[key]/scale, prev_results_compl[key]/scale,
            print ""

    return prev_results

#def codon_transitions(tagDict, annodb, al, dna, count_syn=False):
    ## Adapted from GeneSNPDict code.
    #rc = dict(A='T', C='G', G='C', T='A')
    #moves = []
    #for pos, tag, snp in tagDict.get_snps():
        #try:
            #geneSNP = al[dna[pos - 1:pos]].keys()[0]
        #except IndexError:
            #pass
        #else:
            #if geneSNP.orientation < 0:
                #geneSNP = -geneSNP
            #gene = geneSNP.id
            #ipos = geneSNP.start % 3
            #codonStart = geneSNP.start - ipos
            #codon = str(geneSNP.path.sequence[codonStart:codonStart + 3])
            #b = snp.alt # substitution letter
            #if geneSNP.sequence.orientation < 0: # must complement
                #b = rc[b.upper()]
            #codonAlt = codon[:ipos] + b + codon[ipos + 1:]
            #snp.aaRef = sequtil.translate_orf(codon)
            #snp.aaAlt = sequtil.translate_orf(codonAlt)
            #moves.append((snp.aaRef, snp.aaAlt))
    #return moves

def codon_transition_counts():
    (annodb, al, dna, tagDict, gsd) = load_data()
    l = codon_transitions(tagDict, annodb, al, dna, count_syn=True)
    d = {'null':0}
    for a,b in l:
        if a == b:
            d['null'] += 1
            continue
        try:
            d[b] += 1
        except KeyError:
            d[b] = 1
    print d

def spontaneous_mutations():
    d = {('AT','GC'): 0.633, ('GC', 'AT'): 0.326, ('AT','TA'): 0., ('GC', 'TA'): 0.0247, ('AT', 'CG'): 0.011, ('GC', 'CG'): 0.0055}
    return d

def empirical_mutations():
    d = {('AT','GC'): 89./4099, ('GC', 'AT'): 3960.4099, ('AT','TA'): 3./4099, ('GC', 'TA'): 3./4099, ('AT', 'CG'): 19./4099, ('GC', 'CG'): 25./4099}
    return d

def aa_entropy(mutation_dist, codon_dist=None, mu=1.):
    """Computes entropy of amino acid distribution resulting from applying mutation_dist to each of the possible codons for mutation rate mu"""
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
        codon_transitions[codon] = []
        for i in [0,1,2]:
            for k,v in mutation_dist.items():
                for j in [0,1]:
                    if codon[i] == k[0][j]:
                        result_codon = list(codon)
                        result_codon[i] = k[1][j]
                        result_codon = "".join(result_codon)
                        # 1/3 assumes mutation is equally likely at each position in the codon.
                        codon_transitions[codon].append((result_codon, mu * 1./3 * normed_mutations[k[0]][k[1]]))
        codon_transitions[codon].append((codon, 1.-mu))



    ## Entropy for amino acid distribution for each codon.
    entropies = {}
    for codon in codons:
        #print codon,
        aa_dist = {}
        for (new_codon, p) in codon_transitions[codon]:
            new_aa = sequtil.translate_orf(new_codon)
            try:
                aa_dist[new_aa] += p
            except KeyError:
                aa_dist[new_aa] = p
        entropies[codon] = entropy(normalize_dict(aa_dist))

    # Weighted local entropies.
    total_entropy = 0
    for codon in codons:
        try:
            total_entropy += codon_dist[codon] * entropies[codon]
        except KeyError:
            # If codon_dist lacks a codon, the contribution to the some is zero, so just pass.
            pass
    return total_entropy


    ## Total entropy over all codons, weighted by codon_dist
    #aa_results = {}
    #for codon in codons:
        #for (new_codon, p) in codon_transitions[codon]:
            #new_aa = sequtil.translate_orf(new_codon)
            #try:
                #aa_results[new_aa] += codon_dist[codon] * p
            #except KeyError:
                #aa_results[new_aa] = codon_dist[codon] * p
    #return entropy(normalize_dict(aa_results))

def empirical_codon_distribution():
    annodb, al, dna = read_genome_annots("NC_000913.gbk")
    codons = {}
    for k in annodb.keys():
        seq = str(annodb[k].sequence)
        if not len(seq) % 3:
            continue
        for i in range(0,len(seq), 3):
            codon = seq[i:i+3]
            try:
                codons[codon] += 1
            except KeyError:
                codons[codon] = 1
    return normalize_dict(codons)        

def NTG_codon_hit_distribution():
    # Adapted from GeneSNPDict code.
    (annodb, al, dna, tagDict, gsd) = load_data()
    rc = dict(A='T', C='G', G='C', T='A')
    codons = {}
    for pos, tag, snp in tagDict.get_snps():
        try:
            geneSNP = al[dna[pos - 1:pos]].keys()[0]
        except IndexError:
            pass
        else:
            if geneSNP.orientation < 0:
                geneSNP = -geneSNP
            gene = geneSNP.id
            ipos = geneSNP.start % 3
            codonStart = geneSNP.start - ipos
            codon = str(geneSNP.path.sequence[codonStart:codonStart + 3])
            try:
                codons[codon] += 1
            except KeyError:
                codons[codon] = 1.          
    return normalize_dict(codons)

def GC_bias_codon_dist(gc_content=None):
    if not gc_content:
        annodb, al, dna = read_genome_annots("NC_000913.gbk")
        counts = {'GC':0., 'AT':0.}
        # Option to restrict to coding regions
        #for k in annodb.keys():
            #seq = str(annodb[k].sequence)
            #if not len(seq) % 3:
                #continue
            #for s in seq:
                #if s in 'GC':
                    #counts['GC'] += 1
                #else:
                    #counts['AT'] += 1
        for s in str(dna):
            if s in 'GC':
                counts['GC'] += 1
            else:
                counts['AT'] += 1
        GC_dist = normalize_dict(counts)
    else:
        GC_dist = {'GC':gc_content, 'AT':1-gc_content}
    codon_dist = {}
    for a in 'ATCG':
        for b in 'ATCG':
            for c in 'ATCG':
                weight = 1.
                codon = a + b+ c
                for x in [a,b,c]:
                    for key in GC_dist.keys():
                        if x in key:
                            weight *= GC_dist[key]
                            break
                codon_dist[codon] = weight
    return normalize_dict(codon_dist)

def empirical_context_bias():
    (annodb, al, dna, tagDict, gsd) = load_data()
    mutations = list()
    # Find and classify mutations.
    for k, v in gsd.items():
        for i, j in v:
            ori, new = j.ref, j.alt
            if j.qual < 99:
                continue
            # Find surrounding bases.
            pos = j.pos
            prev = str(dna[pos-2:pos-1])
            next = str(dna[pos:pos+1])
            mutations.append((ori, new, prev, next))
            
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

    #ntg_context_bias = dict()
    #for prev in 'ATCG':
        #for orig in 'ATCG':
            #for next in 'ATCG':
                #d = {}
                #for new in 'ATCG':
                    #if new == orig:
                        #continue
                    #d[new] = float(len([(ori_, new_, prev_, next_) for (ori_, new_, prev_, next_) in mutations if (ori_, new_, prev_, next_) == (orig, new, prev, next)]))                
                    #d[new] += float(len([(ori_, new_, prev_, next_) for (ori_, new_, prev_, next_) in mutations if (ori_, new_, prev_, next_) == (complement[orig], new, complement[next], complement[prev])]))                    
                #if not sum(d.values()):
                    #ntg_context_bias[(prev, orig, next)] = d
                #else:
                    #ntg_context_bias[(prev, orig, next)] = normalize_dict(d)
    #return ntg_context_bias

    ntg_context_bias = dict()
    for prev in 'ATCG':
        for orig in 'ATCG':
            d = {}
            for new in 'ATCG':
                if new == orig:
                    continue
                d[new] = float(len([(ori_, new_, prev_, next_) for (ori_, new_, prev_, next_) in mutations if (ori_, new_, prev_,) == (orig, new, prev)]))                
                #d[new] += float(len([(ori_, new_, prev_, next_) for (ori_, new_, prev_, next_) in mutations if (ori_, new_, prev_) == (complement[orig], complement[new], complement[next])]))                    
            if not sum(d.values()):
                ntg_context_bias[(prev, orig)] = d
            else:
                ntg_context_bias[(prev, orig)] = normalize_dict(d)
    return ntg_context_bias


#def aa_entropy_context_bias(context_dist=None):
    #(annodb, al, dna, tagDict, gsd) = load_data()

    ## Compute positional entropies based on context bias.
    ## Iterate over coding regions, over codons, computing amino acid products for codon based on context
    #dna = str(dna)
    #entropies = []

    #for k in annodb.keys():
        #handle = annodb[k].sequence
        #if handle.orientation < 0:
            #seq = str(dna[handle.start+1:handle.stop-1])
        #else:
            #seq = str(dna[handle.start-1:handle.stop+1])

        #if not (len(seq) % 3 == 2):
            #continue
        #for i in range(1, len(seq) - 3, 3):
            #d = {}
            #codon = seq[i:i+3]
            #for j in range(i, i+3):
                #prev = seq[j-1]
                #next = seq[j+1]
                #orig = seq[j]
                #for new in 'ACTG':
                    #if new == orig:
                        #continue                    
                    #m = ntg_context_bias[(prev, orig)][new]
                    #result_codon = list(codon)
                    #result_codon[j-i] = new
                    #result_codon = "".join(result_codon)
                    #try:
                        #d[result_codon] += m
                    #except KeyError:
                        #d[result_codon] = m
            #if not sum(d.values()):
                #entropies.append(0)
                #continue
            #e = entropy(normalize_dict(d))
            #entropies.append(e)
            ##print len(seq), i, e
    #return sum(entropies) / len(entropies)

def gc_content_plot(mutation_dist, N=100):
    points = []
    for i in range(1, N):
        x = float(i) / N
        GCcd = GC_bias_codon_dist(x)
        e = aa_entropy(mutation_dist=mutation_dist, codon_dist=GCcd)
        points.append((i,e))
    
    #pyplot.clf()
    pyplot.plot([x for (x,y) in points], [y for (x,y) in points])


#def get_all_snps(annodb, al, dna, tagDict):
    
    #for v in tagDict.get_snps():
        #_, _, j = v
        #ori, new = j.ref, j.alt
        #pos = j.pos
        #prev = str(dna[pos-2:pos-1])
        #next = str(dna[pos:pos+1])
        ##if j.orientation < 0:
            ##prev, next = next, prev
        ##mutations.append((ori, new, prev, next))

    
    
    #mutations = list()
    ## Find and classify mutations.
    #for k, v in gsd.items():
        #for i, j in v:
            #ori, new = j.ref, j.alt
            #if j.qual < 99:
                #continue
            ## Find surrounding bases.
            #pos = j.pos
            #prev = str(dna[pos-2:pos-1])
            #next = str(dna[pos:pos+1])
            #mutations.append((ori, new, prev, next))
            
    #complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}


def num_strains_per_library(tag):
    if tag in ['ACAGTG', 'CTTGTA']:
        return 4
    return 3

def rows_to_rst_list_table(rows):
    rst_str = ""
    for row in rows:
        rst_str += "* - %s \n" % row[0]
        for elem in row[1:]:
            rst_str += "  - %s \n" % elem
    return rst_str

def mutations_per_strain():
    (annodb, al, dna, tag_dict, gsd) = load_data()
    tags = tag_dict.keys()
    mutation_types = [('AT','GC'), ('GC', 'AT'), ('AT','TA'), ('GC', 'TA'), ('AT', 'CG'), ('GC', 'CG')]
    mutation_map = dict()
    mutations = dict(zip(tags, [dict() for _ in range(len(tags))]))
    for k in mutation_types:
        a, b = k
        mutation_map[(a[0], b[0])] = k
        mutation_map[(a[1], b[1])] = k
        for tag in tags:
            mutations[tag][k] = []
    #prev_counts = dict(zip('ACGT', [0,0,0,0]))
    #next_counts = dict(zip('ACGT', [0,0,0,0]))
    from collections import defaultdict
    positions = defaultdict(list)
    # Find and classify mutations
    for pos, tag, snp in tag_dict.get_snps():
        j = snp
        ori, new = j.ref, j.alt
        pos = j.pos
        prev = str(dna[pos-2:pos-1])
        next = str(dna[pos:pos+1])
        #prev_counts[prev] += 1
        #next_counts[next] += 1
        mutations[tag][mutation_map[(ori, new)]].append((prev, next))
        positions[pos].append(tag)
    # Print table by tag for latex
    rows = []
    row = ['Library']
    row.extend(tags)
    rows.append(row)
    for k in mutation_types:
        row = ["-".join(k)]
        for tag in tags:
            v = len(mutations[tag][k])/float(num_strains_per_library(tag))
            row.append(round(v,2))
        mean = numpy.mean(map(float, row[1:]))
        std = numpy.std(map(float, row[1:]))
        row.extend([round(mean,2), round(std,2)])
        rows.append(row)
    # Transition transversion ratio
    row = ['TR/TV ratio']
    for tag in tags:
        ratio = float( len(mutations[tag][('GC', 'AT')]) + len(mutations[tag][('AT', 'GC')])) / float( len(mutations[tag][('GC', 'TA')]) + len(mutations[tag][('GC', 'CG')]) + len(mutations[tag][('AT', 'TA')]) + len(mutations[tag][('AT', 'CG')]))
        row.append(round(ratio,2))
    mean = numpy.mean(map(float, row[1:]))
    std = numpy.std(map(float, row[1:]))
    row.extend([round(mean,2), round(std,2)])
    rows.append(row)
    # Transition ratio.
    row = ['Transition ratio']
    for tag in tags:
        ratio = float( len(mutations[tag][('GC', 'AT')])) / len(mutations[tag][('AT', 'GC')])
        row.append(round(ratio,2))
    mean = numpy.mean(map(float, row[1:]))
    std = numpy.std(map(float, row[1:]))
    row.extend([round(mean,2), round(std,2)])
    rows.append(row)

    rst = rows_to_rst_list_table(rows)
    print rst


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


#def empirical_codon_distribution():
    #annodb, al, dna = read_genome_annots("NC_000913.gbk")
    #codons = {}
    #for k in annodb.keys():
        #seq = str(annodb[k].sequence)
        #if not len(seq) % 3:
            #continue
        #for i in range(0,len(seq), 3):
            #codon = seq[i:i+3]
            #try:
                #codons[codon] += 1
            #except KeyError:
                #codons[codon] = 1
    #return normalize_dict(codons)  




    ### Count codons.
    ##annodb, al, dna = read_genome_annots("NC_000913.gbk")
    ##codons = {}
    ##for k in annodb.keys():
        ##seq = str(annodb[k].sequence)
        ##for i in range(0,len(seq), 3):
            ##codon = seq[i:i+3]
            ##try:
                ##codons[codon] += 1
            ##except KeyError:
                ##codons[codon] = 1

def count_codons(seq=None):
    # Count codons.
    codons = {}
    for i in range(0,len(seq), 3):
        codon = seq[i:i+3]
        try:
            codons[codon] += 1
        except KeyError:
            codons[codon] = 1
    return codons

def count_sites(codons, d=sites_count_dict()):
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

if __name__ == '__main__':
    #annodb, al, dna = read_genome_annots("NC_000913.gbk")
    #seq = str(dna)
    #codons = count_codons(seq)
    #print count_sites(codons)

    #codon_transitions = codon_weighting(spontaneous_mutations())
    #for k,v in d.items():
        #print k, v
    #exit()


    (annodb, al, dna, tagDict, gsd) = load_data()
    rc = dict(A='T', C='G', G='C', T='A')

    # Count sites for each gene
    site_counts = {}
    for k in annodb.keys():
        seq = str(annodb[k].sequence)
        codons = count_codons(seq)
        (Na, Ns) = count_sites(codons)
        site_counts[k] = (len(seq), Na, Ns)
    #for x in results:
        #print x

    # Count snps for each gene
    snp_counts = {}
    for pos, tag, snp in tagDict.get_snps():
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

    import rpy2.robjects as robjects
    robjects.r.options(warn=-1) # supress warnings
    fisher_test = robjects.r['fisher.test']
    matrix = robjects.r['matrix']

    sum_ga = 0
    sum_gs = 0
    sum_na = 0
    sum_ns = 0
    for gene, counts in snp_counts.items():
        Ga, Gs = counts
        sum_ga += Ga
        sum_gs += Gs
    for gene, counts in site_counts.items():
        l, Na, Ns = counts
        sum_na += Na
        sum_ns += Ns
    print 'global', (float(sum_ga) / sum_na) / (float(sum_gs) / sum_ns), sum_ga, sum_na, sum_gs, sum_ns, len(dna),
    Ga, Gs, Na, Ns = sum_ga, sum_gs, sum_na, sum_ns
    Ga_exp = Na
    Gs_exp = Ns
    v = robjects.IntVector([Ga, Gs, int(Ga_exp+0.5), int(Gs_exp+0.5)]) 
    fet = fisher_test(matrix(v, nrow=2), alternative="greater")
    print Ga_exp, Gs_exp, fet[0][0], fet[1][0], fet[1][1] # p-value, confidence interval (95%)


    for gene, counts in snp_counts.items():
        Ga, Gs = counts
        l, Na, Ns = site_counts[gene]
        print gene,
        if Gs > 0:
            print (float(Ga) / Na) / (float(Gs) / Ns),
        else:
            print '----',
        print Ga, Na, Gs, Ns, l,
        #Ga_exp = (Ga + Gs) * Na / (Na + Ns)
        #Gs_exp = (Ga + Gs) * Ns / (Na + Ns)
        Ga_exp = Na
        Gs_exp = Ns
        
        v = robjects.IntVector([Ga, Gs, int(Ga_exp+0.5), int(Gs_exp+0.5)]) 
        fet = fisher_test(matrix(v, nrow=2), alternative="greater")
        print Ga_exp, Gs_exp, fet[0][0], fet[1][0], fet[1][1] # p-value, confidence interval (95%)

    

    

    #mutations_per_strain()
    #exit()
        
    #for tag in tags:
        #for k in mutation_types:
            #print k, len(mutations[tag][k])
        #print float(len(mutations[tag][('GC', 'AT')])) / len(mutations[tag][('AT', 'GC')])
    #for k, vs in positions.items():
        #if len(vs) > 2:
            #print k, len(vs)
    


    
    #ntg_context_bias = empirical_context_bias()
    #print ntg_context_bias
    #print aa_entropy_context_bias(context_dist=ntg_context_bias)
    
    ## Local entropy calculations
    #ecd = empirical_codon_distribution()
    #NTGcd = NTG_codon_hit_distribution()
    #GCcd = GC_bias_codon_dist()
    #for (name, mutation_dist) in [('Spontaneous', spontaneous_mutations()), ('NTG derived', empirical_mutations())]:
        #print name, 'mutation spectrum'
        #for (name2, codon_dist) in [('Uniform', None), ('Genome Dist', ecd), ('GC bias', GCcd), ('NTG hits', NTGcd)]:
            #print ' ', name2,
            #print aa_entropy(mutation_dist=mutation_dist, codon_dist=codon_dist)

    ## GC bias plot
    #for (name, mutation_dist) in [('Spontaneous', spontaneous_mutations()), ('NTG derived', empirical_mutations())]:
        #gc_content_plot(mutation_dist)
    #pyplot.ylabel('Entropy')
    #pyplot.xlabel('GC content %')
    #pyplot.title('blue: spontaneous -- grean: NTG')
    #pyplot.show()
    
    #NTG_bias()
    
    # Locational entropy
    


