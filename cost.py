
def get_nstrain(tag):
    'return #strains in each pool'
    if tag == 'ACAGTG' or tag == 'CTTGTA':
        return 4
    else:
        return 3

def sum_nstrain(tags):
    return sum([get_nstrain(tag) for tag in tags])

def sum_nhit(besthits, trueHits=('acrB', 'marC', 'acrA')):
    return sum([1 for t in besthits if t[1] in trueHits])

def nstrain_nhit_dict(subsetDict, nbest=9, results=None):
    if results is None:
        results = {}
    for tags,hits in subsetDict.items():
        nstrain = sum_nstrain(tags)
        nhit = sum_nhit(hits[:nbest])
        results.setdefault(nstrain, []).append(nhit)
    return results

def nstrain_nhit_means(d):
    l = []
    for nstrain, nhits in d.items():
        l.append((nstrain, sum(nhits) / float(len(nhits))))
    l.sort()
    return l

