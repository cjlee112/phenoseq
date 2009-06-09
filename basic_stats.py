
import random

def sample_genes(ngene=4200, nmut=100):
    """Generate counts of genes with 1, 2 or more mutations,
    treating each gene as equally likely."""
    d = {}
    for i in range(nmut):
        k = random.randint(0, ngene)
        d[k] = d.get(k, 0) + 1
    c = {}
    for i in d.values():
        c[i] = c.get(i, 0) + 1
    return c

