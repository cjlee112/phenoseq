
def get_tag(filename):
    return filename.split('.')[0].split('_')[-1]

def get_merge_name(tag, bamlist):
    l = [s.split('_')[2] for s in bamlist]
    l.sort()
    return '_'.join([tag] + l) + '.bam'

def generate_sublists(bamlist):
    'generate leave-one-out jackknife samples'
    for i in range(len(bamlist)):
        yield bamlist[:i] + bamlist[i + 1:]
