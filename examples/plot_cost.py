from matplotlib import pyplot
try:
    from phenoseq import simulate
except ImportError:
    import util
    util.add_root_path()
    from phenoseq import simulate

# run like this to have interactive control over resulting plot:
# python -i plot_cost.py

class FunctionArgs(object):
    'wrapper for passing additional args to a map function'
    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs
    def __call__(self, v):
        return self.func(v, *self.args, **self.kwargs)


def onelib_cost_fig(nstrains=range(3, 34, 3),
                    xlabel=r'experiment cost (\$)',
                    ylabel='average number real targets identified',
                    plotargs={}, mapFunc=map, **kwargs):
    'plot experiment cost vs. yield over a range of nstrain values'
    l = []
    costs = []
    func = FunctionArgs(simulate.min_cost, **kwargs)
    for cost, nlib, cov, y in mapFunc(func, nstrains):
        l.append(y)
        costs.append(cost)
    pyplot.plot(costs, l, **plotargs)
    pyplot.xlabel(xlabel)
    pyplot.ylabel(ylabel)
    return l

def onelib_model_cost_fig(*args, **kwargs):
    onelib_cost_fig(plotargs=dict(marker='+'), **kwargs)
    onelib_cost_fig(nmut=20, plotargs=dict(marker='o'), **kwargs)
    onelib_cost_fig(nmut=100, plotargs=dict(marker='^'), **kwargs)

if __name__ == '__main__':
    from multiprocessing import Pool
    pool = Pool() # use all available CPUs
    onelib_model_cost_fig(mapFunc=pool.map)
    pyplot.show()
