import json

fn = 'DOPS/20190723_006.dat'

d = json.load(open(fn, 'r'))

# d['data'] =  list containing data
# d['data']['pvals'] = second harmonic offset, 'perrs' is error bar
d['params'] = parameters
