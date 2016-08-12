import json
from copy import deepcopy

blind = True
micro = True

def getYields(analysis,mode,mass):
    with open('jsons/{0}/{1}/{2}.json'.format(analysis,mode,mass)) as data_file:    
        data = json.load(data_file)
    return data

class YieldsTex:
    table='''
\\subsection{{{tex} yields}}
{{\\tiny
\\begin{{longtable}}{{ll|r|rr|r}}
\\hline
$\\mHpmpm$ [\\GeV] & Channel & $N_{{exp}}$ & $N_{{AP}}$ & $N_{{PP}}$ & $N_{{obs}}$ \\\\
\\hline
\\endhead
{rows}
\\end{{longtable}}
}}
'''


    row = '''
{mass:4} & ${reco:18}$ & {expected:25} & {ap:25} & {pp:25} & {observed:4} \\\\'''


benchmarks = {
    #'ee100': {'tex': '100\\% $\\Phi \\rightarrow \\Pe\\Pe$',},
    #'em100': {'tex': '100\\% $\\Phi \\rightarrow \\Pe\\Pgm$',},
    #'mm100': {'tex': '100\\% $\\Phi \\rightarrow \\Pgm\\Pgm$',},
    #'et100': {'tex': '100\\% $\\Phi \\rightarrow \\Pe\\Pgt$',},
    #'mt100': {'tex': '100\\% $\\Phi \\rightarrow \\Pgm\\Pgt$',},
    #'tt100': {'tex': '100\\% $\\Phi \\rightarrow \\Pgt\\Pgt$',},
    #'BP1'  : {'tex': 'Benchmark 1',},
    #'BP2'  : {'tex': 'Benchmark 2',},
    #'BP3'  : {'tex': 'Benchmark 3',},
    #'BP4'  : {'tex': 'Benchmark 4',},
    'ee100': {'tex': '100\\% electron-electron',},
    'em100': {'tex': '100\\% electron-muon',},
    'mm100': {'tex': '100\\% muon-muon',},
    'et100': {'tex': '100\\% electron-tau',},
    'mt100': {'tex': '100\\% muon-tau',},
    'tt100': {'tex': '100\\% tau-tau',},
    'BP1'  : {'tex': 'Benchmark 1',},
    'BP2'  : {'tex': 'Benchmark 2',},
    'BP3'  : {'tex': 'Benchmark 3',},
    'BP4'  : {'tex': 'Benchmark 4',},

}

bps = ['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4']
masses = [200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]
flavorMap = {
    'e': '\\Pe',
    'm': '\\Pgm',
    't': '\\Pgt',
}

cats = ['I', 'II', 'III', 'IV', 'V', 'VI']
catTex = {
    'Hpp3l' : {
        'I'   : '\\lll',
        'II'  : '\\llt',
        'III' : '\\ltl',
        'IV'  : '\\ltt',
        'V'   : '\\ttl',
        'VI'  : '\\ttt',
    },
    'Hpp4l' : {
        'I'   : '\\llll',
        'II'  : '\\lllt',
        'III' : '\\lltt',
        'IV'  : '\\ltlt',
        'V'   : '\\lttt',
        'VI'  : '\\tttt',
    },
}
catChans = {
    'Hpp3l' : {
        'I'   : ['eee','eem','eme','emm','mme','mmm'],
        'II'  : ['eet','emt','mmt'],
        'III' : ['ete','etm','mte','mtm'],
        'IV'  : ['ett','mtt'],
        'V'   : ['tte','ttm'],
        'VI'  : ['ttt'],
    },
    'Hpp4l' : {
        'I'   : ['eeee','eeem','eemm','emee','emem','emmm','mmee','mmem','mmmm'],
        'II'  : ['eeet','eemt','emet','emmt','mmet','mmmt','etee','etem','etmm','mtee','mtem','mtmm'],
        'III' : ['eett','emtt','mmtt','ttee','ttem','ttmm'],
        'IV'  : ['etet','etmt','mtet','mtmt'],
        'V'   : ['ettt','mttt','ttet','ttmt'],
        'VI'  : ['tttt'],
    },
}

def latex_float(f):
    f = '{0:g}'.format(f)
    f = f.split('e')
    if len(f)==1: return '${0}$'.format(*f)
    return '${0} \\times 10^{{{1}}}$'.format(f[0],int(f[1]))

print '\\section{Yields\\label{sec:yields}}'
print ''


for bp in bps:
    rowstring = ''
    for mass in masses:
         hpp3l = getYields('Hpp3l',bp,mass)
         hpp4l = getYields('Hpp4l',bp,mass)
         first = True
         if micro:
             rowargs = {'mass':mass if first else '', 'reco':'3\\ell', 'expected':0., 'ap':0., 'pp':0., 'observed': '---' if blind else 0}
             include = False
             for chan in sorted(hpp3l):
                 rowargs['expected'] += hpp3l[chan]['expected']
                 rowargs['ap'] += hpp3l[chan]['ap']
                 rowargs['pp'] += hpp3l[chan]['pp']
                 if not blind: rowargs['observed'] += int(hpp3l[chan]['observed'])
                 first = False
                 include = True
             rowargs['expected'] = latex_float(rowargs['expected'])
             rowargs['ap'] = latex_float(rowargs['ap'])
             rowargs['pp'] = latex_float(rowargs['pp']) if rowargs['pp'] else '---'
             if include: rowstring += YieldsTex.row.format(**rowargs)
             rowargs = {'mass':mass if first else '', 'reco':'4\\ell', 'expected':0., 'ap':'---', 'pp':0., 'observed': '---' if blind else 0}
             include = False
             for chan in sorted(hpp4l):
                 rowargs['expected'] += hpp4l[chan]['expected']
                 rowargs['pp'] += hpp4l[chan]['pp']
                 if not blind: rowargs['observed'] += int(hpp3l[chan]['observed'])
                 first = False
                 include = True
             rowargs['expected'] = latex_float(rowargs['expected'])
             rowargs['pp'] = latex_float(rowargs['pp'])
             if include: rowstring += YieldsTex.row.format(**rowargs)
         
         else:
             for cat in cats:
                 rowargs = {'mass':mass if first else '', 'reco':catTex['Hpp3l'][cat], 'expected':0., 'ap':0., 'pp':0., 'observed': '---' if blind else 0}
                 include = False
                 for chan in sorted(hpp3l):
                     if chan not in catChans['Hpp3l'][cat]: continue
                     rowargs['expected'] += hpp3l[chan]['expected']
                     rowargs['ap'] += hpp3l[chan]['ap']
                     rowargs['pp'] += hpp3l[chan]['pp']
                     if not blind: rowargs['observed'] += int(hpp3l[chan]['observed'])
                     first = False
                     include = True
                 rowargs['expected'] = latex_float(rowargs['expected'])
                 rowargs['ap'] = latex_float(rowargs['ap'])
                 rowargs['pp'] = latex_float(rowargs['pp']) if rowargs['pp'] else '---'
                 if include: rowstring += YieldsTex.row.format(**rowargs)
                 rowargs = {'mass':mass if first else '', 'reco':catTex['Hpp4l'][cat], 'expected':0., 'ap':'---', 'pp':0., 'observed': '---' if blind else 0}
                 include = False
                 for chan in sorted(hpp4l):
                     if chan not in catChans['Hpp4l'][cat]: continue
                     rowargs['expected'] += hpp4l[chan]['expected']
                     rowargs['pp'] += hpp4l[chan]['pp']
                     if not blind: rowargs['observed'] += int(hpp3l[chan]['observed'])
                     first = False
                     include = True
                 rowargs['expected'] = latex_float(rowargs['expected'])
                 rowargs['pp'] = latex_float(rowargs['pp'])
                 if include: rowstring += YieldsTex.row.format(**rowargs)
         rowstring += ' \\hline'
    print YieldsTex.table.format(rows=rowstring,**benchmarks[bp])
    print ''
         

