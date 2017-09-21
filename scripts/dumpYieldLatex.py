#!/usr/bin/env python
import os
import sys
import argparse
import json
from copy import deepcopy

def getYields(analysis,mode,mass):
    with open('jsons/{0}/{1}/{2}.json'.format(analysis,mode,mass)) as data_file:    
        data = json.load(data_file)
    return data

class YieldsTex:
    table='''
\\subsection{{{tex} yields}}
{{\\tiny
\\begin{{longtable}}{{ll|r|rrr|r}}
\\hline
$\\mHpmpm$ [\\GeV] & Channel & $N_{{exp}}$ & $N_{{AP}}$ & $N_{{PP}}$ & $\\rm{{PP_{{R}}}}$ & $N_{{obs}}$ \\\\
\\hline
\\endhead
{rows}
\\end{{longtable}}
}}
'''


    row = '''
{mass:4} & ${reco:18}$ & {expected:25} & {ap:25} & {pp:25} & {ppR:25} & {observed:4} \\\\'''


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

def printYields(args):
    print '\\section{Yields\\label{sec:yields}}'
    print ''
    
    
    for bp in bps:
        rowstring = ''
        for mass in masses:
             hpp3l = getYields('Hpp3l',bp,mass)
             hpp4l = getYields('Hpp4l',bp,mass)
             first = True
             if args.verbose==2:
                 for chan in sorted(hpp3l):
                     if '_SB' in chan: continue
                     rowargs = {'mass':mass if first else '', 'reco':''.join([flavorMap[lep] for lep in chan]), 'expected':0., 'ap':0., 'pp':0., 'ppR':0., 'observed': '---' if not args.unblind else 0}
                     include = False
                     rowargs['expected'] += hpp3l[chan]['expected']
                     rowargs['ap'] += hpp3l[chan]['ap']
                     rowargs['pp'] += hpp3l[chan]['pp']
                     rowargs['ppR'] += hpp3l[chan]['ppR']
                     if args.unblind: rowargs['observed'] += int(hpp3l[chan]['observed'])
                     first = False
                     include = True
                     rowargs['expected'] = latex_float(rowargs['expected'])
                     rowargs['ap'] = latex_float(rowargs['ap'])
                     rowargs['pp'] = latex_float(rowargs['pp']) if rowargs['pp'] else '---'
                     rowargs['ppR'] = latex_float(rowargs['ppR']) if rowargs['ppR'] else '---'
                     if include: rowstring += YieldsTex.row.format(**rowargs)
                 for chan in sorted(hpp4l):
                     if '_SB' in chan: continue
                     rowargs = {'mass':mass if first else '', 'reco':''.join([flavorMap[lep] for lep in chan]), 'expected':0., 'ap':'---', 'pp':0., 'ppR':0., 'observed': '---' if not args.unblind else 0}
                     include = False
                     rowargs['expected'] += hpp4l[chan]['expected']
                     rowargs['pp'] += hpp4l[chan]['pp']
                     rowargs['ppR'] += hpp4l[chan]['ppR']
                     if args.unblind: rowargs['observed'] += int(hpp4l[chan]['observed'])
                     first = False
                     include = True
                     rowargs['expected'] = latex_float(rowargs['expected'])
                     rowargs['pp'] = latex_float(rowargs['pp'])
                     rowargs['ppR'] = latex_float(rowargs['ppR'])
                     if include: rowstring += YieldsTex.row.format(**rowargs)
             elif args.verbose==1:
                 for cat in cats:
                     rowargs = {'mass':mass if first else '', 'reco':catTex['Hpp3l'][cat], 'expected':0., 'ap':0., 'pp':0., 'ppR':0., 'observed': '---' if not args.unblind else 0}
                     include = False
                     for chan in sorted(hpp3l):
                         if chan not in catChans['Hpp3l'][cat]: continue
                         rowargs['expected'] += hpp3l[chan]['expected']
                         rowargs['ap'] += hpp3l[chan]['ap']
                         rowargs['pp'] += hpp3l[chan]['pp']
                         rowargs['ppR'] += hpp3l[chan]['ppR']
                         if args.unblind: rowargs['observed'] += int(hpp3l[chan]['observed'])
                         first = False
                         include = True
                     rowargs['expected'] = latex_float(rowargs['expected'])
                     rowargs['ap'] = latex_float(rowargs['ap'])
                     rowargs['pp'] = latex_float(rowargs['pp']) if rowargs['pp'] else '---'
                     rowargs['ppR'] = latex_float(rowargs['ppR']) if rowargs['ppR'] else '---'
                     if include: rowstring += YieldsTex.row.format(**rowargs)
                 for cat in cats:
                     rowargs = {'mass':mass if first else '', 'reco':catTex['Hpp4l'][cat], 'expected':0., 'ap':'---', 'pp':0., 'ppR':0., 'observed': '---' if not args.unblind else 0}
                     include = False
                     for chan in sorted(hpp4l):
                         if chan not in catChans['Hpp4l'][cat]: continue
                         rowargs['expected'] += hpp4l[chan]['expected']
                         rowargs['pp'] += hpp4l[chan]['pp']
                         rowargs['ppR'] += hpp4l[chan]['ppR']
                         if args.unblind: rowargs['observed'] += int(hpp4l[chan]['observed'])
                         first = False
                         include = True
                     rowargs['expected'] = latex_float(rowargs['expected'])
                     rowargs['pp'] = latex_float(rowargs['pp'])
                     rowargs['ppR'] = latex_float(rowargs['ppR'])
                     if include: rowstring += YieldsTex.row.format(**rowargs)
             else:
                 rowargs = {'mass':mass if first else '', 'reco':'3\\ell', 'expected':0., 'ap':0., 'pp':0., 'ppR':0., 'observed': '---' if not args.unblind else 0}
                 include = False
                 for chan in sorted(hpp3l):
                     rowargs['expected'] += hpp3l[chan]['expected']
                     rowargs['ap'] += hpp3l[chan]['ap']
                     rowargs['pp'] += hpp3l[chan]['pp']
                     rowargs['ppR'] += hpp3l[chan]['ppR']
                     if args.unblind: rowargs['observed'] += int(hpp3l[chan]['observed'])
                     first = False
                     include = True
                 rowargs['expected'] = latex_float(rowargs['expected'])
                 rowargs['ap'] = latex_float(rowargs['ap'])
                 rowargs['pp'] = latex_float(rowargs['pp']) if rowargs['pp'] else '---'
                 rowargs['ppR'] = latex_float(rowargs['ppR']) if rowargs['ppR'] else '---'
                 if include: rowstring += YieldsTex.row.format(**rowargs)
                 rowargs = {'mass':mass if first else '', 'reco':'4\\ell', 'expected':0., 'ap':'---', 'pp':0., 'ppR':0., 'observed': '---' if not args.unblind else 0}
                 include = False
                 for chan in sorted(hpp4l):
                     rowargs['expected'] += hpp4l[chan]['expected']
                     rowargs['pp'] += hpp4l[chan]['pp']
                     rowargs['ppR'] += hpp4l[chan]['ppR']
                     if args.unblind: rowargs['observed'] += int(hpp4l[chan]['observed'])
                     first = False
                     include = True
                 rowargs['expected'] = latex_float(rowargs['expected'])
                 rowargs['pp'] = latex_float(rowargs['pp'])
                 rowargs['ppR'] = latex_float(rowargs['ppR'])
                 if include: rowstring += YieldsTex.row.format(**rowargs)
             
             rowstring += ' \\hline'
        print YieldsTex.table.format(rows=rowstring,**benchmarks[bp])
        print ''
         

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Dump latex yield tables')

    parser.add_argument('--unblind', action='store_true', help='Unblind results')
    parser.add_argument('--verbose', type=int, default=0,
        help='Verbosity level: 0 = sum all, 1 = per category, 2 = per channel.'
    )


    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    printYields(args)

if __name__ == "__main__":
    status = main()
    sys.exit(status)

