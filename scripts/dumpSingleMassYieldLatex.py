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
{{\\tiny
\\begin{{longtable}}{{ll|rrr|r|r}}
\\hline
Benchmark & Channel & $N_{{AP}}$ & $N_{{PP}}$ & $\\rm{{PP_{{R}}}}$ & $N_{{exp}}$ & $N_{{obs}}$ \\\\
\\hline
\\endhead
{rows}
\\end{{longtable}}
}}
'''


    row = '''
{bp:40} & ${reco:18}$ & {ap:25} & {pp:25} & {ppR:25} & {expected:25} & {observed:4} \\\\'''


benchmarks = {
    'ee100': {'tex': '100\\% $\\Phi \\rightarrow \\Pe\\Pe$',},
    'em100': {'tex': '100\\% $\\Phi \\rightarrow \\Pe\\Pgm$',},
    'mm100': {'tex': '100\\% $\\Phi \\rightarrow \\Pgm\\Pgm$',},
    'et100': {'tex': '100\\% $\\Phi \\rightarrow \\Pe\\Pgt$',},
    'mt100': {'tex': '100\\% $\\Phi \\rightarrow \\Pgm\\Pgt$',},
    'tt100': {'tex': '100\\% $\\Phi \\rightarrow \\Pgt\\Pgt$',},
    'BP1'  : {'tex': 'Benchmark 1',},
    'BP2'  : {'tex': 'Benchmark 2',},
    'BP3'  : {'tex': 'Benchmark 3',},
    'BP4'  : {'tex': 'Benchmark 4',},
    #'ee100': {'tex': '100\\% electron-electron',},
    #'em100': {'tex': '100\\% electron-muon',},
    #'mm100': {'tex': '100\\% muon-muon',},
    #'et100': {'tex': '100\\% electron-tau',},
    #'mt100': {'tex': '100\\% muon-tau',},
    #'tt100': {'tex': '100\\% tau-tau',},
    #'BP1'  : {'tex': 'Benchmark 1',},
    #'BP2'  : {'tex': 'Benchmark 2',},
    #'BP3'  : {'tex': 'Benchmark 3',},
    #'BP4'  : {'tex': 'Benchmark 4',},

}

bps = ['ee100','em100','mm100','et100','mt100','tt100']#,'BP1','BP2','BP3','BP4']
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

def latex_float(f,e=0):
    f = '{0:.3g}'.format(f)
    f = f.split('e')
    if len(f)==1: 
        if e:
            return '${0} \\pm {1:.3f}$'.format(f[0],e)
        else:
            return '${0}$'.format(*f)
    p = int(f[1])
    if e:
        return '${0} \\pm {1:.3f} \\times 10^{{{2}}}$'.format(f[0],e/(10**p),p)
    else:
        return '${0} \\times 10^{{{1}}}$'.format(f[0],p)

def add(args, rowargs, chan):
    rowargs['expected'] += chan['expected']
    rowargs['expectedErr2'] += (chan['expected']*chan['alphaError']/chan['alpha'])**2 if chan['expected'] else 0.
    if 'ap' in chan: rowargs['ap'] += chan['ap']
    if 'ap' in chan: rowargs['apErr2'] += chan['apError']**2
    rowargs['pp'] += chan['pp']
    rowargs['ppErr2'] += chan['ppError']**2
    rowargs['ppR'] += chan['ppR']
    rowargs['ppRErr2'] += chan['ppRError']**2
    if args.unblind: rowargs['observed'] += int(chan['observed'])
    return rowargs

def fix(rowargs):
    rowargs['expected'] = latex_float(rowargs['expected'],rowargs['expectedErr2']**0.5)
    rowargs['ap'] = latex_float(rowargs['ap'],rowargs['apErr2']**0.5) if rowargs['ap'] else '---'
    rowargs['pp'] = latex_float(rowargs['pp'],rowargs['ppErr2']**0.5) if rowargs['pp'] else '---'
    rowargs['ppR'] = latex_float(rowargs['ppR'],rowargs['ppRErr2']**0.5) if rowargs['ppR'] else '---'
    return rowargs
    

def printYields(args):
    
    mass = 500
    
    rowstring = ''
    for bp in bps:
        hpp3l = getYields('Hpp3l',bp,mass)
        hpp4l = getYields('Hpp4l',bp,mass)
        first = True
        if args.verbose==2:
            for chan in sorted(hpp3l):
                if '_SB' in chan: continue
                rowargs = {'bp':benchmarks[bp]['tex'] if first else '', 'reco':''.join([flavorMap[lep] for lep in chan]), 'expected':0., 'ap':0., 'pp':0., 'ppR':0., 'expectedErr2':0., 'apErr2':0., 'ppErr2':0., 'ppRErr2':0.,  'observed': '---' if not args.unblind else 0}
                rowargs = add(args, rowargs, hpp3l[chan])
                first = False
                include = True
                rowargs = fix(rowargs)
                if include: rowstring += YieldsTex.row.format(**rowargs)
            for chan in sorted(hpp4l):
                if '_SB' in chan: continue
                rowargs = {'bp':benchmarks[bp]['tex'] if first else '', 'reco':''.join([flavorMap[lep] for lep in chan]), 'expected':0., 'ap':0., 'pp':0., 'ppR':0., 'expectedErr2':0., 'apErr2':0., 'ppErr2':0., 'ppRErr2':0., 'observed': '---' if not args.unblind else 0}
                rowargs = add(args, rowargs, hpp4l[chan])
                first = False
                include = True
                rowargs = fix(rowargs)
                if include: rowstring += YieldsTex.row.format(**rowargs)
        elif args.verbose==1:
            for cat in cats:
                rowargs = {'bp':benchmarks[bp]['tex'] if first else '', 'reco':catTex['Hpp3l'][cat], 'expected':0., 'ap':0., 'pp':0., 'ppR':0., 'expectedErr2':0., 'apErr2':0., 'ppErr2':0., 'ppRErr2':0., 'observed': '---' if not args.unblind else 0}
                include = False
                for chan in sorted(hpp3l):
                    if chan not in catChans['Hpp3l'][cat]: continue
                    rowargs = add(args, rowargs, hpp3l[chan])
                    first = False
                    include = True
                rowargs = fix(rowargs)
                if include: rowstring += YieldsTex.row.format(**rowargs)
            for cat in cats:
                rowargs = {'bp':benchmarks[bp]['tex'] if first else '', 'reco':catTex['Hpp4l'][cat], 'expected':0., 'ap':0., 'pp':0., 'ppR':0., 'expectedErr2':0., 'apErr2':0., 'ppErr2':0., 'ppRErr2':0., 'observed': '---' if not args.unblind else 0}
                include = False
                for chan in sorted(hpp4l):
                    if chan not in catChans['Hpp4l'][cat]: continue
                    rowargs = add(args, rowargs, hpp4l[chan])
                    first = False
                    include = True
                rowargs = fix(rowargs)
                if include: rowstring += YieldsTex.row.format(**rowargs)
        else:
            rowargs = {'bp':benchmarks[bp]['tex'] if first else '', 'reco':'3\\ell', 'expected':0., 'ap':0., 'pp':0., 'ppR':0., 'expectedErr2':0., 'apErr2':0., 'ppErr2':0., 'ppRErr2':0., 'observed': '---' if not args.unblind else 0}
            include = False
            for chan in sorted(hpp3l):
                if '_SB' in chan: continue
                rowargs = add(args, rowargs, hpp3l[chan])
                first = False
                include = True
            rowargs = fix(rowargs)
            if include: rowstring += YieldsTex.row.format(**rowargs)
            rowargs = {'bp':benchmarks[bp]['tex'] if first else '', 'reco':'4\\ell', 'expected':0., 'ap':0., 'pp':0., 'ppR':0., 'expectedErr2':0., 'apErr2':0., 'ppErr2':0., 'ppRErr2':0., 'observed': '---' if not args.unblind else 0}
            include = False
            for chan in sorted(hpp4l):
                if '_SB' in chan: continue
                rowargs = add(args, rowargs, hpp4l[chan])
                first = False
                include = True
            rowargs = fix(rowargs)
            if include: rowstring += YieldsTex.row.format(**rowargs)
        
        rowstring += ' \\hline'
    print YieldsTex.table.format(rows=rowstring)
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

