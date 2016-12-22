#!/usr/bin/env python
import os
import sys
import argparse
import json
from copy import deepcopy

def getYields(analysis):
    dataSamples = ['DoubleMuon','SingleMuon','DoubleEG','SingleElectron','MuonEG','Tau']
    data = {}
    for ds in dataSamples:
        with open('jsons/{0}/skims/{1}.json'.format(analysis,ds)) as data_file:    
            allvals = json.load(data_file)
            for key in [k for k in allvals.keys() if k.startswith('default/')]:
                d,c = key.split('/')
                if c not in data: data[c] = 0
                data[c] += allvals[key]['val']
    return data

class YieldsTex:
    table='''
\\begin{{table}}[!htp]
    \\centering
    \\begin{{tabular}}{{|ll|rrr|r|}}

\\hline
Channel & Preselection Yield \\\\ \\hline

{rows}

    \\end{{tabular}}
    \\caption{{Yields after the preselection requirements for each channel.
    \\label{{tab:preyields}}
    }}
\\end{{table}}
'''


    row = '''
{chan:10} & {yield:10} \\\\'''


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
    
    rowstring = ''
    hpp3l = getYields('Hpp3l')
    hpp4l = getYields('Hpp4l')
    if args.verbose==2:
        for chan in sorted(hpp3l):
            rowargs = {'chan':''.join([flavorMap[lep] for lep in chan]), 'yield':int(hpp3l[chan])}
            rowstring += YieldsTex.row.format(**rowargs)
        rowstring += ' \\hline'
        for chan in sorted(hpp4l):
            rowargs = {'chan':''.join([flavorMap[lep] for lep in chan]), 'yield':int(hpp4l[chan])}
            rowstring += YieldsTex.row.format(**rowargs)
    elif args.verbose==1:
        for cat in cats:
            rowargs = {'chan':catTex['Hpp3l'][cat], 'yield':0}
            for chan in sorted(hpp3l):
                if chan not in catChans['Hpp3l'][cat]: continue
                rowargs['yield'] += int(hpp3l[chan])
            rowstring += YieldsTex.row.format(**rowargs)
        rowstring += ' \\hline'
        for cat in cats:
            rowargs = {'chan':catTex['Hpp4l'][cat], 'yield':0}
            for chan in sorted(hpp4l):
                if chan not in catChans['Hpp4l'][cat]: continue
                rowargs['yield'] += int(hpp4l[chan])
            rowstring += YieldsTex.row.format(**rowargs)

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

