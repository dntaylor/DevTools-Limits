#!/usr/bin/env python
import os
import sys
import json
import pickle

from DevTools.Utilities.utilities import readResults

class YieldsTex:
    table = '''
\\begin{{table}}[!htp]
    \\centering
    \\topcaption{{Yields for the signal region for a mass hypothesis of {mass}~\\GeV.
                The uncertainties shown are combined statistical and systematic.}}
    \\begin{{tabular}}{{|ll|rrr|r|}}

\\hline
                                         &         & \\multicolumn{{4}}{{c|}}{{Signal Region}}                   \\\\
Benchmark                                & Channel & Assoc. Prod.  & Pair Prod.    & Background    & Obs. \\\\ \\hline

{rows}

\\hline
    \\end{{tabular}}
    \\label{{tab:yields{mass}}}
\\end{{table}}
'''

    table2 = '''
\\begin{{table}}[!htp]
    \\centering
    \\topcaption{{Yields for the sideband and signal region for a mass hypothesis of {mass}~\\GeV.
                The uncertainties shown are combined statistical and systematic.}}
    \\begin{{tabular}}{{|ll|rrr|r|rrr|r|}}

\\hline
                                         &         & \\multicolumn{{4}}{{c|}}{{Sideband}}                        & \\multicolumn{{4}}{{c|}}{{Signal Region}}                   \\\\
Benchmark                                & Channel & Assoc. Prod.  & Pair Prod.    & Background    & Obs. & Assoc. Prod.  & Pair Prod.    & Background    & Obs. \\\\ \\hline

{rows}

\\hline
    \\end{{tabular}}
    \\label{{tab:yields{mass}}}
\\end{{table}}
'''

    row = '''
{benchmark:40} & {channel:6} & {apSR:40} & {ppSR:40} & {bgSR:40} & ${obsSR:3}$ \\\\'''

    row2 = '''
{benchmark:40} & {channel:6} & {apSB:25} & {ppSB:25} & {bgSB:25} & ${obsSB:3}$ & {apSR:25} & {ppSR:25} & {bgSR:25} & ${obsSR:3}$ \\\\'''


def readResultsTemp(name):
    pfile = '{0}.pkl'.format(name)
    with open(pfile,'rb') as f:
        results = pickle.load(f)
    return results


def expandChannels(channels):
    hpps = {
        'll' : ['ee','em','mm'],
        'el' : ['ee','em'],
        'ml' : ['em','mm'],
        'ee' : ['ee'],
        'em' : ['em'],
        'mm' : ['mm'],
        'lt' : ['et','mt'],
        'et' : ['et'],
        'mt' : ['mt'],
        'tt' : ['tt']
    }
    hms = {
        'l' : ['e','m'],
        'e' : ['e'],
        'm' : ['m'],
        't' : ['t'],
    }
    allChannels = []
    for chan in channels:
        if len(chan) == 3:
            hpp = chan[:2]
            hm = chan[2:]
            for i in hpps[hpp]:
                for j in hms[hm]:
                    allChannels += [i+j]
        if len(chan) == 4:
            hpp = chan[:2]
            hmm = chan[2:]
            for i in hpps[hpp]:
                for j in hpps[hmm]:
                    allChannels += [i+j]
    return allChannels

bps = ['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4']
masses = [200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]
flavorMap = {
    'e': '\\Pe',
    'm': '\\Pgm',
    't': '\\Pgt',
}

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
}



analyses = ['Hpp3lAP','Hpp3lPP','Hpp4l','HppAP','HppPP','HppComb']
analyses = ['HppComb']
modes = ['ee100','em100','et100','mm100','mt100','tt100']
#modes = ['ee100','em100','mm100']
masses = [200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]
#masses = [500]

channels = {
    'ee100' : {'lll': ['eel'], 'llt': ['eet'], 'llll': ['eeee']},
    'em100' : {'lll': ['eml'], 'llt': ['emt'], 'llll': ['emem']},
    'mm100' : {'lll': ['mml'], 'llt': ['mmt'], 'llll': ['mmmm']},
    'et100' : {'lll': ['ell'], 'llt': ['elt'], 'ltl': ['etl'], 'ltt': ['ett'], 'llll': ['elel'], 'lllt': ['elet'], 'ltlt': ['etet']},
    'mt100' : {'lll': ['mll'], 'llt': ['mlt'], 'ltl': ['mtl'], 'ltt': ['mtt'], 'llll': ['mlml'], 'lllt': ['mlmt'], 'ltlt': ['mtmt']},
    'tt100' : {'lll': ['lll'], 'llt': ['llt'], 'ltl': ['ltl'], 'ltt': ['ltt'], 'ttl': ['ttl'], 'ttt': ['ttt'],  'llll': ['llll'], 'lllt': ['lllt'], 'lltt': ['lltt'], 'ltlt': ['ltlt'], 'lttt': ['lttt'], 'tttt': ['tttt']},
}


def latexfloat(val,errUp,errDown):
    f = '{0:.2e}'.format(val)
    f = f.split('e')
    if len(f)==1: # no exponent
        return '${0:.2f}^{{+{1:.2f}}}_{{-{2:.2f}}}$'.format(val,errUp,errDown)
    else:
        v, p = f
        p = int(p)
        return '${0:.2f}^{{+{1:.2f}}}_{{-{2:.2f}}} \\times 10^{{{3}}}$'.format(float(v),errUp/(10**p),errDown/(10**p),p)

def main():
    data = readResultsTemp('limit_uncertainties')
    for analysis in analyses:
        for mass in masses:
            rows = ''
            for mode in modes:
                inputs3l = readResults('Hpp3l','{0}/{1}'.format(mode,mass))
                inputs4l = readResults('Hpp4l','{0}/{1}'.format(mode,mass))
                chans3l = [c for c in sorted(channels[mode]) if len(c)==3]
                chans4l = [c for c in sorted(channels[mode]) if len(c)==4]
                for chan in chans3l + chans4l:
                    recoChans = expandChannels(channels[mode][chan])
                    inputs = inputs3l if len(chan)==3 else inputs4l
                    obsSR = 0
                    obsSB = 0
                    #print mode, chan, recoChans
                    for reco in recoChans:
                        obsSR += inputs[reco]['observed']
                        obsSB += inputs[reco+'_SB']['observed']
                    args = {}
                    args['benchmark'] = benchmarks[mode]['tex'] if chan==chans3l[0] else ''
                    args['obsSB'] = int(obsSB)
                    args['obsSR'] = int(obsSR)
                    args['channel'] = '\\{0}'.format(chan)
                    apSB         = data[analysis][mode][mass][chan]['apSB']['val']
                    apErrUpSB    = data[analysis][mode][mass][chan]['apSB']['errUp']
                    apErrDownSB  = data[analysis][mode][mass][chan]['apSB']['errDown']
                    ppSB         = data[analysis][mode][mass][chan]['ppSB']['val']
                    ppErrUpSB    = data[analysis][mode][mass][chan]['ppSB']['errUp']
                    ppErrDownSB  = data[analysis][mode][mass][chan]['ppSB']['errDown']
                    bgSB         = data[analysis][mode][mass][chan]['bgSB']['val']
                    bgErrUpSB    = data[analysis][mode][mass][chan]['bgSB']['errUp']
                    bgErrDownSB  = data[analysis][mode][mass][chan]['bgSB']['errDown']
                    apSR         = data[analysis][mode][mass][chan]['apSR']['val']
                    apErrUpSR    = data[analysis][mode][mass][chan]['apSR']['errUp']
                    apErrDownSR  = data[analysis][mode][mass][chan]['apSR']['errDown']
                    ppSR         = data[analysis][mode][mass][chan]['ppSR']['val']
                    ppErrUpSR    = data[analysis][mode][mass][chan]['ppSR']['errUp']
                    ppErrDownSR  = data[analysis][mode][mass][chan]['ppSR']['errDown']
                    bgSR         = data[analysis][mode][mass][chan]['bgSR']['val']
                    bgErrUpSR    = data[analysis][mode][mass][chan]['bgSR']['errUp']
                    bgErrDownSR  = data[analysis][mode][mass][chan]['bgSR']['errDown']
                    args['apSB'] = latexfloat(apSB,apErrUpSB,apErrDownSB) if len(chan)==3 else '---'
                    args['ppSB'] = latexfloat(ppSB,ppErrUpSB,ppErrDownSB) if ppSB>1e-9 else '---'
                    args['bgSB'] = latexfloat(bgSB,bgErrUpSB,bgErrDownSB) if bgSB>1e-9 else '$0.00$'
                    args['apSR'] = latexfloat(apSR,apErrUpSR,apErrDownSR) if len(chan)==3 else '---'
                    args['ppSR'] = latexfloat(ppSR,ppErrUpSR,ppErrDownSR) if ppSR>1e-9 else '---'
                    args['bgSR'] = latexfloat(bgSR,bgErrUpSR,bgErrDownSR) if bgSR>1e-9 else '$0.00$'
                    rows += YieldsTex.row.format(**args)
                    if chan==chans4l[-1]: rows += ' \\hline'
                    #print analysis, mode, mass, chan
                    #if len(chan)==3:
                    #    print obsSR, bgSR, ppSR, apSR
                    #    print obsSB, bgSB, ppSB, apSB
                    #else:
                    #    print obsSR, bgSR, ppSR
                    #    print obsSB, bgSB, ppSB
    
            print YieldsTex.table.format(rows=rows,mass=mass)
            print ''


if __name__ == "__main__":
    status = main()
    sys.exit(status)
