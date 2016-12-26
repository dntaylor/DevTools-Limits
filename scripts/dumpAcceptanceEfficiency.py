#!/usr/bin/env python
import os
import sys
import json
import pickle

from DevTools.Utilities.utilities import readResults
from DevTools.Plotter.higgsUtilities import *

class YieldsTex:
    table = '''
\\begin{{landscape}}
\\begin{{table}}[!htp]
    \\centering
    \\tiny
    \\caption{{Acceptance $\\times$ efficiency for the analysis benchmarks for a mass hypothesis of {mass}~\\GeV in units of percent.}}
    \\scalebox{{.95}}{{
    \\begin{{tabular}}{{|l|rrrrrr|rrrrrr|r|}}

\\hline
Generator         & \\multicolumn{{6}}{{c|}}{{$3\\ell$ reconstructed channel}} & \\multicolumn{{6}}{{c|}}{{$4\\ell$ reconstructed channel}} &       \\\\
channel           & \\lll & \\llt & \\ltl & \\ltt & \\ttl & \\ttt            & \\llll & \\lllt & \\lltt & \\ltlt & \\lttt & \\tttt      & Total \\\\ \\hline

{rows}

    \\end{{tabular}}
    }}
    \\label{{tab:acceff{mass}}}
\\end{{table}}
\\end{{landscape}}
'''

    row = '''
{gen:20}       & {lll:6} & {llt:6} & {ltl:6} & {ltt:6} & {ttl:6} & {ttt:6}            & {llll:6}  & {lllt:6}  & {lltt:6}  & {ltlt:6}  & {lttt:6}  & {tttt:6}       & {total:6}  \\\\'''


flavorMap = {
    'e': '\\Pe',
    'm': '\\Pgm',
    't': '\\Pgt',
}

def main():
    analysis = 'tt100'

    masses = [200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]

    genRecoMap = {}
    genRecoMap['Hpp3l'] = getGenRecoChannelMap('Hpp3l')
    genRecoMap['Hpp4l'] = getGenRecoChannelMap('Hpp4l')

    gen3l = sorted(set([x for x in genRecoMap['Hpp3l'] if len(x)==3]))
    gen4l = sorted(set([x for x in genRecoMap['Hpp3l'] if len(x)==4]))

    for mass in masses:
        # get denom
        counts3l = readResults('Hpp3l','genCounts.{mass}'.format(mass=mass))
        counts4l = readResults('Hpp4l','genCounts.{mass}'.format(mass=mass))
        allCounts = {}
        allCounts3l = []
        allCounts4l = []
        # get num
        countkey3l = 'new/allMassWindow/{mass}/hpp2/{reco}/gen_{gen}'
        countkey4l = 'new/allMassWindow/{mass}/hpp2hmm2/{reco}/gen_{gen}'
        countsAP3l = readResults('Hpp3l','skims/HPlusPlusHMinusHTo3L_M-{mass}_TuneCUETP8M1_13TeV_calchep-pythia8'.format(mass=mass))
        countsPP3l = readResults('Hpp3l','skims/HPlusPlusHMinusMinusHTo4L_M-{mass}_TuneCUETP8M1_13TeV_pythia8'.format(mass=mass))
        countsPP4l = readResults('Hpp4l','skims/HPlusPlusHMinusMinusHTo4L_M-{mass}_TuneCUETP8M1_13TeV_pythia8'.format(mass=mass))
        sigCounts = {}
        for gen in gen3l+gen4l:
            genChan = gen
            h1, h2 = gen[:2], gen[2:]
            if h1>h2 and len(gen)==4: genChan = h2+h1
            if genChan not in allCounts:
                allCounts[genChan] = 0
                if len(genChan)==3: allCounts3l += [genChan]
                if len(genChan)==4: allCounts4l += [genChan]
            allCounts[genChan] += counts3l['decay'][gen] if len(gen)==3 else counts4l['decay'][gen]
            if genChan not in sigCounts: sigCounts[genChan] = {}
            for reco in genRecoMap['Hpp3l'][gen]:
                val = 0
                key = countkey3l.format(mass=mass,reco=reco,gen=gen)
                if len(gen)==3 and key in countsAP3l: val = countsAP3l[key]['count']
                if len(gen)==4 and key in countsPP3l: val = countsPP3l[key]['count']
                recoChan = reco.replace('e','l').replace('m','l')
                h1, h2 = recoChan[:2], recoChan[2:]
                if len(reco)==4 and h1>h2: recoChan = h2+h1
                if recoChan not in sigCounts[genChan]: sigCounts[genChan][recoChan] = 0
                sigCounts[genChan][recoChan] += val
            for reco in genRecoMap['Hpp4l'][gen]:
                val = 0
                key = countkey4l.format(mass=mass,reco=reco,gen=gen)
                if len(gen)==4 and key in countsPP4l: val = countsPP4l[key]['count']
                recoChan = reco.replace('e','l').replace('m','l')
                h1, h2 = recoChan[:2], recoChan[2:]
                if len(reco)==4 and h1>h2: recoChan = h2+h1
                if recoChan not in sigCounts[genChan]: sigCounts[genChan][recoChan] = 0
                sigCounts[genChan][recoChan] += val
                
        allChans = ['lll','llt','ltl','ltt','ttl','ttt','llll','lllt','lltt','ltlt','lttt','tttt']
        rows = ''
        for gen in allCounts3l+allCounts4l:
            rowargs = {}
            rowargs['gen'] = '$'+''.join([flavorMap[x] for x in gen])+'$'
            for ac in allChans:
                rowargs[ac] = '---'
            tot = 0
            for c in sigCounts[gen]:
                rowargs[c] = '${0:.1f}$'.format(float(sigCounts[gen][c])/allCounts[gen]*100)
                tot += sigCounts[gen][c]
            rowargs['total'] = '${0:.1f}$'.format(float(tot)/allCounts[gen]*100)
            rows += YieldsTex.row.format(**rowargs)
            if gen in [allCounts3l[-1], allCounts4l[-1]]: rows += ' \\hline'

        print YieldsTex.table.format(rows=rows,mass=mass)
        print ''


if __name__ == "__main__":
    status = main()
    sys.exit(status)

