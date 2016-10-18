#!/usr/bin/env python
import json
from copy import deepcopy

with open('jsons/Limits/values.json') as data_file:    
    data = json.load(data_file)

class LimitsTex:
    tablePAS = '''
\\begin{{table}}[!htp]
    \\centering
    \\topcaption{{Observed (expected) 95\\% CL limits for associated (AP) and pair production (PP) and the combined limit.}}
    \\begin{{tabular}}{{| c | c | c | c |}}
        \\hline
        Benchmark & AP [GeV] & PP [GeV] & Combined [GeV] \\\\ \\hline
{rows}
        \\hline
    \\end{{tabular}}
    \\label{{tab:limits}}
\\end{{table}}
'''


    table = '''
\\begin{{table}}[!htp]
    \\centering
    \\topcaption{{Expected 95\\% CL limits for associated (AP) and pair production (PP) and the combined limit.}}
    \\begin{{tabular}}{{| c | c | c | c |}}
        \\hline
        Benchmark & AP [GeV] & PP [GeV] & Combined [GeV] \\\\ \\hline
{rows}
        \\hline
    \\end{{tabular}}
    \\label{{tab:exp_limits}}
\\end{{table}}
'''

    table2 = '''
\\begin{{table}}[!htp]
    \\centering
    \\topcaption{{Expected 95\\% CL limits for three and four lepton final states.}}
    \\begin{{tabular}}{{| c | c | c | c |}}
        \\hline
        Benchmark & $3\\ell$ AP [GeV] & $3\\ell$ PP [GeV] & $4\\ell$ PP [GeV] \\\\ \\hline
{rows}
        \\hline
    \\end{{tabular}}
    \\label{{tab:exp_limits2}}
\\end{{table}}
'''

    tableObs = '''
\\begin{{table}}[!htp]
    \\centering
    \\topcaption{{Observed 95\\% CL limits for associated (AP) and pair production (PP) and the combined limit.}}
    \\begin{{tabular}}{{| c | c | c | c |}}
        \\hline
        Benchmark & AP [GeV] & PP [GeV] & Combined [GeV] \\\\ \\hline
{rows}
        \\hline
    \\end{{tabular}}
    \\label{{tab:obs_limits}}
\\end{{table}}
'''

    tableObs2 = '''
\\begin{{table}}[!htp]
    \\centering
    \\topcaption{{Observed 95\\% CL limits for three and four lepton final states.}}
    \\begin{{tabular}}{{| c | c | c | c |}}
        \\hline
        Benchmark & $3\\ell$ AP [GeV] & $3\\ell$ PP [GeV] & $4\\ell$ PP [GeV] \\\\ \\hline
{rows}
        \\hline
    \\end{{tabular}}
    \\label{{tab:obs_limits2}}
\\end{{table}}
'''


    rowPAS = '''
        {tex:40} & {HppAP} ({HppAPExp}) & {HppPP} ({HppPPExp}) & {HppComb} ({HppCombExp}) \\\\'''

    row = '''
        {tex:40} & {HppAP} & {HppPP} & {HppComb} \\\\'''

    row2 = '''
        {tex:40} & {Hpp3lAP} & {Hpp3lPP} & {Hpp4l} \\\\'''

benchmarks = {
    'ee100': {'tex': '100\\% $\\Hpmpm \\rightarrow \\Pe\\Pe$',},
    'em100': {'tex': '100\\% $\\Hpmpm \\rightarrow \\Pe\\Pgm$',},
    'mm100': {'tex': '100\\% $\\Hpmpm \\rightarrow \\Pgm\\Pgm$',},
    'et100': {'tex': '100\\% $\\Hpmpm \\rightarrow \\Pe\\Pgt$',},
    'mt100': {'tex': '100\\% $\\Hpmpm \\rightarrow \\Pgm\\Pgt$',},
    'tt100': {'tex': '100\\% $\\Hpmpm \\rightarrow \\Pgt\\Pgt$',},
    'BP1'  : {'tex': 'Benchmark 1',},
    'BP2'  : {'tex': 'Benchmark 2',},
    'BP3'  : {'tex': 'Benchmark 3',},
    'BP4'  : {'tex': 'Benchmark 4',},

}
benchmarksObs = deepcopy(benchmarks)
benchmarksPAS = deepcopy(benchmarks)

bps = ['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4']


rowstringPAS = ''
rowstring = ''
rowstring2 = ''
rowstringObs = ''
rowstringObs2 = ''
for bp in bps:
    for lim in data[bp]:
        benchmarks[bp][lim] = data[bp][lim][2]
        benchmarksObs[bp][lim] = data[bp][lim][5]
        benchmarksPAS[bp][lim] = data[bp][lim][5]
        benchmarksPAS[bp][lim+'Exp'] = data[bp][lim][2]
    rowstring += LimitsTex.row.format(**benchmarks[bp])
    rowstring2 += LimitsTex.row2.format(**benchmarks[bp])
    rowstringObs += LimitsTex.row.format(**benchmarksObs[bp])
    rowstringObs2 += LimitsTex.row2.format(**benchmarksObs[bp])
    rowstringPAS += LimitsTex.rowPAS.format(**benchmarksPAS[bp])
    if bp=='tt100':
        rowstring += ' \\hline \\hline'
        rowstring2 += ' \\hline \\hline'
        rowstringObs += ' \\hline \\hline'
        rowstringObs2 += ' \\hline \\hline'
        rowstringPAS += ' \\hline \\hline'

print LimitsTex.table2.format(rows=rowstring2)
print ''
print LimitsTex.table.format(rows=rowstring)
print ''
print LimitsTex.tableObs2.format(rows=rowstringObs2)
print ''
print LimitsTex.tableObs.format(rows=rowstringObs)
print ''
print LimitsTex.tablePAS.format(rows=rowstringPAS)
