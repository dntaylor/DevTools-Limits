import json

blind = True

with open('jsons/Limits/values.json') as data_file:    
    data = json.load(data_file)

class LimitsTex:
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

bps = ['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4']


rowstring = ''
rowstring2 = ''
for bp in bps:
    for lim in data[bp]:
        benchmarks[bp][lim] = data[bp][lim][2 if blind else 5]
    rowstring += LimitsTex.row.format(**benchmarks[bp])
    rowstring2 += LimitsTex.row2.format(**benchmarks[bp])
    if bp=='tt100':
        rowstring += ' \\hline \\hline'
        rowstring2 += ' \\hline \\hline'

print LimitsTex.table.format(rows=rowstring)
print ''
print LimitsTex.table2.format(rows=rowstring2)
