import argparse
parser = argparse.ArgumentParser(description='Produce source for spline interpolation.')
parser.add_argument('template_file',  type=str, help='Name of the template file')
parser.add_argument('filename_out', type=str, help='Name of the generated file')
parser.add_argument('-o','--order', default=4, type=int, choices=range(6), help='Maximum interpolation order')

args = parser.parse_args()

filename_in = args.template_file
filename_out = args.filename_out

max_order = args.order

def print_expr(symbs, inds=[]):
    if len(symbs) == 0:
        return 'C[{}]'.format(str.join(',',['i{}+{}'.format(i,k) for i,k in enumerate(inds)]))
    else:
        h = symbs[0]
        q = symbs[1:]
        exprs = [  '{}_{}*({})'.format(h,i,print_expr(q,inds + [i])) for i in range(4)]
        return str.join( ' + ', exprs )


values = []
dvalues = []
for order in range(max_order+1):
    expr = print_expr( ['Phi{}'.format(i) for i in range(order)] )
    values.append( expr )
    dv = []
    for i in range(order):
        args =  ['Phi{}'.format(h) for h in range(order)]
        args[i] = 'dPhi{}'.format(i)
        dexpr = print_expr( args )
        dv.append(dexpr)
    dvalues.append(dv)


import tempita

with file(filename_in) as f:
    txt = f.read()

s = tempita.sub(txt,values=values,dvalues=dvalues,max_order=max_order)

with file(filename_out,'w') as f:
    f.write(s)
