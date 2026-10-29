import pandas as pd, numpy as np
# Prints every simulation number quoted in Sections 3.4, 5.1, 5.4, the abstract and the
# Discussion, from the summaries written by simulation/02_make_tables_figures.R.
# Run from the repository root:  python3 scripts/r500_text_values.py
R='simulation/results/'
d=pd.read_csv(R+'sim_summary.csv'); cv=pd.read_csv(R+'sim_coverage_summary.csv'); cg=pd.read_csv(R+'sim_convergence_summary.csv')
S=list(d.scenario.unique()); HET=[s for s in S if s!='No heterogeneity']
def g(s,m,r,c): return d[(d.scenario==s)&(d.method==m)&(d.rule==r)][c].iloc[0]
P=['MCP','SCAD','LASSO','ENet']
print('MCC_se MCP-H range', d[d.method=='MCP-H'].MCC_se.min(), d[d.method=='MCP-H'].MCC_se.max())
for st in 'HSP':
    v=d[(d.rule=='min')&d.method.isin([f'LASSO-{st}',f'ENet-{st}'])].FDR; print('convex FDR min',st,v.min(),v.max())
v=d[(d.rule=='min')&d.method.str.startswith('MCP-')].FDR; print('MCP FDR min all structures',v.min(),v.max())
v=d[(d.rule=='1se')&d.method.str.match('(LASSO|ENet)-')].FDR; print('convex FDR 1se',v.min(),v.max())
# SCAD less sparse than MCP: FDR SCAD>MCP for each structure/scenario/rule
cnt=0;tot=0
for s in S:
  for st in 'HSP':
    for r in ['1se','min']:
      tot+=1; cnt+= g(s,f'SCAD-{st}',r,'FDR')>g(s,f'MCP-{st}',r,'FDR')
print('SCAD FDR > MCP FDR',cnt,'/',tot)
for p in P:
  for r in ['1se','min']:
    c=sum(g(s,f'{p}-S',r,'MCC')>g(s,f'{p}-P',r,'MCC') for s in S); print('strat>pool MCC',p,r,c)
for p in ['LASSO','ENet']:
  for r in ['1se','min']:
    print('pool>strat MCC',p,r,sum(g(s,f'{p}-P',r,'MCC')>g(s,f'{p}-S',r,'MCC') for s in S),
          'pool FDR<strat',sum(g(s,f'{p}-P',r,'FDR')<g(s,f'{p}-S',r,'FDR') for s in S))
for s in S: print('MCP 1se S-P gain',s, round(g(s,'MCP-S','1se','MCC')-g(s,'MCP-P','1se','MCC'),3))
print('base MCP S/P 1se',g('Base','MCP-S','1se','MCC'),g('Base','MCP-P','1se','MCC'),'strong',g('Strong heterogeneity','MCP-S','1se','MCC'),g('Strong heterogeneity','MCP-P','1se','MCC'))
# hierarchy, 1se: best MCC (ties included)
for r in ['1se','min']:
  best=[]; 
  for s in S:
    x=d[(d.scenario==s)&(d.rule==r)]; mx=x.MCC.max()
    if g(s,'MCP-H',r,'MCC')>=mx-1e-12: best.append(s)
  print('MCP-H best',r,len(best),[s for s in S if s not in best])
for s in S:
  m=g(s,'MCP-H','1se','MCC')-g(s,'MCP-S','1se','MCC'); se=np.hypot(g(s,'MCP-H','1se','MCC_se'),g(s,'MCP-S','1se','MCC_se'))
  m2=g(s,'MCP-H','min','MCC')-max(g(s,'MCP-S','min','MCC'),g(s,'MCP-P','min','MCC'))
  se2=np.hypot(g(s,'MCP-H','min','MCC_se'),g(s,'MCP-S','min','MCC_se') if g(s,'MCP-S','min','MCC')>=g(s,'MCP-P','min','MCC') else g(s,'MCP-P','min','MCC_se'))
  print(f'{s:25s} 1se H-S {m:+.3f} se {se:.4f} {"*" if abs(m)>2*se else ""} | min H-best(S,P) {m2:+.3f} se {se2:.4f} {"*" if abs(m2)>2*se2 else ""}  H {g(s,"MCP-H","min","MCC"):.3f}')
print('base 1se MCP-H',g('Base','MCP-H','1se','MCC'),g('Base','MCP-H','1se','exact'))
for s in ['Strong heterogeneity','p = 500','Dense heterogeneity','No heterogeneity']:
  print(s,'1se H S P SCAD-H', [round(g(s,m,'1se','MCC'),3) for m in ['MCP-H','MCP-S','MCP-P','SCAD-H']])
for p in ['SCAD','LASSO','ENet']:
  dif=[g(s,f'{p}-H','1se','MCC')-g(s,f'{p}-S','1se','MCC') for s in S]
  print('H>S 1se',p,sum(x>0 for x in dif),'maxabs',round(max(map(abs,dif)),3))
# MSE theta ratios S/H
cnt=0
for s in HET:
  for p in P:
    for r in ['1se','min']:
      cnt+= g(s,f'{p}-H',r,'MSE_theta')<g(s,f'{p}-S',r,'MSE_theta')
print('H<S MSE in',cnt,'of',len(HET)*8)
def rat(s,p,r): return g(s,f'{p}-S',r,'MSE_theta')/g(s,f'{p}-H',r,'MSE_theta')
other=[s for s in HET if s not in ['Base','Strong heterogeneity','Small studies, p > n_k']]
for r in ['1se','min']:
  print('MCP ratio',r,'base',round(rat('Base','MCP',r),2),'other',round(min(rat(s,'MCP',r) for s in other),2),round(max(rat(s,'MCP',r) for s in other),2),'strong',round(rat('Strong heterogeneity','MCP',r),2),'small',round(rat('Small studies, p > n_k','MCP',r),2), 'max all non-small', round(max(rat(s,'MCP',r) for s in HET if s!='Small studies, p > n_k'),2),'min nonsmall',round(min(rat(s,'MCP',r) for s in HET if s!='Small studies, p > n_k'),2))
  for p in ['LASSO','ENet']: print(' max',p,r,round(max(rat(s,p,r) for s in HET),2))
for p in ['MCP','SCAD']:
  print('nohet ratio',p,[round(rat('No heterogeneity',p,r),2) for r in ['1se','min']])
print('nohet MCC MCP H/S', [ (g('No heterogeneity','MCP-H',r,'MCC'),g('No heterogeneity','MCP-S',r,'MCC')) for r in ['1se','min']])
print('nohet n_dev MCP-H', [g('No heterogeneity','MCP-H',r,'n_dev') for r in ['1se','min']])
print('small ratios', sorted(round(rat('Small studies, p > n_k',p,r),2) for p in P for r in ['1se','min']))
# near cancel base
for m in ['MCP-H','MCP-S']:
  print('nc',m,[ (round(g('Base',m,r,'err_theta_nc'),3), round(g('Base',m,r,'err_theta_nc_se'),3)) for r in ['1se','min']])
print('n_near_cancel base',g('Base','MCP-H','1se','n_near_cancel'))
# deviations
print('base dev MCP-H', [(round(g('Base','MCP-H',r,'TPR_dev'),3),round(g('Base','MCP-H',r,'FDR_dev'),3)) for r in ['1se','min']])
for r in ['1se','min']:
  x=d[(d.method=='MCP-H')&(d.rule==r)&(d.scenario.isin(HET))]
  print('dev MCP-H',r,'TPR',x.TPR_dev.min(),x.TPR_dev.max(),'FDR',x.FDR_dev.min(),x.FDR_dev.max(), x.loc[x.FDR_dev.idxmax(),'scenario'], x.loc[x.TPR_dev.idxmax(),'scenario'])
print('strong dev', [g('Strong heterogeneity','MCP-H',r,'TPR_dev') for r in ['1se','min']])
print('base convex dev 1se', [(g('Base',m,'1se','TPR_dev'),g('Base',m,'1se','FDR_dev')) for m in ['LASSO-H','ENet-H']])
# timing
print('time base MCP-H', g('Base','MCP-H','1se','seconds'), 'MCP-S', g('Base','MCP-S','1se','seconds'))
# convergence
print('fits',cg.fits.sum(),'iter',cg.fits_with_iter_limit.sum(), cg.fits_with_iter_limit.sum()/cg.fits.sum())
print('boundary mean', ((cg.boundary_min+cg.boundary_1se)/2).mean())
