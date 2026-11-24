import numpy as np
def load(n):
    import os; d=np.genfromtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)),n+'.csv'),delimiter=',',skip_header=1); w,a=d[:,0],d[:,1].copy()
    for lo,hi in [(484,488),(630,641),(653,659)]:          # lamp/filter-change artefacts
        m=(w>=lo)&(w<=hi); k=((w>=lo-6)&(w<lo))|((w>hi)&(w<=hi+6))
        a[m]=np.polyval(np.polyfit(w[k],a[k],2),w[m])
    base=a[(w>=700)&(w<=800)].mean()
    return w,a-base,base
