"""Render computed corrected outputs. No benchmark substitutes for a current curve.
Only the explicitly labelled Figure 2 historical markers read historical arrays.
"""
from pathlib import Path
import os,json
HERE=Path(__file__).resolve().parent
os.environ['MPLCONFIGDIR']=str(HERE/'outputs/.mplconfig')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,to_rgb
from matplotlib.patches import Rectangle,PathPatch
from matplotlib.path import Path as PlotPath
from matplotlib.collections import PatchCollection
OUT=HERE/'outputs';COLORS=['#2166ac','#e08a1e','#2ca25f','#7a4fb0','#c0392b','#202020']
plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','pdf.fonttype':42,'ps.fonttype':42,'legend.fontsize':8.5,'axes.grid':False})
def read(n):return json.loads((OUT/f'fig{n}.json').read_text())
def draw(ax,rows,x,y,**style):
    if rows:ax.plot([r[x] for r in rows],[r[y] for r in rows],**style)
def segments(ax,rows,x,y,color,label,diagnostics=True):
    # Break at failed/masked cases; opacity indicates structural admissibility.
    for good in ([True,False] if diagnostics else [True]):
        run=[];first=True
        for r in rows+[None]:
            keep=r is not None and r.get('status')=='computed' and bool(r.get('quantitative'))==good and y in r
            if keep:run.append(r)
            elif run:
                c=color if good else tuple(.25*np.array(to_rgb(color))+.75)
                draw(ax,run,x,y,color=c,lw=1.7 if good else 1.1,ls='-' if good else '--',label=label if first and good else None);run=[];first=False

def render(n):
    d=read(n)
    if n==2:
        fig,(a,b)=plt.subplots(1,2,figsize=(11,4.3))
        for c,F in zip(COLORS,[1.,1.5,2.]):draw(a,[r for r in d if r['kind']=='forced' and r['F0']==F],'frequency','energy_ratio',color=c,label=f'$F_0={F:g}$ N')
        rows=[r for r in d if r['kind']=='rigid'];draw(b,rows,'natural_frequency','wingbeat_frequency',color=COLORS[4],label='corrected model')
        h=json.loads((HERE/'reference/fig2_original_source.json').read_text())[1][0];b.plot(h['x'],h['y'],'ko',ms=3,label='historical model samples')
        x=[r['natural_frequency'] for r in rows];b.plot(x,x,'--',color='.5',label='resonance')
        a.set(xlabel='flap frequency (Hz)',ylabel='flexible / rigid input energy',title='(a) Diagnostic prescribed-forcing benchmark');b.set(xlabel='natural frequency (Hz)',ylabel='wingbeat frequency (Hz)',title='(b)');a.legend();b.legend()
    elif n==3:
        tr=np.load(OUT/'cycles/fig3_transient.npz');cy=np.load(OUT/'cycles/fig3.npz');fig=plt.figure(figsize=(11,4.5));a=fig.add_subplot(121);b=fig.add_subplot(122,projection='3d');t,y=tr['t'],tr['y'];mask=t<=105
        for i,label,c in [(0,r'$\theta^*$',COLORS[0]),(2,r'$q^*$',COLORS[1]),(4,r'$f^*$',COLORS[3])]:a.plot(t[mask],y[i,mask],label=label,color=c)
        b.plot(y[0,mask][::8],y[2,mask][::8],y[4,mask][::8],color='.7',lw=.8);b.plot(cy['theta'],cy['q'],cy['muscle_state'],color=COLORS[0]);b.view_init(22,-58)
        a.set(xlabel=r'$t^*$',ylabel='amplitude',title='(a)');b.set(xlabel=r'$\theta^*$',ylabel=r'$q^*$',zlabel=r'$f^*$',title='(b)');a.legend()
    elif n==4:
        fig,(a,b)=plt.subplots(1,2,figsize=(11.5,4.6));on=[r for r in d if r['kind']=='onset'];draw(a,on[:90],'r_w','mu_H',color=COLORS[0]);a.fill_between([r['r_w'] for r in on[:90]],0,[r['mu_H'] for r in on[:90]],color=COLORS[0],alpha=.10,label='stable rest');a.axhline(0,color=COLORS[4],ls=':')
        for rw,c in zip([1.3,2.,3.5],COLORS):
            draw(b,sorted([r for r in d if r['kind']=='branch' and r['r_w']==rw],key=lambda r:r['mu_hat']),'mu_hat','stroke_amp',color=c,label=f'$r_w={rw:g}$')
            onset=next(r['mu_H'] for r in on if r['r_w']==rw)
            b.axvline(onset,color=c,ls='--',lw=1)
        a.set(xlabel=r'$r_w$',ylabel=r'$\hat\mu_H$',ylim=(-.02,1.85),title='(a)');b.set(xlabel=r'$\hat\mu$',ylabel=r'$\theta_A$',ylim=(-.003,.125),title='(b)');b.legend(loc='upper right')
    elif n==5:
        fig,a=plt.subplots(figsize=(6.7,4.8))
        for r,c in zip(d,COLORS):
            cy=np.load(OUT/'cycles'/('fig5_'+r['key']+'.npz'));a.plot(cy['theta'],r['mu_hat']*cy['muscle_state'],color=COLORS[4] if r['key']=='boundary' else c,ls=':' if r['key']=='boundary' else '-',label='rigid' if r['key']=='rigid' else f'$r_w={r["r_w"]:.4f}$')
            if r['key']=='boundary':a.fill(cy['theta'],r['mu_hat']*cy['muscle_state'],color=COLORS[4],alpha=.10)
        a.set(xlabel=r'$\theta^*$',ylabel=r'$\hat\mu f^*$');a.legend()
    elif n==6:
        fig,(a,b)=plt.subplots(1,2,figsize=(11.5,4.6))
        for rw,c in zip([1.4,2.,2.5,3.,5.5,7.5],COLORS):segments(a,sorted([r for r in d if r['panel']=='a' and r['r_w']==rw],key=lambda r:r['mu_hat']),'mu_hat','stroke_amp',c,f'$r_w={rw:g}$')
        for t,c in zip([3.,4.4,5.83,8.],COLORS):draw(b,sorted([r for r in d if r['panel']=='b' and r['t0_ms']==t],key=lambda r:r['r_w']),'r_w','stroke_amp',color=c,marker='o',ms=3,label=f'$t_o={t:g}$ ms')
        a.axhline(1,color='.5',ls='--');a.set(xlabel=r'$\hat\mu$',ylabel=r'$\theta_A$',xlim=(0,90),ylim=(0,1.2));b.set(xlabel=r'$r_w$',ylabel=r'$\theta_{\rm sup}$');a.set_title('(a)',loc='left');b.set_title('(b)',loc='left');a.legend();b.legend()
    elif n==7:
        fig,(a,b)=plt.subplots(1,2,figsize=(11,4.3))
        for nw,c in zip([.3,.84,1.3,3.],COLORS):
            rr=sorted([r for r in d if r['N_W']==nw],key=lambda r:r['r_w'])
            for ax,k in [(a,'delta_phi_net_pp'),(b,'delta_phi_aero_pp')]:draw(ax,rr,'r_w',k,color=c,label=f'$N_W={nw:g}$')
        a.set_title('(a)',loc='left');b.set_title('(b)',loc='left')
        for ax,k in [(a,'net'),(b,'aero')]:ax.axhline(0,color='.6',lw=.7);ax.set(xlabel=r'$r_w$',ylabel=fr'$\Delta\varphi_{{{k}}}$ (percentage points)');ax.legend()
    elif n==8:
        fig,a=plt.subplots(figsize=(7.5,4.8))
        for t,c in zip([.2,.5,.8,1.],COLORS):segments(a,sorted([r for r in d if r['target']==t],key=lambda r:r['r_w']),'r_w','G_log',c,rf'$\theta_A={t:g}$')
        for b,c in zip(json.loads((OUT/'fig8_boundaries.json').read_text()),COLORS):a.plot(b['reported_rw'],b['point']['G_log'],'o',mfc='white',mec=c,ms=5)
        a.set(xlabel=r'$r_w$',ylabel=r'$G_A^{\log}$',xlim=(1.45,7.62));a.legend()
    elif n==9:
        fig,a=plt.subplots(figsize=(8.2,4.7));b=a.twinx()
        for tag,c in zip(['rigid','rw2p5','rw2'],COLORS):
            rr=sorted([r for r in d if r['tag']==tag],key=lambda r:r['t0']);rr=[dict(r,t0_ms=1000*r['t0']) for r in rr]
            label={'rigid':'rigid', 'rw2p5':r'$r_w=2.5$', 'rw2':r'$r_w=2.0$'}[tag]
            segments(a,rr,'t0_ms','Omega',c,label,diagnostics=False);segments(b,rr,'t0_ms','G_log',c,label,diagnostics=False)
            supported=[r for r in rr if r.get('status')=='computed' and r.get('quantitative')]
            if tag!='rigid' and supported:
                endpoint=supported[0]
                for ax,field in [(a,'Omega'),(b,'G_log')]:ax.plot(endpoint['t0_ms'],endpoint[field],marker='o',mfc='white',mec=c,ls='none',ms=5,zorder=5)
        for line in b.lines:line.set_linestyle('--')
        a.axvline(5.83,color='.5',ls=':',label='selected baseline');a.set(xlabel=r'$t_o$ (ms)',ylabel=r'$\Omega$');b.set_ylabel(r'$\widehat{G}_f^{\,t_o}$');fig.legend(*a.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,1.),ncol=4)
        from matplotlib.lines import Line2D
        a.legend(handles=[Line2D([0],[0],color='.25',ls='-',label=r'$\Omega$'),Line2D([0],[0],color='.25',ls='--',label=r'$\widehat{G}_f^{\,t_o}$')],loc='lower left')
    elif n==10:
        fig,(a,b)=plt.subplots(1,2,figsize=(10.5,4.3));rr=sorted([r for r in d if r['key']!='rigid'],key=lambda r:r['r_w'])
        for c,k,label in zip(COLORS,['W_plus_over_rigid','peak_muscle_force_over_rigid','peak_positive_muscle_power_over_rigid'],['positive work','peak force','peak positive power']):draw(a,rr,'r_w',k,color=c,marker='o',ms=3,label=label)
        for c,k in zip(COLORS,['phi_net','phi_aero']):b.plot([r['r_w'] for r in rr],[100*r[k] for r in rr],color=c,marker='o',ms=3,label={'phi_net':r'$\varphi_{\rm net}$','phi_aero':r'$\varphi_{\rm aero}$'}[k])
        for ax in [a,b]:ax.axvline(rr[0]['r_w'],color='.6',ls=':',lw=.8);ax.set_xlabel(r'$r_w$');ax.legend()
        a.set_ylabel('demand / matched rigid demand');b.set_ylabel('fraction of positive work (%)');a.set_title('(a)',loc='left');b.set_title('(b)',loc='left')
    elif n==11:
        fig,a=plt.subplots(figsize=(7.3,5.2));valid=[r for r in d if r.get('status')=='computed' and r.get('quantitative')];nws=sorted(set(r['N_W'] for r in valid));edges=np.r_[nws[0]/np.sqrt(nws[1]/nws[0]),np.sqrt(np.array(nws[:-1])*nws[1:]),nws[-1]*np.sqrt(nws[-1]/nws[-2])];patches=[];values=[]
        for j,nw in enumerate(nws):
            rr=sorted([r for r in valid if r['N_W']==nw],key=lambda r:r['r_w']);x=np.array([r['r_w'] for r in rr]);xe=np.r_[x[0],(x[:-1]+x[1:])/2,5.6]
            for i,r in enumerate(rr):patches.append(Rectangle((xe[i],edges[j]),xe[i+1]-xe[i],edges[j+1]-edges[j]));values.append(r['Pi_flex'])
        pc=PatchCollection(patches,cmap='magma',norm=LogNorm(vmin=.28,vmax=max(values)),edgecolor='none');pc.set_array(np.array(values));a.add_collection(pc)
        bb=json.loads((OUT/'fig11_boundaries.json').read_text());xs=[r['reported_rw'] for r in bb];yy=[r['N_W'] for r in bb];ys=np.geomspace(min(yy),max(yy),300);xx=np.interp(np.log(ys),np.log(yy),xs);verts=[(5.6,ys[0])]+list(zip(xx,ys))+[(5.6,ys[-1]),(5.6,ys[0])];pc.set_clip_path(PathPatch(PlotPath(verts),transform=a.transData));a.plot(xs,yy,'k-');a.set(xlabel=r'$r_w$',ylabel=r'$N_W$',yscale='log',xlim=(1.3,5.6),ylim=(.285,3.15));fig.colorbar(pc,ax=a,label=r'$\Pi_{\rm flex}$')
    else:raise ValueError(n)
    fig.tight_layout(rect=(0,0,1,.90) if n==9 else (0,0,1,1));dest=OUT/'figures';dest.mkdir(parents=True,exist_ok=True)
    for ext in ['pdf','png']:fig.savefig(dest/f'fig{n}.{ext}',bbox_inches='tight',dpi=300)
    plt.close(fig)
