"""Split-fidelity benchmark: bounded-error base + exact residual vs direct lossless."""
import sys, os, json, lzma, hashlib, collections, warnings
import numpy as np, pydicom, imagecodecs as ic
warnings.filterwarnings('ignore')
NEARS=[1,2,3,5]
SERIES={ # label: (dir, modality, series description prefix)
 'CT chest 2.0mm (uncompressed src)':('ohif/dcm/scoord3d-and-scoord/scoord-bounding-box','CT',None),
 'CT whole-body 5.0mm':('ohif/dcm/Juno','CT',None),
 'CT 3.75mm (PET/CT)':('ohif/dcm/acrin','CT',None),
 'MR MPRAGE 1.33mm':('ohif/dcm/FreeSurferSeg','MR','MPRAGE'),
 'MR T1w 0.67mm':('mr/dicoms','MR',None),
 'MR T1 3mm 512px':('ohif/dcm/MisterMr','MR','COR T-1'),
 'DX chest radiograph':('ohif/dcm/overlay','DX',None),
}
def noise_sd(a):
    d=(a[:,1:].astype(np.int64)-a[:,:-1]).ravel()
    return 1.4826*np.median(np.abs(d-np.median(d)))/np.sqrt(2)
def run(root,label,spec):
    d,mod,desc=spec; tot=collections.Counter(); n=0; maxerr={k:0 for k in NEARS}; ns=[]
    for dp,_,fs in os.walk(os.path.join(root,d)):
        for f in sorted(fs):
            p=os.path.join(dp,f)
            try: ds=pydicom.dcmread(p)
            except Exception: continue
            if ds.get('Modality')!=mod or 'PixelData' not in ds: continue
            if desc and not str(ds.get('SeriesDescription','')).startswith(desc): continue
            a=ds.pixel_array
            if a.ndim!=2: continue
            off=int(a.min()); u=(a.astype(np.int64)-off).astype(np.uint16)   # unsigned view, offset stored
            n+=1; tot['raw']+=a.size*2; ns.append(noise_sd(u))
            tot['JPEG-LS']+=len(ic.jpegls_encode(u))
            tot['JPEG2000']+=len(ic.jpeg2k_encode(u,level=0,reversible=True))
            tot['HTJ2K']+=len(ic.htj2k_encode(u)) if hasattr(ic,'htj2k_encode') else 0
            jx=ic.jpegxl_encode(u,lossless=True,effort=7); tot['JPEG-XL']+=len(jx)
            assert np.array_equal(ic.jpegxl_decode(jx),u)
            for k in NEARS:
                base=ic.jpegls_encode(u,level=k)
                b=ic.jpegls_decode(base).astype(np.int64)
                r=u.astype(np.int64)-b
                e=int(np.abs(r).max()); maxerr[k]=max(maxerr[k],e); assert e<=k
                r8=(r+k).astype(np.uint8)
                res=min(len(ic.jpegxl_encode(r8,lossless=True,effort=7)), len(lzma.compress(r8.tobytes(),preset=6)), len(ic.jpegls_encode(r8)))
                rec=(b+r).astype(np.uint16)                       # reconstruct
                assert hashlib.sha256(rec.tobytes()).digest()==hashlib.sha256(u.tobytes()).digest()
                tot[f'base{k}']+=len(base); tot[f'res{k}']+=res
    out={'label':label,'slices':n,'noise_sd':float(np.median(ns)),'raw_MB':tot['raw']/1e6}
    for c in ['JPEG-LS','JPEG2000','HTJ2K','JPEG-XL']: out[c]=tot['raw']/tot[c] if tot[c] else None
    best=min(tot[c] for c in ['JPEG-LS','JPEG2000','HTJ2K','JPEG-XL'] if tot[c])
    for k in NEARS:
        out[f'n{k}']={'base_ratio':tot['raw']/tot[f'base{k}'],'residual_ratio':tot['raw']/tot[f'res{k}'],
            'total_ratio':tot['raw']/(tot[f'base{k}']+tot[f'res{k}']),'overhead_vs_best_lossless_pct':100*((tot[f'base{k}']+tot[f'res{k}'])/best-1),
            'overhead_vs_jpegls_pct':100*((tot[f'base{k}']+tot[f'res{k}'])/tot['JPEG-LS']-1),'max_err':maxerr[k]}
    return out
if __name__=='__main__':
    root=sys.argv[1]; res=[]
    for l,s in SERIES.items():
        o=run(root,l,s); res.append(o); print(json.dumps(o),flush=True)
    json.dump(res,open(sys.argv[2],'w'),indent=1)
