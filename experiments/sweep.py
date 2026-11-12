"""Sweep the error bound: how large must the bound be before the base layer gets small?"""
import sys, os, json, lzma, warnings, collections
import numpy as np, pydicom, imagecodecs as ic
warnings.filterwarnings('ignore')
from bench import SERIES
NEARS=[2,4,8,12,16,24]
def frames(root,spec,step):
    d,mod,desc=spec; i=0
    for dp,_,fs in os.walk(os.path.join(root,d)):
        for f in sorted(fs):
            try: ds=pydicom.dcmread(os.path.join(dp,f))
            except Exception: continue
            if ds.get('Modality')!=mod or 'PixelData' not in ds: continue
            if desc and not str(ds.get('SeriesDescription','')).startswith(desc): continue
            i+=1
            if i%step: continue
            a=ds.pixel_array
            if a.ndim==2: yield (a.astype(np.int64)-int(a.min())).astype(np.uint16)
root=sys.argv[1]; out=[]
for label in ['CT chest 2.0mm (uncompressed src)','CT whole-body 5.0mm','CT 3.75mm (PET/CT)','MR MPRAGE 1.33mm','MR T1 3mm 512px']:
    t=collections.Counter()
    for u in frames(root,SERIES[label],4):
        t['raw']+=u.size*2; t['ll']+=len(ic.jpegxl_encode(u,lossless=True,effort=5))
        for k in NEARS:
            # A: JPEG-LS near-lossless base
            ba=ic.jpegls_encode(u,level=k); b=ic.jpegls_decode(ba).astype(np.int64)
            r=u.astype(np.int64)-b; assert np.abs(r).max()<=k
            ra=len(ic.jpegxl_encode((r+k).astype(np.uint8),lossless=True,effort=5))
            # B: uniform quantiser (step 2k+1) coded losslessly with JPEG XL
            q=2*k+1; qi=((u.astype(np.int64)+k)//q); bq=qi*q; r2=u.astype(np.int64)-bq; assert np.abs(r2).max()<=k
            bb=len(ic.jpegxl_encode(qi.astype(np.uint16),lossless=True,effort=5))
            rb=min(len(ic.jpegxl_encode((r2+k).astype(np.uint8),lossless=True,effort=5)),len(lzma.compress((r2+k).astype(np.uint8).tobytes())))
            assert np.array_equal(bq+r2,u)
            t[f'A{k}b']+=len(ba); t[f'A{k}r']+=ra; t[f'B{k}b']+=bb; t[f'B{k}r']+=rb
    row={'label':label,'lossless':round(t['raw']/t['ll'],2)}
    for k in NEARS:
        for m in 'AB':
            row[f'{m}{k}']={'base':round(t['raw']/t[f'{m}{k}b'],1),'res':round(t['raw']/t[f'{m}{k}r'],1),'ovh%':round(100*((t[f'{m}{k}b']+t[f'{m}{k}r'])/t['ll']-1),1)}
    out.append(row); print(json.dumps(row),flush=True)
json.dump(out,open(sys.argv[2],'w'),indent=1)
