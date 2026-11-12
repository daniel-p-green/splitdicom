import sys, os, tempfile, warnings, collections
warnings.filterwarnings('ignore')
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import pydicom, splitdicom as sd
from bench import SERIES
root=sys.argv[1]; tmp=tempfile.mkdtemp(); ok=bad=0; fails=collections.Counter(); independent=0
for label,(d,mod,desc) in SERIES.items():
    for dp,_,fs in os.walk(os.path.join(root,d)):
        for f in sorted(fs):
            p=os.path.join(dp,f)
            try: ds=pydicom.dcmread(p,stop_before_pixels=True)
            except Exception: continue
            if ds.get('Modality')!=mod or 'Rows' not in ds: continue
            if desc and not str(ds.get('SeriesDescription','')).startswith(desc): continue
            b,r,o=[os.path.join(tmp,x) for x in ('b.dcm','r.sdr','o.dcm')]
            try:
                sd.encode(p,b,r,2); sd.decode(b,r,o); diffs=sd.verify(p,o)
                # independent check: base opens through pydicom's own decoder path and honours the bound
                if ok%50==0:
                    import numpy as np
                    bo=pydicom.dcmread(b); ba=bo.pixel_array.astype('int64'); oa=pydicom.dcmread(p).pixel_array.astype('int64')
                    off=int(oa.min()) if oa.min()<0 else 0
                    assert abs((ba+off)-oa).max()<=2; independent+=1
            except SystemExit as e: diffs=[str(e)]
            except Exception as e: diffs=[type(e).__name__+': '+str(e)[:80]]
            if diffs: bad+=1; fails[(label,str(diffs)[:120])]+=1
            else: ok+=1
print('round-trip OK:',ok,'failed:',bad,'| base opened by pydicom decoder and bound held on',independent,'samples')
for k,v in fails.items(): print(v,k)
