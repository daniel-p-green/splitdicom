"""Build NIfTI volumes: original and bounded-error base layers (same codec as splitdicom.py)."""
import sys, os, collections, warnings, json
import numpy as np, pydicom, nibabel as nib, imagecodecs as ic
warnings.filterwarnings('ignore')
SER={'chest2mm':'ohif/dcm/scoord3d-and-scoord/scoord-bounding-box','wb5mm':'ohif/dcm/Juno','petct375':'ohif/dcm/acrin'}
NEARS=[2,4,8,16]
root,out=sys.argv[1],sys.argv[2]; os.makedirs(out,exist_ok=True); info={}
for name,d in SER.items():
    groups=collections.defaultdict(list)
    for dp,_,fs in os.walk(os.path.join(root,d)):
        for f in fs:
            try: ds=pydicom.dcmread(os.path.join(dp,f))
            except Exception: continue
            if ds.get('Modality')=='CT' and 'PixelData' in ds and 'ImagePositionPatient' in ds: groups[ds.SeriesInstanceUID].append(ds)
    sl=max(groups.values(),key=len)
    iop=np.array(sl[0].ImageOrientationPatient,float); rc,cc=iop[:3],iop[3:]; nrm=np.cross(rc,cc)
    sl.sort(key=lambda s: float(np.dot(nrm,np.array(s.ImagePositionPatient,float))))
    ps=[float(x) for x in sl[0].PixelSpacing]; p0=np.array(sl[0].ImagePositionPatient,float); p1=np.array(sl[-1].ImagePositionPatient,float)
    aff=np.eye(4); aff[:3,0]=rc*ps[1]; aff[:3,1]=cc*ps[0]; aff[:3,2]=(p1-p0)/(len(sl)-1); aff[:3,3]=p0; aff[:2,:]*=-1
    vols={k:[] for k in ['orig']+NEARS}; raw=0; size=collections.Counter()
    for s in sl:
        a=s.pixel_array; sl_,ic_=float(s.get('RescaleSlope',1)),float(s.get('RescaleIntercept',0))
        off=int(a.min()) if a.min()<0 else 0; u=(a.astype(np.int64)-off).astype(np.uint16)
        vols['orig'].append(np.round(a*sl_+ic_).astype(np.int16)); raw+=a.size*2
        for k in NEARS:
            cs=ic.jpegls_encode(u,level=k); b=ic.jpegls_decode(cs).astype(np.int64)+off
            assert np.abs(b-a).max()<=k; size[k]+=len(cs)
            vols[k].append(np.round(b*sl_+ic_).astype(np.int16))
    for k,v in vols.items():
        vol=np.stack(v).transpose(2,1,0)
        nib.save(nib.Nifti1Image(vol,aff),os.path.join(out,f'{name}_{k}.nii.gz'))
    info[name]={'slices':len(sl),'shape':list(vol.shape),'spacing':[ps[1],ps[0],float(np.linalg.norm(aff[:3,2]))],'desc':str(sl[0].get('SeriesDescription')),
                'base_ratio':{k:round(raw/size[k],2) for k in NEARS},'kernel':str(sl[0].get('ConvolutionKernel')),'slope':float(sl[0].get('RescaleSlope',1))}
    print(name,info[name],flush=True)
json.dump(info,open(os.path.join(out,'info.json'),'w'),indent=1)
