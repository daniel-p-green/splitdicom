"""Compare segmentations and HU measurements: base layers vs original."""
import sys, json, numpy as np, nibabel as nib
from totalsegmentator.map_to_binary import class_map
segdir=sys.argv[1]; tag=sys.argv[2]; names=class_map['total']; MINVOX=1000
out={}
for n in ['chest2mm','wbtorso','petct375']:
    try: so=np.asarray(nib.load(f'{segdir}/{n}_orig.nii.gz').dataobj)
    except Exception: continue
    vo=np.asarray(nib.load(f'vols/{n}_orig.nii.gz').dataobj).astype(np.float64)
    cnt=np.bincount(so.ravel(),minlength=120); labs=[l for l in range(1,len(cnt)) if cnt[l]>=MINVOX]
    out[n]={'labels':len(labs)}
    for k in ['rand1','2','4','8','16']:
        try: sk=np.asarray(nib.load(f'{segdir}/{n}_{k}.nii.gz').dataobj)
        except Exception: continue
        vk=np.asarray(nib.load(f'vols/{n}_{k}.nii.gz').dataobj).astype(np.float64)
        ck=np.bincount(sk.ravel(),minlength=120); d={}; vol={}; hu={}
        for l in labs:
            mo=so==l; inter=np.count_nonzero(mo&(sk==l)); d[l]=2*inter/(cnt[l]+ck[l]); vol[l]=100*(ck[l]-cnt[l])/cnt[l]
            hu[l]=abs(vk[mo].mean()-vo[mo].mean())
        dv=np.array(list(d.values())); worst=min(d,key=d.get)
        out[n][k]={'dice_median':round(float(np.median(dv)),4),'dice_p5':round(float(np.percentile(dv,5)),4),'dice_min':round(float(dv.min()),4),
            'worst':names.get(worst,str(worst)),'below_0.95':int((dv<0.95).sum()),'vol_change_max_pct':round(float(np.max(np.abs(list(vol.values())))),2),
            'vol_change_median_pct':round(float(np.median(np.abs(list(vol.values())))),3),
            'voxels_changed_pct':round(100*float(np.count_nonzero(so!=sk))/max(1,np.count_nonzero((so>0)|(sk>0))),3),
            'organ_meanHU_shift_max':round(float(max(hu.values())),3),'organ_meanHU_shift_median':round(float(np.median(list(hu.values()))),3),
            'rmse_HU':round(float(np.sqrt(((vk-vo)**2).mean())),2)}
        print(n,k,out[n][k],flush=True)
json.dump(out,open(f'compare_{tag}.json','w'),indent=1)
