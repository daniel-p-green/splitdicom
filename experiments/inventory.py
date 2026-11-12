import sys, os, collections, pydicom
root=sys.argv[1]
ser=collections.defaultdict(lambda:[0,0,None])
for dp,_,fs in os.walk(root):
    if '.git' in dp: continue
    for f in fs:
        p=os.path.join(dp,f)
        try:
            ds=pydicom.dcmread(p,stop_before_pixels=True)
        except Exception: continue
        if 'Rows' not in ds: continue
        k=(os.path.relpath(dp,root)[:50],ds.get('Modality'),str(ds.file_meta.get('TransferSyntaxUID')),ds.get('Rows'),ds.get('Columns'),ds.get('BitsStored'),ds.get('SamplesPerPixel'),str(ds.get('SeriesDescription'))[:25],str(ds.get('SliceThickness')))
        ser[k][0]+=1; ser[k][1]+=os.path.getsize(p)
for k,v in sorted(ser.items(), key=lambda x:-x[1][1])[:45]:
    print(v[0], f"{v[1]/1e6:.1f}MB", k)
