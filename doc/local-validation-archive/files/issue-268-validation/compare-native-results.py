from pathlib import Path
import csv, json, math
root=Path(__file__).resolve().parent/'native-comparison'
checks=[]
for before in sorted(root.glob('*-base.csv')):
    after=before.with_name(before.name.replace('-base.csv','-candidate.csv'))
    old=list(csv.reader(before.open(), delimiter=';'))
    new=list(csv.reader(after.open(), delimiter=';'))
    assert len(old)==len(new),(before,'row count')
    maximum=0.; count=0
    for a,b in zip(old,new):
        assert len(a)==len(b),(before,'column count')
        for x,y in zip(a,b):
            try: u,v=float(x),float(y)
            except ValueError:
                assert x==y,(before,x,y)
                continue
            assert math.isfinite(u) and math.isfinite(v),(before,x,y)
            relative=abs(u-v)/max(abs(u),abs(v)) if u or v else 0.
            maximum=max(maximum,relative);count+=1
            assert relative<2e-12,(before,x,y,relative)
    assert before.with_suffix('.stdout').read_bytes()==after.with_suffix('.stdout').read_bytes(),(before,'native predictions differ')
    assert before.with_suffix('.stderr').read_bytes()==after.with_suffix('.stderr').read_bytes(),(before,'diagnostics differ')
    checks.append(dict(case=before.stem.removesuffix('-base'),probabilities=count,
                       max_relative_error=maximum,byte_identical=before.read_bytes()==after.read_bytes()))
assert len(checks)==12, len(checks)
(root/'results.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
