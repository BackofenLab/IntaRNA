from pathlib import Path
import xml.etree.ElementTree as ET
import cairosvg

# Reuse only the generator's drawing primitives; never execute its old outputs.
primitives = Path('/home/alex/CodexProjects/IntaRNA/issue-258-validation/step9/build_recursions.py').read_text().split('# 1. Fixed-right-boundary')[0]
exec(primitives)
OUT = Path('/home/alex/CodexProjects/IntaRNA/issue-258-validation/review-diagrams')

s = SVG('IntaRNA.PredictorMfeEns2d.outside.svg',
        'Unseeded pair probabilities · fixed-right outside pass',
        'PredictorMfeEns2d  |  --model=P --mode=M --noSeed  |  reverse the existing hybrid recurrence', 1730)
s.lines(48,185,[
    'A pair i=(i₁,i₂) uses local energy coordinates; RNA 2 is reversed. i ≺ k increases both coordinates.',
    'Fix right pair j. Hⱼ(i) is the forward hybrid weight; Aⱼ(i)=∂Z/∂Hⱼ(i) is its outside adjoint.',
    'M(i) sums complete-interaction weights containing actual pair i, over every right boundary j.'
],18)
s.eq(48,285,'d=BW(E_init);   t(i,k)=BW(E_interLeft(i,k));   i⁺=i+(1,1);   s=t(i,i⁺)',23)
s.box(48,322,1344,170,'1  Initialize the same completed-boundary objective')
s.eq(73,389,'Aⱼ(i) ← B(i,j);       Z += B(i,j)Hⱼ(i)',27)
s.lines(73,431,[
    'B(i,j)=BW(energy.getE(i₁,j₁,i₂,j₂,0)) only when !Z_equal(Hⱼ(i),0) and isValidOutputSite() pass.',
    'Otherwise B=0. Complete ED / GU-end filters and native terminal terms match updateZ().'
],18)
s.box(48,520,1344,180,'2  Reverse i in increasing order; credit its explicit pair',fill='#edf6f2',color=TEAL)
s.eq(73,588,'M(i) += Aⱼ(i)Hⱼ(i)',29)
s.lines(73,631,[
    'Every represented suffix contains its left pair i exactly once. Process only nonzero Aⱼ(i) and Hⱼ(i).',
    'Then propagate Aⱼ(i) through each forward term below; later suffix states receive all enclosing mass.'
],18)
s.box(48,728,645,257,'3a  Lonely pairs allowed')
s.eq(73,795,'Hⱼ(j)=d',25)
s.eq(73,840,'Aⱼ(k) += Aⱼ(i)t(i,k)',25)
s.lines(73,883,[
    'Direct terminal: k=j; its pair is counted at j.',
    'Intermediate k: i ≺ k ≺ j, allowed loop,',
    'and !Z_equal(Hⱼ(k),0) as in the forward sum.'
],18)
s.box(720,728,672,257,'3b  noLP: direct stacking',fill='#fff7ed',color=ORANGE)
s.eq(743,795,'Aⱼ(i⁺) += Aⱼ(i)s      if i⁺ ≠ j',25)
s.eq(743,840,'M(j) += Aⱼ(i)s d      if i⁺ = j',25)
s.lines(743,883,[
    'A noLP terminal has Hⱼ(j)=0.',
    'Its second pair j is implicit in the two-pair',
    'initialization and is credited separately.'
],18)
s.box(48,1014,1344,297,'3c  noLP: stack, then a nonstack loop',fill='#fff7ed',color=ORANGE)
s.rna(96,1115,431,[(0,'i'),(.14,'i⁺'),(.58,'k'),(1,'j')],shade=(.58,1,'Hⱼ(k)'))
s.eq(584,1096,'a = Aⱼ(i)s t(i⁺,k)',25)
s.eq(584,1139,'Aⱼ(k) += a',25)
s.eq(584,1182,'M(i⁺) += a Hⱼ(k)',25)
s.lines(73,1249,[
    'Require complementary i⁺; i⁺ ≺ k ≺ j; k ≠ i⁺+(1,1); the same loop limits and epsilon gate on Hⱼ(k).',
    'Here i⁺ is an actual pair skipped by the DP-state jump. Unpaired loop positions receive no mass.'
],18)
s.box(48,1340,1344,210,'4  Commit successful regions and normalize only once',fill='#edf6f2',color=TEAL)
s.eq(73,1409,'P(i,j) = ∑ M_region(i,j) / ∑ Z_region',27)
s.lines(73,1451,[
    'Check the accumulated objective equals the unchanged Zall; map pairs with energy.getBasePair().',
    'Original target/query orientation is stored in the sink. Z=0 → NA; failed arithmetic → no result.',
    'Conditional on an allowed interaction; no unbound-state weight. No outside buffers when output is disabled.'
],18)
s.text(48,1590,'Exact for the existing forward recurrence, including its numerical-zero gates; no alternative ensemble is substituted.',18)
s.footer('PredictorMfeEns2d.cpp: fillHybridZ(), accumulatePairMasses(), predictRegion()',
         'PredictorMfeEns.cpp: addPartitionContribution(); BasePairProbabilities: addRegion(), finalize()')
s.save()

s = SVG('IntaRNA.PredictorMfeEns2dHeuristic.probabilities.svg',
        'Unseeded heuristic probabilities · candidate flow',
        'PredictorMfeEns2dHeuristic  |  --model=P --mode=H --noSeed  |  approximate retained ensemble',1800)
s.lines(48,185,[
    'Each retained cell C(i) represents one selected pair chain, with hybrid weight vᵢ and right boundary jᵢ.',
    'Pair coordinates increase on both strands internally; RNA 2 is reversed. Set i⁺=i+(1,1).',
    'd=BW(E_init), t(i,k)=BW(E_interLeft(i,k)), and B(i,j)=BW(energy.getE(i₁,j₁,i₂,j₂,0)).'
],18)
s.box(48,294,1344,172,'Approximate means the continuation ensemble is pruned',fill='#fff7ed',color=ORANGE)
s.lines(73,357,[
    'Visit i in decreasing order. Generate candidates from already selected child chains C(k).',
    'Every admitted candidate contributes to Z; only the candidate with the lowest complete site energy',
    'is retained as C(i) for later extensions. Strict improvement wins; native energy rounding and ties are preserved.'
],18)
s.box(48,494,1344,358,'1  Record which new pairs each candidate owns')
s.text(74,559,'Candidate',19,BLUE,bold=True)
s.text(570,559,'Hybrid weight v',19,BLUE,bold=True)
s.text(929,559,'Owned pairs S',19,BLUE,bold=True)
s.text(1190,559,'Child k',19,BLUE,bold=True)
rows=[
    ('LP: initial pair','d','{i}','none'),
    ('LP: loop / stack to k','t(i,k)vₖ','{i}','k'),
    ('noLP: initial stack','s d,  s=t(i,i⁺)','{i, i⁺}','none'),
    ('noLP: direct stacking','s vᵢ⁺','{i}','i⁺'),
    ('noLP: nonstack loop to k','s t(i⁺,k)vₖ','{i, i⁺}','k'),
]
for row,y in zip(rows,range(608,809,44)):
    s.text(74,y,row[0],19)
    s.eq(570,y,row[1],22)
    s.eq(929,y,row[2],22)
    s.eq(1190,y,row[3],22)
s.text(73,830,'Native pair, loop/span and right-end GU constraints apply; noLP loops exclude the adjacent child k=i⁺+(1,1).',17)
s.box(48,882,1344,207,'2  Initialize mass and flow from every admitted candidate',fill='#edf6f2',color=TEAL)
s.eq(73,950,'w = B(i,j)v;     Z += w;     M(p) += w  for p ∈ S;     F(k) += w  if a child exists',23)
s.lines(73,995,[
    'M starts at zero; F is accumulated incoming weight for a retained child chain. This includes unselected candidates.',
    'Use exactly updateZ() admission: !Z_equal(v,0), complete-site GU / ED filters, and finite arithmetic.',
    'The current candidate owns S; every remaining pair belongs to its selected child continuation.'
],18)
s.box(48,1119,1344,249,'3  Propagate incoming mass through selected chains')
s.rna(95,1221,368,[(0,'i'),(.48,'k'),(1,'jₖ')],shade=(.48,1,'C(k)'))
s.lines(521,1195,[
    'Sweep i in increasing order. Let Sᵢ and kᵢ describe C(i).',
    'M(p) += F(i) for every p ∈ Sᵢ;',
    'F(kᵢ) += F(i) when C(i) has a child.',
    'Every edge increases both indices, so all parents arrive first.'
],21,step=36)
s.text(73,1333,'No division by vᵢ is needed: a retained cell represents one chain, so its full incoming weight follows every pair.',18)
s.box(48,1397,1344,223,'4  Publish marginals of the retained candidate ensemble',fill='#edf6f2',color=TEAL)
s.eq(73,1465,'P(i,j) = ∑ M_region(i,j) / ∑ Z_region',27)
s.lines(73,1510,[
    'Map energy coordinates back to original strand order, merge disjoint regions, then validate and normalize.',
    'The denominator must equal the unchanged heuristic Zall. Mark the result approximate, including empty results.',
    'Missing continuation chains are not recovered. These marginals differ from the exact --mode=M ensemble.'
],18)
s.text(48,1660,'Only actual paired nucleotides receive mass; positions merely enclosed by a retained interaction remain unmarked.',18)
s.footer('PredictorMfeEns2dHeuristic.cpp: fillHybridZ(), recordProbabilityCandidate(), commitProbabilityRegion()',
         'BasePairProbabilities: markApproximate(), addRegion(), finalize(); PredictorMfeEns.cpp: updateZ()')
s.save()

for name in records:
    path=DEST/name
    doc=ET.parse(path)
    assert b'\r' not in path.read_bytes()
    for node in doc.iter():
        assert node.tag.rsplit('}',1)[-1] not in ('image','script')
    cairosvg.svg2png(url=str(path),write_to=str(OUT/(name+'.png')))
print('\n'.join(records))
