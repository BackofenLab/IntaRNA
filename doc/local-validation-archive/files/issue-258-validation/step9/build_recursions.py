from pathlib import Path
from html import escape
import json
from PIL import ImageFont
from fontTools.ttLib import TTFont
SERIF_GLYPHS = TTFont("/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf").getBestCmap()

DEST = Path('/home/alex/CodexProjects/IntaRNA/pr-worktrees/issue-258-bpsvg/doc/recursions')
BLUE, TEAL, ORANGE, INK = '#275f8c', '#23796a', '#bd6518', '#203244'
records = []


class SVG:
    def __init__(self, name, title, subtitle, height=1150):
        self.name, self.height = name, height
        self.parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="1440" height="{height}" viewBox="0 0 1440 {height}" role="img" aria-labelledby="title desc">',
                      f'<title id="title">{escape(title)}</title>',
                      f'<desc id="desc">{escape(subtitle)}. Internal recurrence documentation for step 9 of PR 258. Equations describe the current implementation; RNA 2 uses reversed energy coordinates.</desc>',
                      '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 Z" fill="#60798e"/></marker>',
                      '<pattern id="seed-hatch" width="8" height="8" patternUnits="userSpaceOnUse"><rect width="8" height="8" fill="#fff0dc"/><path d="M-2 2 L2 -2 M0 8 L8 0 M6 10 L10 6" stroke="#e0b27a" stroke-width="1"/></pattern></defs>',
                      '<rect width="1440" height="100%" fill="white"/>']
        self.text(48, 48, 'INTARNA  /  RECURRENCE NOTES', 14, BLUE, bold=True)
        self.text(48, 91, title, 31, bold=True)
        self.text(48, 125, subtitle, 18)
        self.line(48, 148, 1392, 148)

    def text(self, x, y, text, size=19, color=INK, anchor='start', bold=False, math=False):
        family = 'DejaVu Serif' if math else 'DejaVu Sans'
        font = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/'+('DejaVuSerif' if math else 'DejaVuSans')+('-Bold' if bold else '')+'.ttf', size)
        width = font.getlength(text)
        left = x - (width/2 if anchor == 'middle' else width if anchor == 'end' else 0)
        assert left >= 16 and left + width <= 1424, (self.name, text, left, width)
        assert 15 <= y <= self.height-12, (self.name, text, y)
        rendered = ''.join(escape(c) if not math or ord(c) in SERIF_GLYPHS else f'<tspan font-family="DejaVu Sans, sans-serif">{escape(c)}</tspan>' for c in text)
        self.parts.append(f'<text x="{x}" y="{y}" font-family="{family}, sans-serif" font-size="{size}" fill="{color}" text-anchor="{anchor}"'+(' font-weight="bold"' if bold else '')+f'>{rendered}</text>')

    def lines(self, x, y, lines, size=19, step=29, **kw):
        for i, line in enumerate(lines):
            self.text(x, y+i*step, line, size, **kw)

    def line(self, x1, y1, x2, y2, color='#c4d0da', width=1.2, arrow=False, dash=False):
        self.parts.append(f'<path d="M{x1} {y1} L{x2} {y2}" fill="none" stroke="{color}" stroke-width="{width}"'+(' marker-end="url(#arrow)"' if arrow else '')+(' stroke-dasharray="5 4"' if dash else '')+'/>')

    def box(self, x, y, w, h, title=None, fill='#f3f7fa', color=BLUE):
        self.parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="{fill}" stroke="#d0dbe3"/>')
        if title:
            self.text(x+22, y+34, title, 21, color, bold=True)

    def eq(self, x, y, text, size=25):
        self.text(x, y, text, size, math=True)

    def rna(self, x, y, w, labels, shade=None, seed=None):
        # Pair positions along antiparallel strands; both internal indices increase.
        if shade:
            a,b,label=shade
            self.parts.append(f'<rect x="{x+w*a}" y="{y}" width="{w*(b-a)}" height="60" fill="#dcebf5"/>')
            self.text(x+w*(a+b)/2, y+37, label, 24, BLUE, anchor='middle', math=True)
        if seed:
            a,b=seed
            self.parts.append(f'<rect x="{x+w*a}" y="{y}" width="{w*(b-a)}" height="60" fill="url(#seed-hatch)"/>')
        self.line(x, y, x+w, y, INK, 2)
        self.line(x, y+60, x+w, y+60, INK, 2)
        for fraction,label in labels:
            px=x+w*fraction
            self.line(px,y,px,y+60,ORANGE if seed and seed[0]<=fraction<=seed[1] else INK,2.4)
            if label:
                self.text(px,y-12,label,17,anchor='middle',math=True)

    def footer(self, *sources):
        self.line(48, self.height-95, 1392, self.height-95)
        self.text(48,self.height-71,'Implementation anchors (src/IntaRNA/) · BW(E)=exp(−E/RT) after native energy-unit conversion',14,BLUE,bold=True)
        for i,source in enumerate(sources):
            self.text(48,self.height-48+21*i,source,13)

    def save(self):
        (DEST/self.name).write_text('\n'.join(self.parts+['</svg>'])+'\n')
        records.append(self.name)


# 1. Fixed-right-boundary unseeded partition.
s=SVG('IntaRNA.PredictorMfeEns2d.svg','Unseeded ensemble · fixed right boundary','PredictorMfeEns2d  |  --model=P --mode=M --noSeed',1190)
s.lines(48,184,['Pair indices i=(i₁,i₂), k and j increase on both strands; RNA 2 is reversed.',
                'i ≺ k means strict increase in both coordinates. d=BW(E_init), t(i,k)=BW(E_interLeft(i,k)).'])
s.eq(48,253,'Hⱼ(i) = hybrid partition for all admitted chains from i to a fixed j',25)
for x,title in [(48,'Initiation'),(500,'Direct pair-to-pair step'),(952,'Recursive continuation')]:
    s.box(x,285,440,225,title)
s.rna(106,364,260,[(.5,'i = j')]);s.eq(268,477,'d',28)
s.rna(553,364,330,[(0,'i'),(1,'j')]);s.eq(607,477,'t(i,j) · d',27)
s.rna(1001,364,330,[(0,'i'),(.4,'k'),(1,'j')],shade=(.4,1,'Hⱼ(k)'));s.eq(1001,477,'t(i,k) · Hⱼ(k)',25)
s.text(48,554,'Lonely pairs allowed',23,BLUE,bold=True)
s.eq(48,599,'Hⱼ(j)=d;   Hⱼ(i)=t(i,j)d + ∑ t(i,k)Hⱼ(k)       (i ≺ k ≺ j)',25)
s.lines(48,639,['All other or noncomplementary states are zero. Sum only loop-size-admissible steps.',
                'For every valid right pair j, fill i in decreasing order within the span limits.'])
s.box(48,697,1344,231,'noLP: the left boundary must stack to i⁺ = i + (1,1)',fill='#edf6f2',color=TEAL)
s.eq(72,776,'Hⱼ(j)=0;   Hⱼ(i)=t(i,i⁺) · [ (d if i⁺=j, else Hⱼ(i⁺)) + ∑ t(i⁺,k)Hⱼ(k) ]',24)
s.lines(72,821,['The sum requires i⁺ ≺ k ≺ j and k ≠ i⁺+(1,1): the direct stack is already in Hⱼ(i⁺).',
                'Require a complementary i⁺. Unreachable suffixes are zero; the adjacent two-pair case contributes t(i,j)d.',
                'Thus initialization is a stack, and each nonstack continuation starts another stacked suffix.'],18)
s.eq(48,978,'Zall += B(i,j) · Hⱼ(i);    E_site = energy.getE(i,j, energy.getE(Hⱼ(i)))',24)
s.lines(48,1017,['B(i,j)=BW(energy.getE(i₁,j₁,i₂,j₂,0)) for admitted boundaries; otherwise B=0.',
                 'This adds full-interval ED and terminal terms once; updateCompleteZ submits each boundary once.'],18)
s.footer('PredictorMfeEns2d.cpp: predict(), fillHybridZ()', 'PredictorMfeEns.cpp: updateCompleteZ(), addPartitionContribution()');s.save()

# 2. One successor retained at each left boundary.
s=SVG('IntaRNA.PredictorMfeEns2dHeuristic.svg','Unseeded ensemble · heuristic continuation','PredictorMfeEns2dHeuristic  |  --model=P --mode=H --noSeed',1140)
s.lines(48,185,['C(i)=(vᵢ,jᵢ) retains one hybrid weight and right boundary for each left pair i.',
                'Use d=BW(E_init), t(i,k)=BW(E_interLeft(i,k)); all pair coordinates increase internally.'])
s.text(48,241,'B(i,j)=BW(energy.getE(i₁,j₁,i₂,j₂,0)) is the complete-boundary Boltzmann factor.',17)
s.box(48,254,1344,247,'Enumerate candidates at i, visiting i in decreasing order')
s.rna(82,349,375,[(0,'i'),(.4,'k'),(1,'jₖ')],shade=(.4,1,'vₖ'))
s.lines(505,319,['Lonely pairs allowed:','  initial candidate (d, i), only if !noGUend or i is not GU',
                  '  continuation (t(i,k) · vₖ, jₖ) for each reachable k ≻ i',
                  '  within loop limits and the complete span limits.'],21,step=38)
s.line(712,510,712,561,arrow=True,width=2)
s.box(48,581,645,178,'Accumulate every considered candidate',fill='#edf6f2',color=TEAL)
s.eq(71,659,'updateZ(i,j,v,true)',26)
s.lines(71,697,['For admitted complete sites: add B(i,j)v to Zall',
                'and v to the boundary map; otherwise discard the update.'],18)
s.box(720,581,672,178,'Retain only the best continuation at i')
s.eq(743,659,'C(i) ← argmin E_full(i,j,v)',25)
s.lines(743,697,['Compare complete energy for each candidate.',
                 'Strict improvement replaces vᵢ and jᵢ.'],18)
s.box(48,795,1344,188,'noLP candidate set',fill='#fff7ed',color=ORANGE)
s.lines(72,859,['Let i⁺=i+(1,1) be complementary and s=t(i,i⁺).',
                'Initialize (s·d, i⁺) only if !noGUend or i⁺ is not GU; direct continuation is (s·vᵢ⁺, jᵢ⁺).',
                'Other candidates are (s·t(i⁺,k)·vₖ, jₖ), with k ≻ i⁺ and k ≠ i⁺+(1,1).',
                'The adjacent continuation is excluded from the loop sum because it was handled directly.'],18)
s.lines(48,1020,['The boundary map sums contributions, but the DP cell keeps one right continuation: this is a restricted ensemble.'],18)
s.footer('PredictorMfeEns2dHeuristic.cpp: fillHybridZ()', 'PredictorMfeEns.cpp: updateZ(), updateOptimaUsingZ()');s.save()

# 3. Seed-only mode is enumeration, not an all-seed-substructure partition.
s=SVG('IntaRNA.PredictorMfeEnsSeedOnly.svg','Seed-only ensemble · retained seed enumeration','PredictorMfeEnsSeedOnly  |  --model=P --mode=S',990)
s.box(48,181,1344,160,'1  Fill the seed handler for the active region')
s.lines(72,246,['Retain the handler-selected seed at each admitted start a; its end is e.',
                'Both endpoints must lie in the region. The handler supplies seed pairs and hybrid energy E_seed(a).',
                'A filled handler with no admitted seeds reports an empty result.'],19)
s.line(720,353,720,402,arrow=True,width=2)
s.box(48,420,1344,240,'2  One complete contribution per retained seed',fill='#fff7ed',color=ORANGE)
s.rna(89,524,391,[(0,'a'),(.25,''),(.5,''),(.75,''),(1,'e')],seed=(0,1))
s.eq(540,522,'H_seed = E_seed(a) + E_init',26)
s.eq(540,566,'E_full = energy.getE(a₁,e₁,a₂,e₂,H_seed)',24)
s.eq(540,605,'Zall += BW(E_full)',23)
s.text(540,637,'Only when needZall is enabled and complete-site filters pass.',18)
s.box(48,702,1344,154,'3  Report and trace the retained seed')
s.lines(72,765,['updateOptima(..., isHybridE=true, incrementZall=true) applies the complete-site filters and ranking.',
                'traceBackSeed() supplies its interior pairs; initiation is charged exactly once.',
                'There is no extension recurrence, no sum over alternative seed structures at the same start, and no bpProb output.'],18)
s.footer('PredictorMfeEnsSeedOnly.cpp: predict(), traceBack()', 'PredictorMfe.cpp: updateOptima(); Predictor.h: updateZall() (inherited update path)');s.save()

# 4. Exact stack-seed adapter and ownership, including legacy dispatch.
s=SVG('IntaRNA.PredictorMfeEns2dSeedExtension.svg','Seeded ensemble · exact stack-seed adapter','PredictorMfeEns2dSeedExtension  |  --model=P --mode=M with stack-only seed capability',1190)
s.box(48,177,645,128,'SeedHandler: guaranteesStackOnlySeeds()')
s.lines(71,238,['true → predictStackSeeds(): unique seeded-chain partition',
                'false → legacy anchor extensions; no pair-probability sink'],18)
s.box(720,177,672,128,'Capability is geometric',fill='#fff7ed',color=ORANGE)
s.lines(743,238,['Computed or explicit seeds; mixed lengths and singletons.',
                 'Preserve handler admission and per-start selection rules.'],18)
s.box(48,344,1344,223,'1  Compile active seed occurrences and eligible left boundaries')
s.rna(95,438,398,[(0,'p'),(.40,'a'),(.57,''),(.72,'b'),(1,'q')],seed=(.4,.72))
s.eq(546,433,'max(0, bₖ − Wₖ + 1) ≤ pₖ ≤ aₖ     (k=1,2)',24)
s.lines(546,478,['Admitted stack a…b must fit both the region and span limits W.',
                 'Enumerate the union of these p-boxes once, then require a valid p.',
                 'At p, clip the working grid to both inclusive spans.'],18)
s.box(48,605,1344,227,'2  Give the suffix kernel one consistent energy objective',fill='#edf6f2',color=TEAL)
s.eq(73,677,'d(p)=BW(E_init);    t(r,q)=BW(E_interLeft(r₁,q₁,r₂,q₂))',25)
s.eq(73,726,'B(p,q)=BW(energy.getE(p₁,q₁,p₂,q₂,0))    if the complete site is admitted',24)
s.lines(73,766,['Invalid / inaccessible / noncomplementary vertices and forbidden transitions cannot contribute.',
                'B owns complete-span ED, native dangles, terminal penalties and energyAdd. Seed energies only admit seeds.'],18)
s.box(48,870,1344,173,'3  Commit a completed region')
s.eq(73,937,'Z_region = ∑ B(p,q) · T₁,p(q)',27)
s.lines(580,915,['Complete boundary callback and outside initialization share exactly the same B.',
                  'updateExactCompleteZ preserves the partition objective.',
                  'Pair masses and optional seed masks map through energy.getBasePair().',
                  'Only a successful computation reaches BasePairProbabilities::addRegion().'],17,step=29)
s.footer('PredictorMfeEns2dSeedExtension.cpp: predict(), predictStackSeeds(), exactBoundaryWeight()', 'StackSeedDomain.cpp; SeededPartitionFunction.cpp; PredictorMfeEns.cpp: updateExactCompleteZ()');s.save()

# 5. Forward state machine used in production.
s=SVG('IntaRNA.SeededPartitionFunction.inside.svg','Exact seeded inside recurrence · suffix automaton','SeededPartitionFunction::compute()  |  fixed left boundary p; visit q in increasing pair order',1475)
s.lines(48,185,['C=max(2, longest admitted seed length). ℓ is the current consecutive stack-suffix length, capped at C.',
                'e(q)=shortest admitted seed ending at q, or 0. A run accepts a seed iff 0 < e(q) ≤ ℓ.'])
for x,w,title,lines in [(48,436,'Uℓ: still seed-free',['ℓ=1…C; no seed has occurred.', 'Stack increments ℓ; a loop resets it.']),
                         (502,436,'S₁: seeded, singleton suffix',['A seed has occurred; current run is length 1.', 'Under noLP, must stack before leaving.']),
                         (956,436,'S₂: seeded, usable suffix',['Under noLP: current run has length ≥2.', 'Otherwise: all seeded suffixes use S₂.'])]:
    s.box(x,250,w,147,title,fill='#f3f7fa' if x==48 else '#fff7ed',color=BLUE if x==48 else ORANGE)
    s.lines(x+20,316,lines,16,step=28)
s.line(266,397,266,422,BLUE,1.6)
s.line(266,422,1174,422,BLUE,1.6)
s.line(1174,422,1174,397,BLUE,1.6,arrow=True)
s.text(720,417,'first seed accepted on a stack',16,BLUE,anchor='middle')
s.text(720,377,'stack → S₂',16,ORANGE,anchor='middle')
s.parts.append('<g transform="translate(0 55)">')
s.box(48,437,1344,221,'1  Initialize and extend every state')
s.eq(73,505,'H(p,c₀)=d(p);   other H=0.    c₀ = accept(p,1) ? (noLP ? S₁ : S₂) : U₁',23)
s.eq(73,553,'Stack r=q−(1,1):     H(q,D(c,stack)) += t(r,q) · H(r,c)',25)
s.lines(73,599,['For Uℓ: let ℓ′=min(C,ℓ+1); D=Uℓ′ unless the ending seed is accepted, then D=S₂.',
                'For S₁ or S₂: D=S₂. Once seeded, a later seed never creates a second copy of the path.'],18)
s.box(48,697,1344,260,'2  Add nonstack arrivals, then aggregate eligible departures',fill='#edf6f2',color=TEAL)
s.eq(73,764,'For each nonstack r→q and b∈{0,1}:  H(q,D(b,loop)) += t(r,q) · T_b(r)',24)
s.eq(73,809,'T₀(q)=∑ H(q,Uℓ)   over ℓ≥(noLP ? 2 : 1)',25)
s.eq(73,854,'T₁(q)=H(q,S₂) + (noLP ? 0 : H(q,S₁))',25)
s.lines(73,891,['T_b(r) is already computed. Nonstack: strict increase within loop limits, excluding r=q−(1,1).',
                'Reset the suffix to 1; remain seeded if b=1, or accept a singleton seed at q if b=0.'],18)
s.box(48,995,1344,174,'3  Complete seeded paths at every q',fill='#fff7ed',color=ORANGE)
s.eq(73,1064,'Z += B(p,q) · T₁(q)',29)
s.lines(532,1049,['Seed-free states never terminate an interaction.',
                  'With noLP, singleton runs cannot terminate or cross a nonstack edge.',
                  'Such states remain available for a later stack.'],18)
s.lines(48,1214,['Loop destinations: seeded → S₁ under noLP, otherwise S₂; unseeded → U₁ unless a singleton seed accepts.',
                 'The Boolean seeded state counts each ordered pair chain once, including paths with several seed witnesses.',
                 'Termination-only GU / ED filters belong to B; they do not discard prefixes that can still extend.'],18)
s.text(48,1305,'d, t and B are the initiation, transition and complete-boundary factors supplied by Weights.',17)
s.parts.append('</g>')
s.footer('SeededPartitionFunction.cpp: accepts, dest, nonstack, h, total; StackSeedDomain::seedEndingLength()', 'Code states: Uℓ ↔ c=ℓ−1; S₁ ↔ c=cap; S₂ ↔ c=cap+1. Checked nonnegative arithmetic throughout.');s.save()

# 6. Reverse differentiation and regional normalization.
s=SVG('IntaRNA.SeededPartitionFunction.outside.svg','Actual base-pair probabilities · outside recurrence','Reverse differentiation of the suffix inside computation; same fixed p and same boundary weights B',1610)
s.lines(48,185,['H and T are forward state/aggregate weights. Their adjoints are A=∂Z/∂H and G=∂Z/∂T.',
                'M(q) is the raw mass of complete chains containing actual pair q; it is not site-coverage mass.'])
s.rna(110,275,710,[(0,'p'),(.4,'r'),(.65,'q'),(1,'terminal')],shade=(0,.4,'H(r,c)'))
s.line(571.5,275,571.5,335,ORANGE,5)
s.lines(920,288,['Reverse the incoming edge r → q.', 'Only actual pair q receives this edge’s mass.'],18)
s.parts.append('<g transform="translate(0 140)">')
s.box(48,248,1344,139,'1  Initialize the reverse objective')
s.eq(73,313,'A(q,S₂)=B(p,q);   A(q,S₁)=B(p,q) only without noLP;   other A and G = 0',23)
s.text(73,354,'Initialize at reachable positive seeded totals; a zero terminal coefficient contributes nothing.',18)
s.box(48,423,1344,197,'2  Reverse q in decreasing order; first scatter aggregate adjoints',fill='#edf6f2',color=TEAL)
s.eq(73,489,'A(q,Uℓ) += G₀(q)   for the same departure-eligible ℓ as in the forward sum',24)
s.eq(73,534,'A(q,S₂) += G₁(q);   A(q,S₁) += G₁(q) only without noLP',24)
s.text(73,576,'This scatter must precede reversing the incoming edges at q.',18)
s.box(48,659,645,274,'3a  Reverse a stack edge r→q')
s.eq(73,724,'a = A(q,D(c,stack)) · t(r,q)',24)
s.eq(73,770,'A(r,c) += a',27)
s.eq(73,816,'M(q) += a · H(r,c)',27)
s.lines(73,860,['Repeat for every forward source state c.',
                'The transition marker belongs to destination q.'],18)
s.box(720,659,672,274,'3b  Reverse a nonstack edge r→q')
s.eq(743,724,'a = A(q,D(b,loop)) · t(r,q)',24)
s.eq(743,770,'G_b(r) += a',27)
s.eq(743,816,'M(q) += a · T_b(r)',27)
s.lines(743,860,['Repeat for b=0 and b=1 using the same forward totals.',
                 'Do not mark the unpaired nucleotides in the loop.'],18)
s.box(48,970,1344,130,'4  Mark the source pair after the reverse sweep',fill='#fff7ed',color=ORANGE)
s.eq(73,1043,'M(p) += A(p,c₀) · d(p)',28)
s.text(650,1043,'Every other pair was marked when its incoming edge was reversed.',18)
s.box(48,1139,1344,179,'5  Merge successful disjoint search regions before normalization')
s.eq(73,1212,'Z_total = ∑ Z_region;    M_total(i,j) = ∑ M_region(i,j);    P(i,j)=M_total(i,j)/Z_total',23)
s.lines(73,1255,['Map regional pairs with energy.getBasePair() into original target/query coordinates, then finalize explicitly.',
                 'Z_total=0 → NA. Failed/cancelled computation → no matrix. No unbound-state weight is added.'],18)
s.text(48,1356,'Seed admission and B are fixed while differentiating. Masses sum across every eligible left boundary p.',18)
s.parts.append('</g>')
s.footer('SeededPartitionFunction.cpp: dh (A), dt (G), mass and initial-pair contribution', 'PredictorMfeEns2dSeedExtension::predictStackSeeds(); BasePairProbabilities::addRegion(), finalize(), probabilities()');s.save()

# 7. Legacy extension kernels remain shared with the heuristic seeded predictor.
s=SVG('IntaRNA.PredictorMfeEns2dSeedExtension.legacy.svg','Legacy seeded extensions · anchor factorization','Bulged-seed fallback in mode M; also the extension kernels used by the seeded mode H predictor',1410)
s.lines(48,185,['For one retained anchor seed a…b, use S=BW(E_seed(a)), d=BW(E_init), t(i,k)=BW(E_interLeft(i,k)).',
                'Lₐ(i) extends left to a (includes initiation); R_b(j) extends right from b (excludes initiation).'])
s.rna(120,278,1180,[(0,'i'),(.33,'a'),(.4,''),(.5,''),(.59,''),(.66,'b'),(1,'j')],shade=(0,.33,'Lₐ(i)'),seed=(.33,.66))
s.parts.append('<rect x="899" y="278" width="401" height="60" fill="#dcebf5"/>')
s.rna(899,278,401,[(0,''),(1,'')],shade=(0,1,'R_b(j)'))
s.eq(530,379,'hybrid contribution = Lₐ(i) · S · R_b(j)',25)
s.box(48,421,1344,178,'Lonely pairs allowed')
s.eq(73,488,'Lₐ(a)=d;   Lₐ(i)=∑ t(i,k)Lₐ(k) − Cₐ(i)     (i ≺ k ≼ a)',25)
s.eq(73,538,'R_b(b)=1;   R_b(j)=∑ R_b(k)t(k,j)              (b ≼ k ≺ j)',25)
s.text(73,577,'Only nonzero, complementary states and loop-size-admissible terms participate. Sums use componentwise bounds.',17)
s.box(48,638,1344,245,'noLP: outward boundary stacks are explicit',fill='#edf6f2',color=TEAL)
s.eq(73,707,'Lₐ(i)=t(i,i⁺) · [Lₐ(i⁺) + ∑ t(i⁺,k)Lₐ(k)] − Cₐ(i)    (i⁺ ≺ k ≼ a)',24)
s.eq(73,762,'R_b(j)=t(j⁻,j) · [1{j⁻=b}R_b(b) + ∑ R_b(k)t(k,j⁻)]    (b ≼ k ≺ j⁻)',24)
s.lines(73,808,['i⁺=i+(1,1), j⁻=j−(1,1); require the extra stacked pair. Anchor bases still have Lₐ(a)=d and R_b(b)=1.',
                'These are the legacy loop bounds. The exact stack-seed path instead uses the suffix-state recurrence.'],18)
s.box(48,923,1344,212,'Left-seed correction Cₐ(i)',fill='#fff7ed',color=ORANGE)
s.lines(73,988,['At a competing seed start i: process it if its end is left of a, or it loop-overlaps the anchor.',
                'If a leftmost loop-overlapping seed h is found, subtract BW(getNonOverlappingEnergy(i,h)) · Lₐ(h).',
                'Otherwise subtract its seed weight times the continuation at its end; noLP adds explicit loop continuations.',
                'Apply the implementation’s overlap/loop guards, then clamp a negative corrected L value to zero.'],18,step=30)
s.eq(48,1185,'updateZ(i,j, Lₐ(i) · S · R_b(j), true)    for each nonzero span-admissible combination',23)
s.lines(48,1231,['Complete-site filters and B(i,j) are applied by updateZ. The anchor’s hybrid seed energy is charged once.',
                 'This legacy subtraction-based path has no actual-pair outside pass; do not substitute it for the exact seeded union.'],18)
s.footer('PredictorMfeEns2dSeedExtension.cpp: predict() legacy branch, fillHybridZ_left(), fillHybridZ_right()', 'PredictorMfeEns2dSeedExtension.cpp: getNonOverlappingEnergy(); PredictorMfeEns.cpp: updateZ()');s.save()

# 8. Restriction to one optimum right extension per anchor seed.
s=SVG('IntaRNA.PredictorMfeEns2dHeuristicSeedExtension.svg','Seeded ensemble · one selected right extension','PredictorMfeEns2dHeuristicSeedExtension  |  --model=P --mode=H',1200)
s.lines(48,185,['For anchor seed a…b: S=BW(E_seed), d=BW(E_init), L=Lₐ(i), R=R_b(j).',
                'L and R use the inherited legacy extension recurrences, even when the seed handler is stack-only.'])
s.box(48,249,1344,198,'1  Compute R, submit true right extensions, and select j*')
s.eq(73,317,'updateZ(a,j, d · S · R_b(j), true)     for j ≻ b',25)
s.eq(73,366,'j* = argmin E_full(a,j, d · S · R_b(j))     over right-end-admissible candidates',24)
s.text(73,411,'updateOptRightZ compares complete energy and the right GU-end restriction; ties retain the first candidate.',18)
s.line(720,458,720,504,arrow=True,width=2)
s.box(48,524,1344,288,'2  Compute corrected L and submit left-side combinations',fill='#edf6f2',color=TEAL)
s.rna(92,621,495,[(0,'i'),(.33,'a'),(.45,''),(.60,'b'),(1,'j*')],shade=(0,.33,'L'),seed=(.33,.6))
s.eq(650,603,'updateZ(i,b, Lₐ(i) · S, true)',25)
s.text(650,643,'Includes the seed-only case i=a.',18)
s.eq(650,695,'updateZ(i,j*, Lₐ(i) · S · R_b(j*), true)',23)
s.lines(650,736,['Only for a true left extension and a true selected right extension,',
                 'with both complete spans within their limits.'],17)
s.box(48,852,1344,185,'Candidate family per anchor',fill='#fff7ed',color=ORANGE)
s.lines(73,916,['All true right-only extensions + all corrected left-only extensions + left extensions paired with just j*.',
                'The second phase runs if noGUend is off, or a valid right optimum exists, or anchor end b is not GU.',
                'Initialize j*=b and E_right_opt=∞. A positive selected R is required before combining both sides.',
                'No full Cartesian product of left and right boundaries is retained, and bpProb / bpsvg are unsupported.'],18)
s.text(48,1082,'Every update is filtered and accumulated by PredictorMfeEns::updateZ; ranked reporting is a separate step.',18)
s.footer('PredictorMfeEns2dHeuristicSeedExtension.cpp: predict(), fillHybridZ_right(), fillHybridZ_left()', 'PredictorMfeEns2dHeuristicSeedExtension.h: updateOptRightZ(); inherited legacy fillHybridZ kernels');s.save()

Path(__file__).with_name('new-svg-files.json').write_text(json.dumps(records,indent=2)+'\n')
print(f'Wrote {len(records)} standalone SVG recurrence depictions')
