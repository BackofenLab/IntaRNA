#!/usr/bin/env python3
"""Create Martin's concise briefing and an unsent email with verified attachments."""
import hashlib,json,shutil
from email.message import EmailMessage
from email.policy import SMTP
from pathlib import Path
from xml.sax.saxutils import escape
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,KeepTogether

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT.parent/'output/pdf';OUT.mkdir(parents=True,exist_ok=True)
DELIVERY=ROOT/'martin';DELIVERY.mkdir(exist_ok=True)
DATA={r['case']:r for r in json.loads((ROOT/'results/summary.json').read_text())}
audit=json.loads((ROOT/'results/audit.json').read_text())
assert audit['status']=='pass' and audit['raw_records']==224
assert all(r['output_equal'] for r in json.loads((ROOT/'profiles/validated-outputs.json').read_text())['checks'])
probe=json.loads((ROOT/'results/stream-cleanup-probe.json').read_text());assert probe['fd_growth']==100
def seconds(n):return DATA[n]['median_s']
def rss(n):return DATA[n]['median_rss_mib']
blue=colors.HexColor('#174D68'); ink=colors.HexColor('#1C2933');muted=colors.HexColor('#52636D')
light=colors.HexColor('#F0F5F7');line=colors.HexColor('#D5DFE4');amber=colors.HexColor('#8B451B')
fonts=Path('/usr/share/fonts/truetype/dejavu')
for name,file in [('Brief','DejaVuSans.ttf'),('Brief-Bold','DejaVuSans-Bold.ttf'),('Brief-Italic','DejaVuSans-Oblique.ttf')]:
    pdfmetrics.registerFont(TTFont(name,str(fonts/file)))
pdfmetrics.registerFontFamily('Brief',normal='Brief',bold='Brief-Bold',italic='Brief-Italic',boldItalic='Brief-Bold')
styles={
 'title':ParagraphStyle('title',fontName='Brief-Bold',fontSize=24,leading=28,textColor=blue,spaceAfter=7),
 'kicker':ParagraphStyle('kicker',fontName='Brief-Bold',fontSize=8,leading=11,textColor=muted,spaceAfter=8),
 'deck':ParagraphStyle('deck',fontName='Brief',fontSize=11,leading=15,textColor=ink,spaceAfter=14),
 'h2':ParagraphStyle('h2',fontName='Brief-Bold',fontSize=12,leading=16,textColor=blue,spaceBefore=9,spaceAfter=6),
 'body':ParagraphStyle('body',fontName='Brief',fontSize=9.2,leading=13,textColor=ink,spaceAfter=5),
 'small':ParagraphStyle('small',fontName='Brief',fontSize=8,leading=11,textColor=muted,spaceAfter=5),
 'cell':ParagraphStyle('cell',fontName='Brief',fontSize=8.5,leading=11.7,textColor=ink),
 'head':ParagraphStyle('head',fontName='Brief-Bold',fontSize=8.5,leading=11.7,textColor=colors.white),
}
def P(text,style='body'):return Paragraph(text,styles[style])
W=A4[0]-88
def grid(rows,widths):
    table=Table([[P(str(v),'head' if i==0 else 'cell') for v in row] for i,row in enumerate(rows)],colWidths=widths,repeatRows=1,hAlign='LEFT')
    table.setStyle(TableStyle([
      ('BACKGROUND',(0,0),(-1,0),blue),('VALIGN',(0,0),(-1,-1),'TOP'),
      ('LEFTPADDING',(0,0),(-1,-1),8),('RIGHTPADDING',(0,0),(-1,-1),8),
      ('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),
      ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,light]),
      ('LINEBELOW',(0,-1),(-1,-1),.5,line)]))
    return table
PDF=OUT/'IntaRNA-profiling-brief-for-Martin.pdf'
doc=SimpleDocTemplate(str(PDF),pagesize=A4,rightMargin=44,leftMargin=44,topMargin=38,bottomMargin=44,
    title='IntaRNA profiling: findings and proposed next steps',author='IntaRNA profiling',subject='Technical briefing for Martin, revision 1b83779')
story=[]
story += [P('ENGINEERING BRIEF  /  19 SEPTEMBER 2026','kicker'),P('IntaRNA profiling','title'),
          P('Findings and proposed next steps for Martin','deck'),
          P('<b>Tested revision:</b> refactor-to-c++23 / 1b8377989ba7. This is a profile of one recorded build; no comparison with an earlier release or another PR was performed.','small')]
story += [P('Recommendation','h2'),P('Fix the confirmed input-stream cleanup defect first. Then review the accessibility-file limit behavior and benchmark narrowly scoped CPU optimizations, separately for default, exact and ensemble modes.')]
story += [P('Confirmed defect: input streams are not cleaned up','h2'),
          P('<b>Reproduction:</b> 100 open/read/cleanup cycles against the freshly built library increased open file descriptors from <b>4 to 104</b>. All four heap traces also retained 9,384 bytes in input-stream allocation stacks.'),
          P('<b>Cause:</b> <font name="Brief-Bold">src/IntaRNA/general.cpp:136</font> casts the input stream to <b>filtering_ostream</b>, although it was created as <b>filtering_istream</b>. The cast fails and cleanup is skipped. Resources accumulate in repeated library use; the OS reclaims them when a CLI process exits. Introduction history was not investigated.','body')]
story += [P('Performance findings worth acting on','h2')]
rows=[['Finding','Measured evidence','Interpretation / boundary'],
 ['Default X/H CPU', 'Loop-energy evaluation: <b>57.6%</b> of 368 CPU samples in the biological GcvB batch.', 'First optimization target: repeated internal-loop energy work. Sample share is not a promised speedup.'],
 ['Batch parallelism', f'48 targets × ChiX: <b>{seconds("throughput_1t"):.3f} → {seconds("throughput_12t"):.3f} s</b> at 1 → 12 threads (<b>5.08×</b>).', 'Peak RSS: 18.7 → 65.3 MiB. Batch throughput; not single-pair acceleration.'],
 ['Accessibility reuse', '<b>0.805 → 0.492 s</b> (39% lower latency); same base-pair list.', 'Matched target/query caps: 149/99. One-time generation: 0.847 s, excluded from cached-call time.'],
 ['Windowing, S/H', '<b>127.5 → 24.1 MiB</b> (-81%); <b>21.38 → 59.71 s</b> (2.79× time).', '3,000 × 800 nt; window 300, overlap 150. Default X/H gets no useful RSS reduction and takes 1.72× longer.']]
story.append(grid(rows,[92,206,W-298]))
story += [P('Proposed sequence','h2'),
 P('<b>1.</b> Isolate the cleanup fix and add repeated plain/gzip stream tests. <b>2.</b> Confirm intended ED-file limit semantics with Martin. <b>3.</b> Agree on a production workload, then measure an energy/scoring optimization against this baseline. Proposed changes have not been implemented or assigned an expected speedup.')]
story.append(PageBreak())
story += [P('EVIDENCE AND INTERPRETATION','kicker'),P('How to read the results','title')]
story += [P('Reproducible baseline','h2'),
 P('<b>56 configurations:</b> 168 measured repetitions plus 56 excluded warmups (224 CLI runs). Three measured repetitions per case, deterministic randomized order, sequential execution. Medians are shown; full ranges and raw records are supplied. CPU sampling, heap tracing and diagnostic logging were separate from clean timings.'),
 P('<b>Host/build:</b> Ryzen 5 7530U (6 cores / 12 hardware threads), Linux 6.8.0-139; GCC 16.1.0, C++23, -O3 with symbols/frame pointers; ViennaRNA 2.7.2, Boost 1.85.0. OPENBLAS_NUM_THREADS=1, OMP_DYNAMIC=FALSE, OMP_PLACES=cores and OMP_PROC_BIND=close. Exact commands, hashes, dependencies and input sequences accompany the report.')]
story += [P('Batch scaling: choosing a thread count','h2')]
rows=[['Threads','Median wall s','Observed range s','Peak RSS MiB','Speedup']]
for n in [1,2,4,6,12]:
    r=DATA[f'throughput_{n}t'];rows.append([n,f"{r['median_s']:.3f}",f"{r['min_s']:.3f}–{r['max_s']:.3f}",f"{r['median_rss_mib']:.1f}",f"{seconds('throughput_1t')/r['median_s']:.2f}×"])
story.append(grid(rows,[55,94,134,124,W-407]))
story += [P('48 evaluated 300-nt targets × ChiX (84 nt); 46 reportable interactions. Six threads supply 88% of twelve-thread throughput at 61% of its RSS. Peak RSS values are medians of per-process high-water marks.','small')]
story += [P('Two interpretation traps','h2'),
 P('<b>Cached ED is not automatically a transparent replacement.</b> The writer emits columns through its maximum; the reader advertises one nucleotide less for dangling-end treatment (AccessibilityFromStream.cpp:90). We therefore compared computed and loaded ED at the same 149/99 limits. Whether this is intended needs review before generalizing the saving.'),
 P('<b>Mode changes are not implementation speedups.</b> Exact X/M is dominated by repeated scoring terms (Boltzmann weights 19.4%, ED lookup 18.5% of samples); ensemble mode adds partition-function recurrence and map costs. Disabling accessibility changes the scientific model and was slower in one test. Windowing and algorithm substitutions need their own output checks.')]
story += [P('Validation and uncertainty','h2'),
 P('All timed cases and the automated evidence audit passed; 4,281 assertions in 36 API tests and all 20 CLI regression cases passed. Base-pair lists matched for 1 vs 12 threads and the matched ED comparison. Window comparisons matched coordinates and energy; full structure equality was not separately checked there.'),
 P('This is one laptop, three measured repetitions per case and a bounded workload suite. No release-wide regression or biological-accuracy claim follows. CPU profiles total 1,182 samples at 100 ms. A collector timer warning and incomplete optimized stack ancestry limit precision, although CPU totals were consistent. Hardware counters were blocked; no cache-miss or IPC claims are made.')]
story += [P('Evidence and reporting conventions','h2'),
 P('Attached evidence archive: start with <b>profiling/REPORT.md</b>. Inspect <b>results/summary.csv</b>, <b>results/audit.json</b>, <b>results/stream-cleanup-probe.json</b>, <b>cases.json</b>, and the reproduction scripts. Original data and hashes are retained.','small'),
 P('Reporting follows full configuration disclosure and explicit warmup/repetition statistics, as illustrated by <link href="https://www.spec.org/cpu2017/Docs/runrules.html" color="#174D68">SPEC reporting rules</link> and the <link href="https://github.com/google/benchmark/blob/main/docs/user_guide.md" color="#174D68">Google Benchmark guide</link>. This custom IntaRNA study does not claim SPEC certification.','small')]
def footer(c,d):
    c.saveState();c.setFont('Brief',7.5);c.setFillColor(muted)
    c.drawString(44,25,'IntaRNA  |  profiling date 2026-09-19  |  revision 1b83779')
    c.drawRightString(A4[0]-44,25,f'{d.page}');c.restoreState()
doc.build(story,onFirstPage=footer,onLaterPages=footer)

subject='IntaRNA profiling (1b83779): confirmed cleanup defect and optimization priorities'
body='''Hi Martin,

We completed a profiling pass on refactor-to-c++23 at commit 1b83779: 56 configurations, 168 measured repetitions plus 56 warmups, and separate CPU and heap profiles. The attached two-page brief gives the findings and proposed next steps; the evidence archive contains commands, inputs, raw results, audits and reproduction scripts.

The main findings are:

- Confirmed input-stream cleanup defect. Against the freshly built library, 100 open/read/cleanup calls increased open file descriptors from 4 to 104. general.cpp:136 checks for filtering_ostream although newInputStream creates filtering_istream. This accumulates resources in repeated library use; the OS reclaims them when a CLI process exits. I suggest treating this as the first, isolated fix.
- Default-mode optimization target. Internal-loop energy evaluation accounted for about 58% of CPU samples in the biological GcvB batch. Exact and ensemble modes have different scoring/partition-function hotspots, so optimization work should be measured separately by mode.
- Measured operating tradeoffs. The 48-target ChiX batch fell from 4.944 s on one thread to 0.974 s on twelve (5.08x), with RSS increasing from 18.7 to 65.3 MiB. In S mode, windowing reduced RSS by 81% but took 2.79x longer; default X mode showed no useful memory saving from windowing on that input.
- ED reuse requires a semantic check. With explicitly matched target/query limits of 149/99, loading ED reduced latency from 0.805 to 0.492 s and preserved the base-pair list. The reader currently lowers the available maximum length by one, so the saving should not be generalized to an unqualified round trip.

All regression tests and the evidence audit passed. These are measurements of one build on one laptop, with three measured repetitions per case; they do not establish a cross-version improvement or a universal speedup. Hardware counters were unavailable, and the sampling limitations are documented.

I propose an isolated cleanup fix first, then a review of the ED-file limit behavior, followed by a focused CPU optimization against an agreed production workload. Could you confirm whether the ED limit reduction is intended and which workload should be our acceptance benchmark?

Best,
Alex
'''
(DELIVERY/'email-to-Martin.md').write_text('Subject: '+subject+'\n\n'+body)
message=EmailMessage(policy=SMTP)
message['Subject']=subject;message['X-Unsent']='1'
message.set_content(body)
message.add_attachment(PDF.read_bytes(),maintype='application',subtype='pdf',filename=PDF.name)
archive=ROOT/'intarna-profile-evidence.zip'
message.add_attachment(archive.read_bytes(),maintype='application',subtype='zip',filename=archive.name)
(DELIVERY/'email-to-Martin.eml').write_bytes(message.as_bytes())
(DELIVERY/'delivery-status.json').write_text(json.dumps(dict(status='prepared_not_sent',
    recipient=None,channel=None,subject=subject,pdf=str(PDF),draft=str(DELIVERY/'email-to-Martin.eml'),
    attachments=[dict(name=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [PDF,archive]],
    reason='Recipient and delivery channel require confirmation; no email integration is connected.'),indent=2)+'\n')
print(PDF)
print(DELIVERY/'email-to-Martin.eml')
