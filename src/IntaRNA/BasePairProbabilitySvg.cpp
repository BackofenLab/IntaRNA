#include "IntaRNA/BasePairProbabilityWriter.h"
#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/Interaction.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace IntaRNA {
namespace {

// One shared scale for matrix cells, accessibility frames and legend swatches.
// Fill colors live exclusively in the SVG stylesheet, allowing later editing.
constexpr std::array<double,8> breaks{0,0.01,0.1,0.25,0.5,0.75,0.9,1};
constexpr double cell=20;

std::string xml(const std::string & value)
{
	std::string escaped;
	for (unsigned char c:value) {
		switch(c) {
			case '&': escaped+="&amp;"; break;
			case '<': escaped+="&lt;"; break;
			case '>': escaped+="&gt;"; break;
			case '"': escaped+="&quot;"; break;
			case '\'': escaped+="&apos;"; break;
			default:
				if (c<32 && c!='\t' && c!='\n' && c!='\r')
					throw std::invalid_argument("bpsvg: sequence name contains an invalid XML character");
				escaped+=c;
		}
	}
	return escaped;
}

size_t colorClass(Z_type probability)
{
	for(size_t k=1;k<breaks.size()-1;++k)
		if(probability<breaks[k]) return k-1;
	return breaks.size()-2;
}

Z_type unpaired(const Accessibility & accessibility,size_t i,Z_type RT)
{
	const E_type ed=accessibility.getED(i,i);
	const Z_type p=E_isINF(ed)?Z_type(0):Z_exp(-E_2_Z(ed)/RT);
	if (!std::isfinite(p) || p<0 || p>1)
		throw std::runtime_error("bpsvg: invalid single-nucleotide unpaired probability");
	return p;
}

} // namespace

std::string
BasePairProbabilityWriter::svgBlock(const BasePairProbabilities & result,const InteractionEnergy & energy,
		Z_type targetRT,Z_type queryRT,const Interaction * mfe)
{
	if (!std::isfinite(targetRT) || targetRT<=0 || !std::isfinite(queryRT) || queryRT<=0)
		throw std::invalid_argument("bpsvg: accessibility energy scales must be positive and finite");
	const auto & at=energy.getAccessibility1();
	const auto & aq=energy.getAccessibility2().getAccessibilityOrigin();
	const auto & t=at.getSequence();
	const auto & q=aq.getSequence();
	const bool empty=result.status()==BasePairProbabilities::Status::empty;
	if (!empty && result.status()!=BasePairProbabilities::Status::nonempty)
		throw std::logic_error("bpsvg output requires successful explicit finalization");
	if (result.rawMasses().size1()!=t.size() || result.rawMasses().size2()!=q.size())
		throw std::invalid_argument("bpsvg output sequence dimensions do not match result");
	const auto p=empty?Matrix<Z_type>():result.probabilities();
	const double w=cell*q.size(), h=cell*t.size();
	// Leave space for names even for very short sequences; keep a fixed nt pitch.
	const double nameWidth=std::max(9.0*q.getId().size()+80,7.0*(t.getId().size()+q.getId().size())+100);
	const double x=std::max(130.0,(std::max(820.0,nameWidth)-w)/2);
	const double y=std::max(120.0,(9.0*t.getId().size()-h)/2+60);
	const double width=2*x+w, height=2*y+h+185;
	// SVG y increases downwards, while target coordinates increase upwards.
	auto targetY=[&](size_t i) { return y+cell*(t.size()-1-i); };
	const bool markMfe=mfe && !mfe->isEmpty();
	size_t tStart=t.size(),tEnd=0,qStart=q.size(),qEnd=0;
	if (markMfe) {
		for (const auto & bp:mfe->basePairs) {
			if (bp.first>=t.size() || bp.second>=q.size())
				throw std::invalid_argument("bpsvg: MFE interaction coordinates exceed sequence dimensions");
			tStart=std::min(tStart,bp.first); tEnd=std::max(tEnd,bp.first);
			qStart=std::min(qStart,bp.second); qEnd=std::max(qEnd,bp.second);
		}
	}
	std::ostringstream out;
	out.imbue(std::locale::classic());
	out<<std::setprecision(std::numeric_limits<Z_type>::max_digits10);
	out<<"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""<<width<<"\" height=\""<<height
		<<"\" viewBox=\"0 0 "<<width<<' '<<height<<"\" role=\"img\">\n"
		<<"<title>Base-pair probabilities: "<<xml(t.getId())<<" / "<<xml(q.getId())<<"</title>\n"
		<<"<desc>Target rows run 5' to 3' upwards and query columns run 5' to 3' to the right. Pair probabilities are conditional on an allowed interaction. "
		<<"The frame shows single-nucleotide intramolecular pairing probabilities (1 - Pu) from the accessibility model.</desc>\n"
		<<R"SVG(<style>
text { font-family: sans-serif; fill: #243247; }
.heading { font-size: 18px; font-weight: bold; }
.axis { font-size: 15px; }
.subtitle { font-size: 13px; }
.index, .legend { font-size: 11px; }
.nt { font-family: monospace; font-size: 13px; pointer-events: none; }
.nt-light { fill: #fff; }
/* Segmented probability gradient: edit these seven rules to recolor the plot. */
.p0 { fill: #FFFFFF; }
.p1 { fill: #DEEAF0; }
.p2 { fill: #AEC6CF; }
.p3 { fill: #78A6C2; }
.p4 { fill: #4682B4; }
.p5 { fill: #24558A; }
.p6 { fill: #002060; }
.na { fill: #e5e5e5; }
.guide { stroke: #b0b0b0; stroke-width: 0.6; pointer-events: none; }
.major { stroke: #777; stroke-width: 1.2; }
.origin { stroke: #888; stroke-width: 1.2; stroke-dasharray: 3 2; }
.border { fill: none; stroke: #aaa; stroke-width: 0.7; pointer-events: none; }
.seed { fill: none; stroke: #c65d00; stroke-width: 1.8; pointer-events: none; }
.mfe { fill: none; stroke: #16823b; stroke-width: 2.5; pointer-events: none; }
</style>
)SVG";
	out<<"<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
	auto text=[&](double tx,double ty,const std::string & value,const std::string & cls,const char * anchor="middle") {
		out<<"<text class=\""<<cls<<"\" x=\""<<tx<<"\" y=\""<<ty<<"\" text-anchor=\""<<anchor<<"\">"<<xml(value)<<"</text>\n";
	};
	text(width/2,30,"Base-pair probabilities","heading");
	if (t.getId()!="target" || q.getId()!="query")
		text(width/2,53,t.getId()+" / "+q.getId(),"subtitle");
	text(x+w/2,y+h+65,q.getId()+" (5' to 3')","axis");
	out<<"<text class=\"axis\" transform=\"translate("<<x-90<<' '<<y+h/2
		<<") rotate(-90)\" text-anchor=\"middle\">"<<xml(t.getId())<<" (5' to 3')</text>\n";
	auto rect=[&](double rx,double ry,Z_type value,bool missing,const std::string & attributes,const std::string & title) {
		out<<"<rect class=\"probability "<<(missing?"na":"p"+std::to_string(colorClass(value)))
			<<"\" x=\""<<rx<<"\" y=\""<<ry<<"\" width=\"20\" height=\"20\" "<<attributes
			<<" data-probability=\"";
		if (missing) out<<"NA"; else out<<value;
		out<<"\"><title>"<<xml(title)<<": ";
		if (missing) out<<"NA (empty interaction ensemble)";
		else {
			if (value<0.001) out<<std::scientific<<std::setprecision(2);
			else out<<std::fixed<<std::setprecision(3);
			out<<value<<std::defaultfloat<<std::setprecision(std::numeric_limits<Z_type>::max_digits10);
		}
		out<<"</title></rect>\n";
	};
	out<<"<g class=\"matrix\">\n";
	for(size_t i=0;i<t.size();++i) for(size_t j=0;j<q.size();++j) {
		const auto ti=std::to_string(t.getInOutIndex(i)), qj=std::to_string(q.getInOutIndex(j));
		rect(x+cell*j,targetY(i),empty?0:p(i,j),empty,
			"data-type=\"base-pair\" data-target=\""+ti+"\" data-query=\""+qj+"\"",
			"Base-pair probability ("+t.getId()+" "+ti+", "+q.getId()+" "+qj+")");
	}
	out<<"</g>\n<g class=\"accessibility-frame\">\n";
	for(bool target:{false,true}) {
		const auto & seq=target?t:q;
		const auto & acc=target?at:aq;
		for(size_t i=0;i<seq.size();++i) {
			const Z_type value=1-unpaired(acc,i,target?targetRT:queryRT);
			const auto idx=std::to_string(seq.getInOutIndex(i));
			for(bool far:{false,true}) {
				const double rx=target?(far?x+w:x-cell):x+cell*i;
				const double ry=target?targetY(i):(far?y+h:y-cell);
				rect(rx,ry,value,false,"data-type=\"intramolecular-pairing\" data-strand=\""+std::string(target?"target":"query")
					+"\" data-index=\""+idx+"\"","Intramolecular pairing probability ("+seq.getId()+" "+idx+")");
				text(rx+cell/2,ry+14,seq.asString().substr(i,1),value>=.75?"nt nt-light":"nt");
			}
		}
	}
	out<<"</g>\n<g class=\"guides\">\n";
	for(bool target:{false,true}) {
		const auto & seq=target?t:q;
		for(size_t i=0;i<seq.size();++i) {
			const auto idx=seq.getInOutIndex(i);
			const bool origin=idx==1 && i>0 && seq.getInOutIndex(i-1)==-1;
			if (idx%10!=0 && idx!=1) continue;
			// Negative multiples: before the labelled nt (-11|-10); positive:
			// after it (10|11). At the missing zero: -1|+1.
			const double offset=cell*(i+(idx>0 && !origin?1:0));
			if (idx%10==0 || origin) out<<"<line class=\"guide"<<(origin?" origin":idx%50==0?" major":"")
				<<"\" data-strand=\""<<(target?"target":"query")<<"\" data-index=\""<<idx
				<<"\" x1=\""<<(target?x-cell:x+offset)<<"\" y1=\""<<(target?y+h-offset:y-cell)
				<<"\" x2=\""<<(target?x+w+cell:x+offset)<<"\" y2=\""<<(target?y+h-offset:y+h+cell)<<"\"/>\n";
			const std::string label=idx==1?"+1":std::to_string(idx);
			if(target) {
				text(x-26,targetY(i)+cell/2+4,label,"index","end");
				text(x+w+26,targetY(i)+cell/2+4,label,"index","start");
			} else {
				text(x+cell*(i+.5),y-28,label,"index");
				text(x+cell*(i+.5),y+h+38,label,"index");
			}
		}
	}
	out<<"</g>\n<rect class=\"border\" x=\""<<x<<"\" y=\""<<y<<"\" width=\""<<w<<"\" height=\""<<h<<"\"/>\n";
	bool annotated=false;
	if(result.collectsSeedPairs()) {
		out<<"<g class=\"seed-outlines\">\n";
		for(size_t i=0;i<t.size();++i) for(size_t j=0;j<q.size();++j) if(result.seedPairs()(i,j)) {
			annotated=true;
			out<<"<rect class=\"seed\" data-target=\""<<t.getInOutIndex(i)<<"\" data-query=\""<<q.getInOutIndex(j)
				<<"\" x=\""<<x+cell*j+1<<"\" y=\""<<targetY(i)+1<<"\" width=\"18\" height=\"18\"/>\n";
		}
		out<<"</g>\n";
	}
	if (markMfe) {
		out<<"<rect class=\"mfe\" data-target-start=\""<<t.getInOutIndex(tStart)
			<<"\" data-target-end=\""<<t.getInOutIndex(tEnd)
			<<"\" data-query-start=\""<<q.getInOutIndex(qStart)
			<<"\" data-query-end=\""<<q.getInOutIndex(qEnd)
			<<"\" x=\""<<x+cell*qStart<<"\" y=\""<<targetY(tEnd)
			<<"\" width=\""<<cell*(qEnd-qStart+1)<<"\" height=\""<<cell*(tEnd-tStart+1)<<"\"/>\n";
	}
	const double legendX=(width-7*86)/2, legendY=y+h+110;
	text(width/2,legendY-14,"Probability ranges (intermolecular base pairs and intramolecular pairing)","legend");
	out<<"<g class=\"probability-legend\">\n";
	for(size_t k=0;k<breaks.size()-1;++k) {
		out<<"<rect class=\"p"<<k<<"\" x=\""<<legendX+k*86<<"\" y=\""<<legendY<<"\" width=\"86\" height=\"16\"/>\n";
		std::ostringstream range;
		range.imbue(std::locale::classic());
		range<<'['<<breaks[k]<<", "<<breaks[k+1]<<(k+2==breaks.size()?"]":")");
		text(legendX+(k+.5)*86,legendY+32,range.str(),"legend");
	}
	out<<"</g>\n";
	if (annotated) {
		out<<"<rect class=\"seed\" x=\""<<legendX<<"\" y=\""<<legendY+48<<"\" width=\"18\" height=\"18\"/>\n";
		text(legendX+28,legendY+61,"Orange outline: pair in an admitted seed within a searched region","legend","start");
	}
	if (markMfe) {
		const double mfeLegendY=legendY+(annotated?76:48);
		out<<"<rect class=\"mfe\" x=\""<<legendX<<"\" y=\""<<mfeLegendY<<"\" width=\"18\" height=\"18\"/>\n";
		text(legendX+28,mfeLegendY+13,"Green outline: minimum-free-energy interaction region","legend","start");
	}
	if(empty) text(width/2,legendY+90,"Gray / NA: empty interaction ensemble; accessibility remains defined","legend");
	out<<"</svg>\n";
	// Transfer the potentially large document instead of copying its buffer.
	return std::move(out).str();
}

void BasePairProbabilityWriter::writeSvg(std::ostream & out,const BasePairProbabilities & result,const InteractionEnergy & energy,
		Z_type targetRT,Z_type queryRT,const Interaction * mfe)
{
	emit(out,svgBlock(result,energy,targetRT,queryRT,mfe));
}

void BasePairProbabilityWriter::writeSvgFile(const std::string & filename,const BasePairProbabilities & result,const InteractionEnergy & energy,
		Z_type targetRT,Z_type queryRT,const Interaction * mfe)
{
	emitFile(filename,svgBlock(result,energy,targetRT,queryRT,mfe));
}

} // namespace IntaRNA
