#include "IntaRNA/BasePairProbabilityWriter.h"
#include "IntaRNA/InteractionEnergy.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>

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

Z_type unpaired(const Accessibility & accessibility,size_t i,const InteractionEnergy & energy)
{
	const E_type ed=accessibility.getED(i,i);
	const Z_type p=E_isINF(ed)?Z_type(0):energy.getBoltzmannWeight(ed);
	if (!std::isfinite(p) || p<0 || p>1)
		throw std::runtime_error("bpsvg: invalid single-nucleotide unpaired probability");
	return p;
}

} // namespace

std::string
BasePairProbabilityWriter::svgBlock(const BasePairProbabilities & result,const InteractionEnergy & energy)
{
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
	const double x=std::max(130.0,(std::max(820.0,9.0*q.getId().size()+80)-w)/2);
	const double y=std::max(120.0,(9.0*t.getId().size()-h)/2+60);
	const double width=2*x+w, height=2*y+h+160;
	std::ostringstream out;
	out.imbue(std::locale::classic());
	out<<std::setprecision(std::numeric_limits<Z_type>::max_digits10);
	out<<"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""<<width<<"\" height=\""<<height
		<<"\" viewBox=\"0 0 "<<width<<' '<<height<<"\" role=\"img\">\n"
		<<"<title>Base-pair probabilities: "<<xml(t.getId())<<" / "<<xml(q.getId())<<"</title>\n"
		<<"<desc>Target rows and query columns run 5' to 3'. Pair probabilities are conditional on an allowed interaction. "
		<<"The frame shows single-nucleotide unpaired probabilities from the accessibility model.</desc>\n"
		<<R"SVG(<style>
text { font-family: sans-serif; fill: #243247; }
.heading { font-size: 18px; font-weight: bold; }
.axis { font-size: 15px; }
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
</style>
)SVG";
	out<<"<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
	auto text=[&](double tx,double ty,const std::string & value,const std::string & cls,const char * anchor="middle") {
		out<<"<text class=\""<<cls<<"\" x=\""<<tx<<"\" y=\""<<ty<<"\" text-anchor=\""<<anchor<<"\">"<<xml(value)<<"</text>\n";
	};
	text(width/2,30,"Base-pair probabilities","heading");
	text(x+w/2,y-65,q.getId()+" (5' to 3')","axis");
	out<<"<text class=\"axis\" transform=\"translate("<<x-90<<' '<<y+h/2
		<<") rotate(-90)\" text-anchor=\"middle\">"<<xml(t.getId())<<" (5' to 3')</text>\n";
	auto rect=[&](double rx,double ry,Z_type value,bool missing,const std::string & attributes,const std::string & title) {
		out<<"<rect class=\"probability "<<(missing?"na":"p"+std::to_string(colorClass(value)))
			<<"\" x=\""<<rx<<"\" y=\""<<ry<<"\" width=\"20\" height=\"20\" "<<attributes
			<<" data-probability=\"";
		if (missing) out<<"NA"; else out<<value;
		out<<"\"><title>"<<xml(title)<<": ";
		if (missing) out<<"NA (empty interaction ensemble)"; else out<<value;
		out<<"</title></rect>\n";
	};
	out<<"<g class=\"matrix\">\n";
	for(size_t i=0;i<t.size();++i) for(size_t j=0;j<q.size();++j) {
		const auto ti=std::to_string(t.getInOutIndex(i)), qj=std::to_string(q.getInOutIndex(j));
		rect(x+cell*j,y+cell*i,empty?0:p(i,j),empty,
			"data-type=\"base-pair\" data-target=\""+ti+"\" data-query=\""+qj+"\"",
			"Base-pair probability (target "+ti+", query "+qj+")");
	}
	out<<"</g>\n<g class=\"accessibility-frame\">\n";
	for(bool target:{false,true}) {
		const auto & seq=target?t:q;
		const auto & acc=target?at:aq;
		for(size_t i=0;i<seq.size();++i) {
			const Z_type value=unpaired(acc,i,energy);
			const auto idx=std::to_string(seq.getInOutIndex(i));
			for(bool far:{false,true}) {
				const double rx=target?(far?x+w:x-cell):x+cell*i;
				const double ry=target?y+cell*i:(far?y+h:y-cell);
				rect(rx,ry,value,false,"data-type=\"unpaired\" data-strand=\""+std::string(target?"target":"query")
					+"\" data-index=\""+idx+"\"","Unpaired probability ("+std::string(target?"target ":"query ")+idx+")");
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
				<<"\" x1=\""<<(target?x-cell:x+offset)<<"\" y1=\""<<(target?y+offset:y-cell)
				<<"\" x2=\""<<(target?x+w+cell:x+offset)<<"\" y2=\""<<(target?y+offset:y+h+cell)<<"\"/>\n";
			const std::string label=idx==1?"+1":std::to_string(idx);
			if(target) {
				text(x-26,y+cell*(i+.5)+4,label,"index","end");
				text(x+w+26,y+cell*(i+.5)+4,label,"index","start");
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
				<<"\" x=\""<<x+cell*j+1<<"\" y=\""<<y+cell*i+1<<"\" width=\"18\" height=\"18\"/>\n";
		}
		out<<"</g>\n";
	}
	const double legendX=(width-7*86)/2, legendY=y+h+85;
	text(width/2,legendY-14,"Probability ranges (base pairs and single-nucleotide accessibility)","legend");
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
	if(empty) text(width/2,legendY+90,"Gray / NA: empty interaction ensemble; accessibility remains defined","legend");
	out<<"</svg>\n";
	return out.str();
}

void BasePairProbabilityWriter::writeSvg(std::ostream & out,const BasePairProbabilities & result,const InteractionEnergy & energy)
{
	emit(out,svgBlock(result,energy));
}

void BasePairProbabilityWriter::writeSvgFile(const std::string & filename,const BasePairProbabilities & result,const InteractionEnergy & energy)
{
	emitFile(filename,svgBlock(result,energy));
}

} // namespace IntaRNA
