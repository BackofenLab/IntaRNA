// Standalone diagnostic for IntaRNA #265. This is a stress probe, not proof of thread safety.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <omp.h>
extern "C" {
#include <ViennaRNA/vrna_config.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/LPfold.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/constraints/basic.h>
#include <ViennaRNA/params/io.h>
}

struct Case {
 std::string sequence;
 std::string constraint;
 vrna_md_t model;
 int window;
};
struct Result {
 int length, window, status = 0, callbacks = 0;
 bool invalid = false;
 std::vector<double> values;
 std::vector<unsigned char> seen;
 Result(int n, int w): length(n), window(w), values((n+1)*(w+1)), seen(values.size()) {}
};
using Compound = std::unique_ptr<vrna_fold_compound_t, decltype(&vrna_fold_compound_free)>;
Compound compound(const Case &c) {
 vrna_md_t model = c.model;
 Compound fc(vrna_fold_compound(c.sequence.c_str(), &model, VRNA_OPTION_PF | VRNA_OPTION_WINDOW), &vrna_fold_compound_free);
 if (!fc) { std::cerr << "fold compound allocation failed\n"; std::exit(3); }
 if (!c.constraint.empty()) vrna_constraints_add(fc.get(), c.constraint.c_str(), VRNA_CONSTRAINT_DB_DEFAULT | VRNA_CONSTRAINT_DB_ENFORCE_BP);
 return fc;
}
void capture(FLT_OR_DBL *pr, int pr_size, int j, int max, unsigned int type, void *data) {
 auto &r = *static_cast<Result *>(data);
 if (!(type & VRNA_PROBS_WINDOW_UP)) return;
 ++r.callbacks;
 if (j < 1 || j > r.length) { r.invalid = true; return; }
 for (int l = 1; l <= std::min({j, pr_size, max, r.window}); ++l) {
  const size_t pos = j*(r.window+1)+l;
  r.values[pos] = pr[l]; r.seen[pos] = 1;
  if (!std::isfinite(pr[l]) || pr[l] < 0 || pr[l] > 1.00000001) r.invalid = true;
 }
}
void calculate(vrna_fold_compound_t *fc, const Case &c, Result &r) {
 r.status = vrna_probs_window(fc, c.window, VRNA_PROBS_WINDOW_UP, capture, &r);
}
std::vector<Case> cases() {
 std::vector<Case> cs;
 uint32_t rng = 0x265;
 for (int id=0; id<24; ++id) {
  Case c;
  const int n = 77 + id*13;
  for (int i=0; i<n; ++i) { rng ^= rng<<13; rng ^= rng>>17; rng ^= rng<<5; c.sequence += "ACGU"[rng%4]; }
  // Include stable stem motifs among pseudorandom backgrounds.
  if (id%4==0) c.sequence.replace(10, 20, "GGGGGGGAAAAAACCCCCCC");
  c.window = std::min(n, 32+(id%7)*16);
  vrna_md_set_default(&c.model);
  c.model.window_size = c.window;
  c.model.max_bp_span = (id%3==0 ? c.window-8 : c.window);
  c.model.temperature = 20 + (id%5)*8.5;
  c.model.noLP = id%2;
  c.model.noGUclosure = (id/2)%2;
  c.model.dangles = 2;
  c.model.gquad = 0;
  c.model.compute_bpp = 1;
  if (id%3==1) { c.constraint.assign(n, '.'); c.constraint[0] = c.constraint[n-1] = 'x'; }
  cs.push_back(c);
 }
 return cs;
}
int main(int argc, char **argv) {
 const int threads = argc>1 ? std::atoi(argv[1]) : 4;
 const int repetitions = argc>2 ? std::atoi(argv[2]) : 10;
 const bool lifecycle = argc>3 && std::string(argv[3])=="lifecycle";
 const bool lock = argc>4 && std::string(argv[4])=="lock";
 const bool setupLock = argc>4 && std::string(argv[4])=="setup-lock";
 if (threads<1 || repetitions<1) return 2;
 omp_set_dynamic(0); omp_set_num_threads(threads);
 vrna_params_load_RNA_Turner2004();
 const auto inputs = cases();
 std::vector<Result> reference;
 for (const auto &c: inputs) {
  reference.emplace_back(c.sequence.size(), c.window);
  auto fc=compound(c); calculate(fc.get(), c, reference.back());
  if (!reference.back().status || reference.back().invalid) { std::cerr<<"invalid serial baseline\n"; return 3; }
 }
 size_t compared = 0, mismatches = 0;
 const auto start = std::chrono::steady_clock::now();
 for (int rep=0; rep<repetitions; ++rep) {
  std::vector<Result> observed;
  std::vector<Compound> compounds;
  for (const auto &c: inputs) {
   observed.emplace_back(c.sequence.size(), c.window);
   compounds.emplace_back(lifecycle ? nullptr : compound(c).release(), &vrna_fold_compound_free);
  }
  #pragma omp parallel for schedule(dynamic, 1)
  for (int id=0; id<static_cast<int>(inputs.size()); ++id) {
   if (lifecycle) {
    if (setupLock) {
     #pragma omp critical(probe_setup)
     compounds[id] = compound(inputs[id]);
    } else compounds[id] = compound(inputs[id]);
   }
   if (lock) {
    #pragma omp critical(probe_probs_window)
    calculate(compounds[id].get(), inputs[id], observed[id]);
   } else calculate(compounds[id].get(), inputs[id], observed[id]);
   if (lifecycle) compounds[id].reset();
  }
  for (size_t id=0; id<inputs.size(); ++id) {
   const auto &a = reference[id]; const auto &b = observed[id];
   if (a.status != b.status || b.invalid || a.callbacks != b.callbacks || a.seen != b.seen) {
    std::cerr<<"metadata mismatch repetition="<<rep<<" case="<<id<<" status="<<b.status<<" invalid="<<b.invalid<<" callbacks="<<b.callbacks<<'\n'; ++mismatches;
   }
   for (size_t i=0; i<a.values.size(); ++i) if (a.seen[i]) {
    ++compared;
    if (a.values[i]!=b.values[i]) {
     if (mismatches<10) std::cerr<<std::setprecision(17)<<"probability mismatch repetition="<<rep<<" case="<<id<<" j="<<i/(a.window+1)<<" l="<<i%(a.window+1)<<" serial="<<a.values[i]<<" parallel="<<b.values[i]<<'\n';
     ++mismatches;
    }
   }
  }
  if (mismatches) break;
 }
 const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 std::cout<<"ViennaRNA="<<VRNA_VERSION<<" threads="<<threads<<" repetitions="<<repetitions<<" cases="<<inputs.size()<<" mode="<<(lifecycle?"lifecycle":"window-only")<<" lock="<<lock<<" setup_lock="<<setupLock<<" compared="<<compared<<" mismatches="<<mismatches<<" elapsed_seconds="<<elapsed<<'\n';
 return mismatches?1:0;
}
