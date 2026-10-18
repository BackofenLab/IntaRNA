#include <IntaRNA/Matrix.h>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <string>

using Key = std::pair<std::size_t, std::size_t>;
std::uint64_t checks = 0;
void require(bool ok) { ++checks; if (!ok) throw std::runtime_error("audit mismatch at check " + std::to_string(checks)); }
struct Oracle {
    int kind;
    std::size_t rows=0, cols=0, upper=0;
    std::map<Key,int> cells;
    bool stored(std::size_t i, std::size_t j) const {
        return kind==0 || (i<=j && (kind==1 || j-i<=upper));
    }
    int get(std::size_t i, std::size_t j) const {
        auto pos=cells.find({i,j});
        return pos==cells.end() ? 0 : pos->second;
    }
    void reshape(std::size_t r, std::size_t c, std::size_t u, bool preserve) {
        rows=r; cols=c; upper=u;
        for (auto it=cells.begin();it!=cells.end();) {
            const auto [i,j]=it->first;
            if (!preserve || i>=r || j>=c || !stored(i,j)) it=cells.erase(it);
            else ++it;
        }
    }
};
template<class M> void verify(const M& m, const Oracle& o) {
    require(m.size1()==o.rows); require(m.size2()==o.cols);
    for(std::size_t i=0;i<o.rows;++i)
        for(std::size_t j=0;j<o.cols;++j) require(m(i,j)==o.get(i,j));
}
template<class M> void fill(M& m, Oracle& o, std::mt19937& rng) {
    for(std::size_t i=0;i<o.rows;++i)
        for(std::size_t j=0;j<o.cols;++j)
            if(o.stored(i,j)) m(i,j)=o.cells[{i,j}]=int(rng()%2001)-1000;
}
void resize(IntaRNA::Matrix<int>& m, const Oracle& o, bool preserve) { m.resize(o.rows,o.cols,preserve); }
void resize(IntaRNA::UpperTriangularMatrix<int>& m, const Oracle& o, bool preserve) { m.resize(o.rows,o.cols,preserve); }
void resize(IntaRNA::UpperBandedMatrix<int>& m, const Oracle& o, bool preserve) { m.resize(o.rows,o.cols,0,o.upper,preserve); }
template<class M> void audit(int kind) {
    std::mt19937 rng(246+kind);
    M m;
    Oracle oracle{kind};
    for(int step=0;step<10000;++step) {
        switch(rng()%7) {
            case 0: case 1: {
                auto r=rng()%41, c=kind==1 ? r : rng()%41, u=rng()%61;
                bool preserve=rng()%2;
                oracle.reshape(r,c,u,preserve);
                resize(m,oracle,preserve);
                if(!preserve) fill(m,oracle,rng);
                break;
            }
            case 2: fill(m,oracle,rng); break;
            case 3: m.clear(); oracle.cells.clear(); break;
            case 4: {
                M copied(m); verify(copied,oracle);
                Oracle changed=oracle;
                fill(copied,changed,rng);
                verify(m,oracle); verify(copied,changed);
                M assigned; assigned=copied;
                copied.clear(); verify(assigned,changed);
                m=assigned; oracle=changed;
                break;
            }
            case 5: {
                M moved(std::move(m)); verify(moved,oracle);
                require(m.size1()==0 && m.size2()==0);
                resize(m,oracle,false); m.clear();
                verify(moved,oracle);
                m=std::move(moved);
                require(moved.size1()==0 && moved.size2()==0);
                auto* alias=&m; m=std::move(*alias);
                break;
            }
            case 6: {
                Oracle other{kind};
                auto r=rng()%41, c=kind==1 ? r : rng()%41;
                other.reshape(r,c,rng()%61,false);
                M swapped; resize(swapped,other,false); fill(swapped,other,rng);
                m.swap(swapped); verify(swapped,oracle); oracle=other;
                break;
            }
        }
        verify(m,oracle);
    }
}
int main() {
    audit<IntaRNA::Matrix<int>>(0);
    audit<IntaRNA::UpperTriangularMatrix<int>>(1);
    audit<IntaRNA::UpperBandedMatrix<int>>(2);
    std::cout << "Passed 30000 seeded operations, " << checks << " independent reference checks\n";
}
