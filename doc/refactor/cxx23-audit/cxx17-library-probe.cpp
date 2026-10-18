#include <unordered_map>
int main() { std::unordered_map<int,double> p; auto [it,ok]=p.try_emplace(1,2.0); return it->second==2.0 && ok ? 0 : 1; }
