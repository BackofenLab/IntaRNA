#include <filesystem>
#include <iostream>
#include <string>
namespace IntaRNA {
std::istream* newInputStream(const std::string&);
void deleteInputStream(std::istream*&);
}
std::size_t descriptor_count() {
    std::size_t count=0;
    for (const auto& entry: std::filesystem::directory_iterator("/proc/self/fd")) {
        (void)entry; ++count;
    }
    return count;
}
int main(int argc,char** argv) {
    if (argc!=2) return 2;
    auto before=descriptor_count();
    for (int i=0;i<100;++i) {
        auto* input=IntaRNA::newInputStream(argv[1]);
        if (!input) return 3;
        std::string line;
        std::getline(*input,line);
        IntaRNA::deleteInputStream(input);
        if (input!=nullptr) return 4;
    }
    auto after=descriptor_count();
    std::cout << "{\"iterations\":100,\"fd_before\":" << before
              << ",\"fd_after\":" << after << ",\"fd_growth\":" << after-before << "}\n";
}
