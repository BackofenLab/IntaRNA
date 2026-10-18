#include <cstddef>
#include <sstream>
#include <string>
#include <string_view>
int main() {
std::string text{"IntaRNA"};
const bool has_rna = text.contains("RNA") && text.contains('R');
text.resize_and_overwrite(8, [](char * data, std::size_t size) noexcept {
	for (std::size_t i = 0; i < size; ++i)
		data[i] = 'x';
	return size;
});
std::stringstream stream;
stream << text;
const std::string_view stream_view = stream.view();
std::stringstream sink;
sink << stream.view();
(void)has_rna;
(void)stream_view;
}
