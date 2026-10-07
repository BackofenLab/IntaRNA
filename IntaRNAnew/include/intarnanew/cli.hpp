#pragma once

#include "intarnanew/config.hpp"

#include <expected>
#include <span>
#include <string>
#include <string_view>

namespace intarnanew {

enum class OptionGroup { query, target, seed, shape, interaction, helix, output, general };
enum class OptionValueMode { none, required, optionalBoolean };
enum class OptionSupport { implemented, compatibilityOnly, unavailable };

struct OptionSpec {
    std::string_view longName;
    char shortName{};
    OptionGroup group{OptionGroup::general};
    OptionValueMode valueMode{OptionValueMode::required};
    std::string_view valueName;
    std::string_view defaultValue;
    std::string_view description;
    OptionSupport support{OptionSupport::implemented};
    bool repeatable{};
    bool basic{};
};

class Cli {
public:
    [[nodiscard]] static auto parse(std::span<const std::string_view> arguments)
        -> std::expected<Config, std::string>;
    [[nodiscard]] static auto parse(
        std::span<const std::string_view> arguments,
        std::string_view invocationName) -> std::expected<Config, std::string>;
    [[nodiscard]] static auto optionRegistry() noexcept -> std::span<const OptionSpec>;
    [[nodiscard]] static auto help(bool full) -> std::string;
    [[nodiscard]] static auto version() -> std::string;
};

} // namespace intarnanew
