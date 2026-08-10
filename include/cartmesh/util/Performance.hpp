#pragma once

#include <cstdint>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

namespace cartmesh {

[[nodiscard]] inline std::uint64_t peak_rss_bytes() noexcept {
#if defined(__unix__) || defined(__APPLE__)
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return 0;
    }
#if defined(__APPLE__)
    return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
    return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024ULL;
#endif
#else
    return 0;
#endif
}

} // 命名空间 cartmesh
