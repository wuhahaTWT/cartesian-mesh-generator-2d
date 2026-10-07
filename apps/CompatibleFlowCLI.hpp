#pragma once
#include <optional>
// No selector means the existing CLI path. A compatible request is handled
// completely here, including errors, rather than silently ignoring old flags.
std::optional<int> tryCompatibleFlowCLI2D(int argc,char** argv);
