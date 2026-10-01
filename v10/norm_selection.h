#pragma once
#include "format.h"
#include <sstream>
namespace hic10 {
inline void select_norms(const std::string &names, bool &vc, bool &vcs, bool &scale, bool &selected) {
    if (!selected) { vc = vcs = scale = false; selected = true; }
    std::istringstream input(names);
    std::string name;
    check(!names.empty() && names.back() != ',', "empty normalization selection");
    while (std::getline(input, name, ',')) {
        if (name == "VC") vc = true;
        else if (name == "VC_SQRT") vcs = true;
        else if (name == "SCALE") scale = true;
        else check(false, "unknown computed normalization " + name + "; use --vectors to import another type");
    }
}
}
