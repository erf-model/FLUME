#include "FLUME_Source.H"
#include "FLUME_Friction.H"

#include <AMReX_Print.H>

namespace flume {

std::vector<std::unique_ptr<Source>>
make_sources (const Inputs& in, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    std::vector<std::unique_ptr<Source>> sources;
    if (in.friction_manning > amrex::Real(0.0)) {
        sources.push_back(std::make_unique<Friction>(ba, dm, in.friction_manning, in.gravity,
                                                     in.dry_depth));
    }
    std::string names;
    for (const auto& s : sources) {
        if (!names.empty()) { names += ", "; }
        names += s->name();
    }
    amrex::Print() << "FLUME: sources: " << (names.empty() ? std::string("none") : names) << "\n";
    return sources;
}

} // namespace flume
