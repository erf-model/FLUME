#include "FLUME_Checkpoint.H"
#include "FLUME_WaterBalance.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>

namespace flume {

namespace {

constexpr const char* magic = "FLUME_checkpoint";
constexpr int version = 1;

bool valid_key (const std::string& key)
{
    return !key.empty() && std::all_of(key.begin(), key.end(), [] (char c) {
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9')
               || c == '_' || c == '/';
    });
}

std::string join_doubles (const std::vector<double>& v)
{
    std::ostringstream s;
    s << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const double x : v) { s << " " << x; }
    return s.str();
}

} // namespace

// --- CheckpointWriter --------------------------------------------------------

CheckpointWriter::CheckpointWriter (std::string dir) : m_dir(std::move(dir))
{
    amrex::UtilCreateCleanDirectory(m_dir, true);
    put_ints("version", {version});
    put_ints("nprocs", {amrex::ParallelDescriptor::NProcs()});
    put_ints("real_bytes", {static_cast<std::int64_t>(sizeof(amrex::Real))});
}

void
CheckpointWriter::add (const std::string& key, const std::string& line)
{
    if (!valid_key(key)) { amrex::Abort("FLUME: checkpoint key '" + key + "' is not valid"); }
    if (std::find(m_keys.begin(), m_keys.end(), key) != m_keys.end()) {
        amrex::Abort("FLUME: checkpoint key '" + key + "' written twice");
    }
    m_keys.push_back(key);
    m_lines.push_back(key + line);
}

void
CheckpointWriter::put (const std::string& key, const std::vector<double>& values)
{
    add(key, join_doubles(values));
}

void
CheckpointWriter::put_ints (const std::string& key, const std::vector<std::int64_t>& values)
{
    std::ostringstream s;
    for (const std::int64_t x : values) { s << " " << x; }
    add(key, s.str());
}

void
CheckpointWriter::put_words (const std::string& key, const std::vector<std::string>& words)
{
    std::string line;
    for (const auto& w : words) {
        if (w.empty() || w.find_first_of(" \t\n\r") != std::string::npos) {
            amrex::Abort("FLUME: checkpoint word '" + w + "' of key '" + key + "' is not valid");
        }
        line += " " + w;
    }
    add(key, line);
}

void
CheckpointWriter::put_per_rank (const std::string& key, const std::vector<double>& values)
{
    const int np = amrex::ParallelDescriptor::NProcs();
    const int io = amrex::ParallelDescriptor::IOProcessorNumber();
    std::vector<double> all(values.size() * static_cast<std::size_t>(np));
    if (!values.empty()) {
        amrex::ParallelDescriptor::Gather(values.data(), values.size(), all.data(), values.size(),
                                          io);
    }
    // Rank r's values are entries r n to r n + n - 1. Only the I/O rank holds
    // them; the others write nothing, so their line does not matter.
    add(key, join_doubles(all));
}

void
CheckpointWriter::put_multifab (const std::string& name, const amrex::MultiFab& mf)
{
    if (!valid_key(name) || name.find('/') != std::string::npos) {
        amrex::Abort("FLUME: checkpoint field name '" + name + "' is not valid");
    }
    // The valid cells only: the ghost cells are filled from them.
    amrex::MultiFab valid(mf.boxArray(), mf.DistributionMap(), mf.nComp(), 0);
    amrex::MultiFab::Copy(valid, mf, 0, 0, mf.nComp(), 0);
    amrex::VisMF::Write(valid, m_dir + "/" + name);
    put_ints("multifab/" + name, {mf.nComp()});
}

void
CheckpointWriter::finish ()
{
    if (amrex::ParallelDescriptor::IOProcessor()) {
        const std::string path = m_dir + "/Header";
        std::ofstream f(path);
        if (!f) { amrex::Abort("FLUME: cannot open " + path + " for writing"); }
        f << magic << "\n";
        for (const auto& l : m_lines) { f << l << "\n"; }
        f.close();
        if (!f) { amrex::Abort("FLUME: cannot write " + path); }
    }
    amrex::ParallelDescriptor::Barrier();
}

// --- CheckpointReader --------------------------------------------------------

CheckpointReader::CheckpointReader (std::string dir) : m_dir(std::move(dir))
{
    const std::string path = m_dir + "/Header";
    int exists = 0;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        exists = amrex::FileExists(path) ? 1 : 0;
    }
    amrex::ParallelDescriptor::Bcast(&exists, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    if (exists == 0) { fail("no file " + path + ": not a FLUME checkpoint"); }

    amrex::Vector<char> buf;
    amrex::ParallelDescriptor::ReadAndBcastFile(path, buf);
    std::istringstream in(std::string(buf.dataPtr()));
    std::string line;
    if (!std::getline(in, line) || line != magic) {
        fail("the first line of " + path + " is not '" + magic + "': not a FLUME checkpoint");
    }
    while (std::getline(in, line)) {
        std::istringstream ls(line);
        std::string key;
        if (!(ls >> key)) { continue; }
        std::vector<std::string> v;
        std::string t;
        while (ls >> t) { v.push_back(t); }
        if (!m_entries.emplace(key, std::move(v)).second) {
            fail("the key " + key + " appears twice in " + path);
        }
    }
    const std::int64_t ver = get_ints("version", 1)[0];
    if (ver != version) {
        fail("written by checkpoint version " + std::to_string(ver) + "; this build reads "
             + std::to_string(version));
    }
    const std::int64_t real_bytes = get_ints("real_bytes", 1)[0];
    const auto own_bytes = static_cast<std::int64_t>(sizeof(amrex::Real));
    if (real_bytes != own_bytes) {
        auto prec = [] (std::int64_t b) {
            return b == 8 ? std::string("double") : b == 4 ? std::string("single")
                                                           : std::to_string(b) + "-byte";
        };
        fail("written by a " + prec(real_bytes) + "-precision build; this build is "
             + prec(own_bytes) + " precision");
    }
    const std::int64_t np = get_ints("nprocs", 1)[0];
    if (np < 1 || np > std::numeric_limits<int>::max()) { fail("nprocs is not valid"); }
    m_nprocs = static_cast<int>(np);
}

void
CheckpointReader::fail (const std::string& msg) const
{
    amrex::Abort("flume.restart: " + m_dir + ": " + msg);
    std::abort();  // amrex::Abort does not return; it is not declared [[noreturn]]
}

bool
CheckpointReader::has (const std::string& key) const
{
    return m_entries.count(key) > 0;
}

const std::vector<std::string>&
CheckpointReader::tokens (const std::string& key, std::size_t n) const
{
    const auto it = m_entries.find(key);
    if (it == m_entries.end()) { fail("the checkpoint holds no " + key); }
    if (it->second.size() != n) {
        fail(key + " holds " + std::to_string(it->second.size()) + " value(s), not "
             + std::to_string(n));
    }
    return it->second;
}

std::vector<double>
CheckpointReader::get (const std::string& key, std::size_t n) const
{
    std::vector<double> v;
    for (const auto& t : tokens(key, n)) {
        errno = 0;
        char* end = nullptr;
        const double x = std::strtod(t.c_str(), &end);
        // strtod may set ERANGE for a subnormal number, which it reads
        // exactly; only an overflow is an error.
        if (end == t.c_str() || *end != '\0' || (errno == ERANGE && std::isinf(x))) {
            fail(key + " holds '" + t + "', which is not a number");
        }
        v.push_back(x);
    }
    return v;
}

std::vector<std::int64_t>
CheckpointReader::get_ints (const std::string& key, std::size_t n) const
{
    std::vector<std::int64_t> v;
    for (const auto& t : tokens(key, n)) {
        errno = 0;
        char* end = nullptr;
        const std::int64_t x = std::strtoll(t.c_str(), &end, 10);
        if (end == t.c_str() || *end != '\0' || errno == ERANGE) {
            fail(key + " holds '" + t + "', which is not an integer");
        }
        v.push_back(x);
    }
    return v;
}

std::vector<std::string>
CheckpointReader::get_words (const std::string& key) const
{
    const auto it = m_entries.find(key);
    if (it == m_entries.end()) { fail("the checkpoint holds no " + key); }
    return it->second;
}

std::vector<double>
CheckpointReader::get_per_rank (const std::string& key, std::size_t n) const
{
    const std::vector<double> all = get(key, n * static_cast<std::size_t>(m_nprocs));
    const int np = amrex::ParallelDescriptor::NProcs();
    const int me = amrex::ParallelDescriptor::MyProc();
    std::vector<double> v(n, 0.0);
    if (np == m_nprocs) {
        for (std::size_t i = 0; i < n; ++i) { v[i] = all[static_cast<std::size_t>(me) * n + i]; }
    } else if (amrex::ParallelDescriptor::IOProcessor()) {
        for (std::size_t i = 0; i < n; ++i) {
            double carry = 0.0;
            for (int r = 0; r < m_nprocs; ++r) {
                compensated_add(v[i], carry, all[static_cast<std::size_t>(r) * n + i]);
            }
        }
    }
    return v;
}

void
CheckpointReader::get_multifab (const std::string& name, amrex::MultiFab& dst) const
{
    const std::int64_t ncomp = get_ints("multifab/" + name, 1)[0];
    if (ncomp != dst.nComp()) {
        fail(name + " has " + std::to_string(ncomp) + " component(s), not "
             + std::to_string(dst.nComp()));
    }
    amrex::MultiFab src;
    amrex::VisMF::Read(src, m_dir + "/" + name);
    const amrex::BoxArray& sba = src.boxArray();
    const amrex::BoxArray& dba = dst.boxArray();
    if (src.nComp() != dst.nComp() || sba.minimalBox() != dba.minimalBox()
        || sba.numPts() != dba.numPts() || sba.ixType() != dba.ixType()) {
        std::ostringstream s;
        s << name << " covers the cells " << sba.minimalBox() << ", not those of the inputs, "
          << dba.minimalBox() << " (flume.n_cell)";
        fail(s.str());
    }
    dst.ParallelCopy(src, 0, 0, dst.nComp(), 0, 0);
}

// --- Logs ----------------------------------------------------------------------

bool
rewind_log (const std::string& path, int step, bool keep_step)
{
    std::vector<std::string> kept;
    {
        std::ifstream in(path);
        if (!in) { return false; }
        std::size_t fields = 0;  // of a data line; 0 until "# step ..." is read
        std::string line;
        while (std::getline(in, line)) {
            std::istringstream ls(line);
            std::vector<std::string> words;
            std::string w;
            while (ls >> w) { words.push_back(w); }
            if (!line.empty() && line[0] == '#') {
                kept.push_back(line);
                if (words.size() > 2 && words[0] == "#" && words[1] == "step") {
                    fields = words.size() - 1;
                }
                continue;
            }
            if (words.empty() || (fields > 0 && words.size() != fields)) { continue; }
            errno = 0;
            char* end = nullptr;
            const long long s = std::strtoll(words[0].c_str(), &end, 10);
            if (*end != '\0' || errno == ERANGE) { continue; }
            if (s < step || (keep_step && s == step)) { kept.push_back(line); }
        }
    }
    const std::string tmp = path + ".rewind";
    {
        std::ofstream out(tmp, std::ios::trunc);
        if (!out) { amrex::Abort("FLUME: cannot open " + tmp + " for writing"); }
        for (const auto& l : kept) { out << l << "\n"; }
        out.close();
        if (!out) { amrex::Abort("FLUME: cannot write " + tmp); }
    }
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        amrex::Abort("FLUME: cannot rename " + tmp + " to " + path);
    }
    return true;
}

} // namespace flume
