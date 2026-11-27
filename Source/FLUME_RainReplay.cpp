#include "FLUME_RainReplay.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <system_error>
#include <utility>
#include <vector>

namespace flume {

namespace {

//! `x` with 17 significant digits, which read back to the same double.
std::string
digits17 (double x)
{
    std::ostringstream s;
    s.imbue(std::locale::classic());
    s << std::setprecision(std::numeric_limits<double>::max_digits10) << x;
    return s.str();
}

[[noreturn]] void
fail (const std::string& path, const std::string& msg)
{
    amrex::Abort("flume.rain.replay_file: " + path + ": " + msg);
    std::abort();  // amrex::Abort does not return; it is not declared [[noreturn]]
}

} // namespace

RainReplay::RainReplay (std::string index_path, const amrex::Geometry& geom,
                        const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
    : m_index(std::move(index_path)), m_geom(geom), m_ba(ba), m_dm(dm), m_held(ba, dm, 1, 0)
{
    const std::filesystem::path p(m_index);
    m_dir = p.parent_path().string();
    m_stem = p.stem().string();
    const std::string file = p.filename().string();
    if (m_stem.empty() || file.empty() || file == "." || file == "..") {
        fail(m_index, "not the name of a file (the grids are named after it)");
    }
    // A line of the index is a time and a name: the name holds no blank and
    // no '#', which starts a comment (FLUME_RainField.H).
    if (file.find_first_of(" \t#") != std::string::npos) {
        fail(m_index, "its name holds a blank or a '#', which a line of the index cannot hold");
    }
    m_held.setVal(amrex::Real(0.0));
    amrex::Print() << "FLUME: rain: writing the rain of the host model for a replay to "
                   << m_index << "\n";
}

void
RainReplay::start (double t)
{
    m_held.setVal(amrex::Real(0.0));
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::error_code ec;
        if (!m_dir.empty()) {
            std::filesystem::create_directories(m_dir, ec);
            if (ec) { fail(m_index, "cannot make the directory " + m_dir + ": " + ec.message()); }
        }
        if (std::filesystem::is_directory(m_index, ec)) { fail(m_index, "is a directory"); }
    }
    amrex::ParallelDescriptor::Barrier();
    amrex::MultiFab zero(m_ba, m_dm, 1, 0);
    zero.setVal(amrex::Real(0.0));
    write_grid(zero, zero_name());
    m_grids = 0;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        const std::string line = digits17(t) + " " + zero_name() + "\n";
        std::ofstream out(m_index, std::ios::binary | std::ios::trunc);
        out << line;
        out.close();
        if (!out) { fail(m_index, "cannot be written"); }
        m_close_at = 0;
        m_close_time = t;
        m_close_step = 0.0;
        m_size = line.size();
    }
}

void
RainReplay::resume (const CheckpointReader& chk)
{
    if (!chk.has("multifab/rain_replay_held") || !chk.has("rain_replay_index")
        || !chk.has("rain_replay_close")) {
        chk.fail("it was written without flume.rain.replay_file; a run that writes the rain"
                 " for a replay writes it from its start");
    }
    chk.get_multifab("rain_replay_held", m_held);
    const std::vector<std::int64_t> ints = chk.get_ints("rain_replay_index", 3);
    const std::vector<double> close = chk.get("rain_replay_close", 2);
    if (ints[0] < 0 || ints[2] < 0 || ints[1] <= ints[2] || !std::isfinite(close[0])
        || !(close[1] >= 0.0)) {
        chk.fail("rain_replay_index or rain_replay_close is not valid");
    }
    m_grids = ints[0];
    const auto size = static_cast<std::uintmax_t>(ints[1]);
    const auto close_at = static_cast<std::uintmax_t>(ints[2]);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream in(m_index, std::ios::binary);
        if (!in) {
            fail(m_index, "cannot be opened: a run that continues from a checkpoint cuts back"
                          " and goes on with the index of the run before it");
        }
        // The index as the run wrote it up to the checkpoint: `m_grids` lines
        // of rain and lines of no rain in its first `close_at` bytes, then
        // its closing line, at the time of the checkpoint's closing line.
        // The run may since have written over the closing line (a step that
        // joined it) and gone on; that is cut off, and the closing line is
        // written again. So is a line cut short by a run that was stopped.
        std::string text(static_cast<std::size_t>(close_at), '\0');
        in.read(text.data(), static_cast<std::streamsize>(close_at));
        bool same = static_cast<std::uintmax_t>(in.gcount()) == close_at
                    && (close_at == 0 || text.back() == '\n');
        std::string next;
        same = same && static_cast<bool>(std::getline(in, next)) && !in.eof();
        in.close();
        const std::string closing = digits17(close[0]) + " " + zero_name() + "\n";
        if (same) {
            // The line at close_at: the closing line, or the line of the
            // next grid of rain at the same time.
            std::ostringstream num;
            num << std::setw(6) << std::setfill('0') << (m_grids + 1);
            const std::string time = digits17(close[0]) + " ";
            same = close_at + closing.size() == size
                   && (next + "\n" == closing
                       || next == time + m_stem + "_" + num.str() + ".asc");
        }
        long long grids = 0;
        std::size_t pos = 0;
        while (same && pos < text.size()) {
            const std::size_t eol = text.find('\n', pos);
            std::istringstream s(text.substr(pos, eol - pos));
            s.imbue(std::locale::classic());
            double time = 0.0;
            std::string name;
            std::string extra;
            if (!(s >> time >> name) || (s >> extra)) {
                same = false;
            } else if (name != zero_name()) {
                ++grids;
            }
            pos = eol + 1;
        }
        if (!same || grids != m_grids) {
            std::ostringstream msg;
            msg << std::setprecision(17) << "does not hold what the run wrote up to the checkpoint ("
                << m_grids << " grids of rain, closed at " << close[0]
                << " s): a run that continues from a checkpoint goes on with the index of the"
                   " run that wrote the checkpoint";
            fail(m_index, msg.str());
        }
        write_at(close_at, closing);
        m_size = size;
        m_close_at = close_at;
        m_close_time = close[0];
        m_close_step = close[1];
    }
    // The grid of no rain, as the run before wrote it.
    amrex::MultiFab zero(m_ba, m_dm, 1, 0);
    zero.setVal(amrex::Real(0.0));
    write_grid(zero, zero_name());
}

void
RainReplay::record (const amrex::MultiFab& depth, double t0, double t1)
{
    AMREX_ASSERT(depth.boxArray() == m_ba && depth.DistributionMap() == m_dm);
    // contains_nan() is collective, so every rank stops.
    if (depth.contains_nan(0, 1)) {
        fail(m_index, "the rain of the host step is not a number in a cell");
    }
    // What is written: the depth and what is held back, where it is above 0;
    // the rest is held back.
    amrex::MultiFab w(m_ba, m_dm, 1, 0);
    {
        const auto& d = depth.const_arrays();
        const auto& h = m_held.arrays();
        const auto& o = w.arrays();
        amrex::ParallelFor(w, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
        {
            const amrex::Real x = d[b](i, j, k) + h[b](i, j, k);
            o[b](i, j, k) = amrex::max(x, amrex::Real(0.0));
            h[b](i, j, k) = amrex::min(x, amrex::Real(0.0));
        });
        amrex::Gpu::streamSynchronize();
    }
    // max() is collective; a step of no rain adds no line.
    if (!(w.max(0) > amrex::Real(0.0))) { return; }

    // The rate in mm/h, in double: depth / (t1 - t0) * 3.6e6.
    const double scale = 3.6e6 / (t1 - t0);
    std::ostringstream num;
    num << std::setw(6) << std::setfill('0') << (m_grids + 1);
    const std::string name = m_stem + "_" + num.str() + ".asc";
    write_grid(w, name, scale);
    ++m_grids;

    if (amrex::ParallelDescriptor::IOProcessor()) {
        // A step that starts where the last one with rain ended takes the
        // place of its closing line. The tolerance is the engine's
        // (follow_clock): in single precision the host hands over its step
        // rounded to float, so the end of the last step and the start of this
        // one can differ by that rounding.
        const double tol = 1.0e-6 * std::max(t1 - t0, m_close_step)
                           + 4.0 * std::numeric_limits<double>::epsilon() * std::abs(t0);
        if (t0 < m_close_time - tol) {
            std::ostringstream msg;
            msg << std::setprecision(17) << "the host step from " << t0
                << " s starts before the end of the last step with rain, " << m_close_time
                << " s";
            fail(m_index, msg.str());
        }
        const bool joined = t0 <= m_close_time + tol;
        const std::uintmax_t base = joined ? m_close_at : m_size;
        const std::string rain = digits17(joined ? m_close_time : t0) + " " + name + "\n";
        const std::string close = digits17(t1) + " " + zero_name() + "\n";
        write_at(base, rain + close);
        m_close_at = base + rain.size();
        m_close_time = t1;
        m_close_step = t1 - t0;
        m_size = m_close_at + close.size();
    }
}

void
RainReplay::write_checkpoint (CheckpointWriter& chk) const
{
    chk.put_multifab("rain_replay_held", m_held);
    // The state of the index, kept on the I/O rank.
    std::int64_t ints[2] = {static_cast<std::int64_t>(m_size),
                            static_cast<std::int64_t>(m_close_at)};
    double close[2] = {m_close_time, m_close_step};
    const int io = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(ints, 2, io);
    amrex::ParallelDescriptor::Bcast(close, 2, io);
    chk.put_ints("rain_replay_index", {static_cast<std::int64_t>(m_grids), ints[0], ints[1]});
    chk.put("rain_replay_close", {close[0], close[1]});
}

void
RainReplay::write_grid (const amrex::MultiFab& depth, const std::string& name, double scale) const
{
    // The whole grid on the I/O rank, in memory the host can read.
    const amrex::Box& domain = m_geom.Domain();
    const int io = amrex::ParallelDescriptor::IOProcessorNumber();
    const amrex::BoxArray one(domain);
    const amrex::DistributionMapping one_dm(amrex::Vector<int>{io});
    amrex::MultiFab all(one, one_dm, 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    all.ParallelCopy(depth, 0, 0, 1);
    amrex::Gpu::streamSynchronize();
    if (!amrex::ParallelDescriptor::IOProcessor()) { return; }

    const amrex::Real* v = all[0].dataPtr();
    const int nx = domain.length(0);
    const int ny = domain.length(1);
    const double dx = static_cast<double>(m_geom.CellSize(0));
    const double dy = static_cast<double>(m_geom.CellSize(1));
    const std::string path = m_dir.empty() ? name : m_dir + "/" + name;
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) { fail(m_index, "cannot write the grid " + path); }
    std::ostringstream s;
    s.imbue(std::locale::classic());
    s << "ncols " << nx << "\nnrows " << ny << "\nxllcorner "
      << digits17(static_cast<double>(m_geom.ProbLo(0))) << "\nyllcorner "
      << digits17(static_cast<double>(m_geom.ProbLo(1))) << "\n";
    if (dx == dy) {
        s << "cellsize " << digits17(dx) << "\n";
    } else {
        s << "dx " << digits17(dx) << "\ndy " << digits17(dy) << "\n";
    }
    out << s.str();
    // Rows from the top (the largest y) down, one at a time.
    for (int j = ny - 1; j >= 0; --j) {
        std::ostringstream row;
        row.imbue(std::locale::classic());
        row << std::setprecision(std::numeric_limits<double>::max_digits10);
        for (int i = 0; i < nx; ++i) {
            const double x = static_cast<double>(v[static_cast<std::size_t>(i)
                                                   + static_cast<std::size_t>(nx)
                                                         * static_cast<std::size_t>(j)])
                             * scale;
            row << (i > 0 ? " " : "") << x;
        }
        row << "\n";
        out << row.str();
    }
    out.close();
    if (!out) { fail(m_index, "cannot write the grid " + path); }
}

void
RainReplay::write_at (std::uintmax_t at, const std::string& lines) const
{
    // Write over the end of the file in place: the file never loses its
    // closing line between two steps (a step with rain writes its line over
    // the closing line before it, and a new closing line after it).
    std::error_code ec;
    const std::uintmax_t now = std::filesystem::file_size(m_index, ec);
    if (ec) { fail(m_index, "cannot be read: " + ec.message()); }
    if (now < at) { fail(m_index, "has been cut short by something other than this run"); }
    std::fstream f(m_index, std::ios::in | std::ios::out | std::ios::binary);
    f.seekp(static_cast<std::streamoff>(at));
    f << lines;
    f.close();
    if (!f) { fail(m_index, "cannot be written"); }
    if (now > at + lines.size()) {
        std::filesystem::resize_file(m_index, at + lines.size(), ec);
        if (ec) { fail(m_index, "cannot be cut back: " + ec.message()); }
    }
}

} // namespace flume
