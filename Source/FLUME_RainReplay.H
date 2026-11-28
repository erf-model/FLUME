#ifndef FLUME_RAINREPLAY_H_
#define FLUME_RAINREPLAY_H_

#include "FLUME_Checkpoint.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <cstdint>
#include <string>

namespace flume {

/**
 * Writes the rain a host model hands over (flume.rain.type = microphysics)
 * as gridded rain (flume.rain.type = gridded, FLUME_RainField.H), so that a
 * stand-alone run can read it back without the host (flume.rain.replay_file).
 *
 * The index file (the path of flume.rain.replay_file) holds lines of a time
 * in s and the name of an ESRI ASCII grid of the rain rate in mm/h, in the
 * index file's directory. Each host step with rain from t0 to t1 adds the
 * line "t0 <stem>_<n>.asc", where the grid is the depth of the step over
 * t1 - t0 and n counts the grids from 000001. The index always ends with the
 * line "t <stem>_zero.asc", a grid of no rain at the end of the last step
 * with rain (or at the start of the run): read with
 * flume.rain.interpolation = step, the gridded rain gives the host's rain
 * step by step and no rain after it, for as long as the stand-alone run goes
 * on. A step of no rain adds no line, and the closing line of the step
 * before it then stands for the dry spell.
 *
 * The index is written in place: a step with rain writes its line over the
 * closing line before it and a new closing line after it, so a run that is
 * stopped between two steps leaves an index that ends with its closing line.
 * The checkpoint holds the state of the index (its size, its closing line and
 * the grids of rain in it); a restart cuts the index back to that state and
 * stops if the index does not hold it.
 *
 * A negative depth (a correction from a finer level of the host) is not
 * written in its step: it is held per cell and taken from the rain of the
 * steps that follow. The coupled Rain source (FLUME_Rain.H) takes such a
 * deficit from the depth at the next step with rain; a rate cannot be
 * negative, so the replay takes it from that rain only, and the two agree
 * once the rain has made up the deficit. What is held goes into the
 * checkpoint.
 *
 * Rates are written with 17 significant digits. The stand-alone run reads
 * them back to within 2 units of round-off of the host's depth per step.
 */
class RainReplay
{
public:
    //! For the cells of `geom` (one box, the union of `ba`) in the boxes
    //! `ba` on the ranks `dm`, writing the index file `index_path`. Aborts,
    //! naming flume.rain.replay_file, when the path does not end in a file
    //! name, or the name holds a blank or a '#'.
    RainReplay (std::string index_path, const amrex::Geometry& geom, const amrex::BoxArray& ba,
                const amrex::DistributionMapping& dm);

    //! A run that starts at time `t`, s: write the grid of no rain and an
    //! index of the one line "t <stem>_zero.asc" (an existing index is
    //! replaced). Collective.
    void start (double t);

    /**
     * A run that continues from the checkpoint `chk`: read what is held back
     * and the state of the index, and cut the index back to what the run had
     * written when it wrote the checkpoint. Aborts, naming
     * flume.rain.replay_file, when the index does not exist or does not
     * begin with what the run had written, and naming flume.restart when the
     * checkpoint was written without flume.rain.replay_file. Collective.
     */
    void resume (const CheckpointReader& chk);

    /**
     * The host step from `t0` to `t1` (s, t1 > t0) with the depth of rain
     * `depth`, m, on the cells of `ba` (component 0, valid cells). Writes a
     * grid and the index when a cell has rain. Collective.
     */
    void record (const amrex::MultiFab& depth, double t0, double t1);

    //! What is held back, as rain_replay_held, and the state of the index,
    //! as rain_replay_index (grids, size, start of the closing line) and
    //! rain_replay_close (time of the closing line, length of the last step
    //! with rain). Collective.
    void write_checkpoint (CheckpointWriter& chk) const;

    //! The name of the grid of no rain, <stem>_zero.asc.
    [[nodiscard]] std::string zero_name () const { return m_stem + "_zero.asc"; }

private:
    //! Write `depth` times `scale` (component 0, valid cells) as the grid
    //! `name` in the index file's directory. Collective; the I/O rank writes.
    void write_grid (const amrex::MultiFab& depth, const std::string& name,
                     double scale = 1.0) const;

    //! Write `lines` over the index from byte `at` on, and cut it after
    //! them. I/O rank.
    void write_at (std::uintmax_t at, const std::string& lines) const;

    std::string m_index;       //!< the index file
    std::string m_dir;         //!< its directory ("" for the current one)
    std::string m_stem;        //!< its name without the extension
    amrex::Geometry m_geom;
    amrex::BoxArray m_ba;
    amrex::DistributionMapping m_dm;
    amrex::MultiFab m_held;    //!< the negative depth held back, m (<= 0)
    // On the I/O rank only (m_grids on every rank):
    std::uintmax_t m_size = 0;        //!< the index file's size, bytes
    std::uintmax_t m_close_at = 0;    //!< where its closing line starts, bytes
    double m_close_time = 0.0;        //!< the time of its closing line, s
    double m_close_step = 0.0;        //!< the length of the last step with rain, s
    long long m_grids = 0;            //!< the grids of rain written so far
};

} // namespace flume

#endif
