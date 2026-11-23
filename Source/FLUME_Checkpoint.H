#ifndef FLUME_CHECKPOINT_H_
#define FLUME_CHECKPOINT_H_

#include <AMReX_MultiFab.H>

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace flume {

/**
 * Writes a checkpoint: a directory that holds a text file `Header` and the
 * MultiFabs, each written with amrex::VisMF in the precision of the build.
 *
 * The header holds one entry per line, a key and its values. Doubles are
 * written with 17 significant digits, which read back to the same bits.
 * Every rank calls every method in the same order (some are collective).
 * Keys are letters, digits, `_` and `/`; a key may be written once.
 */
class CheckpointWriter
{
public:
    //! Make the directory `dir` (an existing one is renamed, as amrex does
    //! for plotfiles). Collective.
    explicit CheckpointWriter (std::string dir);

    //! Values that are the same on every rank.
    void put (const std::string& key, const std::vector<double>& values);
    void put_ints (const std::string& key, const std::vector<std::int64_t>& values);
    void put_words (const std::string& key, const std::vector<std::string>& words);

    //! Values that differ between ranks (`values` holds this rank's, the same
    //! number on every rank): the header keeps those of every rank. Collective.
    //! They are read back by rank (get_per_rank) and only their sum over the
    //! ranks is kept on another rank count: use it only for values that enter
    //! the results as a sum over the ranks.
    void put_per_rank (const std::string& key, const std::vector<double>& values);

    //! The valid cells of `mf`, in the file <dir>/<name>. Collective.
    void put_multifab (const std::string& name, const amrex::MultiFab& mf);

    //! Write the header on the I/O rank. Collective.
    void finish ();

    [[nodiscard]] const std::string& dir () const { return m_dir; }

private:
    void add (const std::string& key, const std::string& line);

    std::string m_dir;
    std::vector<std::string> m_keys;
    std::vector<std::string> m_lines;
};

/**
 * Reads a checkpoint written by CheckpointWriter. Every rank reads the
 * header (the I/O rank reads the file and sends it). The constructor aborts
 * on a key given twice; every method aborts, naming flume.restart and the
 * checkpoint, when a key is missing or holds a different number of values
 * than asked for.
 */
class CheckpointReader
{
public:
    explicit CheckpointReader (std::string dir);

    [[nodiscard]] bool has (const std::string& key) const;
    [[nodiscard]] std::vector<double> get (const std::string& key, std::size_t n) const;
    [[nodiscard]] std::vector<std::int64_t> get_ints (const std::string& key,
                                                      std::size_t n) const;
    [[nodiscard]] std::vector<std::string> get_words (const std::string& key) const;

    /**
     * This rank's values of a put_per_rank key, `n` of them. On the rank
     * count the checkpoint was written with, rank r gets what rank r wrote.
     * On another rank count the I/O rank gets the sum over the ranks of each
     * value (compensated) and the other ranks get 0: the sum over the ranks is
     * kept to round-off, not to the bit.
     */
    [[nodiscard]] std::vector<double> get_per_rank (const std::string& key, std::size_t n) const;

    /**
     * Fill the valid cells of `dst` from the file <dir>/<name>, on any boxes
     * and ranks: the file must cover the same cells of the domain with the
     * same number of components. The ghost cells of `dst` are not set.
     * Collective.
     */
    void get_multifab (const std::string& name, amrex::MultiFab& dst) const;

    //! The rank count the checkpoint was written with.
    [[nodiscard]] int nprocs () const { return m_nprocs; }

    [[nodiscard]] const std::string& dir () const { return m_dir; }

    //! Abort with "flume.restart: <dir>: <msg>".
    [[noreturn]] void fail (const std::string& msg) const;

private:
    [[nodiscard]] const std::vector<std::string>& tokens (const std::string& key,
                                                          std::size_t n) const;

    std::string m_dir;
    std::map<std::string, std::vector<std::string>> m_entries;
    int m_nprocs = 1;
};

/**
 * Cut the log file `path` back to step `step` and leave it to be appended
 * to. Kept, in order: the comment lines (starting with `#`), and the lines
 * whose first number, the step, is below `step`, or equal to it with
 * `keep_step`, and that have as many fields as the comment line that names
 * the columns ("# step ..."), so that a line cut short by a run that was
 * stopped is removed. The file is written anew under another name and then
 * renamed over the old one. Returns false when the file does not exist.
 * I/O rank only.
 */
bool rewind_log (const std::string& path, int step, bool keep_step);

} // namespace flume

#endif
