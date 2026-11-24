#ifndef FLUME_TIMETABLE_H_
#define FLUME_TIMETABLE_H_

#include "FLUME_Inputs.H"

#include <string>
#include <vector>

namespace flume {

/**
 * A value that is a function of time, given by a table: times
 * t_0 < t_1 < ... < t_N (s) and values v_0 ... v_N (>= 0). With
 * RainInterpolation::step the value is v_k from t_k to t_{k+1}; with
 * RainInterpolation::linear it is linear in time between the entries. Before
 * t_0 the value is v_0, and from t_N on it is v_N, in both modes.
 *
 * The rate of a hyetograph (FLUME_Forcing.H) and the discharge of an inflow
 * hydrograph (FLUME_Inflow.H) are such tables. The integral over an interval
 * is I(t1) - I(t0), with I the integral of the value from t_0, exact for the
 * piecewise-constant or piecewise-linear value.
 */
class TimeTable
{
public:
    //! `time` strictly increasing, `value` >= 0 and finite, the same size, at
    //! least one entry.
    TimeTable (std::vector<double> time, std::vector<double> value, RainInterpolation mode);

    //! The value at time t.
    [[nodiscard]] double value (double t) const;
    //! The integral of the value from t_0 to t (below 0 before t_0).
    [[nodiscard]] double integral (double t) const;
    //! The largest value at any time from t0 to t1 (t1 >= t0; t1 may be
    //! +infinity).
    [[nodiscard]] double max_value (double t0, double t1) const;

    [[nodiscard]] const std::vector<double>& times () const { return m_time; }
    [[nodiscard]] const std::vector<double>& values () const { return m_value; }

private:
    std::vector<double> m_time;
    std::vector<double> m_value;
    std::vector<double> m_cumulative;  //!< integral from t_0 to t_k
    RainInterpolation m_mode = RainInterpolation::step;
};

/**
 * Read a table file of the input `key` (for example flume.rain.file): lines
 * of two numbers, a time in s and a `quantity` in `unit` (for example "rain
 * rate" and "mm/h"), with the times strictly increasing; `#` starts a comment
 * that runs to the end of the line, and blank lines are skipped. The I/O
 * rank reads the file and sends the table to every rank. Every rank aborts
 * with the same message, naming the input, the file and the line, when the
 * file cannot be opened or is not a regular file, a line does not hold
 * exactly two finite numbers, a value is below 0, a time is not above the one
 * before, or the file holds no entry. Returns the times (s) and the values
 * (as written). Collective.
 */
void read_time_table (const std::string& key, const std::string& path, const std::string& quantity,
                      const std::string& unit, std::vector<double>& time,
                      std::vector<double>& value);

} // namespace flume

#endif
