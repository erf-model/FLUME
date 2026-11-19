Testing
=======

Running the tests
-----------------

.. code-block:: bash

   ctest --test-dir build-rel --output-on-failure

Every test runs in its own directory, ``<build>/Tests/run/<name>``, which is
removed and made again at the start of the test. Nothing is written into the
source tree. ``Tests/run_case.py`` runs the executable and then the check; the
test fails if either fails.

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Label
     - Tests
   * - ``regression``
     - Runs of the executable whose output is checked against pinned numbers.
   * - ``abort``
     - Runs that must stop with a non-zero exit code and a message naming the
       input.

Select a label with ``ctest -L <label>``.

Current tests
-------------

``skeleton_np1``, ``skeleton_np2``, ``skeleton_np4``
   ``Tests/Skeleton/inputs`` on 1, 2 and 4 MPI ranks: 48 by 32 cells of 2 m on
   [-10, 86] by [5, 69] m, in 6 boxes of at most 16 cells, written to
   ``skel00000``. Every optional key has a value that is not its default.
   ``Tests/check_plotfile.py`` reads the plotfile's headers and checks the
   variable name, the domain, the cell size, that the boxes cover every cell
   once, that the depth is zero in every box, and that the boxes were written
   to as many data files as the run has ranks. Without MPI only the 1-rank
   test is registered.

   The check compares the cell size with a relative tolerance of 1e-12 and the
   corners with 1e-12 of the domain extent. The pinned numbers are exact in
   single and double precision. The depth is read from the minimum and maximum
   in the plotfile's header, not from the data files.

``skeleton_defaults``
   ``Tests/Skeleton/inputs_defaults`` gives only the required keys: 128 by 65
   cells. The check pins the defaults: lower corner at the origin, 4 boxes
   (``flume.max_grid_size = 64``) and the plotfile name ``plt00000``.

``abort_*``
   One test per check in :doc:`Inputs`: missing keys, wrong numbers of values,
   values outside the range (including integers that would wrap), corners that
   are not finite, not ordered, too large, or too close for the cell size, an
   unknown key in the ``flume`` namespace and an unused key in another
   namespace. The wrong-type values that AMReX's parser stops are not tested.

Adding a test
-------------

``Tests/CMakeLists.txt`` has two helpers.

.. code-block:: cmake

   flume_add_regtest(<name> NPROCS <n> INPUTS <file> [ARGS <key=value>...]
                     CHECK <command and its arguments>)

   flume_add_abort_test(<name> EXPECT <regular expression> ARGS <key=value>...)

The check command runs in the test's directory. A new test is first shown to
fail against the defect it guards or against a one-line change to the code,
and only then committed.

Continuous integration
----------------------

Three workflows run on every pull request to ``development``.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Workflow
     - Jobs
   * - ``ci.yml``
     - GCC and Clang with MPI, GCC with OpenMP, and GCC Debug with AMReX's
       testing mode: build with warnings as errors, then every test. No MPI
       with single precision: build and the 1-rank tests. CUDA, HIP and SYCL:
       compile only, because the hosted runners have no GPU.
   * - ``style.yml``
     - ``codespell``, no tabs, no trailing whitespace, and ``actionlint`` on
       the workflow files.
   * - ``docs.yml``
     - This documentation, with warnings as errors.
