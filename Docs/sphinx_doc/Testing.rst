Testing
=======

Running the tests
-----------------

.. code-block:: bash

   ctest --test-dir build-rel --output-on-failure
   ctest --test-dir build-rel -L regression
   ctest --test-dir build-rel -R dam_break

Every test runs in its own directory, ``<build>/Tests/run/<name>``, made
again at the start of the test. Nothing is written into the source tree.
``Tests/run_case.py`` runs the executable and then the check; the test fails
if either fails. Without MPI only the 1-rank tests are registered.

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Label
     - Tests
   * - ``regression``
     - Runs whose output is checked against closed forms or pinned numbers.
   * - ``abort``
     - Runs that must stop with a non-zero exit code and a message naming
       the input.

:doc:`Verification` gives the analytical cases and their errors;
:doc:`TestCatalog` lists every test.

Adding a test
-------------

``Tests/CMakeLists.txt`` has two helpers:

.. code-block:: cmake

   flume_add_regtest(<name> NPROCS <n> INPUTS <file> [ARGS <key=value>...]
                     [OUTPUT <regular expression>]
                     CHECK <command and its arguments>)

   flume_add_abort_test(<name> EXPECT <regular expression> ARGS <key=value>...)

The check command runs in the test's directory. With ``OUTPUT``, the output
of the run must also match the regular expression. A new test is first shown
to fail against the defect it guards, or against a one-line change to the
code, and only then committed.

A case with its own files (decks, data, check script) lives in one folder
under ``Tests`` with a ``README.md``, as ``Tests/Rain`` and ``Tests/Raster``
do.

Continuous integration
----------------------

Three workflows run on every pull request to ``development``.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Workflow
     - Jobs
   * - ``ci.yml``
     - GCC and Clang with MPI, and GCC Debug with AMReX's testing mode:
       build with warnings as errors, then every test. No MPI with single
       precision: build and the 1-rank tests. CUDA, HIP and SYCL: compile
       only.
   * - ``style.yml``
     - ``codespell``, no tabs, no trailing whitespace, ``actionlint``.
   * - ``docs.yml``
     - This documentation, with warnings as errors.
