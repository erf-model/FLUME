Building
========

Requirements
------------

* CMake 3.25 or newer.
* A C++20 compiler.
* MPI (optional; on by default).
* Python 3 for the tests (or ``-DFLUME_ENABLE_TESTS=OFF``).

Getting the code
----------------

FLUME needs two submodules: ERF, and the AMReX that ERF pins. The stand-alone
build compiles no ERF code and needs none of ERF's other submodules.

.. code-block:: bash

   git clone https://github.com/erf-model/FLUME.git
   cd FLUME
   git submodule update --init Submodules/ERF
   git -C Submodules/ERF submodule update --init --depth 1 Submodules/AMReX

``git clone --recursive`` also works; it fetches every submodule of ERF.

Configure, build and test
-------------------------

.. code-block:: bash

   cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
   cmake --build build-rel -j 8
   ctest --test-dir build-rel --output-on-failure

The executable is ``build-rel/flume``. The default build type is ``Release``.
Build directories named ``build-<something>`` are ignored by git.

The flood grid is two-dimensional. AMReX is built in three dimensions with
boxes one cell thick, so the same code links with ERF.

Building ERF with FLUME
-----------------------

With ``-DFLUME_WITH_ERF=ON``, FLUME adds ERF as a CMake subproject, the way
ERF's own ``Tests/ParentProject`` and GALE do. One configure builds ``flume``
and ERF's ``erf_exec``, and both link one AMReX, built with ERF's AMReX
options. ERF needs no submodule other than AMReX.

.. code-block:: bash

   cmake -S . -B build-erf -DCMAKE_BUILD_TYPE=Release -DFLUME_WITH_ERF=ON
   cmake --build build-erf -j 8
   ctest --test-dir build-erf --output-on-failure

``erf_exec`` is in ``build-erf/ERF/Exec``. ERF is built with its flood
coupling (``ERF_ENABLE_FLOOD=ON``), and FLUME is compiled into ``erf_exec`` as
the flood engine ``flume`` (:ref:`sec:InsideERF`). FLUME adds
``Source/ERFEngine`` to ERF's ``ERF_FLOOD_ENGINE_DIRS``; directories of other
engines given in that cache entry are kept.

FLUME's options set ERF's: ``FLUME_MPI`` sets ``ERF_ENABLE_MPI``,
``FLUME_PRECISION`` sets ``ERF_PRECISION`` and ``FLUME_GPU_BACKEND`` sets
``ERF_ENABLE_CUDA``, ``ERF_ENABLE_HIP`` or ``ERF_ENABLE_SYCL``. OpenMP stays
off. An ERF cache entry that disagrees with both the FLUME option and the
value FLUME wrote at the last configure, for example
``-DERF_ENABLE_MPI=OFF`` with ``FLUME_MPI=ON`` on a first configure, stops
the configure step and names the FLUME option to set instead. A value equal
to the one FLUME wrote at the last configure is taken as FLUME's own and is
replaced when the FLUME option changes. The stand-alone build
(``FLUME_WITH_ERF=OFF``, the default) compiles no ERF code.

ERF's other options keep ERF's defaults. ERF's own tests
(``ERF_ENABLE_TESTS=ON``) also need ERF's ``googletest`` submodule, and add
ERF's suite, with ERF's labels, to FLUME's ``ctest``.

With CUDA, give ``CMAKE_CUDA_ARCHITECTURES``: FLUME enables CUDA before ERF
is added, so ERF's detection of the architecture on Cray machines is skipped. ERF's AMReX options
also make a device lambda that captures ``this`` a compile error.

The continuous integration builds and tests this layout with gcc, MPI and
double precision on the CPU. GPU builds with ERF are not tested there.

CMake options
-------------

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Option
     - Default
     - Effect
   * - ``FLUME_MPI``
     - ``ON``
     - Build with MPI.
   * - ``FLUME_PRECISION``
     - ``DOUBLE``
     - ``DOUBLE`` or ``SINGLE``.
   * - ``FLUME_GPU_BACKEND``
     - ``NONE``
     - ``NONE``, ``CUDA``, ``HIP`` or ``SYCL``. HIP and SYCL need the vendor
       compiler as ``CMAKE_CXX_COMPILER``; HIP also needs ``AMReX_AMD_ARCH``.
       For CUDA on a machine without a GPU, give ``CMAKE_CUDA_ARCHITECTURES``.
   * - ``FLUME_TESTING``
     - ``OFF``
     - AMReX's assertions and testing mode: new arrays are filled with
       signalling NaNs.
   * - ``FLUME_ENABLE_TESTS``
     - ``ON``
     - Register the tests with CTest. Needs Python 3.
   * - ``FLUME_WERROR``
     - ``OFF``
     - Warnings in FLUME's sources are errors. CPU builds only.
   * - ``FLUME_MPIEXEC_PREFLAGS``
     - empty
     - Extra ``mpiexec`` flags for the tests, for example
       ``--oversubscribe``.
   * - ``FLUME_WITH_ERF``
     - ``OFF``
     - Also build ERF (``erf_exec``) from ``FLUME_ERF_DIR`` as a subproject,
       on one AMReX. See `Building ERF with FLUME`_.
   * - ``FLUME_ERF_DIR``
     - ``Submodules/ERF``
     - ERF source tree. AMReX is taken from its ``Submodules/AMReX``. A
       folder without AMReX stops the configure step. With
       ``FLUME_WITH_ERF``, a folder that is not an ERF tree with the flood
       coupling (erf-model/ERF#5375 or later) also stops it.

A value outside the list of ``FLUME_PRECISION`` or ``FLUME_GPU_BACKEND`` stops
the configure step. FLUME builds without OpenMP; ``FLUME_OMP=ON`` stops the
configure step.

Moving the ERF submodule
------------------------

.. code-block:: bash

   git -C Submodules/ERF fetch --tags
   git -C Submodules/ERF checkout <tag-or-commit>
   git -C Submodules/ERF submodule update --init --depth 1 Submodules/AMReX
   git add Submodules/ERF

AMReX moves with ERF: rebuild and run the tests afterwards.

Building the documentation
--------------------------

.. code-block:: bash

   python -m pip install -r Docs/sphinx_doc/requirements.txt
   sphinx-build -W --keep-going -b html Docs/sphinx_doc Docs/sphinx_doc/_build/html

Warnings are errors, here and in the continuous integration.
