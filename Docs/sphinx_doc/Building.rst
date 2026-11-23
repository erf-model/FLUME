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
boxes one cell thick, so the same code can later link with ERF.

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
