Building
========

Requirements
------------

* CMake 3.25 or newer (required by the pinned AMReX).
* A C++20 compiler.
* MPI (optional; on by default).
* Python 3 for the tests. Without it the configure step stops, unless
  ``-DFLUME_ENABLE_TESTS=OFF`` is given.

Getting the code
----------------

FLUME needs two submodules: ERF, and the AMReX that ERF pins. ERF's other
submodules (RRTMGP, Noah-MP, WW3, EKAT, GoogleTest) are not used by the
stand-alone build and need not be fetched.

.. code-block:: bash

   git clone https://github.com/erf-model/FLUME.git
   cd FLUME
   git submodule update --init Submodules/ERF
   git -C Submodules/ERF submodule update --init --depth 1 Submodules/AMReX

``git clone --recursive`` also works. It fetches every submodule of ERF.

FLUME has no AMReX pin of its own. AMReX is built from
``Submodules/ERF/Submodules/AMReX``. The stand-alone build compiles no ERF code.

Configure, build and test
-------------------------

.. code-block:: bash

   cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
   cmake --build build-rel -j 8
   ctest --test-dir build-rel --output-on-failure

The executable is ``build-rel/flume``. Without ``CMAKE_BUILD_TYPE`` the build
is ``Release``. Name each build directory
``build-<something>``; those names are ignored by git.

The flood grid is two-dimensional. AMReX is built in three dimensions and the
boxes are one cell thick, so the same code can later link with ERF.

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
   * - ``FLUME_OMP``
     - ``OFF``
     - Build with OpenMP. Cannot be combined with a GPU backend.
   * - ``FLUME_PRECISION``
     - ``DOUBLE``
     - Floating-point precision: ``DOUBLE`` or ``SINGLE``.
   * - ``FLUME_GPU_BACKEND``
     - ``NONE``
     - GPU backend: ``NONE``, ``CUDA``, ``HIP`` or ``SYCL``. HIP and SYCL need
       the vendor compiler as ``CMAKE_CXX_COMPILER``; HIP also needs
       ``AMReX_AMD_ARCH``. For CUDA on a machine without a GPU, give
       ``CMAKE_CUDA_ARCHITECTURES``: AMReX's default is the GPU of the build
       machine.
   * - ``FLUME_TESTING``
     - ``OFF``
     - Turn on AMReX's assertions and its testing mode, in which new arrays are
       filled with signalling NaNs.
   * - ``FLUME_ENABLE_TESTS``
     - ``ON``
     - Register the tests with CTest. Needs Python 3.
   * - ``FLUME_WERROR``
     - ``OFF``
     - Treat warnings in FLUME's own sources as errors. CPU builds only: with
       a GPU backend FLUME adds no warning flags and this option has no
       effect. The continuous integration turns it on.
   * - ``FLUME_MPIEXEC_PREFLAGS``
     - empty
     - Extra ``mpiexec`` flags for the tests, for example ``--oversubscribe``.
       With MPI the tests need a launcher (``MPIEXEC_EXECUTABLE``).

A value of ``FLUME_PRECISION`` or ``FLUME_GPU_BACKEND`` outside its list stops
the configure step.

Moving the ERF submodule
------------------------

.. code-block:: bash

   git -C Submodules/ERF fetch --tags
   git -C Submodules/ERF checkout <tag-or-commit>
   git -C Submodules/ERF submodule update --init --depth 1 Submodules/AMReX
   git add Submodules/ERF

AMReX moves with ERF, so rebuild and run the tests after moving the pin.

Building the documentation
--------------------------

.. code-block:: bash

   python -m pip install -r Docs/sphinx_doc/requirements.txt
   sphinx-build -W --keep-going -b html Docs/sphinx_doc Docs/sphinx_doc/_build/html

Warnings are errors, here and in the continuous integration.
