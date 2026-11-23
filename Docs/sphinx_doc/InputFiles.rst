Input files
===========

A relative file name is taken from the directory the run starts in. A name
with spaces goes in double quotes. Each message about a file names the key,
the file and the position of the fault.

.. _sec:RasterFiles:

Raster files
------------

``flume.bed.file``, ``flume.buildings.file`` and
``flume.friction.manning_file`` name ESRI ASCII grids (``.asc``), which QGIS
and ArcGIS write:

.. code-block:: text

   ncols         4
   nrows         3
   xllcorner     500000.0
   yllcorner     4100000.0
   cellsize      10.0
   NODATA_value  -9999
   12.1 12.4 12.9 13.3
   11.8 12.0 12.5 12.8
   11.6 11.7 12.1 12.4

.. list-table:: Header keywords (any case, any order, each once)
   :header-rows: 1
   :widths: 34 66

   * - Keyword
     - Meaning
   * - ``ncols``, ``nrows``
     - Whole numbers from 1 to 2\ :sup:`30`.
   * - ``xllcorner``, ``yllcorner``
     - Lower left corner; or ``xllcenter``, ``yllcenter``, the centre of
       the lower left cell.
   * - ``cellsize``
     - Cell size, greater than 0; or ``dx`` and ``dy`` for cells that are
       not square.
   * - ``NODATA_value``
     - Optional marker of missing data.

Rules:

* The values follow, ``ncols`` per row, the northern row first: exactly
  ``ncols`` times ``nrows`` finite numbers.
* The raster is on the coordinates of the run, in metres, and holds the
  whole domain. FLUME does not reproject or change units.
* A cell that reads a ``NODATA_value`` cell stops the run.
* The bed is interpolated bilinearly between raster cell centres
  (:ref:`sec:BedRaster`). Buildings and the Manning coefficient take the
  raster cell that holds the cell centre; a centre on a raster cell edge
  takes the cell above it.
* On a periodic domain the raster is not wrapped.
* The file may be much larger than the domain: only the part over the
  domain is kept. One rank reads the file, so it must fit in the memory of
  one rank.

.. _sec:HyetographFiles:

Hyetograph files
----------------

``flume.rain.file`` holds lines of a time in s and a rain rate in mm/h:

.. code-block:: text

   # time (s)  rate (mm/h)
   0      36
   600    72   # the rate doubles
   1200   0

Rules:

* ``#`` starts a comment. Blank lines are skipped.
* Each line holds two finite numbers; the rate is at least 0; the times
  increase strictly. At least one line. Times may be negative.
* Before the first time the rate is the first one, and after the last time
  the last one. A storm that starts after the run starts needs a first line
  of rate 0; a storm that ends needs a last line of rate 0. The run prints a
  note when the first time is after 0 with a rate above 0.

.. _sec:RainIndexFiles:

Gridded rain index files
------------------------

With ``flume.rain.type = gridded``, ``flume.rain.file`` holds lines of a
time in s and the path of an ESRI ASCII grid of the rain rate in mm/h:

.. code-block:: text

   # time (s)  grid of the rate (mm/h)
   0      radar_0000.asc
   300    radar_0300.asc   # five minutes later
   600    radar_0600.asc

Rules:

* ``#`` starts a comment. Blank lines are skipped. Each line holds exactly
  two items, so a grid path holds no space and no ``#``. The times are
  finite and increase strictly. At least one line.
* A relative grid path is taken from the directory of the index file.
* Each grid has its own cells. It must hold the whole domain, hold no
  ``NODATA_value`` in a cell that shares area with the domain, and hold no
  rate below 0 anywhere. Every grid is read and checked at start-up; a
  failed check names the key, the index file, the entry and the time.
* Before the first time the first grid holds, and after the last time the
  last one.
* During the run a grid is read again when a step needs it. The grids a
  later step can use are kept (one with ``step``, two with ``linear``); a
  grid only the current step uses is read and freed at once. A grid of no
  rain is not read again. The files must not change during the run.

.. _sec:Precision:

Single precision and large coordinates
--------------------------------------

The domain corners are stored in the precision of the build. In single
precision a corner at UTM-size coordinates is rounded to the nearest
representable number, 0.5 m apart at 4.4e6 m, and the cell size changes with
it. Use double precision, or corners near the origin, when this matters.

Raster files are read in double precision, and the position of each cell in
the raster is formed from differences of corners in double precision, so map
coordinates cost no digits in the sampling.

An elevation :math:`z` is held to about :math:`6 \times 10^{-8}\,|z|` in
single precision (0.24 mm at 4000 m). A depth taken as the difference of two
elevations is no finer than that.
