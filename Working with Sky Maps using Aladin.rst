
Working with Sky Maps using Aladin Desktop
==========================================


The tutorial focuses on some basic strategies for working with
gravitational-wave sky localizations in the context of electromagnetic
follow-up activities. Here we propose the usage of `Aladin
Desktop <https://aladin.u-strasbg.fr/java/nph-aladin.pl?frame=downloading>`__.
The following main topics are addressed.

1. Sky map visualizations and crebible regions.
2. Access to existing catalogs.
3. Sky maps comparisons.

MOC and GW sky localizations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The contours of a GW sky localization, which enclose a given percentage
of the total probability, are constructed using a water-filling
algorithm: the pixels from most probable to least are ranked, and summed
up to get a fixed level of probability. The enclosed area within a given
probability level contour of a GW sky map can be effectively described
through the Multi-Order Coverage (MOC) method. It is a standard of the
Virtual Observatory which provides a multi-scale mapping based on
HEALPix sky tessellation. Basically, the algorithm maps irregular and
complex sky regions into hierarchically grouped predefined cells. Each
MOC cell is defined by two numbers: the hierarchy level (HEALPIX ORDER)
and the pixel index (HEALPIX NPIX).The NUNIQ scheme defines an algorithm
for packing an (ORDER, NPIX) pair into a single integer for compactness:

.. math:: uniq = 4\times 4^{(order)} + npix

 MOCs are serialized as ``FITS`` or ``JSON`` files. The MOC resolution
is determined by the map resolution parameter :math:`N_{side}`, which is
used for defining the resolution of the grid. For a typical
:math:`N_{side}` = 512, the MOC order resolution is :math:`order` = 9;
(:math:`N_{side}` = 2\ :math:`^{order}`). The MOC maps make database
queries for retrieving objects, logical operation ( such as union,
intersection, subtraction, difference) extremely simple and fast even
for very complex sky regions. If databases are adapted to support MOC
based queries, such as
`VizieR <http://vizier.u-strasbg.fr/viz-bin/VizieR>`__, they offer a
useful method allowing any support of sky region query.

cite

[1] Pérez, Fernando and Granger, Brian E.. 2007. *IPython: a System for
Interactive Scientific Computing*. `URL <http://ipython.org>`__





1. Running Aladin Desktop
-------------------------

`Aladin
Desktop <https://aladin.u-strasbg.fr/java/nph-aladin.pl?frame=downloading>`__
is developed in Java. As any Java tool, Aladin Desktop requires a `Java
Virtual Machine <https://www.java.com/en/>`__ on your machine. Download
the ``Aladin.jar`` file from
`here <https://aladin.u-strasbg.fr/java/Aladin.jar>`__ and execute it
from a terminal by typing

::

                            $ java -Xmx2g -jar Aladin.jar           

The flag ``-Xmx< ammount of memory >`` specifies the maximum memory
allocation pool for a Java Virtual Machine (JVM). Here 2GB of memory is
allocated. For GW sky localizations with ``NSIDE = 2048``, increase the
memory allocated up to 3GB, ``-Xmx3g``.

2. Loading a GW sky localization
--------------------------------

You can *copy&paste* the sky map location from the
`GraceDB <https://gracedb.ligo.org/>`__ or drag the local file in the
main Aladin window. Aladin recognizes only the sky maps with extention
**.fits.gz**. For a deep discussion about the skymap formats supported
from LVC, see the Section 3 of the User Guide.

3. Building a Credible Region
-----------------------------

The sequence of the Aladin GUI (Graphical User Interface) commands to
create a credible region at a defined confidence level is reported
below. From the main menu press
``→ Coverage → Generate a MOC based on... → The current probability skymap → MOC generation window``

The ``MOC generation`` window requires two mandadory parameters 1) the
*probability sky map* and 2) the *threshold*. The probability skymap
entry - that means the GW sky localization in our contex - can be
selected from the drop down menu. They are the GW sky localization
previously loaded and showed in the Aladin Stack. The threshold input
represents the percentage of the credible region passed in decimal
(ranging from 0 to 1). Press the ``CREATE`` button to generate the
resulting credible region. The credible region is created and loaded in
the Aladin Stack. The credible region obtained so far are decoded in the
Multi Order Coverage map (MOC) and each new level can be independently
used.

4. *Properties* Window
----------------------

The associated ``Properties`` windows allows to change the *drawing
methods* in *perimeter* in order to simultaneously visualize multiple
confidence levels. To open the ``Properties`` window, right click on the
selected plan in the Aladin stack. This operation facilitates tiling
operations by telescopes monitoring the highest probability areas. The
enclosed sky area in square degrees and the percentage of the sky
coverage are quoted for each credible region either **i)** by leaving
the cursor on the corresponding plan loaded in the Aladin stack or
**ii)** by opening the associated Properties windows.

You can overlap a large dataset of image backgrounds provided by the
`HiPS list aggregator <https://aladin.unistra.fr/hips/list>`__ or you
can generate your own HiPS from image/cube data. For doing this, from
the main menu press
\`\ ``→ Tool → Generate a HIPS based on... → An image collections (FITS, JPEG, PNG)...``

5. Querying and Filtering a Galaxy Catalog
------------------------------------------

Singer et al. (2016) [4] discuss a fast algorithm for obtaining a
three-dimensional probability estimates of sky location and luminosity
distance from observations of binary compact object mergers with
Advanced LIGO and Virgo. Combining the reconstructed gravitational wave
volumes with positions and redshifts of possible host galaxies provides
a manageable list of sky location targets to search for the
electromagnetic counterpart of the gravitational wave signal. These
tasks can efficiently be performed in the Aladin Desktop using the data
collections tree and the filter methods as follows.

``Aladin data collections tree → Select → click on the catalog item → in the popup menu check → by region & MOC``

Now we can filter the galaxy catalog...

``Catalog → Create a filter → Properties window → Advanced mode→ Or enter your filter definition:``

An example about the Aladin filter using as galaxy selection the
marginal distance posterior distribution integrated over the whole sky
is reported below:
``${Dist} > DISTMEAN-DISTSTB && ${Dist} < DISTMEAN+DISTSTB {draw}``. The
posterior mean distance (Mpc) and the posterior standard deviation of
distance (Mpc) are reported in the fits file header with the keywords
``DISTMEAN`` and ``DISTSTB``.

6. Sky map comparisons
----------------------

References
----------






.. raw:: html

   <!--bibitex
   @Article{PER-GRA:2007,
     Author    = {P\'erez, Fernando and Granger, Brian E.},
     Title     = {{IP}ython: a System for Interactive Scientific Computing},
     Journal   = {Computing in Science and Engineering},
     Volume    = {9},
     Number    = {3},
     Pages     = {21--29},
     month     = may,
     year      = 2007,
     url       = "https://ipython.org",
     ISSN      = "1521-9615",
     doi       = {10.1109/MCSE.2007.53},
     publisher = {IEEE Computer Society},
   }
   -->

