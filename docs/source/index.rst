.. pyro documentation main file, created by
   sphinx-quickstart on Mon Dec 25 18:42:54 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

***************************
Simple Astronomy Animations
***************************

A collection of astronomy animations produced for my intro astronomy classes.  These are
all written in python (using matplotlib) and the videos are rendered using FFmpeg.

* :ref:`earthmoon` : demonstrations of the orbit and phases of the moon and Earth's seasons.

* :ref:`ellipses` : animations showing the properties of ellipses

* :ref:`orbits` : demonstrations of different orbits a satellite can make around Earth.

* :ref:`solarsystem` : Kepler's laws, parallax, retrograde motion, and more

* :ref:`binary` : properties of binary star systems and exoplanet systems

* :ref:`thermo` : blackbody spectrum, random walk, thermal motion

* :ref:`waves` : properties of waves, including Doppler effect

* :ref:`nuclear` : nuclear physics concepts

.. note::

   All the source code for these animations is available on github: https://github.com/zingale/astro_animations

.. note::

   Videos are rendered with `ffmpeg <https://www.ffmpeg.org/>`_ using:

   .. code:: bash

      ffmpeg -framerate 15 -f image2 -pattern_type glob -i "*.png" -vcodec mpeg4 -c:v libx264 -crf 20 -pix_fmt yuv420p -movflags +faststart movie.mp4


.. toctree::
   :maxdepth: 1
   :caption: Animation Library
   :hidden:

   earth_moon
   ellipse
   mechanics
   solar_system
   binary
   thermo
   waves
   nuclear
