:maintainers:
  `andrewtarzia <https://github.com/andrewtarzia/>`_

Overview
========

A public code base for using the `stk <https://stk.readthedocs.io/en/stable/>`_
ecosystem to construct and optimise cage molecules.

Each directory in ``src/makecages`` corresponds to a publication.


Installation per project
========================

Each project will have an ``environment.yml`` representing the ``conda`` or ``mamba`` environment 
used to run the script.

.. code-block:: bash

  mamba env create -f environment.yml

Then, you must install the library to access the modules of ``makecages``. So within the home directory of this repostory, run:

.. code-block:: bash

  pip install -e .


Studies
=======

- Added 04/09/2026: `study_3_azo`, associated with the paper `Inorganic Chemistry Frontiers , 2026, DOI: 10.1039/D6QI01683A <https://pubs.rsc.org/qi/article/doi/10.1039/D6QI01683A/1298566/Light-and-Metal-Controlled-Structural?searchresult=1>`_

