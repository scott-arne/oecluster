Python API Reference
====================

The main Python guide is :doc:`../python-api`. This page adds generated
reference material for the importable Python modules.

Top-Level Package
-----------------

The public API lives directly on the top-level :mod:`oecluster` package:
distance computation (:func:`oecluster.pdist`, :func:`oecluster.cdist`),
clustering (:func:`oecluster.butina`, :func:`oecluster.dbscan`,
:func:`oecluster.hdbscan`, :func:`oecluster.agglomerative`,
:func:`oecluster.k_medoids`, :func:`oecluster.bitbirch`,
:func:`oecluster.murcko`), representative
selection, quality reporting, and the supporting option, result, and
comparison classes.

.. automodule:: oecluster
   :members:
   :undoc-members:

CLI Wrapper
-----------

.. automodule:: oecluster._cli
   :members:
   :undoc-members:
