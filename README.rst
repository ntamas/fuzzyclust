Algorithm for fuzzy community detection in complex networks
===========================================================

:Author: Tamas Nepusz
:Version: 0.1

This is a possible implementation of the fuzzy community detection algorithm
published in [1]_. It is not the same as the one we used in [1]_ - that copy
was lost in a hard drive crash some time ago, but it should work similarly.

The algorithm is mostly suitable for small and mid-size networks due to its
computational cost and the shape of the goal function landscape that has
many local minima. The algorithm tries its best to climb out of those.

Requirements
------------

- igraph_ 0.6 or later (the C library of course).

- ``gcc`` and ``make`` of course.

.. _igraph: http://igraph.org/c/

Usage
-----

See the manpage in ``doc/fuzzyclust.1`` for more information.

Bugs, questions
---------------

Send them to me here on GitHub or to ``nepusz at hal elte hu``.

.. [1] Nepusz T, Petróczi A, Négyessy L, Bazsó F: Fuzzy communities and
       the concept of bridgeness in complex networks. *Physical Review E*
       77:016107, 2008.