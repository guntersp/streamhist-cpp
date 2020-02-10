# streamhist-cpp
header only streamhist port in C++



Overview
========

This project is C++ port of the streamhist library written in Python and 
an implementation of the streaming, one-pass histograms
described in Ben-Haim's `Streaming Parallel Decision
Trees <http://jmlr.org/papers/volume11/ben-haim10a/ben-haim10a.pdf>`__.
The histograms act as an approximation of the underlying dataset. The
histogram bins do not have a preset size, so as values stream into the
histogram, bins are dynamically added and merged as needed. One
particularly nice feature of streaming histograms is that they can be
used to approximate quantiles without sorting (or even individually
storing) values. Additionally, they can be used for learning,
visualization, discretization, or analysis. The histograms may be built
independently and merged, making them convenient for parallel and
distributed algorithms.

This ``C++`` version is a port of the ``Python`` `streamhist` library and its of the algorithm combines ideas and code from
`BigML <https://bigml.com>`__'s `Streaming Histograms for
Clojure/Java <https://github.com/bigmlcom/histogram>`__ and
`VividCortex <https://vividcortex.com>`__'s `Streaming approximate
histograms in Go <https://github.com/VividCortex/gohistogram>`__.



License
=======

* Copyright © 2015 Carson Farmer carsonfarmer@gmail.com
* Copyright © 2013 VividCortex
* All rights reserved. MIT Licensed.
* Copyright © 2013 BigML
* Licensed under the Apache License, Version 2.0


