# globalETAS
globally indexed ETAS (and other) earthquake models. the global indexing should be applicable to numerous applications.

This repository contains all the base code for globalETAS. Most of the salient code is in the globalETAS.py module; there are some scripts, sample implmentations, and similar collections in various *.py modules. Some support codes are contained in one or more "submodule". These repositories can be found as independent repositories alongside globalETAS; they have been imported as submodules for this project (see clone instructions below).

In addition to *.py modules, there are a bunch of *.ipynb "jupyter notebooks," the modern answer to ipython notebooks. These do/will include some unit tests, demo scripts, and sample implementations/reports for several significant earthquakes over the past several years. Perticularly, the newer notebooks can serve as a nice template for a new earthquake of interest; as with everything else, these are a work in progress.

Some things to look into: There may be issues with notebooks and multi-processing. I believe I have found instances where mpp data do not appear to pipe back correctly to the parent process, but this needs to be confirmed.

Cloning:
 > git clone --recursive https://github.com/markyoder/globalETAS.git
 
 If you forget the --recursive, you won't get the submodules. fix it like:
 
> git submodule update --init --recursive
which tells git to run two commands:

> git submodule init
... and subsequently, 
> git submodule update

We're working on making this happen as smoothly as possible, but it may still be necessary to install a few extra bits. We recommend using Anaconda Python 3.x; for those stubborn amongst us who insist to run on their system Python, you an probably just replace "conda" installations with "pip":

- on a fresh linux install... stuff we have to do besides just clone this (NOTE: you may or may not need to specify the chanel, `-c conda-forge`:
   - `pip install geopy`
   - `conda install basemap`
   - `conda install -c conda-forge basemap-data-hires`
   - `pip install geographiclib`
   - `conda install -c ioos rtree`
      ** newer versions of Anaconda (and other Python 3.6) may require an updated rtree library; try:
     conda install -c conda-forge rtree
      ** in general, you can search for an appropriate installation with:
          anaconda search -t conda rtree

Managing dependencies remains a work in progress. Please contact the primary author(s) with any questions or problems; we will do our best to identify an demploy appropriate dependency libraries, code around them, etc. If rtree libraries become problematic, we will push development of other indexing objects, such as the (under development) 'bindex' model.

