# globalETAS
globally indexed ETAS (and other) earthquake models. the global indexing should be applicable to numerous applications.

This repository contains all the base code for globalETAS. Most of the salient code is in the globalETAS.py module; there are some scripts, sample implmentations, and similar collections in various *.py modules. Some support codes are contained in one or more "submodule". These repositories can be found as independent repositories alongside globalETAS; they have been imported as submodules for this project (see clone instructions below).

In addition to *.py modules, there are a bunch of *.ipynb "jupyter notebooks," the modern answer to ipython notebooks. These do/will include some unit tests, demo scripts, and sample implementations/reports for several significant earthquakes over the past several years. Perticularly, the newer notebooks can serve as a nice template for a new earthquake of interest; as with everything else, these are a work in progress.

TODOs: An incomplete list
- There may be issues with notebooks and multi-processing. I believe I have found instances where mpp data do not appear to pipe back correctly to the parent process, but this needs to be confirmed. Most folks familiar with Python MPP and notebooks will say it's a twitchy combination. In particular, during development, if your processes are crashing, your notbook is probably crashing too. Be prepared to restart your service.
- Catalogs: 
 - There seems to be no good solution right now. I used to use a web API hack of the ANSS catalog, but that was abandoned for USGS's ComCat, which was awesome for a couple of years (just `conda install libcomcat`), but now seems to be not terribly well supported, so it does not easily install on newer Python, and particularly in combination with some other packages. 
 - So you can install comcat it en an `environment`, but again, there can be some difficulty playing with other packages. Also, it's really big, so may not play well in user-HPC enviornments (where `$HOME` spase is limited).
 - I wrote a little hack to the comcat web API which works really well, except that it is limited to small subcatalogs, so a major TODO is to write a composite wrapper to break large queries (how do we know -- in advance, that a query is large?) into small queries and then put them together.

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

