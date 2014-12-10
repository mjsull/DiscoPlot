DiscoPlot
=========

DiscoPlot allows the user to quickly identify genomic rearrangements, misassemblies and sequencing artefacts by providing a scalable method for visualising large sections of the genome. It reads single-end or paired read alignments in SAM, BAM or standard BLAST tab format and creates a scatter plot of opaque crosses representing the alignments to a reference. DiscoPlot is freely available (under a GPL license) for download (Mac OS X, Unix and Windows) at https://mjsull.github.io/DiscoPlot.

.. image:: https://pypip.in/version/DiscoPlot/badge.svg
        :target: https://pypi.python.org/pypi/DiscoPlot/
        :alt: Latest Version

.. image:: https://pypip.in/download/DiscoPlot/badge.svg
        :target: https://pypi.python.org/pypi/DiscoPlot/
        :alt: Downloads

.. image:: https://travis-ci.org/mscook/DiscoPlot.svg?branch=master
        :target: https://travis-ci.org/mjsull/DiscoPlot
        :alt: Build status

.. image:: https://landscape.io/github/mjsull/DiscoPlot/master/landscape.png
        :target: https://landscape.io/github/mjsull/DiscoPlot/master
        :alt: Code Health

Documentation
-------------

Please use this README.rst as the core DiscoPlot user documentation. 


News
----

**4/12/14:** Version 0.5.0 released.


Citation
--------

Cite this Github repository if you use DiscoPlot to generate figures 
for publications:: 

    SULLIVAN MJ, BEATSON SA^. 
    DiscoPlot - Visualising discordant readss.
    https://github.com/mjsull/DiscoPlot.

TODO:
 * Sam compatibility
 * Print selected read names, alignments or sequences


Installation
------------

DiscoPlot is a commandline application. If you're not familiar with the 
commandline we recommend you ask local IT support to help you install it.

We now test DiscoPlots builds on both Linux (Ubuntu >= 12.04) and Windows 8 system. 

You will need to install/have installed:
    * python >= 2.7 (**Python 3 is not supported**)

To automatically generate BLAST aligments (For long read DiscoPlots) using DiscoPlot you will need to install/have installed:
    * ncbiblast >= 2.2.27
    
You can check these are installed by::
    
    $ python --version
    $ blastn -version

Installation of python or blastn (without a package manager) is beyond the 
scope of this document.

If you have both python and blastn you need to (if not already present) 
install pip_.

You can check if pip_ exists with::

    $ which pip

If you get a "not found", please read the `pip installation instructions`_. 

**If you already have pip we do suggest you upgrade it.** We are using version 
1.5.6 at the time of writing this document. 

You can upgrade pip_ like this::

    $ pip install --upgrade pip


The following python libraries_ should be installed (automatically) if you follow 
the installation instructions detailed below.

We use the following python libraries_:
    * numpy >= 1.6.1
    * matplotlib >= 1.1.0
    * pysam >= 0.8.1

Pysam is only required for generating DiscoPlots with BAM files. SAM compatability has been included to allow windows users to generate DiscoPlots. PySam will not install on Windows, don't bother trying (or if you've succeeded please let me know how).

These libraries will also have dependencies (i.e. atlas, lapack, fortran 
compilers, freetype and png). **These most likely won't be installed on 
your computer. Please install these before attempting the installation.**

Linux (Ubuntu)
~~~~~~~~~~~~~~

SeqFindr uses 3rd party packages that are extremely important for scientific 
computing but are notoriously difficult to install. While *pip install * 
*--user DiscoPlot* may work we recommend you install these 3rd party packages 
using apt-get.

Run::

    $ sudo apt-get install python-numpy python-matplotlib python-dev libatlas-dev liblapack-dev gfortran libfreetype6-dev libfreetype6 libpng-dev 

Now pip_ install SeqFindr::
    
    $ pip install --user DiscoPlot

We use the --user option of pip_ to put SeqFindr in: /home/$USER/.local/bin/
You need to add this location to you ~/.bash_profile. 

Add SeqFindr to your path::

    $ echo 'export PATH=$PATH:/home/$USER/.local/bin/' >> ~/.bash_profile

Finally install BLAST+::

    $ sudo apt-get install ncbi-blast+ 

**Test it:**

Run::
    
    $ SeqFindr -h 
    $ python -c 'import SeqFindr; print SeqFindr'


MacOSX (Mavericks)
~~~~~~~~~~~~~~~~~~

**You'll need to have the equivalents of python-dev libatlas-dev liblapack-dev 
gfortran libfreetype6-dev libfreetype6 & libpng-dev installed.** We had no 
problems installing SeqFindr on a recently acquired OSX Mavericks machine 
using the homebrew package manager.

The installed packages on this machine via::

    $ brew list 

Are available at this gist_.

pip install DiscoPlot::
    
    $ pip install --user DiscoPlot

We use the --user option of pip_ to put SeqFindr in: /home/$USER/.local/bin/
You need to add this location to you ~/.bash_profile. 

Add DiscoPlot to your path::

    $ echo 'export PATH=$PATH:/home/$USER/.local/bin/' >> ~/.bash_profile

Finally install BLAST+::

    $ sudo brew install blast 

**Test it:**

Run::
    
    $ DiscoPlot -h 
    $ python -c 'import DiscoPlot; print DiscoPlotr'


Upgrading DiscoPlot
~~~~~~~~~~~~~~~~~~

You can upgrade like this::
    
    pip install --upgrade DiscoPlot


**Please regularly check back to make sure you're running the most recent 
DiscoPlot version.**



Example figure produced by DiscoPlot
-----------------------------------

DiscoPlot of a mock genome. A mock genome was created by adding genomic rearrangements to the chromosome of E. coli str. UTI89.  Paired-end reads generated from the mock genome (query) with GemSim (ref) and mapped back to UTI89 (reference). The first ~500 Kbp were then visualised using DiscoPlot.

.. image:: https://raw.github.com/mscook/SeqFindr/master/example/CU_fimbriae.png
    :alt: SeqFindr CU fimbriae genes image
    :align: center

DiscoPlots of common structural variants. Each box shows a common genomic rearrangement represented by a DiscoPlot. Rows A and B were created using 100 bp long paired-end reads with an insert size of 300bp. Rows C and D were created using single-end reads with an average length of 1000bp. For each box the rearrangement in the sequenced genome is listed, followed by the scale of the gridlines in brackets. A1,  C1: 300 bp deletion (400 bp). A2, C2: 300 bp insertion (400 bp). A3, C3: 300 bp inversion (400 bp). A4, C4: 300 bp sequence translocated 50 Kbp upstream (10 Kbp). B1, D1: 3000 bp deletion (1000 bp). B2, D2: 3000 bp insertion (500 bp). B3, D3: 3000 bp inversion (1000 bp). B4, D4: 3000 bp sequence translocated 50 Kbp upstream (10 Kbp). C1) 

.. image:: https://raw.github.com/mscook/SeqFindr/master/example/CU_fimbriae.png
    :alt: SeqFindr CU fimbriae genes image
    :align: center


Tutorial
--------

We provide a script_ to run all the examples. **Note:** We have changed the 
color generation code. As a consequence the background colors will be 
different when running this yourself. The results will not change.

Navigate to the SeqFindr/example directory (from git clone). The following files should be present:
    * A database file called *Antibiotic_markers.fa* 
    * An ordering file called *dummy.order* (-i option)
    * An assemblies directory containing *strain1.fa, strain2.fa and strain3.fa*
    * A consensus directory containing *strain1.fa, strain2.fa and strain3.fa*
      (-m option)

**Note:** the assembly and consensus directories contain:
    * the same number of files (3 each)
    * there is a 1-1 filename mapping (strain1.fa, strain2.fa, strain3.fa == 
      strain1.fa, strain2.fa, strain3.fa)
    * there are only fasta files. If you wish to include complete genomes 
      either download the genomes in fasta format OR convert the Genbank or 
      EMBL files to fasta format. 

The toy assemblies and consensuses were generated such that:
    * **strain1** was missing: 70-shv86, 70-ctx143 and 70-aac3(IV)380 with 
      mis-assembly of 70-aphA(1)1310 & 70-tem8674
    * **strain2** was missing: 70-oxa(7)295, 70-pse(4)348 70-ctx143, 
      70-aadA1588, 70-aadB1778 and 70-aacC(2)200
    * **strain2** was missing 70-shv86, 70-ctx143 and 70-aac3(IV)380 with 
      mis-assembly of 70-aphA(1)1310, 70-tem8674 and 70-aadA1588


Running all the examples at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Something like this::

    $ # Assuming you git cloned, python setup.py install
    $ cd SeqFindr/example
    $ ./run_examples.sh
    $ # See directories run1/ run2/ run3/ run4/


Run 1 - Looking at only assemblies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Command::

    SeqFindr Antibiotic_markers.fa assemblies/ -o run1 -l 

.. image:: https://raw.github.com/mscook/SeqFindr/master/example/run1_small.png
    :alt: run1
    :align: center


Link to full size run1_.


Run 2 - Combining assembly and mapping consensus data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Command::

    SeqFindr Antibiotic_markers.fa assemblies/ -m consensus/ -o run2 -l

.. image:: https://raw.github.com/mscook/SeqFindr/master/example/run2_small.png
    :alt: run2
    :align: center


Link to full size run2_.


Run 3 - Combining assembly and mapping consensus data with differentiation between hits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Command::

    SeqFindr Antibiotic_markers.fa assemblies/ -m consensus/ -o run3 -l -r

.. image:: https://raw.github.com/mscook/SeqFindr/master/example/run3_small.png
    :alt: run3
    :align: center


Link to full size run3_.


The clustering dendrogram looks like this:

.. image:: https://raw.github.com/mscook/SeqFindr/master/example/dendrogram_run3_small.png
    :alt: run3 dendrogram
    :align: center


Link to full size dendrogram_.


Run 4 - Combining assembly and mapping consensus data with defined ordering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Note:** the ordering file is defined using the option *--index_file*. The 
ordering file **must** contain the same number of strains as the assemblies 
directory and the strain names must agree (TODO - add a script to flag issues).

Command::

    SeqFindr Antibiotic_markers.fa assemblies/ -m consensus/ -o run4 -l -r --index_file dummy.order

.. image:: https://raw.github.com/mscook/SeqFindr/master/example/run4_small.png
    :alt: run4
    :align: center


Link to full size run4_.


How to generate mapping consensus data
--------------------------------------

**We strongly recommend that you use mapping consensus data.** It minimises 
the effects of missassembly and collapsed repeats.

We use Nesoni_. We use the database file (in multi-fasta format) as the 
reference for mapping. Nesoni_ has no issues with multifasta files as 
references (BWA will treat them as separate chromosomes). 
The workflow is something like this::

    $ nesoni make-reference myref ref-sequences.fa
    $ # for each strain
    $ #     nesoni analyse-sample: mysample myref pairs: reads1.fastq reads2.fastq
    $ #     extract the consensus.fa file


For those of you using a cluster running PBSPro see:
https://github.com/mscook/SeqFindr_nesoni
This is a script that generates a job array, submits and cleans up the
mapping results ready for input to SeqFindr.

The output from the described workflow and SeqFindr_nesoni is a consensus.fa 
file which we term the mapping consensus. This file is a multi-fasta file of 
the consensus base calls relative to the database sequences.

Caveats: 
    * you will probably want to allow multi-mapping reads (giving *--monogamous
      no --random yes* to nesoni consensus) (this is default for
      SeqFindr_nesoni), 
    * The (poor) alignment of reads at the start and the end of the database 
      genes can result in N base calls. This can result in downstream false 
      negatives.

**SeqFindr now provides a solution to minimise the effects of poor mapping at 
the start and end of the given sequences.** 

The SeqFindr option is -s or --STRIP::

    -s STRIP, --strip STRIP Strip the 1st and last N bases of mapping consensuses & database [default = 10]

By default this strips the 1st and last 10 bases from the mapping consensuses. 
We have had good results with this value. Feel free to experiment with 
different values (say, -s 0, -s 5, -s 10, -s 15). Please see image-compare_ 
a script we developed to compare the effects of different values of -s on the 
resultant figures. 


SeqFindr usage options
----------------------

See the help listing_. You can get this yourself with::

    $ SeqFindr -h


Future
------

Please see the TODO_ for future SeqFindr project directions.





.. _pip: http://www.pip-installer.org/en/latest/
.. _libraries: https://github.com/mscook/SeqFindr/blob/master/requirements.txt
.. _image-compare: https://github.com/mscook/image-compare
.. _listing: https://github.com/mscook/SeqFindr/blob/master/HELP.rst
.. _changelog: https://github.com/mscook/SeqFindr/blob/master/CHANGES.rst
.. _TODO:  https://github.com/mscook/SeqFindr/blob/master/TODO.rst
.. _script: https://raw.github.com/mscook/SeqFindr/master/example/run_examples.sh
.. _run1: https://raw.github.com/mscook/SeqFindr/master/example/run1.png
.. _run2: https://raw.github.com/mscook/SeqFindr/master/example/run2.png
.. _run3: https://raw.github.com/mscook/SeqFindr/master/example/run3.png
.. _dendrogram: https://raw.github.com/mscook/SeqFindr/master/example/dendrogram_run3.png
.. _run4: https://raw.github.com/mscook/SeqFindr/master/example/run4.png
.. _Nesoni: http://www.vicbioinformatics.com/software.nesoni.shtml
.. _SeqFindr documentation: http://seqfindr.rtfd.org
.. _SeqFindr official site: http://mscook.github.io/SeqFindR/
.. _gist: https://gist.github.com/mscook/ef7499fc9d2138f17c7f
.. _pip installation instructions: http://pip.readthedocs.org/en/latest/installing.html
