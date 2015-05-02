DiscoPlot
=========

.. image:: https://pypip.in/version/DiscoPlot/badge.svg
        :target: https://pypi.python.org/pypi/DiscoPlot/
        :alt: Latest Version

.. image:: https://pypip.in/download/DiscoPlot/badge.svg
        :target: https://pypi.python.org/pypi/DiscoPlot/
        :alt: Downloads

.. image:: https://travis-ci.org/BeatsonLab-MicrobialGenomics/DiscoPlot.svg?branch=master
        :target: https://travis-ci.org/BeatsonLab-MicrobialGenomics/DiscoPlot
        :alt: Build status


DiscoPlot allows the user to quickly identify genomic rearrangements, 
misassemblies and sequencing artefacts by providing a scalable method for 
visualising large sections of the genome. It reads single-end or paired read 
alignments in SAM, BAM or standard BLAST tab format and creates a scatter 
plot of opaque crosses representing the alignments to a reference. 
DiscoPlot is freely available (under a GPL license) for download (Mac OS 
X, Unix and Windows) at: 
https://github.com/BeatsonLab-MicrobialGenomics/DiscoPlot/releases.

.. image:: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_3.png
        :target: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_3.png
        :alt: DiscoPlot figure

**DiscoPlot of a mock genome.** A mock genome was created by adding genomic 
rearrangements to the chromosome of E. coli str. UTI89.  Paired-end reads 
generated from the mock genome (query) with GemSim and mapped back to UTI89 
(reference). The first ~500 Kbp were then visualised using DiscoPlot.


Documentation
-------------

Please use this README.rst as the core DiscoPlot user documentation. 


Citation
--------

Cite this Github repository if you use DiscoPlot to generate figures 
for publications:: 

    SULLIVAN MJ & BEATSON SA^. 
    DiscoPlot - Visualising discordant reads.
    https://github.com/BeatsonLab-MicrobialGenomics/DiscoPlot.


TODO
----

On the roadmap:
    * Sam compatibility
    * Print selected read names, alignments or sequences


Installation
------------

DiscoPlot is a commandline application. If you're not familiar with the 
commandline we recommend you ask local IT support to help you install it.

You will need to install/have installed:
    * python >= 2.7 (**Python 3 is not supported**)

To automatically generate BLAST aligments (For long read DiscoPlots) using DiscoPlot you will need to install/have installed:
    * ncbiblast+ >= 2.2.27
    
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
    * numpy >= 1.8.1
    * matplotlib >= 1.3.1
    * pysam >= 0.8.1

Pysam is only required for generating DiscoPlots with BAM files. SAM 
compatability has been included to allow windows users to generate 
DiscoPlots. PySam will not install on Windows, don't bother trying (or if 
you've succeeded please let me know how).


Linux (Ubuntu)
~~~~~~~~~~~~~~

Discoplot uses 3rd party packages that are extremely important for scientific 
computing but may be difficult to install. While *pip install * 
*--user DiscoPlot* may work we recommend you install these 3rd party packages 
using apt-get.

Run::

    $ sudo apt-get install python-numpy python-matplotlib 

Now pip_ install DiscoPlot::
    
    $ pip install --user DiscoPlot

We use the --user option of pip_ to put DiscoPlot in: /home/$USER/.local/bin/
You need to add this location to you ~/.bash_profile. 

Add DiscoPlot to your path::

    $ echo 'export PATH=$PATH:/home/$USER/.local/bin/' >> ~/.bash_profile

Finally install BLAST+::

    $ sudo apt-get install ncbi-blast+ 


MacOSX (Mavericks)
~~~~~~~~~~~~~~~~~~

**You'll need to have the equivalents of python-dev libatlas-dev liblapack-dev 
gfortran libfreetype6-dev libfreetype6 & libpng-dev installed.** We had no 
problems installing DiscoPlot on a recently acquired OSX Mavericks machine 
using the homebrew package manager.

The installed packages on this machine via::

    $ brew list 

Are available at this gist_.

pip install DiscoPlot::
    
    $ pip install --user DiscoPlot

We use the --user option of pip_ to put DiscoPlot in: /home/$USER/.local/bin/
You need to add this location to you ~/.bash_profile. 

Add DiscoPlot to your path::

    $ echo 'export PATH=$PATH:/home/$USER/.local/bin/' >> ~/.bash_profile

Finally install BLAST+::

    $ sudo brew install blast 


Windows
~~~~~~~
Download and install numpy and matplotlib.
To make this process easier you can download a distribution of python with matplotlib and numpy already installed
such as anaconda_.

pip install DiscoPlot::
    
    $ pip install DiscoPlot

Finally download and install BLAST_.


Testing DiscoPlot Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run::
    
    $ DiscoPlot -h 
    $ python -c 'import DiscoPlot; print DiscoPlot'





Upgrading DiscoPlot
~~~~~~~~~~~~~~~~~~~

You can upgrade like this::
    
    pip install --upgrade DiscoPlot

**Please regularly check back to make sure you're running the most recent 
DiscoPlot version.**


Example of figures produced by DiscoPlot
----------------------------------------

.. image:: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_3.png
        :target: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_3.png
        :alt: DiscoPlot figure
        :align: center

**DiscoPlot of a mock genome.** A mock genome was created by adding genomic 
rearrangements to the chromosome of E. coli str. UTI89.  Paired-end reads 
generated from the mock genome (query) with GemSim (ref) and mapped back to 
UTI89 (reference). The first ~500 Kbp were then visualised using DiscoPlot.

.. image:: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_4.png
    :target: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_4.png
    :alt: DiscoPlots of structural variants
    :align: center

**DiscoPlots of common structural variants.** Each box shows a common genomic 
rearrangement represented by a DiscoPlot. Rows A and B were created using 
100 bp long paired-end reads with an insert size of 300bp. Rows C and D were 
created using single-end reads with an average length of 1000bp. 
For each box the rearrangement in the sequenced genome is listed, followed by 
the scale of the gridlines in brackets.
A1,  C1: 300 bp deletion (400 bp).
A2, C2: 300 bp insertion (400 bp).
A3, C3: 300 bp inversion (400 bp).
A4, C4: 300 bp sequence translocated 50 Kbp upstream (10 Kbp). 
B1, D1: 3000 bp deletion (1000 bp). 
B2, D2: 3000 bp insertion (500 bp).
B3, D3: 3000 bp inversion (1000 bp). 
B4, D4: 3000 bp sequence translocated 50 Kbp upstream (10 Kbp). C1) 

.. image:: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_5.png
    :target: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_5.png
    :alt: DiscoPlot of E. coli genome
    :align: center

**The dynamic nature of the genome of Escherichia coli str. UTI89.** 
Discoplot of paired-end reads from a clonal culture of UTI89 mapped back 
to the published reference chromosome and plasmid. Coordinates from
0 to 5,065,741 represent the chromosome of E. coli UTI89, 
coordinates ≥ 5,066,000 represent the plasmid of E. coli UTI89


.. image:: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_6.png
    :target: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_6.png
    :alt: DiscoPlot of E. coli genome
    :align: center


Tutorials
---------

**Coming Soon**


Commands
--------


usage: DiscoPlot.py [-h] [-r READ_FILE] [-ref REFERENCE_FILE] [-bam BAM_FILE]
                    [-sam SAM_FILE] [-hm HEATMAP] [-B GEN_BLAST]
                    [-b BLAST_FILE] [-o OUTPUT_FILE] [-s SIZE] [-bin BIN_SIZE]
                    [-g GAP] [-sub SUBSECTION [SUBSECTION ...]]
                    [-wb WRITE_READS [WRITE_READS ...]] [-c MIN_HITS]
                    [-m MAX_HITS] [-dpi IMAGE_QUALITY] [-i MIN_IDENT]
                    [-l MIN_LENGTH] [-d UNMAPPED] [-a ALPHA] [-a2 ALPHA2]
                    [-mc M_COUNT] [-ms M_SIZE] [-log] [-sw] [-nl] [-ng] [-na]
                    [-split SPLIT_GRAPH [SPLIT_GRAPH ...]]
                    [-hl HIGHLIGHT [HIGHLIGHT ...]] [-mw MARKER_EDGE_WIDTH]

DiscoPlot.py - Visualising discordant reads.

USAGE: DiscoPlot.py -bam bamfile.bam -o output_file.bmp -size 5000
          Create a bmp file from a bamfile of paired-end reads with 5000 bins
       DiscoPlot.py -r reads.fa -B blast_prefix -r reference -o output_file.png -bin 10000
          Create a png file using reads.fa aligned to the reference, automatically generate blast file. Use a bin size of 10,000bp.

In paired read mode DiscoPlot must be provided with a BAM or SAM file.
In Single read mode DiscoPlit must be provided with a alignment file (BLAST tab delimited format) or reads and a reference (in FASTA format).

-bin (size of bins in bp) or -s (size of bins) must be specified.

``
optional arguments:
  -h, --help            show this help message and exit
  -r READ_FILE, --read_file READ_FILE
                        read file - provide DiscoPlot with a read file to
                        BLAST (long read mode).
  -ref REFERENCE_FILE, --reference_file REFERENCE_FILE
                        Reference file - Reference for generating long reads
                        alignments.
  -bam BAM_FILE, --bam_file BAM_FILE
                        bam file - paired read mode. (Requires pysam).
  -sam SAM_FILE, --sam_file SAM_FILE
                        sam file - paired read mode. (pysam not required)
  -hm HEATMAP, --heatmap HEATMAP
                        Heatmap file - provide DiscoPlot with custom generated
                        heatmap.
  -B GEN_BLAST, --gen_blast GEN_BLAST
                        Generate blast files, use argument as prefix for
                        output.
  -b BLAST_FILE, --blast_file BLAST_FILE
                        Provide DiscoPlot with alignment file (long read mode)
                        (BLAST tab delimited file - output format 6)
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        output file [gif/bmp/png]
  -s SIZE, --size SIZE  Number of bins.
  -bin BIN_SIZE, --bin_size BIN_SIZE
                        Bin size (in bp)
  -g GAP, --gap GAP     Gap size - gap size between entries in reference.
  -sub SUBSECTION [SUBSECTION ...], --subsection SUBSECTION [SUBSECTION ...]
                        Only display subection of genome [ref]/[min_cutoff
                        max_cutoff]/[ref min_cutoff max_cutoff]
  -wb WRITE_READS [WRITE_READS ...], --write_reads WRITE_READS [WRITE_READS ...]
                        Write reads in rectangle to bam/sam [x1 y1 x2 y2
                        out.bam]
  -c MIN_HITS, --min_hits MIN_HITS
                        Only show bins with more than this number of hits.
  -m MAX_HITS, --max_hits MAX_HITS
                        Only show bins with less hits than this.
  -dpi IMAGE_QUALITY, --image_quality IMAGE_QUALITY
                        Image quality (in DPI)
  -i MIN_IDENT, --min_ident MIN_IDENT
                        Min. idenity of hits to draw (long read mode).
  -l MIN_LENGTH, --min_length MIN_LENGTH
                        Min. length of hits to draw (long read mode).
  -d UNMAPPED, --unmapped UNMAPPED
                        Unmapped bases on edge for RMaD to consider read
                        partially unmapped.
  -a ALPHA, --alpha ALPHA
                        Transparency of mapped read markers
  -a2 ALPHA2, --alpha2 ALPHA2
                        Transparency of unmapped read markers
  -mc M_COUNT, --m_count M_COUNT
                        The count of a bin to be used as the median value for
                        calculating the size of the dot [auto]
  -ms M_SIZE, --m_size M_SIZE
                        Set the width (in bins) of a marker with a median
                        count.
  -log, --log           Log10 bin counts. (For data with highly variable
                        coverage).
  -sw, --switch         Draw most common (inverted/direct) hits first.
  -nl, --no_legend      Don't create legend.
  -ng, --no_gridlines   Don't draw gridlines.
  -na, --no_label       No axis labels.
  -split SPLIT_GRAPH [SPLIT_GRAPH ...], --split_graph SPLIT_GRAPH [SPLIT_GRAPH ...]
                        Show multiple subsections of graph [start1 stop1
                        start2 stop2 etc.] or [ref1 start1 stop1 ref2 start2
                        stop2 etc.]
  -hl HIGHLIGHT [HIGHLIGHT ...], --highlight HIGHLIGHT [HIGHLIGHT ...]
                        Highlight subsections of graph [alpha start1 stop1
                        start2 stop2 etc.] or [alphref1 start1 stop1 ref2
                        start2 stop2 etc.]
  -mw MARKER_EDGE_WIDTH, --marker_edge_width MARKER_EDGE_WIDTH
                        Marker width (default is roughly 20x bin size)

Thanks for using DiscoPlot.py
``




.. _pip: http://www.pip-installer.org/en/latest/
.. _libraries: https://github.com/BeatsonLab-MicrobialGenomics/DiscoPlot/blob/master/requirements.txt
.. _gist: https://gist.github.com/mscook/ef7499fc9d2138f17c7f
.. _pip installation instructions: http://pip.readthedocs.org/en/latest/installing.html
.. _anaconda: http://continuum.io/downloads
.. _BLAST: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
