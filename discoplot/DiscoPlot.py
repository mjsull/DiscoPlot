#!/usr/bin/env python

# DiscoPlot: identify genomic rearrangements, misassemblies and sequencing
# artefacts in NGS data
# Copyright (C) 2013-2015 Mitchell Sullivan
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Mitchell Sullivan
# mjsull@gmail.com
# School of Chemistry & Molecular Biosciences
# The University of Queensland
# Brisbane, QLD 4072.
# Australia

__title__ = 'DiscoPlot'
__version__ = '1.0.2'
__description__ = ("DiscoPlot: identify genomic rearrangements, misassemblies "
                   "and sequencing artefacts in NGS data")
__author__ = 'Mitchell Sullivan'
__license__ = 'GPLv3'
__author_email__ = "mjsull@gmail.com"
__url__ = 'https://github.com/BeatsonLab-MicrobialGenomics/DiscoPlot'

import argparse
import numpy
import sys
import subprocess


def read_sbam(args):
    import pysam
    if not args.bam_file is None:
        sam = pysam.Samfile(args.bam_file, 'rb')
    elif not args.sam_file:
        sam = pysam.Samfile(args.sam_file)
    global refpos
    global cuta
    global cutb
    cuta = 0
    cutb = float('inf')
    refpos = {}
    if not args.subsection is None:
        if len(args.subsection) == 1:
            refpos[args.subsection[0]] = 0
            totallength = None
            for i in range(0, len(sam.references)):
                if sam.references[i] == args.subsection[0]:
                    totallength = sam.lengths[i]
            if totallength is None:
                sys.stderr.write('Selected reference not found.')
                sys.exit()
        elif len(args.subsection) == 2:
            refpos[sam.references[0]] = 0
            cuta = int(args.subsection[0])
            cutb = int(args.subsection[1])
            totallength = cutb - cuta
        elif len(args.subsection) == 3:
            refpos[args.subsection[0]] = 0
            cuta = int(args.subsection[1])
            cutb = int(args.subsection[2])
            totallength = cutb - cuta
        else:
            sys.stderr.write('Too many arguments given for subsection')
            sys.exit()
        if args.bin_size is None:
            args.bin_size = totallength / args.size + 1
        else:
            args.size = totallength / args.bin_size + 1
    else:
        references = sam.references
        reflengths = sam.lengths
        currpos = 0
        if args.bin_size is None:
            args.bin_size = sum(reflengths) / (args.size - (len(reflengths) -1) * (args.gap + 1)) + 1
        else:
            args.size = sum(map(lambda x: x/args.bin_size, reflengths)) + (len(reflengths) -1) * args.gap + 1
        for i in range(len(references)):
            refpos[references[i]] = currpos
            currpos += reflengths[i] / args.bin_size + args.gap
    global invgrid, dirgrid, unmapped_for, unmapped_rev
    unmapped_rev = {}
    unmapped_for = {}
    invgrid = {}
    dirgrid = {}
    for read in sam.fetch():
        ref = sam.getrname(read.tid)
        if ref in refpos:
            if read.is_read1:
                if cuta <= read.pos <= cutb:
                    pos1 = (read.pos - cuta) / args.bin_size + refpos[ref]
                    if read.mate_is_unmapped:
                        if read.is_reverse:
                            if pos1 in unmapped_rev:
                                unmapped_rev[pos1] += 1
                            else:
                                unmapped_rev[pos1] = 1
                        else:
                            if pos1 in unmapped_for:
                                unmapped_for[pos1] += 1
                            else:
                                unmapped_for[pos1] = 1
                    else:
                        mref = sam.getrname(read.rnext)
                        if mref in refpos:
                            if cuta <= read.pnext <= cutb:
                                pos2 = (read.pnext - cuta) / args.bin_size + refpos[mref]
                                if read.is_reverse:
                                    if read.mate_is_reverse:
                                        if pos1 < pos2:
                                            if pos2 in dirgrid and pos1 in dirgrid[pos2]:
                                                dirgrid[pos2][pos1] += 1
                                            elif pos2 in dirgrid:
                                                dirgrid[pos2][pos1] = 1
                                            else:
                                                dirgrid[pos2] = {pos1:1}
                                        else:
                                            if pos1 in dirgrid and pos2 in dirgrid[pos1]:
                                                dirgrid[pos1][pos2] += 1
                                            elif pos1 in dirgrid:
                                                dirgrid[pos1][pos2] = 1
                                            else:
                                                dirgrid[pos1] = {pos2:1}
                                    else:
                                        if pos2 in invgrid and pos1 in invgrid[pos2]:
                                            invgrid[pos2][pos1] += 1
                                        elif pos2 in invgrid:
                                            invgrid[pos2][pos1] = 1
                                        else:
                                            invgrid[pos2] = {pos1:1}
                                else:
                                    if read.mate_is_reverse:
                                        if pos1 in invgrid and pos2 in invgrid[pos1]:
                                            invgrid[pos1][pos2] += 1
                                        elif pos1 in invgrid:
                                            invgrid[pos1][pos2] = 1
                                        else:
                                            invgrid[pos1] = {pos2:1}
                                    else:
                                        if pos1 < pos2:
                                            if pos1 in dirgrid and pos2 in dirgrid[pos1]:
                                                dirgrid[pos1][pos2] += 1
                                            elif pos1 in dirgrid:
                                                dirgrid[pos1][pos2] = 1
                                            else:
                                                dirgrid[pos1] = {pos2:1}
                                        else:
                                            if pos2 in dirgrid and pos1 in dirgrid[pos2]:
                                                dirgrid[pos2][pos1] += 1
                                            elif pos2 in dirgrid:
                                                dirgrid[pos2][pos1] = 1
                                            else:
                                                dirgrid[pos2] = {pos1:1}
            else:
                if read.mate_is_unmapped:
                    ref = sam.getrname(read.tid)
                    if ref in refpos:
                        if cuta <= read.pos <= cutb:
                            pos = (read.pos - cuta) / args.bin_size + refpos[ref]
                            if read.is_reverse:
                                if pos in unmapped_rev:
                                    unmapped_rev[pos] += 1
                                else:
                                    unmapped_rev[pos] = 1
                            else:
                                if pos in unmapped_for:
                                    unmapped_for[pos] += 1
                                else:
                                    unmapped_for[pos] = 1


def read_sing(args):
    readlen = None
    if not args.read_file is None:
        reads = open(args.read_file)
        first = True
        getfq = 0
        readlen = {}
        for line in reads:
            if first:
                first = False
                if line.startswith('@'):
                    getfq = 2
                name = line.rstrip()[1:]
                seq = ''
            elif line.startswith('>'):
                readlen[name] = len(seq)
                name = line.rstrip()[1:]
                seq = ''
            elif getfq == 0:
                seq += line.rstrip()
            elif getfq == 1:
                readlen[name] = len(seq)
                name = line.rstrip()
                seq = ''
            elif getfq == 2:
                seq += line.rstrip()
                getfq = 3
            elif getfq == 3:
                getfq = 4
            elif getfq == 4:
                getfq = 1
        readlen[name] = len(seq)
    if not args.reference_file is None:
        ref = open(args.reference_file)
        first = True
        references = []
        reflengths = []
        for line in ref:
            if line.startswith('>'):
                if first:
                    first = False
                else:
                    references.append(name)
                    reflengths.append(len(seq))
                name = line.rstrip()[1:]
                seq = ''
            else:
                seq += line
        references.append(name)
        reflengths.append(len(seq))
    else:
        blast = open(args.blast_file)
        refdict = {}
        for line in blast:
            if line.split()[1] in refdict:
                if max([int(line.split()[8]), int(line.split()[9])]) > refdict[line.split()[1]]:
                    refdict[line.split()[1]] = max([int(line.split()[8]), int(line.split()[9])])
            else:
                refdict[line.split()[1]] = max([int(line.split()[8]), int(line.split()[9])])
        blast.close()
        references = []
        reflengths = []
        for i in refdict:
            references.append(i)
            reflengths.append(refdict[i])
    cuta = 0
    cutb = float('inf')
    refpos = {}
    if not args.subsection is None:
        if len(args.subsection) == 1:
            refpos[args.subsection[0]] = 0
            totallength = None
            for i in range(0, len(references)):
                if references[i] == args.subsection[0]:
                    totallength = reflengths[i]
            if totallength is None:
                sys.stderr.write('Selected reference not found.')
                sys.exit()
        elif len(args.subsection) == 2:
            refpos[references[0]] = 0
            cuta = int(args.subsection[0])
            cutb = int(args.subsection[1])
            totallength = cutb - cuta
        elif len(args.subsection) == 3:
            refpos[args.subsection[0]] = 0
            cuta = int(args.subsection[0])
            cutb = int(args.subsection[1])
            totallength = cutb - cuta
        else:
            sys.stderr.write('Too many arguments given for subsection')
            sys.exit()
        if args.bin_size is None:
            args.bin_size = totallength / args.size
        else:
            args.size = totallength / args.bin_size
    else:
        currpos = 0
        if args.bin_size is None:
            args.bin_size = sum(reflengths) / (args.size - (len(reflengths) -1) * (args.gap + 1))
        else:
            args.size = sum(map(lambda x: x/args.bin_size, reflengths)) + (len(reflengths) -1) * args.gap
        for i in range(len(references)):
            refpos[references[i]] = currpos
            currpos += reflengths[i] / args.bin_size + args.gap
    global invgrid, dirgrid, unmapped_for, unmapped_rev
    unmapped_rev = {}
    unmapped_for = {}
    invgrid = {}
    dirgrid = {}
    blast = open(args.blast_file)
    lastquery = ''
    hits = []
    for line in blast:
        query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
        qstart, qstop, rstart, rstop, length, mm, indel = map(int, [qstart, qstop, rstart, rstop, length, mm, indel])
        if query != lastquery and lastquery != '':
            hits.sort(reverse=True)
            newhits = [hits[0]]
            qtaken = set()
            for i in range(hits[2], hits[3] + 1):
                qtaken.add(i)
            for i in hits[1:]:
                if i[:-3] == newhits[-1][:-3]:
                    newhits.append(i)
                else:
                    getit = False
                    for j in range(hits[2], hits[3] + 1):
                        if not j in qtaken:
                            getit = True
                            qtaken.add(j)
                    if getit:
                        newhits.append(i)
            anchor = None
            revseq = None
            for i in newhits:
                bitscore, length, qstart, qstop, rstart, rstop, subject = i
                if anchor is None:
                    if rstart < rstop:
                        anchor = rstart
                        revseq = False
                    else:
                        anchor = rstop
                        revseq = True
                    if min(qtaken) >= args.unmapped:
                        if revseq:
                            if anchor in unmapped_for:
                                unmapped_for[anchor] += 1
                            else:
                                unmapped_for[anchor] = 1
                        else:
                            if anchor in unmapped_rev:
                                unmapped_rev[anchor] += 1
                            else:
                                unmapped_rev[anchor] = 1
                    if max(qtaken) <= readlen[lastquery] - args.unmapped:
                        if revseq:
                            if anchor in unmapped_rev:
                                unmapped_rev[anchor] += 1
                            else:
                                unmapped_rev[anchor] = 1
                        else:
                            if anchor in unmapped_for:
                                unmapped_for[anchor] += 1
                            else:
                                unmapped_for[anchor] = 1
                lastxpos = None
                lastypos = None
                oldstart, oldstop = qstart, qstop
                if revseq:
                    rstart, rstop = rstop, rstart
                    qstart = readlen[lastquery] - qstop
                    qstop = readlen[lastquery] - oldstart
                for j in range(qstart, qstop):
                    xpos = refpos[subject] + (anchor + j - cuta) / args.bin_size
                    ypos = refpos[subject] + (rstart + int(((j - qstart) * 1.0 / (qstop - qstart)) * (rstop - rstart))) / args.bin_size
                    if xpos != lastxpos or ypos != lastypos:
                        if rstart < rstop:
                            if xpos in dirgrid:
                                if ypos in dirgrid[xpos]:
                                    dirgrid[xpos][ypos] += 1
                                else:
                                    dirgrid[xpos][ypos] = 1
                            else:
                                dirgrid[xpos] = {ypos:1}
                        else:
                            if xpos in invgrid:
                                if ypos in invgrid[xpos]:
                                    invgrid[xpos][ypos] += 1
                                else:
                                    invgrid[xpos][ypos] = 1
                            else:
                                invgrid[xpos] = {ypos:1}
                    lastxpos, lastypos = xpos, ypos

        if ident >= args.min_ident and length >= args.min_length and subject in refpos and ((cuta <= rstart <= cutb) or (cuta <= rstop <= cutb)):
            hits.append((float(bitscore), length, qstart, qstop, rstart, rstop, subject))
        lastquery = query


def generate_blast(args):
    subprocess.Popen('makeblastdb -dbtype nucl -out ' + args.gen_blast + '.db -in ' +
                     args.reference_file, shell=True, stdout=subprocess.PIPE).wait()
    subprocess.Popen('blastn -db ' + args.gen_blast + '.db -outfmt 6 -query ' +
                      args.read_file + ' -out ' + args.gen_blast + '.out', shell=True).wait()
    args.blast_file = args.gen_blast + '.out'


def draw_dotplot(args):
    global refpos
    global cuta
    global cutb
    vals1, vals2 = [], []
    for i in invgrid:
        for j in invgrid[i]:
            vals1.append(invgrid[i][j])
            vals2.append(invgrid[i][j])
    for i in dirgrid:
        for j in dirgrid[i]:
            vals1.append(dirgrid[i][j])
            vals2.append(dirgrid[i][j])
    vals2 = numpy.array(vals2)
    for i in unmapped_rev:
        vals1.append(unmapped_rev[i])
    for i in unmapped_for:
        vals1.append(unmapped_for[i])
    vals1 = numpy.array(vals1)
    med = numpy.median(vals2)
    numvals = numpy.size(vals1)
    sizemod = 2000.0 / args.size / med
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, aspect='equal')
    x = numpy.zeros(numvals, dtype='u4')
    y = numpy.zeros(numvals, dtype='u4')
    sizes = numpy.zeros(numvals, dtype='f4')
    colours = numpy.array(['x' for i in range(numvals)])
    count = 0
    for i in dirgrid:
        for j in dirgrid[i]:
            if args.max_hits >= dirgrid[i][j] >= args.min_hits:
                x[count] = i * args.bin_size + cuta
                y[count] = j * args.bin_size + cuta
                sizes[count] = dirgrid[i][j] * sizemod
                colours[count] = 'r'
                count += 1
    for i in invgrid:
        for j in invgrid[i]:
            if args.max_hits >= invgrid[i][j] >= args.min_hits:
                x[count] = i * args.bin_size + cuta
                y[count] = j * args.bin_size + cuta
                sizes[count] = invgrid[i][j] * sizemod
                colours[count] = 'b'
                count += 1
    for i in unmapped_for:
        if args.max_hits >= unmapped_for[i] >= args.min_hits:
            x[count] = cuta
            y[count] = i * args.bin_size + cuta
            sizes[count] = unmapped_for[i] * sizemod
            colours[count] = 'g'
            count += 1
    for i in unmapped_rev:
        if args.max_hits >= unmapped_rev[i] >= args.min_hits:
            x[count] = i * args.bin_size + cuta
            y[count] = cuta
            sizes[count] = unmapped_rev[i] * sizemod
            colours[count] = 'g'
            count += 1
    count1, count2, count3 = 0, 0, 0
    for i in colours:
        if i == 'b':
            count1 += 1
        elif i == 'r':
            count2 += 1
        elif i == 'g':
            count3 += 1
    ax.scatter(x, y, s=sizes, c=colours, edgecolor='none', alpha=0.3)
    sizes = []
    names = []
    for i in [10, 25, 50, 75, 90]:
        sizes.append(numpy.percentile(vals2, i))
        names.append(str(i) + '% Normal ' + str(sizes[-1]))
    names.append('50% Inverted ' + str(sizes[2]))
    a = plt.scatter(-100, -100, s=sizes[2] * sizemod, c='b', edgecolor='none', alpha=0.3)
    b = plt.scatter(-100, -100, s=sizes[0] * sizemod, c='r', edgecolor='none', alpha=0.3)
    c = plt.scatter(-100, -100, s=sizes[1] * sizemod, c='r', edgecolor='none', alpha=0.3)
    d = plt.scatter(-100, -100, s=sizes[2] * sizemod, c='r', edgecolor='none', alpha=0.3)
    e = plt.scatter(-100, -100, s=sizes[3] * sizemod, c='r', edgecolor='none', alpha=0.3)
    f = plt.scatter(-100, -100, s=sizes[4] * sizemod, c='r', edgecolor='none', alpha=0.3)
    leg = ax.legend([b, c, d, e, f, a], names, loc=4)
    leg.draggable(state=True)
    for i in refpos:
        if not refpos[i] == 0:
            ax.axhspan(refpos[i] * args.bin_size, refpos[i] * args.bin_size - args.gap * args.bin_size, facecolor='g', alpha=0.3)
            ax.axvspan(refpos[i] * args.bin_size, refpos[i] * args.bin_size - args.gap * args.bin_size, facecolor='g', alpha=0.3)
    if cutb == float('inf'):
        cutb = args.size * args.bin_size + cuta
    plt.xlim([cuta - args.bin_size * 10, cutb])
    plt.ylim([cuta - args.bin_size * 10, cutb])
    plt.grid(True)
    if not args.output_file is None:
        plt.savefig(args.output_file, dpi=args.image_quality)
    else:
        plt.show()



parser = argparse.ArgumentParser(prog='DiscoPlot', formatter_class=argparse.RawDescriptionHelpFormatter, description='''
DiscoPlot - read mapping visualisation in the large

USAGE: DiscoPlot -bam bamfile.bam -o output_file.bmp -size 5000
          Create a bmp file from a bamfile of paired-end reads with a width and height of 5000px
       DiscoPlot -r reads.fa -B blast_prefix -r reference -o output_file.png -bin bin_size
          Create a png file from reads.fa, generate blast file. Image size will be reference length / bin_size
''', epilog="Thanks for using DiscoPlot")
parser.add_argument('-r', '--read_file', action='store', default=None, help='read file')
parser.add_argument('-ref', '--reference_file', action='store', default=None, help='reference file')
parser.add_argument('-bam', '--bam_file', action='store', default=None, help='bam file')
parser.add_argument('-sam', '--sam_file', action='store', default=None, help='sam file')
parser.add_argument('-B', '--gen_blast', action='store', default=None, help='Generate blast files, use argument as prefix for output.')
parser.add_argument('-b', '--blast_file', action='store', default=None, help='Blast file (output format 6)')
parser.add_argument('-o', '--output_file', action='store', default=None, help='output file [gif/bmp/png]')
parser.add_argument('-s', '--size', action='store', type=int, default=None, help='Number of bins')
parser.add_argument('-bin', '--bin_size', action='store', type=int, default=None, help='Bin size (in bp)')
parser.add_argument('-g', '--gap', action='store', type=int, default=5, help='Gap size')
parser.add_argument('-sub', '--subsection', nargs='+', action='store', default=None, help='Only display subection of genome [ref]/[min_cutoff max_cutoff]/[ref min_cutoff max_cutoff]')
parser.add_argument('-c', '--min_hits', action='store', type=int, default=1, help='Min hits to be shown')
parser.add_argument('-m', '--max_hits', action='store', type=float, default=float('inf'), help='Bins with more hits than this will be skipped.')
parser.add_argument('-dpi', '--image_quality', action='store', type=int, default=1600, help='Image quality (in DPI)')


args = parser.parse_args()
if args.size is None and args.bin_size is None:
    sys.stderr.write('Please give a image size or bin size.')
    sys.exit()

if not args.gen_blast is None:
    if args.reference_file is None:
        sys.stderr.write('Please provide a reference file')
        sys.exit()
    if args.read_file is None:
        sys.stderr.write('Please provide a read file (FASTA)')
        sys.exit()
    generate_blast(args)

if not args.output_file is None:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt

if not args.size is None and not args.bin_size is None:
    sys.stderr.write('Only provide bin size or image size, not both.')
    sys.exit()
if not args.sam_file is None or not args.bam_file is None:
    read_sbam(args)
elif args.blast_file is None:
    sys.stderr.write('Please either generate or provide a BLAST comparison')
    sys.exit()
else:
    read_sing(args)
draw_dotplot(args)
