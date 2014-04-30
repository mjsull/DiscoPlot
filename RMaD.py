import argparse
import numpy
import pysam
import sys
import subprocess
import random
import matplotlib.pyplot as plt
import scipy.sparse as sparse

def read_sbam(args):
    if not args.bam_file is None:
        sam = pysam.Samfile(args.bam_file, 'rb')
    else:
        sam = pysam.Samfile(args.sam_file)
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
        totallength = totallength / args.bin_size
        print totallength
    else:
        references = sam.references
        reflengths = sam.lengths
        currpos = 0
        if args.bin_size is None:
            args.bin_size = sum(reflengths) / (args.size - (len(reflengths) -1) * (args.gap + 1))
        else:
            args.size = sum(map(lambda x: x/args.bin_size, reflengths)) + (len(reflengths) -1) * args.gap
        for i in range(len(references)):
            refpos[references[i]] = currpos
            currpos += reflengths[i] / args.bin_size + args.gap
        totallength = currpos - args.gap + 1
    global invgrid, dirgrid, unmapped_for, unmapped_rev
    unmapped_rev = numpy.zeros(totallength, dtype='u2')
    unmapped_for = numpy.zeros(totallength, dtype='u2')
    invgrid = numpy.zeros((totallength, totallength), dtype='u2')
    dirgrid = numpy.zeros((totallength, totallength), dtype='u2')
    for read in sam.fetch():
        ref = sam.getrname(read.tid)
        if ref in refpos:
            if read.is_read1:
                if cuta <= read.pos <= cutb:
                    pos1 = (read.pos - cuta) / args.bin_size + refpos[ref]
                    if read.mate_is_unmapped:
                        if read.is_reverse:
                            unmapped_rev[pos1] += 1
                        else:
                            unmapped_for[pos1] += 1
                    else:
                        mref = sam.getrname(read.rnext)
                        if mref in refpos:
                            if cuta <= read.pnext <= cutb:
                                pos2 = (read.pnext - cuta) / args.bin_size + refpos[mref]
                                if read.is_reverse:
                                    if read.mate_is_reverse:
                                        if pos1 < pos2:
                                            dirgrid[pos2][pos1] += 1
                                        else:
                                            dirgrid[pos1][pos2] += 1
                                    else:
                                        invgrid[pos2][pos1] += 1
                                else:
                                    if read.mate_is_reverse:
                                        invgrid[pos1][pos2] += 1
                                    else:
                                        if pos1 < pos2:
                                            dirgrid[pos1][pos2] += 1
                                        else:
                                            dirgrid[pos2][pos1] += 1
            else:
                if read.mate_is_unmapped:
                    ref = sam.getrname(read.tid)
                    if ref in refpos:
                        if cuta <= read.pos <= cutb:
                            pos = (read.pos - cuta) / args.bin_size + refpos[ref]
                            if read.is_reverse:
                                unmapped_rev[pos] += 1
                            else:
                                unmapped_for[pos] += 1


def read_sing(args):
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
    global refpos
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
        totallength = totallength / args.bin_size
        print totallength
    else:
        currpos = 0
        if args.bin_size is None:
            args.bin_size = sum(reflengths) / (args.size - (len(reflengths) -1) * (args.gap + 1))
        else:
            args.size = sum(map(lambda x: x/args.bin_size, reflengths)) + (len(reflengths) -1) * args.gap
        for i in range(len(references)):
            refpos[references[i]] = currpos
            currpos += reflengths[i] / args.bin_size + args.gap
        totallength = currpos - args.gap + 1
    global invgrid, dirgrid, unmapped_for, unmapped_rev
    unmapped_rev = numpy.zeros(totallength, dtype='u2')
    unmapped_for = numpy.zeros(totallength, dtype='u2')
    # invgrid = sparse.lil_matrix((totallength, totallength), dtype='u2')
    # dirgrid = sparse.lil_matrix((totallength, totallength), dtype='u2')
    invgrid = numpy.zeros((totallength, totallength), dtype='u2')
    dirgrid = numpy.zeros((totallength, totallength), dtype='u2')
    blast = open(args.blast_file)
    lastquery = ''
    hits = []
    for line in blast:
        query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
        qstart, qstop, rstart, rstop, length, mm, indel = map(int, [qstart, qstop, rstart, rstop, length, mm, indel])
        if query != lastquery and lastquery != '':
            hits.sort(reverse=True)
            finalhits = []
            qhits = set()
            lasthit = [None]
            randomhits = []
            for i in hits:
                if not lasthit[0] is None and i[:-3] == lasthit[0][:-3]:
                    lasthit.append(i)
                else:
                    if lasthit != [None]:
                        randomhits.append(random.choice(lasthit))
                    lasthit = [i]
            if lasthit != [None]:
                randomhits.append(random.choice(lasthit))
            for i in randomhits:
                getit = False
                for j in range(i[2], i[3] + 1):
                    if not j in qhits:
                        qhits.add(j)
                        getit = True
                if getit:
                    finalhits.append((i[2], i[3], i[4], i[5], i[6]))
            finalhits.sort()
            anchor = min([finalhits[0][2], finalhits[0][3]]) - finalhits[0][0] + 1
            anchorref = finalhits[0][4]
            anchororient = finalhits[0][2] < finalhits[0][3]
            if finalhits[0][0] >= 30:
                pos = (finalhits[0][2] - cuta) / args.bin_size + refpos[finalhits[0][4]]
                unmapped_rev[pos] += 1
            for i in finalhits:
                lastpos = None
                if i[4] in refpos:
                    if anchororient:
                        if i[2] < i[3]:
                            for j in range(i[0], i[1] + 1):
                                if (cuta <= i[2] + j <= cutb) and (cuta <= anchor + j <= cutb):
                                    xpos = (anchor + j - cuta) / args.bin_size + refpos[anchorref]
                                    ypos = (i[2] + j - cuta) / args.bin_size + refpos[i[4]]
                                    dirgrid[xpos,ypos] += 1
                                    lastpos = (xpos, ypos)
                        else:
                            for j in range(i[0], i[1] + 1):
                                if (cuta <= i[3] - j <= cutb) and (cuta <= anchor + j <= cutb):
                                    xpos = (anchor + j - cuta) / args.bin_size + refpos[anchorref]
                                    ypos = (i[3] - j - cuta) / args.bin_size + refpos[i[4]]
                                    invgrid[xpos,ypos] += 1
                                    lastpos = (xpos, ypos)
                    else:
                        if i[2] < i[3]:
                            for j in range(i[0], i[1] + 1):
                                if (cuta <= i[2] + j <= cutb) and (cuta <= anchor - j <= cutb):
                                    xpos = (anchor - j - cuta) / args.bin_size + refpos[anchorref]
                                    ypos = (i[2] + j - cuta) / args.bin_size + refpos[i[4]]
                                    invgrid[xpos,ypos] += 1
                                    lastpos = (xpos, ypos)
                        else:
                            for j in range(i[0], i[1] + 1):
                                if (cuta <= i[3] - j <= cutb) and (cuta <= anchor - j <= cutb):
                                    xpos = (anchor - j - cuta) / args.bin_size + refpos[anchorref]
                                    ypos = (i[3] - j - cuta) / args.bin_size + refpos[i[4]]
                                    dirgrid[xpos,ypos] += 1
                                    lastpos = (xpos, ypos)

            hits = []
        hits.append((float(bitscore), length, qstart, qstop, rstart, rstop, subject))
        lastquery = query


def generate_blast(args):
    subprocess.Popen('makeblastdb -dbtype nucl -out ' + args.gen_blast + '.db -in ' +
                     args.reference_file, shell=True, stdout=subprocess.PIPE).wait()
    subprocess.Popen('blastn -db ' + args.gen_blast + '.db -outfmt 6 -query ' +
                      args.read_file + ' -out ' + args.gen_blast + '.out', shell=True).wait()
    args.blast_file = args.gen_blast + '.out'


def draw_dotplot(args):
    vals = numpy.concatenate((numpy.extract(invgrid >= args.min_hits, invgrid), numpy.extract(dirgrid >= args.min_hits, dirgrid)))#,
#                             numpy.extract(unmapped_for >= args.min_hits, unmapped_for), numpy.extract(unmapped_rev >= args.min_hits, unmapped_rev)))
    med = numpy.median(vals)
    numvals = numpy.size(vals)
    sizemod = 100.0 / med
    gridsize = numpy.shape(dirgrid)[0]
    print 'ding'
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    x = numpy.zeros(numvals, dtype='u4')
    y = numpy.zeros(numvals, dtype='u4')
    sizes = numpy.zeros(numvals, dtype='u4')
    colours = numpy.array(['x' for i in range(numvals)])
    count = 0
    for i in range(gridsize):
        for j in range(gridsize):
            if dirgrid[i,j] >= args.min_hits:
                x[count] = i * args.bin_size
                y[count] = j * args.bin_size
                sizes[count] = dirgrid[i,j] * sizemod
                colours[count] = 'r'
                # h = int((dirgrid[i][j] * sizemod) ** 0.5)
                # ax.add_artist(Rectangle(xy=(i * args.bin_size, j * args.bin_size),
                #   color='red', width=h, height=h, alpha=0.3))
                count += 1
            if invgrid[i,j] >= args.min_hits:
                x[count] = i * args.bin_size
                y[count] = j * args.bin_size
                sizes[count] = dirgrid[i,j] * sizemod
                colours[count] = 'b'
                # h = int((dirgrid[i][j] * sizemod) ** 0.5)
                # ax.add_artist(Rectangle(xy=(i * args.bin_size, j * args.bin_size),
                #   color='blue', width=h, height=h, alpha=0.3))
                count += 1
    ax.scatter(x, y, s=sizes, c=colours, alpha=0.3)
    ax.legend(loc="bottom right")
    for i in refpos:
        if not refpos[i] == 0:
            ax.axhspan(refpos[i] * args.bin_size, refpos[i] * args.bin_size - args.gap * args.bin_size, facecolor='g', alpha=0.3)
            ax.axvspan(refpos[i] * args.bin_size, refpos[i] * args.bin_size - args.gap * args.bin_size, facecolor='g', alpha=0.3)
    print 'dong'
    plt.xlim([0, gridsize * args.bin_size])
    plt.ylim([0, gridsize * args.bin_size])
    plt.show()



def draw_grid(args):
    vals = numpy.concatenate((numpy.extract(invgrid >= args.min_hits, invgrid), numpy.extract(dirgrid >= args.min_hits, dirgrid),
                             numpy.extract(unmapped_for >= args.min_hits, unmapped_for), numpy.extract(unmapped_rev >= args.min_hits, unmapped_rev)))
    gridsize = numpy.shape(invgrid)[0]
    im = Image.new("RGB", (gridsize, gridsize + 50), "white")
    pix = im.load()
    scale = set()
    if not args.scale is None:
        for i in range(0, gridsize):
            if i * args.bin_size / int(args.scale) != (i + 1) * args.bin_size / int(args.scale):
                scale.add(i)
    percentiles = []
    for i in range(5, 96, 5):
        percentiles.append(numpy.percentile(vals, i))
    refposit = set()
    for i in refpos:
        if refpos[i] != 0:
            for j in range(refpos[i] - args.gap, refpos[i]):
                refposit.add(j)
    for i in range(gridsize):
        for j in range(gridsize):
            if dirgrid[i][j] >= args.min_hits or invgrid[i][j] >= args.min_hits:
                numhits = dirgrid[i][j] + invgrid[i][j]
                if numhits < percentiles[0]:
                    xrat = 0.05
                elif numhits < percentiles[1]:
                    xrat = 0.1
                elif numhits < percentiles[2]:
                    xrat = 0.15
                elif numhits < percentiles[3]:
                    xrat = 0.2
                elif numhits < percentiles[4]:
                    xrat = 0.25
                elif numhits < percentiles[5]:
                    xrat = 0.3
                elif numhits < percentiles[6]:
                    xrat = 0.35
                elif numhits < percentiles[7]:
                    xrat = 0.4
                elif numhits < percentiles[8]:
                    xrat = 0.45
                elif numhits < percentiles[9]:
                    xrat = 0.5
                elif numhits < percentiles[10]:
                    xrat = 0.55
                elif numhits < percentiles[11]:
                    xrat = 0.6
                elif numhits < percentiles[12]:
                    xrat = 0.65
                elif numhits < percentiles[13]:
                    xrat = 0.7
                elif numhits < percentiles[14]:
                    xrat = 0.75
                elif numhits < percentiles[15]:
                    xrat = 0.8
                elif numhits < percentiles[16]:
                    xrat = 0.85
                elif numhits < percentiles[17]:
                    xrat = 0.9
                elif numhits < percentiles[18]:
                    xrat = 0.95
                else:
                    xrat = 1
                yrat = invgrid[i][j] * 1.0 / (dirgrid[i][j] + invgrid[i][j])
                if xrat > 1.0:
                    xrat = 1.0
                if yrat > 1.0:
                    yrat = 1.0
                tr = 255
                tb = 0
                tg = 255 - xrat * 255
                br = 0
                bb = xrat * 255
                bg = 255 - xrat * 255
                r = int(tr * yrat + br * (1-yrat))
                g = int(tg * yrat + bg * (1-yrat))
                b = int(tb * yrat + bb * (1-yrat))
                pix[i, gridsize - j - 1] = (r,g,b)
            else:
                if i in scale or j in scale:
                    pix[i, gridsize - j - 1] = (0, 0, 0)
                elif i in refposit or j in refposit:
                    pix[i, gridsize - j - 1] = (100, 100, 100)

    for i in range(20):
        for j in range(20):
            xrat = i * 5.0 / 100
            yrat = j * 5.0 / 100
            tr = 255
            tb = 0
            tg = 255 - xrat * 255
            br = 0
            bb = xrat * 255
            bg = 255 - xrat * 255
            r = int(tr * yrat + br * (1- yrat))
            g = int(tg * yrat + bg * (1- yrat))
            b = int(tb * yrat + bb * (1- yrat))
            pix[5 + i * 2, gridsize + 5 + j * 2] = (r,g,b)
            pix[6 + i * 2, gridsize + 5 + j * 2] = (r,g,b)
            pix[5 + i * 2, gridsize + 6 + j * 2] = (r,g,b)
            pix[6 + i * 2, gridsize + 6 + j * 2] = (r,g,b)
    im.save(args.output_file)


parser = argparse.ArgumentParser(prog='coif.py', formatter_class=argparse.RawDescriptionHelpFormatter, description='''
RMaD.py - read mapping visualisation in the large

USAGE: RMaD.py -bam bamfile.bam -o output_file.bmp -size 5000
          Create a bmp file from a bamfile of paired-end reads with a width and height of 5000px
       RMaD.py -r reads.fa -B blast_prefix -r reference -o output_file.png -bin bin_size
          Create a png file from reads.fa, generate blast file. Image size will be reference length / bin_size
''', epilog="Thanks for using RMaD.py")
parser.add_argument('-r', '--read_file', action='store', default=None, help='read file')
parser.add_argument('-ref', '--reference_file', action='store', default=None, help='reference file')
parser.add_argument('-bam', '--bam_file', action='store', default=None, help='bam file')
parser.add_argument('-sam', '--sam_file', action='store', default=None, help='sam file')
parser.add_argument('-B', '--gen_blast', action='store', default=None, help='Generate blast files, use argument as prefix for output.')
parser.add_argument('-b', '--blast_file', action='store', default=None, help='Blast file (output format 6)')
parser.add_argument('-o', '--output_file', action='store', required=True, help='output file [gif/bmp/png]')
parser.add_argument('-s', '--size', action='store', type=int, default=None, help='Image size')
parser.add_argument('-bin', '--bin_size', action='store', type=int, default=None, help='Bin size')
parser.add_argument('-g', '--gap', action='store', type=int, default=5, help='Gap size')
parser.add_argument('-sub', '--subsection', nargs='+', action='store', default=None, help='Only display subection of genome [ref min_cutoff max_cutoff')
parser.add_argument('-m', '--scale', action='store', default=None, help='Draw scale lines [bp]')
parser.add_argument('-c', '--min_hits', action='store', type=int, default=1, help='Min hits to be shown')


args = parser.parse_args()
if args.size is None and args.bin_size is None:
    sys.stderr.write('Please give a image size or bin size.')
    sys.exit()
if not args.size is None and not args.bin_size is None:
    sys.stderr.write('Only provide bin size or image size, not both.')
    sys.exit()
if not args.sam_file is None or not args.bam_file is None:
    read_sbam(args)
elif not args.gen_blast is None:
    if args.reference_file is None:
        sys.stderr.write('Please provide a reference file')
        sys.exit()
    if args.read_file is None:
        sys.stderr.write('Please provide a read file (FASTA)')
        sys.exit()
    generate_blast(args)
if args.blast_file is None:
    sys.stderr.write('Please either generate or provide a BLAST comparison')
    sys.exit()
else:
    read_sing(args)
draw_dotplot(args)