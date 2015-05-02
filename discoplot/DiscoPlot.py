import argparse
import numpy
import sys
import subprocess
import random
from collections import namedtuple as ntup

class readSam:
    def __init__(self, sam_file):
        self.header = ''
        self.references = []
        self.lengths = []
        self.read = ntup('blast', 'pos pnext rname rnext is_reverse mate_is_reverse is_read1 is_unmapped mate_is_unmapped line')
        self.sam = open(sam_file)
        line = self.sam.readline()
        lastline = None
        while line.startswith('@'):
            self.header += line
            if line.startswith('@SQ'):
                for i in line.split():
                    if i.startswith('SN:'):
                        self.references.append(i[3:])
                    elif i.startswith('LN:'):
                        self.lengths.append(int(i[3:]))
            lastline = line
            line = self.sam.readline()
        self.sam.seek(0)
        if not lastline is None:
            getit = True
            while getit:
                line = self.sam.readline()
                if line == lastline:
                    getit = False

    def __iter__(self):
        return self

    def next(self):
        line = self.sam.readline()
        if line == '':
            raise StopIteration
        name, flag, rname, pos, mapq, cigar, rnext, pnext = line.split()[:8]
        if rnext == '=' or rnext == '*':
            rnext = rname
        flag = bin(int(flag)).zfill(12)
        read = self.read(int(pos), int(pnext), rname, rnext, flag[-5] == '1', flag[-6] == '1', flag[-7] == '1', flag[-3] == '1', flag[-4] == '1', line)
        return read

class writeSam:
    def __init__(self, samfile, header):
        self.out = open(samfile, 'w')
        self.out.write(header)

    def write(self, read):
        self.out.write(read.line)



def read_sbam(args):
    try:
        import pysam
        havepysam = True
    except ImportError:
        havepysam = False
    if not args.bam_file is None:
        sam = pysam.Samfile(args.bam_file, 'rb')
    elif not args.sam_file is None and havepysam:
        sam = pysam.Samfile(args.sam_file)
    elif not args.sam_file is None:
        sam = readSam(args.sam_file)
    else:
        sys.stderr.write('Please install pysam to read bam files (pysam not needed for sam files).')
        return
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
    if not args.write_reads is None:
        if havepysam:
            if args.bam_file is None:
                newsam = pysam.Samfile(args.write_reads[4], 'w', template=sam)
            else:
                newsam = pysam.Samfile(args.write_reads[4], 'wb', template=sam)
        else:
            newsam = writeSam(args.write_reads[4], sam.header)
        if len(args.write_reads) == 5:
            refpos[sam.references[0]] = 0
            cutw = int(args.write_reads[0])
            cutx = int(args.write_reads[1])
            cuty = int(args.write_reads[2])
            cutz = int(args.write_reads[3])
    for read in sam:
        if not args.write_reads is None:
            pos1 = read.pos
            pos2 = read.pnext
            if read.is_reverse:
                if read.mate_is_reverse:
                    if pos1 <= pos2:
                        if cutw <= pos2 <= cutx and cuty <= pos1 <= cutz:
                            newsam.write(read)
                    else:
                        if cutw <= pos1 <= cutx and cuty <= pos2 <= cutz:
                            newsam.write(read)
                else:
                    if cutw <= pos2 <= cutx and cuty <= pos1 <= cutz:
                        newsam.write(read)

            else:
                if read.mate_is_reverse:
                    if cutw <= pos1 <= cutx and cuty <= pos2 <= cutz:
                        newsam.write(read)
                else:
                    if pos1 <= pos2:
                        if cutw <= pos1 <= cutx and cuty <= pos2 <= cutz:
                            newsam.write(read)
                    else:
                        if cutw <= pos2 <= cutx and cuty <= pos1 <= cutz:
                            newsam.write(read)
        if havepysam:
            if read.tid >= 0:
                ref = sam.getrname(read.tid)
            else:
                ref = ''
        else:
            ref = read.rname
        if ref in refpos:
            if read.is_read1:
                if cuta <= read.pos <= cutb and not read.is_unmapped:
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
                        if havepysam:
                            mref = sam.getrname(read.rnext)
                        else:
                            mref = read.rnext
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
                    if havepysam:
                        if read.tid >= 0:
                            ref = sam.getrname(read.tid)
                        else:
                            ref = ''
                    else:
                        ref = read.rname
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
    global cuta
    global cutb
    global refpos
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
                name = line.split()[0][1:]
                seq = ''
            elif line.startswith('>'):
                readlen[name] = len(seq)
                name = line.split()[0][1:]
                seq = ''
            elif getfq == 0:
                seq += line.rstrip()
            elif getfq == 1:
                readlen[name] = len(seq)
                name = line.split()[0][1:]
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
                name = line.split()[0][1:]
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
    maxrstart = 0
    for line in blast:
        query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
        qstart, qstop, rstart, rstop, length, mm, indel = map(int, [qstart, qstop, rstart, rstop, length, mm, indel])
        if rstart >= maxrstart:
            maxrstart = rstart
        if query != lastquery and lastquery != '' and hits != []:
            hits.sort(reverse=True)
            newhits = [hits[0]]
            qtaken = set()
            for i in range(hits[0][2], hits[0][3] + 1):
                qtaken.add(i)
            hitsizes = set()
            hitsizes.add(hits[0][:4])
            for i in hits[1:]:
                if i[:-3] == newhits[-1][:-3]:
                    newhits.append(i)
                    hitsizes.add(i[:4])
                else:
                    getit = False
                    for j in range(i[2], i[3] + 1):
                        if not j in qtaken:
                            getit = True
                            qtaken.add(j)
                    if getit:
                        newhits.append(i)
                        hitsizes.add(i[:4])
            if len(hitsizes) == 1 and len(newhits) != 1:
                newhits = [random.choice(newhits)]
            anchor = None
            revseq = None
            newhits2 = []
            for i in newhits:
                 if subject in refpos and ((cuta <= i[4] <= cutb) or (cuta <= i[5] <= cutb)):
                     newhits2.append(i)
            for i in newhits2:
                bitscore2, length2, qstart2, qstop2, rstart2, rstop2, subject2 = i
                if anchor is None:
                    if rstart2 < rstop2:
                        anchor = refpos[subject2] * args.bin_size + rstart2 - qstart2
                        revseq = False
                    else:
                        anchor = refpos[subject2] * args.bin_size + rstop2 - (readlen[lastquery] - qstop2)
                        revseq = True
                    if min(qtaken) >= args.unmapped and min(qtaken) == qstart2:
                        if revseq:
                            if (anchor + readlen[lastquery] - qstart2 - cuta)/args.bin_size + refpos[subject2] in unmapped_for:
                                unmapped_for[(anchor + readlen[lastquery] - qstart2 - cuta)/args.bin_size + refpos[subject2]] += 1
                            else:
                                unmapped_for[(anchor + readlen[lastquery] - qstart2 - cuta)/args.bin_size + refpos[subject2]] = 1
                        else:
                            if (anchor + qstart2 - cuta)/args.bin_size in unmapped_rev:
                                unmapped_rev[(anchor + qstart2 - cuta)/args.bin_size] += 1
                            else:
                                unmapped_rev[(anchor + qstart2 - cuta)/args.bin_size] = 1
                    if max(qtaken) <= readlen[lastquery] - args.unmapped and max(qtaken) == qstop2:
                        if revseq:
                            if (anchor + readlen[lastquery] - qstop2 - cuta)/args.bin_size + refpos[subject2] in unmapped_rev:
                                unmapped_rev[(anchor + readlen[lastquery] - qstop2 - cuta)/args.bin_size + refpos[subject2]] += 1
                            else:
                                unmapped_rev[(anchor + readlen[lastquery] - qstop2 - cuta)/args.bin_size + refpos[subject2]] = 1
                        else:
                            if (anchor + qstop2 - cuta)/args.bin_size in unmapped_for:
                                unmapped_for[(anchor + qstop2 - cuta)/args.bin_size] += 1
                            else:
                                unmapped_for[(anchor + qstop2 - cuta)/args.bin_size] = 1
                lastxpos = None
                lastypos = None
                oldstart, oldstop = qstart2, qstop2
                if revseq:
                    rstart2, rstop2 = rstop2, rstart2
                    qstart2 = readlen[lastquery] - qstop2
                    qstop2 = readlen[lastquery] - oldstart
                for j in range(qstart2, qstop2):
                    xpos = (anchor + j - cuta) / args.bin_size
                    ypos = refpos[subject2] + ((rstart2 + int(((j - qstart2) * 1.0 / (qstop2 - qstart2)) * (rstop2 - rstart2))) - cuta) / args.bin_size
                    if xpos != lastxpos or ypos != lastypos:
                        if rstart2 < rstop2:
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
        if query != lastquery:
            hits = []
        if ident >= args.min_ident and length >= args.min_length:
            hits.append((float(bitscore), length, qstart, qstop, rstart, rstop, subject))
        lastquery = query

def read_hist(args):
    histFile = open(args.heatmap)
    global unmapped_for
    global unmapped_rev
    global dirgrid
    global invgrid
    global cuta
    global cutb
    global refpos
    unmapped_for, unmapped_rev, dirgrid, invgrid = {}, {}, {}, {}
    cuta = 0
    cutb = float('inf')
    refpos = {}
    header = True
    args.bin_size = None
    for line in histFile:
        if not line.startswith('#'):
            if header:
                header = False
                headpos = []
                args.size = len(line.split())
                for i in line.split():
                    headpos.append(int(i.split(':')[1].split('-')[0]))
                args.bin_size = headpos[1] - headpos[0]
            else:
                name = line.split()[0]
                vals = line.split()[1:]
                pos = int(name.split(':')[1].split('-')[0])
                for i in range(len(vals)):
                    if vals[i] != 0.0:
                        if pos/args.bin_size in dirgrid:
                            dirgrid[pos/args.bin_size][headpos[i]/args.bin_size] = float(vals[i])
                        else:
                            dirgrid[pos/args.bin_size] = {headpos[i]/args.bin_size:float(vals[i])}

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
    numvals = 0
    vals = []
    diagdict = {}
    for i in invgrid:
        for j in invgrid[i]:
            if args.max_hits >= invgrid[i][j] >= args.min_hits:
                numvals += 1
                if i - j in diagdict:
                    diagdict[i-j].append(invgrid[i][j])
                else:
                    diagdict[i-j] = [invgrid[i][j]]
    thesum = 0
    for i in diagdict:
        thesum += sum(diagdict[i])
    for i in dirgrid:
        for j in dirgrid[i]:
            if args.max_hits >= dirgrid[i][j] >= args.min_hits:
                numvals += 1
                if i - j in diagdict:
                    diagdict[i-j].append(dirgrid[i][j])
                else:
                    diagdict[i-j] = [dirgrid[i][j]]
    thesum2 = 0
    for i in diagdict:
        thesum2 += sum(diagdict[i])
    if thesum2 / 2 < thesum:
        args.switch = not args.switch
    maxsum = 0
    for i in diagdict:
        if sum(diagdict[i]) > maxsum:
            maxsum = sum(diagdict[i])
            vals = diagdict[i]
    numvals2 = 0
    for i in unmapped_rev:
        if args.max_hits >= unmapped_rev[i] >= args.min_hits:
            numvals2 += 1
    for i in unmapped_for:
        if args.max_hits >= unmapped_for[i] >= args.min_hits:
            numvals2 += 1
    if args.log:
        vals = map(numpy.log10, vals)
    if args.m_count != -1:
        med = args.m_count
    else:
        med = numpy.median(vals)
    sizemod = (864.0 / args.size * args.m_size) ** 2 / med
    themew = 864.0 / args.size * args.marker_edge_width
    if args.split_graph is None:
        ax = plt.subplot(aspect=1)
    else:
        if args.split_graph[0].isdigit():
            start = True
            starts = []
            stops = []
            for i in args.split_graph:
                if start:
                    starts.append(int(i))
                else:
                    stops.append(int(i))
                start = not start
        else:
            count = 0
            starts = []
            stops = []
            for i in args.split_graph():
                if count % 3 == 0:
                    name = i
                elif count % 3 == 1:
                    starts.append(refpos[name] * args.bin_size + int(i))
                else:
                    stops.append(refpos[name] * args.bin_size + int(i))
        widths = [a - b for a, b in zip(stops, starts)]
        heights = widths[::-1]
        gs = gridspec.GridSpec(len(starts), len(starts), width_ratios=widths, height_ratios=heights)
        axgrid = [[None for i in range(len(starts))] for i in range(len(starts))]
        for i in range(len(starts) * len(starts)):
            if i % len(starts) == 0:
                axgrid[i%len(starts)][i/len(starts)] = plt.subplot(gs[i], aspect=1)
            else:
                axgrid[i%len(starts)][i/len(starts)] = plt.subplot(gs[i], aspect=1)#, sharey=axgrid[0][i/len(starts)])
    if not args.highlight is None:
        hstarts = []
        hstops = []
        halpha = float(args.highlight[0])
        if args.highlight[1].isdigit():
            start = True
            for i in args.highlight[1:]:
                if start:
                    hstarts.append(int(i))
                else:
                    hstops.append(int(i))
                start = not start
        else:
            count = 0
            for i in args.highlight[1:]:
                if count % 3 == 0:
                    name = i
                elif count % 3 == 1:
                    hstarts.append(refpos[name] * args.bin_size + int(i))
                else:
                    hstops.append(refpos[name] * args.bin_size + int(i))
        if args.split_graph is None:
            for i in range(len(hstarts)):
                ax.axhspan(hstarts[i], hstops[i], facecolor='g', alpha=halpha)
                ax.axvspan(hstarts[i], hstops[i], facecolor='g', alpha=halpha)

    x = numpy.zeros(numvals, dtype='u4')
    y = numpy.zeros(numvals, dtype='u4')
    sizes = numpy.zeros(numvals, dtype='f4')
    colours = numpy.array(['x' for i in range(numvals)])
    count = 0
    if args.switch:
        for i in invgrid:
            for j in invgrid[i]:
                if args.max_hits >= invgrid[i][j] >= args.min_hits:
                    x[count] = i * args.bin_size + cuta
                    y[count] = j * args.bin_size + cuta
                    if args.log:
                        sizes[count] = numpy.log10(invgrid[i][j]) * sizemod
                    else:
                        sizes[count] = invgrid[i][j] * sizemod
                    colours[count] = 'b'
                    count += 1
    for i in dirgrid:
        for j in dirgrid[i]:
            if args.max_hits >= dirgrid[i][j] >= args.min_hits:
                x[count] = i * args.bin_size + cuta
                y[count] = j * args.bin_size + cuta
                if args.log:
                    sizes[count] = numpy.log10(dirgrid[i][j]) * sizemod
                else:
                    sizes[count] = dirgrid[i][j] * sizemod
                colours[count] = 'r'
                count += 1
    if not args.switch:
        for i in invgrid:
            for j in invgrid[i]:
                if args.max_hits >= invgrid[i][j] >= args.min_hits:
                    x[count] = i * args.bin_size + cuta
                    y[count] = j * args.bin_size + cuta
                    if args.log:
                        sizes[count] = numpy.log10(invgrid[i][j]) * sizemod
                    else:
                        sizes[count] = invgrid[i][j] * sizemod
                    colours[count] = 'b'
                    count += 1
    if args.split_graph is None:
        ax.scatter(x, y, s=sizes, c=colours, alpha=args.alpha, marker='x', lw=themew)
    else:
        for i in range(len(starts)):
            for j in range(len(starts)):
                axgrid[i][j].scatter(x, y, s=sizes, c=colours, alpha=args.alpha, marker='x', lw=themew)
    count = 0
    x = numpy.zeros(numvals2, dtype='u4')
    y = numpy.zeros(numvals2, dtype='u4')
    sizes = numpy.zeros(numvals2, dtype='f4')
    colours = numpy.array(['x' for i in range(numvals2)])
    for i in unmapped_for:
        if args.max_hits >= unmapped_for[i] >= args.min_hits:
            x[count] = cuta
            y[count] = i * args.bin_size + cuta
            if args.log:
                sizes[count] = numpy.log10(unmapped_for[i]) * sizemod
            else:
                sizes[count] = unmapped_for[i] * sizemod
            colours[count] = 'g'
            count += 1
    for i in unmapped_rev:
        if args.max_hits >= unmapped_rev[i] >= args.min_hits:
            x[count] = i * args.bin_size + cuta
            y[count] = cuta
            if args.log:
                sizes[count] = numpy.log10(unmapped_rev[i]) * sizemod
            else:
                sizes[count] = unmapped_rev[i] * sizemod
            colours[count] = 'g'
            count += 1
    if args.split_graph is None:
        ax.scatter(x, y, s=sizes, c=colours, alpha=args.alpha2, marker='+', lw=themew)
    else:
        for i in range(len(starts)):
            for j in range(len(starts)):
                axgrid[i][j].scatter(x, y, s=sizes, c=colours, alpha=args.alpha2, marker='+', lw=themew)
    sizes = []
    names = []
    for i in [10, 25, 50, 75, 90]:
        if args.log:
            sizes.append(10**numpy.percentile(vals, i))
        else:
            sizes.append(numpy.percentile(vals, i))
        names.append(str(i) + '% Normal ' + str(sizes[-1]))
    names.append('50% Inverted ' + str(sizes[2]))
    if args.log:
        a = plt.scatter([], [], s=numpy.log10(sizes[2]) * sizemod, c='b', marker='x', lw=themew)
        b = plt.scatter([], [], s=numpy.log10(sizes[0]) * sizemod, c='r', marker='x', lw=themew)
        c = plt.scatter([], [], s=numpy.log10(sizes[1]) * sizemod, c='r', marker='x', lw=themew)
        d = plt.scatter([], [], s=numpy.log10(sizes[2]) * sizemod, c='r', marker='x', lw=themew)
        e = plt.scatter([], [], s=numpy.log10(sizes[3]) * sizemod, c='r', marker='x', lw=themew)
        f = plt.scatter([], [], s=numpy.log10(sizes[4]) * sizemod, c='r', marker='x', lw=themew)
    else:
        a = plt.scatter([], [], s=sizes[2] * sizemod, c='b', marker='x', lw=themew)
        b = plt.scatter([], [], s=sizes[0] * sizemod, c='r', marker='x', lw=themew)
        c = plt.scatter([], [], s=sizes[1] * sizemod, c='r', marker='x', lw=themew)
        d = plt.scatter([], [], s=sizes[2] * sizemod, c='r', marker='x', lw=themew)
        e = plt.scatter([], [], s=sizes[3] * sizemod, c='r', marker='x', lw=themew)
        f = plt.scatter([], [], s=sizes[4] * sizemod, c='r', marker='x', lw=themew)
    if args.split_graph is None:
        if args.no_legend:
            leg = ax.legend([b, c, d, e, f, a], names, loc=4)
            leg.draggable(state=True)
        for i in refpos:
            if not refpos[i] == 0:
                ax.axhspan(refpos[i] * args.bin_size, refpos[i] * args.bin_size - args.gap * args.bin_size, facecolor='g', alpha=args.alpha)
                ax.axvspan(refpos[i] * args.bin_size, refpos[i] * args.bin_size - args.gap * args.bin_size, facecolor='g', alpha=args.alpha)
    if cutb == float('inf'):
        cutb = args.size * args.bin_size + cuta
    if args.split_graph is None:
        plt.xlim([cuta - args.bin_size, cutb])
        plt.ylim([cuta - args.bin_size, cutb])
        if args.no_gridlines:
            plt.grid(True)
        if args.no_label:
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
    else:
        gs.update(wspace=0.05, hspace=0.05)#, left=0.05, top=0.85, bottom=0.05)
        for i in range(len(starts)):
            for j in range(len(starts)):
                axgrid[i][j].set_xlim((starts[i], stops[i]))
                axgrid[i][j].set_ylim((starts[-1-j], stops[-1-j]))
                if j != len(starts) - 1:
                    axgrid[i][j].set_xticklabels([])
                if i != 0:
                    axgrid[i][j].set_yticklabels([])
                if args.no_gridlines:
                    axgrid[i][j].grid(True)
                if args.no_label:
                    axgrid[i][j].xaxis.set_visible(False)
                    axgrid[i][j].yaxis.set_visible(False)
    if not args.output_file is None:
        fig = plt.gcf()
        fig.set_size_inches(12,12)
        plt.savefig(args.output_file, dpi=args.image_quality)
    else:
        plt.show()



parser = argparse.ArgumentParser(prog='DiscoPlot.py', formatter_class=argparse.RawDescriptionHelpFormatter, description='''
DiscoPlot.py - read mapping visualisation in the large

USAGE: DiscoPlot.py -bam bamfile.bam -o output_file.bmp -size 5000
          Create a bmp file from a bamfile of paired-end reads with 5000 bins
       DiscoPlot.py -r reads.fa -B blast_prefix -r reference -o output_file.png -bin 10000
          Create a png file using reads.fa aligned to the reference, automatically generate blast file. Use a bin size of 10,000bp.
''', epilog="Thanks for using DiscoPlot.py")
parser.add_argument('-r', '--read_file', action='store', default=None, help='read file - provide DiscoPlot with a read file to BLAST (long read mode).')
parser.add_argument('-ref', '--reference_file', action='store', default=None, help='Reference file - Reference for generating long reads alignments.')
parser.add_argument('-bam', '--bam_file', action='store', default=None, help='bam file - paired read mode. (Requires pysam).')
parser.add_argument('-sam', '--sam_file', action='store', default=None, help='sam file - paired read mode. (pysam not required)')
parser.add_argument('-hm', '--heatmap', action='store', default=None, help='Heatmap file - provide DiscoPlot with custom generated heatmap.')
parser.add_argument('-B', '--gen_blast', action='store', default=None, help='Generate blast files, use argument as prefix for output.')
parser.add_argument('-b', '--blast_file', action='store', default=None, help='Provide DiscoPlot with alignment file (long read mode) (BLAST tab delimited file - output format 6)')
parser.add_argument('-o', '--output_file', action='store', default=None, help='output file [gif/bmp/png]')
parser.add_argument('-s', '--size', action='store', type=int, default=None, help='Number of bins.')
parser.add_argument('-bin', '--bin_size', action='store', type=int, default=None, help='Bin size (in bp)')
parser.add_argument('-g', '--gap', action='store', type=int, default=5, help='Gap size - gap size between entries in reference.')
parser.add_argument('-sub', '--subsection', nargs='+', action='store', default=None, help='Only display subection of genome [ref]/[min_cutoff max_cutoff]/[ref min_cutoff max_cutoff]')
parser.add_argument('-wb', '--write_reads', nargs='+', action='store', default=None, help='Write reads in rectangle to bam/sam [x1 y1 x2 y2 out.bam]')
parser.add_argument('-c', '--min_hits', action='store', type=int, default=1, help='Only show bins with more than this number of hits.')
parser.add_argument('-m', '--max_hits', action='store', type=float, default=float('inf'), help='Only show bins with less hits than this.')
parser.add_argument('-dpi', '--image_quality', action='store', type=int, default=1600, help='Image quality (in DPI)')
parser.add_argument('-i', '--min_ident', action='store', type=float, default=85.0, help='Min. idenity of hits to draw (long read mode).')
parser.add_argument('-l', '--min_length', action='store', type=int, default=50, help='Min. length of hits to draw (long read mode).')
parser.add_argument('-d', '--unmapped', action='store', type=int, default=100, help='Unmapped bases on edge for RMaD to consider read partially unmapped.')
parser.add_argument('-a', '--alpha', action='store', type=float, default=0.1, help='Transparency of mapped read markers')
parser.add_argument('-a2', '--alpha2', action='store', type=float, default=0.8, help='Transparency of unmapped read markers')
parser.add_argument('-mc', '--m_count', action='store', type=int, default=-1, help='The count of a bin to be used as the median value for calculating the size of the dot [auto]')
parser.add_argument('-ms', '--m_size', action='store', type=float, default=20, help='Set the width (in bins) of a marker with a median count.')
parser.add_argument('-log', '--log', action='store_true', default=False, help='Log10 bin counts. (For data with highly variable coverage).')
parser.add_argument('-sw', '--switch', action='store_true', default=False, help='Draw most common (inverted/direct) hits first.')
parser.add_argument('-nl', '--no_legend', action='store_false', default=True, help='Don\'t create legend.')
parser.add_argument('-ng', '--no_gridlines', action='store_false', default=True, help='Don\'t draw gridlines.')
parser.add_argument('-na', '--no_label', action='store_true', default=False, help='No axis labels.')
parser.add_argument('-split', '--split_graph', nargs='+', action='store', default=None, help='Show multiple subsections of graph [start1 stop1 start2 stop2 etc.] or [ref1 start1 stop1 ref2 start2 stop2 etc.]')
parser.add_argument('-hl', '--highlight', nargs='+', action='store', default=None, help='Highlight subsections of graph [alpha start1 stop1 start2 stop2 etc.] or [alphref1 start1 stop1 ref2 start2 stop2 etc.]')
parser.add_argument('-mw', '--marker_edge_width', action='store', type=int, default=20, help='Marker width (default is roughly 20x bin size)')



args = parser.parse_args()
if args.size is None and args.bin_size is None and args.heatmap is None:
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
import matplotlib.gridspec as gridspec

if not args.size is None and not args.bin_size is None:
    sys.stderr.write('Only provide bin size or image size, not both.')
    sys.exit()
if not args.sam_file is None or not args.bam_file is None:
    read_sbam(args)
elif not args.heatmap is None:
    read_hist(args)
elif args.blast_file is None:
    sys.stderr.write('Please either generate or provide a BLAST comparison')
    sys.exit()
else:
    read_sing(args)
if args.write_reads is None:
    draw_dotplot(args)