__author__ = 'mjsul_000'

def readSam(sam_file):
    self.header = ''
    self.references = []
    self.lengths = []
    blast = ntup('qname rname pos mapq cigar rnext pnext tlen seq qual')
    with open(sam_file) as sam:
        for line in sam:
            if line.startswith('@'):
                self.header += line
                if line.startswith('SQ'):
                    for i in line.split():
                        if i.startswith('SN:'):
                            self.references.append(i[3:])
                        elif i.startswith('LN:')
                            self.lengths.append(int(i[3:]))
            yield line
