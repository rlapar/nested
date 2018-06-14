#!/usr/bin/env python3

from io import StringIO

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline

from nested.config.config import config


class Domain(object):
    """Class representing domain parsed from blastx output.

    Atrributes:
        evalue (float): evalue (same as blastx output)
        frame (tuple): frame (same as blastx output)
        score (float): score (same as blastx output)
        title (str): title (same as blastx output)
        location (list): location [from, to] (from > to for negative strand)
        type (str): type of domain (AP, GAG, INT, ...)
    """

    def __init__(self, evalue=None, frame=None, score=None, title=None, location=None, domainType=None):
        self.evalue = evalue
        self.frame = frame
        self.score = score
        self.title = title
        self.location = location
        self.type = domainType

    def __str__(self):
        strlen = 15
        lines = ['{{title = {a}{b},'.format(a=self.title[:strlen], b='...' if len(self.title) > strlen else ''),
                 ' type = {},'.format(self.type),
                 ' location = {},'.format(self.location),
                 ' evalue = {},'.format(self.evalue),
                 ' frame = {},'.format(self.frame),
                 ' score = {}}}'.format(self.score)]
        return '\n'.join(lines)


def run_blastx(sequence):
    """Run blastx and get the list of found domains

    Arguments:
        sequence (Bio.Seq.Seq): sequence

    Returns:
        list[Domain]: list of matched domains
    """
    domains = []
    blastx_cline = NcbiblastxCommandline(**config['blastx']['args'])
    xml_out, stderr = blastx_cline(stdin=str(sequence))
    blast_records = NCBIXML.parse(StringIO(xml_out))
    try:
        blast_record = next(blast_records).alignments
    except ValueError:
        return domains

    if blast_record:
        for alignment in blast_record:
            for hsp in alignment.hsps:
                domain = Domain(evalue=hsp.expect,
                                frame=hsp.frame,
                                score=hsp.score,
                                title=alignment.title)

                if hsp.frame[0] < 0:  # complementary strand
                    domain.location = [hsp.query_start, hsp.query_start - 3 * hsp.align_length]
                else:
                    domain.location = [hsp.query_start, hsp.query_start + 3 * hsp.align_length]

                domain.type = alignment.title.split(' ')[1].split('_')[0]
                domains.append(domain)

    return domains
