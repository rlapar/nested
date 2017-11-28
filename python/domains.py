from StringIO import StringIO

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline

from python import config

def runBlastx(genes):
    #print 'Running blastx...'
    for g in genes:
        genes[g]['domains'] = findDomains(genes[g])
        
    return genes

def findDomains(gene):
    domains = {'INT':[],'GAG':[],'AP':[],'RNaseH':[],'RT':[],'XX':[]}
    blastx_cline = NcbiblastxCommandline(**config.blastx_args)
    
    xml_out, stderr = blastx_cline(stdin=str(gene['sequence']))
    blast_records = NCBIXML.parse(StringIO(xml_out))
    try:
        blast_record = next(blast_records).alignments
    except ValueError:
        return domains
    #print len(blast_record)
    if blast_record:
        for alignment in blast_record:
            for hsp in alignment.hsps:
                domain = {
                    'evalue': hsp.expect,
                    'frame': hsp.frame,
                    'score': hsp.score,
                    'title': alignment.title
                }
                if hsp.frame[0] < 0: #complementary strand
                    domain['location'] = [hsp.query_start, hsp.query_start - 3 * hsp.align_length]
                else: #forward strand
                    domain['location'] =[hsp.query_start, hsp.query_start + 3 * hsp.align_length]               
                
                domain_type = alignment.title.split(' ')[1].split('_')[0]
                domain['type'] = domain_type
                if domain_type not in domains.keys():
                    domains['XX'].append(domain)
                else:
                    domains[domain_type].append(domain)
        
    return domains