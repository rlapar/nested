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
    blastx_cline = NcbiblastxCommandline(db=config.blastx_gydb_protein_db, outfmt=5, 
        num_threads=config.blastx_args['num_threads'], dbsize=config.blastx_args['dbsize'], 
        word_size=config.blastx_args['word_size'], evalue=config.blastx_args['evalue'])
    
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
                if domain_type not in domains.keys():
                    domains['XX'].append(domain)
                else:
                    domains[domain_type].append(domain)
        
    return domains