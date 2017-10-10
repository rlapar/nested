from StringIO import StringIO

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline

from python import config

def findDomains(gene):
    blastx_cline = NcbiblastxCommandline(db=config.blastx_gydb_protein_db, outfmt=5, 
        num_threads=config.blastx_args['num_threads'], dbsize=config.blastx_args['dbsize'], 
        word_size=config.blastx_args['word_size'], evalue=config.blastx_args['evalue'])
    
    xml_out, stderr = blastx_cline(stdin=str(gene['sequence']))
    
    blast_records = NCBIXML.parse(StringIO(xml_out))
    blast_record = next(blast_records).alignments
    #print len(blast_record)
    domains = []
    if blast_record:
        for alignment in blast_record:
            for hsp in alignment.hsps:
                domain = {
                    'E-value': hsp.expect,
                    'Frame': hsp.frame,
                    'Score': hsp.score,
                    'Title': alignment.title,
                    'Type': alignment.title.split(' ')[1]
                }
                if hsp.frame[0] < 0: #complementary strand
                    domain['Location'] = [hsp.query_start, hsp.query_start - 3 * hsp.align_length]
                else: #forward strand
                    domain['Location'] =[hsp.query_start, hsp.query_start + 3 * hsp.align_length]               
                
                domains.append(domain)
        
    return domains