#!/usr/bin/env python3

from nested.core.te import run_ltr_finder
from nested.core.domain import run_blastx
from nested.scoring.graph import Graph
from nested.utils import intervals

class Gene(object):
	"""Class representing a genereic gene in the program. Used in recursive calls, contain information about the gene domains and found TE's

	Attributes:
		seqid (str): id of sequence 
		sequence (Bio.Seq.Seq): sequence
		te_list (list[TE]): list of pairs of LTR's found from LTR finder
		domain_list (list[Domain]): list of found domains
	
	"""
	def __init__(self, seqid=None, sequence=None):
		self.seqid = seqid
		self.sequence = sequence
		self.te_list = run_ltr_finder(seqid, sequence)
		self.domain_list = run_blastx(sequence)

	def __str__(self):
		strlen = 15
		lines = ['{{id = {},'.format(self.seqid),
				 ' sequence = {a}{b},'.format(a=self.sequence[:strlen], b='...' if len(self.sequence) > strlen else ''),
				 ' teList.size = {},'.format(len(self.te_list)),
				 ' domainList.size = {}}}'.format(len(self.domain_list))]
		return '\n'.join(lines)

	def get_best_candidate(self):
		"""Evaluate all LTR pairs and return best scored pair, or None if no pair found

		Returns:
			TE: best evaluated pair
		"""
		scores = []
		for te in self.te_list:
			scores.append(self._evaluate_te(te))

		if not scores:
			return None

		return self.te_list[scores.index(max(scores))]

	def _evaluate_te(self, te):
		#evaluate LTR pair using graph, set up score and features
		graph = Graph(te, self.domain_list)
		te.score, te.features = graph.get_score()
		te.score /= float(intervals.length(te.location))
		return te.score