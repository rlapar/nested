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
		scores (list[float]): list of evaluations of each TE
	
	"""
	def __init__(self, seqid=None, sequence=None):
		self.seqid = seqid
		self.sequence = sequence
		self.te_list = run_ltr_finder(seqid, sequence)
		self.domain_list = run_blastx(sequence)
		self.scores = []
		self._evaluate_te_list()

	def __str__(self):
		strlen = 15
		lines = ['{{id = {a}{b},'.format(a=self.seqid[:strlen], b='...' if len(self.sequence) > strlen else ''),
				 ' sequence = {a}{b},'.format(a=self.sequence[:strlen], b='...' if len(self.sequence) > strlen else ''),
				 ' te_list.size = {},'.format(len(self.te_list)),
				 ' domain_list.size = {}}}'.format(len(self.domain_list))]
		return '\n'.join(lines)

	def get_best_candidate(self):
		"""Evaluate all LTR pairs and return best scored pair, or None if no pair found

		Returns:
			TE: best evaluated pair
		"""
		if not self.scores:
			return None
		return self.te_list[self.scores.index(max(self.scores))]

	def get_candidates_above_threshold(self, threshold=0.15):
		candidates = []
		for i in range(len(self.te_list)):
			if self.scores[i] >= threshold:
				candidates.append(self.te_list[i])
		return candidates

	def _evaluate_te_list(self):
		for te in self.te_list:
			self.scores.append(self._evaluate_te(te))

	def _evaluate_te(self, te):
		#evaluate LTR pair using graph, set up score and features
		graph = Graph(te, self.domain_list)
		te.score, te.features = graph.get_score()
		te.score /= float(intervals.length(te.location))
		return te.score