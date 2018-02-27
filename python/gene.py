#!/usr/bin/env python3

import copy

from python import te as pythonTe
from python import domain as pythonDomain
from python.transposonGraph import TransposonGraph

class Gene(object):
	"""Class representing a genereic gene in the program. Used in recursive calls, contain information about the gene domains and found TE's

	Attributes:
		seqid (str): id of sequence 
		sequence (Bio.Seq.Seq): sequence
		teList (list[TE]): list of pairs of LTR's found from LTR finder
		domainList (list[Domain]): list of found domains
	
	"""
	def __init__(self, seqid=None, sequence=None):
		self.seqid = seqid
		self.sequence = sequence
		self.teList = pythonTe.runLtrFinder(seqid, sequence)
		self.domainList = pythonDomain.runBlastx(sequence)

	def __str__(self):
		strlen = 15
		lines = ['{{id = {},'.format(self.seqid),
				 ' sequence = {a}{b},'.format(a=self.sequence[:strlen], b='...' if len(self.sequence) > strlen else ''),
				 ' teList.size = {},'.format(len(self.teList)),
				 ' domainList.size = {}}}'.format(len(self.domainList))]
		return '\n'.join(lines)

	def getBestCandidate(self):
		"""Evaluate all LTR pairs and return best scored pair

		Returns:
			TE: best evaluated pair
		"""
		scores = []
		for te in self.teList:
			scores.append(self._evaluateTe(te))

		if not scores:
			return None

		return self.teList[scores.index(max(scores))]		

	def _evaluateTe(self, te): #evaluate LTR pair using graph, set up score and features
		graph = TransposonGraph(te, self.domainList)
		te.score, te.features = graph.getScore()
		te.score /= float(te.location[1] - te.location[0])
		return te.score
