#!/usr/bin/env python3

from python import te as pythonTe
from python import domain as pythonDomain

class Gene(object):
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
		scores = []
		for te in self.teList:
			scores.append(self._evaluateTe(te))

		if not scores:
			return None

		return self.teList[scores.index(max(scores))]		

	def _evaluateTe(self, te):
		#TODO scores and features (domains)
		#TODO genegraph

		#graph = Genegraph()
		#score, features = graph.getScore()

		te.score = 1 / float(te.location[1] - te.location[0])
		te.features = None

		return te.score