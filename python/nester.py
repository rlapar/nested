#!/usr/bin/env python3

from python import intervals
from python.gene import Gene


class Nester(object):
	def __init__(self, sequence):
		self.seqid = sequence.id
		self.sequence = sequence.seq
		self.nestedList = []
		self._findNestedTransposonList()

	def _findNestedTransposonList(self):
		self.nestedList = self._getUnexpandedTransposonList(self.sequence)
		self._expandTransposonList()

	def _expandTransposonList(self):
		for i in reversed(range(len(self.nestedList) - 1)):
			for j in range(i + 1, len(self.nestedList)):
				self.nestedList[j].location = intervals.expand(self.nestedList[i].location, self.nestedList[j].location)
				for domain in self.nestedList[j].features['domains']:
					domain.location = intervals.expand(self.nestedList[i].location, domain.location)
				self.nestedList[j].features['ppt'] = intervals.expand(self.nestedList[i].location, self.nestedList[j].features['ppt'])
				self.nestedList[j].features['pbs'] = intervals.expand(self.nestedList[i].location, self.nestedList[j].features['pbs'])

	def _getUnexpandedTransposonList(self, sequence):
		gene = Gene(self.seqid, sequence)
		bestCandidate = gene.getBestCandidate()
		
		if not bestCandidate:
			return []

		nested = [bestCandidate]

		#crop TE and call recursivelly
		croppedSequence = sequence[:(bestCandidate.location[0] - 1)] + sequence[(bestCandidate.location[1] + 1):]
		nested += self._getUnexpandedTransposonList(croppedSequence)

		return nested

		
		


