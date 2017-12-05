#!/usr/bin/env python3

from python import intervals
from python.gene import Gene


class Nester(object):
	def __init__(self, sequence):
		self._seqid = sequence.id
		self._sequence = sequence.seq
		self.nestedList = []
		self._findNestedTransposonList()

	def _findNestedTransposonList(self):
		self.nestedList = self._getUnexpandedTransposonList(self._sequence)
		self._expandTransposonList()

	def _expandTransposonList(self):
		for i in reversed(range(len(self.nestedList) - 1)):
			for j in range(i + 1, len(self.nestedList)):
				self.nestedList[j].location = intervals.expandInterval(self.nestedList[i].location, self.nestedList[j].location)
				#TODO expand domains
				#TODO expand ppt, pbs

	def _getUnexpandedTransposonList(self, sequence):
		gene = Gene(self._seqid, sequence)
		bestCandidate = gene.getBestCandidate()
		
		if not bestCandidate:
			return []

		nested = [bestCandidate]

		#crop TE and call recursivelly
		croppedSequence = sequence[:(bestCandidate.location[0] - 1)] + sequence[(bestCandidate.location[1] + 1):]
		nested += self._getUnexpandedTransposonList(croppedSequence)

		return nested

		
		


