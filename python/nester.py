#!/usr/bin/env python3

from python import intervals
from python.gene import Gene


class Nester(object):
	"""Class represent nesting in sequence, recursivelly find best evaluated transposon and crop it until no new transposons found
	
	Attributes:
		seqid (str): id of sequence
		sequence (Bio.Seq.Seq): sequence
		nestedList (list[TE]): expanded list of nested transposons, the lower the index the sooner TE has been cropped
	"""
	def __init__(self, sequence):
		self.seqid = sequence.id
		self.sequence = sequence.seq
		self.nestedList = []
		self._findNestedTransposonList()

	def _findNestedTransposonList(self): #find list of nested transposons
		self.nestedList = self._getUnexpandedTransposonList(self.sequence)
		self._expandTransposonList()

	def _expandTransposonList(self): #backwards expanding of intervals according to previously found and cropped elements
		for i in reversed(range(len(self.nestedList) - 1)):
			for j in range(i + 1, len(self.nestedList)):
				self.nestedList[j].location = intervals.expand(self.nestedList[i].location, self.nestedList[j].location)
				for domain in self.nestedList[j].features['domains']:
					domain.location = intervals.expand(self.nestedList[i].location, domain.location)
				self.nestedList[j].features['ppt'] = intervals.expand(self.nestedList[i].location, self.nestedList[j].features['ppt'])
				self.nestedList[j].features['pbs'] = intervals.expand(self.nestedList[i].location, self.nestedList[j].features['pbs'])

	def _getUnexpandedTransposonList(self, sequence): #recursivelly find and crop best evaluated transposon, return unexpanded list of found transposons
		gene = Gene(self.seqid, sequence)
		bestCandidate = gene.getBestCandidate()
		
		if not bestCandidate:
			return []

		nested = [bestCandidate]

		#crop TE and call recursivelly
		croppedSequence = sequence[:(bestCandidate.location[0] - 1)] + sequence[(bestCandidate.location[1] + 1):]
		nested += self._getUnexpandedTransposonList(croppedSequence)

		return nested

		
		


