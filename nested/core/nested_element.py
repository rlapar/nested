#!/usr/bin/env python3

class NestedElement(object):
    """Transfer class used for output.

    Attributes:
        id (str): id of sequence
        sequence (Bio.Seq.Seq): sequence
        nested_list (list[TE]): expanded list of nested transposons, the lower the index the sooner TE has been cropped
    """
    def __init__(self, id, sequence, nested_list):
        self.id = id
        self.sequence = sequence
        self.nested_list = nested_list

    def __str__(self):
        strlen = 15
        lines = ['{{id = {},'.format(self.id),
                 'sequence = {a}{b},'.format(a=self.sequence[:strlen], b='...' if len(self.sequence) > strlen else ''),
                 'nested_list = {}}}'.format(self.nested_list)]
        return '\n'.join(lines)