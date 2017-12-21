"""Expand target interval by source interval
	- if source[0] is before target interval, then just shift target interval
	- if source[0] is in target interval then, expand the target interval

Arguments:
	source(list): interval to use for expansion
	target(list): interval to expand

Returns:
	list: expanded target interval
"""
def expand(source, target):
	length = (source[1] - source[0])
	if source[0] <= target[0]:
		target[0] += length
	if source[0] <= target[1]:
		target[1] += length
	return target

"""Compare intervals a, b
	
Arguments:
	a(list)
	b(list)

Returns:
	int: -1 if a[1]<b[0], 0 if a[1]=b[0], 1 otherwise
"""
def compare(a, b):
	return -1 if a[1] < b[0] else (0 if a[1] == b[0] else 1)

#a contains b as subinterval
"""Checks if a contains b as subinterval

Arguments:
	a(list): interval to check
	b(list): query subinterval

Returs:
	bool: True if b is complete subinterval of a
"""
def contains(a, b):
	return a[0] <= b[0] and a[1] >= b[1]

#remove a from b
"""Remove a from b

Arguments:
	a(list): interval to remove
	b(list): interval to be removed from

Returs:
	list: b with a removed
"""
def remove(a, b):
	return [[b[0],a[0]], [a[1], b[1]]]

#get intervals of a cropped by b
"""Crop intervals B from a

Arguments:
	a(list): interval to be cropped
	B(list): list of intervals to remove from a
"""
def crop(a, B):
	if not len(B):
		return [a]

	#sort children
	B.sort(key=lambda x: x[0])

	cropped = [[a[0], B[0][0]]]
	for i in range(len(B) - 1):
		cropped.append([B[i][1], B[i+1][0]])
	cropped.append([B[-1][1], a[1]])	

	return cropped
