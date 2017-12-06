def expand(source, target):
	length = (source[1] - source[0])
	if source[0] <= target[0]:
		target[0] += length
	if source[0] <= target[1]:
		target[1] += length
	return target

def compare(a, b):
	return -1 if a[1] < b[0] else (0 if a[1] == b[0] else 1)

#a contains b as subinterval
def contains(a, b):
	return a[0] <= b[0] and a[1] >= b[1]

#remove a from b
def remove(a, b):
	return [[b[0],a[0]], [a[1], b[1]]]

#get intervals of a cropped by b
def crop(a, b):
	if not len(b):
		return [a]

	#sort children
	b.sort(key=lambda x: x[0])

	cropped = [[a[0], b[0][0]]]
	for i in range(len(b) - 1):
		cropped.append([b[i][1], b[i+1][0]])
	cropped.append([b[-1][1], a[1]])	

	return cropped
