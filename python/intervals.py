def expandInterval(source, target):
	length = (source[1] - source[0])
	if source[0] <= target[0]:
		target[0] += length
	if source[0] <= target[1]:
		target[1] += length
	return target