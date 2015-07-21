#!/usr/bin/python


fin1 = open('edge_kmerFreq.txt', 'r')
fin2 = open('predict.txt', 'r')

edge_freq = {}
count = 0
try:
	while 1:
		line = fin1.next().strip()
		count += 1
		(a, b, c) = line.split()
		a = int(a)
		b = int(b)
		c = int(c)
		
		if (a,b) in edge_freq:
			print 'repeat'
		else:
			edge_freq[(a, b)] = c
except:
	print 'read end'

print 'count' + str(count)
fin1.close()

count1 = 0
count2 = 0

repeat_to_uniq = 0
uniq_to_repeat = 0

try:
	while 1:
		line = fin2.next().strip()
			
#		print 'predict ....' + line

		(a, b, c) = line.split()
		a = int(a)
		b = int(b)
		c = int(round(float(c),0))
		if (a, b) in edge_freq:
			if edge_freq[(a, b)] == c:
				count1 += 1
			else:
				print str((a, b)) + '\t' + str(edge_freq[(a, b)]) + '\t' + str(c)
				if edge_freq[(a, b)] == 1:
					uniq_to_repeat += 1
				if c == 1:
					repeat_to_uniq += 1
				count2 += 1
				
#		else:
#			print 'help !!!!!!!!!!!!!!!!'

except:
	fin1.close()
	fin2.close()
	print str(float(count1)/(count2 + count1))
	print 'uniq_to_repeat = ' + str(uniq_to_repeat)
	print 'repeat_to_uniq = ' + str(repeat_to_uniq)


