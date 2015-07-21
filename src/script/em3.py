#!/usr/bin/python
import math
import sys
import random

dataClass = []
mu1 = 0
delta1 = 10.0 #insert size variation
mu2 = 0
delta2 = 10.0
MAXITERATION = 10000
class1 = 1
class2 = 1
w1 = 0.0
w2 = 0.0
EPISILON = 1.0

def Gaussian1(x):
	global w1
	global delta1
	global mu1
	if delta1 == 0:
		return (math.log(w1))
	##print 'w1 = ' + str(w1) + '\t' + str(mu1) + '\t' + str(delta1) + '\t' + str(x)
	if  ((math.exp(-((((x - mu1)**2))/(2*(delta1**2)))))) == 0:
		return (math.log(w1)) + math.log( ((math.sqrt(2*((math.pi))))*delta1)**(-1) * ((math.exp(-((((mu1 + 5 * delta1 - mu1)**2))/(2*(delta1**2)))))))
	#return (math.log(w1)) + math.log( ((math.sqrt(2*((math.pi))))*delta1)**(-1) * ((math.exp(-((((x - mu1)**2))/(2*(delta1**2)))))))
	return (math.log(w1)) + math.log( ((math.sqrt(2*((math.pi))))*delta1)**(-1)) -  (((-((((x - mu1)**2))/(2*(delta1**2))))))
			

def Gaussian2(x):
	global w2
	global delta2
	global mu2
	if delta2 == 0:
		return (math.log(w2))
	##print 'w2 = ' + str(w2) + '\t' + str(mu2) + '\t' + str(delta2) + '\t' + str(x)
	if  ((math.exp(-((((x - mu2)**2))/(2*(delta2**2)))))) == 0:
		return (math.log(w2)) + math.log( ((math.sqrt(2*((math.pi))))*delta2)**(-1) * ((math.exp(-((((mu2 + 5*delta2 - mu2)**2))/(2*(delta2**2)))))))
	#return (math.log(w2)) +math.log( ((math.sqrt(2*((math.pi))))*delta2)**(-1) * ((math.exp(-((((x - mu2)**2))/(2*(delta2**2)))))))
	return (math.log(w2)) + math.log( ((math.sqrt(2*((math.pi))))*delta2)**(-1)) - ((math.exp(-((((x - mu2)**2))/(2*(delta2**2))))))


def count(data):
	global sum1
	global sum2
	global mu1
	global mu2
	global delta1
	global delta2
	global class1
	global class2
	global dataClass
	global w1
	global w2
	sum1 = 0.0
	sum2 = 0.0
	squareSum1 = 0.0
	squareSum2 = 0.0
	class1 = 1
	class2 = 1
	for i in range(len(data)):
		if dataClass[i] == 1:
			class1 += 1
			sum1 += data[i]
		if dataClass[i] == 2:
			class2 += 1
			sum2 += data[i]
	#print str(class1)
	#print str(class2)
	mu1 = (sum1*(1.0))/class1
	mu2 = (sum2*(1.0))/class2
	w1 = (class1*1.0)/(class1+class2)
	w2 = (class2*1.0)/(class1+class2)
	for i in range(len(data)):
		if dataClass[i] == 1:
			squareSum1 += (data[i]-mu1)**2
		if dataClass[i] == 2:
			squareSum2 += (data[i]-mu2)**2
			
def em3(data):
	global sum1
	global sum2
	global mu1
	global mu2
	global delta1
	global delta2
	global class1
	global class2
	global dataClass
	global w1
	global w2
	class1 = 1
	class2 = 1
	mu1 = 0
	mu2 = 0
	delta1 = 10.0 
	delta2 = 10.0
	MAXITERATION = 10000
	class1 = 1
	class2 = 1
	w1 = 0.0
	w2 = 0.0
	EPISILON = 1.0
	dataClass = []
	for i in range(len(data)):
		dataClass.append(random.randint(1,2))

	for i in range(len(data)):
		if dataClass[i] == 1:
			class1 += 1
		if dataClass[i] == 2:
			class2 += 1

	expectation = [() for i in range(len(data))]
	for i in range(MAXITERATION):
		c1 = 0
		c2 = 0
		for i in range(len(data)):
			count(data)
			x = data[i]
			below = 0
			below += Gaussian1(x) + Gaussian2(x)
			upper1 = Gaussian1(x)
			upper2 = Gaussian2(x)
			expect1 = upper1/below
			expect2 = upper2/below
			pair = (expect1,expect2)
			if abs(expect1) <= abs(expect2) and dataClass[i] == 2:
				dataClass[i] = 1
				c1 += 1
			if abs(expect1) > abs(expect2) and dataClass[i] == 1:
				c2 += 1
				dataClass[i] = 2
		if c1 <= 1 and c2 <= 1:
			break

	sum1 = 0
	sum2 = 0
	class1 = 0
	class2 = 0
	squareSum1 = 0.0
	squareSum2 = 0.0

	#out = open('pointClass1','w')
	for i in range(len(data)):
		if dataClass[i] == 1:
			sum1 += data[i]
			class1 += 1
			#out.write(str(data[i]) + '\t' + str(dataClass[i]) + '\n')
	#out.close()

	#out = open('pointClass2','a')
	for i in range(len(data)):
		if dataClass[i] == 2:
			sum2 += data[i]
			class2 += 1
			#out.write(str(data[i]) + '\t' + str(dataClass[i]) + '\n')
	#out.close()

	for i in range(len(data)):
		if dataClass[i] == 1:
			squareSum1 += (data[i]-mu1)**2
		if dataClass[i] == 2:
			squareSum2 += (data[i]-mu2)**2
	if class1 > 0:
		delta1 = math.sqrt(squareSum1/class1)
		mu1 = int(sum1*1.0/class1)
		#print 'class 1 mu = ' + str(sum1*1.0/class1)
		#print 'class 1 delta1 = ' + str(delta1)
	if class2 > 0:
		delta2 = math.sqrt(squareSum2/class2)
		mu2 = int(sum2*1.0/class2)
		#print 'class 2 mu = ' + str(sum2*1.0/class2)
		#print 'class 2 delta2 = ' + str(delta2)
