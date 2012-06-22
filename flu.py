from Bio import SeqIO
import numpy as np

class sequence(object):
    def __init__(self, name, place, year, source, seq):
        self.id = name
        self.place = place
        self.year = year
        self.source = source
	self.seq = seq

seqlist = []

for rec in SeqIO.parse('H5N HA.fst', 'fasta'):

	info = rec.description.split('/')

	if '' in info:
		info.remove('')

	if len(info) > 4:
		animal = info[1]	
		place = info[2]
		
		
		if ' ' in info[4]:
			year = info[4].split(' ')[1]
		elif ' ' in info[5]:
			year = info[5].split(' ')[1]
		else:
			year = '?'		

		seqlist.append(sequence(str(rec.id), place, year, animal, str(rec.seq)))		



def KLdivergence(vec1,vec2):

	xent = 0.0
	for p,q in zip(vec1,vec2):
		if not (p == 0):
			xent = xent + p*np.log(p/q)
	return xent



def entropy(vec):

	ent = 0.0
	for p in vec:
		if not (p == 0):
			ent = ent - p*np.log(p)
	return ent



def initM(freq):

	f = open('flumatrix.txt','r')
	rows = [line.split() for line in f]
	f.close()

	Q = np.array([[0.0 for i in range(20)] for j in range(20)])

	for row in rows:
		if not (row == []):
			i = rows.index(row)
			j = 0
			for word in row:
				Q[i][j] = float(word)*freq[j]
				Q[j][i] = float(word)*freq[i]
				j = j + 1

	for i in range(20):
		sum = np.sum(Q[i])
		Q[i][i] = -sum

	evals, U = np.linalg.eig(Q)
	Ut = np.linalg.inv(U)
	L = np.diag(evals)
	
	return U,evals,Ut



def expQ(freq,time):

	expLt = np.diag(np.exp(evals*time))
	M = np.dot(np.dot(U,expLt),Ut)
	return np.dot(M,freq)



g = np.array([float(num) for num in '0.0470718	0.0509102	0.0742143	0.0478596	0.0250216	0.0333036	0.0545874	0.0763734	0.0199642	0.0671336	0.0714981	0.0567845	0.0181507	0.0304961	0.0506561	0.0884091	0.0743386	0.0185237	0.0314741	0.0632292'.split()])
f=g/np.sum(g)
#f = np.array([ 0.05105 , 0.04219 , 0.08457 , 0.04646 , 0.02647 , 0.03385 , 0.07029 , 0.06871 , 0.02502 , 0.07154 , 0.08383 , 0.06506 , 0.02428 , 0.0329 , 0.03358 , 0.07621 , 0.04835 , 0.01764 , 0.04183 , 0.05618])

[U,evals,Ut] = initM(f)

