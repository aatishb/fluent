import numpy as np
import flu
import matplotlib.pyplot as plt

AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def count_AA(seq):

	dict = {}

	for word in AA:
		dict[word] = 0

	count = 0.0	
	for word in seq:
		if word in AA:		
			dict[word] += 1
			count += 1

	if not (count == 0):
		return [dict[word]/count for word in AA]

	else:
		return [dict[word] for word in AA]


def deltat(virarray):

	meandist = []

	print "Calculating mean deltaT"
	for year1 in years:
		for year2 in years:
			if years.index(year2) == years.index(year1) + 1:

				i = years.index(year1)
			
				dist = []

				for word1 in virarray[i]:
					for word2 in virarray[i+1]:
						label1 = str(word1.place)+'/'+str(word1.year)+'/'+str(word1.id)
						label2 = str(word2.place)+'/'+str(word2.year)+'/'+str(word2.id)
						dist.append(tree.distance(label1,label2))

				meandist.append(np.mean(dist))
	print "done"

	return meandist

def plotstuff(virarray, time):

	# Create array of frequencies at each site. freqarray : years, sites, freq of amino acids
	freqarray = np.array([[count_AA(word) for word in zip(*[vir.seq for vir in virpool])] for virpool in virarray], dtype=float)

	#for year in years:
	#	i = years.index(year)
	#	print year, len(virarray[i]), np.sum([flu.entropy(word) for word in freqarray[i]])

	print "\n"

	# Next, create and populate 3 arrays, for relative entropy, probability of fitting neutral model, and selection
	prob = [[] for year in years[0:-1]]
	relent = [[] for year in years[0:-1]]
	sel = [[] for year in years[0:-1]]

	for year1 in years:
		for year2 in years:
			if years.index(year2) == years.index(year1) + 1:
				i = years.index(year1)

				#calculate relative entropy pairs of sites across years
				xent = []
				for (site1,site2) in zip(freqarray[i+1],freqarray[i]):
					if not (np.sum(site1)*np.sum(site2) == 0):
						xent.append(flu.KLdivergence(site1,flu.expQ(site2,time[i])))
					else:
						xent.append(0.0)

				N = len(virarray[i+1])
				print N
				relent[i] = np.array(xent,dtype=float)
				prob[i] = np.exp(N*np.array(xent,dtype=float))
				sel[i] = np.less_equal(prob[i],0.001)
			
	print "\n"


	print "Generating figures"
	sites = range(588)

	plt.clf()
	plt.subplots_adjust(bottom=0)
	plt.imshow(sel,interpolation='nearest', cmap=plt.cm.Reds,aspect='auto')
	plt.yticks(range(len(years[0:-1])), years[0:-1])
	plt.savefig("heatmap.png",bbox_inches='tight',dpi=100)


	plt.clf()
	plt.subplots_adjust(bottom=-0.5)

	for i in range(len(years)-1):
		ax = plt.subplot(str(511+i))
		ax.set_yscale('log')
		ax.plot(sites, prob[i], label=years[i])
	plt.savefig("prob.png",bbox_inches='tight',dpi=100)


	plt.clf()
	plt.subplots_adjust(bottom=-0.5)

	for i in range(len(years)-1):
		ax = plt.subplot(str(511+i))
		ax.set_yscale('linear')
		ax.plot(sites, relent[i], label=years[i])
		ax.plot(sites, (np.log(0.01)/len(virarray[i+1]))*np.ones(588))
		ax.plot(sites, (np.log(0.001)/len(virarray[i+1]))*np.ones(588))
	plt.savefig("relent.png",bbox_inches='tight',dpi=100)

	print "Done"


from Bio import Phylo
tree = Phylo.read('egyptH5NHA.phy_phyml_tree.txt', 'newick')
nodes = tree.get_terminals()

clade1 = Phylo.read('clade1.txt', 'newick')
c1nodes = [word.name for word in clade1.get_terminals()]
c1ids = [word.name.split('/')[2] for word in clade1.get_terminals()]

clade2 = Phylo.read('clade2.txt', 'newick')
c2nodes = [word.name for word in clade2.get_terminals()]
c2ids = [word.name.split('/')[2] for word in clade2.get_terminals()]


years = ['2006','2007','2008','2009','2010','2011']
virarray  = [[word for word in flu.seqlist if (word.place == "Egypt" and word.year == year)] for year in years]


time = [0.0143905347475, 0.039710268946, 0.0457343136228, 0.0451385507375, 0.0425961707733]
timec1 = [0.013992550923,0.040050920329,0.0492249534222,0.0624247846509,0.0557367576217]
timec2 = [0.0144531601812,0.0208246072785,0.0251183491289,0.0228981542124,0.0199896902198]

#print deltat(virarray)
plotstuff(virarray, time)



