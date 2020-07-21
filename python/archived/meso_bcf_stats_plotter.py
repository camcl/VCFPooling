import csv
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
# sns.set()
import json
from General import BCFstatsUtils as butils
from General.Style import ind_palette
from General.Style import geno_palette,type_palette
from General.Style import fwidth_small
from General.Style import fwidth_med
from General.Style import  lw, if_title, if_overlap_legend,if_ind_legend, if_geno_legend, kept_alpha, dropped_alpha,legfontsize, lw2,retained_col, fheight_med,fheight_small, lw3, fheight_med_smaller, lw4,dashstyle
import matplotlib.ticker as mticker
from General.Style import fdpi
from matplotlib import colors
import General.Style
import General.Style
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap


def plot_nrds_by_subset(dir):

	with open(dir + 'all_tables.json', 'r') as fp:
		all_tables = json.load(fp)

	# inds = ["ans17", "car16", "LBK", "Loschbour", "ne1", "ne1CP"]
	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]

	print(inds)

	# filters = ["all", "filtered", "thets"]
	filters = ["filtered"]

	subsets = ["kept", "dropped"]

	fig, ax = plt.subplots(figsize=(8,9))
	plt.suptitle("Genotype discordance by genotype in HC\n (% of genotypes wrong in imputed)\n ________________________________________\nkept                          |                          dropped")
	counter = 0
	for filter in filters:
		for subset in subsets:
			counter += 1
			all_nrds = []

			# plt.subplots()
			for ind in inds:

				nrd_dict = all_tables[ind][filter][subset][0]["# NRDs"]
				keys = ["Ref/Ref discordance", "Ref/Alt discordance", "Alt/Alt discordance","NRD"]
				this_nrds = []
				for k in keys:
					this_nrds.append(nrd_dict[k])

				all_nrds.append(this_nrds)
			print(all_nrds)


			ax = plt.subplot(3,2,counter)

			x = [0,10,20,30]


			plt.bar(map(lambda y : y-3, x),all_nrds[0], width = 1, label=inds[0])
			plt.bar(map(lambda y : y-2, x),all_nrds[1], width = 1, label=inds[1])
			plt.bar(map(lambda y : y-1, x),all_nrds[2], width = 1, label=inds[2])
			plt.bar(map(lambda y : y-0, x),all_nrds[3], width = 1, label=inds[3])
			plt.bar(map(lambda y : y+1, x),all_nrds[4], width = 1, label=inds[4])
			# plt.bar(map(lambda y : y+2, x),all_nrds[5], width = 1, label=inds[5])

			if counter == 1:
				# plt.title(subset)
				plt.ylabel("all markers")

				plt.legend(fontsize=7.5)
			if counter == 3:
			# 	plt.title(subset)
				plt.ylabel("quality filtered")

			if counter == 5:
			# 	plt.title(subset)
				plt.ylabel("trusted heterozygotes")

			# xlabels = ["R/R discordance", "R/A discordance", "A/A discordance","Non-Ref Discordance"]
			xlabels = ["Hom-Ref", "Het", "Hom-Alt","Non-Ref"]

			plt.xticks(x,xlabels)
			plt.ylim(ymax=20)
			# plt.xticks(rotation=-15)



	plt.legend()
	plt.savefig(dir+"Genotype_discordance_by_subset.pdf")
	#plt.show()


def plot_nrds_ind(dir):

	with open(dir + 'all_tables.json', 'r') as fp:
		all_tables = json.load(fp)


	inds = ["ans17", "car16", "LBK", "Loschbour", "ne1", "ne1CP"]
	# inds = ["ans17", "car16", "LBK", "Loschbour", "ne1"]

	print(inds)

	filters = ["all", "filtered", "thets"]
	subsets = ["kept", "dropped"]


	for ind in inds:

		fig, ax = plt.subplots(figsize=(4,9))
		plt.tight_layout()
		plt.suptitle("Genotype discordance "+ind)
		plt.subplots_adjust(top=0.93)
		counter = 0

		for filter in filters:
			counter += 1
			plt.subplot(3,1,counter)

			all_nrds = []
			for subset in subsets:

				nrd_dict = all_tables[ind][filter][subset][0]["# NRDs"]
				keys = ["Ref/Ref discordance", "Ref/Alt discordance", "Alt/Alt discordance","NRD"]
				this_nrds = []
				for k in keys:
					this_nrds.append(nrd_dict[k])

				all_nrds.append(this_nrds)
			print(all_nrds)


			x = [0,10,20,30]


			plt.bar(map(lambda y : y-1, x),all_nrds[0], width = 1, label=subsets[0], color="black")
			plt.bar(map(lambda y : y-0, x),all_nrds[1], width = 1, label=subsets[1], color="grey")

			if counter == 1:
				# plt.title(subset)
				plt.title("All Markers")

				plt.legend(fontsize=7.5)
			if counter == 2:
			# 	plt.title(subset)
				plt.title("Quality Filtered")

			if counter == 3:
			# 	plt.title(subset)
				plt.title("Trusted Heterozygotes")

			# xlabels = ["R/R discordance", "R/A discordance", "A/A discordance","Non-Ref Discordance"]
			xlabels = ["Hom-Ref", "Het", "Hom-Alt","Non-Ref"]

			plt.xticks(x,xlabels)
			plt.ylim(ymax=11)
			# plt.xticks(rotation=-15)


		# plt.legend()
		plt.savefig(dir+"Genotype_discordance_"+ind+".pdf")
		# plt.show()


#Plot discordance and what types of errors there are in each category for files in dir
def  plot_nrds_gcts_stacked(dir, cov, refset):

	with open(dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json', 'r') as fp:
		all_tables = json.load(fp)



	lims = [7,7,11,10,12,20]
	# lims = [2,2,2,2,2,2]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]


	subsets = ["kept", "dropped"]

	for ind in inds:

		fig, ax = plt.subplots(figsize=(5,4))
		# plt.tight_layout()
		plt.suptitle("Genotype discordance by genotype in HC - "+ind)
		plt.subplots_adjust(top=0.93)
		counter = 0

		counter += 1
		ax = plt.subplot(1,1,counter)

		# for this individual and filter : kept and dropped gcts
		all_gcts_bottom = []
		all_gcts_top = []

		for subset in subsets:

			gct_dict = all_tables[ind][subset][1]["# GCTs"]
			keys = ["Ref/Ref discordance", "Ref/Alt discordance", "Alt/Alt discordance","NRD"]
			keys = [u'RR Hom -> RA Het',  u'RR Hom -> AA Hom',u'RA Het -> RR Hom', u'RA Het -> AA Hom',  u'AA Hom -> RR Hom',  u'AA Hom -> RA Het']


			this_gcts_bottom = []
			this_gcts_top = []

			# RR
			num_RR_RA = gct_dict['RR Hom -> RA Het']
			num_RR_AA = gct_dict['RR Hom -> AA Hom']
			num_RR_RR = gct_dict['RR Hom -> RR Hom']

			print("RR-RA: " + str(num_RR_RA))
			print("RR-AA: " + str(num_RR_AA))

			RR_wrongs= int(num_RR_RA) + int(num_RR_AA)
			RR_all= int(num_RR_RA) + int(num_RR_AA) + int(num_RR_RR)


			if RR_wrongs > 0:
				perc_RR_RA =  float(num_RR_RA) / RR_wrongs
				perc_RR_AA =  float(num_RR_AA) / RR_wrongs
				RR_disc = 100* float(RR_wrongs) / RR_all

			else:
				perc_RR_RA = 0.0
				perc_RR_AA = 0.0
				RR_disc = 0.0

			perc_RR_RA = perc_RR_RA * RR_disc
			perc_RR_AA = perc_RR_AA * RR_disc


			print(perc_RR_RA)
			print(perc_RR_AA)

			# RA
			num_RA_RR = gct_dict['RA Het -> RR Hom']
			num_RA_AA = gct_dict['RA Het -> AA Hom']
			num_RA_RA = gct_dict['RA Het -> RA Het']


			RA_wrongs= int(num_RA_RR) + int(num_RA_AA)
			RA_all= int(num_RA_RR) + int(num_RA_AA) + int(num_RA_RA)
			RA_disc = 100* float(RA_wrongs) / RA_all


			perc_RA_RR =  float(num_RA_RR) / RA_wrongs
			perc_RA_AA =  float(num_RA_AA) / RA_wrongs


			perc_RA_RR = perc_RA_RR * RA_disc
			perc_RA_AA = perc_RA_AA * RA_disc


			# AA
			num_AA_RR = gct_dict['AA Hom -> RR Hom']
			num_AA_RA = gct_dict['AA Hom -> RA Het']
			num_AA_AA = gct_dict['AA Hom -> AA Hom']

			AA_wrongs= int(num_AA_RR) + int(num_AA_RA)
			AA_all= int(num_AA_RR) + int(num_AA_RA) + int(num_AA_AA)

			if AA_wrongs > 0:
				perc_AA_RR =  float(num_AA_RR) / AA_wrongs
				perc_AA_RA =  float(num_AA_RA) / AA_wrongs
				AA_disc = 100* float(AA_wrongs) / AA_all

			else:
				perc_AA_RR = 0
				perc_AA_RA = 0
				AA_disc = 0.0

			perc_AA_RR = perc_AA_RR * AA_disc
			perc_AA_RA = perc_AA_RA * AA_disc


			this_gcts_bottom = [perc_RR_RA, perc_RA_RR,perc_AA_RR]
			this_gcts_top = [perc_RR_AA, perc_RA_AA,perc_AA_RA]



			nrd_dict = all_tables[ind][subset][0]["# NRDs"]
			keys = ["Ref/Ref discordance", "Ref/Alt discordance", "Alt/Alt discordance","NRD"]
			this_nrds = []
			for k in keys:
				this_nrds.append(nrd_dict[k])

			this_gcts_bottom.append(nrd_dict["NRD"])

			all_gcts_bottom.append(this_gcts_bottom)
			all_gcts_top.append(this_gcts_top)

			RRcol = "steelblue"
			RAcol = "darkmagenta"
			AAcol = "green"

			bottom_texts = ["RA", "RR", "RR"]
			bottom_colors = [RAcol, RRcol, RRcol, "black"]

			top_texts = ["AA", "AA", "RA"]
			top_colors = [AAcol, AAcol, RAcol, "black"]

			print(all_gcts_bottom)
			print(all_gcts_top)

			x_bottom = [0,10,20,30]
			x_top = [0,10,20]

			# x = range(len(keys))
			# x = map(lambda y : 10*y, x)

			print(x_bottom)
			width=3.0

			# ax.text(1, 2, "agg", ha='center',va='center')

			# for

			# plt.text(x_top, all_gcts_bottom[0][:-1], bottom_texts)
			# ax.text(x_top, map(lambda y: y+5 , all_gcts_top[0][:-1]), top_texts)


		#kept
		#label=subsets[0]
		plt.bar(map(lambda y : y-(width/2), x_bottom),all_gcts_bottom[0], width = width, color=bottom_colors)
		plt.bar(map(lambda y : y-(width/2), x_top),all_gcts_top[0], bottom=all_gcts_bottom[0][:-1], width = width, color=top_colors)

		#dropped
		plt.bar(map(lambda y : y+(width/2), x_bottom),all_gcts_bottom[1], width = width, color=bottom_colors, alpha=0.5)
		plt.bar(map(lambda y : y+(width/2), x_top),all_gcts_top[1], bottom=all_gcts_bottom[1][:-1], width = width, color=top_colors, alpha=0.5)

		adj=1.5
		locs_kept = map(lambda y : y-(width/2), x_top)
		for i in range(len(all_gcts_top[0])):
			plt.text(locs_kept[i]-adj, all_gcts_bottom[0][i]+all_gcts_top[0][i]+0.05,'{0:.2f}'.format(all_gcts_bottom[0][i]+all_gcts_top[0][i]),size=7)

		locs_dropped = map(lambda y : y+(width/2), x_top)
		for i in range(len(all_gcts_top[1])):
			plt.text(locs_dropped[i]-adj, all_gcts_bottom[1][i]+all_gcts_top[1][i]+0.05,'{0:.2f}'.format(all_gcts_bottom[1][i]+all_gcts_top[1][i]),size=7)




		if counter == 1:
			# plt.title(subset)
			# plt.title("All Markers")
			plt.bar([0],[0], label="RR", color=RRcol)
			plt.bar([0],[0], label="RA", color=RAcol)
			plt.bar([0],[0], label="AA", color=AAcol)
			plt.legend(fontsize=8)
		if counter == 2:
			plt.title("Quality Filtered")

		if counter == 3:
			plt.title("Trusted Heterozygotes")

		xlabels = ["Hom-Ref", "Het", "Hom-Alt","Non-Ref"]

		plt.ylabel("% wrong genotype in imputed")
		plt.xticks(x_bottom,xlabels)
		plt.ylim(ymax=lims[inds.index(ind)])
		# plt.xticks(rotation=-15)


		# plt.legend()
		plt.savefig(dir+'_'+refset+'/'+"Genotype_discordance_stacked"+str(cov)+"_"+refset+"_"+ind+".pdf")
		# plt.show()

#Plot discordance in each true geno category for files in dir IN SAME FILE
def  plot_gcts_allinds(dir, cov, refset):

	with open(dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json', 'r') as fp:
		all_tables = json.load(fp)



	lims = [7,7,11,10,12,20]
	# lims = [2,2,2,2,2,2]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]
	ind_offsets = [-2, -1, 0, 1, 2]


	subsets = ["kept", "dropped"]

	fig, ax = plt.subplots(figsize=(fwidth_med, fheight_med))
	# plt.tight_layout()
	# plt.suptitle("Genotype discordance by genotype in HC")
	# plt.subplots_adjust(top=0.93)
	# ax = plt.subplot(1,1,1)


	for ind in inds:

		counter = 0

		counter += 1

		# for this individual and filter : kept and dropped gcts
		# all_gcts_bottom = []
		# all_gcts_top = []
		all_gcts = []

		for subset in subsets:

			gct_dict = all_tables[ind][subset][1]["# GCTs"]

			# RR
			num_RR_RA = gct_dict['RR Hom -> RA Het']
			num_RR_AA = gct_dict['RR Hom -> AA Hom']
			num_RR_RR = gct_dict['RR Hom -> RR Hom']

			print("RR-RA: " + str(num_RR_RA))
			print("RR-AA: " + str(num_RR_AA))

			RR_wrongs= int(num_RR_RA) + int(num_RR_AA)
			RR_all= int(num_RR_RA) + int(num_RR_AA) + int(num_RR_RR)


			if RR_wrongs > 0:
				perc_RR_RA =  float(num_RR_RA) / RR_wrongs
				perc_RR_AA =  float(num_RR_AA) / RR_wrongs
				RR_disc = 100* float(RR_wrongs) / RR_all

			else:
				perc_RR_RA = 0.0
				perc_RR_AA = 0.0
				RR_disc = 0.0

			perc_RR_RA = perc_RR_RA * RR_disc
			perc_RR_AA = perc_RR_AA * RR_disc


			print(perc_RR_RA)
			print(perc_RR_AA)

			# RA
			num_RA_RR = gct_dict['RA Het -> RR Hom']
			num_RA_AA = gct_dict['RA Het -> AA Hom']
			num_RA_RA = gct_dict['RA Het -> RA Het']


			RA_wrongs= int(num_RA_RR) + int(num_RA_AA)
			RA_all= int(num_RA_RR) + int(num_RA_AA) + int(num_RA_RA)
			RA_disc = 100* float(RA_wrongs) / RA_all


			perc_RA_RR =  float(num_RA_RR) / RA_wrongs
			perc_RA_AA =  float(num_RA_AA) / RA_wrongs


			perc_RA_RR = perc_RA_RR * RA_disc
			perc_RA_AA = perc_RA_AA * RA_disc


			# AA
			num_AA_RR = gct_dict['AA Hom -> RR Hom']
			num_AA_RA = gct_dict['AA Hom -> RA Het']
			num_AA_AA = gct_dict['AA Hom -> AA Hom']

			AA_wrongs= int(num_AA_RR) + int(num_AA_RA)
			AA_all= int(num_AA_RR) + int(num_AA_RA) + int(num_AA_AA)

			if AA_wrongs > 0:
				perc_AA_RR =  float(num_AA_RR) / AA_wrongs
				perc_AA_RA =  float(num_AA_RA) / AA_wrongs
				AA_disc = 100* float(AA_wrongs) / AA_all

			else:
				perc_AA_RR = 0
				perc_AA_RA = 0
				AA_disc = 0.0

			perc_AA_RR = perc_AA_RR * AA_disc
			perc_AA_RA = perc_AA_RA * AA_disc


			# this_gcts_bottom = [perc_RR_RA, perc_RA_RR,perc_AA_RR]
			# this_gcts_top = [perc_RR_AA, perc_RA_AA,perc_AA_RA]

			this_gcts = [perc_RR_RA + perc_RR_AA, perc_RA_RR +perc_RA_AA , perc_AA_RR + perc_AA_RA]



			nrd_dict = all_tables[ind][subset][0]["# NRDs"]
			keys = ["Ref/Ref discordance", "Ref/Alt discordance", "Alt/Alt discordance","NRD"]
			this_nrds = []
			for k in keys:
				this_nrds.append(nrd_dict[k])

			# this_gcts_bottom.append(nrd_dict["NRD"])

			all_gcts.append(this_gcts)



			# print(all_gcts_bottom)
			# print(all_gcts_top)

			x_bottom = [0,15,30]
			x_top = [0,15,30]

			# x = range(len(keys))
			# x = map(lambda y : 10*y, x)

			print(x_bottom)
			width=1.0

			# ax.text(1, 2, "agg", ha='center',va='center')

			# for

			# plt.text(x_top, all_gcts_bottom[0][:-1], bottom_texts)
			# ax.text(x_top, map(lambda y: y+5 , all_gcts_top[0][:-1]), top_texts)



		RRcol = "steelblue"
		RAcol = "darkmagenta"
		AAcol = "green"

		RRh = ""
		RAh = "/"
		AAh = "."


		bottom_texts = ["RA", "RR", "RR"]
		bottom_colors = [RAcol, RRcol, RRcol, "black"]
		bottom_hatches = [RAh, RRh, RRh]

		top_texts = ["AA", "AA", "RA"]
		top_colors = [AAcol, AAcol, RAcol, "black"]
		top_hatches = [AAh, AAh, RAh]

				#kept
		label=subsets[0]
		bsb = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 - (width/2), x_bottom),all_gcts[0], width = width, color=ind_palette[inds.index(ind)], alpha=0.8)

		# for b,h in zip(bsb, bottom_hatches):
			# b.set_edgecolor(ind_palette[inds.index(ind)])
			# b.set_alpha(0.1)
			# b.set_hatch(h)

		# for b,h in zip(bst, top_hatches):
			# b.set_edgecolor(ind_palette[inds.index(ind)])
			# b.set_edgecolor(ind_palette[inds.index(ind)])
			# b.set_hatch(h)

		#dropped
		bsb = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 + (width/2), x_bottom),all_gcts[1], width = width, color=ind_palette[inds.index(ind)], alpha=0.5)

		#
		# for b,h in zip(bsb, bottom_hatches):
		# 	b.set_edgecolor(ind_palette[inds.index(ind)])
		# 	b.set_hatch("/")
		#
		# for b,h in zip(bst, top_hatches):
		# 	b.set_edgecolor(ind_palette[inds.index(ind)])
		# 	b.set_hatch("/")


		# adj=1.5
		# locs_kept = map(lambda y : y-(width/2), x_top)
		# for i in range(len(all_gcts_top[0])):
		# 	plt.text(locs_kept[i]-adj, all_gcts_bottom[0][i]+all_gcts_top[0][i]+0.05,'{0:.2f}'.format(all_gcts_bottom[0][i]+all_gcts_top[0][i]),size=7)
		#
		# locs_dropped = map(lambda y : y+(width/2), x_top)
		# for i in range(len(all_gcts_top[1])):
		# 	plt.text(locs_dropped[i]-adj, all_gcts_bottom[1][i]+all_gcts_top[1][i]+0.05,'{0:.2f}'.format(all_gcts_bottom[1][i]+all_gcts_top[1][i]),size=7)
		#



		# if counter == 1:
		# 	# plt.title(subset)
		# 	# plt.title("All Markers")
		# 	plt.bar([0],[0], label="RR", color=RRcol)
		# 	plt.bar([0],[0], label="RA", color=RAcol)
		# 	plt.bar([0],[0], label="AA", color=AAcol)
		# 	plt.legend(fontsize=8)
		# if counter == 2:
		# 	plt.title("Quality Filtered")
		#
		# if counter == 3:
		# 	plt.title("Trusted Heterozygotes")

	xlabels = ["Hom-Ref", "Het", "Hom-Alt","Non-Ref"]

	plt.ylabel("Genotype discordance")
	plt.xlabel("True genotype")

	plt.xticks(x_bottom,xlabels)
	plt.ylim(ymax=1)

	# plt.xticks(rotation=-15)


	# plt.legend()
	print(dir+'_'+refset+'/'+"Genotype_discordance_allinds"+str(cov)+"_"+refset+".pdf")
	plt.savefig("Genotype_discordance_allinds"+str(cov)+"_"+refset+".svg", dpi=fdpi, bbox_inches="tight")
	plt.show()
	plt.close()

#Plot discordance in each true geno category for files in dir IN SAME FILE
def  plot_gcts_allinds_by_overlap(dir, cov, refset):

	with open(dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json', 'r') as fp:
		all_tables = json.load(fp)



	lims = [7,7,11,10,12,20]
	# lims = [2,2,2,2,2,2]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]
	ind_offsets = [-2, -1, 0, 1, 2]


	subsets = ["kept", "dropped"]

	fig, ax = plt.subplots(figsize=(fwidth_med, fheight_med))
	# plt.tight_layout()
	# plt.suptitle("Genotype discordance by genotype in HC")
	# plt.subplots_adjust(top=0.93)
	# ax = plt.subplot(1,1,1)


	for ind in inds:

		counter = 0

		counter += 1

		# for this individual and filter : kept and dropped gcts
		# all_gcts_bottom = []
		# all_gcts_top = []
		all_gcts = []

		for subset in subsets:

			gct_dict = all_tables[ind][subset][1]["# GCTs"]

			# RR
			num_RR_RA = gct_dict['RR Hom -> RA Het']
			num_RR_AA = gct_dict['RR Hom -> AA Hom']
			num_RR_RR = gct_dict['RR Hom -> RR Hom']

			print("RR-RA: " + str(num_RR_RA))
			print("RR-AA: " + str(num_RR_AA))

			RR_wrongs= int(num_RR_RA) + int(num_RR_AA)
			RR_all= int(num_RR_RA) + int(num_RR_AA) + int(num_RR_RR)


			if RR_wrongs > 0:
				perc_RR_RA =  float(num_RR_RA) / RR_wrongs
				perc_RR_AA =  float(num_RR_AA) / RR_wrongs
				RR_disc = 100* float(RR_wrongs) / RR_all

			else:
				perc_RR_RA = 0.0
				perc_RR_AA = 0.0
				RR_disc = 0.0

			perc_RR_RA = perc_RR_RA * RR_disc
			perc_RR_AA = perc_RR_AA * RR_disc


			print(perc_RR_RA)
			print(perc_RR_AA)

			# RA
			num_RA_RR = gct_dict['RA Het -> RR Hom']
			num_RA_AA = gct_dict['RA Het -> AA Hom']
			num_RA_RA = gct_dict['RA Het -> RA Het']


			RA_wrongs= int(num_RA_RR) + int(num_RA_AA)
			RA_all= int(num_RA_RR) + int(num_RA_AA) + int(num_RA_RA)
			RA_disc = 100* float(RA_wrongs) / RA_all


			perc_RA_RR =  float(num_RA_RR) / RA_wrongs
			perc_RA_AA =  float(num_RA_AA) / RA_wrongs


			perc_RA_RR = perc_RA_RR * RA_disc
			perc_RA_AA = perc_RA_AA * RA_disc


			# AA
			num_AA_RR = gct_dict['AA Hom -> RR Hom']
			num_AA_RA = gct_dict['AA Hom -> RA Het']
			num_AA_AA = gct_dict['AA Hom -> AA Hom']

			AA_wrongs= int(num_AA_RR) + int(num_AA_RA)
			AA_all= int(num_AA_RR) + int(num_AA_RA) + int(num_AA_AA)

			if AA_wrongs > 0:
				perc_AA_RR =  float(num_AA_RR) / AA_wrongs
				perc_AA_RA =  float(num_AA_RA) / AA_wrongs
				AA_disc = 100* float(AA_wrongs) / AA_all

			else:
				perc_AA_RR = 0
				perc_AA_RA = 0
				AA_disc = 0.0

			perc_AA_RR = perc_AA_RR * AA_disc
			perc_AA_RA = perc_AA_RA * AA_disc


			# this_gcts_bottom = [perc_RR_RA, perc_RA_RR,perc_AA_RR]
			# this_gcts_top = [perc_RR_AA, perc_RA_AA,perc_AA_RA]

			this_gcts = [perc_RR_RA + perc_RR_AA, perc_RA_RR +perc_RA_AA , perc_AA_RR + perc_AA_RA]



			nrd_dict = all_tables[ind][subset][0]["# NRDs"]
			keys = ["Ref/Ref discordance", "Ref/Alt discordance", "Alt/Alt discordance","NRD"]
			this_nrds = []
			for k in keys:
				this_nrds.append(nrd_dict[k])

			# this_gcts_bottom.append(nrd_dict["NRD"])

			all_gcts.append(this_gcts)



			# print(all_gcts_bottom)
			# print(all_gcts_top)

			x_bottom = [0,15,30]
			x_top = [0,15,30]

			# x = range(len(keys))
			# x = map(lambda y : 10*y, x)

			print(x_bottom)
			width=1.0

			# ax.text(1, 2, "agg", ha='center',va='center')

			# for

			# plt.text(x_top, all_gcts_bottom[0][:-1], bottom_texts)
			# ax.text(x_top, map(lambda y: y+5 , all_gcts_top[0][:-1]), top_texts)



		RRcol = "steelblue"
		RAcol = "darkmagenta"
		AAcol = "green"

		RRh = ""
		RAh = "/"
		AAh = "."


		bottom_texts = ["RA", "RR", "RR"]
		bottom_colors = [RAcol, RRcol, RRcol, "black"]
		bottom_hatches = [RAh, RRh, RRh]

		top_texts = ["AA", "AA", "RA"]
		top_colors = [AAcol, AAcol, RAcol, "black"]
		top_hatches = [AAh, AAh, RAh]

				#kept
		label=subsets[0]
		bsb = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)] - (width/2), x_bottom),all_gcts[0], width = width, color=ind_palette[inds.index(ind)], alpha=0.8, label=ind)

		# for b,h in zip(bsb, bottom_hatches):
			# b.set_edgecolor(ind_palette[inds.index(ind)])
			# b.set_alpha(0.1)
			# b.set_hatch(h)

		# for b,h in zip(bst, top_hatches):
			# b.set_edgecolor(ind_palette[inds.index(ind)])
			# b.set_edgecolor(ind_palette[inds.index(ind)])
			# b.set_hatch(h)

		#dropped
		bsb = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]+5 + (width/2), x_bottom),all_gcts[1], width = width, color=ind_palette[inds.index(ind)], alpha=0.5, label=ind)

		#
		# for b,h in zip(bsb, bottom_hatches):
		# 	b.set_edgecolor(ind_palette[inds.index(ind)])
		# 	b.set_hatch("/")
		#
		# for b,h in zip(bst, top_hatches):
		# 	b.set_edgecolor(ind_palette[inds.index(ind)])
		# 	b.set_hatch("/")


		# adj=1.5
		# locs_kept = map(lambda y : y-(width/2), x_top)
		# for i in range(len(all_gcts_top[0])):
		# 	plt.text(locs_kept[i]-adj, all_gcts_bottom[0][i]+all_gcts_top[0][i]+0.05,'{0:.2f}'.format(all_gcts_bottom[0][i]+all_gcts_top[0][i]),size=7)
		#
		# locs_dropped = map(lambda y : y+(width/2), x_top)
		# for i in range(len(all_gcts_top[1])):
		# 	plt.text(locs_dropped[i]-adj, all_gcts_bottom[1][i]+all_gcts_top[1][i]+0.05,'{0:.2f}'.format(all_gcts_bottom[1][i]+all_gcts_top[1][i]),size=7)
		#



		# if counter == 1:
		# 	# plt.title(subset)
		# 	# plt.title("All Markers")
		# 	plt.bar([0],[0], label="RR", color=RRcol)
		# 	plt.bar([0],[0], label="RA", color=RAcol)
		# 	plt.bar([0],[0], label="AA", color=AAcol)
		# 	plt.legend(fontsize=8)
		# if counter == 2:
		# 	plt.title("Quality Filtered")
		#
		# if counter == 3:
		# 	plt.title("Trusted Heterozygotes")

	xlabels = ["Hom-Ref", "Het", "Hom-Alt","Non-Ref"]

	plt.ylabel("Genotype discordance")
	plt.xlabel("True genotype")

	plt.xticks([i+2.5 for i in x_bottom],xlabels)
	plt.ylim(ymax=1)

	# plt.xticks(rotation=-15)


	# plt.legend()
	# print(dir+'_'+refset+'/'+"Genotype_discordance_allinds"+str(cov)+"_"+refset+".pdf")
	plt.savefig("Genotype_discordance_allinds_byoverlap"+str(cov)+"_"+refset+".pdf", dpi=fdpi, bbox_inches="tight")

	plt.savefig("Genotype_discordance_allinds_byoverlap"+str(cov)+"_"+refset+".svg", dpi=fdpi, bbox_inches="tight")
	plt.show()
	plt.close()



#Plot bar graph of concordance in each true geno category for files in dir IN SAME FILE
def  plot_concordance_per_geno_by_ind(dir, cov, refset, ref_title, ymin=0.8, resdir="", grid=False, wid=6, heig=4):


	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]
	ind_offsets = [-2, -1, 0, 1, 2]


	subsets = ["kept", "dropped"]

	fig, ax = plt.subplots(figsize=(wid, heig))

	if grid:
		plt.grid(b=True, which='major',axis='y', linewidth=lw3)

	for ind in inds:

		counter = 0

		counter += 1

		all_gcts = []

		for subset in subsets:

			concordance_RR = get_concordance_from_gcss(dir,refset,cov,ind,subset,['RR Hom'])
			concordance_RA = get_concordance_from_gcss(dir,refset,cov,ind,subset,['RA Het'])
			concordance_AA = get_concordance_from_gcss(dir,refset,cov,ind,subset,['AA Hom'])


			this_gcts = [concordance_RR, concordance_RA, concordance_AA]

			all_gcts.append(this_gcts)

			x_bottom = [0,15,30]

			width=1.0


		#kept
		label=subsets[0]
		bsb = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)] - (width/2), x_bottom),all_gcts[0], width = width, color=ind_palette[inds.index(ind)], alpha=kept_alpha, label=ind)
		# bsb = sns.barplot(map(lambda y : y + ind_offsets[inds.index(ind)] - (width/2), x_bottom),all_gcts[0], color=ind_palette[inds.index(ind)], alpha=kept_alpha, label=ind)

		#dropped
		bsb = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]+5 + (width/2), x_bottom),all_gcts[1], width = width, color=ind_palette[inds.index(ind)], alpha=dropped_alpha, label=ind)
		# bsb = sns.barplot(map(lambda y : y + ind_offsets[inds.index(ind)]+5 + (width/2), x_bottom),all_gcts[1],color=ind_palette[inds.index(ind)], alpha=dropped_alpha, label=ind)

	xlabels = ["RR", "RA", "AA"]

	plt.ylabel("Genotype concordance")
	plt.xlabel("Genotype in high-quality sample")

	plt.xticks([i+2.5 for i in x_bottom],xlabels)
	plt.ylim(ymax=1.0, ymin=ymin)

	plt.setp(ax.spines.values(), linewidth=0)

	if if_title:
		plt.title("Coverage: "+ str(cov) +" reference: "+ ref_title)



	plt.tight_layout(pad=0)
	plt.savefig(resdir+"conc_per_geno_by_ind_"+str(cov)+"_"+ref_title+".pdf", dpi=fdpi)


	# plt.savefig(resdir+"conc_per_geno_by_ind_"+str(cov)+"_"+ref_title+".pdf", dpi=fdpi, bbox_inches="tight")



	# plt.show()
	plt.close()


#Plot discordance and what types of errors there are in each category for files in dir IN SAME FILE
def  plot_nrds_gcts_stacked_allinds(dir, cov, refset):

	with open(dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json', 'r') as fp:
		all_tables = json.load(fp)



	lims = [7,7,11,10,12,20]
	# lims = [2,2,2,2,2,2]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]
	ind_offsets = [-2, -1, 0, 1, 2]


	subsets = ["kept", "dropped"]

	fig, ax = plt.subplots(figsize=(fwidth_small, fheight))
	# plt.tight_layout()
	# plt.suptitle("Genotype discordance by genotype in HC")
	# plt.subplots_adjust(top=0.93)
	# ax = plt.subplot(1,1,1)


	for ind in inds:

		counter = 0

		counter += 1

		# for this individual and filter : kept and dropped gcts
		all_gcts_bottom = []
		all_gcts_top = []

		for subset in subsets:

			gct_dict = all_tables[ind][subset][1]["# GCTs"]
			keys = ["Ref/Ref discordance", "Ref/Alt discordance", "Alt/Alt discordance","NRD"]
			keys = [u'RR Hom -> RA Het',  u'RR Hom -> AA Hom',u'RA Het -> RR Hom', u'RA Het -> AA Hom',  u'AA Hom -> RR Hom',  u'AA Hom -> RA Het']


			this_gcts_bottom = []
			this_gcts_top = []

			# RR
			num_RR_RA = gct_dict['RR Hom -> RA Het']
			num_RR_AA = gct_dict['RR Hom -> AA Hom']
			num_RR_RR = gct_dict['RR Hom -> RR Hom']

			print("RR-RA: " + str(num_RR_RA))
			print("RR-AA: " + str(num_RR_AA))

			RR_wrongs= int(num_RR_RA) + int(num_RR_AA)
			RR_all= int(num_RR_RA) + int(num_RR_AA) + int(num_RR_RR)


			if RR_wrongs > 0:
				perc_RR_RA =  float(num_RR_RA) / RR_wrongs
				perc_RR_AA =  float(num_RR_AA) / RR_wrongs
				RR_disc = 100* float(RR_wrongs) / RR_all

			else:
				perc_RR_RA = 0.0
				perc_RR_AA = 0.0
				RR_disc = 0.0

			perc_RR_RA = perc_RR_RA * RR_disc
			perc_RR_AA = perc_RR_AA * RR_disc


			print(perc_RR_RA)
			print(perc_RR_AA)

			# RA
			num_RA_RR = gct_dict['RA Het -> RR Hom']
			num_RA_AA = gct_dict['RA Het -> AA Hom']
			num_RA_RA = gct_dict['RA Het -> RA Het']


			RA_wrongs= int(num_RA_RR) + int(num_RA_AA)
			RA_all= int(num_RA_RR) + int(num_RA_AA) + int(num_RA_RA)
			RA_disc = 100* float(RA_wrongs) / RA_all


			perc_RA_RR =  float(num_RA_RR) / RA_wrongs
			perc_RA_AA =  float(num_RA_AA) / RA_wrongs


			perc_RA_RR = perc_RA_RR * RA_disc
			perc_RA_AA = perc_RA_AA * RA_disc


			# AA
			num_AA_RR = gct_dict['AA Hom -> RR Hom']
			num_AA_RA = gct_dict['AA Hom -> RA Het']
			num_AA_AA = gct_dict['AA Hom -> AA Hom']

			AA_wrongs= int(num_AA_RR) + int(num_AA_RA)
			AA_all= int(num_AA_RR) + int(num_AA_RA) + int(num_AA_AA)

			if AA_wrongs > 0:
				perc_AA_RR =  float(num_AA_RR) / AA_wrongs
				perc_AA_RA =  float(num_AA_RA) / AA_wrongs
				AA_disc = 100* float(AA_wrongs) / AA_all

			else:
				perc_AA_RR = 0
				perc_AA_RA = 0
				AA_disc = 0.0

			perc_AA_RR = perc_AA_RR * AA_disc
			perc_AA_RA = perc_AA_RA * AA_disc


			this_gcts_bottom = [perc_RR_RA, perc_RA_RR,perc_AA_RR]
			this_gcts_top = [perc_RR_AA, perc_RA_AA,perc_AA_RA]



			nrd_dict = all_tables[ind][subset][0]["# NRDs"]
			keys = ["Ref/Ref discordance", "Ref/Alt discordance", "Alt/Alt discordance","NRD"]
			this_nrds = []
			for k in keys:
				this_nrds.append(nrd_dict[k])

			# this_gcts_bottom.append(nrd_dict["NRD"])

			all_gcts_bottom.append(this_gcts_bottom)
			all_gcts_top.append(this_gcts_top)



			print(all_gcts_bottom)
			print(all_gcts_top)

			x_bottom = [0,15,30]
			x_top = [0,15,30]

			# x = range(len(keys))
			# x = map(lambda y : 10*y, x)

			print(x_bottom)
			width=1.0

			# ax.text(1, 2, "agg", ha='center',va='center')

			# for

			# plt.text(x_top, all_gcts_bottom[0][:-1], bottom_texts)
			# ax.text(x_top, map(lambda y: y+5 , all_gcts_top[0][:-1]), top_texts)



		RRcol = "steelblue"
		RAcol = "darkmagenta"
		AAcol = "green"

		RRh = ""
		RAh = "/"
		AAh = "."


		bottom_texts = ["RA", "RR", "RR"]
		bottom_colors = [RAcol, RRcol, RRcol, "black"]
		bottom_hatches = [RAh, RRh, RRh]

		top_texts = ["AA", "AA", "RA"]
		top_colors = [AAcol, AAcol, RAcol, "black"]
		top_hatches = [AAh, AAh, RAh]

		#kept
		#label=subsets[0]
		# plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 - (width/2), x_bottom),all_gcts_bottom[0], width = width, color=ind_palette[inds.index(ind)])
		# plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 - (width/2), x_top),all_gcts_top[0], bottom=all_gcts_bottom[0], width = width, color=ind_palette[inds.index(ind)], alpha=0.8)
		#
		# #dropped
		# plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 + (width/2), x_bottom),all_gcts_bottom[1], width = width, color=ind_palette[inds.index(ind)], alpha=0.4)
		# plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 + (width/2), x_top),all_gcts_top[1], bottom=all_gcts_bottom[1], width = width, color=ind_palette[inds.index(ind)], alpha=0.2)

				#kept
		label=subsets[0]
		bsb = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 - (width/2), x_bottom),all_gcts_bottom[0], width = width, color=ind_palette[inds.index(ind)], alpha=0.8)
		bst = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 - (width/2), x_top),all_gcts_top[0], bottom=all_gcts_bottom[0], width = width, color=ind_palette[inds.index(ind)], alpha=0.8)

		# for b,h in zip(bsb, bottom_hatches):
			# b.set_edgecolor(ind_palette[inds.index(ind)])
			# b.set_alpha(0.1)
			# b.set_hatch(h)

		# for b,h in zip(bst, top_hatches):
			# b.set_edgecolor(ind_palette[inds.index(ind)])
			# b.set_edgecolor(ind_palette[inds.index(ind)])
			# b.set_hatch(h)

		#dropped
		bsb = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 + (width/2), x_bottom),all_gcts_bottom[1], width = width, color=ind_palette[inds.index(ind)], alpha=0.5)
		bst = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 + (width/2), x_top),all_gcts_top[1], bottom=all_gcts_bottom[1], width = width, color=ind_palette[inds.index(ind)], alpha=0.5)
		#
		# for b,h in zip(bsb, bottom_hatches):
		# 	b.set_edgecolor(ind_palette[inds.index(ind)])
		# 	b.set_hatch("/")
		#
		# for b,h in zip(bst, top_hatches):
		# 	b.set_edgecolor(ind_palette[inds.index(ind)])
		# 	b.set_hatch("/")


		# adj=1.5
		# locs_kept = map(lambda y : y-(width/2), x_top)
		# for i in range(len(all_gcts_top[0])):
		# 	plt.text(locs_kept[i]-adj, all_gcts_bottom[0][i]+all_gcts_top[0][i]+0.05,'{0:.2f}'.format(all_gcts_bottom[0][i]+all_gcts_top[0][i]),size=7)
		#
		# locs_dropped = map(lambda y : y+(width/2), x_top)
		# for i in range(len(all_gcts_top[1])):
		# 	plt.text(locs_dropped[i]-adj, all_gcts_bottom[1][i]+all_gcts_top[1][i]+0.05,'{0:.2f}'.format(all_gcts_bottom[1][i]+all_gcts_top[1][i]),size=7)
		#



		# if counter == 1:
		# 	# plt.title(subset)
		# 	# plt.title("All Markers")
		# 	plt.bar([0],[0], label="RR", color=RRcol)
		# 	plt.bar([0],[0], label="RA", color=RAcol)
		# 	plt.bar([0],[0], label="AA", color=AAcol)
		# 	plt.legend(fontsize=8)
		# if counter == 2:
		# 	plt.title("Quality Filtered")
		#
		# if counter == 3:
		# 	plt.title("Trusted Heterozygotes")

	xlabels = ["Hom-Ref", "Het", "Hom-Alt","Non-Ref"]

	plt.ylabel("Genotype discordance")
	plt.xlabel("True genotype")

	plt.xticks(x_bottom,xlabels)
	plt.ylim(ymax=1)

	# plt.xticks(rotation=-15)


	# plt.legend()
	print(dir+'_'+refset+'/'+"Genotype_discordance_stacked_allinds"+str(cov)+"_"+refset+".pdf")
	plt.savefig(dir+'_'+refset+'/'+"Genotype_discordance_stacked_allinds"+str(cov)+"_"+refset+".pdf", dpi=fdpi, bbox_inches="tight")
	plt.show()
	plt.close()



#Plot NRD for FULL, EUR, EUR ind
def  plot_nrds_allinds_compref(cov):


	lims = [7,7,11,10,12,20]
	# lims = [2,2,2,2,2,2]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]
	ind_offsets = [-2, -1, 0, 1, 2]


	subsets = ["kept", "dropped"]

	fig, ax = plt.subplots(figsize=(fwidth_small, fheight))
	# plt.tight_layout()
	# plt.suptitle("Genotype discordance by genotype in HC")

	ax = plt.subplot(1,1,1)


	for ind in inds:

		counter = 0
		counter += 1

		all_nrds_this_ind = []

		for subset in subsets:

			nrds_this_subset = []


			with open('stats_rmts_61_5_impnofilter_FULL/all_stats_'+str(cov)+'_FULL.json', 'r') as fp:
				all_tables = json.load(fp)
				nrd_dict = all_tables[ind][subset][0]["# NRDs"]
				nrds_this_subset.append(nrd_dict["NRD"])


			with open('stats_rmts_61_5_impnofilter_EUR/all_stats_'+str(cov)+'_EUR.json', 'r') as fp:
				all_tables = json.load(fp)
				nrd_dict = all_tables[ind][subset][0]["# NRDs"]
				nrds_this_subset.append(nrd_dict["NRD"])

			with open('stats_rmts_ind_EUR/all_stats_'+str(cov)+'_EUR.json', 'r') as fp:
				all_tables = json.load(fp)
				nrd_dict = all_tables[ind][subset][0]["# NRDs"]
				nrds_this_subset.append(nrd_dict["NRD"])


			all_nrds_this_ind.append(nrds_this_subset)

			x_bottom = [0,15,30]
			x_top = [0,15,30]
			print(x_bottom)
			width=1.0


			RRcol = "steelblue"
			RAcol = "darkmagenta"
			AAcol = "green"

			RRh = ""
			RAh = "/"
			AAh = "."


			bottom_texts = ["RA", "RR", "RR"]
			tcols = [RAcol, RRcol, RRcol, "red", "pink"]
			bottom_hatches = [RAh, RRh, RRh]

			top_texts = ["AA", "AA", "RA"]
			top_colors = [AAcol, AAcol, RAcol, "black"]
			top_hatches = [AAh, AAh, RAh]

			#kept
			#label=subsets[0]
			# plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 - (width/2), x_bottom),all_gcts_bottom[0], width = width, color=ind_palette[inds.index(ind)])
			# plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 - (width/2), x_top),all_gcts_top[0], bottom=all_gcts_bottom[0], width = width, color=ind_palette[inds.index(ind)], alpha=0.8)
			#
			# #dropped
			# plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 + (width/2), x_bottom),all_gcts_bottom[1], width = width, color=ind_palette[inds.index(ind)], alpha=0.4)
			# plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 + (width/2), x_top),all_gcts_top[1], bottom=all_gcts_bottom[1], width = width, color=ind_palette[inds.index(ind)], alpha=0.2)

		#kept
		label=subsets[0]
		bsl = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 - (width/2), x_bottom),all_nrds_this_ind[0], width = width, color=ind_palette[inds.index(ind)], alpha=0.8)
		#dropped
		bsb = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 + (width/2), x_bottom),all_nrds_this_ind[1], width = width, color=ind_palette[inds.index(ind)], alpha=0.5)

		# for b in bsb:
			# b.set_alpha(0.0)
			# b.set_edgecolor(ind_palette[inds.index(ind)])
			# b.set_edgecolor(matplotlib.colors.ColorConverter.to_rgba(ind_palette[inds.index(ind)], alpha=0.0))
			# # newcol = matplotlib.colors.ColorConverter.to_rgba('black', alpha=.1)
			# # newcol = matplotlib.colors.ColorConverter.to_rgba(matplotlib.colors.to_rgba(ind_palette[inds.index(ind)], alpha=0.5), alpha=0.5)
			# # newcol2 = matplotlib.colors.ColorConverter.to_rgba(matplotlib.colors.to_rgba(ind_palette[inds.index(ind)]))
			# #
			# # print(newcol)
			# # print(newcol2)
			# #
			# # print(matplotlib.colors.to_hex(newcol))
			# # print(matplotlib.colors.to_hex(newcol2))
			# #
			# # print(ind_palette[inds.index(ind)])
			# b.set_edgecolor(newcol)

			# b.set_hatch('/')




	xlabels = ["FULL + ANCIENT", "EUR + ANCIENT", "EUR + SINGLE"]
	plt.xlabel("Imputation configuration")
	plt.ylabel("NRD")
	plt.xticks(x_bottom,xlabels)
	plt.ylim(ymax=7)

	# plt.xticks(rotation=-15)


	# plt.legend()
	print("NRD_imp_config"+str(cov)+".pdf")
	plt.savefig("NRD_imp_config"+str(cov)+".svg", dpi=fdpi, bbox_inches="tight")
	# plt.show()
	# plt.close()

#Plot NRD for FULL, EUR, EUR ind
def  plot_nrds_allinds_compref_by_overlap(cov):


	lims = [7,7,11,10,12,20]
	# lims = [2,2,2,2,2,2]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]
	ind_offsets = [-2, -1, 0, 1, 2]


	subsets = ["kept", "dropped"]

	fig, ax = plt.subplots(figsize=(fwidth_med, fheight_med))

	# plt.tight_layout()
	# plt.suptitle("Genotype discordance by genotype in HC")

	ax = plt.subplot(1,1,1)


	for ind in inds:

		counter = 0
		counter += 1

		all_nrds_this_ind = []

		for subset in subsets:

			nrds_this_subset = []


			with open('stats_rmts_61_5_impnofilter_FULL/all_stats_'+str(cov)+'_FULL.json', 'r') as fp:
				all_tables = json.load(fp)
				nrd_dict = all_tables[ind][subset][0]["# NRDs"]
				nrds_this_subset.append(nrd_dict["NRD"])


			with open('stats_rmts_61_5_impnofilter_EUR/all_stats_'+str(cov)+'_EUR.json', 'r') as fp:
				all_tables = json.load(fp)
				nrd_dict = all_tables[ind][subset][0]["# NRDs"]
				nrds_this_subset.append(nrd_dict["NRD"])

			with open('stats_rmts_ind_EUR/all_stats_'+str(cov)+'_EUR.json', 'r') as fp:
				all_tables = json.load(fp)
				nrd_dict = all_tables[ind][subset][0]["# NRDs"]
				nrds_this_subset.append(nrd_dict["NRD"])


			all_nrds_this_ind.append(nrds_this_subset)

			x_bottom = [0,15,30]
			x_top = [0,15,30]
			print(x_bottom)
			width=1.0


			RRcol = "steelblue"
			RAcol = "darkmagenta"
			AAcol = "green"

			RRh = ""
			RAh = "/"
			AAh = "."


			bottom_texts = ["RA", "RR", "RR"]
			tcols = [RAcol, RRcol, RRcol, "red", "pink"]
			bottom_hatches = [RAh, RRh, RRh]

			top_texts = ["AA", "AA", "RA"]
			top_colors = [AAcol, AAcol, RAcol, "black"]
			top_hatches = [AAh, AAh, RAh]

			#kept
			#label=subsets[0]
			# plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 - (width/2), x_bottom),all_gcts_bottom[0], width = width, color=ind_palette[inds.index(ind)])
			# plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 - (width/2), x_top),all_gcts_top[0], bottom=all_gcts_bottom[0], width = width, color=ind_palette[inds.index(ind)], alpha=0.8)
			#
			# #dropped
			# plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 + (width/2), x_bottom),all_gcts_bottom[1], width = width, color=ind_palette[inds.index(ind)], alpha=0.4)
			# plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)]*2 + (width/2), x_top),all_gcts_top[1], bottom=all_gcts_bottom[1], width = width, color=ind_palette[inds.index(ind)], alpha=0.2)

		#kept
		label=subsets[0]
		bsl = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)] - (width/2), x_bottom),all_nrds_this_ind[0], width = width, color=ind_palette[inds.index(ind)], alpha=0.8)
		#dropped
		bsb = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)] + 5 + (width/2), x_bottom),all_nrds_this_ind[1], width = width, color=ind_palette[inds.index(ind)], alpha=0.5)

		# for b in bsb:
			# b.set_alpha(0.0)
			# b.set_edgecolor(ind_palette[inds.index(ind)])
			# b.set_edgecolor(matplotlib.colors.ColorConverter.to_rgba(ind_palette[inds.index(ind)], alpha=0.0))
			# # newcol = matplotlib.colors.ColorConverter.to_rgba('black', alpha=.1)
			# # newcol = matplotlib.colors.ColorConverter.to_rgba(matplotlib.colors.to_rgba(ind_palette[inds.index(ind)], alpha=0.5), alpha=0.5)
			# # newcol2 = matplotlib.colors.ColorConverter.to_rgba(matplotlib.colors.to_rgba(ind_palette[inds.index(ind)]))
			# #
			# # print(newcol)
			# # print(newcol2)
			# #
			# # print(matplotlib.colors.to_hex(newcol))
			# # print(matplotlib.colors.to_hex(newcol2))
			# #
			# # print(ind_palette[inds.index(ind)])
			# b.set_edgecolor(newcol)

			# b.set_hatch('/')




	xlabels = ["FULL + ANCIENT", "EUR + ANCIENT", "EUR + SINGLE"]
	plt.xlabel("Imputation configuration")
	plt.ylabel("NRD")
	plt.xticks([2.5,17.5,32.5], xlabels)
	plt.ylim(ymax=7)

	# plt.xticks(rotation=-15)


	# plt.legend()
	print("NRD_imp_config"+str(cov)+".pdf")
	plt.savefig("NRD_imp_config_byoverlap_"+str(cov)+".svg", dpi=fdpi, bbox_inches="tight")
	plt.savefig("NRD_imp_config_byoverlap_"+str(cov)+".pdf", dpi=fdpi, bbox_inches="tight")

	plt.show()
	# plt.close()



#Plot concordance for FULL, EUR, EUR ind
def  plot_concordance_per_refset_by_ind(cov, genos, ymin=0.99, ymax=1.0,resdir="", grid=False,wid=6, heig=4,if_ind_legend=False):

	if genos == ['RR Hom','RA Het','AA Hom']:
		genos_title = "all_genos"

	if genos == ['RA Het']:
		genos_title = "hets"

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]
	ind_offsets = [-2, -1, 0, 1, 2]

	subsets = ["kept", "dropped"]

	fig, ax = plt.subplots(figsize=(wid, heig))

	# ax = plt.subplot(1,1,1)
	if grid:
		plt.grid(b=True, which='major',axis='y', linewidth=lw3)
		plt.grid(b=True, which='minor',axis='y', linewidth=lw3)


	for ind in inds:

		counter = 0
		counter += 1

		all_concs_this_ind = []

		for subset in subsets:

			concs_this_subset = []


			concordance = get_concordance_from_gcss('stats_rmts_ind_nofilter_anno','EUR',cov,ind,subset,genos)
			concs_this_subset.append(concordance)

			concordance = get_concordance_from_gcss('stats_rmts_61_5_impnofilter_anno','EUR',cov,ind,subset,genos)
			concs_this_subset.append(concordance)

			concordance = get_concordance_from_gcss('stats_rmts_61_5_impnofilter_anno','FULL',cov,ind,subset,genos)
			concs_this_subset.append(concordance)


			all_concs_this_ind.append(concs_this_subset)

			x_bottom = [0,15,30]
			width=1.0


		#kept
		label=subsets[0]
		bsl = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)] - (width/2), x_bottom),all_concs_this_ind[0], width = width, color=ind_palette[inds.index(ind)], alpha=kept_alpha)
		#dropped
		bsb = plt.bar(map(lambda y : y + ind_offsets[inds.index(ind)] + 5 + (width/2), x_bottom),all_concs_this_ind[1], width = width, color=ind_palette[inds.index(ind)], alpha=dropped_alpha)

		############## text above bars ###############
		# #
		# adj=0.0004
		# # locs_kept = map(lambda y : y-(width/2), x_bottom)
		# rot = 45
		# for i in range(len(all_concs_this_ind[0])):
		# 	# plt.text(locs_kept[i]-adj, all_concs_this_ind[0][i]+0.05,'{0:.2f}'.format(all_concs_this_ind[0][i]),size=70)
		# 	plt.text(map(lambda y : y + ind_offsets[inds.index(ind)] - (width/2), x_bottom)[i], all_concs_this_ind[0][i]+ adj,'{0:.4f}'.format(all_concs_this_ind[0][i]) ,size=7, rotation=rot)
		#
		# locs_dropped = map(lambda y : y+(width/2), x_bottom)
		# for i in range(len(all_concs_this_ind[1])):
		# 	# plt.text(locs_dropped[i]-adj, all_concs_this_ind[1][i]+0.05,'{0:.2f}'.format(all_concs_this_ind[1][i]),size=70)
		# 	plt.text(map(lambda y : y + ind_offsets[inds.index(ind)] + 5 + (width/2), x_bottom)[i], all_concs_this_ind[1][i]+ adj,'{0:.4f}'.format(all_concs_this_ind[1][i]) ,size=7, rotation=rot)


	xlabels = ["1", "2", "3"]
	plt.xlabel("Imputation configuration")
	plt.ylabel("Genotype concordance")
	plt.ylim(ymin=ymin, ymax=ymax)
	plt.xticks([2.5,17.5,32.5], xlabels)
	plt.setp(ax.spines.values(), linewidth=0)

	# legend_properties = {'weight':650}
	if if_ind_legend:
		# color legend  bbox_to_anchor=(0.0,0.0,1.1,1.0)
		for ind in range(len(inds)):
				# plt.scatter([-1],[-1], color=ind_palette[ind], s=0, label=inds[ind])
				plt.bar([-1],[-1], color=ind_palette[ind], width=0.1, label=inds[ind])


		# leg = plt.legend(fontsize=4.5, loc=1,borderpad=0.1,handlelength=0,bbox_to_anchor=(1.011,1.025), labelspacing=0.2,framealpha=1.0)
		# leg.get_frame().set_linewidth(lw3)
		# leg_texts = leg.get_texts() # list of matplotlib Text instances.
		#
		# c = 0
		# for leg_text in leg_texts:
		# 	leg_text.set_color(ind_palette[c])
		# 	leg_text.set_fontweight("semibold")
		# 	leg_text.set_position((-1.9,0))
		# 	c+=1


		leg = plt.legend(fontsize=4.5, ncol=5, loc=9, borderpad=0.6, labelspacing=0.8,framealpha=1.0,handlelength=3, mode="expand",bbox_to_anchor=(0.0,1.17,1.0,0.0),handletextpad=0.5)
		leg.get_frame().set_linewidth(lw3)
		leg_texts = leg.get_texts() # list of matplotlib Text instances.

		c = 0
		for leg_text in leg_texts:
			# leg_text.set_color(ind_palette[c])
			# leg_text.set_fontweight("semibold")
			# leg_text.set_position((-1.9,0))
			c+=1

		for leg_patch in leg.get_patches():
			leg_patch.set_linewidth(0)




	if if_title:
		plt.title("cov: " + str(cov) + " genotyped: " + genos_title)

	 # (left, bottom, right, top)
	plt.tight_layout(pad=0,rect=[0.0,0.0,1.0,0.9])
	plt.savefig(resdir+"conc_per_ref_by_ind_"+str(cov)+"_"+genos_title+".pdf", dpi=fdpi)


	# plt.savefig(resdir+"conc_per_ref_by_ind_"+str(cov)+"_"+genos_title+".pdf", dpi=fdpi, bbox_inches="tight")
	# plt.show()
	plt.close()


#Plot discordance and what types of errors there are in each category for files in dir
def plot_nrds_gcts_stacked_comp2(dir1,dir2,title1,title2):

	inds = ["ans17", "car16", "LBK", "Loschbour", "ne1", "ne1CP"]
	lims = [7,7,11,10,12,20]
	# inds = ["ans17", "car16", "LBK", "Loschbour", "ne1"]

	# inds = ["LBK"]


	print(inds)

	filters = ["all", "filtered", "thets"]
	# filters = ["all"]


	subsets = ["kept", "dropped"]


	for ind in inds:

		fig, ax = plt.subplots(figsize=(9,12))
		# plt.tight_layout()
		plt.suptitle("Genotype discordance by genotype in HC - "+ind)
		plt.subplots_adjust(top=0.93)


		# Dir 1 : left hand side of plot, ODD indices

		with open(dir1 + 'all_tables.json', 'r') as fp:
			all_tables = json.load(fp)
		counter = -1
		for filter in filters:
			counter += 2
			ax = plt.subplot(3,2,counter)

			# for this individual and filter : kept and dropped gcts
			all_gcts_bottom = []
			all_gcts_top = []

			for subset in subsets:

				gct_dict = all_tables[ind][filter][subset][1]["# GCTs"]
				keys = ["Ref/Ref discordance", "Ref/Alt discordance", "Alt/Alt discordance","NRD"]
				keys = [u'RR Hom -> RA Het',  u'RR Hom -> AA Hom',u'RA Het -> RR Hom', u'RA Het -> AA Hom',  u'AA Hom -> RR Hom',  u'AA Hom -> RA Het']


				this_gcts_bottom = []
				this_gcts_top = []

				# for k in keys:
				# 	this_gcts_bottom.append(gct_dict[k])

				# RR
				num_RR_RA = gct_dict['RR Hom -> RA Het']
				num_RR_AA = gct_dict['RR Hom -> AA Hom']
				num_RR_RR = gct_dict['RR Hom -> RR Hom']

				print("RR-RA: " + str(num_RR_RA))
				print("RR-AA: " + str(num_RR_AA))

				RR_wrongs= int(num_RR_RA) + int(num_RR_AA)
				RR_all= int(num_RR_RA) + int(num_RR_AA) + int(num_RR_RR)


				if RR_wrongs > 0:
					perc_RR_RA =  float(num_RR_RA) / RR_wrongs
					perc_RR_AA =  float(num_RR_AA) / RR_wrongs
					RR_disc = 100* float(RR_wrongs) / RR_all

				else:
					perc_RR_RA = 0.0
					perc_RR_AA = 0.0
					RR_disc = 0.0

				perc_RR_RA = perc_RR_RA * RR_disc
				perc_RR_AA = perc_RR_AA * RR_disc


				print(perc_RR_RA)
				print(perc_RR_AA)

				# RA
				num_RA_RR = gct_dict['RA Het -> RR Hom']
				num_RA_AA = gct_dict['RA Het -> AA Hom']
				num_RA_RA = gct_dict['RA Het -> RA Het']


				RA_wrongs= int(num_RA_RR) + int(num_RA_AA)
				RA_all= int(num_RA_RR) + int(num_RA_AA) + int(num_RA_RA)
				RA_disc = 100* float(RA_wrongs) / RA_all


				perc_RA_RR =  float(num_RA_RR) / RA_wrongs
				perc_RA_AA =  float(num_RA_AA) / RA_wrongs


				perc_RA_RR = perc_RA_RR * RA_disc
				perc_RA_AA = perc_RA_AA * RA_disc


				# AA
				num_AA_RR = gct_dict['AA Hom -> RR Hom']
				num_AA_RA = gct_dict['AA Hom -> RA Het']
				num_AA_AA = gct_dict['AA Hom -> AA Hom']

				AA_wrongs= int(num_AA_RR) + int(num_AA_RA)
				AA_all= int(num_AA_RR) + int(num_AA_RA) + int(num_AA_AA)

				if AA_wrongs > 0:
					perc_AA_RR =  float(num_AA_RR) / AA_wrongs
					perc_AA_RA =  float(num_AA_RA) / AA_wrongs
					AA_disc = 100* float(AA_wrongs) / AA_all

				else:
					perc_AA_RR = 0
					perc_AA_RA = 0
					AA_disc = 0.0

				perc_AA_RR = perc_AA_RR * AA_disc
				perc_AA_RA = perc_AA_RA * AA_disc


				this_gcts_bottom = [perc_RR_RA, perc_RA_RR,perc_AA_RR]
				this_gcts_top = [perc_RR_AA, perc_RA_AA,perc_AA_RA]



				nrd_dict = all_tables[ind][filter][subset][0]["# NRDs"]
				keys = ["Ref/Ref discordance", "Ref/Alt discordance", "Alt/Alt discordance","NRD"]
				this_nrds = []
				for k in keys:
					this_nrds.append(nrd_dict[k])

				this_gcts_bottom.append(nrd_dict["NRD"])

				all_gcts_bottom.append(this_gcts_bottom)
				all_gcts_top.append(this_gcts_top)

			RRcol = "steelblue"
			RAcol = "darkmagenta"
			AAcol = "green"

			bottom_texts = ["RA", "RR", "RR"]
			bottom_colors = [RAcol, RRcol, RRcol, "black"]

			top_texts = ["AA", "AA", "RA"]
			top_colors = [AAcol, AAcol, RAcol, "black"]

			print(all_gcts_bottom)
			print(all_gcts_top)

			x_bottom = [0,10,20,30]
			x_top = [0,10,20]

			# x = range(len(keys))
			# x = map(lambda y : 10*y, x)

			print(x_bottom)
			width=3.0

			# ax.text(1, 2, "agg", ha='center',va='center')

			# for

			# plt.text(x_top, all_gcts_bottom[0][:-1], bottom_texts)
			# ax.text(x_top, map(lambda y: y+5 , all_gcts_top[0][:-1]), top_texts)

			#kept
			#label=subsets[0]
			plt.bar(map(lambda y : y-(width/2), x_bottom),all_gcts_bottom[0], width = width, color=bottom_colors)
			plt.bar(map(lambda y : y-(width/2), x_top),all_gcts_top[0], bottom=all_gcts_bottom[0][:-1], width = width, color=top_colors)

			#dropped
			plt.bar(map(lambda y : y+(width/2), x_bottom),all_gcts_bottom[1], width = width, color=bottom_colors, alpha=0.5)
			plt.bar(map(lambda y : y+(width/2), x_top),all_gcts_top[1], bottom=all_gcts_bottom[1][:-1], width = width, color=top_colors, alpha=0.5)

			adj=1.5
			locs_kept = map(lambda y : y-(width/2), x_top)
			for i in range(len(all_gcts_top[0])):
				plt.text(locs_kept[i]-adj, all_gcts_bottom[0][i]+all_gcts_top[0][i]+0.05,'{0:.2f}'.format(all_gcts_bottom[0][i]+all_gcts_top[0][i]),size=7)

			locs_dropped = map(lambda y : y+(width/2), x_top)
			for i in range(len(all_gcts_top[1])):
				plt.text(locs_dropped[i]-adj, all_gcts_bottom[1][i]+all_gcts_top[1][i]+0.05,'{0:.2f}'.format(all_gcts_bottom[1][i]+all_gcts_top[1][i]),size=7)



			if counter == 1:
				# plt.title(subset)
				plt.title(title1 + "\nAll Markers")
				plt.bar([0],[0], label="RR", color=RRcol)
				plt.bar([0],[0], label="RA", color=RAcol)
				plt.bar([0],[0], label="AA", color=AAcol)
				plt.legend(fontsize=8)
			if counter == 3:
			# 	plt.title(subset)

				plt.title("Quality Filtered")

			if counter == 5:
			# 	plt.title(subset)
				plt.title("Trusted Heterozygotes")

			# xlabels = ["R/R discordance", "R/A discordance", "A/A discordance","Non-Ref Discordance"]
			xlabels = ["Hom-Ref", "Het", "Hom-Alt","Non-Ref"]

			# xlabels = keys
			# print(gct_dict.keys())

			plt.ylabel("% wrong genotype in imputed")
			plt.xticks(x_bottom,xlabels,size=7)
			plt.ylim(ymax=lims[inds.index(ind)])
			# plt.xticks(rotation=-15)


		# Dir 2 : right hand side of plot, EVEN indices

		with open(dir2 + 'all_tables.json', 'r') as fp:
			all_tables = json.load(fp)
		counter = 0
		for filter in filters:
			counter += 2
			ax = plt.subplot(3,2,counter)

			# for this individual and filter : kept and dropped gcts
			all_gcts_bottom = []
			all_gcts_top = []

			for subset in subsets:

				gct_dict = all_tables[ind][filter][subset][1]["# GCTs"]
				keys = ["Ref/Ref discordance", "Ref/Alt discordance", "Alt/Alt discordance","NRD"]
				keys = [u'RR Hom -> RA Het',  u'RR Hom -> AA Hom',u'RA Het -> RR Hom', u'RA Het -> AA Hom',  u'AA Hom -> RR Hom',  u'AA Hom -> RA Het']


				this_gcts_bottom = []
				this_gcts_top = []

				# for k in keys:
				# 	this_gcts_bottom.append(gct_dict[k])

				# RR
				num_RR_RA = gct_dict['RR Hom -> RA Het']
				num_RR_AA = gct_dict['RR Hom -> AA Hom']
				num_RR_RR = gct_dict['RR Hom -> RR Hom']

				print("RR-RA: " + str(num_RR_RA))
				print("RR-AA: " + str(num_RR_AA))

				RR_wrongs= int(num_RR_RA) + int(num_RR_AA)
				RR_all= int(num_RR_RA) + int(num_RR_AA) + int(num_RR_RR)


				if RR_wrongs > 0:
					perc_RR_RA =  float(num_RR_RA) / RR_wrongs
					perc_RR_AA =  float(num_RR_AA) / RR_wrongs
					RR_disc = 100* float(RR_wrongs) / RR_all

				else:
					perc_RR_RA = 0.0
					perc_RR_AA = 0.0
					RR_disc = 0.0

				perc_RR_RA = perc_RR_RA * RR_disc
				perc_RR_AA = perc_RR_AA * RR_disc


				print(perc_RR_RA)
				print(perc_RR_AA)

				# RA
				num_RA_RR = gct_dict['RA Het -> RR Hom']
				num_RA_AA = gct_dict['RA Het -> AA Hom']
				num_RA_RA = gct_dict['RA Het -> RA Het']


				RA_wrongs= int(num_RA_RR) + int(num_RA_AA)
				RA_all= int(num_RA_RR) + int(num_RA_AA) + int(num_RA_RA)
				RA_disc = 100* float(RA_wrongs) / RA_all


				perc_RA_RR =  float(num_RA_RR) / RA_wrongs
				perc_RA_AA =  float(num_RA_AA) / RA_wrongs


				perc_RA_RR = perc_RA_RR * RA_disc
				perc_RA_AA = perc_RA_AA * RA_disc


				# AA
				num_AA_RR = gct_dict['AA Hom -> RR Hom']
				num_AA_RA = gct_dict['AA Hom -> RA Het']
				num_AA_AA = gct_dict['AA Hom -> AA Hom']

				AA_wrongs= int(num_AA_RR) + int(num_AA_RA)
				AA_all= int(num_AA_RR) + int(num_AA_RA) + int(num_AA_AA)

				if AA_wrongs > 0:
					perc_AA_RR =  float(num_AA_RR) / AA_wrongs
					perc_AA_RA =  float(num_AA_RA) / AA_wrongs
					AA_disc = 100* float(AA_wrongs) / AA_all

				else:
					perc_AA_RR = 0
					perc_AA_RA = 0
					AA_disc = 0.0

				perc_AA_RR = perc_AA_RR * AA_disc
				perc_AA_RA = perc_AA_RA * AA_disc


				this_gcts_bottom = [perc_RR_RA, perc_RA_RR,perc_AA_RR]
				this_gcts_top = [perc_RR_AA, perc_RA_AA,perc_AA_RA]



				nrd_dict = all_tables[ind][filter][subset][0]["# NRDs"]
				keys = ["Ref/Ref discordance", "Ref/Alt discordance", "Alt/Alt discordance","NRD"]
				this_nrds = []
				for k in keys:
					this_nrds.append(nrd_dict[k])

				this_gcts_bottom.append(nrd_dict["NRD"])

				all_gcts_bottom.append(this_gcts_bottom)
				all_gcts_top.append(this_gcts_top)

			RRcol = "steelblue"
			RAcol = "darkmagenta"
			AAcol = "green"

			bottom_texts = ["RA", "RR", "RR"]
			bottom_colors = [RAcol, RRcol, RRcol, "black"]

			top_texts = ["AA", "AA", "RA"]
			top_colors = [AAcol, AAcol, RAcol, "black"]

			print(all_gcts_bottom)
			print(all_gcts_top)

			x_bottom = [0,10,20,30]
			x_top = [0,10,20]

			# x = range(len(keys))
			# x = map(lambda y : 10*y, x)

			print(x_bottom)
			width=3.0

			# ax.text(1, 2, "agg", ha='center',va='center')

			# for

			# plt.text(x_top, all_gcts_bottom[0][:-1], bottom_texts)
			# ax.text(x_top, map(lambda y: y+5 , all_gcts_top[0][:-1]), top_texts)


			#kept
			#label=subsets[0]
			plt.bar(map(lambda y : y-(width/2), x_bottom),all_gcts_bottom[0], width = width, color=bottom_colors)
			plt.bar(map(lambda y : y-(width/2), x_top),all_gcts_top[0], bottom=all_gcts_bottom[0][:-1], width = width, color=top_colors)

			#dropped
			plt.bar(map(lambda y : y+(width/2), x_bottom),all_gcts_bottom[1], width = width, color=bottom_colors, alpha=0.5)
			plt.bar(map(lambda y : y+(width/2), x_top),all_gcts_top[1], bottom=all_gcts_bottom[1][:-1], width = width, color=top_colors, alpha=0.5)


			adj=1.5
			locs_kept = map(lambda y : y-(width/2), x_top)
			for i in range(len(all_gcts_top[0])):
				plt.text(locs_kept[i]-adj, all_gcts_bottom[0][i]+all_gcts_top[0][i]+0.05,'{0:.2f}'.format(all_gcts_bottom[0][i]+all_gcts_top[0][i]),size=7)

			locs_dropped = map(lambda y : y+(width/2), x_top)
			for i in range(len(all_gcts_top[1])):
				plt.text(locs_dropped[i]-adj, all_gcts_bottom[1][i]+all_gcts_top[1][i]+0.05,'{0:.2f}'.format(all_gcts_bottom[1][i]+all_gcts_top[1][i]),size=7)



			if counter == 2:
				# plt.title(subset)
				plt.title(title2 + "\nAll Markers")
				plt.bar([0],[0], label="RR", color=RRcol)
				plt.bar([0],[0], label="RA", color=RAcol)
				plt.bar([0],[0], label="AA", color=AAcol)
				plt.legend(fontsize=8)
			if counter == 4:
			# 	plt.title(subset)

				plt.title("Quality Filtered")

			if counter == 6:
			# 	plt.title(subset)
				plt.title("Trusted Heterozygotes")

			# xlabels = ["R/R discordance", "R/A discordance", "A/A discordance","Non-Ref Discordance"]
			xlabels = ["Hom-Ref", "Het", "Hom-Alt","Non-Ref"]

			# xlabels = keys
			# print(gct_dict.keys())

			plt.ylabel("% wrong genotype in imputed")
			plt.xticks(x_bottom,xlabels,size=7)
			plt.ylim(ymax=lims[inds.index(ind)])
			# plt.xticks(rotation=-15)


		# plt.legend()
		plt.savefig(dir1+"Genotype_discordance_stacked_"+title1+"_"+title2+"_"+ind+".pdf")
		# plt.show()
#

#Plot number of markers in the different coverages for the different subsets
def  plot_nummarkers_per_cov(dir, refset):

	# marker_nums[ind][0][c] : num markers for ind for KEPT markers at cov c
	# marker_nums[ind][1][c] : num markers for ind for DROPPED markers at cov c

	marker_nums=[[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]
	# inds = ["ans17", "sf12", "ne1"]


	subsets = ["kept", "dropped"]
	covs = [0.25, 0.5, 0.75, 1,  1.25, 1.5, 1.75]

	for cov in covs:
		print(dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json')
		with open(dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json', 'r') as fp:
			all_tables = json.load(fp)

		for ind in range(len(inds)):

			for subset in range(len(subsets)):

				num_markers = 0
				gcss_dict = all_tables[inds[ind]][subsets[subset]][2]["# GCsS"]
				num_markers += float(gcss_dict['AA Hom matches'])
				num_markers += float(gcss_dict['AA Hom mismatches'])

				num_markers += float(gcss_dict['RA Het matches'])
				num_markers += float(gcss_dict['RA Het mismatches'])

				num_markers += float(gcss_dict['RR Hom matches'])
				num_markers += float(gcss_dict['RR Hom mismatches'])

				marker_nums[ind][subset].append(num_markers)



	for ind in range(len(inds)):
		plt.plot(covs, marker_nums[ind][0], linestyle="solid", color=ind_palette[ind], label=inds[ind])
		plt.plot(covs, marker_nums[ind][1], linestyle="dashed", dashes =dashstyle, color=ind_palette[ind])
		plt.plot(covs, np.add(marker_nums[ind][1],marker_nums[ind][0]) , linewidth=4, color=ind_palette[ind])

		plt.xticks(covs)
		plt.legend()

	plt.show()
	# plt.savefig(dir+'_'+refset+'/'+"NRDs_per_cov_"+refset+".pdf")



#Plot Non-reference discordance for 5 individuals for coverages  on x-axis
def  plot_nrds_per_cov(dir, refset):

	# nrds[ind][0][c] : non-reference discrodance for ind for KEPT markers at cov c
	# nrds[ind][1][c] : non-reference discrodance for ind for DROPPED markers at cov c

	nrds=[[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]
	# inds = ["ans17", "sf12", "ne1"]

	# plt.figure(figsize=(fwidth, fheight))
	plt.figure(figsize=(fwidth_med, fheight_med))
	subsets = ["kept", "dropped"]
	covs = [0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75,2]

	for cov in covs:
		print("Reading " + dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json')
		with open(dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json', 'r') as fp:
			all_tables = json.load(fp)

		for ind in range(len(inds)):

			for subset in range(len(subsets)):

				nrd_dict = all_tables[inds[ind]][subsets[subset]][0]["# NRDs"]
				keys = ["Ref/Ref discordance", "Ref/Alt discordance", "Alt/Alt discordance","NRD"]
				this_nrd = nrd_dict["NRD"]
				nrds[ind][subset].append(this_nrd)



	for ind in range(len(inds)):
		plt.plot(covs, nrds[ind][0], linestyle="solid", color=ind_palette[ind], linewidth=1.5)
		plt.plot(covs, nrds[ind][1], linestyle="dashed", dashes =dashstyle, color=ind_palette[ind], linewidth=1.5)

	plt.xticks([0.1, 0.5, 1, 1.5,2])
	# plt.legend()
	# plt.title("Non-reference discordance of imputed")
	plt.ylim(ymax=18, ymin=0)
	plt.xlim(xmax=2, xmin=0.1)

	plt.plot([-1],[-1], linestyle="solid", color="black", label="overlapping", linewidth=1.5)
	plt.plot([-1],[-1], linestyle="dashed", dashes =dashstyle, color="black",label="non-overlapping", linewidth=1.5)
	plt.legend()

	plt.xlabel("Coverage (x)")
	plt.ylabel("NRD")

	# plt.savefig("NRDs_per_cov_"+refset+".svg", dpi=fdpi, bbox_inches="tight")
	# plt.savefig("NRDs_per_cov_"+refset+".pdf", dpi=fdpi, bbox_inches="tight")

	plt.show()
	plt.close()




#Plot concordance for 5 individuals for coverages  on x-axis
def  plot_concordance_per_cov(dir, refset, genos, ymin=0.9, resdir="", grid=False, axis=None):


	concordances=[[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]

	# plt.figure(figsize=(fwidth_med, fheight_med))
	# plt.figure(figsize=(fwidth_small, fheight_small))


	subsets = ["kept", "dropped"]
	covs = [0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75,2]
	majors = [0.0, 0.5, 1.0, 1.5, 2.0]
	minors = [0.25, 0.75, 1.25, 1.75]

	axis.set_xticks(majors)
	axis.set_xticks(minors, minor = True)

	axis.set_yticks([])
	axis.set_yticks([], minor = True)


	if grid:
		axis.grid(b=True, which='major', linewidth=lw2)
		axis.grid(b=True, which='minor', linewidth=lw3)

	for cov in covs:

		for ind in range(len(inds)):

			for subset in range(len(subsets)):
				concordance = get_concordance_from_gcss(dir,refset,cov,inds[ind],subsets[subset],genos)
				concordances[ind][subset].append(concordance)

	for ind in range(len(inds)):
		# if axis:
		# 	sns.lineplot(covs, concordances[ind][0], color=ind_palette[ind], linewidth=lw)
		# 	sns.lineplot(covs, concordances[ind][1], dashes=[(2,1)], color=ind_palette[ind], linewidth=lw)


		# else:
		plt.plot(covs, concordances[ind][0], linestyle="solid", color=ind_palette[ind], linewidth=lw, marker="o", markersize=1)#, markeredgecolor="black",markeredgewidth=0.09)
		plt.plot(covs, concordances[ind][1], linestyle="dashed", dashes =dashstyle, color=ind_palette[ind], linewidth=lw, marker="v",  markersize=1)#, markeredgecolor="black",markeredgewidth=0.09)

		print("MAX= " + str(max(concordances[ind][0])))


	# axis.set_xticks([0.0, 0.5, 1, 1.5,2])
	# plt.legend()
	# plt.title("Non-reference discordance of imputed")
	# axis.set_ylim(ymax=1.0, ymin=ymin)
	axis.set_xlim(xmax=2.0, xmin=0.0)

	if genos == ['RR Hom','RA Het','AA Hom']:
		genos_title = "all_genos"

	if genos == ['RA Het']:
		genos_title = "hets"


	# if if_overlap_legend:
	# 	# overlappingness legend
	# 	plt.plot([-1],[-1], linestyle="solid", color="black", label="overlapping", linewidth=lw)
	# 	plt.plot([-1],[-1], linestyle="dashed", dashes =dashstyle, color="black",label="non-overlapping", linewidth=lw)
	# 	plt.legend(loc=4)

	# if if_ind_legend:
	# 	# color legend
	# 	for ind in range(len(inds)):
	# 			plt.plot([-1],[-1], linestyle="solid", color=ind_palette[ind], label=inds[ind] ,linewidth=lw)
	# 	plt.legend()


	# if if_title:
	# 	plt.title("Concordance for genotypes : " +str(genos_title))


	# axis.xlabel("Coverage (x)")
	# plt.title("Concordance compared to high coverage for genotypes : " +str(genos))

	# plt.savefig(resdir+"conc_per_cov_"+refset+"_"+str(genos_title)+".pdf", dpi=fdpi, bbox_inches="tight")
	# plt.show()
	# plt.close()




#Plot concordance for 1 individuals for coverages  on x-axis
def  plot_concordance_per_cov_by_geno(dir, refset, ind, ymin=0.8):

	#concordances[geno][subset][cov]
	concordances=[[[],[]],[[],[]],[[],[]]]


	plt.figure(figsize=(fwidth_med, fheight_med))
	subsets = ["kept", "dropped"]
	genos = ['RR Hom','RA Het','AA Hom']
	covs = [0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75,2]


	for cov in covs:

		for geno in genos:

			for subset in range(len(subsets)):
				concordance = get_concordance_from_gcss(dir,refset,cov,ind,subsets[subset],[geno])
				concordances[genos.index(geno)][subset].append(concordance)

	for geno in range(len(genos)):
		plt.plot(covs, concordances[geno][0], linestyle="solid", color=geno_palette[geno], linewidth=lw)
		plt.plot(covs, concordances[geno][1], linestyle="dashed", dashes =dashstyle, color=geno_palette[geno], linewidth=lw)

	plt.xticks([0.1, 0.5, 1, 1.5,2])
	# plt.legend()
	# plt.title("Non-reference discordance of imputed")
	plt.ylim(ymax=1.0, ymin=ymin)
	plt.xlim(xmax=2.0, xmin=0.0)



	if if_overlap_legend:
		# overlappingness legend
		plt.plot([-1],[-1], linestyle="solid", color="black", label="overlapping", linewidth=lw)
		plt.plot([-1],[-1], linestyle="dashed", dashes =dashstyle, color="black",label="non-overlapping", linewidth=lw)
		plt.legend()

	if if_geno_legend:
		plt.plot([-1],[-1], linestyle="solid", color=geno_palette[0], label="RR", linewidth=lw)
		plt.plot([-1],[-1], linestyle="solid", color=geno_palette[1],label="RA", linewidth=lw)
		plt.plot([-1],[-1], linestyle="solid", color=geno_palette[2], label="AA", linewidth=lw)
		plt.legend()


	if if_title:
		plt.title("Concordance for : " + ind)



	plt.xlabel("Coverage (x)")
	plt.ylabel("Concordance")
	# plt.title("Concordance compared to high coverage for genotypes : " +str(genos))



	plt.savefig("conc_per_cov_"+refset+"_by_geno_"+ind+".pdf", dpi=fdpi, bbox_inches="tight")
	plt.show()
	plt.close()


#Plot concordance for 5 individuals for coverages  on x-axis
# from comparison to low coverage data
def  plot_lc_concordance_per_cov(dir, refset, genos):

	concordances=[[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]

	plt.figure(figsize=(fwidth_med, fheight_med))
	covs = [0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75,2]
	covs = [0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75]

	for cov in covs:

		for ind in range(len(inds)):

			concordance = get_concordance_from_gcss(dir,refset,cov,inds[ind],"all",genos,lc=True)
			concordances[ind][0].append(concordance)

	for ind in range(len(inds)):
		plt.plot(covs, concordances[ind][0], linestyle="solid", color=ind_palette[ind], linewidth=lw)

	plt.xticks([0.1, 0.5, 1, 1.5,2])
	# plt.legend()
	# plt.title("Non-reference discordance of imputed")
	plt.ylim(ymax=1.0, ymin=0.9)
	plt.xlim(xmax=2.0, xmin=0.0)


	plt.xlabel("Coverage (x)")
	plt.ylabel("Concordance")
	plt.title("Concordance compared to low coverage for genotypes : " +str(genos))
	plt.savefig("conc_per_cov_lc_"+refset+"_"+str(genos)+".pdf", dpi=fdpi, bbox_inches="tight")
	plt.show()
	plt.close()

#
# Get total entropy aggregated over MAF bins for genotypes genos
# genos e.g. [0,1,2]
def get_avg_entropy(e_table, genos):
	allele_freqs = np.array(e_table)[:,0]

	#aggregated over MAF bins for cross-checking againt reported GcSs
	entropy_agg= 0.0
	num_markers_agg = 0.0

	for a in range(len(allele_freqs)):

		entropy_this_bin = 0.0
		num_markers_this_bin = 0.0

		for geno in genos:

			entropy_this_geno = e_table[a][1+geno]
			num_markers_this_geno = e_table[a][4+geno]

			entropy_this_bin += entropy_this_geno*num_markers_this_geno
			num_markers_this_bin += num_markers_this_geno


		average_entropy_this_bin = entropy_this_bin / num_markers_this_bin
		# print("---af: "+str(allele_freqs[a]) + " num markers= " + str(num_markers_this_bin) + "avg entropy= " + str(average_entropy_this_bin))

		entropy_agg += average_entropy_this_bin * num_markers_this_bin
		num_markers_agg += num_markers_this_bin


	avg_entropy_agg = float(entropy_agg / num_markers_agg)

	return avg_entropy_agg



#Plot entropy for 5 individuals for coverages  on x-axis
def  plot_entropy_per_cov(datadir, statsdir, refset, genos, ymax=0.6):

	entropies=[[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]

	plt.figure(figsize=(fwidth_med, fheight_med))
	subsets = ["kept", "dropped"]
	covs = [0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75,2]

	for cov in covs:

		for ind in range(len(inds)):

			for subset in range(len(subsets)):
				res_file = datadir+ statsdir + "/entropy_HCimp"+impset+"_"+str(cov)+"x_"+subsets[subset]+"_"+refset+"_"+inds[ind]+".output"
				e_table,tot_entropy = load_entropy_results(res_file)
				avg_entropy = get_avg_entropy(e_table, genos)
				print(str(cov) +" " + str(inds[ind]) + " " + subsets[subset] +" " + str(genos) +" estimated = " + str(avg_entropy) + " counted = " + str(tot_entropy))
				if len(genos) == 3:
					entropies[ind][subset].append(tot_entropy)
				else:
					entropies[ind][subset].append(avg_entropy)

	for ind in range(len(inds)):
		plt.plot(covs, entropies[ind][0], linestyle="solid", color=ind_palette[ind], linewidth=lw)
		plt.plot(covs, entropies[ind][1], linestyle="dashed", dashes =dashstyle, color=ind_palette[ind], linewidth=lw)

	plt.xticks([0.1, 0.5, 1, 1.5,2])
	# plt.legend()
	# plt.title("Non-reference discordance of imputed")
	plt.ylim(ymax=ymax, ymin=0.0)
	plt.xlim(xmax=2.0, xmin=0.0)

	plt.plot([-1],[-1], linestyle="solid", color="black", label="overlapping", linewidth=lw)
	plt.plot([-1],[-1], linestyle="dashed", dashes =dashstyle, color="black",label="non-overlapping", linewidth=lw)
	plt.legend()

	plt.xlabel("Coverage (x)")
	plt.ylabel("Entropy")
	plt.title("Entropy compared to high coverage for genotypes : " + str(genos))
	plt.savefig("entr_per_cov_"+refset+"_"+str(genos)+".pdf", dpi=fdpi, bbox_inches="tight")
	plt.show()
	plt.close()



# Plot entropy for 5 individuals for coverages  on x-axis
# For comparison to low cov data
def  plot_lc_entropy_per_cov(datadir, statsdir, refset, genos, type="entropy"):

	entropies=[[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]

	plt.figure(figsize=(fwidth_med, fheight_med))
	covs = [0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75,2]

	for cov in covs:

		for ind in range(len(inds)):

			res_file = datadir+ statsdir + "/"+type+"_lcimp"+impset+"_"+str(cov)+"x_"+refset+"_"+inds[ind]+".output"
			e_table,tot_entropy = load_entropy_results(res_file)
			avg_entropy = get_avg_entropy(e_table, genos)
			print(str(cov) +" " + str(inds[ind]) + " LC " + str(genos) +" estimated = " + str(avg_entropy) + " counted = " + str(tot_entropy))
			if len(genos) == 3:
				entropies[ind][0].append(tot_entropy)
			else:
				entropies[ind][0].append(avg_entropy)

	for ind in range(len(inds)):
		plt.plot(covs, entropies[ind][0], linestyle="solid", color=ind_palette[ind], linewidth=lw)

	plt.xticks([0.1, 0.5, 1, 1.5,2])

	plt.ylim(ymin=0.0)
	plt.xlim(xmax=2.0, xmin=0.0)


	plt.xlabel("Coverage (x)")
	plt.ylabel("Entropy")
	plt.title("Entropy compared to low coverage for genotypes : " + str(genos))
	plt.savefig("entr_per_cov_lc_"+refset+"_"+str(genos)+".pdf", dpi=fdpi, bbox_inches="tight")
	plt.show()
	plt.close()


#Plot error (discordance, GCsAF) for 5 individuals for MAF  on x-axis
#
def  plot_discordance_per_maf(dir, refset, cov):

	# errors[ind][0][c] : discordance for ind for KEPT markers at maf c
	# errors[ind][1][c] : discordance for ind for DROPPED markers at maf c

	errors=[[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]

	plt.figure(figsize=(fwidth_med, fheight_med))
	subsets = ["kept", "dropped"]



	print("Reading " + dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json')
	with open(dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json', 'r') as fp:
		all_tables = json.load(fp)

	for ind in range(len(inds)):

		for subset in range(len(subsets)):
			print("keys: " + str(all_tables[inds[ind]][subsets[subset]][3].keys()))
			error_dict = all_tables[inds[ind]][subsets[subset]][3]["# GCsAF"]

			allele_freqs = map(float,error_dict['allele frequency'])
			this_row = []
			ALL_correct = 0
			ALL = 0
			for a in range(len(allele_freqs)):

				num_AA_match = float(error_dict['AA Hom matches'][a])
				num_AA_mismatch = float(error_dict['AA Hom mismatches'][a])

				num_RA_match = float(error_dict['RA Het matches'][a])
				num_RA_mismatch = float(error_dict['RA Het mismatches'][a])

				num_RR_match = float(error_dict['RR Hom matches'][a])
				num_RR_mismatch = float(error_dict['RR Hom mismatches'][a])

				discordance = (num_AA_mismatch + num_RA_mismatch + num_RR_mismatch) / float(num_AA_match + num_RA_match + num_RR_match +num_AA_mismatch + num_RA_mismatch + num_RR_mismatch)
				# concordance = (num_AA_match + num_RA_match + num_RR_match) / float(num_AA_match + num_RA_match + num_RR_match +num_AA_mismatch + num_RA_mismatch + num_RR_mismatch)

				print("Allele frequency: " +str(allele_freqs[a]) + "  discordance = " + str(discordance))
				this_row.append(str(discordance))

				ALL_correct +=  num_AA_match + num_RA_match + num_RR_match
				ALL += num_AA_match + num_RA_match + num_RR_match +num_AA_mismatch + num_RA_mismatch + num_RR_mismatch

			concordance = ALL_correct / float(ALL)

			print(inds[ind] + " , " + subsets[subset] + " : concordance = " + str(concordance))
			errors[ind][subset] = this_row



	for ind in range(len(inds)):
		plt.plot(allele_freqs, errors[ind][0], linestyle="solid", color=ind_palette[ind], linewidth=1.5)
		plt.plot(allele_freqs, errors[ind][1], linestyle="dashed", dashes =dashstyle, color=ind_palette[ind], linewidth=1.5)

	plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
	# plt.legend()
	# plt.title("Non-reference discordance of imputed")
	plt.ylim(ymin=0.5, ymax=1.0)
	plt.xlim(xmax=0.5, xmin=0.0)

	plt.plot([-1],[-1], linestyle="solid", color="black", label="overlapping", linewidth=1.5)
	plt.plot([-1],[-1], linestyle="dashed", dashes =dashstyle, color="black",label="non-overlapping", linewidth=1.5)
	plt.legend()
	plt.title(cov)
	plt.xlabel("MAF")
	plt.ylabel("Genotype discordance")

	# plt.savefig("NRDs_per_cov_"+refset+".svg", dpi=fdpi, bbox_inches="tight")
	# plt.savefig("NRDs_per_cov_"+refset+".pdf", dpi=fdpi, bbox_inches="tight")

	plt.show()
	plt.close()



#Plot error (discordance, GCsAF) for 5 individuals for MAF  on x-axis
#
def  plot_concordance_per_maf_by_geno(dir, refset, cov, ind, ymin=0.5):

	# errors[genotype][0][c] : discordance for genotype for KEPT markers at maf c
	# errors[genotype][1][c] : discordance for genotype for DROPPED markers at maf c

	errors=[[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]


	plt.figure(figsize=(fwidth_med, fheight_med))
	subsets = ["kept", "dropped"]



	print("Reading " + dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json')
	with open(dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json', 'r') as fp:
		all_tables = json.load(fp)


	for subset in range(len(subsets)):
		print("keys: " + str(all_tables[ind][subsets[subset]][3].keys()))
		error_dict = all_tables[ind][subsets[subset]][3]["# GCsAF"]

		allele_freqs = map(float,error_dict['allele frequency'])
		this_AA = []
		this_RA = []
		this_RR = []

		for a in range(len(allele_freqs)):

			num_AA_match = float(error_dict['AA Hom matches'][a])
			num_AA_mismatch = float(error_dict['AA Hom mismatches'][a])
			AA_concordance = float(num_AA_match / (num_AA_match+num_AA_mismatch))

			num_RA_match = float(error_dict['RA Het matches'][a])
			num_RA_mismatch = float(error_dict['RA Het mismatches'][a])
			RA_concordance = float(num_RA_match / (num_RA_match+num_RA_mismatch))

			num_RR_match = float(error_dict['RR Hom matches'][a])
			num_RR_mismatch = float(error_dict['RR Hom mismatches'][a])
			RR_concordance = float(num_RR_match / (num_RR_match+num_RR_mismatch))

			this_AA.append(AA_concordance)
			this_RA.append(RA_concordance)
			this_RR.append(RR_concordance)


		errors[0][subset] = this_RR
		errors[1][subset] = this_RA
		errors[2][subset] = this_AA


	for geno in range(3):
		plt.plot(allele_freqs, errors[geno][0], linestyle="solid", color=geno_palette[geno], linewidth=1.5)
		plt.plot(allele_freqs, errors[geno][1], linestyle="dashed", dashes =dashstyle, color=geno_palette[geno], linewidth=1.5)

	plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
	# plt.legend()
	# plt.title("Non-reference discordance of imputed")
	plt.ylim(ymin=ymin, ymax=1.0)
	plt.xlim(xmin=0.0, xmax=0.5)

	plt.plot([-1],[-1], linestyle="solid", color="black", label="overlapping", linewidth=1.5)
	plt.plot([-1],[-1], linestyle="dashed", dashes =dashstyle, color="black",label="non-overlapping", linewidth=1.5)

	plt.legend()
	plt.title(ind + "   " + str(cov))
	plt.xlabel("MAF")
	plt.ylabel("Genotype concordance")

	plt.savefig("conc_per_MAF_by_geno/"+dir+"_"+refset+"/conc_per_MAF_by_geno_"+refset+"_"+str(cov)+"_"+ind +".pdf", dpi=fdpi, bbox_inches="tight")

	# plt.show()
	plt.close()




#Plot error (discordance, GCsAF) for 5 individuals for MAF  on x-axis
#
def  plot_entropy_per_maf_by_geno(datadir,statsdir, refset, cov, ind, ymax=12):

	# errors[genotype][0][c] : discordance for genotype for KEPT markers at maf c
	# errors[genotype][1][c] : discordance for genotype for DROPPED markers at maf c

	impset = statsdir.split("imp")[-1]

	errors=[[[],[]],[[],[]],[[],[]]]


	plt.figure(figsize=(fwidth_med, fheight_med))
	subsets = ["kept", "dropped"]


	for subset in range(len(subsets)):


		res_file = datadir+ statsdir + "/entropy_HCimp"+impset+"_"+str(cov)+"x_"+subsets[subset]+"_"+refset+"_"+ind+".output"
		e_table,tot_entropy = load_entropy_results(res_file)

		allele_freqs = np.array(e_table)[:,0]

		this_conc = []

		#aggregated over MAF bins for cross-checking
		entropy_agg= 0.0
		num_markers_agg = 0.0

		this_RR = []
		this_RA = []
		this_AA = []

		for a in range(len(allele_freqs)):

			entropy_this_bin = 0.0
			num_markers_this_bin = 0.0

			for geno in [0,1,2]:

				entropy_this_geno = e_table[a][1+geno]
				num_markers_this_geno = e_table[a][4+geno]

				errors[geno][subset].append(entropy_this_geno)

				entropy_this_bin += entropy_this_geno*num_markers_this_geno
				num_markers_this_bin += num_markers_this_geno


			average_entropy_this_bin = entropy_this_bin / num_markers_this_bin
			print("---af: "+str(allele_freqs[a]) + " num markers= " + str(num_markers_this_bin) + "avg entropy= " + str(average_entropy_this_bin))

			this_conc.append(average_entropy_this_bin)

			entropy_agg += average_entropy_this_bin * num_markers_this_bin

	for geno in range(3):
		plt.plot(allele_freqs, errors[geno][0], linestyle="solid", color=geno_palette[geno], linewidth=1.5)
		plt.plot(allele_freqs, errors[geno][1], linestyle="dashed", dashes =dashstyle, color=geno_palette[geno], linewidth=1.5)

	plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
	# plt.legend()
	# plt.title("Non-reference discordance of imputed")
	plt.ylim(ymin=0.0, ymax=ymax)
	plt.xlim(xmin=0.0, xmax=0.5)

	plt.plot([-1],[-1], linestyle="solid", color="black", label="overlapping", linewidth=1.5)
	plt.plot([-1],[-1], linestyle="dashed", dashes =dashstyle, color="black",label="non-overlapping", linewidth=1.5)
	plt.legend()
	plt.title(ind + "   " + str(cov))
	plt.xlabel("MAF")
	plt.ylabel("Entropy")

	plt.savefig("entr_per_MAF_by_geno/"+statsdir+"_"+refset+"/entr_per_MAF_by_geno_"+refset+"_"+str(cov)+"_"+ind +".pdf", dpi=fdpi, bbox_inches="tight")

	# plt.show()
	plt.close()


def rebin(a, shape):
	return np.convolve(a, np.ones((shape,))/shape, mode='valid')


#Plot concordance for 5 individuals for MAF  on x-axis for genotypes genos
# genos is a list that may contain at least 1 of: 'RR Hom' 'Ra Het' 'AA Hom'
#
def plot_concordance_per_maf_by_ind(datadir, refset, cov, genos, ymin=0.95, ymax=1.0, y2min=0.6, ref_title=None, metric="concordance", retained=False,
									rebin_shape=None, resdir="", show_y=True, grid=True, make_figure=True, input_ax=None, tot_markers=False):

	# errors[ind][0][c] : discordance for ind for KEPT markers at maf c
	# errors[ind][1][c] : discordance for ind for DROPPED markers at maf c

	if ref_title == None:
		ref_title = ref

	errors=[[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]

	if make_figure:
		fig, ax1 = plt.subplots(figsize=(fwidth_small, fheight_small))

	else:
		ax1 = input_ax

	if retained:
		ax2 = ax1.twinx()
		ax2.set_ylim(ymin=y2min, ymax=0.9)

	if tot_markers:
		ax2 = ax1.twinx()
		# ax2.set_ylim(ymin=y2min, ymax=0.9)


	# plt.rc('axes.formatter', useoffset=False)

	subsets = ["kept", "dropped"]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]



	print("GCsAF: reading " + datadir + '_' + refset + '/all_stats_' + str(cov) + '_' + refset + '.json')
	with open(datadir+ '_'+refset+ '/all_stats_'+str(cov)+ '_'+refset+ '.json', 'r') as fp:
		all_tables = json.load(fp)

	for ind in range(len(inds)):
		for subset in range(len(subsets)):

			marker_counts_by_af_bin=[]


			error_dict = all_tables[inds[ind]][subsets[subset]][3]["# GCsAF"]
			allele_freqs = map(float,error_dict['allele frequency'])[:-1]

			this_conc = []

			#aggregated over MAF bins for cross-checking againt reported GcSs
			num_match_agg = 0
			num_mismatch_agg = 0

			print("Allele frequency mbin marker counts for " + inds[ind] + " " + subsets[subset])

			for a in range(len(allele_freqs)):

				num_match = 0
				num_mismatch = 0

				for geno in genos:
					num_match += float(error_dict[geno+' matches'][a])
					num_mismatch += float(error_dict[geno+' mismatches'][a])

				print(str(num_match+num_mismatch))
				marker_counts_by_af_bin.append(num_match+num_mismatch)


				concordance = float(num_match / (num_match+num_mismatch))
				this_conc.append(concordance)

				num_match_agg += num_match
				num_mismatch_agg += num_mismatch


			with open("marker_counts/"+ref_title+"_"+inds[ind] + "_" + subsets[subset]+"_marker_counts_by_af", "w") as mfile:
				for c in marker_counts_by_af_bin:
						mfile.write(str(c))
						mfile.write('\n')

			concordance_agg = float(num_match_agg / (num_match_agg+num_mismatch_agg))
			concordance_gcss = get_concordance_from_gcss(datadir, refset, cov, inds[ind], subsets[subset], genos)

			# print("Total num markers = " + str(num_match_agg+num_mismatch_agg) + "   concordance = " + str(concordance_agg) + "  gcss concordance = " + str(concordance_gcss))

			print(inds[ind] + " " + subsets[subset] + " " + str(cov) + " " + str(genos))
			print("GCsAF concordance = " + str(concordance_agg))
			print("---- match: " + str(num_match_agg))
			print("---- mismatch: " + str(num_mismatch_agg))
			print("---- total: " + str(num_match_agg + num_mismatch_agg))


			assert abs(concordance_gcss-concordance_agg) < 0.0005

			errors[ind][subset] = this_conc

		if metric == "discordance":
			# turn concordance into discordance
			errors[ind][0] = map(lambda x : 1 - x , errors[ind][0])
			errors[ind][1] = map(lambda x : 1 - x , errors[ind][1])

			if inds[ind] == "sf12":
				print("FIRST DISCORDNACES SF12 before rebin" )
				print(errors[ind][0][0:5])
				print(errors[ind][1][0:5])

	print("Allele freqs: " + str(allele_freqs))

	ax1.set_xticks([0.001, 0.1])

	if rebin_shape:
		allele_freqs = rebin(np.array(allele_freqs), rebin_shape)
		print("before rebin ")
		print(np.array(errors).shape)
		for ind in range(len(inds)):
			errors[ind][0] = rebin(np.array(errors[ind][0]), rebin_shape)
		for ind in range(len(inds)):
			errors[ind][1] = rebin(np.array(errors[ind][1]), rebin_shape)

		print("after rebin")
		print(np.array(errors).shape)
		# purely for plotting purposes: replacing 0.015 with 0.01 to get log axis to look comprehensible
		allele_freqs[0] = 0.01

		print("FIRST DISCORDNACES SF12 AFTER rebin" )
		print(errors[1][0][0:5])
		print(errors[1][1][0:5])


	print("Allele freqs: " + str(allele_freqs))
	ax1.set_xlim(xmin=min(allele_freqs), xmax = max(allele_freqs))




	# ax1.plot([-1],[-1], linestyle="solid", color="black", label="overlapping", linewidth=lw)
	# ax1.plot([-1],[-1], linestyle="dashed", dashes =dashstyle, color="black",label="non-overlapping", linewidth=lw)


	if genos == ['RR Hom','RA Het','AA Hom']:
		genos_title = "all_genos"

	if genos == ['RA Het']:
		genos_title = "hets"

	if rebin_shape:
		genos_title += "_rebin_"+str(rebin_shape)

	if show_y:
		ax1.set_ylabel("Genotype " + metric)
		plt.setp(ax1.get_yticklabels(), visible=True)
		ax1.tick_params(axis='x', which='major', pad=3)

	else:
		plt.setp(ax1.get_yticklabels(), visible=False)


	ax1.set_xlabel("MAF")


	################### add percentage markers kept ############

	if retained:
		# get this datas marker counts
		this_marker_counts = []
		for ind in inds:
			this_ind_marker_counts = [0 for c in marker_counts_by_af_bin]
			for subset in subsets:
				print("marker_counts/"+ref_title+"_"+ind + "_" + subset+"_marker_counts_by_af")

				with open("marker_counts/"+ref_title+"_"+ind + "_" + subset+"_marker_counts_by_af", "r") as thisfile:
					for c in range(len(marker_counts_by_af_bin)):
						marker_count = float(thisfile.readline()[:-1])
						this_ind_marker_counts[c] += marker_count
						# print("..................."+str(marker_count))


			this_marker_counts.append(this_ind_marker_counts)

		print("this marker counts: ")
		print(np.array(this_marker_counts))

		# get corresponding unfiltered marker counts
		unfiltered_marker_counts = []
		for ind in inds:
			this_ind_marker_counts = [0 for c in marker_counts_by_af_bin]
			for subset in subsets:
				unfiltered_title = ref_title.split('_')[0]
				print("marker_counts/"+unfiltered_title+"_"+ind + "_" + subset+"_marker_counts_by_af")
				with open("marker_counts/"+unfiltered_title+"_"+ind + "_" + subset+"_marker_counts_by_af", "r") as thisfile:
					for c in range(len(marker_counts_by_af_bin)):
						marker_count = float(thisfile.readline()[:-1])
						this_ind_marker_counts[c] += marker_count
						# print("..................."+str(marker_count))


			unfiltered_marker_counts.append(this_ind_marker_counts)

		print("unfiltered marker counts: ")
		print(np.array(unfiltered_marker_counts))



		#########################################################

		diffs = list(np.divide(this_marker_counts, unfiltered_marker_counts))

		if rebin_shape:
			for ind in range(len(inds)):
				diffs[ind]=rebin(diffs[ind], rebin_shape)

		diffs=np.array(diffs)
		print(diffs.shape)

		diffs_avg = np.average(diffs, axis=0)
		diffs_std = np.std(diffs, axis=0)

		print("DIFFS AVG")
		print(diffs_avg)
		print("DIFFS STD")
		print(diffs_std)
		print("DIFFS MIN")
		print(np.min(diffs, axis=0))
		print("DIFFS MAX")
		print(np.max(diffs, axis=0))


		ax2.plot(allele_freqs, diffs_avg, color=retained_col, linewidth=lw4, linestyle="solid", alpha=0.6)
		# ax2.fill_between(allele_freqs, diffs_avg+diffs_std, diffs_avg-diffs_std, facecolor=retained_col, alpha=0.1)
		ax2.fill_between(allele_freqs, np.max(diffs, axis=0), np.min(diffs, axis=0), facecolor=retained_col, alpha=0.07)


		# ax2.tick_params(axis='y', colors=retained_col)


		# for ind in range(len(inds)):
		# 	# ax2.plot(allele_freqs, diffs[ind], color="#e6e6e6", linewidth=lw)
		# 	# ax2.plot(allele_freqs, diffs[ind], color=ind_palette[ind], linewidth=lw, linestyle="--",dashes=(2,5))
		#
		# 	# ax2.plot(allele_freqs, diffs[ind], color="gray", linewidth=lw/2)
		# 	# ax2.plot(allele_freqs, diffs[ind], color=ind_palette[ind], linewidth=lw, linestyle="--", dashes=(3,5), alpha=1.0)
		#
		# 	ax2.plot(allele_freqs, diffs[ind], color=ind_palette[ind], linewidth=lw*1.4, linestyle="solid")
		# 	ax2.plot(allele_freqs, diffs[ind], color="#e6e6e6", linewidth=lw*0.7)
		#




		plt.setp(ax2.get_yticklabels(), visible=True)

		ax2.set_ylabel("Retained markers", rotation=270, labelpad=10)#, color=retained_col)
		ax2.set_yticks([0.6, 0.7, 0.8, 0.9])


	if tot_markers:
		# get this datas marker counts
		this_marker_counts = []
		for ind in inds:
			this_ind_marker_counts = [0 for c in marker_counts_by_af_bin]
			for subset in subsets:
				print("marker_counts/"+ref_title+"_"+ind + "_" + subset+"_marker_counts_by_af")

				with open("marker_counts/"+ref_title+"_"+ind + "_" + subset+"_marker_counts_by_af", "r") as thisfile:
					for c in range(len(marker_counts_by_af_bin)):
						marker_count = float(thisfile.readline()[:-1])
						this_ind_marker_counts[c] += marker_count
						# print("..................."+str(marker_count))


			this_marker_counts.append(this_ind_marker_counts)


		if rebin_shape:
			for ind in range(len(inds)):
				this_marker_counts[ind]=rebin(this_marker_counts[ind], rebin_shape)

		counts_avg = np.average(this_marker_counts, axis=0)



		ax2.plot(allele_freqs, counts_avg, color=retained_col, linewidth=lw4, linestyle="solid", alpha=0.6)
		# ax2.fill_between(allele_freqs, diffs_avg+diffs_std, diffs_avg-diffs_std, facecolor=retained_col, alpha=0.1)
		ax2.fill_between(allele_freqs, np.max(this_marker_counts, axis=0), np.min(this_marker_counts, axis=0), facecolor=retained_col, alpha=0.07)

		plt.setp(ax2.get_yticklabels(), visible=True)
		ax2.set_ylabel("Total markers", rotation=270, labelpad=10)#, color=retained_col)
		# ax2.set_yticks([0.6, 0.7, 0.8, 0.9])





	if metric == "discordance":
		for ind in range(len(inds)):
			ax1.loglog(allele_freqs, errors[ind][0], linestyle="solid", color=ind_palette[ind], linewidth=lw)
			ax1.loglog(allele_freqs, errors[ind][1], linestyle="dashed", dashes =dashstyle, color=ind_palette[ind], linewidth=lw)

		ax1.set_ylim(ymin=ymin, ymax=ymax)
		if ymin == 0.001:
			ax1.set_yticks([0.001,0.01, 0.1])

		if ymin == 0.01:
			ax1.set_yticks([0.01, 0.1])






		# ax1.set_ylim(ymin=0.0, ymax=0.5)

	if metric == "concordance":
		print("Plotting concordance : ")
		for ind in range(len(inds)):
			print(errors[ind][0])
			print(errors[ind][1])

			ax1.semilogx(allele_freqs, errors[ind][0], linestyle="solid", color=ind_palette[ind], linewidth=lw)
			ax1.semilogx(allele_freqs, errors[ind][1], linestyle="dashed", dashes =dashstyle, color=ind_palette[ind], linewidth=lw)

		ax1.set_ylim(ymin=ymin, ymax=ymax)



	# ax1.xaxis.set_major_formatter(mticker.ScalarFormatter())
	# ax1.xaxis.get_major_formatter().set_scientific(True)
	# ax1.xaxis.get_major_formatter().set_useOffset(False)
	# ax1.xaxis.set_minor_formatter(mticker.ScalarFormatter())


	# ax1.yaxis.set_major_formatter(mticker.ScalarFormatter())
	# ax1.yaxis.get_major_formatter().set_scientific(False)
	# ax1.yaxis.get_major_formatter().set_useOffset(False)
	# ax1.yaxis.set_minor_formatter(mticker.ScalarFormatter())
	#

	# ax1.get_yaxis().set_major_locator(mticker.AutoLocator())
	# ax1.get_yaxis().set_minor_locator(mticker.AutoMinorLocator())

	# plt.rcParams['ytick.major.pad'] = 0


	if grid:
		print("grid?")
		ax1.grid(b=True, which='major', linewidth=lw2)
		ax1.grid(b=True, which='minor', linewidth=lw3, linestyle="--")

	if retained:
		ax2.grid(b=False, which='major')
		plt.setp(ax2.get_yticklabels(), visible=True)
		plt.setp(ax1.get_yticklabels(), visible=True)

	# else:
	# 	ax1.grid(b=True, which='major', linewidth=0.5)


	# plt.yticks([0.02, 0.1])

	# plt.yticks([0.1])


	if if_ind_legend:
		# color legend
		for ind in range(len(inds)):
				plt.plot([-1],[-1], linestyle="solid", color=ind_palette[ind], label=inds[ind] ,linewidth=lw)
		plt.legend(loc=4)


	if if_title:
		plt.title("genos: " + str(genos_title) + "   cov: " + str(cov) + "   ref: " + ref_title)


	# plt.rcParams['xtick.major.size'] = 0
	# plt.rcParams['ytick.major.size'] = 0
	# plt.rcParams['xtick.minor.size'] = 0
	# plt.rcParams['ytick.minor.size'] = 0

	ax1.tick_params(axis=u'both', which=u'both',length=0)

	if retained:
		ax2.tick_params(axis=u'both', which=u'both',length=0)

	plt.rcParams['xtick.major.pad'] =  3
	plt.rcParams['ytick.major.pad'] =  2


	if make_figure:

		plt.tight_layout(pad=0)
		plt.savefig(resdir+metric+"_per_MAF_by_ind_"+ref_title+"_"+str(cov)+"_"+str(genos_title)+".pdf", dpi=fdpi)


		# plt.savefig(resdir+metric+"_per_MAF_by_ind_"+ref_title+"_"+str(cov)+"_"+str(genos_title)+".pdf", dpi=fdpi, bbox_inches="tight")

		# plt.show()
		plt.close()




def load_entropy_results(res_file):
	print("Reading from " + str(res_file))
	with open(res_file,'r') as resfile:
		true_file = resfile.readline()
		comp_file = resfile.readline()
		af_bins = resfile.readline()
		af_bins = resfile.readline()
		af_bins = resfile.readline()
		af_bins = resfile.readline()
		af_bins = resfile.readline()

		header = resfile.readline()

		entropy_table = []

		for line in resfile:
			if line.startswith("#"):
				break
			t =  map(float,line.split(','))
			entropy_table.append(t)

		tot_line = resfile.next()
		total_entropy = float(tot_line.split(' ')[-1])
	print(entropy_table)
	return entropy_table,total_entropy


def load_compref_results(res_file):
	print("Reading from " + str(res_file))
	with open(res_file,'r') as resfile:
		true_file = resfile.readline()
		comp_file = resfile.readline()
		af_bins = resfile.readline()
		af_bins = resfile.readline()
		af_bins = resfile.readline()
		af_bins = resfile.readline()
		af_bins = resfile.readline()

		header = resfile.readline()



		comp_table = []

		for line in resfile:
			if line.startswith("#"):
				break
			t =  map(float,line.split(','))
			comp_table.append(t)

		all_same_line = resfile.next()
		all_same = float(all_same_line.split(' ')[-1])

		all_markers_line = resfile.next()
		all_markers = float(all_markers_line.split(' ')[-1])
		total_refcomp = all_same/all_markers

	# print("Comp Table: ")
	# print(comp_table)
	# print("Total refcomp: " + str(total_refcomp))
	return comp_table,total_refcomp


# a,b = load_entropy_results('/home/kristiina/reslink/concordance/stats_rmts_61_5_impnofilter_anno/entropy_HCimp_anno_1x_dropped_FULL_ans17.output')
# print a
# print b
# exit()



#Plot concordance for 5 individuals for MAF  on x-axis for genotypes genos
# genos is a list that may contain at least 1 of: 'RR Hom' 'Ra Het' 'AA Hom'
#
def  plot_lc_concordance_per_maf_by_ind(dir, refset, cov, genos, ymin=0.95):

	# errors[ind][0][c] : discordance for ind for markers at maf c

	errors=[[[]],[[]],[[]],[[]],[[]]]


	plt.figure(figsize=(fwidth_med, fheight_med))

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]


	print("GCsAF: reading " + dir+'_'+refset+'/all_stats_lc_'+str(cov)+'_'+refset+'.json')
	with open(dir+'_'+refset+'/all_stats_lc_'+str(cov)+'_'+refset+'.json', 'r') as fp:
		all_tables = json.load(fp)

	for ind in range(len(inds)):
		# print("keys: " + str(all_tables[inds[ind]][subsets[subset]][3].keys()))
		print(all_tables[inds[ind]].keys())
		error_dict = all_tables[inds[ind]]["all"][3]["# GCsAF"]

		allele_freqs = map(float,error_dict['allele frequency'])[:-1]
		this_conc = []

		#aggregated over MAF bins for cross-checking againt reported GcSs
		num_match_agg = 0
		num_mismatch_agg = 0

		for a in range(len(allele_freqs)):

			num_match = 0
			num_mismatch = 0

			for geno in genos:
				num_match += float(error_dict[geno+' matches'][a])
				num_mismatch += float(error_dict[geno+' mismatches'][a])

			print("---af: "+str(allele_freqs[a]) + " num markers= " + str(float(error_dict[geno+' matches'][a])+float(error_dict[geno+' mismatches'][a])))
			concordance = float(num_match / (num_match+num_mismatch))
			this_conc.append(concordance)

			num_match_agg += num_match
			num_mismatch_agg += num_mismatch


		concordance_agg = float(num_match_agg / (num_match_agg+num_mismatch_agg))
		# concordance_gcss = get_concordance_from_gcss(dir, refset,cov,inds[ind],subsets[subset],genos)

		# print("Total num markers = " + str(num_match_agg+num_mismatch_agg) + "   concordance = " + str(concordance_agg) + "  gcss concordance = " + str(concordance_gcss))

		print(inds[ind] + " " + " " + str(cov) + " " + str(genos))
		print("GCsAF concordance = " + str(concordance_agg))
		print("---- match: " + str(num_match_agg))
		print("---- mismatch: " + str(num_mismatch_agg))
		print("---- total: " + str(num_match_agg + num_mismatch_agg))


		# assert abs(concordance_gcss-concordance_agg) < 0.0005

		errors[ind][0] = this_conc

	for ind in range(len(inds)):
		plt.plot(allele_freqs, errors[ind][0], linestyle="solid", color=ind_palette[ind], linewidth=lw)
		# plt.plot(allele_freqs, errors[ind][1], linestyle="dashed", dashes =dashstyle, color=ind_palette[ind], linewidth=lw)
		print(allele_freqs)
	# plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
	# plt.legend()
	# plt.title("Non-reference discordance of imputed")
	plt.ylim(ymin=ymin, ymax=1.0)
	plt.xlim(xmin=0.0, xmax=0.5)

	plt.plot([-1],[-1], linestyle="solid", color="black", label="overlapping", linewidth=lw)
	# plt.plot([-1],[-1], linestyle="dashed", dashes =dashstyle, color="black",label="non-overlapping", linewidth=lw)
	plt.legend()
	plt.title(str(genos) + "   " + str(cov))
	plt.xlabel("MAF")
	plt.ylabel("Genotype concordance")


	plt.savefig(dir+'_'+refset+"/conc_per_MAF_by_ind/conc_lc_per_MAF_by_ind_"+refset+"_"+str(cov)+"_"+str(genos) +".pdf", dpi=fdpi, bbox_inches="tight")

	plt.show()
	plt.close()

# Plot entropy for 5 individuals for MAF  on x-axis for genotypes genos
# genos is a list that may contain at least 1 of: 0,1,2
#
def  plot_entropy_per_maf_by_ind(datadir, statsdir, refset, cov, genos, ymax):

	# errors[ind][0][c] : entropy for ind for KEPT markers at maf c
	# errors[ind][1][c] : entropy for ind for DROPPED markers at maf c

	errors=[[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]


	plt.figure(figsize=(fwidth_med, fheight_med))
	subsets = ["kept", "dropped"]
	# subsets = ["dropped"]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]
	# inds = ["ans17"]


	for ind in range(len(inds)):
		for subset in range(len(subsets)):

			# res_file="~/Projects/PythonUtils/General/test.output"
			res_file = datadir + statsdir + "/entropy_HCimp"+impset+"_"+str(cov)+"x_"+subsets[subset]+"_"+refset+"_"+inds[ind]+".output"

			e_table,tot_entropy = load_entropy_results(res_file)

			allele_freqs = np.array(e_table)[:,0]


			this_conc = []

			#aggregated over MAF bins for cross-checking againt reported GcSs
			entropy_agg= 0.0
			num_markers_agg = 0.0

			for a in range(len(allele_freqs)):

				entropy_this_bin = 0.0
				num_markers_this_bin = 0.0

				for geno in genos:

					entropy_this_geno = e_table[a][1+geno]
					num_markers_this_geno = e_table[a][4+geno]

					entropy_this_bin += entropy_this_geno*num_markers_this_geno
					num_markers_this_bin += num_markers_this_geno


				average_entropy_this_bin = entropy_this_bin / num_markers_this_bin
				print("---af: "+str(allele_freqs[a]) + " num markers= " + str(num_markers_this_bin) + "avg entropy= " + str(average_entropy_this_bin))

				this_conc.append(average_entropy_this_bin)

				entropy_agg += average_entropy_this_bin * num_markers_this_bin
				num_markers_agg += num_markers_this_bin


			avg_entropy_agg = float(entropy_agg / num_markers_agg)

			# print("Total num markers = " + str(_+num_markers_agg) + "   concordance = " + str(avg_entropy_agg) + "  gcss concordance = " + str(concordance_gcss))

			print(inds[ind] + " " + subsets[subset] + " " + str(cov) + " " + str(genos))
			print("Estimated entropy " + str(avg_entropy_agg))
			print("Counted entropy " + str(tot_entropy))


			errors[ind][subset] = this_conc

	for ind in range(len(inds)):
		plt.plot(allele_freqs, errors[ind][0], linestyle="solid", color=ind_palette[ind], linewidth=lw)
		plt.plot(allele_freqs, errors[ind][1], linestyle="dashed", dashes =dashstyle, color=ind_palette[ind], linewidth=lw)
		print(allele_freqs)
	plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
	# plt.legend()
	# plt.title("Non-reference discordance of imputed")
	plt.ylim(ymin=0.0, ymax=ymax)
	plt.xlim(xmin=0.0, xmax=0.5)

	plt.plot([-1],[-1], linestyle="solid", color="black", label="overlapping", linewidth=lw)
	plt.plot([-1],[-1], linestyle="dashed", dashes =dashstyle, color="black",label="non-overlapping", linewidth=lw)
	plt.legend()
	plt.title(str(genos) + "   " + str(cov))
	plt.xlabel("MAF")
	plt.ylabel("Entropy")


	# plt.savefig(statsdir+'_'+refset+"/conc_per_MAF_by_ind/entr_per_MAF_by_ind_"+refset+"_"+str(cov)+"_"+str(genos) +".pdf", dpi=fdpi, bbox_inches="tight")

	plt.show()
	plt.close()


# Plot fraction of genotypes that concofmr with the references most frequent, per MAF, for ind
# separated by kept, dropped
def  plot_compref_per_maf(datadir, statsdir, filter, annotype, refset, cov, ind, ymin):

	# errors[type][0][c] : refcomp for type for KEPT markers at maf c
	# errors[type][1][c] : refcomp for type for DROPPED markers at maf c

	# type = HC or imputed

	errors=[[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]


	plt.figure(figsize=(fwidth_med, fheight_med))
	subsets = ["kept", "dropped"]
	# subsets = ["dropped"]

	types = ["HC", "imp"+filter]

	for type in range(len(types)):
		for subset in range(len(subsets)):

			res_file = datadir + statsdir + "/compref_"+types[type]+"_"+annotype+"_"+str(cov)+"x_"+subsets[subset]+"_"+refset+"_"+ind + ".output"

			comp_table,total_refcomp = load_compref_results(res_file)

			allele_freqs = np.array(comp_table)[:,0]


			this_conc = []

			#aggregated over MAF bins for cross-checking againt reported total compref
			num_same_agg= 0.0
			num_markers_agg = 0.0

			for a in range(len(allele_freqs)):

				num_same_this_bin = comp_table[a][1]
				num_markers_this_bin = comp_table[a][2]

				if num_markers_this_bin > 0:
					average_compref_this_bin = num_same_this_bin / num_markers_this_bin
				else:
					average_compref_this_bin = 0

				print("---af: "+str(allele_freqs[a]) + " same: "+str(num_same_this_bin)+ " markers: " + str(num_markers_this_bin)+ " compref = " + str(average_compref_this_bin))

				this_conc.append(average_compref_this_bin)

				num_same_agg += num_same_this_bin
				num_markers_agg += num_markers_this_bin


			compref_agg = float(num_same_agg / num_markers_agg)

			print(types[type] + " " + subsets[subset] + " " + str(cov))
			print("Estimated compref " + str(compref_agg))
			print("Counted compref " + str(total_refcomp))


			errors[type][subset] = this_conc

	for type in range(len(types)):
		plt.plot(allele_freqs, errors[type][0], linestyle="solid", color=type_palette[type], linewidth=lw, label=subsets[0]+ " " + types[type])
		plt.plot(allele_freqs, errors[type][1], linestyle="dashed", dashes =dashstyle, color=type_palette[type], linewidth=lw, label=subsets[1]+ " " + types[type])
		print(allele_freqs)
	# plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])

	plt.ylim(ymin=ymin, ymax=1.0)
	plt.xlim(xmin=0.0, xmax=1.0)

	# plt.plot([-1],[-1], linestyle="solid", color="black", label="overlapping", linewidth=lw)
	# plt.plot([-1],[-1], linestyle="dashed", dashes =dashstyle, color="black",label="non-overlapping", linewidth=lw)

	# for type in range(len(types)):
	# 	plt.plot([-1],[-1], linestyle="solid", color=type_palette[type], label=types[type], linewidth=lw)


	plt.legend()
	plt.title(str(ind) + "   " + str(cov))
	plt.ylabel("Reference concordance")

	# plt.xlabel("MAF")
	# plt.savefig("compref_per_MAF_"+refset+"_"+str(cov)+"_"+str(ind) +".pdf", dpi=fdpi, bbox_inches="tight")

	plt.xlabel("AF")
	plt.savefig("compref_per_AAF_"+refset+"_"+str(cov)+"_"+str(ind) +".pdf", dpi=fdpi, bbox_inches="tight")

	plt.show()
	plt.close()

# Plot fraction of genotypes that concofmr with the references most frequent, per MAF, for ind
# kept, dropped together
def  plot_compref_per_maf_nosplit(datadir, statsdir, filter, annotype, refset, cov, ind):

	# comprefs[type][c] : refcomp for type at maf c
	# type = HC or imputed

	comprefs=[[[]],[[]]]

	# plt.figure(	figsize=(fwidthmed, fheight_med))
	subsets = ["kept", "dropped"]
	# subsets = ["dropped"]

	types = ["HC", "imp"+filter]

	for type in range(len(types)):

		# just for checking how many af bins
		res_file = datadir + statsdir + "/compref_"+types[type]+"_"+annotype+"_"+str(cov)+"x_"+subsets[0]+"_"+refset+"_"+ind + ".output"
		comp_table,total_refcomp = load_compref_results(res_file)
		allele_freqs = np.array(comp_table)[:,0]


		num_same_this_type = [0 for a in range(len(allele_freqs))]
		num_markers_this_type = [0 for a in range(len(allele_freqs))]
		this_type_conc = []


		for subset in range(len(subsets)):

			res_file = datadir + statsdir + "/compref_"+types[type]+"_"+annotype+"_"+str(cov)+"x_"+subsets[subset]+"_"+refset+"_"+ind + ".output"
			comp_table,total_refcomp = load_compref_results(res_file)
			allele_freqs = np.array(comp_table)[:,0]



			#aggregated over MAF bins for cross-checking againt reported total compref
			num_same_agg= 0.0
			num_markers_agg = 0.0

			for a in range(len(allele_freqs)):

				num_same_this_bin = comp_table[a][1]
				num_markers_this_bin = comp_table[a][2]

				num_same_this_type[a] += num_same_this_bin
				num_markers_this_type[a] += num_markers_this_bin

				num_same_agg += num_same_this_bin
				num_markers_agg += num_markers_this_bin


			compref_agg = float(num_same_agg / num_markers_agg)
			print("Average compref for " + types[type] +" " + subsets[subset] +" calulated: " + str(compref_agg) +" from file: " + str(total_refcomp))

		for a in range(len(num_same_this_type)):
			if num_markers_this_type[a] == 0:
				this_type_conc.append(1)
			else:
				this_type_conc.append(float(num_same_this_type[a]) / float(num_markers_this_type[a]))

		comprefs[type] = this_type_conc

	diff = []
	for a in range(len(allele_freqs)):
		diff.append(comprefs[0][a] - comprefs[1][a])

	# for type in range(len(types)):
	# 	plt.plot(allele_freqs, comprefs[type], linestyle="solid", color=type_palette[type], linewidth=lw, label=types[type])

	plt.plot(allele_freqs, diff, label = cov)

	#
	# # plt.ylim(ymin=ymin, ymax=1.0)
	#
	# plt.legend()
	# plt.title(str(ind) + "   " + str(cov))
	# plt.ylabel("Reference concordance")
	#
	# plt.xlabel("MAF")
	# plt.savefig("compref_per_MAF_nosplit"+refset+"_"+str(cov)+"_"+str(ind) +".pdf", dpi=fdpi, bbox_inches="tight")
	# plt.xlim(xmin=0.0, xmax=0.5)
	#
	# # plt.xlabel("AF")
	# # plt.savefig("compref_per_AAF_notsplit_"+refset+"_"+str(cov)+"_"+str(ind) +".pdf", dpi=fdpi, bbox_inches="tight")
	# # plt.xlim(xmin=0.0, xmax=1.0)
	#
	#
	# plt.show()
	# plt.close()



# Plot fraction of genotypes that concofmr with the references most frequent, per MAF, for ind
# kept, dropped together
def  get_compref_per_maf_nosplit(datadir, statsdir, filter, annotype, refset, cov, ind):

	# comprefs[type][c] : refcomp for type at maf c
	# type = HC or imputed

	comprefs=[[[]],[[]]]

	# plt.figure(	figsize=(fwidthmed, fheight_med))
	subsets = ["kept", "dropped"]
	# subsets = ["dropped"]

	types = ["HC", "imp"+filter]

	for type in range(len(types)):

		# just for checking how many af bins
		res_file = datadir + statsdir + "/compref_"+types[type]+"_"+annotype+"_"+str(cov)+"x_"+subsets[0]+"_"+refset+"_"+ind + ".output"
		comp_table,total_refcomp = load_compref_results(res_file)
		allele_freqs = np.array(comp_table)[:,0]


		num_same_this_type = [0 for a in range(len(allele_freqs))]
		num_markers_this_type = [0 for a in range(len(allele_freqs))]
		this_type_conc = []


		for subset in range(len(subsets)):

			res_file = datadir + statsdir + "/compref_"+types[type]+"_"+annotype+"_"+str(cov)+"x_"+subsets[subset]+"_"+refset+"_"+ind + ".output"
			comp_table,total_refcomp = load_compref_results(res_file)
			allele_freqs = np.array(comp_table)[:,0]



			#aggregated over MAF bins for cross-checking againt reported total compref
			num_same_agg= 0.0
			num_markers_agg = 0.0

			for a in range(len(allele_freqs)):

				num_same_this_bin = comp_table[a][1]
				num_markers_this_bin = comp_table[a][2]

				num_same_this_type[a] += num_same_this_bin
				num_markers_this_type[a] += num_markers_this_bin

				num_same_agg += num_same_this_bin
				num_markers_agg += num_markers_this_bin


			compref_agg = float(num_same_agg / num_markers_agg)
			print("Average compref for " + types[type] +" " + subsets[subset] +" calulated: " + str(compref_agg) +" from file: " + str(total_refcomp))

		for a in range(len(num_same_this_type)):
			if num_markers_this_type[a] == 0:
				this_type_conc.append(1)
			else:
				this_type_conc.append(float(num_same_this_type[a]) / float(num_markers_this_type[a]))

		comprefs[type] = this_type_conc

	diff = []
	print(ind + " " + str(cov))
	print("HC aff")
	print(comprefs[0])
	print("imp aff")

	print(comprefs[1])

	for a in range(len(allele_freqs)):
		diff.append(comprefs[0][a] - comprefs[1][a])

	# num:markers_this_type is num markers for imputed data

	return allele_freqs,diff, num_markers_this_type




# Plot fraction of genotypes that concofmr with the references most frequent, per MAF, for ind
def  plot_HCcompref_per_maf(datadir, statsdir, filter, annotype, refset, cov, ind, ymin, genos=[0,1,2] , byminor=False):

	# errors[type][0][c] : refcomp for type for KEPT markers at maf c
	# errors[type][1][c] : refcomp for type for DROPPED markers at maf c

	# type = HC or imputed

	errors=[[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]


	plt.figure(figsize=(fwidth_med, fheight_med))
	subsets = ["kept", "dropped"]
	# subsets = ["dropped"]

	types = ["HC"]

	byminor_str = ""
	if byminor:
		byminor_str = "_byminor_"



	for type in range(len(types)):
		for subset in range(len(subsets)):

			res_file = datadir + statsdir + "/HC"+byminor_str+"compref_"+types[type]+"_"+annotype+"_"+str(cov)+"x_"+subsets[subset]+"_"+refset+"_"+ind + ".output"

			comp_table,total_refcomp = load_compref_results(res_file)

			# print(comp_table)
			# exit()
			allele_freqs = np.array(comp_table)[:,0]


			this_conc = []

			#aggregated over MAF bins for cross-checking againt reported total compref
			num_same_agg= 0.0
			num_markers_agg = 0.0

			for a in range(len(allele_freqs)):

				num_same_this_bin = []
				num_markers_this_bin = []

				for geno in genos:
					num_same_this_bin_this_geno = comp_table[a][1 + geno]
					num_markers_this_bin_this_geno = comp_table[a][4 + geno]
					num_same_this_bin.append(num_same_this_bin_this_geno)
					num_markers_this_bin.append(num_markers_this_bin_this_geno)

				average_compref_this_bin = float(np.sum(num_same_this_bin)) / np.sum(num_markers_this_bin)
				print("---af: "+str(allele_freqs[a])  + "average compref = " + str(average_compref_this_bin))

				this_conc.append(average_compref_this_bin)

				num_same_agg += np.sum(num_same_this_bin)
				num_markers_agg += np.sum(num_markers_this_bin)


			compref_agg = float(num_same_agg) / num_markers_agg

			print(types[type] + " " + subsets[subset] + " " + str(cov))
			print("Estimated compref " + str(compref_agg))
			print("Counted compref " + str(total_refcomp))


			errors[type][subset] = this_conc

	for type in range(len(types)):
		plt.plot(allele_freqs, errors[type][0], linestyle="solid", color=type_palette[type], linewidth=lw)
		plt.plot(allele_freqs, errors[type][1], linestyle="dashed", dashes =dashstyle, color=type_palette[type], linewidth=lw)
		print(allele_freqs)

	xmax = 1.0
	# plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])

	plt.ylim(ymin=ymin, ymax=1.0)
	plt.xlim(xmin=0.0, xmax=xmax)

	# plt.plot([-1],[-1], linestyle="solid", color="black", label="overlapping", linewidth=lw)
	# plt.plot([-1],[-1], linestyle="dashed", dashes =dashstyle, color="black",label="non-overlapping", linewidth=lw)

	for type in range(len(types)):
		plt.plot([-1],[-1], linestyle="solid", color=type_palette[type], label=types[type], linewidth=lw)


	plt.legend()
	plt.title(str(ind) + "   " + str(cov))
	plt.xlabel("MAF")
	plt.ylabel("Reference concordance")


	plt.savefig("HC"+byminor_str+"compref_per_MAF_"+refset+"_"+str(cov)+"_"+str(ind) +"_"+str(genos)+".pdf", dpi=fdpi, bbox_inches="tight")

	plt.show()
	plt.close()

# Plot entropy for 5 individuals for MAF  on x-axis for genotypes genos
# genos is a list that may contain at least 1 of: 0,1,2
# entropy between low cov files and imputed

def  plot_lc_entropy_per_maf_by_ind(datadir, statsdir, refset, cov, genos, ymax):

	# errors[ind][0][c] : entropy for ind for KEPT markers at maf c
	# errors[ind][1][c] : entropy for ind for DROPPED markers at maf c

	errors=[[[]],[[]],[[]],[[]],[[]]]


	plt.figure(figsize=(fwidth_med, fheight_med))
	# subsets = ["dropped"]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]
	# inds = ["ans17"]


	for ind in range(len(inds)):

		# res_file="~/Projects/PythonUtils/General/test.output"
		res_file = datadir + statsdir + "/entropy_lcimp"+impset+"_"+str(cov)+"x_"+refset+"_"+inds[ind]+".output"

		e_table,tot_entropy =load_entropy_results(res_file)

		allele_freqs = np.array(e_table)[:,0]


		this_conc = []

		#aggregated over MAF bins for cross-checking againt reported GcSs
		entropy_agg= 0.0
		num_markers_agg = 0.0

		for a in range(len(allele_freqs)):

			entropy_this_bin = 0.0
			num_markers_this_bin = 0.0

			for geno in genos:

				entropy_this_geno = e_table[a][1+geno]
				num_markers_this_geno = e_table[a][4+geno]

				entropy_this_bin += entropy_this_geno*num_markers_this_geno
				num_markers_this_bin += num_markers_this_geno


			average_entropy_this_bin = entropy_this_bin / num_markers_this_bin
			print("---af: "+str(allele_freqs[a]) + " num markers= " + str(num_markers_this_bin) + "avg entropy= " + str(average_entropy_this_bin))

			this_conc.append(average_entropy_this_bin)

			entropy_agg += average_entropy_this_bin * num_markers_this_bin
			num_markers_agg += num_markers_this_bin


		avg_entropy_agg = float(entropy_agg / num_markers_agg)

		# print("Total num markers = " + str(_+num_markers_agg) + "   concordance = " + str(avg_entropy_agg) + "  gcss concordance = " + str(concordance_gcss))

		print(inds[ind] + " "  + " " + str(cov) + " " + str(genos))
		print("Estimated entropy " + str(avg_entropy_agg))
		print("Counted entropy " + str(tot_entropy))


		errors[ind][0] = this_conc


	for ind in range(len(inds)):
		plt.plot(allele_freqs, errors[ind][0], linestyle="solid", color=ind_palette[ind], linewidth=lw)
		# plt.plot(allele_freqs, errors[ind][1], linestyle="dashed", dashes =dashstyle, color=ind_palette[ind], linewidth=lw)
		print(allele_freqs)
	plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
	# plt.legend()
	# plt.title("Non-reference discordance of imputed")
	plt.ylim(ymin=0.0, ymax=ymax)
	plt.xlim(xmin=0.0, xmax=0.5)

	plt.plot([-1],[-1], linestyle="solid", color="black", label="overlapping", linewidth=lw)
	plt.plot([-1],[-1], linestyle="dashed", dashes =dashstyle, color="black",label="non-overlapping", linewidth=lw)
	plt.legend()
	plt.title(str(genos) + "   " + str(cov) + " LC")
	plt.xlabel("MAF")
	plt.ylabel("Entropy")


	plt.savefig(statsdir+'_'+refset+"/conc_per_MAF_by_ind/entr_lc_per_MAF_by_ind_"+refset+"_"+str(cov)+"_"+str(genos) +".pdf", dpi=fdpi, bbox_inches="tight")

	plt.show()
	plt.close()


#
# Calculate overall genotype concordance based on stats from the GCsS for genotypes genos
# genos is a list that may contain at least 1 of: 'RR Hom' 'Ra Het' 'AA Hom'
# lc indicates if low cov comparison should be used
#
def get_concordance_from_gcss(dir, refset, cov, ind, subset, genos, lc=False):

	if lc:
		print("GcSs: reading " + dir+'_'+refset+'/all_stats_lc_'+str(cov)+'_'+refset+'.json')
		with open(dir+'_'+refset+'/all_stats_lc_'+str(cov)+'_'+refset+'.json', 'r') as fp:
			all_tables = json.load(fp)

	else:
		print("GcSs: reading " + dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json')
		with open(dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json', 'r') as fp:
			all_tables = json.load(fp)

	gct_dict = all_tables[ind][subset][2]["# GCsS"]

	num_match = 0
	num_mismatch = 0

	for geno in genos:
		print(float(gct_dict[geno+' matches']))
		num_match += float(gct_dict[geno+' matches'])
		num_mismatch += float(gct_dict[geno+' mismatches'])



	concordance = num_match / (num_match + num_mismatch)

	print(ind + " " + subset + " " + str(cov) + " " + str(genos))
	print("GcSs concordance = " + str(concordance))
	print("---- match: " + str(num_match))
	print("---- mismatch: " + str(num_mismatch))
	print("---- total: " + str(num_match + num_mismatch))

	return concordance


#calculate overall genotype concordance based on stats from the GCTs
def calc_concordance_from_gcts(dir, cov, refset):

	with open(dir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+'.json', 'r') as fp:
		all_tables = json.load(fp)

	table=[["Sample", "Kept", "Dropped", "All"]]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]

	subsets = ["kept", "dropped"]


	for ind in inds:


		# for this individual : kept and dropped gcts

		ALL_wrongs = 0
		ALL = 0

		this_row = [ind]

		for subset in subsets:

			gct_dict = all_tables[ind][subset][1]["# GCTs"]

			# markers that are RR in True file

			num_RR_RA = gct_dict['RR Hom -> RA Het']
			num_RR_AA = gct_dict['RR Hom -> AA Hom']
			num_RR_RR = gct_dict['RR Hom -> RR Hom']


			RR_wrongs= int(num_RR_RA) + int(num_RR_AA)
			RR_all= int(num_RR_RA) + int(num_RR_AA) + int(num_RR_RR)

			# RA
			num_RA_RR = gct_dict['RA Het -> RR Hom']
			num_RA_AA = gct_dict['RA Het -> AA Hom']
			num_RA_RA = gct_dict['RA Het -> RA Het']


			RA_wrongs= int(num_RA_RR) + int(num_RA_AA)
			RA_all= int(num_RA_RR) + int(num_RA_AA) + int(num_RA_RA)


			# AA
			num_AA_RR = gct_dict['AA Hom -> RR Hom']
			num_AA_RA = gct_dict['AA Hom -> RA Het']
			num_AA_AA = gct_dict['AA Hom -> AA Hom']

			AA_wrongs= int(num_AA_RR) + int(num_AA_RA)
			AA_all= int(num_AA_RR) + int(num_AA_RA) + int(num_AA_AA)

			# Total genotype concordance
			ALL_wrongs_subset = RR_wrongs + RA_wrongs + AA_wrongs
			ALL_subset = RR_all + RA_all + AA_all



			ALL_wrongs += ALL_wrongs_subset
			ALL += ALL_subset

			print("----sample:" + ind)
			print("----subset:" + subset)

			print("----all wrongs: " + str(ALL_wrongs_subset))
			print("----all sites : " + str(ALL_subset))
			conc = 1.0 - float(ALL_wrongs_subset)/float(ALL_subset)
			print("----concordance: " + str(conc))

			this_row.append(str(conc))

		print("sample:" + ind)
		print("total:")
		print("all wrongs: " + str(ALL_wrongs))
		print("all sites: " + str(ALL))
		conc = 1.0 - float(ALL_wrongs)/float(ALL)

		print("concordance: " + str(conc))
		this_row.append(str(conc))
		table.append(this_row)

	with open(dir +"_"+refset+ "/"+ "concordance_"+str(cov)+"_"+refset+".csv", 'w') as outfile:
		for l in table:
			outfile.write(",".join(l)+"\n")





def plot_gcts_ind(dir):

	with open(dir + 'all_tables.json', 'r') as fp:
		all_tables = json.load(fp)


	inds = ["ans17", "car16", "LBK", "Loschbour", "ne1", "ne1CP"]
	# inds = ["ans17", "car16", "LBK", "Loschbour", "ne1"]

	print(inds)

	filters = ["all", "filtered", "thets"]
	subsets = ["kept", "dropped"]


	for ind in inds:

		fig, ax = plt.subplots(figsize=(5,9))
		# plt.tight_layout()
		plt.suptitle("Genotype errors "+ind)
		plt.subplots_adjust(top=0.93)
		counter = 0

		for filter in filters:
			counter += 1
			plt.subplot(3,1,counter)

			all_gcts = []
			for subset in subsets:

				gct_dict = all_tables[ind][filter][subset][1]["# GCTs"]
				keys = ["Ref/Ref discordance", "Ref/Alt discordance", "Alt/Alt discordance","NRD"]
				keys = [u'RR Hom -> RA Het',  u'RR Hom -> AA Hom',u'RA Het -> RR Hom', u'RA Het -> AA Hom',  u'AA Hom -> RR Hom',  u'AA Hom -> RA Het']
				this_gcts = []
				for k in keys:
					this_gcts.append(gct_dict[k])

				all_gcts.append(this_gcts)
			print(all_gcts)

			x = [0,10,20,30]
			x = range(len(keys))
			x = map(lambda y : 10*y, x)
			print(x)

			plt.bar(map(lambda y : y-1, x),all_gcts[0], width = 1, label=subsets[0], color="black")
			plt.bar(map(lambda y : y-0, x),all_gcts[1], width = 1, label=subsets[1], color="grey")

			if counter == 1:
				# plt.title(subset)
				plt.title("All Markers")

				plt.legend(fontsize=7.5)
			if counter == 2:
			# 	plt.title(subset)
				plt.title("Quality Filtered")

			if counter == 3:
			# 	plt.title(subset)
				plt.title("Trusted Heterozygotes")

			# xlabels = ["R/R discordance", "R/A discordance", "A/A discordance","Non-Ref Discordance"]
			# xlabels = ["Hom-Ref", "Het", "Hom-Alt","Non-Ref"]
			xlabels = keys
			# print(gct_dict.keys())


			plt.xticks(x,xlabels,size=4)
			# plt.ylim(ymax=11)
			plt.xticks(rotation=-15)


		# plt.legend()
		plt.savefig(dir+"Genotype_errors_"+ind+".pdf")
		# plt.show()


def plot_table_transposed(table_rows, filename, title=""):
	fig, ax = plt.subplots()
	ax.xaxis.set_visible(False)
	ax.yaxis.set_visible(False)
	table = _plot_table_transposed(ax,table_rows)
	table.set_fontsize(10)
	plt.title(title)
	plt.show()
	# plt.savefig(filename+".pdf")


def _plot_table_transposed(ax, cellText, id=None):

	if id == None:
		cellTextt = np.transpose(cellText)[2:]

	if id == "GCTs":
		cellTextt = np.transpose(cellText)[7:10]
	if id == "NRDs":
		cellTextt = np.transpose(cellText)[2:]


	table = ax.table(cellText = cellTextt,loc='center')
	return table


def plot_table(table_rows, filename):

	colnames = table_rows[0]
	table_rows = table_rows[1:]
	fig, ax = plt.subplots()
	ax.xaxis.set_visible(False)
	ax.yaxis.set_visible(False)
	table = ax.table(cellText = table_rows,colLabels = colnames, loc='center')
	plt.savefig(filename+".pdf")


# plot table 'id' for individual 'sample'
# for
# the three filter levels: all, quality filtered, thets only
# for
# the two subsets HC_all HC_intersected in all cases
def plot_table_ind(id, ind, datadir,resdir,ref=""):

	if not ref == "":
		ref = "_" + ref

	filters = ["", "_filtered", "_thets"]
	subsets = ["only", "intersect"]
	titles = ["All markers :       dropped", "kept","Filtered :       dropped", "kept", "Trusted hets :      dropped", "kept"]

	# table_rows_matrix = [[],[],[]]

	fig, ax = plt.subplots()
	plt.suptitle("Genotype discordance " + ind)
	# ax.xaxis.set_visible(False)
	# ax.yaxis.set_visible(False)
	counter = 0

	for f in range(len(filters)):
		for s in range(len(subsets)):
			counter += 1
			print(datadir+"stats_HCimp_1x_"+subsets[s]+filters[f]+"_"+ind+".output")
			table_rows = get_table(datadir+"stats_HCimp_1x_"+subsets[s]+filters[f]+ref+"_"+ind+".output", id)
			# table_rows = []

			print(table_rows)
			csvfile = open(id+"_"+ind+"_"+filters[f]+"_"+subsets[s]+".csv", 'w')
			writer = csv.writer(csvfile, delimiter=',')
			writer.writerows(table_rows)

			ax = plt.subplot("32"+str(counter))
			plt.title(titles[counter-1])
			ax.xaxis.set_visible(False)
			ax.yaxis.set_visible(False)
			# tester = [[1,2],[3,4],[5,6]]
			# _plot_table_transposed(ax,tester)
			_plot_table_transposed(ax,table_rows,id=id)
	plt.savefig(resdir + id + "_" + ind +".pdf")
	# plt.show()
	plt.close()


# plot table 'id' for 'filter'
# for
# all individuals
# for
# the two subsets HC_all HC_intersected in all cases
def plot_all_inds(filter, id, datadir,resdir, ref=""):

	if not ref == "":
		ref = "_" + ref
	subsets = ["only", "intersect"]
	inds = ["ans17", "car16", "LBK", "Loschbour", "ne1", "ne1CP"]
	titles = [inds[0] + " :       dropped", "kept",
			  inds[1] +" :       dropped",  "kept",
			  inds[2] +" :      dropped", "kept",
			  inds[3] +" :      dropped", "kept",
			  inds[4] +" :      dropped", "kept",
  			  inds[5] +" :      dropped", "kept"]


	# table_rows_matrix = [[],[],[]]

	fig, ax = plt.subplots()
	plt.suptitle("Genotype concordance table " + filter)
	# ax.xaxis.set_visible(False)
	# ax.yaxis.set_visible(False)
	counter = 0
	for i in range(len(inds)):
		for s in range(len(subsets)):
			counter += 1
			print(datadir+"stats_HCimp_1x_"+subsets[s]+filter+ref+"_"+inds[i]+".output")
			table_rows = butils.get_table(datadir+"stats_HCimp_1x_"+subsets[s]+filter+ref+"_"+inds[i]+".output", id)


			csvfile = open(id+"_"+inds[i]+"_"+filter+"_"+subsets[s]+".csv", 'w')
			writer = csv.writer(csvfile, delimiter=',')
			writer.writerows(table_rows)

			print(table_rows)
			ax = plt.subplot(6,2,counter)
			plt.title(titles[counter-1])
			ax.xaxis.set_visible(False)
			ax.yaxis.set_visible(False)
			# tester = [[1,2],[3,4],[5,6]]
			# _plot_table_transposed(ax,tester)
			_plot_table_transposed(ax,table_rows)
	plt.savefig(resdir + id + filter+ "_all_inds_ne1CP.pdf")
	# plt.show()
	plt.close()


def plot_avg_compref_diffs(resdir, grid=False, tot_markers=False):
	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]

	linewidths = [lw, lw, lw, lw*0.8, lw*0.8, lw*0.7, lw*0.7, lw*0.7, lw*0.7]

	linestyles = ["solid","dashed", "solid","dashed", "solid","dashed", "solid","dashed", "solid"]
	linestyles = ["solid" for i in range(9)]

	colors=["black", "red"]

	# fig = plt.figure(figsize=(fwidth_small, fheight_med_smaller))

	fig,ax = plt.subplots(figsize=(fwidth_small, fheight_med_smaller))

	if tot_markers:
		ax2 = ax.twinx()



	combos = np.array(np.meshgrid(linewidths,colors,linestyles)).T.reshape(-1,3)

	print(combos)
	colors = plt.cm.viridis_r(np.linspace(0,1,9))
	# colors = plt.cm.coolwarm(np.linspace(0,1,9))
	# colors = sns.color_palette(sns.color_palette("gray", 9))
	colors = sns.cubehelix_palette(9, start=.5, rot=-.75)


	# cols = sns.color_palette(sns.color_palette("RdBu_r", 13))
	# colors= []
	# for i in range(4):
	# 	colors.append(cols[i])
	#
	# colors.append([0,0,0,0.8])
	#
	# for i in range(4):
	# 	colors.append(cols[-(4-i)])
	colors = colors[::-1]




	bounds = covs
	# cmap=sns.color_palette("gray")
	cmap= LinearSegmentedColormap.from_list(
        "cols", colors, N=9)


	Z = [[0,0],[0,0]]
	levels = covs+[2.25]
	CS3 = plt.contourf(Z, levels, cmap=cmap)
	# plt.clf()

	# print(colors)
	# colors[4] = [0,0,0,0.8]

	# colors = sns.color_palette(sns.color_palette("GnBu_d", 9))
	# colors = sns.color_palette(sns.color_palette("YlOrRd", 9))
	# colors = sns.color_palette(sns.color_palette("YlOrBr", 11))[2:]


	if grid:
		plt.grid(b=True, which='major', linewidth=lw2)
		plt.grid(b=True, which='minor', linewidth=lw3, linestyle="--")

	lines = []

	for cov in covs:
		diffs_inds = []
		for ind in inds:
			freqs, diffs, num_markers_imp = get_compref_per_maf_nosplit(datadir, 'compref_rmts_61_5_impnofilter_annofull', "nofilter","annofull", "EUR", cov, ind)
			diffs_inds.append(diffs)


		if cov == 1.0:
			num_markers_imp_1x = num_markers_imp

		print(diffs_inds)
		print("averaging " + str(cov))
		diffs_avg = np.average(diffs_inds, axis=0)
		print(diffs_avg.shape)
		print(diffs_avg)
		print(freqs.shape)



		# freqs = freqs[0:-1]
		# diffs_avg = diffs_avg[0:-1]
		#
		# freqs = rebin(freqs, 1)
		# diffs_avg = rebin(diffs_avg, 1)



		# plt.plot(freqs, diffs_avg, label=str(cov), linewidth=float(combos[covs.index(cov)][0]),linestyle=combos[covs.index(cov)][2], color= colors[covs.index(cov)])
		lin = ax.plot(freqs[0:-1], diffs_avg[0:-1], label="{0:.2f}".format(cov), linewidth=linewidths[covs.index(cov)],linestyle=linestyles[covs.index(cov)], color= colors[covs.index(cov)])
		lines.append(lin)


	if tot_markers:
		print("Num markers 1x")
		print(num_markers_imp_1x)
		ax2.plot(freqs, num_markers_imp_1x, color= "black")

	else:
		sns.despine()
		covticks = [i+0.125 for i in covs]
		covticks[0] = 0.175

		cbar = plt.colorbar(CS3,ticks=covticks)


		cbar.ax.tick_params(length=0, pad=2)
		cbar.outline.set_visible(False)
		# cbar.dividers.set_linewidth(0)


		cbar.ax.set_yticklabels(['{0:.2f}'.format(cov) for cov in covs])

	if not grid:
		ax.plot(freqs, np.zeros(len(freqs)),linewidth=lw/2, color="gray", linestyle="dotted")


	# legend=plt.legend(title="Coverage (x)")
	# plt.setp(legend.get_title(),fontsize=legfontsize)
	# legend.get_frame().set_linewidth(lw2)



	ax.set_ylabel("Reference affinity difference")
	plt.xlabel("MAF")
	print(freqs[0:-1])
	# plt.ticklabel_format(style='sci',scilimits=(0,0),axis='y')

	plt.xlim(xmin=min(freqs), xmax=0.5)

	ax.set_ylim(ymax=0.005, ymin=-0.02)
	ax.set_ylim(ymax=0.02, ymin=-0.02)

	# (left, bottom, right, top)
	plt.tight_layout(w_pad=0, h_pad=0, rect=[-0.045,-0.09,1.05,1.07])
	plt.savefig(resdir+"compref_diff_avginds.pdf", dpi=fdpi)

	# plt.savefig(resdir+"compref_diff_avginds_w_markers.pdf", dpi=fdpi)

	plt.show()
	plt.close()




# given a directory with bcf stats files from mesolithic data
# write NRDs, GCTs and GCsS to a json file in
# statsdir: the name of the folder where results are. they will be written to a folder with the same name in this project dir
# datadir: the location of the statsdir folder
#
# for a given coverage and refernce
# eg write_stats_to_json(/media/kristiina/My Passport/Data/Mesolithic/results/concordance/, stats_rmts_61_5/):
#
def write_stats_to_json(datadir, statsdir, cov, refset, impset, individual=""):

	all_stats = {
	"ans17" : {"dropped" : [],
				 "kept": []
			 },
	"LBK" : {"dropped" : [],
				 "kept": []
			 },
	"Loschbour" : {"dropped" : [],
				 "kept": []
			 },
	"sf12" : {"dropped" : [],
				 "kept": []
			 },
	"ne1" : {"dropped" : [],
				 "kept": []
			 },
	"ne1CP" : {"dropped" : [],
				 "kept": []
			 },
	}

	subsets = ["dropped", "kept"]
	subsets_tags = ["dropped", "kept"]

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]


	for ind in inds:
		for s in range(len(subsets)):
			statsfile = datadir + statsdir + "/stats_HCimp"+impset+"_"+str(cov)+"x_"+subsets[s]+"_"+refset+"_"+ind+".output"

			print("Doing statfle " + statsfile)

			# statsfile = datadir + statsdir + "/stats_HCimpnofilter_"+str(cov)+"x_"+subsets[s]+"_"+refset+"_"+ind+".output"

			print("doing NRDs")
			table_rows = butils.get_table(statsfile, "NRDs")
			nrd_dict = butils.table_to_dict(table_rows)
			# print(nrd_dict)
			print("doing GCTs")
			table_rows = butils.get_table(statsfile, "GCTs")
			gct_dict = butils.table_to_dict(table_rows)
			# print(gct_dict)
			print("doing GCsS")
			table_rows = butils.get_table(statsfile, "GCsS")
			gcs_dict = butils.table_to_dict(table_rows)
			# print(gcs_dict)
			print("doing GCsAF")
			table_rows = butils.get_table(statsfile, "GCsAF")
			gcsAF_dict = butils.table_to_dict(table_rows)
			print(gcsAF_dict)



			all_stats[ind][subsets_tags[s]].append(nrd_dict)
			all_stats[ind][subsets_tags[s]].append(gct_dict)
			all_stats[ind][subsets_tags[s]].append(gcs_dict)
			all_stats[ind][subsets_tags[s]].append(gcsAF_dict)

	print("writing to " + statsdir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+individual+'.json')
	with open(statsdir+'_'+refset+'/all_stats_'+str(cov)+'_'+refset+individual+'.json', 'w') as fp:
	# with open(statsdir+'_'+refset+'/all_stats_impnofilter_'+str(cov)+'_'+refset+'.json', 'w') as fp:
		json.dump(all_stats, fp, sort_keys=True,indent=4)




# given a directory with bcf stats files from mesolithic data
# write NRDs, GCTs and GCsS to a json file in
# statsdir: the name of the folder where results are. they will be written to a folder with the same name in this project dir
# datadir: the location of the statsdir folder
#
# for a given coverage and refernce
# eg write_stats_to_json(/media/kristiina/My Passport/Data/Mesolithic/results/concordance/, stats_rmts_61_5/):
#
def write_lc_stats_to_json(datadir, statsdir, cov, refset, impset, individual=""):

	all_stats = {
	"ans17" : {"all" : []
			 },
	"LBK" : {"all" : []
			 },
	"Loschbour" : {"all" : []
			 },
	"sf12" : {"all" : []
			 },
	"ne1" : {"all" : []
			 },
	"ne1CP" : {"all" : []
			 },
	}

	inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]


	for ind in inds:
		statsfile = datadir + statsdir + "/stats_lcimp"+impset+"_"+str(cov)+"x_"+refset+"_"+ind+".output"

		print("Doing statfle " + statsfile)

		# statsfile = datadir + statsdir + "/stats_HCimpnofilter_"+str(cov)+"x_"+subsets[s]+"_"+refset+"_"+ind+".output"

		print("doing NRDs")
		table_rows = butils.get_table(statsfile, "NRDs")
		nrd_dict = butils.table_to_dict(table_rows)
		# print(nrd_dict)
		print("doing GCTs")
		table_rows = butils.get_table(statsfile, "GCTs")
		gct_dict = butils.table_to_dict(table_rows)
		# print(gct_dict)
		print("doing GCsS")
		table_rows = butils.get_table(statsfile, "GCsS")
		gcs_dict = butils.table_to_dict(table_rows)
		# print(gcs_dict)
		print("doing GCsAF")
		table_rows = butils.get_table(statsfile, "GCsAF")
		gcsAF_dict = butils.table_to_dict(table_rows)
		print(gcsAF_dict)



		all_stats[ind]["all"].append(nrd_dict)
		all_stats[ind]["all"].append(gct_dict)
		all_stats[ind]["all"].append(gcs_dict)
		all_stats[ind]["all"].append(gcsAF_dict)

	print("writing to " + statsdir+'_'+refset+'/all_stats_lc_'+str(cov)+'_'+refset+individual+'.json')
	with open(statsdir+'_'+refset+'/all_stats_lc_'+str(cov)+'_'+refset+individual+'.json', 'w') as fp:
	# with open(statsdir+'_'+refset+'/all_stats_impnofilter_'+str(cov)+'_'+refset+'.json', 'w') as fp:
		json.dump(all_stats, fp, sort_keys=True,indent=4)



# all_tables = {
# "ans17" : {"all" : {
# 			 "dropped" : [],
# 			 "kept": []
# 		 },
# 		 "filtered" : {
# 			 "dropped" : [],
# 			 "kept": []
# 		 },
# 		 "thets" : {
# 			"dropped" : [],
# 			 "kept": []
#
# 		 }},
# "car16" : {"all" : {
# 			 "dropped" : [],
# 			 "kept": []
# 		 },
# 		 "filtered" : {
# 			 "dropped" : [],
# 			 "kept": []
# 		 },
# 		 "thets" : {
# 			"dropped" : [],
# 			 "kept": []
#
# 		 }},
# "LBK" : {"all" : {
# 			 "dropped" : [],
# 			 "kept": []
# 		 },
# 		 "filtered" : {
# 			 "dropped" : [],
# 			 "kept": []
# 		 },
# 		 "thets" : {
# 			"dropped" : [],
# 			 "kept": []
#
# 		 }},
# "Loschbour" : {"all" : {
# 			 "dropped" : [],
# 			 "kept": []
# 		 },
# 		 "filtered" : {
# 			 "dropped" : [],
# 			 "kept": []
# 		 },
# 		 "thets" : {
# 			"dropped" : [],
# 			 "kept": []
#
# 		 }},
# "ne1" : {"all" : {
# 			 "dropped" : [],
# 			 "kept": []
# 		 },
# 		 "filtered" : {
# 			 "dropped" : [],
# 			 "kept": []
# 		 },
# 		 "thets" : {
# 			"dropped" : [],
# 			 "kept": []
#
# 		 }},
# "ne1CP" : {"all" : {
# 			 "dropped" : [],
# 			 "kept": []
# 		 },
# 		 "filtered" : {
# 			 "dropped" : [],
# 			 "kept": []
# 		 },
# 		 "thets" : {
# 			"dropped" : [],
# 			 "kept": []
#
# 		 }},
# }





# These are the rmts_lc_61_5 runs also, but with old post-imputation quality filtering.
#full_1x_dir = 'stats_1x_FULL/'
# eur_1x_dir = 'stats_1x_EUR/'


cov=1
ref="FULL"
impset="nofilter"
impset=""

ref='EUR'
# impset="nofilter"
impset="nofilter_anno"
# impset="_anno"
impset="nofilter_annofull"
impset="nofilter_aafanno"


datadir = '/home/kristiina/reslink/concordance/'
# statsdir = 'stats_rmts_61_5'
# statsdir = 'stats_rmts_ind'


##### IMPUTATION CONFIGURATION - INDIV or ANCIENT ##############3
statsdir = 'stats_rmts_61_5_imp'+impset
# statsdir = 'stats_rmts_ind_'+impset



covs = [0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75 ,2]
# covs = [1]


inds = ["ans17", "sf12", "LBK", "Loschbour", "ne1"]
#
# for cov in covs:
# 	write_stats_to_json(datadir,statsdir,cov,ref, impset)
	# write_lc_stats_to_json(datadir,statsdir,cov,ref, impset)









################################### metrics per MAF #####################################################

# for ind in inds:
# 	plot_compref_per_maf(datadir, 'compref_rmts_61_5_impnofilter_annofull', "nofilter","annofull", "EUR", 1, ind, 0.25)

# inds = ["Loschbour", "ne1"]
# for ind in inds:
# 	plot_compref_per_maf(datadir, 'compref_rmts_61_5_impnofilter_annofull', "nofilter","annofull", "FULL", 1, ind, 0.25)


# plot_compref_per_maf(datadir, 'compref_rmts_61_5_impnofilter_anno', "nofilter","anno", "EUR", 1, "sf12", 0.25)

# plot_compref_per_maf(datadir, 'compref_rmts_61_5_impnofilter_aafanno', "nofilter","aafanno", "EUR", 1, "ans17", 0.25)

# plot_compref_per_maf_nosplit(datadir, 'compref_rmts_61_5_impnofilter_aafanno', "nofilter","aafanno", "EUR", 1, "ans17", 0.25)


############## compref mer MAF, all covs in one plot: one plot per ind ###########################
# covs = [0.1, 0.25, 0.5, 1]
# covs = [2]




# plt.figure(	figsize=(fwidthmed, fheight_med))
# ind = "ne1"
#
# for cov in covs:
# 	plot_compref_per_maf_nosplit(datadir, 'compref_rmts_61_5_impnofilter_annofull', "nofilter","annofull", "EUR", cov, ind)
# # plot_compref_per_maf_nosplit(datadir, 'compref_rmts_61_5_impnofilter_annofull', "nofilter","annofull", "FULL", 1, "ans17")
#
# plt.legend()
# plt.ylabel("Reference concordance")
#
# plt.xlabel("MAF")
#
# plt.xlim(xmin=0.0, xmax=0.5)
# plt.ylim( ymax=0.01)
#
# plt.savefig("compref_diff"+ind+".pdf", bbox_inches="tight")

# ###############



# plt.show()
# plt.close()
#################################### HC compref #################################################################

inds = ["ans17"]
# for ind in inds:
# 	plot_HCcompref_per_maf(datadir, 'compref_rmts_61_5_impnofilter_annofull', "nofilter","annofull", "EUR", 1, ind, 0.0, genos=[0], byminor=True)
# 	plot_HCcompref_per_maf(datadir, 'compref_rmts_61_5_impnofilter_annofull', "nofilter","annofull", "EUR", 1, ind, 0.0, genos=[1], byminor=True)
# 	plot_HCcompref_per_maf(datadir, 'compref_rmts_61_5_impnofilter_annofull', "nofilter","annofull", "EUR", 1, ind, 0.0, genos=[2], byminor=True)

# for ind in inds:
	# plot_HCcompref_per_maf(datadir, 'compref_rmts_61_5_impnofilter_aafanno', "nofilter","aafanno", "EUR", 1, ind, 0.0, genos=[0], byminor=False)
	# plot_HCcompref_per_maf(datadir, 'compref_rmts_61_5_impnofilter_aafanno', "nofilter","aafanno", "EUR", 1, ind, 0.0, genos=[1], byminor=False)
	# plot_HCcompref_per_maf(datadir, 'compref_rmts_61_5_impnofilter_aafanno', "nofilter","aafanno", "EUR", 1, ind, 0.0, genos=[2], byminor=False)



cov = 1
# plot_entropy_per_maf_by_ind(datadir, statsdir, ref, cov, [1], ymax=1.0)
# plot_entropy_per_maf_by_ind(datadir, statsdir, ref, cov, [0,1,2], ymax=8)

# plot_lc_entropy_per_maf_by_ind(datadir, statsdir, ref, cov, [1], ymax=1.0)
# plot_lc_entropy_per_maf_by_ind(datadir, statsdir, ref, cov, [0,1,2], ymax=8)

# plot_lc_concordance_per_maf_by_ind(statsdir, ref, cov, ['RR Hom','RA Het','AA Hom'], ymin=0.0)
# plot_concordance_per_maf_by_ind(statsdir, ref, cov, ['RR Hom','RA Het','AA Hom'], ymin=0.95)

# plot_lc_concordance_per_maf_by_ind(statsdir, ref, cov, ['RA Het'], ymin=0.0)
# plot_concordance_per_maf_by_ind(statsdir, ref, cov, ['RA Het'], ymin=0.4)


#####################################################################################3




# for ind in inds:
# 	for cov in covs:
# 		plot_concordance_per_maf_by_geno('stats_rmts_61_5_impnofilter_anno', 'EUR', cov, ind)
# 		# plot_concordance_per_maf_by_geno('stats_rmts_61_5_imp_anno', 'EUR', cov, ind)



# for ind in inds:
# 	plot_concordance_per_maf_by_geno('stats_rmts_61_5_impnofilter_anno', 'FULL', 1, ind, ymin=0.7)
# 	plot_concordance_per_maf_by_geno('stats_rmts_61_5_imp_anno', 'FULL', 1, ind, ymin=0.8)


# for ind in inds:
# 	plot_entropy_per_maf_by_geno(datadir,'stats_rmts_61_5_impnofilter_anno', 'FULL', 1, ind, ymax=3)
	# plot_entropy_per_maf_by_geno(datadir,'stats_rmts_61_5_imp_anno', 'FULL', 1, ind)




# for cov in covs:
# 	plot_lc_concordance_per_maf_by_ind(statsdir, ref, cov, ['RA Het'], ymin=0.5)
# 	plot_lc_concordance_per_maf_by_ind(statsdir, ref, cov, ['RR Hom','RA Het','AA Hom'], ymin=0.0)


# for cov in covs:
	# plot_concordance_per_maf_by_ind(statsdir, ref, cov, ['RA Het'], ymin=0.5)
	# plot_concordance_per_maf_by_ind(statsdir, ref, cov, ['RR Hom','RA Het','AA Hom'], ymin=0.0)



# for cov in covs:
	# plot_discordance_per_maf(statsdir, 'EUR', cov)

################################### metrics per cov #####################################################

# article
# plot_concordance_per_cov(statsdir, 'EUR', ['RR Hom','RA Het','AA Hom'], ymin=0.98)
# plot_concordance_per_cov(statsdir, 'EUR', ['RA Het'], ymin=0.8)


# plot_entropy_per_cov(datadir, statsdir,'EUR', [0,1,2], ymax=0.1)
# plot_entropy_per_cov(datadir, statsdir,'EUR', [1])

# plot_lc_concordance_per_cov(statsdir, 'EUR', ['RR Hom','RA Het','AA Hom'])
# plot_lc_concordance_per_cov(statsdir, 'EUR', ['RA Het'])


# plot_lc_entropy_per_cov(datadir, statsdir,'EUR', [0,1,2], type="entropy")
# plot_lc_entropy_per_cov(datadir, statsdir,'EUR', [1], type="entropy")


# plot_lc_entropy_per_cov(datadir, statsdir,'EUR', [0,1,2], type="entropynew")
# plot_lc_entropy_per_cov(datadir, statsdir,'EUR', [1], type="entropynew")
#
# for ind in inds:
# 	plot_concordance_per_cov_by_geno(statsdir, ref, ind, ymin=0.8)


# plot_nrds_per_cov(statsdir, 'EUR')

###################################### metrics per ref set ########################################################

# plot_concordance_per_refset_by_ind(1,  ['RR Hom','RA Het','AA Hom'], ymin=0.99)
# plot_concordance_per_refset_by_ind(1,  ['RA Het'], ymin=0.9)




# plot_gcts_allinds('stats_rmts_61_5_imp',1,"FULL")
# plot_gcts_allinds_by_overlap('stats_rmts_61_5_imp',1,"FULL")
# plot_nrds_allinds_compref(1)
# plot_nrds_allinds_compref_by_overlap(1)


#################################### metrics per geno ########################33

#nonfiltered
# plot_concordance_per_geno_by_ind('stats_rmts_61_5_impnofilter_anno', cov, 'FULL', ref_title="ancient+FULL", ymin=0.95)
# plot_concordance_per_geno_by_ind('stats_rmts_61_5_impnofilter_anno', cov, 'EUR', ref_title="ancient+EUR")
# plot_concordance_per_geno_by_ind('stats_rmts_ind_nofilter_anno', cov, 'EUR', ref_title="individual+EUR")


# plot_concordance_per_geno_by_ind('stats_rmts_61_5_imp_anno', cov, 'FULL', ref_title="ancient+FULL_filtered", ymin=0.99)


#
#  plot_nummarkers_per_cov(statsdir, ref)

# for the ones phased with 61+5 samples
individual=""
#
# for the ones phased alone # These rewrote the stats-files for EUR 1x 61+5 ones - but the json files are for the 61+5 ones
# individual="_individual"


# plot_nrds_gcts_stacked_allinds('stats_rmts_61_5_imp',1,"FULL")




# plot_nrds_by_subset(full_1x_dir)
# plot_nrds_ind(full_1x_dir)
# plot_gcts_ind(full_1x_dir)
# 	plot_nrds_gcts_stacked(statsdir,cov,ref)



# plot_nrds_gcts_stacked_comp2(full_1x_dir, eur_1x_dir,"FULL","EUR")



	# calc_concordance_from_gcts(statsdir, cov, ref)
	# calc_concordance_from_gcss(statsdir, cov, ref)

###### table plotting ########
# for ind in ["ans17", "sf12", "LBK", "Loschbour", "ne1"]:
# 	plot_table_ind("NRDs",ind, full_1x_dir,full_1x_dir, ref="FULL")


# filters = ["", "_filtered", "_thets"]
# for f in filters:
# 	plot_all_inds(f, "GCTs",datadir+eur_1x_dir,eur_1x_dir,ref="EUR")


### overall concordance stats ###


