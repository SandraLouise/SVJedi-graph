#!/usr/bin/env python3

import sys
import numpy as np

def main(infiles):
	''' Return contingency table '''
	dict_truth = {}
	dict_prediction = {}
	
	counter = 0

	#Truth file
	with open(infiles[0]) as truthfile:
		for line in truthfile:
			if not line.startswith('#'):
				chrom, start, __, ___, ____, _____, ______, info, _______, geno, *_ = line.rstrip('\n').split('\t')
				
				genotype = geno.split(':')[0]
				length = info.split('SVLEN=')[1].split(';')[0]
				if length == 'False' or length == '0':
					end = info.split(';END=')[1].split(';')[0]
					length = int(end) - int(start)
				
				else:
					length = abs(int(length))
				
				if genotype == '0/0' or genotype == "0|0" : genotype = '0'
				elif genotype == '0/1' or genotype == '0|1' or genotype == '1|0' : genotype = '1'
				elif genotype == '1/1' or genotype == "1|1" : genotype = '2'
				else: continue
				
				#add variant in truth dictionary
				dict_truth[chrom + '_' + start + '-' + str(length)] = genotype
	
	#Prediction file
	with open(infiles[1]) as predictionfile:
		for line in predictionfile:
			if not line.startswith('#') and len(line.split('\t'))>2:
				bis_chrom, bis_start, __, ___, ____, _____, ______, bis_info, _______, bis_genotype = line.rstrip('\n').split('\t')
				
				bis_length = bis_info.split('SVLEN=')[1].split(';')[0]
				if bis_length == 'False' or bis_length == '0':
					bis_end = bis_info.split('END=')[1].split(';')[0]
					bis_length = int(bis_end) - int(bis_start)
				else:
					bis_length = abs(int(bis_length))
				
				sv_id = bis_chrom + '_' + bis_start + '-' + str(bis_length)
				
				genotype_predicted = bis_genotype.split(':')[0]
				
				if genotype_predicted == '0/0' : genotype_predicted = '0'
				elif genotype_predicted == '0/1' : genotype_predicted = '1'
				elif genotype_predicted == '1/1' : genotype_predicted = '2'
				elif genotype_predicted == './.' : genotype_predicted = '3'
				
				#add variant i prediction dictionary
				dict_prediction[sv_id] =  genotype_predicted
	
	FP, FN = 0, 0
		
	a = np.array([0] * 12).reshape(3, 4)
	for truth_sv in list(dict_truth.keys()):
		if truth_sv in list(dict_prediction.keys()):
			a[int(dict_truth[truth_sv]), int(dict_prediction[truth_sv])] += 1
			
			if dict_prediction[truth_sv] == '3': FN += 1
			elif dict_truth[truth_sv] != dict_prediction[truth_sv]: 
				FP += 1
				print('FP = ' + truth_sv)	

	TP = sum([a[i][i] for i in range(3)])
	
	print('---------------------')
	print('Table of contingency:')
	print(a)
	print()	
	print('Genotyping accuracy: ' 	+ str(round(TP / (TP + FP) * 100, 1)))
	print('Genotyping rate: ' 		+ str(round((TP + FP) / (TP + FP + FN) * 100, 1)))
	print('Number of predicted SV: ' + str(len(list(dict_prediction.keys())) - FN))
	print('Number of unpredicted SV: ' + str(FN))
	print('Number of known SVs: ' + str(len(list(dict_truth.keys()))))


if __name__ == '__main__':
	main(sys.argv[1:])
