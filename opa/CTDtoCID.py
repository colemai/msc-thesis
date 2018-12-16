#!/bin/bash/python3

import pdb
from sys import argv
import urllib.request
import pickle


def get_input_IDs (file_path):
	"""
	Input: Txt file path, one ID on each line
	Output: 
	"""
	with open (file_path, 'r') as file_object:
		lines = file_object.readlines()
		lines = [i.strip() for i in lines] # strip 
		return(lines)

def create_map_ctd_cid (ctd_ids):
	"""
	Input: LIST ctd_ids
	Output: DICT ctd ids as keys and cid ids as values
	"""
	map_ctd_cid = {}
	for item in ctd_ids:
		url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceid/Comparative%20Toxicogenomics%20Database/" + item + "/cids/TXT/"
		try:
			map_ctd_cid[item] = urllib.request.urlopen(url).read().strip()
		except Exception:
			print('Exception caught: ', Exception)
			continue  # or you could use 'continue'
	return map_ctd_cid
		
		

def save_obj(obj, name):
    with open('./'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
	input_IDs = get_input_IDs(argv[1])
	mapper = create_map_ctd_cid(input_IDs)
	save_obj(mapper, 'ctd_cid_map')
	print(map_ctd_cid)

