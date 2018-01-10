from libtbx import easy_pickle
import datetime
import os
from pymongo import MongoClient

client=MongoClient("localhost",27017)
db=client.pdb_database
coll=db.pdb_feature


  
workdir="/home/yanting/feature_set_pickle/"
for pickle_files in os.listdir(workdir):
  pickle_file=os.path.join(workdir,pickle_files)
  result=easy_pickle.load(pickle_file)[0]
  try:
    data={}
    data['pdb_code']=pickle_files[:4]
    data['data_type']=result.data_type
    data['space_group']=result.space_group_symbol
    data['unit_cell']=list(result.unit_cell)
    data['resolution']=result.resolution
    data['number_of_atoms']=result.number_of_atoms
    data['number_of_atoms_super_sphere']=result.number_of_atoms_super_sphere
    data['ligand_name_and_num']=result.ligands.ligand_name_and_num
    data['alt_conf_frac']=result.occupancies.alt_conf_frac
    data['alt_loc_dist']=result.occupancies.alt_loc_dist
    data['symmetry_ss_bonds']=result.symmetry_ss_bonds
    data['fraction_of_nonH_incomplete']=result.fraction_of_nonH_incomplete
    coll.insert(data)
  except Exception,e:
    print pickle_files[:4]
