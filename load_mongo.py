from mongokit import Document,Connection
from libtbx import easy_pickle
import datetime
import os

connection=Connection("localhost",27017)

@connection.register
class Insert_pdb_data(Document):
  """feature through the whole PDB files"""
  __collection__='pdb_feature'
  __database__='pdb_database'
  structure={
     'pdb_code':basestring,
     'data_type':basestring,
     'space_group':basestring,
     'unit_cell':[float],
     'resolution':float,
     'number_of_atoms':int,
     'number_of_atoms_super_sphere':int,
     'ligand':{'ligand_name_and_num':dict
              }, 
     'altloc':{'alt_conf_frac':float,
               'alt_loc_dist':dict,
              },
     'symmetry_ss_bonds':list,
     'fraction_of_nonH_incomplete':float
  }
  
workdir="/home/yanting/feature_set_pickle/"
for pickle_files in os.listdir(workdir):
  insert_pdb=connection.Insert_pdb_data()
  pickle_file=os.path.join(workdir,pickle_files)
  result=easy_pickle.load(pickle_file)
  try:
    insert_pdb['pdb_code']=pickle_files[:4]
    insert_pdb['data_type']=result.data_type
    insert_pdb['space_group']=result.space_group_symbol
    insert_pdb['unit_cell']=result.unit_cell
    insert_pdb['resolution']=result.resolution
    insert_pdb['number_of_atoms']=result.number_of_atoms
    insert_pdb['number_of_atoms_super_sphere']=result.number_of_atoms_super_sphere
    insert_pdb['ligand']['ligand_name_and_num']=result.ligands.ligand_name_and_num
    insert_pdb['altloc']['alt_conf_frac']=result.occupancies.alt_conf_frac
    insert_pdb['altloc']['alt_loc_dist']=result.occupancies.alt_loc_dist
    insert_pdb['symmetry_ss_bonds']=result.symmetry_ss_bonds
    insert_pdb['fraction_of_nonH_incomplete']=result.fraction_of_nonH_incomplete
    insert_pdb.save()
  except Exception,e:
    print "Failed",e
    print pickle_file
    print "**********************************************************************"
