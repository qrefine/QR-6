from __future__ import division
import iotbx.pdb
import collections
import math
import os
import pickle
from libtbx import group_args
from iotbx.pdb import remark_2_interpretation
from qrefine.super_cell import expand
from mmtbx.monomer_library import server
from mmtbx.building import extend_sidechains
from libtbx import easy_pickle
from libtbx.easy_mp import parallel_map
from libtbx import Auto

mon_lib_server = server.server()
pdb_dir = "/home/yanting/pdb/pdb/"
structure_factors_dir = "/home/yanting/pdb/structure_factors"

"""Altlocs"""

def get_altloc_counts(pdb_hierarchy):
  number_of_residues = 0
  number_of_alt_confs = 0
  alt_loc_dist = collections.Counter ( )
  for rg in pdb_hierarchy.residue_groups():
    number_of_residues += 1
    n_confs = len(rg.conformers())
    if (n_confs > 1):
      number_of_alt_confs += 1
      alt_loc_dist[n_confs] += 1
    alt_conf_frac = number_of_alt_confs * 100. / number_of_residues
    # for key, value in zip(alt_loc_dist.keys(), alt_loc_dist.values()):
    #    print key, value
    return group_args(
      alt_conf_frac=alt_conf_frac,
      alt_loc_dist=alt_loc_dist)

"""Ligand counts""" """Metals and ions"""

def get_non_standard_items(pdb_hierarchy):
  result = collections.Counter()
  get_class = iotbx.pdb.common_residue_names_get_class
  ignore = [
    "common_amino_acid",
    "modified_amino_acid",
    "common_rna_dna",
    "modified_rna_dna",
    "ccp4_mon_lib_rna_dna",
    "common_water"]

  for rg in pdb_hierarchy.residue_groups():
    for urn in rg.unique_resnames():
      if (not get_class(urn) in ignore):
        result[urn] += 1

    # print ",".join(["%s:%s"%(k,v) for k,v in zip(result.keys(), result.values())])
  return group_args(
      ligand_name_and_num=result
  )

"""S-S"""

def find_ss_across_symmetry(super_cell, ss_bond_length_cutoff=2.5):
  ph_ss = super_cell.ph_super_sphere
  result = []
  cys_master = []
  cys_copies = []
  def fill_it(container, chain):
    for rg in chain.residue_groups():
      for ag in rg.atom_groups():
        if(ag.resname == "CYS"):
          container.append(ag)
  for chain in ph_ss.chains():
    if(len(chain.id.strip()) == 1):
      fill_it(container = cys_master, chain = chain)
    else:
      fill_it(container = cys_copies, chain = chain)
  for master_ag in cys_master:
    for master_atom in master_ag.atoms():
      if(master_atom.element.strip().upper() == "S"):
        master_S = master_atom
    if(master_S is not None):
      for copies_ag in cys_copies:
        copy_S = None
        for copy_atom in copies_ag.atoms():
          if(copy_atom.element.strip().upper() == "S"):
            copy_S = copy_atom
        if(copy_S is not None):
          r1 = master_S.xyz
          r2 = copy_S.xyz
          dist = math.sqrt(
            (r1[0]-r2[0])**2 + (r1[1]-r2[1])**2 + (r1[2]-r2[2])**2)
          if(dist < ss_bond_length_cutoff):
            result.append(master_S.quote())
  return result

"""Metals (identity and counts), ions"""

def get_resolution(pdb_inp):
  resolution = None
  resolutions = iotbx.pdb.remark_2_interpretation.extract_resolution(
    pdb_inp.extract_remark_iii_records(2))
  if(resolutions is not None):
    resolution = resolutions[0]
  return resolution

"""Fraction of non-H atom-incomplete residues"""
def complete_model(pdb_hierarchy):
  number_of_residues = len(list(pdb_hierarchy.residue_groups()))
  n_changed = extend_sidechains.extend_protein_model(
    pdb_hierarchy,
    mon_lib_server,
    add_hydrogens=False)
  fraction_of_nonH_incomplete = n_changed * 100. /number_of_residues
  return fraction_of_nonH_incomplete

def run(file_name):
  pdb_inp = iotbx.pdb.input(file_name=file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  n_atoms = pdb_hierarchy.atoms().size()
  if(n_atoms > 10000): return None
  fraction_of_nonH_incomplete = complete_model(pdb_hierarchy=pdb_hierarchy)
  crystal_symmetry = pdb_inp.crystal_symmetry()
  resolution = get_resolution(pdb_inp = pdb_inp)
  super_cell = expand(
    pdb_hierarchy    = pdb_hierarchy,
    crystal_symmetry = pdb_inp.crystal_symmetry())
  symmetry_ss_bonds  = find_ss_across_symmetry(super_cell = super_cell)
  result_occupancies = get_altloc_counts(pdb_hierarchy=pdb_hierarchy)
  ligands = get_non_standard_items(pdb_hierarchy=pdb_hierarchy)
  return group_args(
    number_of_atoms             = pdb_hierarchy.atoms().size(),
    number_of_atoms_super_sphere= super_cell.ph_super_sphere.atoms().size(),
    occupancies                 = result_occupancies,
    unit_cell                   = crystal_symmetry.unit_cell().parameters(),
    space_group_symbol          = crystal_symmetry.space_group().type().lookup_symbol(),
    resolution                  = resolution,
    data_type                   = pdb_inp.get_experiment_type(),
    ligands                     = ligands,
    symmetry_ss_bonds           = symmetry_ss_bonds,
    fraction_of_nonH_incomplete = fraction_of_nonH_incomplete)
def dump_pickle(file_name):
  try:
      pdb_code = os.path.basename(file_name)[3:7]
      result = run(file_name)
      if(result is not None):
        easy_pickle.dump(pdb_code+".pkl", result)
  except KeyboardInterrupt:raise 
  except Exception, e:
      print "FAILED:",file_name 
      print str(e)
      print "-"*79


if __name__ == '__main__':
  # PDB model files
  path = "/net/cci/pdb_mirror/pdb/"
  of = open("".join([path,"INDEX"]),"r")
  files = ["".join([path,f]).strip() for f in of.readlines()]
  of.close()
  # PDB reflection data files (list of corresponding codes)
  dpath = "/net/cci/pdb_mirror/structure_factors/"
  of = open("".join([dpath,"INDEX"]),"r")
  dfiles = [
    os.path.basename("".join([path,f]).strip())[1:5] for f in of.readlines()]
  of.close()
  #
  for f in files:
    pdb_code = os.path.basename(f)[3:7]
    if(pdb_code in dfiles):
      try:
        result = run(file_name=f)
        if(result is not None):
          easy_pickle.dump(pdb_code+".pkl", result)
      except KeyboardInterrupt: raise
      except Exception, e:
        print "FAILED:", f
        print str(e)
        print "-"*79