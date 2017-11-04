from __future__ import division
import iotbx.pdb
import collections
from libtbx import group_args
from iotbx.pdb import remark_2_interpretation
from qrefine.super_cell import expand

"""Altlocs"""


def get_altloc_counts(pdb_hierarchy):
    number_of_residues = 0
    number_of_alt_confs = 0
    alt_loc_dist = collections.Counter()
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


def find_ss_across_symmetry(super_cell):
  ph_ss = super_cell.ph_super_sphere
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
      
      
  print "cys_master:", len(cys_master)
  print "cys_copies:", len(cys_copies)
  STOP()
          
          #print dir(ag)
          #STOP()
        #for conf in rg.conformers():
        #  print dir(conf)
        #  STOP()
          
#      for atom in chain.atoms().extract_element():
#        if atom =="S":
#            print dir(atom)
#
#    else:
#      print "general"
#
  
  STOP()


"""Metals (identity and counts), ions"""

def get_resolution(pdb_inp):
  resolution = None
  resolutions = iotbx.pdb.remark_2_interpretation.extract_resolution(
    pdb_inp.extract_remark_iii_records(2))
  if(resolutions is not None):
    resolution = resolutions[0]
  return resolution

def run(file_name):
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    """Crystal symmetry"""
    crystal_symmetry = pdb_inp.crystal_symmetry()
    """ Data resolution"""
    resolution = get_resolution(pdb_inp = pdb_inp)
    """Number of atoms(super sphere)"""
    super_cell = expand(
      pdb_hierarchy=pdb_hierarchy,
      crystal_symmetry=pdb_inp.crystal_symmetry())
    """Find SS across symmetry"""
    find_ss_across_symmetry(super_cell = super_cell)
    result_occupancies = get_altloc_counts(pdb_hierarchy=pdb_hierarchy)
    ligands = get_non_standard_items(pdb_hierarchy=pdb_hierarchy)
    return group_args(
      number_of_atoms=pdb_hierarchy.atoms().size(),
      number_of_atoms_super_sphere=super_cell.ph_super_sphere.atoms().size(),
      occupancies=result_occupancies,
      unit_cell=crystal_symmetry.unit_cell().parameters(),
      space_group_symbol=crystal_symmetry.space_group().type().lookup_symbol(),
      resolution=resolution,
      data_type=pdb_inp.get_experiment_type(),
      ligands=ligands)


if __name__ == '__main__':
    result = run(file_name="1f8t.pdb")
    print result
