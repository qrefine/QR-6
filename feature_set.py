from __future__ import division
import iotbx.pdb
import collections
from libtbx import group_args
from iotbx.pdb import remark_2_interpretation
from qrefine.super_cell import expand

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
  #for key, value in zip(alt_loc_dist.keys(), alt_loc_dist.values()):
  #    print key, value
  return group_args(
    alt_conf_frac = alt_conf_frac,
    alt_loc_dist  = alt_loc_dist)
    
def get_non_standard_items(pdb_hierarchy):
  result = collections.Counter()
  if(sel.count(True)>0):
    ph = ph.select(sel)
    for rg in ph.residue_groups():
      for urn in rg.unique_resnames():
        result[urn]+=1
    print ",".join(["%s:%s"%(k,v) for k,v in zip(result.keys(), result.values())])

def run(file_name):
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    """Metals (identity and counts), ions"""

    """Crystal symmetry"""
    crystal_symmetry = pdb_inp.crystal_symmetry()
    """ Data resolution"""
    resolution = iotbx.pdb.remark_2_interpretation.extract_resolution(
      pdb_inp.extract_remark_iii_records(2))[0]
    """Number of atoms(super sphere)"""
    super_cell = expand(
        pdb_hierarchy=pdb_hierarchy,
        crystal_symmetry=pdb_inp.crystal_symmetry())
    """Altlocs"""
    result_occupancies = get_altloc_counts(pdb_hierarchy = pdb_hierarchy)
    return group_args(
      number_of_atoms              = pdb_hierarchy.atoms().size(),
      number_of_atoms_super_sphere = super_cell.ph_super_sphere.atoms().size(),
      occupancies        = result_occupancies,
      unit_cell          = crystal_symmetry.unit_cell().parameters(),
      space_group_symbol = crystal_symmetry.space_group().type().lookup_symbol(),
      resolution         = resolution,
      data_type          = pdb_inp.get_experiment_type())


if __name__ == '__main__':
    result = run(file_name="3dtj.pdb")
    print result
