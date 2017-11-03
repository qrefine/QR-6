from __future__ import division
import iotbx.pdb
import collections
from libtbx import group_args
from iotbx.pdb import remark_2_interpretation
from qrefine.super_cell import expand


def run(file_name):
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    #print dir(pdb_inp)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    """Metals (identity and counts), ions"""

    """X-ray data available"""
    print pdb_inp.get_experiment_type()

    """Extract_authors"""
    try:
        pdb_inp.extract_authors()
    except IndexError:
        print "None"
    else:
        print pdb_inp.extract_authors()

    """ Data resolution"""
    resolution = iotbx.pdb.remark_2_interpretation.extract_resolution(
      pdb_inp.extract_remark_iii_records(2))
    print resolution[0]

    """Crystal symmetry information"""
    symmetry = pdb_inp.crystal_symmetry()
    print symmetry.space_group_info()
    # symmetry = pdb_inp.crystal_symmetry_from_cryst1()
    # space_group = symmetry.space_group()
    # print space_group.type().lookup_symbol()
    # print symmetry.unit_cell().parameters()

    """Number of atoms(super sphere)"""
    # atom_count = pdb_inp.model_atom_counts()
    # print list(atom_count)
    # atom_count = pdb_inp.atoms()
    # print atom_count.size()

    super_cell = expand(
        pdb_hierarchy=pdb_hierarchy,
        crystal_symmetry=pdb_inp.crystal_symmetry()
    )
    ph_ss = super_cell.ph_super_sphere

    print pdb_hierarchy.atoms().size()
    print super_cell.ph_super_cell.atoms().size()
    print ph_ss.atoms().size()

    """Altlocs"""
    number_of_residues = 0
    number_of_alt_confs = 0
    alt_loc_dist = collections.Counter()
    for rg in pdb_hierarchy.residue_groups():
        number_of_residues += 1
        n_confs = len(rg.conformers())
        if (n_confs > 1):
            number_of_alt_confs += 1
        alt_loc_dist[n_confs] += 1
    print number_of_residues
    print number_of_alt_confs
    alt_conf_frac = number_of_alt_confs * 100. / number_of_residues
    print alt_conf_frac
    print "alt_loc_dist:", alt_loc_dist
    for key, value in zip(alt_loc_dist.keys(), alt_loc_dist.values()):
        print key, value
    return group_args(
      #atom_count=atom_count,
        alt_conf_frac = alt_conf_frac,
        alt_loc_dist  = alt_loc_dist)


if __name__ == '__main__':
    result = run(file_name="3dtj.pdb")
    print result
    print result.alt_conf_frac
    print result.alt_loc_dist