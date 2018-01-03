# QR-6

### 1. Run script over entire PDB.
       phenix.python feature_set.py & (change the pdb_path & sf_path & nprocs)
##### get the result ------- pdb_filtered.pkl
##### load the pdb_filetered.pkl file get 80000+ pdb pickle file 
##### each pdb_code.pkl contain the information 
     [group_args
       data_type                      : X-RAY DIFFRACTION
       fraction_of_nonH_incomplete    : 0.0
       ligands                        : group_args
       ligand_name_and_num            : Counter()
       number_of_atoms                : 2288
       number_of_atoms_super_sphere   : 4248
       occupancies                    : group_args
       alt_conf_frac                  : 0.0
      alt_loc_dist                   : Counter()
      resolution                     : 4.0
      space_group_symbol             : P 1
      symmetry_ss_bonds              : []
      unit_cell                      : (51.603, 51.675, 51.797, 109.94, 108.48, 110.02)]
     
### 2. Populate a mongoDB based on pkl files.
     phenix.python load_mongo.py &
### 3. Loop over interesting systems and assert that qr.finalize works.

