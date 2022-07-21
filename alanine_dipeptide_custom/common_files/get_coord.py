import mdtraj
import numpy
import os.path

#dist_parent = mdtraj.compute_distances(parent, [[0,1]], periodic=True)
#dist_traj = mdtraj.compute_distances(traj, [[0,1]], periodic=True)
#dist = numpy.append(dist_parent,dist_traj)
#d_arr = numpy.asarray(dist)
#d_arr = d_arr*10
#numpy.savetxt("dist.dat", d_arr)

# Topology and other paths
topology_path = os.path.expandvars('$WEST_SIM_ROOT/common_files/diala_nowat_ff14_hmr.prmtop')
traj_path = os.path.expandvars('$WEST_CURRENT_SEG_DATA_REF/seg.nc')
parent_path = os.path.expandvars('$WEST_PARENT_DATA_REF/seg.ncrst')
ref_path = os.path.expandvars('$WEST_SIM_ROOT/common_files/diala_nowat_eq2.pdb')
atom_slice = numpy.asarray(range(0, 22))
ref_slice = atom_slice

parent_traj = mdtraj.load(parent_path, top=topology_path, atom_indices=atom_slice)
seg_traj = mdtraj.load(traj_path, top=topology_path, atom_indices=atom_slice)
ref_file = mdtraj.load(ref_path, top=topology_path, atom_indices=atom_slice)

# Start loading and stuff
full_traj = parent_traj.join(seg_traj)
#full_traj = full_traj.atom_slice(atom_indices=atom_slice)
full_traj = full_traj.superpose(ref_file, atom_indices=ref_slice)
all_coords = full_traj._xyz * 10

print(all_coords.shape)
numpy.save('coord.npy',all_coords)

#full_traj.save('coord.pdb')

#numpy.savetxt('coord.dat', [all_coords]) 
