# from mpi4py import MPI #equivalent to the use of MPI_init() in C
# import matplotlib.pyplot as plt
# import numpy as np
# from pysemtools.io.ppymech.neksuite import pynekread, pynekwrite
# from pysemtools.datatypes.msh import Mesh
# from pysemtools.datatypes.field import FieldRegistry

# # Get mpi info
# comm = MPI.COMM_WORLD

# msh = Mesh(comm, create_connectivity=True)
# fld = FieldRegistry(comm)
# fname_mesh = '/home/sachinbm/cpu-neko-tf/neko-0.9.1/examples/turb_channel/field0.f01001'

# fname_1 = '/home/sachinbm/cpu-neko-tf/neko-0.9.1/examples/turb_channel/tau0.f01001'
# pynekread(fname_1, comm, data_dtype=np.double, msh=msh, fld=fld)

# print(fld.registry.keys())
# print(fld.registry['p'].shape)
# print(fld.t)

# print(fld.registry['p'])


# fname = lambda i: f'/home/sachinbm/cpu-neko-tf/neko-0.9.1/examples/turb_channel/field0.f{str(i).zfill(5)}'
# pynekread(fname_mesh, comm, data_dtype=np.double, msh=msh)
# for i in range(1001, 1010):
    # pynekread(fname(i), comm, data_dtype=np.double, fld=fld)

# --------------------------------------------------------------------------------------------------------------

from mpi4py import MPI
import matplotlib.pyplot as plt
import numpy as np
from pysemtools.io.ppymech.neksuite import pynekread, pynekwrite
from pysemtools.datatypes.msh import Mesh
from pysemtools.datatypes.field import FieldRegistry
from pysemtools.interpolation.mesh_to_mesh import PMapper

comm = MPI.COMM_WORLD
msh = Mesh(comm, create_connectivity=True)

# Function to generate filenames
fname_field = lambda i: f'/home/sachinbm/cpu-neko-tf/neko-develop/examples/turb_channel/field0.f{str(i).zfill(5)}'
fname_tau = lambda i: f'/home/sachinbm/cpu-neko-tf/neko-develop/examples/turb_channel/tau0.f{str(i).zfill(5)}'
fname_rew = lambda i: f'/home/sachinbm/cpu-neko-tf/neko-develop/examples/turb_channel/rl_data0.f{str(i).zfill(5)}'

# Read mesh from first field file
pynekread(fname_field(1001), comm, msh=msh)

# Collect tau pressure data and time values
tau_list = []
time_values = []

for i in range(1001, 1011):
    # Create a new FieldRegistry for each file to avoid conflicts
    fld_temp = FieldRegistry(comm)
    
    # Read the tau file into the temporary field registry
    pynekread(fname_tau(i), comm, msh=msh, fld=fld_temp)
    
    print(f"At time: {fld_temp.t}")
    
    # Store time and pressure data
    time_values.append(fld_temp.t)
    tau_p = fld_temp.registry['p'].copy()
    tau_list.append(tau_p)

# Convert to numpy arrays
tau_array = np.array(tau_list)
time_values = np.array(time_values)

# Calculate time intervals for weighted averaging
# For time-weighted averaging, we need to determine the time intervals
# Method 1: Use intervals between consecutive time points
time_intervals = np.diff(time_values)
# Add the first interval (assume same as second interval)
time_intervals = np.concatenate([[time_intervals[0]], time_intervals])

print(f"Time values: {time_values}")
print(f"Time intervals: {time_intervals}")

# Weighted time average
# Weight each field by its time interval
# Example: (10, 1000, 8, 8, 8) * (10, 1, 1, 1, 1) → (10, 1000, 8, 8, 8)
weighted_sum = np.sum(tau_array * time_intervals.reshape(-1, 1, 1, 1, 1), axis=0)
total_time = np.sum(time_intervals)
tau_weighted_mean = weighted_sum / total_time

# Compare with simple arithmetic mean
tau_simple_mean = np.mean(tau_array, axis=0)

print(f"\nWeighted time average range: [{tau_weighted_mean.min():.6f}, {tau_weighted_mean.max():.6f}]")
print(f"Simple arithmetic mean range: [{tau_simple_mean.min():.6f}, {tau_simple_mean.max():.6f}]")

# Check the difference
difference = tau_weighted_mean - tau_simple_mean
print(f"Max difference between weighted and simple average: {np.abs(difference).max():.8f}")

# Weighted standard deviation
weighted_mean_expanded = tau_weighted_mean.reshape(1, *tau_weighted_mean.shape)
weighted_variance = np.sum(time_intervals.reshape(-1, 1, 1, 1, 1) * 
                          (tau_array - weighted_mean_expanded)**2, axis=0) / total_time
tau_weighted_std = np.sqrt(weighted_variance)

print(f"Weighted standard deviation range: [{tau_weighted_std.min():.6f}, {tau_weighted_std.max():.6f}]")


# Create new FieldRegistry objects for the computed statistics
fld_mean = FieldRegistry(comm)
fld_std = FieldRegistry(comm)

# Add the computed mean field to the registry
# You need to add it as a field that can be written
fld_mean.add_field(comm, field_name='tau_mean', field=tau_weighted_mean)
# print(fld_mean.registry['tau_mean'])

# Add the computed std field to the registry  
fld_std.add_field(comm, field_name='tau_std', field=tau_weighted_std)

# Write the mean field to a file
output_mean_file = '/home/sachinbm/cpu-neko-tf/neko-develop/examples/turb_channel/statistics/tau_mean0.f00001'
pynekwrite(output_mean_file, comm, msh=msh, fld=fld_mean) # , istep=0, write_mesh=True

# Write the std field to a file
output_std_file = '/home/sachinbm/cpu-neko-tf/neko-develop/examples/turb_channel/statistics/tau_std0.f00001'
pynekwrite(output_std_file, comm, msh=msh, fld=fld_std) # , istep=0, write_mesh=True

print(f"Written mean field to: {output_mean_file}")
print(f"Written std field to: {output_std_file}")