import numpy as np
import sys

dft_file = "dft_energy_simple.txt"
nn_file = "nn_energy_simple.txt"
if len(sys.argv) == 3:
	dft_file, nn_file = sys.argv[1:]
# use 'awk ${print $5}' to get energy from grep-result
dft_data = []
nn_data = []

def get_rmse(nn_data, dft_data):
	return np.sqrt(np.mean((nn_data-dft_data) ** 2))

with open(dft_file, 'r') as fo:
	for line in fo:
		line_num = float(line.strip())
		dft_data.append(line_num)
with open(nn_file, 'r') as fo:
	for line in fo:
		line_num = float(line.strip())
		nn_data.append(line_num)

dft_data = np.array(dft_data, dtype=float)
nn_data = np.array(nn_data, dtype=float)

rmse = get_rmse(nn_data, dft_data)
print(f"RMSE for NN/DFT is {rmse} eV")




