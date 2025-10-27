import numpy as np
import os

#set wd to file location
os.chdir(os.path.dirname(os.path.abspath(__file__)))

filename = "data/sample_finalreport.txt"
chunk_size = 500

data_chunks = []

with open(filename, 'r') as f:
    header = f.readline().strip().split('\t')

    while True:

        lines = [next(f, None) for _ in range(chunk_size)]
        lines = [line for line in lines if line is not None]
        
        if not lines:
            break

        chunk_array = np.array([line.strip().split('\t') for line in lines], dtype=str)
        
        
        data_chunks.append(chunk_array)


data = np.vstack(data_chunks)


print("Data shape:", data.shape)
print("First row:", data[0])
