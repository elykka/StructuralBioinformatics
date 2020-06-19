import os
from os.path import join, isfile
import re

# put path to edges files
edges = '/home/elykka/Documents/CMC/data/vhl/edges'

# put path where to save file with all paths
data_path = '/home/elykka/Documents/CMC/data'

list_names = [f for f in os.listdir(edges) if isfile(join(edges, f))]
list_names.sort(key=lambda f: int(re.sub('\D', '', f)))

protein_name = edges.split('/')[-2]

file = open(join(data_path, protein_name + '.txt'), 'w+')
for name in list_names:
    file.write(join(edges, name) + '\n')

file.close()