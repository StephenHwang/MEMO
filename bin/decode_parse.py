import numpy as np

parse_file_path = 'HG01358.1_w_rc.fa.parse'
dict_file_path = 'dict.txt'                     # dict.txt made from "gzip -cdfq HG01358.1_w_rc.fa.dict > dict.txt"

# read parse
parse = np.fromfile(parse_file_path, dtype='uint32')

# read dict
#   - removing all odd chars (keeping '\x01' for splitting)
#   - view chars by set(dict_txt)
dict_txt = open(dict_file_path, 'r').read()
dict_txt = dict_txt.replace('\x00','').replace('\x02','')
dict_list = dict_txt.split('\x01')

# now reconstruct string (v0)
# phrases = [ dict_list[phrase_idx - 1] for phrase_idx in parse]
# print( ''.join(phrases) )

# now reconstruct string (v1)
print(dict_list[parse[0] - 1], end='')
for phrase_idx in parse[1:]:
    print(dict_list[phrase_idx - 1][10:], end='')
print()

