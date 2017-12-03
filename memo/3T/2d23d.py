import re
s = '''
{
 {i0^2 ke[0] + j0^2 ke[1] + 2 i0 j0 ke[2], 
  j0 (j1 ke[1] + i1 ke[2]) + i0 (i1 ke[0] + j1 ke[2]), 
  j0 (j2 ke[1] + i2 ke[2]) + i0 (i2 ke[0] + j2 ke[2]), 
  i0^2 ke[5] + j0^2 ke[8] + i0 j0 (ke[4] + ke[9]), 
  j0 (i1 ke[4] + j1 ke[8]) + i0 (i1 ke[5] + j1 ke[9]), 
  j0 (i2 ke[4] + j2 ke[8]) + i0 (i2 ke[5] + j2 ke[9]), 
  i0^2 ke[14] + j0^2 ke[19] + i0 j0 (ke[13] + ke[20]), 
  j0 (i1 ke[13] + j1 ke[19]) + i0 (i1 ke[14] + j1 ke[20]), 
  j0 (i2 ke[13] + j2 ke[19]) + i0 (i2 ke[14] + j2 ke[20])},
 {0, i1^2 ke[0] + j1^2 ke[1] + 2 i1 j1 ke[2], 
  j1 (j2 ke[1] + i2 ke[2]) + i1 (i2 ke[0] + j2 ke[2]), 
  i0 (j1 ke[4] + i1 ke[5]) + j0 (j1 ke[8] + i1 ke[9]), 
  i1^2 ke[5] + j1^2 ke[8] + i1 j1 (ke[4] + ke[9]), 
  j1 (i2 ke[4] + j2 ke[8]) + i1 (i2 ke[5] + j2 ke[9]), 
  i0 (j1 ke[13] + i1 ke[14]) + j0 (j1 ke[19] + i1 ke[20]), 
  i1^2 ke[14] + j1^2 ke[19] + i1 j1 (ke[13] + ke[20]), 
  j1 (i2 ke[13] + j2 ke[19]) + i1 (i2 ke[14] + j2 ke[20])},
 {0, 0, i2^2 ke[0] + j2^2 ke[1] + 2 i2 j2 ke[2], 
  i0 (j2 ke[4] + i2 ke[5]) + j0 (j2 ke[8] + i2 ke[9]), 
  i1 (j2 ke[4] + i2 ke[5]) + j1 (j2 ke[8] + i2 ke[9]), 
  i2^2 ke[5] + j2^2 ke[8] + i2 j2 (ke[4] + ke[9]), 
  i0 (j2 ke[13] + i2 ke[14]) + j0 (j2 ke[19] + i2 ke[20]), 
  i1 (j2 ke[13] + i2 ke[14]) + j1 (j2 ke[19] + i2 ke[20]), 
  i2^2 ke[14] + j2^2 ke[19] + i2 j2 (ke[13] + ke[20])},
 {0, 0, 0, i0^2 ke[3] + j0^2 ke[6] + 2 i0 j0 ke[7], 
  j0 (j1 ke[6] + i1 ke[7]) + i0 (i1 ke[3] + j1 ke[7]), 
  j0 (j2 ke[6] + i2 ke[7]) + i0 (i2 ke[3] + j2 ke[7]), 
  i0^2 ke[12] + j0^2 ke[17] + i0 j0 (ke[11] + ke[18]), 
  j0 (i1 ke[11] + j1 ke[17]) + i0 (i1 ke[12] + j1 ke[18]), 
  j0 (i2 ke[11] + j2 ke[17]) + i0 (i2 ke[12] + j2 ke[18])},
 {0, 0, 0, 0, i1^2 ke[3] + j1^2 ke[6] + 2 i1 j1 ke[7], 
  j1 (j2 ke[6] + i2 ke[7]) + i1 (i2 ke[3] + j2 ke[7]), 
  i0 (j1 ke[11] + i1 ke[12]) + j0 (j1 ke[17] + i1 ke[18]), 
  i1^2 ke[12] + j1^2 ke[17] + i1 j1 (ke[11] + ke[18]), 
  j1 (i2 ke[11] + j2 ke[17]) + i1 (i2 ke[12] + j2 ke[18])},
 {0, 0, 0, 0, 0, i2^2 ke[3] + j2^2 ke[6] + 2 i2 j2 ke[7], 
  i0 (j2 ke[11] + i2 ke[12]) + j0 (j2 ke[17] + i2 ke[18]), 
  i1 (j2 ke[11] + i2 ke[12]) + j1 (j2 ke[17] + i2 ke[18]), 
  i2^2 ke[12] + j2^2 ke[17] + i2 j2 (ke[11] + ke[18])},
 {0, 0, 0, 0, 0, 0, i0^2 ke[10] + j0^2 ke[15] + 2 i0 j0 ke[16], 
  j0 (j1 ke[15] + i1 ke[16]) + i0 (i1 ke[10] + j1 ke[16]), 
  j0 (j2 ke[15] + i2 ke[16]) + i0 (i2 ke[10] + j2 ke[16])},
 {0, 0, 0, 0, 0, 0, 0, i1^2 ke[10] + j1^2 ke[15] + 2 i1 j1 ke[16], 
  j1 (j2 ke[15] + i2 ke[16]) + i1 (i2 ke[10] + j2 ke[16])},
 {0, 0, 0, 0, 0, 0, 0, 0, i2^2 ke[10] + j2^2 ke[15] + 2 i2 j2 ke[16]}
}
'''

s = s.replace('\n', '').replace('{', '').replace('}', '')
s = re.sub(r'ke\[(\d+)\]', r'ke\1', s)

s = re.sub(r'(\w+) ([\(\w]+)', r'\1 * \2', s)
s = re.sub(r'(\w+) ([\(\w]+)', r'\1 * \2', s)
s = re.sub(r'([ij]\d)\^2', r'\1 * \1', s)

s = re.sub(r'([ij])(\d)', r'\1[\2]', s)
s = re.sub(r'ke(\d+)', r'ke[\1]', s)
items = s.split(',')

tab = []
for i in range(9):
    tab.append(items[9*i:9*(i+1)])

def trim(s):
    while s[0] == ' ':
        s = s[1:]
    return s

count = 0
for j in range(9):
    for i in range(j, -1, -1):
        print(f'    Matrix[{count}] = {trim(tab[i][j])};')
        count += 1

