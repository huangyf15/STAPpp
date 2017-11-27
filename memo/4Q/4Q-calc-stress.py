import re

raw = '''
{CB[0] de[0] + CB[1] de[2] + CB[2] de[4] + CB[3] de[6] + 
  v (CB[4] de[1] + CB[5] de[3] + CB[6] de[5] + CB[7] de[7]), 
 CB[4] de[1] + CB[5] de[3] + CB[6] de[5] + 
  v (CB[0] de[0] + CB[1] de[2] + CB[2] de[4] + CB[3] de[6]) + 
  CB[7] de[7], 
 d33 (CB[4] de[0] + CB[0] de[1] + CB[5] de[2] + CB[1] de[3] + 
    CB[6] de[4] + CB[2] de[5] + CB[7] de[6] + CB[3] de[7])}
'''

raw = re.sub(r'[\n\{\}]', '', raw)
raw = re.sub('\s+', ' ', raw)
raw = re.sub('^ ', '', raw)
# raw = re.sub(r'\] \(', '] * (', raw)
raw = re.sub(r'CB\[(\d)\]\^2', r'CB[\1] * CB[\1]', raw)
raw = re.sub(r'\] (\w)', r'] * \1', raw)
raw = re.sub(r'v ', r'v * ', raw)
raw = re.sub(r'd33 ', r'd33 * ', raw)
raw = re.sub('CB', 'B', raw)
raws = raw.split(', ')
# print(raws)
# K = [raws[8*i:8*(i+1)] for i in range(12)]

for i in range(3):
    print('    stress[%d] = cof * (%s);'%(i, raws[i]))

# count = 0
# for j in range(8):
#     for i in range(0, j+1):
#         print('    ke[%d] += cof * (%s);'%(count, K[i][j]))
#         count += 1
