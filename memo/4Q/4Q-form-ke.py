import re
raw = '''
{
 {CB[0]^2 + d33 CB[4]^2, d33 CB[0] CB[4] + v CB[0] CB[4], 
  CB[0] CB[1] + d33 CB[4] CB[5], d33 CB[1] CB[4] + v CB[0] CB[5], 
  CB[0] CB[2] + d33 CB[4] CB[6], d33 CB[2] CB[4] + v CB[0] CB[6], 
  CB[0] CB[3] + d33 CB[4] CB[7], d33 CB[3] CB[4] + v CB[0] CB[7]},
 {0, d33 CB[0]^2 + CB[4]^2, v CB[1] CB[4] + d33 CB[0] CB[5], 
  d33 CB[0] CB[1] + CB[4] CB[5], v CB[2] CB[4] + d33 CB[0] CB[6], 
  d33 CB[0] CB[2] + CB[4] CB[6], v CB[3] CB[4] + d33 CB[0] CB[7], 
  d33 CB[0] CB[3] + CB[4] CB[7]},
 {0, 0, CB[1]^2 + d33 CB[5]^2, d33 CB[1] CB[5] + v CB[1] CB[5], 
  CB[1] CB[2] + d33 CB[5] CB[6], d33 CB[2] CB[5] + v CB[1] CB[6], 
  CB[1] CB[3] + d33 CB[5] CB[7], d33 CB[3] CB[5] + v CB[1] CB[7]},
 {0, 0, 0, d33 CB[1]^2 + CB[5]^2, v CB[2] CB[5] + d33 CB[1] CB[6], 
  d33 CB[1] CB[2] + CB[5] CB[6], v CB[3] CB[5] + d33 CB[1] CB[7], 
  d33 CB[1] CB[3] + CB[5] CB[7]},
 {0, 0, 0, 0, CB[2]^2 + d33 CB[6]^2, d33 CB[2] CB[6] + v CB[2] CB[6], 
  CB[2] CB[3] + d33 CB[6] CB[7], d33 CB[3] CB[6] + v CB[2] CB[7]},
 {0, 0, 0, 0, 0, d33 CB[2]^2 + CB[6]^2, 
  v CB[3] CB[6] + d33 CB[2] CB[7], d33 CB[2] CB[3] + CB[6] CB[7]},
 {0, 0, 0, 0, 0, 0, CB[3]^2 + d33 CB[7]^2, 
  d33 CB[3] CB[7] + v CB[3] CB[7]},
 {0, 0, 0, 0, 0, 0, 0, d33 CB[3]^2 + CB[7]^2}
}
'''

raw = re.sub(r'[\n\{\}]', '', raw)
raw = re.sub('\s+', ' ', raw)
raw = re.sub('^ ', '', raw)
# raw = re.sub(r'\] \(', '] * (', raw)
raw = re.sub(r'CB\[(\d)\]\^2', r'CB[\1] * CB[\1]', raw)
raw = re.sub(r'\] C', '] * C', raw)
raw = re.sub(r'v C', r'v * C', raw)
raw = re.sub('d33 C', 'd33 * C', raw)
raw = re.sub('CB', 'B', raw)
raws = raw.split(', ')
# print(raws)
K = [raws[8*i:8*(i+1)] for i in range(12)]

count = 0
for j in range(8):
    for i in range(0, j+1):
        print('    ke[%d] += cof * (%s);'%(count, K[i][j]))
        count += 1
        # print(K[i][j], end='\t')
    # print()
    
