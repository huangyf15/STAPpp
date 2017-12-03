import re

s = '''
{
 {-(-1 + v) x32^2 + 2 y23^2, (1 + v) x32 y23, 
  x13 (x32 - v x32) + 2 y23 y31, 
  2 v x13 y23 + x32 y31 - v x32 y31, -(-1 + v) x21 x32, 2 v x21 y23},
 {0, 2 x32^2 - (-1 + v) y23^2, x13 (y23 - v y23) + 2 v x32 y31, 
  2 x13 x32 - (-1 + v) y23 y31, -(-1 + v) x21 y23, 2 x21 x32},
 {0, 0, -(-1 + v) x13^2 + 2 y31^2, (1 + v) x13 y31, -(-1 + v) x13 x21,
   2 v x21 y31},
 {0, 0, 0, 2 x13^2 - (-1 + v) y31^2, -(-1 + v) x21 y31, 2 x13 x21},
 {0, 0, 0, 0, -(-1 + v) x21^2, 0},
 {0, 0, 0, 0, 0, 2 x21^2}
}
'''

s = s.replace('\n', '')
s = s.replace('{', '').replace('}', '')
items = s.split(',')

tab = [[] for i in range(6)]
for i in range(6):
    for j in range(6):
        item = items[i*6 + j]
        item = item.replace('-(-1 + v)', '+ va') # va = 1 - v
        item = item.replace('- (-1 + v)', '+ va') # va = 1 - v
        if item[:2] == ' +':
            item = item[2:]
        item = re.sub(r'([\w\)]) (\w)', r'\1 * \2', item)
        item = re.sub(r'([\w\)]) (\w)', r'\1 * \2', item)
        item = re.sub(r'(\w+)\^2', r'(\1 * \1)', item)
        tab[i].append(item)

count = 0
for j in range(6):
    for i in range(j, -1, -1):
        print(f'    ke[{count}] = cof * ({tab[i][j]});')
        count += 1
