import sys, os
from collections import defaultdict


def read_set(file_name):
  result = defaultdict(set)
  line_counter = 0
  
  for line in open(file_name):
    line_counter += 1
    splited = map(int, line.split(" "))
    result[(splited[0], splited[1])].add(int(splited[2]))
 
  size = sum(map(len, result.values()))
  if line_counter != size:
    print "Be aware collision in files", line_counter, size

  return result




map1 = read_set(sys.argv[1])
map2 = read_set(sys.argv[2])

out_dict = sys.argv[3]

true_file = open(os.path.join(out_dict, "true_diags.txt"), "w")
false_file = open(os.path.join(out_dict, "false_diags.txt"), "w")
missing_file = open(os.path.join(out_dict, "missing_diags.txt"), "w")

true_counter = 0
false_counter = 0

seen = set()
for key in map1.keys():
  for v1 in map1[key]:
    not_find = True
    for vt in range(v1 - 6, v1 + 6):
      if vt in map2[key] and (key, vt) not in seen:
        true_counter += 1
        not_find = False
        seen.add((key, vt))
        #map2[key].remove(vt)
        to_file = true_file
        break
    if not_find:
      false_counter += 1
      to_file = false_file
    
    to_file.write(str(key[0]) + " " + str(key[1]) + " " + str(v1) + "\n")
seen = set()      
missing_counter = 0
for key in map2.keys():
  for v1 in map2[key]:
    not_find = True
    for vt in range(v1 - 6, v1 + 6):
      if vt in map1[key] and (key, vt) not in seen:
        not_find = False
        seen.add((key, vt))
        break
    if not_find:
      missing_counter += 1
      missing_file.write(str(key[0]) + " " + str(key[1]) + " " + str(v1) + "\n")

true_file.close()
false_file.close()
missing_file.close()


print "true_counter", true_counter
print "false_counter", false_counter
print "missing_counter", missing_counter

