import sys

if len(sys.argv) != 3:
  print ("<input fasta file with contigs> <output>")

fin = open(sys.argv[1])
fout = open(sys.argv[2], "w")
lens = []
line = fin.readline()
while line:
  all_seq = []
  seq = fin.readline()
  length = 0
  while seq and seq[0] != '>':
    all_seq.append(seq)
    length += len(seq.strip())
    seq = fin.readline()
  lens.append((length, (line, all_seq)))
  line = seq
lens.sort()
for item in lens:
  (length, (id_line, all_seq)) = item
  fout.write(id_line.strip() + " " + str(length) + "\n")
  for seq in all_seq:
    fout.write(seq)
fout.close()
