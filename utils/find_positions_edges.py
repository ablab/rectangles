import sys

if len(sys.argv) != 3:
  print "<input nucmer file> <output>"
  exit(1)

CONTIG = "CONTIG"
UNALIGNED_CONTIG = "This contig is unaligned."
REAL_ALIGN = "Real Alignment"
TOTAL_ALIGN = "One align captures most of this contig"
SKIPPING = "Skipping"
MISASSEMBl = "Exte"
MARKING = "Marking as"
AMBIGUOUS = "Ambiguous Ali"

nucmer = open(sys.argv[1], "r")
out = open(sys.argv[2], "w")

infos = []
skipping = []
for line in nucmer:
  if CONTIG in line:     
    for line1 in skipping:
      out.write(line1)
      out.write("\n")
   
    skipping = []
    contig_name = line.split("_")[1]
  
  elif REAL_ALIGN in line or TOTAL_ALIGN in line or AMBIGUOUS in line:
    if contig_name != None:
      out.write('contig ' + contig_name + "\n")
      contig_name = None

    ##Example:  One align captures most of this contig: 2660551 2724138 | 1 63588 | 63588 63588 | 100.0 | gi|48994873|gb|U00096.2|_Escherichia_coli_str._K-12_substr._MG1655__complete_genome NODE_235_length_63588_cov_218.868
    lst = line.split(":")[1].split("|")
    gen_align = lst[0].strip().split(" ")
    contig_align = lst[1].strip().split(" ")
    gen_align_begin = int(gen_align[0])
    gen_align_end = int(gen_align[1])
    contig_align_begin = int(contig_align[0])
    contig_align_end = int(contig_align[1])
    
    out.write(lst[0].strip() + " "  + lst[1].strip())
    out.write("\n")
    
  elif SKIPPING in line:
    ##Example:      Skipping [4612278][4612357] redundant alignment 2 4612278 4612357 | 34 113 | 80 80 | 97.5 | gi|48994873|gb|U00096.2|_Escherichia_coli_str._K-12_substr._MG1655__complete_genome NODE_353_length_42569_cov_220.145     
    lst = line.split("|")
    tmp = lst[0].strip().split(" ")
    skip_line = "Skip " + tmp[-2] + " " + tmp[-1] + lst[1] 
    skipping.append(skip_line)
  elif MARKING in line:
    lst = line.split(":")[1].split("|")
    skipping.append("Marking as ambiguous " + lst[0].strip() + " "  + lst[1].strip())
out.close()

