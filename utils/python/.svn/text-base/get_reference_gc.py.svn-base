from pyfasta import Fasta
f = Fasta('/xchip/cga2/aaron/static/hg19.fasta',key_fn=lambda key: key.split()[0])
csv = open("/xchip/cga2/aaron/copy_number/data/targets/hg19/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.csv","r")
output = open("target_gc_content.txt","w")
output.write("target\tgc\tlength\n")

header = csv.readline()
zeros = 0
for line in csv:
    sp = line.strip().split(",")
    # print len(sp)
    # print sorted(f.keys())
    sequence = f[sp[1]][int(sp[2]):int(sp[3])]
    gc = 0
    for s in sequence:
        if s == "G" or s == "g" or s == "C" or s == "c":
            gc += 1
    if len(sequence) < 1:
        print "target " + sp[0] + " is less than one"
        zeros += 1
    else:
        output.write(sp[0] + "\t" + str(float(gc)/float(len(sequence))) + "\t" + str(len(sequence)) + "\n")
print "zeros = " + str(zeros)

