import random

bases = ['A', 'T', 'G', 'C']
D = [100, 200, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000]

for d in D:
    seq = ""
    for _ in range(d):
        seq += random.choice(bases)
    with open(str(d)+"/"+str(d)+"_reference.fasta", 'w') as f:
        f.write(">"+str(d)+" reference\n")
        f.write(seq+"\n")
        
print "References ready, creating backbones:"
for d in D:
    # read reference sequence
    with open(str(d)+"/"+str(d)+"_reference.fasta", 'r') as f:
        seq = f.readlines()[1].strip()

    # change 15% of the sequence   
    indices = random.sample(range(len(seq)), int(0.15*len(seq)))
    for i in indices:
        new = seq[i]
        while new==seq[i]:
            new = random.choice(bases)
        seq = seq[:i]+new+seq[i+1:]
    # write new sequence as backbone 
    with open(str(d)+"/"+str(d)+"_backbone.fasta", 'w') as f:
        f.write(">"+str(d)+"_backbone\n")
        f.write(seq+"\n")
    print d, "done"

print "Backbones ready, creating reads:"

p = 0.1
FASTQ = " !\"#$%&'()*+,-./0123456789:;" #lower end of quality marks

for d in D:
    with open(str(d)+'/'+str(d)+"_reference.fasta", 'r') as f:
        seq = f.readlines()[1].strip()
    with open(str(d)+'/'+str(d)+"_reads.fastq", 'w') as f:
        pass
    
    reads = []
    L = int(p*len(seq))
    N = min(100, len(seq)-L-1)
    k = 1
    # pick N random places to start read, create L-length reads
    for i in random.sample(range(len(seq)-L), N):
        temp = seq[i:i+L]

        # random quality
        quality = ""
        for _ in range(len(temp)):
            quality += random.choice(FASTQ[len(FASTQ)/2:]) # higher quality to correct bases
            
        # 85% accuracy on read
        indices = random.sample(range(len(temp)), int(0.15*len(temp)))
        for j in indices:
            new = temp[j]
            while new==temp[j]:
                new = random.choice(bases)
            temp = temp[:j]+new+temp[j+1:]
            
            # lower quality to changed bases
            new_q = random.choice(FASTQ[:len(FASTQ)/2])
            quality = quality[:j]+new_q+quality[j+1:]
            
        with open(str(d)+'/'+str(d)+"_reads.fastq", 'a') as f:
            f.write('@'+str(k)+'\n'+temp+'\n'+'+'+'\n'+quality+'\n')
        k+=1
    print d, "done"

    
