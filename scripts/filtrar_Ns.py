import sys

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

def read_fasta(f):
    name, seq = None, []
    for line in open(f):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield name, "".join(seq)

total = kept = 0
with open(output_fasta, "w") as out:
    for h, s in read_fasta(input_fasta):
        total += 1
        if len(s) > 0 and s.upper().count("N") / len(s) < 0.20:
            kept += 1
            out.write(h + "\n")
            for i in range(0, len(s), 80):
                out.write(s[i:i+80] + "\n")

print(f"Retenidas: {kept}/{total}")
