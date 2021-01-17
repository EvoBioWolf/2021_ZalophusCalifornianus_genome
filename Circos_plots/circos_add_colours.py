fname = "final.bundle"
karyotype = "canFam3_karyotype.txt"
outname = "final_with_colours.txt"



d = {}
with open(karyotype, "r") as f:
    for line in f:
       (key, val) = line.split()[3], line.split()[6]
       d[(key)] = val

with open(fname, "r") as g:
	with open (outname, "w") as o:
		for line in g:
			outline = line.strip() + ",color=" + d[line.split()[0]]
			o.write(outline + "\n")
