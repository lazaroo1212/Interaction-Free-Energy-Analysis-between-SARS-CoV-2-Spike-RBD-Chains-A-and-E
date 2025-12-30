file1 = "BioPhysics/Energy_analysis_project/NACCESS/NACCESS/7EKF_cmip.asa"
file2 = "BioPhysics/Energy_analysis_project/NACCESS/NACCESS/7EKF_cmip_chain_A.asa"
file3 = "BioPhysics/Energy_analysis_project/NACCESS/NACCESS/7EKF_cmip_chain_B.asa"
outfile = "BioPhysics/Energy_analysis_project/NACCESS/NACCESS/7EKF_cmip_interaction.asa"

def parse_cols(line):
	return [x for x in line.strip().split() if x]

with open(file1) as f1, open(file2) as f2, open(file3) as f3:
	lines1 = f1.readlines()
	lines2 = f2.readlines()
	lines3 = f3.readlines()

with open(outfile, "w") as out:
	# First: file1 - file2
	c = 0
	for line1, line2 in zip(lines1, lines2):
		cols1 = parse_cols(line1)
		cols2 = parse_cols(line2)

		diff = float(cols1[9]) - float(cols2[9])
		cols1[9] = f"{diff:.4f}"

		out.write("\t".join(cols1) + "\n")
		c+=1

	# Second: file1 - file3
	for line1, line3 in zip(lines1[c:], lines3):
		cols1 = parse_cols(line1)
		cols3 = parse_cols(line3)

		diff = float(cols1[9]) - float(cols3[9])
		cols1[9] = f"{diff:.4f}"

		out.write("\t".join(cols1) + "\n")
		c+=1
