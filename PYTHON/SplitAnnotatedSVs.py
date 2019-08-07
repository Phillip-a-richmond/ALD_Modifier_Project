import sys, os, argparse

#####################################
# Author:   Phillip Richmond        #
# Contact:  prichmond@cmmt.ubc.ca   #
# Open source GNU licensing         #
#####################################

# This will be hard coded for ALD, given the ordering of samples below

##################
### Initialize ###
##################
def GetArgs():
	
	#Read in your arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("-I","--Infile",help="Input file of SVs from AnnotSV",required=True)
	parser.add_argument("-T","--Type",help="Type of SVs in file to guage format: ERDS, SMOOVE, MELT",required=True)
	args = parser.parse_args()
	return args	


# Given two genotypes and a family ID, write the line to one of 4 files:
# Family#_DominantProtective.tsv
# Family#_RecessiveProtective.tsv
# Family#_DominantDamaging.tsv
# Family#_RecessiveDamaging.tsv

def WriteFamilyDiscordance(family,cald_GT,amn_GT,line):
	DominantProtective_outfile = open("%s_DominantProtective.tsv"%family,'a')
	RecessiveProtective_outfile = open("%s_RecessiveProtective.tsv"%family,'a')
	DominantDamaging_outfile = open("%s_DominantDamaging.tsv"%family,'a')
	RecessiveDamaging_outfile = open("%s_RecessiveDamaging.tsv"%family,'a')
	Shared_outfile = open("%s_Shared.tsv"%family,'a')
	if cald_GT == amn_GT:
		Shared_outfile.write(line)
	elif cald_GT == '0/1' and amn_GT == '0/0':
		DominantDamaging_outfile.write(line)
	elif cald_GT == '1/1' and ( (amn_GT == '0/0') or (amn_GT == '0/1') ):
		RecessiveDamaging_outfile.write(line)
	elif amn_GT == '0/1' and cald_GT == '0/0':
		DominantProtective_outfile.write(line)
	elif amn_GT == '1/1' and ( (cald_GT == '0/0') or (cald_GT == '0/1') ):
		RecessiveProtective_outfile.write(line)

#AnnotSV ID	SV chrom	SV start	SV end	SV length	SV typeREF	ALT	FORMAT	ALD010_BWAmem	ALD011_BWAmem	ALD026_BWAmem	ALD027_BWAmem	ALD036_BWAmem	ALD041_BWAmem	ALD042_BWAmem	ALD048_BWAmem	ALD049_BWAmem	ALD058_BWAmem	ALD059_BWAmem	ALD065_BWAmem	AnnotSV type	Gene name	NM	CDS length	tx length	location	intersectStart	intersectEnd	DGV_GAIN_IDs	DGV_GAIN_n_samples_with_SV	DGV_GAIN_n_samples_tested	DGV_GAIN_Frequency	DGV_LOSS_IDs	DGV_LOSS_n_samples_with_SV	DGV_LOSS_n_samples_tested	DGV_LOSS_Frequency	GD_ID	GD_AN	GD_N_HET	GD_N_HOMALT	GD_AF	GD_POPMAX_AF	GD_ID_others	DDD_SV	DDD_DUP_n_samples_with_SV	DDD_DUP_Frequency	DDD_DEL_n_samples_with_SV	DDD_DEL_Frequency	1000g_event	1000g_AF	1000g_max_AF	IMH_ID	IMH_AF	IMH_ID_others	promoters	dbVar_event	dbVar_variant	dbVar_status	TADcoordinates	ENCODEexperiments	GCcontent_left	GCcontent_rightRepeats_coord_left	Repeats_type_left	Repeats_coord_right	Repeats_type_right	ACMG	HI_CGscore	TriS_CGscore	DDD_status	DDD_modDDD_consequence	DDD_disease	DDD_pmids	HI_DDDpercent	synZ_ExAC	misZ_ExAC	pLI_ExAC	delZ_ExAC	dupZ_ExAC	cnvZ_ExAC	morbidGenes	morbidGenesCandidates	Mim Number	Phenotypes	Inheritance	AnnotSV ranking

#1_869471_870235_DEL	1	869471	870235	764	DEL	N	<DEL>	GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB:DHFC:DHFFC:DHBFC:DHSP	0/0:60:0:-0,-6,-20:20:20:0:20:0:7:0:0:13:0:0:0.3:0.290323:0.3:0	0/1:28:194.88:-20,-1,-4:16:7:8:7:8:0:0:0:7:8:0.53:0.0625:0.0645161:0.0588235:9	0/0:60:0:-0,-6,-20:21:21:0:20:0:5:0:0:15:0:0:0.322581:0.30303:0.30303:2	0/1:10:384.16:-40,-1,-2:23:7:15:7:15:0:0:0:7:15:0.68:0.125:0.133333:0.121212:10/0:51:0:-0,-5,-17:17:17:0:17:0:5:0:0:12:0:0:0.235294:0.235294:0.216216:1	1/1:34:826.3:-84,-5,-1:37:6:30:6:30:0:0:0:6:30:0.83:0:0:0:31	1/1:55:856.29:-86,-6,-1:34:3:30:3:30:0:0:0:3:30:0.91:0:0:0:31	1/1:9:531.38:-55,-2,-2:27:6:20:6:20:0:0:0:6:20:0.77:0:0:0:21	0/1:79:431.71:-44,-1,-9:37:18:18:18:18:0:0:0:18:18:0.5:0.155556:0.162791:0.152174:21	1/1:72:1190.8:-120,-8,-1:48:5:42:5:42:0:0:0:5:42:0.89:0:0:0:43	1/1:24:423.16:-43,-3,-1:18:2:15:2:15:0:0:0:2:15:0.88:0:0:0:16	0/1:13:357.02:-37,-1,-3:22:7:14:7:14:0:0:0:7:14:0.67:0.15625:0.172414:0.172414:15	full	SAMD11			-1		0	0	-1					-1	-1	gnomAD_v2_DUP_1_62			-1		-1	-1	-1		-1								0.715000	0.740000					75.37	-3.08373522872216	-3.77554068437523		-1.30231557829469	-2.53125239294469	-2.47162789907313			2

def ProduceFamilialDiscordantSVs(infilename,filetype):
	infile = open(infilename,'r')

	# Our headerline 
	headerline=infile.readline()
	# split into columns
	headercols = headerline.strip('\n').split('\t')
	print headercols[9:21]
	
	for line in infile:
		# I'll initialize a new variant dict, where the key is the column number and the 
		varDict = {}
		cols = line.strip('\n').split('\t')
		# Now I need to loop through the brother pairs and split the variants out based on their genotypes
		ALD010 = cols[9].split(':')[0]
		ALD011 = cols[10].split(':')[0]
		ALD026 = cols[11].split(':')[0]
		ALD027 = cols[12].split(':')[0]
		ALD036 = cols[13].split(':')[0]
		ALD041 = cols[14].split(':')[0]
		ALD042 = cols[15].split(':')[0]
		ALD048 = cols[16].split(':')[0]
		ALD049 = cols[17].split(':')[0]
		ALD058 = cols[18].split(':')[0]
		ALD059 = cols[19].split(':')[0]
		ALD065 = cols[20].split(':')[0]

		# Family 1
		WriteFamilyDiscordance('Family1',ALD010,ALD011,line)
		WriteFamilyDiscordance('Family2',ALD026,ALD027,line)
		WriteFamilyDiscordance('Family3',ALD036,ALD065,line)
		WriteFamilyDiscordance('Family4',ALD042,ALD041,line)
		WriteFamilyDiscordance('Family5',ALD049,ALD048,line)
		WriteFamilyDiscordance('Family6',ALD058,ALD059,line)

def Main():
	print "That'll do pig"
	ARGS = GetArgs()
	print "Working with this file: %s"%ARGS.Infile
	

	ProduceFamilialDiscordantSVs(ARGS.Infile,ARGS.Type)

if __name__=="__main__":
	Main()


