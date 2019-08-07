import sys, os

# The purpose of this script is to take in a directory of txt files from the UpSet tool (set directory), and output the number of variants in each of the 6 possibilities:
# 1. 6/6 overlap
# 2. 5/6 overlap
# 3. 4/6 overlap
# 4. 3/6 overlap
# 5. 2/6 overlap
# 6. 1/6 overlap


def CountLinesInFile(infilename):
	infile = open(infilename,'r')
	lines=infile.readlines()
	varCount=len(lines)
	return varCount


def ReadDirIntoArray(dirname):
	fileList=os.listdir(dirname)
	print(fileList)
	return fileList

# This function will return a set number from a filename.
# File names look like this:
#000011_Family5_Family6.txt   
def DefineSetNumberFromFilename(filename):
	cols=filename.split('_')
	vector=cols[0]
	OneCount = 0
	for i in range(len(vector)):
		if vector[i]=='1':
			OneCount+=1
	label = str(OneCount)
	return label

# Takes the list of files and prints out our dictionary
def ArrayOfFiles2Dict(fileList,dirName):
	# initialize a dictionary, where the keys are overlaps (6 - 6/6 overlaps, 5: 5/6 overlaps...etc)
	FileDict={'1':0,'2':0,'3':0,'4':0,'5':0,'6':0}
	for each in fileList:
		if each[-3:] != 'txt':
			continue
		varCount=CountLinesInFile('%s%s'%(dirName,each))
		label = DefineSetNumberFromFilename(each)
		FileDict[label]+=varCount
	return FileDict

def PrintDict(FileDict):
	for i in range(1,7):
		Intersect=str(i)
		print("%s\t%d"%(Intersect,FileDict[Intersect]))

def Main():
	print("DominantDamaging")
	dirName='/Users/philliprichmond/Dropbox/Grad_School/ALD/MEI/All_VarLevel_MEIid_DominantDamaging_INTERVENE_UPSET/sets/'
	fileList = ReadDirIntoArray('/Users/philliprichmond/Dropbox/Grad_School/ALD/MEI/All_VarLevel_MEIid_DominantDamaging_INTERVENE_UPSET/sets/')
	FileDict = ArrayOfFiles2Dict(fileList,dirName)
	PrintDict(FileDict)

	print("RecessiveDamaging")
	dirName='/Users/philliprichmond/Dropbox/Grad_School/ALD/MEI/All_VarLevel_MEIid_RecessiveDamaging_INTERVENE_UPSET/sets/'
	fileList = ReadDirIntoArray('/Users/philliprichmond/Dropbox/Grad_School/ALD/MEI/All_VarLevel_MEIid_RecessiveDamaging_INTERVENE_UPSET/sets/')
	FileDict = ArrayOfFiles2Dict(fileList,dirName)
	PrintDict(FileDict)

	print("DominantProtective")
	dirName='/Users/philliprichmond/Dropbox/Grad_School/ALD/MEI/All_VarLevel_MEIid_DominantProtective_INTERVENE_UPSET/sets/'	
	fileList = ReadDirIntoArray('/Users/philliprichmond/Dropbox/Grad_School/ALD/MEI/All_VarLevel_MEIid_DominantProtective_INTERVENE_UPSET/sets/')
	FileDict = ArrayOfFiles2Dict(fileList,dirName)
	PrintDict(FileDict)

	print("RecessiveProtective")
	dirName='/Users/philliprichmond/Dropbox/Grad_School/ALD/MEI/All_VarLevel_MEIid_RecessiveProtective_INTERVENE_UPSET/sets/'
	fileList = ReadDirIntoArray('/Users/philliprichmond/Dropbox/Grad_School/ALD/MEI/All_VarLevel_MEIid_RecessiveProtective_INTERVENE_UPSET/sets/')
	FileDict = ArrayOfFiles2Dict(fileList,dirName)
	PrintDict(FileDict)


	return

if __name__=="__main__":
	Main()




