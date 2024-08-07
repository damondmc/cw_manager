import re

match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')

def readMemoryUsage(filePath):
    try:
        file = open(filePath, 'r+')
        lines = file.readlines()
        #Same thing as before but with time instead of template count
        for line in lines:
            if 'completion-loop' in line:
                strings = line.split(',')
                memory = re.findall(match_number, strings[3])[0]
                break
        file.close()
        return float(memory)
    except:
        print(filePath)

def readTemplateCount(filePath):
    try:
        file = open(filePath, 'r+')
        lines = file.readlines()
        #Same thing as before but with time instead of template count
        for line in lines:
            if 'Number of semicoherent templates' in line:
                strings = line.split('=')
                nTemp = re.findall(match_number, strings[1])[0]
                break
        file.close()
        return int(nTemp) 
    except:
        print(filePath)

def readRunTime(filePath):
    try:
        file = open(filePath, 'r+')
        lines = file.readlines()
        #Same thing as before but with time instead of template count
        for line in lines:
            if 'completion-loop' in line:
                strings = line.split(',')
                memory = re.findall(match_number, strings[1])[0]
                break
        file.close()
        return float(memory)
    except:
        print(filePath)

def readEstimatedUpperStrainLimit(filePath):
    try:
        file = open(filePath, 'r+')
        lines = file.readlines()
        if 'DONE' in lines[-1]:
            strings = lines[-2].split('=')
            h0 = re.findall(match_number, strings[2])[0]
        else:
            print('The prgroam is not finished.')
            strings = lines[-1].split('=')
            h0 = re.findall(match_number, strings[2])[0]
        file.close()
        return float(h0)
    except:
        print(filePath)
