import re

class UnrenderedMathML:
    def __init__(self,replacement=None,oldPattern=None,newPattern=None):
        if replacement is not None:
            patterns = re.sub("\s","",replacement).split("->")
        elif oldPattern is not None and newPattern is not None:
            patterns = [re.sub("\s","",oldPattern),re.sub("\s","",newPattern)]
        else:
            print("UnrenderedMathML requires either a replacement instuction or both the old and new patterns.")
            return(None)
        self.name = re.match("(\w+)\(",patterns[0]).group(1)
        self.oldPattern = patterns[0]
        self.newPattern = patterns[1]
        self.oldArgs = re.findall("#\d+",self.oldPattern)
        self.newArgs = re.findall("#\d+",self.newPattern)
        self.reTest = re.compile(re.sub("\\\\#\d+",".*?",re.escape(self.oldPattern)))
    def identify(self,formula):
        formula = re.sub("\s","",formula)
        return(len(self.reTest.findall(formula)) is not 0)
    def replace(self,formula):
        formula = re.sub("\s","",formula)
        numParamEntries = len(self.oldArgs)
        brackets = []
        toReturn = []
        argVec = []
        toReplace = []
        
        for ee in re.finditer("%s\("%self.name,formula):    
            started = False
            pos1 = ee.span()[1]
            pos2 = len(self.name)+1
            args = []
            char = 0
            startArg = pos1
            while pos2 is not len(self.oldPattern):
                if formula[pos1]==self.oldPattern[pos2]:                    
                    pos1 += 1
                    pos2 += 1
                else:
                    startArg=pos1
                    pos2 +=1
                    while self.oldPattern[pos2].isdigit():
                        pos2 += 1                   
                    level=0
                    
                    while formula[pos1]!=self.oldPattern[pos2] or level>0:
                        if formula[pos1]=="(":
                            level +=1
                        elif formula[pos1]==")":
                            level -= 1
                        pos1 += 1
                    args.append(formula[startArg:pos1])
            toReplace.append(formula[ee.span()[0]:pos1])
            argVec.append(args)
        for i in range(len(argVec)):
            newFormula = self.newPattern
            for j in range(len(argVec[i])):
                newFormula = newFormula.replace(self.oldArgs[j],argVec[i][j])
            formula = formula.replace(toReplace[i],newFormula)
        return(formula)
