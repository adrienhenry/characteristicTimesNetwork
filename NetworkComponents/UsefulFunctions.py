import urllib3
import re
import sympy
import numpy as np
import  bioservices
from scipy.integrate import ode

def getMetInfo(metName=None,metID=None,pH=7.5,ionic=0.25,refConc=0.001,method=None):
    http = urllib3.PoolManager()
    isCompoundID = False
    if metID is not None :
         metIDs = [metID]
         isCompoundID = True
    else:
        r = http.request('GET', 'http://equilibrator.weizmann.ac.il/search?query=%s'%metName)
        r.status
        metIDs = unique(re.findall("compoundId=(C\d{5})",str(r.data)))
    met_results = []
    for met in metIDs:
        try:
            if method is None or method is "kegg":
                infos = get_KeggInfo(met)
            elif method is "chebi":
                infos = get_chebiInfo(met)                
                
            if metName is not None:
                if ( re.sub(",|-|\(|\)","",metName).upper() in [re.sub(",|-|\(|\)","",nn).upper() for nn in infos["synonyms"]]):
                    pass
                else:
                    met_results.append(infos)
            elif metID is not None:
                met_results = [infos]
        except:
            continue
    ## Get ready to return a metabolite
    if len(met_results)==0:
        print("The program did not find any corresponding metabolite.")
        return(None)
    elif len(met_results)==1:        
        return(met_results[0])
    else:
        print("Multiple compounds matched the search:")
        while True:
            try:
                for i in range(len(met_results)):
                    print("[%d]\t%s"%(i,met_results[i]["synonyms"][0]))
                answer = input("Choose one or exit\n")
                
                if answer.upper() == "EXIT":
                    return(None)
                answer = int(answer)
                assert(answer>=0 and answer<len(met_results))
                
                met_results = met_results[answer]
                break
            except:
                print("You returned an invalid answer.")
                continue    
    return(met_results)

def unique(vec):
    temp = []
    for val in vec:
        if val not in temp:
            temp.append(val)
    return(temp)

    return({"formula":formula,"synonyms":None,"energy":{"value":None,"error":None}})

def get_eQuilibratorInfo(compoundID,pH=7.5,ionic=0.25,refConc=0.001):
    http = urllib3.PoolManager()
    r = http.request('GET', 'http://equilibrator.weizmann.ac.il/compound?compoundId=%s&ph=%.1f&ionic_strength=%.2f&reactantsConcentrationPrefactor=%.3f'%(compoundID,pH,ionic,refConc))
    formula = re.search("Formula.*?colspan=\"100%\">(.*?)(</sub>)?</td>",str(r.data)).group(1)
    formula = re.findall("([A-Z](?:[a-z]+)?)(\d+)?",re.sub("</?sub>","",formula))
    temp = {}
    for xx in formula:
        if xx[1] is "":
            temp[xx[0]] = 1
        else:
            temp[xx[0]] = int(xx[1])
            
    formula = temp        
    synonyms = re.search("Common names.*?<td colspan=\"100%\">(.*?)</td>",str(r.data)).group(1)
    synonyms = re.sub("\\\\t|\\\\n|\s","",synonyms).split(";")
            
    energy = re.search("G\\\\'<sup>m.*?<strong>(.*?)</strong>.*?plusmn; (\d+(\.\d+))",str(r.data))
    energy = {"value":float(energy.group(1)),"error":float(energy.group(2))}
    return({"formula":formula,"synonyms":synonyms,"energy":energy})       

def get_KeggInfo(compoundID):
    s = bioservices.KEGG()
    data = s.get("cpd:%s"%compoundID,parse=True)
    try:
        data = s.get("cpd:%s"%compoundID,parse=True)
        formula = data["FORMULA"] 
    except:
        
        print("get_chebi could not find the KEGG %s."%compoundID)
        return(None)   
    formula = data["FORMULA"]          
    formula = decomposeChemicalFormula(formula)
    synonyms = [xx.replace(";","") for xx in data["NAME"]]    
    return({"formula":formula,"synonyms":synonyms,"energy":{"value":None,"error":None}})  

def get_chebiInfo(chebiID):
    s = bioservices.ChEBI()
    try:
        formula = s.getCompleteEntity(chebiID).Formulae[0].data
    except:
        print("get_chebi could not find the ChEBI %s."%chebiID)
        return(None)        
    formula = decomposeChemicalFormula(formula)
    return({"formula":formula,"synonyms":None,"energy":{"value":None,"error":None}})  


def decomposeChemicalFormula(formula):
    formula = re.findall("([A-Z](?:[a-z]+)?)(\d+)?",re.sub("</?sub>","",formula))
    temp = {}
    for xx in formula:
        if xx[1] is "":
            temp[xx[0]] = 1
        else:
            temp[xx[0]] = int(xx[1])
    return(temp)

def findForwarBackward(formula):
    fragments = findBrackets(formula)
    
    for i in range(len(fragments)):
        #formula = re.sub("(\w)"+re.escape(fragments[0]),r"\1_fragmentNb%d_"%i,formula)
        formula = formula.replace(fragments[i],"_fragmentNb%d_"%i)
    variables = unique(re.findall("_?[A-Za-z]\w*",formula))
    command = ",".join(variables)+" = sympy.symbols(\"" + ",".join(variables) + "\")"
    exec(command)
    formula = sympy.expand(formula)
    
    backward = []
    forward = []
    additiveBlocks = sympy.srepr(formula)
    if additiveBlocks[0:3] == "Add":
        additiveBlocks = sympy.sympify(additiveBlocks[4:len(additiveBlocks)-1])
    else:
        additiveBlocks = [sympy.sympify(additiveBlocks)]
    for ee in additiveBlocks:
        ee = str(sympy.simplify(ee))
        ee = str(ee)
        if ee[0] == "-":
            ee = re.sub("^-","",ee)
            backward.append(ee)        
        else:
            forward.append(ee)
    backward = "+".join(backward)
    forward = "+".join(forward)
    for i in range(len(fragments)):
        backward = backward.replace("_fragmentNb%d_"%i,fragments[i])
        forward = forward.replace("_fragmentNb%d_"%i,fragments[i])

        
        # backward = re.sub("(?<!\w)fragmentNb%d(?!\w)"%i,fragments[i],backward)
        # forward = re.sub("(?<!\w)fragmentNb%d(?!\w)"%i,fragments[i],forward)
    return(forward,backward)

def findBrackets(formula):
    brackets = []
    toReturn = []
    test = []
    for ee in re.finditer("\(",formula):
        bracket = {"start":ee.span()[0],"end":None ,"text" : None }        
        brackets.append(bracket)
    lastEnding = -1
    for bracket in brackets:        
        level = 1                   
        position = bracket["start"]+1

        
        while level>0:
            if formula[position] == "(":
                level += 1
            elif formula[position] == ")":
                level -= 1
            position += 1                          
        bracket["end"] = position
        bracket["text"] = formula[bracket["start"]:bracket["end"]]
        if re.search("-",bracket["text"]) is not None:
            continue
        elif bracket["start"]<lastEnding:
            continue
        else:
            toReturn.append(bracket["text"])
            lastEnding = bracket["end"]
    return(toReturn)
    
def findSS(conc,param,deriv,time_max,verbose=False):
    convergence_ok = False
    r = ode(deriv).set_integrator('lsoda',rtol=1e-6,nsteps=1000)
    r.set_initial_value(conc, 0).set_f_params(param)
    dt = 0.001
    tmax = 5
    out = [np.array(conc)]
    out_deriv = [deriv(0,conc,param)]
    out_t = [0]
    i = 0
    while True:
        i+=1
        r.integrate(r.t+dt)
        out.append(r.y)        
        out_deriv.append(deriv(r.t,r.y,param))
        out_t.append(r.t)
        
        if r.t>time_max and time_max!=0:
            convergence_ok = False
            break
        if (np.all(np.abs(out_deriv[-1][out[-1]>0]/out[-1][out[-1]>0])<1e-5)):
            second_deriv = np.abs(out_deriv[-1])-np.abs(out_deriv[-2])
            if(np.max(second_deriv)<1e-8):
                convergence_ok = True
                break
            else:
                dt *= 2
        else:       
            dt *= 2
        if i%1==0 and verbose:
            second_deriv = np.abs(out_deriv[-1])-np.abs(out_deriv[-2])
            print("%.9f"%(np.max(np.abs(out_deriv[-1][out[-1]>0]/out[-1][out[-1]>0]))),"\t",np.max(second_deriv))
    t_tot = 0
    return(convergence_ok,out_t,np.array(out[-1]))


def findFunctions(formula):
    brackets = []
    toReturn = []
    test = []
    for ee in re.finditer("\w+\(",formula):
        bracket = {"start":ee.span()[0],"end":None ,"text" : None ,"inner":ee.span()[1],"level":0}     
        brackets.append(bracket)
    lastEnding = -1
    for bracket in brackets:
        level = 1            
        position = bracket["inner"]
        
        while level>0:
            if formula[position] == "(":
                level += 1
            elif formula[position] == ")":
                level -= 1
            position += 1        
        bracket["end"] = position
        bracket["text"] = formula[bracket["start"]:bracket["end"]]
    for bracket in brackets:
        level = 0
        for bb in brackets:                        
            if bracket["start"]>bb["start"] and bracket["end"]<bb["end"]:
                level += 1
        bracket["level"] = level
    sortby = np.argsort([bb["level"] for bb in brackets])[::-1]
    toReturn = [brackets[i]["text"] for i in sortby]  
    return(toReturn)
    
def findArguments(formula):
    formula = formula.replace(" ","")
    args = []
    formula = re.search("\((.*)\)",formula).group(1)
    level = 0
    start = 0
    for i in range(len(formula)):
        if formula[i] == "(":level+=1
        if formula[i] == ")":level-=1
        if formula[i] == "," and level==0:
            args.append(formula[start:i])
            start = i+1
    args.append(formula[start:len(formula)])
    return(args)
