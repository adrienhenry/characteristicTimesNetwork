#! /usr/bin/python3
import os
import re
import libsbml
import math
from .UsefulFunctions import *
from .UnrenderedMathML import *
import numpy as np
import pandas as pd
from scipy import linalg
#from getMetInfo import *

""" == Metabolic Network =="""

class Network:
    def __init__(self,name):
        
        self.name = name
        self.metabolites = []
        self.reactions = []
        self.compartments = []        
        self.functionDefinitions = []
        self.parameters = []
        self.globalParameters = []
        self.derivatives = None
        self.rates = None
        self.JacobianPerturbation = None
        self.JacobianTracer = None
    def add_compartment(self,name,volume):
        comp = 0
        try:
            i = [xx.name for xx in self.compartments].index(name)
            comp = self.compartments[i]
        except ValueError:
            self.compartments.append(Compartment(name,volume))
            comp = self.compartments[-1]
            comp.whereami.append(self)
        return(comp)
            
    def addMetabolite(self,name,compartment="cell",quantity=None,concentration=None):
        comp = self.add_compartment(compartment,None)
        vol = comp.volume
        if quantity is None:
            quantity = concentration*vol        
        met = comp.get_metabolite(name)
        if met is None:
            met = Metabolite(name,comp,quantity)        
        assert (quantity == met.quantity) , "%s is already in the network with a different quantity."%met.name
        comp.metabolites.append(met)
        self.metabolites.append(met)        
        met.whereami.append(comp)
        met.whereami.append(self)
        return(met)

    def addReaction(self,name,metabolite_names,stoichiometry,kineticLaw=""):        
        index = search_index(name,[self.reactions[i] for i in range(len(self.reactions))])
        try:
            assert(index is None)
        except:
            print("%s is already in the network."%name)
        reac = Reaction(name)
        reac.stoichiometry = stoichiometry
        for mm in range(len(metabolite_names)):
            index = search_index(metabolite_names[mm],[met.name for met in self.metabolites])
            if index is None:
                try:
                    assert(len(self.compartments)==1)
                except:
                    print("%s is not in the network and there is ambiguity about the its compartment."%name)
                self.addMetabolite(metabolite_names[mm],self.compartments[0],quantity=None,concentration=1)
                index = -1            
            metabolite = self.metabolites[index]
            metabolite.reactions.append(reac)
            metabolite.stoichiometry.append(stoichiometry[mm])
            if stoichiometry[mm]<0:
                reac.substrates.append(metabolite)
            elif stoichiometry[mm]==0:
                reac.effectors.append(metabolite)
            elif stoichiometry[mm]>0:
                reac.products.append(metabolite)
            reac.metabolites.append(metabolite)
        
        reac.network = self
        reac.kineticLaw = kineticLaw
        self.reactions.append(reac)
        return(reac)
    
    def addFunctionDefinition(self,name,formula,arguments):        
        func = FunctionDefinition(name,formula,arguments)
        func.whereami.append(self)
        self.functionDefinitions.append(func)
        return(func)

    def addParameter(self,name,value,reaction=None):        
        param = Parameter(name,value,self)
        if reaction is not None:
            try:
                assert(search_index(name,[par.name for par in reaction.parameters]) is  None)
            except:
                print("Parameter %s is defined multiple times for reaction %s."%(name,reaction.name))
            reaction.parameters.append(param)            
            param.reactions.append(reaction)
        self.parameters.append(param)
        return(param)
    
    def readSBML(self,filename):
        document = libsbml.readSBML(filename)        
        if (document.getNumErrors() > 0):           
            document.printErrors();
            return(1);
        model = document.getModel();
        
        if (model == None):
            print("The model is empty." + "\n");
            return(1);

        if model.getNumCompartments()==0:

            print("There is no compartment in the model." + "\n");
            return(1);

        for cc in range(model.getNumCompartments()):
            try:
                self.add_compartment(model.getCompartment(cc).getId(),model.getCompartment(cc).getSize())
                
            except:
               
                print("Could not add the compartment number %d."%cc)
                return(1);

        if model.getNumSpecies()==0:
            print("There is no species in the model." + "\n");
            return(1);
        for mm in range(model.getNumSpecies()):
            try:
                met_name = model.getSpecies(mm).getId()
                isConstant = model.getSpecies(mm).getConstant()                
                compart = model.getSpecies(mm).getCompartment()
                quantity = model.getSpecies(mm).getInitialAmount()
                if quantity ==0 :
                    concentration = model.getSpecies(mm).getInitialConcentration()                    
                    met = self.addMetabolite(met_name,compart,concentration=concentration)
                else:
                    met = self.addMetabolite(met_name,compart,quantity=quantity)
                if model.getSpecies(mm).getName() is not "":
                    met.longName = model.getSpecies(mm).getName()
                met.isConstant = isConstant
                if isConstant:                    
                    met.exchangeCarbon = False                
                annotation = model.getSpecies(mm).getAnnotationString()
                compoundID = re.search("C\d{5}",annotation)
                chebiID = re.search("CHEBI.(\d{5})",annotation)
                if compoundID is not None:
                    compoundID = compoundID.group(0)
                    met.compoundID = compoundID
                if chebiID is not None:
                    met.chebiID = chebiID.group(1)
                    
                
            except:
                import pdb;pdb.set_trace()
                print("Could not add the metabolite number %d."%mm)
                return(1);

        for pp in range(model.getListOfParameters().size()):
            try:
                param = model.getListOfParameters().get(pp)
                pname = param.getName()
                pid = param.getId()
                if pname == "":
                        pname = pid
                param = self.addParameter(pname,param.getValue(),None)
                param.id = pid
                self.globalParameters.append(param)
            except:
                print("Could not add the global parameter number %d."%pp)
        if model.getNumReactions()==0:
            print("There is no species in the model." + "\n");
            return(1);

        this_dir, this_filename = os.path.split(__file__)
        DATA_PATH = os.path.join(this_dir, "knownUnrenderedMathML")
        knownUnrenderedMathML = []
        for line in open(DATA_PATH).read().split("\n"):
            if line is not "":
                knownUnrenderedMathML.append(UnrenderedMathML(line))            
        for ff in range(model.getNumFunctionDefinitions()):
            isComputable = True
            try:
                fd = model.getFunctionDefinition(ff)
                math = fd.getMath()
                arguments = [(math.getChild(n)).getName() for n in range ( math.getNumChildren()-1)] 
                math = math.getChild(math.getNumChildren() - 1)
                formula = libsbml.formulaToString(math)
                nameFd = fd.getId()
                for unknown in knownUnrenderedMathML:
                    if unknown.identify(formula):
                        formula = unknown.replace(formula)
                func = self.addFunctionDefinition(nameFd,formula,arguments)
                isComputable = bool(func.isComputable())
                # if not isComputable:
                #    # import pdb;pdb.set_trace()
                # assert(isComputable)
            except:
                if not bool(func.isComputable()):
                    print("Function number %s is not computable."%ff)
                print("Could not add the function definition number %d."%ff)
                return(1);
            
        for rr in range(model.getNumReactions()):            
            theSBMLreac = model.getReaction(rr)
            reac_name = theSBMLreac.getId()
            metabolite_names = []
            stoichiometry = []            
            try:
                for mm in range(theSBMLreac.getNumReactants()):
                    metabolite_names.append(theSBMLreac.getReactant(mm).getSpecies())
                    stoichiometry.append(-theSBMLreac.getReactant(mm).getStoichiometry())
                for mm in range(theSBMLreac.getNumModifiers()):
                    metabolite_names.append(theSBMLreac.getModifier(mm).getSpecies())
                    stoichiometry.append(0)
                for mm in range(theSBMLreac.getNumProducts()):
                    metabolite_names.append(theSBMLreac.getProduct(mm).getSpecies())
                    stoichiometry.append(theSBMLreac.getProduct(mm).getStoichiometry())
                kineticLaw = theSBMLreac.getKineticLaw()
                kineticLawMath = libsbml.formulaToString(kineticLaw.getMath())
                for unknown in knownUnrenderedMathML:
                    if unknown.identify(kineticLawMath):
                        kineticLawMath = unknown.replace(kineticLawMath)
                reac = self.addReaction(reac_name,metabolite_names,stoichiometry,kineticLawMath)
                reac.longName = theSBMLreac.getName()
                for pp in range(kineticLaw.getListOfParameters().size()):                    
                    param = kineticLaw.getListOfParameters().get(pp)
                    pname = param.getName()
                    if pname == "":
                        pname = param.getId()
                    self.addParameter(pname,param.getValue(),reac)
                
                
                funcDef = re.findall("(\w+)\([,\s\w]+?\)",kineticLawMath)
                for ff in funcDef:
                    func_result = self.functionDefinitions[getNames(self.functionDefinitions).index(ff)]
                    reac.functionDefinitions.append(func_result)
                    func_result.reactions.append(reac)
            except:
                print("Could not add the reaction number %d."%rr)
                return(1);
        
            
    def searchMetabolicInformations(self,withEnergy=False):
        ind = 0
        for met in self.metabolites:
                print("%d/%d\t%s(%s)"%(ind,len(self.metabolites)-1,met.longName,met.compoundID))
                ind+=1
                infos = None
                print(met.compoundID)

                if met.chebiID != "":
                    infos = getMetInfo(metID=met.chebiID,method="chebi")
                elif met.compoundID != "":
                    infos = getMetInfo(metID=met.compoundID,method="kegg")
                elif met.longName is not "":
                    infos = getMetInfo(met.longName)
                else:
                    infos = getMetInfo(name)
                if infos is None:
                    answer = None
                    while answer not in ["Y","N"]:
                        answer = input("Do you want to enter manually %s's informations ?(y/n)\n"%met.longName).upper()
                    if answer == "Y":
                        
                        if(withEnergy):
                            flag = True
                            while flag:                            
                                try:
                                    answer = input("Gibbs formation energy:\n")
                                    if answer == "None":
                                        met.gibbsEnergy = None
                                        met.energyError = None
                                    else:
                                        met.gibbsEnergy = float(answer)
                                        met.energyError = float(input("Energy error:\n"))
                                    flag = False
                                except:
                                    flag = True
                                    print("Energy information must be expressed by floats")
                        flag = True
                        while flag:
                            try:
                                formula = input("Metabolite's formula:\n")
                                formula = decomposeChemicalFormula(formula)
                                met.formula = formula
                                flag = False
                            except:
                                flag = True
                                print("The formula was not correct.")
                    continue
                else:
                    met.synonyms = infos["synonyms"]
                    met.gibbsEnergy = infos["energy"]["value"]
                    met.energyError = infos["energy"]["error"]
                    met.formula = infos["formula"]
                print(met.formula)
                print("Done")
    def testCarbonBalance(self):
        for reac in self.reactions:
            carb = np.array(reac.carbons())*reac.stoichiometry
            carb = np.array([met.isConstant is False for met in reac.metabolites])*carb
            reac.carbonBalance = np.sum(carb)
    def writeInformations(self,file):
        ff = open(file, 'w')
        ff.write("name,C,exchangeCarbon\n")
        for met in self.metabolites:            
            if "C" not in met.formula:
                met.formula = {"C":0}
            if met.exchangeCarbon is None:
                met.exchangeCarbon = True
            longName = met.longName
            longName = longName.replace(",","-")
            ff.write("%s,%s,%d,%r\n"%(met.name,longName,met.formula["C"],met.exchangeCarbon))
        ff.close()
    def readInformations(self,file):
        data = pd.read_csv(file,index_col=0)        
        for met in self.metabolites:
            try:
                met.formula = {"C":data.loc[met.name].C}
                met.exchangeCarbon = data.loc[met.name].exchangeCarbon
            except:
                print("The function could not find the informations for %s in the file %s."%(met.name,file))
    def updateParamIndexes(self):
        i = 0
        for param in self.parameters:            
            representation = "param[%d]"%i
            param.index = i
            param.representation = representation
            i += 1
        for comp in self.compartments:
            representation = "param[%d]"%i
            comp.index = i
            comp.representation = representation
            i += 1
    def updateConcIndexes(self):
        i=0
        for met in self.metabolites:
            representation = "conc[%d]"%i
            met.index = i
            met.representation = representation
            i += 1
    def updateReacIndexes(self):
        i = 0
        for reac in self.reactions:
            reac.editReaction()
            reac.index = i
            reac.representation = "v[%d]"%i
            reac.editedKineticLaw = "v[%d] = %s"%(i,reac.editedKineticLaw)
            i += 1
            
    def updateNetwork(self):
        self.updateParamIndexes()
        self.updateConcIndexes()
        self.updateReacIndexes()
        for met in self.metabolites:
            met.editMetDynamics()
            
    def separateForwardBackwardFluxes(self):
        functionsToSplit = [ff for ff in self.functionDefinitions]
        indd =0
        for ff in functionsToSplit:
            indd+=1
            oldName = ff.name
            oldFormula = ff.formula
            vp,vm = ff.splitForwardBackward()
            try:
                assert((vp is not None) and (vm is not None))
            except:
                import pdb;pdb.set_trace()
                print("The program did not manage to split separate function %s."%ff.name)
                continue
            if (vm is "") and (vp is not ""):
                ff.name = ff.name + "_p"
                ff.formula = vp
            if (vp is "") and (vm is not ""):
                ff.name = ff.name + "_m"
                ff.formula = vm
            if (vp is not "") and (vm is not ""):
                ff.name = ff.name + "_p"
                ff.formula = vp
                fm = self.addFunctionDefinition(oldName + "_m",vm,ff.arguments)
                ff.reverseFunction = fm
                fm.reverseFunction = ff
                ff.reactions = []
                
        reacToSplit = [rr for rr in self.reactions]
        for reac in reacToSplit:
            reac.functionDefinitions = []
            functions = findFunctions(reac.kineticLaw)
            newKinetic = reac.kineticLaw
            funcList = []
            for func in functions:                
                oldName = re.search(".*(?=\()",func).group(0)
                fragment = ""
                ffp = reac.searchFunction(oldName+"_p",warning=False)
                if ffp is not None:
                    fragment +=  re.sub("(?<!\w)%s(?!\w)"%re.escape(oldName),"%s_p"%oldName,func)
                    funcList.append(ffp)
                ffm = reac.searchFunction(oldName+"_m",warning=False)
                if ffm is not None:
                    fragment += "-" + re.sub("(?<!\w)%s(?!\w)"%re.escape(oldName),"%s_m"%oldName,func)
                    funcList.append(ffp)                
                newKinetic = re.sub("(?<!\w)%s(?!\w)"%re.escape(func),"(%s)"%fragment,newKinetic)
                
            nkp,nkm = findForwarBackward(newKinetic)            
            if (nkm is "") and (nkp is not ""):
                reac.name = reac.name + "_p"
                reac.kineticLaw = nkp
                reac.functionDefinitions = funcList
            if (nkp is "") and (nkm is not ""):
                reac.name = reac.name + "_m"
                reac.kineticLaw = nkm
                reac.functionDefinitions = funcList
            if (nkp is not "") and (nkm is not ""):                
                reac.name = reac.name + "_p"
                reac.kineticLaw = nkp
                reacm = self.addReaction(reac.name + "_m",getNames(reac.metabolites),[-1*xx for xx in reac.stoichiometry],nkm)
                reac.reverseReaction = reacm
                reacm.reverseReaction = reac
                for pp in reac.parameters:reacm.parameters.append(pp)
                for func in funcList:
                    if re.search(func.name,reac.kineticLaw) is not None:
                        reac.functionDefinitions.append(func)
                    if re.search(func.name,reacm.kineticLaw) is not None:
                        reacm.functionDefinitions.append(func)
                                                                    
    def generateDerivatives(self):
        funcText = "def derivatives(t,conc,param):"
        funcText += "\n\t" + "v = np.zeros(%d)"%len(self.reactions)
        
        for reac in self.reactions:
            funcText += "\n\t" + reac.editedKineticLaw
        funcText += "\n\tdconc = np.zeros(%d)"%len(self.metabolites)
        for met in self.metabolites:
            funcText += "\n\t" + met.editedKineticLaw
        funcText += "\n\treturn(dconc)"
        derivatives = None
        ns = {}
        try:
            exec(funcText,globals(),ns)
        except:import pdb;pdb.set_trace()
        self.derivatives = ns["derivatives"]
    def generateRates(self):
        funcText = "def rates(t,conc,param):"
        funcText += "\n\t" + "v = np.zeros(%d)"%len(self.reactions)
        for reac in self.reactions:
            funcText += "\n\t" + reac.editedKineticLaw
        funcText += "\n\treturn(v)"
        rates = None
        ns = {}
        exec(funcText,globals(),ns)    
        self.rates = ns["rates"]
        
    def getStoichiometryMatrix(self):
        S = np.zeros([len(self.metabolites),len(self.reactions)])
        for reac in self.reactions:
            for i in range(len(reac.metabolites)):
                met = reac.metabolites[i]
                S[met.index,reac.index] = reac.stoichiometry[i]
        return(S)


    def generateTracerJacobian_2(self,maxTime=10000,conc0=None,param=None,verbose=False): 
        if param is None:
            param = self.getParameterValues()
        if conc0 is None:
            conc0 = self.getMetaboliteValues()
            convergence_ok = False
            while not convergence_ok:
                convergence_ok,out_t,conc0 = findSS(conc0,param,self.derivatives,maxTime,verbose=False)
        concSS = conc0
        vSS = self.rates(out_t,conc0,param)
        try:
            assert(convergence_ok)
        except:
            print("The model did no manage to find a steady state in a simulation of %f."%maxTime)
            return(None)
        nMet = len(self.metabolites)
        J = np.zeros([nMet,nMet])
        
        for met in self.metabolites:
            k = []
            if met.exchangeCarbon:
                for reac in met.reactions:
                    if reac.carbonBalance!=0:
                        continue
                    cells = reac.getTracerMatrixComponents(met)                
                    if cells is None:
                        continue
                    k.append(cells)
                    for cc in cells:                   
                        J[met.index,cc[0]] += eval(cc[1])
            #if met.longName=="F6P":import pdb;pdb.set_trace()
        self.JacobianTracer = J
        return(J)
    
    def generateTracerJacobian(self,maxTime=10000,conc0=None,param=None,verbose=False,reaction=None): 
        if param is None:
            param = self.getParameterValues()
        if conc0 is None:
            conc0 = self.getMetaboliteValues()
            convergence_ok = False
            while not convergence_ok:
                convergence_ok,out_t,conc0 = findSS(conc0,param,self.derivatives,maxTime,verbose=False)
        concSS = conc0
        vSS = self.rates(out_t,conc0,param)
        try:
            assert(convergence_ok)
        except:
            print("The model did no manage to find a steady state in a simulation of %f."%maxTime)
            return(None)
        nMet = len(self.metabolites)
        J = np.zeros([nMet,nMet])
        if reaction is None:
            listReactions = [reac for reac in self.reactions]
        else:
            listReactions = [self.reactions[getNames(self.reactions).index(reaction)]]
        for reac in listReactions:
            carb = np.array(reac.carbons())
            carb = reac.stoichiometry*carb
            sub = np.array(reac.metabolites)[carb<0]
            prod = np.array(reac.metabolites)[carb>0]
            SSsub = np.array(reac.stoichiometry)[carb<0]
            SSprod = np.array(reac.stoichiometry)[carb>0]
            nbC = np.sum(carb[carb>0])
            if verbose:
                print("+".join(["%.1f%s(%d)"%(SSsub[i],sub[i].longName,sub[i].formula["C"]) for i in range(len(sub))]),end=" ")
                print("->",end=" ")
                print("+".join(["%.1f%s(%d)"%(SSprod[i],prod[i].longName,prod[i].formula["C"]) for i in range(len(prod))]))
            if np.sum(carb)>0:
                continue
            for i in range(len(sub)):
                ss = sub[i]
                if not ss.exchangeCarbon:
                    continue
                F = vSS[reac.index]/concSS[ss.index]
                J[ss.index,ss.index] -= F
                for j in range(len(prod)):
                    pp = prod[j]
                    if not pp.exchangeCarbon:
                        continue
                    J[pp.index,ss.index] += ss.formula["C"]*SSprod[j]*F/nbC
        self.JacobianTracer = J
        return(J)
            
    def generatePerturbationJacobian(self,maxTime=1e7,conc0=None,param=None,verbose=False):
        if param is None:
            param = self.getParameterValues()
        if conc0 is None:
            conc0 = self.getMetaboliteValues()
            convergence_ok = False
            while not convergence_ok:
                convergence_ok,out_t,conc0 = findSS(conc0,param,self.derivatives,maxTime,verbose=False)
        concSS = conc0
        vSS = self.rates(out_t,concSS,param)
        try:
            assert(convergence_ok)
        except:
            print("The model did no manage to find a steady state in a simulation of %f."%maxTime)
            return(None)
        nMet = len(self.metabolites)
        J = []
        dx = concSS*0.0001
        #for i in range(nMet):
        for j in range(nMet):
            conc = np.array(concSS)
            conc[j] += dx[j]
            derivp = self.derivatives(0,conc,param)
            conc = np.array(concSS)
            conc[j] -= dx[j]
            derivm = self.derivatives(0,conc,param)
            J.append((derivp-derivm)/(2*dx[j]))
        J = np.array(J)
        self.JacobianPerturbation = J
        return(J)
    
    def computeCharacteristicTimes(self,what,norm_type=1,met_to_perturb=[],method="inverseJacobian",verbose=False):
        if what == "tracer":
            J = np.array(self.JacobianTracer)
        elif what == "perturbation":
            J = np.array(self.JacobianPerturbation)
        else:
            print("There is no such system as %s."%what)
            return(None)
        nonNull = np.arange(len(J))[np.apply_along_axis(np.linalg.norm,1,J)>0]
        nonNull =np.arange(len(J))[self.derivatives(0,self.getMetaboliteValues(),self.getParameterValues())!=0]
        if met_to_perturb==[]:
            met_to_perturb = np.array(nonNull)
        J = J[nonNull,:][:,nonNull]
        
        ev = np.linalg.eig(J)[0]
        # if np.any(np.iscomplex(ev)) and method=="inverseJacobian":
        #     print("The method inverseJacobian is not suited for the problem the jacobian has complex eigen values.")
        #     return(None,None)        
        ev = np.real(ev)
        ev = ev[ev!=0]
        tau = -1/np.max(ev)
        Ts = None
        if method == "inverseJacobian":                        
            inverse = np.linalg.pinv(J)       
            Ts = np.apply_along_axis(lambda x:np.linalg.norm(x,1),1,inverse)
            return(tau,np.max(Ts))        
        elif method=="integration":
            dt =  tau
            mat_dt = linalg.expm(J*dt)            
            Ts = []
            for i in range(len(J)):
                dconc = np.zeros(len(J))
                dconc[i] = 1
                norm_num = [0]
                norm = [np.linalg.norm(dconc,norm_type)]
                tt = [0]
                count = 0
                while True:
                    tt.append(tt[-1] + dt)
                    dconc = np.dot(mat_dt,dconc)            
                    norm.append(np.linalg.norm(dconc,norm_type))
                    norm_num.append(0.5*(norm[-1]+norm[-2])*(tt[-1]-tt[-2]))
                    
                    if (norm[-1][norm[0]>0]/norm[0][norm[0]>0])< 1e-10:
                        break                    
                    elif verbose and count%100==0:
                        print((norm[-1]/norm[0]))
                    elif np.isnan((norm[-1]/norm[0])):
                        import pdb;pdb.set_trace()
                        print(norm[0])
                    count+=1
                Ts.append(np.sum(norm_num)/norm[0])
        T = np.max(Ts)
        return(tau,T)

    def getMetaboliteValues(self,type="concentration"):
        values = None
        if type == "concentration":
            values = np.array([met.quantity/met.compartment.volume for met in self.metabolites])
        elif type == "quantity":
            values = np.array([met.quantity for met in self.metabolites])
        return(values)
    def getParameterValues(self):
        values = np.append(np.array([param.value for param in self.parameters]),np.array([comp.volume for comp in self.compartments]))
        #print(np.array([comp.volume for comp in self.compartments]))
        return(values)

    
    
class Compartment:
    def __init__(self,name,volume=1):
        self.name = name
        self.unit="liter"
        self.volume = volume
        self.metabolites = []
        self.whereami = []
        self.index = None
        self.representation = None
    def get_metabolite(self,name):
        try:
            i = [xx.name for xx in self.metabolites].index(name)
            return(self.metabolites[i])
        except ValueError:
            return(None)

    def delete(self):
        while len(self.metabolites)>0:
            self.metabolites[0].delete()
        for parent in self.whereami:
            try:
                index = [xx.name for xx in parent.compartments].index(self.name)                
                del parent.compartments[index]
            except ValueError:
                print("Warning: %s has alread been removed from %s"%(self.name,parent.name))
                      
class Metabolite:
    def __init__(self,name,compartment,quantity):
        self.name = name
        self.longName = ""
        self.quantity = quantity
        self.compoundID = ""
        self.chebiID = ""
        self.gibbsEnergy = None
        self.energyError = None
        self.synonymous = []
        self.formula = {}
        self.compartment = compartment
        self.reactions = []
        self.tag = []
        self.stoichiometry = []
        self.whereami = []
        self.index = None
        self.representation = None
        self.editedKineticLaw = None
        self.exchangeCarbon = None
        self.isConstant = False
    def editMetDynamics(self):
        if self.isConstant:
           self.editedKineticLaw = "dconc[%d] = 0"%self.index
        else:
            reacContrib = []
            for rr in range(len(self.reactions)):
                if self.stoichiometry[rr]!=0:
                    reacContrib.append("+%.2f*%s"%(self.stoichiometry[rr],self.reactions[rr].representation))
            if len(reacContrib)==0:
                reacContrib = ["0"]
            reacContrib = "".join(reacContrib)
            reacContrib = re.sub("\+-"," - ",reacContrib)
            reacContrib = re.sub("\+"," + ",reacContrib)
            self.editedKineticLaw = "dconc[%d] = %s"%(self.index,reacContrib)

    def delete(self):
        while len(self.reactions)>0:
            self.reactions[0].delete()
        for parent in self.whereami:            
            parent.metabolites.remove(self)

class Reaction:
    def __init__(self,name):
        self.name = name
        self.longName = name
        self.substrates = []
        self.products = []
        self.effectors = []
        self.metabolites = []
        self.stoichiometry = []
        self.parameters = []
        self.network = None
        self.whereami = []
        self.kineticLaw = ""
        self.editedKineticLaw = ""
        self.reverseReaction = None
        self.functionDefinitions = []
        self.carbonBalance = 0
    def delete(self,deleteReverse = False):
        if self.reverseReaction is not None:
            if deleteReverse:
                self.reverseReaction.delete(False)
            else:
                self.reverseReaction.reverseReaction = None
        for met in self.metabolites:
            index = met.reactions.index(self)
            met.stoichiometry.pop(index)
            met.reactions.remove(self)
        for where in self.whereami:
            where.reactions.remove(self)
        
        self.network.reactions.remove(self)
        for ff in self.functionDefinitions:
            ff.reactions.remove(self)

    def carbons(self):
        try:
            return([met.formula["C"] for met in self.metabolites])
        except:
            for met in self.metabolites:
                if "C" not in met.formula:
                    print("There is no information about the number of carbons in %s."%met.name)
            return(None)
        
        
    def searchArgument(self,arg):
        for met in self.metabolites:
            if met.name == arg:
                return(met)
        for pp in self.parameters:
            if pp.name == arg or pp.id == arg:
                return(pp)
        for pp in  self.network.globalParameters:
            if pp.name == arg or pp.id == arg:
                return(pp)
        for cc in self.network.compartments:            
            if cc.name == arg:
                return(cc)
        return(None)
    def searchFunction(self,arg,warning=True):
        for func in self.network.functionDefinitions:
            if func.name == arg:
                return(func)
        if warning:
            print("No function named %s was found."%(arg))
        return(None)
    
    def editReaction(self):
        functions = findFunctions(self.kineticLaw)
        editedKineticLaw = self.kineticLaw
        count = 0
        for func in functions:
            # if self.name == 'vPFK':
            #     import pdb;pdb.set_trace()
            try:
                arguments = findArguments(func)                
                representations = []
                for pp in arguments:
                    temp = self.searchArgument(pp)
                    if temp is None:
                        representations.append(pp)
                    else:
                        representations.append(temp.representation)
                
                funcOb = re.search(".*(?=\()",func).group(0)
                funcOb = self.searchFunction(funcOb)
                fragment  = funcOb.editFunction(representations)
                editedKineticLaw = editedKineticLaw.replace(func,fragment)             
            except:
                print("The program could not edit the kinetic law for reaction %s."%(self.name))
        freeArguments = re.sub("conc\[\d+\]|param\[\d+\]|(?<!\w)\d+(?!\w)","",editedKineticLaw)
        freeArguments = re.findall("(?<!\w)\w+(?!\w)",freeArguments)
        for arg in freeArguments:
            argOb = self.searchArgument(arg)
            try:
                assert(arg is not None)                
            except:
                print("Editing %s was impossible because %s could not be found."%(self.name,arg))
                return(None)
            if argOb is None:
                import pdb;pdb.set_trace()
                print("The parameter %s is unknown."%arg)
            editedKineticLaw = re.sub("(?<!\w)%s(?!\w)"%re.escape(arg),argOb.representation,editedKineticLaw)
        self.editedKineticLaw = editedKineticLaw            
        
    def getTracerMatrixComponents(self,met0):
        mm = self.metabolites.index(met0)
        rr = self.index
        ## met0 substrate
        if self.stoichiometry[mm]<=0:
            nbC = np.array(self.carbons())[np.array(self.stoichiometry)<0]
            if len(nbC)>0:
                nbC = np.sum(nbC)
            else:
                return(None)
            stoich = "%.2f"%self.stoichiometry[mm]
            if stoich == "0.00":
                return(None)
            return([[met0.index,"-vSS[%d]/(concSS[%d])"%(self.index,met0.index)]])
        ## met0 product
        else:        
            nbC = self.carbons()*np.array(self.stoichiometry)
            nbC = nbC[nbC>0]            
            if len(nbC)>0:
                nbC = np.sum(nbC)
            else:
                return(None)
            toReturn = []
            front = self.stoichiometry[mm]/nbC
            for met in np.array(self.metabolites)[np.array(self.stoichiometry)<0]:
                n = met.formula["C"]
                if not met.exchangeCarbon:
                    continue                
                ratio = n*front
                toReturn.append([met.index,"%.15f*vSS[%d]/concSS[%d]"%(ratio,self.index,met.index)])
            return(toReturn)

class FunctionDefinition:  
    def __init__(self,name,formula,arguments):
        self.name = name
        self.formula = formula        
        self.arguments = arguments
        self.reactions = []
        self.whereami = []
        self.reverseFunction = None
    def editFunction(self,args):
        
        if len(args)!=len(self.arguments):
            print("The number of arguments submitted does not fit the formula's number of argument.")
            return(None)
        newFormula = self.formula
        for i in range(len(self.arguments)):            
            newFormula = re.sub("(?<!\w)%s(?!\w)"%self.arguments[i],args[i],newFormula)
        return "(%s)"%newFormula
    
    def isComputable(self,args=None,debug=False):
        newFormula = self.formula
        if args is None:
            args = np.random.rand(len(self.arguments))
        else:
            if len(args)!=len(self.arguments):
                print("The number of arguments submitted does not fit the formula's number of argument.")
                return(0)
        for i in range(len(self.arguments)):
            newFormula = re.sub("(?<!\w)%s(?!\w)"%self.arguments[i],str(args[i]),newFormula)
        if debug is True:
            print(newFormula)
        try:
            tempF = eval(newFormula)
            return(tempF)
        except:
            return(0)
    def splitForwardBackward(self,debug=False):
        vp,vm = findForwarBackward(self.formula)
        numericalVp = vp
        if numericalVp=="":numericalVp = "0"
        numericalVm = vm
        if numericalVm=="":numericalVm = "0"
        numericalF = self.formula
        args = np.random.lognormal(1, 1,len(self.arguments))
        for i in range(len(self.arguments)):
            numericalVp = re.sub("(?<!\w)%s(?!\w)"%self.arguments[i],str(args[i]),numericalVp)
            numericalVm = re.sub("(?<!\w)%s(?!\w)"%self.arguments[i],str(args[i]),numericalVm)
            numericalF = re.sub("(?<!\w)%s(?!\w)"%self.arguments[i],str(args[i]),numericalF)
        try:            
            numericalVp = eval(numericalVp)
            numericalVm = eval(numericalVm)
            numericalF = eval(numericalF)        
            np.testing.assert_almost_equal(numericalF , (numericalVp - numericalVm))
        except:
            if debug:
                print(vp)
                print(numericalVp)
                print(vm)
                print(numericalVm)
                print(numericalF)
            return(None,None)
        return(vp,vm)
    def delete(self):
        for ff in self.whereami:
            ff.functionDefinitions.remove(self)
        for rr in self.reactions:
            rr.functionDefinitions.remove(self)

class Parameter:
    def __init__(self,name,value,network):        
        self.name = name
        self.value = value
        self.id = ""
        self.reactions = []
        self.network = [network]
        self.tag = []
        self.index = None
        self.representation = None
    def delete(self):
        for rr in self.reactions:
            rr.parameters.remove(self)
        for nn in self.network:
            nn.parameters.remove(self)
            index = search_index(self,nn.globalParameters)
            if index is not None:
                nn.globalParameters.remove(self)
        

""" Useful small functions """

def search_index(elmt,the_list):
    try:
        index = the_list.index(elmt)
        return(index)
    except:
        return(None)
       

def unique(vec):
    temp = []
    for val in vec:
        if val not in temp:
            temp.append(val)
    return(temp)

def getNames(vec):
    return([xx.name for xx in vec])
