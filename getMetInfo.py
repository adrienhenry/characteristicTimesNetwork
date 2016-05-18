import re
import urllib3
def getMetInfo(metName,pH=7.5,ionic=0.25,refConc=0.001):
    http = urllib3.PoolManager()
    isCompoundID = False
    if re.search("C\d{5}",metName) is not None :
         metIDs = [metName]
         isCompoundID = True
    else:
        r = http.request('GET', 'http://equilibrator.weizmann.ac.il/search?query=%s'%metName)
        r.status
        metIDs = unique(re.findall("compoundId=(C\d{5})",str(r.data)))
    met_results = []
    for met in metIDs:
        try:
            infos = get_KeggInfo(met)
            #import pdb;pdb.set_trace()
            if (re.sub(",|-|\(|\)","",metName).upper() in [re.sub(",|-|\(|\)","",nn).upper() for nn in infos["synonyms"]]) or isCompoundID:
                met_results = [infos]
                break
            else:
                met_results.append(infos)
                
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
    http = urllib3.PoolManager()
    r = http.request('GET', 'http://rest.kegg.jp/get/cpd:%s'%(compoundID))
    formula = re.search("FORMULA\s+(\w+)",str(r.data)).group(1)
    formula = re.findall("([A-Z](?:[a-z]+)?)(\d+)?",re.sub("</?sub>","",formula))
    temp = {}
    
    for xx in formula:
        if xx[1] is "":
            temp[xx[0]] = 1
        else:
            temp[xx[0]] = int(xx[1])
            
    formula = temp
    synonyms = re.search("(?:NAME)(.*?)\\\\n[A-Z]",str(r.data),).group(1)
    synonyms = re.sub("\s+|\\\\n","",synonyms).split(";")
    return({"formula":formula,"synonyms":synonyms,"energy":{"value":None,"error":None}})  
