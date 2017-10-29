import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import os 
import itertools
import numpy as np
from tqdm import tqdm
PATH=os.getcwd()
def Transition(PrevState,deletion):
    #PrevState['
    NextState={'RGT2':True,
               'YCK1_2':True,
               'GRR1':True,
               'STD1':True,
               'RGT1':True,
               'GLC7':True,
               'REG1':True,
               'SNF1':True,
               'SNF4':True,
               'MIG1':True,
               'MALT':True,
               'GAL2':True,
               'GAL11':True,
               'GAL80':True
    }
    # Nutrient input doesn't change, so update with previous value
    NutrientInputs=['glucose_ext','galactose_ext','maltose_ext']
    for N in NutrientInputs:
        NextState[N]=PrevState[N]
    # Nutrients and metabolites
    NextState['galactose_int'] = PrevState['galactose_ext'] and PrevState['Gal2p']
    NextState['maltose_int'] = ( PrevState['maltose_ext'] and ( ( ( PrevState['MalTp'] ) ) )    )
 
    # Regulatory genes
    NextState['CAT8'] =  not ( (PrevState['Mig1p']) )
    NextState['GAL1'] = ( (PrevState['Gal4p']) and not ( PrevState['Mig1p'] ) ) 
    NextState['GAL3'] =  not ( (PrevState['Mig1p']) )
    NextState['GAL4'] =  not ( (PrevState['Mig1p']) )
    NextState['MALR'] =  not ( (PrevState['Mig1p']) )
    NextState['MIG2'] =  not ( (PrevState['Rgt1p']) )
    NextState['MIG3'] =  not ( (PrevState['Rgt1p']) )
    NextState['MTH1'] =  not ( (PrevState['Mig1p'] and ( ( ( PrevState['Mig2p']) ) )    ) )
    NextState['SIP4'] = (PrevState['Cat8p']) 
    NextState['SNF3']= not (PrevState['Mig1p']) or not (PrevState['Mig2p'])  # Added 
    NextState['Yck1p'] = (PrevState['YCK1_2'] ) 
 
    # Output
    NextState['ACS1'] =(PrevState['Cat8p']) 
    NextState['FBP1'] =(PrevState['Cat8p'])  or(PrevState['Sip4p']) 
    NextState['GAL10'] =(PrevState['Gal4p']) 
    NextState['GAL5'] =(PrevState['Gal4p']) 
    NextState['GAL7'] =(PrevState['Gal4p']) 
    NextState['HXT1'] =  not ((PrevState['Rgt1p'] and (((PrevState['Mth1p']  or PrevState['Std1p'])))))
    NextState['HXT2'] =  not ((PrevState['Rgt1p'])  or(PrevState['Mig1p']) )
    NextState['HXT3'] =  not ((PrevState['Rgt1p'] and (((PrevState['Mth1p'])))))
    NextState['HXT4'] =  not ((PrevState['Mig1p'])  or ( PrevState['Rgt1p'] and (((PrevState['Mth1p'])))))
    NextState['HXT5'] =  not ((PrevState['Rgt1p']))
    NextState['HXT8'] =  not ((PrevState['Rgt1p']))
    NextState['ICL1'] =(PrevState['Cat8p'])  or(PrevState['Sip4p']) 
    NextState['IDP2'] =(PrevState['Cat8p']) 
    NextState['JEN1'] =(PrevState['Cat8p']) 
    NextState['MALS'] = ( (PrevState['MalRp']) and not (PrevState['Mig1p']) ) 
    NextState['MDH2'] =(PrevState['Cat8p'])  or(PrevState['Sip4p']) 
    NextState['MEL1'] = ((PrevState['Gal4p'])  ) or  not ( PrevState['Mig1p'] or PrevState['Gal4p'] ) 
    NextState['MLS1'] =(PrevState['Cat8p'])  or(PrevState['Sip4p']) 
    NextState['PCK1'] =(PrevState['Cat8p']) 
    NextState['SFC1'] =(PrevState['Cat8p']) 
    NextState['SUC2'] =  not ((PrevState['Mig1p'])  or(PrevState['Mig2p']) )
    NextState['CAT2'] =(PrevState['Cat8p'])
    NextState['4orfs'] = not ((PrevState['Rgt1p']) )
 
     # Proteins
    NextState['Cat8p'] = (PrevState['CAT8'] and ( ( (PrevState['Snf1p']) ) )   ) 
    NextState['Gal11p'] = (PrevState['GAL11']) 
    NextState['Gal1p'] = (PrevState['GAL1']) 
    NextState['Gal2p'] = (PrevState['GAL2']) 
    NextState['Gal3p'] = (PrevState['GAL3'] and ( ( (PrevState['galactose_int']) ) )   ) 
    NextState['Gal4p'] = ( (PrevState['GAL4']) and not (PrevState['Gal80p']) ) 
    NextState['Gal80p'] = ( ( (PrevState['GAL80']) and not (PrevState['Gal1p'])  ) and not (PrevState['Gal3p']) ) 
    NextState['Glc7Reg1'] = (PrevState['GLC7'] and ( ( (PrevState['REG1'] and PrevState['glucose_ext'] ) ) )   ) 
    NextState['MalRp'] = (PrevState['MALR']and ( ( (PrevState['maltose_int']) ) )   ) 
    NextState['MalTp'] = (PrevState['MALT']) 
    NextState['Mig1p'] = ( (PrevState['MIG1']) and not (PrevState['Snf1p']) ) 
    NextState['Mig2p'] = (PrevState['MIG2']) 
    NextState['Mig3p'] = ( (PrevState['MIG3']) and not (PrevState['Snf1p']) ) 
    NextState['Mth1p'] = (PrevState['MTH1'] and not PrevState['Yck1p'])  or (PrevState['MTH1'] and not PrevState['Rgt2p'] and not PrevState['Snf3p']) or (PrevState['MTH1'] and not PrevState['SCF_grr1']) #( ( ( ( (PrevState['MTH1']) and not (PrevState['SCF_grr1'])  ) and not (PrevState['Rgt2p'])  ) and not (PrevState['Yck1p'])  ) and not (PrevState['Snf3p']) ) 
    NextState['Rgt1p'] = (PrevState['RGT1'] and ( ( (PrevState['Mth1p'] or PrevState['Std1p']) ) )   ) 
    NextState['Rgt2p'] = (PrevState['glucose_ext'] and ( ( (PrevState['RGT2'])))) 
    NextState['SCF_grr1'] = (PrevState['GRR1']) 
    NextState['Sip4p'] = (PrevState['SIP4'] and ( ( (PrevState['Snf1p']) ) )   ) 
    NextState['Snf1p'] = ( (PrevState['SNF1'] and ( ( (PrevState['SNF4']) ) )    ) and not (PrevState['Glc7Reg1'])) 
    NextState['Snf3p'] = (PrevState['glucose_ext'] and ( ( (PrevState['SNF3'])))) 
    NextState['Std1p'] = (PrevState['STD1'] and not PrevState['Yck1p']) or (PrevState['STD1'] and not PrevState['SCF_grr1']) or (PrevState['STD1'] and not PrevState['Snf3p'] and not PrevState['Rgt2p'])#( ( ( ( (PrevState['STD1']) and not (PrevState['Rgt2p'])) and not (PrevState['SCF_grr1'] )  ) and not (PrevState['Yck1p'])  ) and not (PrevState['Snf3p'] ) )
    if deletion:
        for k in deletion.keys():
            NextState[k]=deletion[k]
    return NextState

def set_initial(modifier_dict):
    InitialValues={'glucose_ext':False,
                   'galactose_ext':False,
                   'maltose_ext':False,
                   'galactose_int':False,
                   'maltose_int':False,
                   'SNF3':False,
                   'RGT2':True,
                   'YCK1_2':True,
                   'GRR1':True,
                   'MTH1':False,
                   'STD1':True,
                   'RGT1':True,
                   'GLC7':True,
                   'REG1':True,
                   'SNF1':True,
                   'SNF4':True,
                   'MIG1':True,
                   'MIG2':False,
                   'MIG3':False,
                   'MALR':False,
                   'MALT':True,
                   'GAL1':False,
                   'GAL2':True,
                   'GAL3':False,
                   'GAL4':False,
                   'GAL11':True,
                   'GAL80':True,
                   'CAT8':False,
                   'SIP4':False,
                   'SUC2':False,
                   'HXT1':False,
                   'HXT2':False,
                   'HXT3':False,
                   'HXT4':False,
                   'HXT5':False,
                   'HXT8':False,
                   '4orfs':False,
                   'MALS':False,
                   'GAL5':False,
                   'GAL7':False,
                   'GAL10':False,
                   'MEL1':False,
                   'ICL1':False,
                   'FBP1':False,
                   'PCK1':False,
                   'MLS1':False,
                   'MDH2':False,
                   'ACS1':False,
                   'SFC1':False,
                   'CAT2':False,
                   'IDP2':False,
                   'JEN1':False,
                   'Snf3p':False,
                   'Rgt2p':False,
                   'Yck1p':False,
                   'SCF_grr1':False,
                   'Mth1p':False,
                   'Std1p':False,
                   'Rgt1p':False,
                   'Glc7Reg1':False,
                   'Snf1p':False,
                   'Mig1p':False,
                   'Mig2p':False,
                   'Mig3p':False,
                   'MalRp':False,
                   'MalTp':False,
                   'Gal1p':False,
                   'Gal2p':False,
                   'Gal3p':False,
                   'Gal4p':False,
                   'Gal11p':False,
                   'Gal80p':False,
                   'Cat8p':False,
                   'Sip4p':False,
            }
    NewInitialValues=InitialValues.copy()
    
    for k in modifier_dict.keys():
        if k not in NewInitialValues.keys():
            print(k, 'not found')
        NewInitialValues[k]=modifier_dict[k]
    return NewInitialValues

def ss_extractor(simulation):
    ss={}
    for k in simulation.keys():
        ss[k]=simulation[k][-1]
    return ss
def plot_as_heatmap(inputdict):
    df=pd.DataFrame.from_dict(inputdict)
    sns.heatmap(df)
    plt.show()

    
def LSS(Initial,NumIter,TransFunc,TransFuncArgs={}):
    States={'glucose_ext':[],'galactose_ext':[],'maltose_ext':[],'galactose_int':[],'maltose_int':[],'SNF3':[],'RGT2':[],'YCK1_2':[],'GRR1':[],'MTH1':[],'STD1':[],'RGT1':[],'GLC7':[],'REG1':[],'SNF1':[],'SNF4':[],'MIG1':[],'MIG2':[],'MIG3':[],'MALR':[],'MALT':[],'GAL1':[],'GAL2':[],'GAL3':[],'GAL4':[],'GAL11':[],'GAL80':[],'CAT8':[],'SIP4':[],'SUC2':[],'HXT1':[],'HXT2':[],'HXT3':[],'HXT4':[],'HXT5':[],'HXT8':[],'4orfs':[],'MALS':[],'GAL5':[],'GAL7':[],'GAL10':[],'MEL1':[],'ICL1':[],'FBP1':[],'PCK1':[],'MLS1':[],'MDH2':[],'ACS1':[],'SFC1':[],'CAT2':[],'IDP2':[],'JEN1':[],'Snf3p':[],'Rgt2p':[],'Yck1p':[],'SCF_grr1':[],'Mth1p':[],'Std1p':[],'Rgt1p':[],'Glc7Reg1':[],'Snf1p':[],'Mig1p':[],'Mig2p':[],'Mig3p':[],'MalRp':[],'MalTp':[],'Gal1p':[],'Gal2p':[],'Gal3p':[],'Gal4p':[],'Gal11p':[],'Gal80p':[],'Cat8p':[],'Sip4p':[]}
    
    #print('Initializing State Tracker...')
    for k in States.keys():
        States[k].append(Initial[k])
    #print('Successfully initialized State Tracker.') 

    CurrentState=Initial
    NextState={}
    if TransFuncArgs=={}:
        print('wt simulation')
    for i in range(0,NumIter):
        NextState=TransFunc(CurrentState,TransFuncArgs)
        for k in States.keys():
            States[k].append(NextState[k])
        CurrentState=NextState

    return States


NumIter=15
#StateTracker=LSS(set_initial({'glucose_ext':True}),NumIter,Transition)
#plot_as_heatmap(StateTracker)


# NutrientInputs=['glucose_ext','galactose_ext','maltose_ext']

# inputmodifierlist=[]
# for N in NutrientInputs:
#     inputmodifierlist.append({N:True})

# SS=[]
# for I in inputmodifierlist:
    
#     StateTracker=LSS(set_initial(I),NumIter,Transition)

#     SS.append(ss_extractor(StateTracker))

#     state_diff=[]

#     for i in range(0,NumIter-1):
#         counter=0
#         for k in StateTracker.keys():
#             if StateTracker[k][i]!= StateTracker[k][i+1]:
#             #print(k)
#                 counter=counter+1
#         state_diff.append(counter)
#     what_am_i_plotting=list(I.keys())
#     plt.plot(state_diff,label=what_am_i_plotting[0])
# plt.legend()
# plt.xlabel('Iteration number')
# plt.ylabel('# of changed states')
# plt.show()


FixedRegulators=['RGT2','YCK1_2','GRR1','STD1','RGT1','GLC7','REG1','SNF1','SNF4','MIG1','MALT','GAL2','GAL11','GAL80']
GeneDeletions={'wt':{},
    'cat8':{'CAT8':False},
    'gal1':{'GAL1':False},
    'gal11':{'GAL11':False},
    'gal2':{'GAL2':False},
    'gal3':{'GAL3':False},
    'gal4':{'GAL4':False},
    'gal80':{'GAL80':False},
    'glc7':{'GLC7':False},
    'grr1':{'GRR1':False},
    'malR':{'MALR':False},
    'malT':{'MALT':False},
    'mig1':{'MIG1':False},
    'mig1mig2':{'MIG1':False,'MIG2':False},
    'mig2':{'MIG2':False},
    'mig3':{'MIG3':False},
    'mth1':{'MTH1':False},
    'reg1':{'REG1':False},
    'rgt1':{'RGT1':False},
    'rgt2':{'RGT2':False},
    'rgt2snf3':{'RGT2':False,'SNF3':False},
    'sip4':{'SIP4':False},
    'snf1':{'SNF1':False},
    'snf1snf4':{'SNF1':False,'SNF4':False},
    'snf3':{'SNF3':False},
    'snf4':{'SNF4':False},
    'std1':{'STD1':False},
    'yck1_2':{'YCK1_2':False}
}
                     
# ###############################################################
# ###### Investigate wt-SS predictions for various nutrient inputs

# WTcomp=pd.read_csv(PATH+'/wt-responses.csv',sep='\t',index_col=0)

# WTcompdict=WTcomp.to_dict()
# MAP={'none':{},
#      'gal':{'galactose_ext':True},
#      'glc+mal':{'glucose_ext':True,
#                 'maltose_ext':True},
#      'gal+ mal':{'galactose_ext':True,
#                  'maltose_ext':True},
#      'mal':{'maltose_ext':True},
#      'all':{'maltose_ext':True,
#             'galactose_ext':True,
#             'glucose_ext':True},
#      'glc':{'glucose_ext':True},
#      'glc + gal':{'glucose_ext':True,
#                   'galactose_ext':True}}
# mismatch={}
# nutrient_spec_SS={}
# for k in MAP.keys():
#     StateTracker=LSS(set_initial(MAP[k]),NumIter,Transition)
#     nutrient_spec_SS[k]=ss_extractor(StateTracker)
#     counter=0
#     mismatch[k]=[]
#     for v in nutrient_spec_SS[k].keys():
#         if nutrient_spec_SS[k][v]!=WTcompdict[k][v]:
#             mismatch[k].append(v)
#             counter=counter+1
#     print("Mismatches in simulation", k, "=",counter)
# #####################################################################


###############################################################
###### Investigate KO-SS predictions for various nutrient inputs
# KOcomparison=pd.read_csv('../KO-responses.tsv',sep='\t',header=None)
# Delnames=[]
# for k in GeneDeletions.keys():
#     Delnames.append(k)


# varnamecol=list(KOcomparison[0])
# Genes=varnamecol[2:-1]

# KOdata={}

# # Construct dictionary from dataframe
# for c in range(1,225):
    
#     coldata=list(KOcomparison[c])
#     name=str(coldata[0])+':'+str(coldata[1])
#     expressionvalues=coldata[2:-1]
#     KOdata[name]={}
#     KOdata[name]['deletion']=str(coldata[0])
#     KOdata[name]['medium']=str(coldata[1])
#     KOdata[name]['expressionstates']={}
#     for i in range(0,len(expressionvalues)):
#         if expressionvalues[i]==0:
#             KOdata[name]['expressionstates'][Genes[i]]=False
#         elif expressionvalues[i]==1:
#             KOdata[name]['expressionstates'][Genes[i]]=True
            


# MediumMAP={'none':{},
#      'gal':{'galactose_ext':True},
#      'glc+mal':{'glucose_ext':True,
#                 'maltose_ext':True},
#      'gal+ mal':{'galactose_ext':True,
#                  'maltose_ext':True},
#      'mal':{'maltose_ext':True},
#      'all':{'maltose_ext':True,
#             'galactose_ext':True,
#             'glucose_ext':True},
#      'glc':{'glucose_ext':True},
#      'glc + gal':{'glucose_ext':True,
#                   'galactose_ext':True}}
# mismatch={}
# condition_spec_SS={}
# print("Del: Nut\t\t#Mismatch\n")

# for k in KOdata.keys():
#     StateTracker=LSS(set_initial(MediumMAP[KOdata[k]['medium']]),NumIter,Transition,GeneDeletions[KOdata[k]['deletion']])
#     condition_spec_SS[k]=ss_extractor(StateTracker)
#     conditioncounter=0
#     mismatch[k]=[]
#     for v in KOdata[k]['expressionstates'].keys():
#         if condition_spec_SS[k][v]!=KOdata[k]['expressionstates'][v]:
#             mismatch[k].append(v)
#             conditioncounter=conditioncounter+1
#     print('Mismatches for sim ',k,':',conditioncounter)

###############################################################

outputs=['ACS1','FBP1','GAL10','GAL5','GAL7','HXT1','HXT2','HXT3','HXT4','HXT5','HXT8','ICL1','IDP2','JEN1','MALS','MDH2','MEL1','MLS1','PCK1','SFC1','SUC2','CAT2','4orfs']

allgenes=['SNF3','RGT2','YCK1_2','GRR1','MTH1','STD1','RGT1','GLC7','REG1','SNF1','SNF4','MIG1','MIG2','MIG3','MALR','MALT','GAL1','GAL2','GAL3','GAL4','GAL11','GAL80','CAT8','SIP4','SUC2','HXT1','HXT2','HXT3','HXT4','HXT5','HXT8','4orfs','MALS','GAL5','GAL7','GAL10','MEL1','ICL1','FBP1','PCK1','MLS1','MDH2','ACS1','SFC1','CAT2','IDP2','JEN1','Snf3p','Rgt2p','Yck1p','SCF_grr1','Mth1p','Std1p','Rgt1p','Glc7Reg1','Snf1p','Mig1p','Mig2p','Mig3p','MalRp','MalTp','Gal1p','Gal2p','Gal3p','Gal4p','Gal11p','Gal80p','Cat8p','Sip4p']
regulatorygenes=[]

for i in range(0,len(allgenes)):
    if allgenes[i] not in outputs:
        regulatorygenes.append(allgenes[i])

AllSingleGeneDeletions={}
for i in range()
SingleGeneDeletions={#'wt':[{}},{
    'cat8':('CAT8',False,'1'),
    'gal1':('GAL1',False,'2'),
    'gal11':('GAL11',False,'3'),
    'gal2':('GAL2',False,'4'),
    'gal3':('GAL3',False,'5'),
    'gal4':('GAL4',False,'6'),
    'gal80':('GAL80',False,'7'),
    'glc7':('GLC7',False,'8'),
    'grr1':('GRR1',False,'9'),
    'malR':('MALR',False,'10'),
    'malT':('MALT',False,'11'),
    'mig1':('MIG1',False,'12'),
    'mig2':('MIG2',False,'13'),
    'mig3':('MIG3',False,'14'),
    'mth1':('MTH1',False,'15'),
    'reg1':('REG1',False,'16'),
    'rgt1':('RGT1',False,'17'),
    'rgt2':('RGT2',False,'18'),
    'sip4':('SIP4',False,'19'),
    'snf1':('SNF1',False,'20'),
    'snf3':('SNF3',False,'21'),
    'snf4':('SNF4',False,'22'),
    'std1':('STD1',False,'23'),
    'yck1_2':('YCK1_2',False,'24')
}

SGD=list(SingleGeneDeletions.keys())
Readout_states=[]
READOUT_VAR='SUC2' # SUC2 readout
delstringlist=[]
MutantName=[]
# {'total':,'viable':}
OrderViabilityStats={1:{'total':0,'viable':0},
                     2:{'total':0,'viable':0},
                     3:{'total':0,'viable':0},
                     4:{'total':0,'viable':0},
                     5:{'total':0,'viable':0}
}

count=1
for j in [1,2,3,4,5]:
    for mutant in tqdm(itertools.combinations(SGD,j)):
        READOUT={}
        delstring={}
        name=mutant[0] # initialize
        idstring=SingleGeneDeletions[mutant[0]][2]
        for i in range(1,len(mutant)):
            name=name+'|'+mutant[i]
            idstring=idstring+','+SingleGeneDeletions[mutant[i]][2]
        MutantName.append(name)
        for m in mutant:
            M=SingleGeneDeletions[m]
            delstring[SingleGeneDeletions[m][0]]=SingleGeneDeletions[m][1]
        delstringlist.append(delstring)
        
        StateTracker=LSS(set_initial({}),NumIter,Transition,delstring)# Initial condition is starvation
        SS=ss_extractor(StateTracker)
        READOUT['name']=name
        READOUT['value']=SS[READOUT_VAR]
        READOUT['id']=idstring
        
        OrderViabilityStats[j]['total']=OrderViabilityStats[j]['total']+1

        if READOUT['value']==True:
            OrderViabilityStats[j]['viable']=OrderViabilityStats[j]['viable']+1
        Readout_states.append(READOUT)
        count=count+1
            
glucRepressionCount=0
gluNonRepressionCount=0

for R in Readout_states:
    if R['value']==True:
        gluNonRepressionCount=gluNonRepressionCount+1
    else:
        glucRepressionCount=glucRepressionCount+1

with open('./mutant-predictions-starvation.txt','w') as out: 
    out.write('#id\tpet_file_id\tname\tspec\tviability\n')
    idnum=1
    for R in tqdm(Readout_states):
        if R['value']==False: # Conditions are reversed because of starvation simulation
            out.write(str(idnum)+'\t'+'NA'+'\t'+R['name']+'\t'+R['id']+'\t'+'inviable'+'\n')
        elif R['value']==True:
            out.write(str(idnum)+'\t'+'NA'+'\t'+R['name']+'\t'+R['id']+'\t'+'viable'+'\n')
        idnum=idnum+1
