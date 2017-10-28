import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import os 
PATH=os.getcwd()
def Transition(PrevState):
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
    NextState['SNF3']= not (PrevState['Mig1p']) and not (PrevState['Mig2p'])  # Added 
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
    NextState['4orfs'] =  not ((PrevState['RGT1']) )
 
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
    NextState['Mth1p'] = ( ( ( ( (PrevState['MTH1']) and not (PrevState['SCF_grr1'])  ) and not (PrevState['Rgt2p'])  ) and not (PrevState['Yck1p'])  ) and not (PrevState['Snf3p']) ) 
    NextState['Rgt1p'] = (PrevState['RGT1'] and ( ( (PrevState['Mth1p'] or PrevState['Std1p']) ) )   ) 
    NextState['Rgt2p'] = (PrevState['glucose_ext'] and ( ( (PrevState['RGT2'])))) 
    NextState['SCF_grr1'] = (PrevState['GRR1']) 
    NextState['Sip4p'] = (PrevState['SIP4'] and ( ( (PrevState['Snf1p']) ) )   ) 
    NextState['Snf1p'] = ( (PrevState['SNF1'] and ( ( (PrevState['SNF4']) ) )    ) and not (PrevState['Glc7Reg1'])) 
    NextState['Snf3p'] = (PrevState['glucose_ext'] and ( ( (PrevState['SNF3'])))) 
    NextState['Std1p'] = ( ( ( ( (PrevState['STD1']) and not (PrevState['Rgt2p'])) and not (PrevState['SCF_grr1'] )  ) and not (PrevState['Yck1p'])  ) and not (PrevState['Snf3p'] ) )
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

    
def LSS(Initial,NumIter,TransFunc):
    States={'glucose_ext':[],'galactose_ext':[],'maltose_ext':[],'galactose_int':[],'maltose_int':[],'SNF3':[],'RGT2':[],'YCK1_2':[],'GRR1':[],'MTH1':[],'STD1':[],'RGT1':[],'GLC7':[],'REG1':[],'SNF1':[],'SNF4':[],'MIG1':[],'MIG2':[],'MIG3':[],'MALR':[],'MALT':[],'GAL1':[],'GAL2':[],'GAL3':[],'GAL4':[],'GAL11':[],'GAL80':[],'CAT8':[],'SIP4':[],'SUC2':[],'HXT1':[],'HXT2':[],'HXT3':[],'HXT4':[],'HXT5':[],'HXT8':[],'4orfs':[],'MALS':[],'GAL5':[],'GAL7':[],'GAL10':[],'MEL1':[],'ICL1':[],'FBP1':[],'PCK1':[],'MLS1':[],'MDH2':[],'ACS1':[],'SFC1':[],'CAT2':[],'IDP2':[],'JEN1':[],'Snf3p':[],'Rgt2p':[],'Yck1p':[],'SCF_grr1':[],'Mth1p':[],'Std1p':[],'Rgt1p':[],'Glc7Reg1':[],'Snf1p':[],'Mig1p':[],'Mig2p':[],'Mig3p':[],'MalRp':[],'MalTp':[],'Gal1p':[],'Gal2p':[],'Gal3p':[],'Gal4p':[],'Gal11p':[],'Gal80p':[],'Cat8p':[],'Sip4p':[]}
    
    print('Initializing State Tracker...')
    for k in States.keys():
        States[k].append(Initial[k])
    print('Successfully initialized State Tracker.') 

    CurrentState=Initial
    NextState={}

    for i in range(0,NumIter):
        NextState=TransFunc(CurrentState)
        for k in States.keys():
            States[k].append(NextState[k])
        CurrentState=NextState

    return States


NumIter=15
StateTracker=LSS(set_initial({'glucose_ext':True}),NumIter,Transition)
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

WTcomp=pd.read_csv(PATH+'/wt-responses.csv',sep='\t',index_col=0)

WTcompdict=WTcomp.to_dict()

# sim_gluc=ss_extractor(StateTracker)
# counter=0
# for k in WTcompdict['STATE'].keys():
#     if sim_gluc[k]!=WTcompdict['STATE'][k]:
#         print(k)
# print(counter)
