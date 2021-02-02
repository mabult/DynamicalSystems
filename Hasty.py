import numpy as np

#Enzymatic Degradation
#------------------------------------------------------------------------------------------------
def function_f(state,p):
    x = (2*state[1] + 4*state[3] + 6*state[4] + 4*state[5] + 8*state[6] + 10*state[7] + 8*state[8]
    + 2*state[10] + 4*state[12] + 6*state[13] + 4*state[14] + 8*state[15] + 10*state[16] + 8*state[17]
    + 2*state[19] + 4*state[21] + 6*state[22] + 4*state[23] + 8*state[24] + 10*state[25] + 8*state[26]   
    + state[30] + state[31] + state[32] + state[33] + state[34] + state[35] + 2*state[36] + 2*state[37] + 4*state[38])
    
    return p[22]/(p[23] + x)


#------------------------------------------------------------------------------------------------

def Hasty(state,t,p):
    # unpack the state vector
    # 39 dimensions ...

    deriv_state = [0.0] * 39

    for x in state:
        x = max(0,x) #just in case...
        
# first 9 states represent the promoter encoding araC and the promoter encoding GFP   
    # state[0] : P(a,0,0) no a2, no r4
    # state[1] : P(a,1,0) 1 a2, no r4
    # state[2] : P(a,L,0) looped DNA, no r4
    # state[3] : P(a,0,1) no a2, 1 r4
    # state[4] : P(a,1,1) 1 a2, 1 r4
    # state[5] : P(a,L,1) looped DNA, 1 r4
    # state[6] : P(a,0,2) no a2, 2 r4
    # state[7] : P(a,1,2) 1 a2, 2 r4
    # state[8] : P(a,L,2) looped DNA, 2 r4
    
    # state[9] : P(r,0,0) no a2, no r4
    # state[10] : P(r,1,0) 1 a2, no r4
    # state[11] : P(r,L,0) looped DNA, no r4
    # state[12] : P(r,0,1) no a2, 1 r4
    # state[13] : P(r,1,1) 1 a2, 1 r4
    # state[14] : P(r,L,1) looped DNA, 1 r4 
    # state[15] : P(r,0,2) no a2, 2 r4
    # state[16] : P(r,1,2) 1 a2, 2 r4
    # state[17] : P(r,L,2) looped DNA, 2 r4

    # state[18] : P(g,0,0) no a2, no r4
    # state[19] : P(g,1,0) 1 a2, no r4
    # state[20] : P(g,L,0) looped DNA, no r4
    # state[21] : P(g,0,1) no a2, 1 r4
    # state[22] : P(g,1,1) 1 a2, 1 r4
    # state[23] : P(g,L,1) looped DNA, 1 r4 
    # state[24] : P(g,0,2) no a2, 2 r4
    # state[25] : P(g,1,2) 1 a2, 2 r4
    # state[26] : P(g,L,2) looped DNA, 2 r4

    # state[27] : m_a mRNA for activator
    # state[28] : m_r mRNA for repressor
    # state[29] : m_g mRNA for GFP    
    # state[30] : a_uf unfolded activator
    # state[31] : r_uf unfolded repressor   
    # state[32] : g_uf unfolded GFP    
    # state[33] : a folded activator
    # state[34] : r folded repressor
    # state[35] : g folded GFP  
    # state[36] : a2 activator in dimeric form
    # state[37] : r2 repressor in dimeric form
    # state[38] : r4 repressor in tetrameric form
    


    # Parameters p (defined outside of function)
    
    # Activation/Repression of Promoters
    # p[0] : k(a) 
    # p[1] : k(-a)
    # p[2] : k(r)
    # p[3] : k(-r)
    # p[4] : k(l)
    # p[5] : k(ul)
    
    # Protein Synthesis (transcription, translation, folding)
    # p[6] : b(a)
    # p[7] : b(r)
    # p[8] : alpha
    # p[9] : t(a)
    # p[10] : t(r)
    # p[11] : k(fa)
    # p[12] : k(fr)
    
    # Dimerisation / tetramerisation of Proteins
    # p[13] : k(da)
    # p[14] : k(-da)
    # p[15] : k(dr)
    # p[16] : k(-dr) 
    # p[17] : k(t)
    # p[18] : k(-t)
    
    # Degradation parameters
    # p[19] : d(a)
    # p[20] : d(r)
    # p[21] : lambda <-> d(a) Relative degradation rate of the activator / repressor 
    # p[23] : c(l) (function f(X))
    # p[24] : epsilon
    
    # Added Parameters (relating to Reporting System)
    # p[26] : k(fg) Maturation rate of the reporter protein
    # p[27] : lambda2 <-> d(g) Relative degradation rate of the reporter / repressor 
    
    # Now ... the reaction rates
    # We are going to list all the reactions and compute their rates
    # the rates will be stored in the list rates
    
    # binding of a2 to first promoter 
    # rRate : P(a,0,0) + a2 -> P(a,1,0) 
    rRate = p[0]*state[36]*state[0]
    deriv_state[0] += -rRate
    deriv_state[36] += -rRate   
    deriv_state[1] += rRate
    # rRate : P(a,0,0) + a2 <- P(a,1,0)     
    rRate = p[1]*state[1]    
    deriv_state[0] += rRate
    deriv_state[36] += rRate   
    deriv_state[1] += -rRate
    # rRate : P(a,0,1) + a2 -> P(a,1,1) 
    rRate = p[0]*state[36]*state[3]
    deriv_state[3] += -rRate
    deriv_state[36] += -rRate   
    deriv_state[4] += rRate
    # rRate : P(a,0,1) + a2 <- P(a,1,1)     
    rRate = p[1]*state[4]
    deriv_state[3] += rRate
    deriv_state[36] += rRate   
    deriv_state[4] += -rRate   
    # rRate : P(a,0,2) + a2 -> P(a,1,2) 
    rRate = p[0]*state[36]*state[6]
    deriv_state[6] += -rRate
    deriv_state[36] += -rRate   
    deriv_state[7] += rRate    
    # rRate : P(a,0,2) + a2 <- P(a,1,2)     
    rRate = p[1]*state[7]
    deriv_state[6] += rRate
    deriv_state[36] += rRate   
    deriv_state[7] += -rRate       
    
    # binding of a2 to second promoter 
    # rRate : P(r,0,0) + a2 -> P(r,1,0) 
    rRate = p[0]*state[36]*state[9]
    deriv_state[9] += -rRate
    deriv_state[36] += -rRate   
    deriv_state[10] += rRate     
    # rRate : P(r,0,0) + a2 <- P(r,1,0)     
    rRate= p[1]*state[10]   
    deriv_state[9] += rRate
    deriv_state[36] += rRate   
    deriv_state[10] += -rRate  
    # rRate : P(r,0,1) + a2 -> P(r,1,1) 
    rRate = p[0]*state[36]*state[12]
    deriv_state[12] += -rRate
    deriv_state[36] += -rRate   
    deriv_state[13] += rRate       
    # rRate : P(r,0,1) + a2 <- P(r,1,1)     
    rRate = p[1]*state[13]
    deriv_state[12] += rRate
    deriv_state[36] += rRate   
    deriv_state[13] += -rRate   
    # rRate : P(r,0,2) + a2 -> P(r,1,2) 
    rRate = p[0]*state[36]*state[15]
    deriv_state[15] += -rRate
    deriv_state[36] += -rRate 
    deriv_state[16] += rRate    
    # rRate : P(r,0,2) + a2 <- P(r,1,2)     
    rRate = p[1]*state[16]    
    deriv_state[15] += rRate
    deriv_state[36] += rRate
    deriv_state[16] += -rRate    
    
    # binding of a2 to third promoter 
    # rRate : P(g,0,0) + a2 -> P(g,1,0) 
    rRate = p[0]*state[36]*state[18]
    deriv_state[18] += -rRate
    deriv_state[36] += -rRate   
    deriv_state[19] += rRate     
    # rRate : P(g,0,0) + a2 <- P(g,1,0)     
    rRate = p[1]*state[19]   
    deriv_state[18] += rRate
    deriv_state[36] += rRate   
    deriv_state[19] += -rRate  
    # rRate : P(g,0,1) + a2 -> P(g,1,1) 
    rRate = p[0]*state[36]*state[21]
    deriv_state[21] += -rRate
    deriv_state[36] += -rRate   
    deriv_state[22] += rRate       
    # rRate : P(g,0,1) + a2 <- P(g,1,1)     
    rRate = p[1]*state[22]
    deriv_state[21] += rRate
    deriv_state[36] += rRate   
    deriv_state[22] += -rRate   
    # rRate : P(g,0,2) + a2 -> P(g,1,2) 
    rRate = p[0]*state[36]*state[24]
    deriv_state[24] += -rRate
    deriv_state[36] += -rRate 
    deriv_state[25] += rRate    
    # rRate : P(g,0,2) + a2 <- P(g,1,2)     
    rRate = p[1]*state[25]    
    deriv_state[24] += rRate
    deriv_state[36] += rRate
    deriv_state[25] += -rRate       
    
    # First binding of r4 to first promoter 
    # rRate : P(a,0,0) + r4 -> P(a,0,1) 
    rRate = 2*p[2]*state[38]*state[0]
    deriv_state[0] += -rRate
    deriv_state[38] += -rRate
    deriv_state[3] += rRate  
    # rRate : P(a,0,0) + r4 <- P(a,0,1)     
    rRate = p[3]*state[3]
    deriv_state[0] += rRate
    deriv_state[38] += rRate
    deriv_state[3] += -rRate  
    # rRate : P(a,1,0) + r4 -> P(a,1,1) 
    rRate = 2*p[2]*state[38]*state[1]
    deriv_state[1] += -rRate
    deriv_state[38] += -rRate
    deriv_state[4] += rRate 
    # rRate : P(a,1,0) + r4 <- P(a,1,1)     
    rRate = p[3]*state[4]
    deriv_state[1] += rRate
    deriv_state[38] += rRate
    deriv_state[4] += -rRate
    # Second binding of r4 to first promoter 
    # rRate : P(a,0,1) + r4 -> P(a,0,2) 
    rRate = p[2]*state[38]*state[3]
    deriv_state[3] += -rRate
    deriv_state[38] += -rRate
    deriv_state[6] += rRate     
    # rRate : P(a,0,1)+ r4 <- P(a,0,2)     
    rRate = 2*p[3]*state[6]
    deriv_state[3] += rRate
    deriv_state[38] += rRate
    deriv_state[6] += -rRate         
    # rRate : P(a,1,1) + r4 -> P(a,1,2)
    rRate = p[2]*state[38]*state[4]
    deriv_state[4] += -rRate
    deriv_state[38] += -rRate
    deriv_state[7] += rRate      
    # rRate : P(a,1,1) + r4 <- P(a,1,2)   
    rRate = 2*p[3]*state[7]
    deriv_state[4] += rRate
    deriv_state[38] += rRate
    deriv_state[7] += -rRate       

    # First binding of r4 to second promoter 
    # rRate : P(r,0,0) + r4 -> P(r,0,1) 
    rRate = 2*p[2]*state[38]*state[9]
    deriv_state[9] += -rRate
    deriv_state[38] += -rRate
    deriv_state[12] += rRate       
    # rRate : P(r,0,0) + r4 <- P(r,0,1)     
    rRate = p[3]*state[12]
    deriv_state[9] += rRate
    deriv_state[38] += rRate
    deriv_state[12] += -rRate       
    # rRate : P(r,1,0) + r4 -> P(r,1,1) 
    rRate = 2*p[2]*state[38]*state[10]
    deriv_state[10] += -rRate
    deriv_state[38] += -rRate
    deriv_state[13] += rRate      
    # rRate : P(r,1,0) + r4 <- P(r,1,1)     
    rRate = p[3]*state[13]
    deriv_state[10] += rRate
    deriv_state[38] += rRate
    deriv_state[13] += -rRate 
    
    # Second binding of r4 to second promoter 
    # rRate : P(r,0,1) + r4 -> P(r,0,2) 
    rRate = p[2]*state[38]*state[12]
    deriv_state[12] += -rRate
    deriv_state[38] += -rRate
    deriv_state[15] += rRate       
    # rRate : P(r,0,1)+ r4 <- P(r,0,2)     
    rRate = 2*p[3]*state[15]
    deriv_state[12] += rRate
    deriv_state[38] += rRate
    deriv_state[15] += -rRate     
    # rRate : P(r,1,1) + r4 -> P(r,1,2)
    rRate = p[2]*state[38]*state[13]
    deriv_state[13] += -rRate
    deriv_state[38] += -rRate
    deriv_state[16] += rRate       
    # rRate : P(r,1,1) + r4 <- P(r,1,2)   
    rRate = 2*p[3]*state[16]
    deriv_state[13] += rRate
    deriv_state[38] += rRate
    deriv_state[16] += -rRate       
    
    # First binding of r4 to third promoter 
    # rRate : P(g,0,0) + r4 -> P(g,0,1) 
    rRate = 2*p[2]*state[38]*state[18]
    deriv_state[18]+=-rRate
    deriv_state[38] += -rRate
    deriv_state[21] += rRate       
    # rRate : P(g,0,0) + r4 <- P(g,0,1)     
    rRate = p[3]*state[21]
    deriv_state[18] += rRate
    deriv_state[38] += rRate
    deriv_state[21] += -rRate       
    # rRate : P(g,1,0) + r4 -> P(g,1,1) 
    rRate = 2*p[2]*state[38]*state[19]
    deriv_state[19] += -rRate
    deriv_state[38] += -rRate
    deriv_state[22] += rRate      
    # rRate : P(g,1,0) + r4 <- P(g,1,1)     
    rRate = p[3]*state[22]
    deriv_state[19] += rRate
    deriv_state[38] += rRate
    deriv_state[22] += -rRate 
    
    # Second binding of r4 to third promoter 
    # rRate : P(g,0,1) + r4 -> P(g,0,2) 
    rRate = p[2]*state[38]*state[21]
    deriv_state[21] += -rRate
    deriv_state[38] += -rRate
    deriv_state[24] += rRate       
    # rRate : P(g,0,1)+ r4 <- P(g,0,2)     
    rRate = 2*p[3]*state[24]
    deriv_state[21] += rRate
    deriv_state[38] += rRate
    deriv_state[24] += -rRate     
    # rRate : P(g,1,1) + r4 -> P(g,1,2)
    rRate = p[2]*state[38]*state[22]
    deriv_state[22] += -rRate
    deriv_state[38] += -rRate
    deriv_state[25] += rRate       
    # rRate : P(g,1,1) + r4 <- P(g,1,2)   
    rRate = 2*p[3]*state[25]
    deriv_state[22] += rRate
    deriv_state[38] += rRate
    deriv_state[25] += -rRate 
    
    
    # Loops on first promoter
    # rRate : P(a,1,2) -> P(a,L,2) +a2
    rRate = p[4]*state[7]
    deriv_state[7] += -rRate
    deriv_state[8] += rRate
    deriv_state[36] += rRate    
    # rRate : P(a,0,2) -> P(a,L,2)
    rRate = p[4]*state[6]  
    deriv_state[6] += -rRate
    deriv_state[8] += rRate
    # rRate : P(a,L,0) -> P(a,0,0)
    # this is odd shouldn't it be 2 instead of 0?????
    # Answer: No as there is a path from P(a,L,2) to P(a,L,0) through degradation
    rRate = p[5]*state[2]
    deriv_state[2] += -rRate
    deriv_state[0] += rRate
    
    # Loops on second promoter
    # rRate : P(r,1,2) -> P(r,L,2) +a2
    rRate = p[4]*state[16]
    deriv_state[16] += -rRate
    deriv_state[17] += rRate
    deriv_state[36] += rRate  
    # rRate : P(r,0,2) -> P(r,L,2)
    rRate = p[4]*state[15] 
    deriv_state[15] += -rRate
    deriv_state[17] += rRate
    # rRate : P(r,L,0) -> P(r,0,0)
    # this is odd shouldn't it be 2 instead of 0?????
    # Answer: No as there is a path from P(r,L,2) to P(r,L,0) through degradation
    rRate = p[5]*state[11]
    deriv_state[11] += -rRate
    deriv_state[9] += rRate    
    
    # Loops on third promoter
    # rRate : P(g,1,2) -> P(g,L,2) +a2
    rRate = p[4]*state[25]
    deriv_state[25] += -rRate
    deriv_state[26] += rRate
    deriv_state[36] += rRate  
    # rRate : P(g,0,2) -> P(g,L,2)
    rRate = p[4]*state[24] 
    deriv_state[24] += -rRate
    deriv_state[26] += rRate
    # rRate : P(g,L,0) -> P(g,0,0)
    # this is odd shouldn't it be 2 instead of 0?????
    # Answer: No as there is a path from P(g,L,2) to P(g,L,0) through degradation
    rRate = p[5]*state[20]
    deriv_state[20] += -rRate
    deriv_state[18] += rRate 
    
    
    # Ok now the more common stuff... 
    # transcription first promoter
    # rRate : P(a,0,0) -> P(a,0,0) +ma
    rRate = p[6]*state[0]  
    deriv_state[27] += rRate      
    # rRate : P(a,1,0) -> P(a,1,0) +ma
    rRate = p[6]*p[8]*state[1]
    deriv_state[27] += rRate      
    # transcription second promoter
    # rRate : P(r,0,0) -> P(r,0,0) +mr
    rRate = p[6]*state[9] 
    deriv_state[28] += rRate      
    # rRate : P(r,1,0) -> P(r,1,0) +mr
    rRate = p[6]*p[8]*state[10] 
    deriv_state[28] += rRate   
    
    # rRate : P(g,0,0) -> P(g,0,0) +mg
    rRate = p[6]*state[18]    
    deriv_state[29] += rRate 
    # rRate : P(g,1,0) -> P(g,1,0) +mg
    rRate = p[6]*p[8]*state[19]
    deriv_state[29] += rRate  
    
    # translation
    # rRate : ma -> auf +ma
    rRate = p[9]*state[27] 
    deriv_state[30] += rRate       
    # rRate : mr -> ruf +mr
    rRate = p[10]*state[28] 
    deriv_state[31] += rRate  
    # rRate : mg -> guf +mg
    rRate = p[9]*state[29]    
    deriv_state[32] += rRate 
    
    # folding
    # rRate : auf -> a
    rRate = p[11]*state[30] 
    deriv_state[33] += rRate  
    deriv_state[30] += -rRate      
    # rRate : ruf -> r
    rRate = p[12]*state[31] 
    deriv_state[34] += rRate  
    deriv_state[31] += -rRate  
    # rRate : guf -> g 
    rRate = p[26]*state[32]  # MODIFIED : PARAMETER ADDED
    deriv_state[32] += -rRate 
    deriv_state[35] += rRate  
    
    # dimerisation and tetramerisation
    # rRate : a+a -> a2
    rRate = p[13]*state[33]*state[33] 
    #rRate=p[13]*state[33]*(state[33]-1)/2 
    deriv_state[36] += rRate  
    deriv_state[33] += -2*rRate      
    # rRate : a+a <- a2
    rRate = p[14]*state[36] 
    deriv_state[36] += -rRate  
    deriv_state[33] += 2*rRate     
    # rRate : r+r -> r2
    rRate = p[15]*state[34]*state[34] 
    #rRate=p[15]*state[34]*(state[34]-1)/2 
    deriv_state[37] += rRate  
    deriv_state[34] += -2*rRate     
    # rRate : r+r <- r2
    rRate = p[16]*state[37]
    deriv_state[37] += -rRate  
    deriv_state[34] += 2*rRate    
     # rRate : r2+r2 -> r4
    rRate = p[17]*state[37]*state[37]
    #rRate=p[17]*state[37]*(state[37]-1)/2 
    deriv_state[38] += rRate  
    deriv_state[37] += -2*rRate  
    # rRate : r2+r2 <- r4
    rRate = p[18]*state[38] 
    deriv_state[38] += -rRate  
    deriv_state[37] += 2*rRate   
    
    
    # Degradation
    # rRate : ma -> 
    rRate = p[19]*state[27]
    deriv_state[27] += -rRate 
    # rRate : mr -> 
    rRate = p[20]*state[28]
    deriv_state[28] += -rRate 
    # rRate : mg -> 
    rRate = p[25]*p[19]*state[29]   
    deriv_state[29] += -rRate 
    
    # Introduce f function (Enzymatic Degradation)
    f_f = function_f(state,p)
    
    # rRate : auf -> 
    rRate = p[21]*f_f*state[30]
    deriv_state[30] += -rRate 
    # rRate : ruf -> 
    rRate = f_f*state[31]
    deriv_state[31] += -rRate 
    # rRate : guf -> 
    rRate= p[27]*f_f*state[32] # MODIFIED : PARAMETER ADDED
    deriv_state[32] += -rRate 
    
    # rRate : a -> 
    rRate = p[21]*f_f*state[33]
    deriv_state[33] += -rRate 
    # rRate : r -> 
    rRate = f_f*state[34]
    deriv_state[34] += -rRate 
    # rRate : g-> 
    rRate = p[27]*f_f*state[35] # MODIFIED : PARAMETER ADDED
    deriv_state[35] += -rRate   
    
    # rRate : a2 -> 
    rRate = p[21]*f_f*state[36]
    deriv_state[36] += -rRate 
    # rRate : r2 -> 
    rRate = f_f*state[37] 
    deriv_state[37] += -rRate 
    # rRate : r4 -> 
    rRate = f_f*state[38]  
    deriv_state[38] += -rRate 

    # Degradation of protein complexes bound to the promoter
    # Degradation of a2
    # rRate : P(a,1,0)->P(a,0,0)
    rRate = f_f*state[1]
    deriv_state[1] += -rRate  
    deriv_state[0] += rRate 
    # rRate : P(a,1,1)->P(a,0,1)
    rRate = f_f*state[4]    
    deriv_state[4] += -rRate  
    deriv_state[3] += rRate 
    # rRate : P(a,1,2)->P(a,0,2)
    rRate = f_f*state[7]   
    deriv_state[7] += -rRate  
    deriv_state[6] += rRate 
    
    # rRate : P(r,1,0)->P(r,0,0)
    rRate = f_f*state[10]  
    deriv_state[10] += -rRate  
    deriv_state[9] += rRate 
    # rRate : P(r,1,1)->P(r,0,1)
    rRate = f_f*state[13]  
    deriv_state[13] += -rRate  
    deriv_state[12] += rRate     
    # rRate : P(r,1,2)->P(r,0,2)
    rRate = f_f*state[16]       
    deriv_state[16] += -rRate  
    deriv_state[15] += rRate 

    # rRate : P(g,1,0)->P(g,0,0)
    rRate = f_f*state[19]  
    deriv_state[19] += -rRate  
    deriv_state[18] += rRate 
    # rRate : P(g,1,1)->P(g,0,1)
    rRate = f_f*state[22]  
    deriv_state[22] += -rRate  
    deriv_state[21] += rRate     
    # rRate : P(g,1,2)->P(g,0,2)
    rRate = f_f*state[25]       
    deriv_state[25] += -rRate  
    deriv_state[24] += rRate
    
    # Degradation of r2
    # rRate : P(a,0,1)->P(a,0,0)
    rRate = f_f*state[3]
    deriv_state[3] += -rRate  
    deriv_state[0] += rRate 
    # rRate : P(a,1,1)->P(a,1,0)
    rRate = f_f*state[4]
    deriv_state[4] += -rRate  
    deriv_state[1] += rRate     
    # rRate : P(a,0,2)->P(a,0,1)
    rRate = 2*f_f*state[6]   
    deriv_state[6] += -rRate  
    deriv_state[3] += rRate 
    # rRate : P(a,1,2)->P(a,1,1)
    rRate = 2*f_f*state[7] 
    deriv_state[7] += -rRate  
    deriv_state[4] += rRate     
    
    # rRate : P(r,0,1)->P(r,0,0)
    rRate = f_f*state[12] 
    deriv_state[12] += -rRate  
    deriv_state[9] += rRate 
    # rRate : P(r,1,1)->P(r,1,0)
    rRate = f_f*state[13]
    deriv_state[13] += -rRate  
    deriv_state[10] += rRate 
    # rRate : P(r,0,2)->P(r,0,1)
    rRate = 2*f_f*state[15]  
    deriv_state[15] += -rRate  
    deriv_state[12] += rRate 
    # rRate : P(r,1,2)->P(r,1,1)
    rRate = 2*f_f*state[16]      
    deriv_state[16] += -rRate  
    deriv_state[13] += rRate 

    # rRate : P(g,0,1)->P(g,0,0)
    rRate = f_f*state[21]
    deriv_state[21] += -rRate  
    deriv_state[18] += rRate 
    # rRate : P(g,1,1)->P(g,1,0)
    rRate = f_f*state[22]
    deriv_state[22] += -rRate  
    deriv_state[19] += rRate 
    # rRate : P(g,0,2)->P(g,0,1)
    rRate = 2*f_f*state[24]   
    deriv_state[24] += -rRate  
    deriv_state[21] += rRate 
    # rRate : P(g,1,2)->P(g,1,1)
    rRate = 2*f_f*state[25] 
    deriv_state[25] += -rRate  
    deriv_state[22] += rRate 
    
    
    # Unlooping the DNA by degrading r4
    # rRate : P(a,L,2)->P(a,L,1)
    rRate = 2*p[24]*f_f*state[8]    
    deriv_state[8] += -rRate  
    deriv_state[5] += rRate 
    # rRate : P(a,L,1)->P(a,L,0)
    rRate = p[24]*f_f*state[5]  
    deriv_state[5] += -rRate  
    deriv_state[2] += rRate 
    # rRate : P(r,L,2)->P(r,L,1)
    rRate = 2*p[24]*f_f*state[17]  
    deriv_state[17] += -rRate  
    deriv_state[14] += rRate 
    # rRate : P(r,L,1)->P(r,L,0)
    rRate = p[24]*f_f*state[14]   
    deriv_state[14] += -rRate  
    deriv_state[11] += rRate   
    # rRate : P(g,L,2)->P(g,L,1)
    rRate = 2*p[24]*f_f*state[26]  
    deriv_state[26] += -rRate  
    deriv_state[23] += rRate 
    # rRate : P(g,L,1)->P(g,L,0)
    rRate = p[24]*f_f*state[23]   
    deriv_state[23] += -rRate  
    deriv_state[20] += rRate     
    
   
    return deriv_state

#------------------------------------------------------------------------------------------------

# kr and ka parameters
def k_r(iptg,A,Cr_min,Cr_max,b_1,kr_1):
    x = A*(Cr_min+(Cr_max-Cr_min)/(1+(iptg/kr_1)**b_1))
    return x

def k_a(ara,iptg,A,Ca_min,Ca_max,c_1,ka_1,b_2,kr_2):
    x = A*(Ca_min+(Ca_max-Ca_min)*((ara**c_1)/(ara**c_1+ka_1**c_1))/(1+(iptg/kr_2)**b_2))
    return x


#------------------------------------------------------------------------------------------------
def initial_parameters():
    p = [0.0] * 28# list with 28 elements (2 new paramters realted to GFP)

    # Activation/Repression of Promoters
    # p[0] : k(a) 
    # See expression in next cell
    # p[1] : k(-a)
    p[1] = 1.8
    # p[2] : k(r)
    # See expression in next cell  
    # p[3] : k(-r)
    p[3] = 1.8    
    # p[4] : k(l)
    p[4] = 0.36
    # p[5] : k(ul)
    p[5] = 0.18

    # Protein Synthesis (transcription, translation, folding)
    # p[6] : b(a)
    p[6] = 0.36
    # p[7] : b(r)
    p[7] = 0.36
    # p[8] : alpha
    p[8] = 20    
    # p[9] : t(a)
    p[9] = 90
    # p[10] : t(r)
    p[10] = 90
    # p[11] : k(fa)
    p[11] = 0.90
    # p[12] : k(fr)
    p[12] = 0.90    
 
    # Dimerisation / tetramerisation of Proteins
    # p[13] : k(da)
    p[13] = 0.018
    # p[14] : k(-da)
    p[14] = 0.00018
    # p[15] : k(dr)
    p[15] = 0.018
    # p[16] : k(-dr) 
    p[16] = 0.00018
    # p[17] : k(t)
    p[17] = 0.018
    # p[18] : k(-t)
    p[18] = 0.00018    

    # Degradation parameters
    # p[19] : d(a)
    p[19] = 0.54 

    # p[20] : d(r)
    p[20] = 0.54
    #p[20] = 0.01
    # p[21] : lambda
    p[21] = 2.5
    # p[22] : gamma (function f(X))
    p[22] = 1080
    # p[23] : c(l) (function f(X))
    p[23] = 0.1
    # p[24] : epsilon
    p[24] = 0.2
    # p[25] addition : ratio of degradation rate of mg/ma
    p[25] = 1
    #p[25]=2.8
    
    # extra parameters 
    # p[26] : k(fg) # folding parameter for G
    p[26] = 0.90

    # p[27] : lambda2 (relative degradation rate / R)
    p[27] = 2.5
      
    return p

# Parameters are identical to vector above, except a few components listed in modifications
# and components 0 and 2 which are computed from the rest
def modify_parameters(parameters,modifications,inducers):
    for m in modifications:
        parameters[int(m[0])] = m[1]
    parameters[0] = k_a(inducers['ara']['value'],inducers['iptg']['value'],parameters[1],0,1,2,2.5,2,1.8)
    parameters[2] = k_r(inducers['iptg']['value'],parameters[3],0.01,0.2,2,0.035)

    
# Setting Up Initial Conditions
# Initial conditions are all zero, except a few components listed in modifications
def initial_conditions(modifications):
    initial = np.zeros(39)
    for m in modifications:
        initial[int(m[0])] = m[1]
    return initial
