
#This function calculates a cord across the Earth.  LLo and LLhi are the minimum and maximum values
# of the parameter along the cord.  LToH converts values of L to atmospheric heights H
def LToHConv(b,Rs,H):
    LToH=0
    if b>1.:
        LToH=lambda(L):(pow(L*L + b*b*Rs*Rs, 0.5)-Rs)
        LLo=-Rs*pow(pow(1+H/Rs,2)-b*b,0.5)
        LHi=Rs*pow(pow(1+H/Rs,2)-b*b,0.5)
    else :
        LToH=lambda(L):(Rs* (pow(b*b+pow(pow(1-b*b,0.5)+L/Rs,2),0.5)-1) )
        LLo=0
        LHi=pow(pow(H+Rs,2)-b*b*Rs*Rs,0.5)-Rs*pow(1-b*b,0.5)
    return [LToH, LLo, LHi]

