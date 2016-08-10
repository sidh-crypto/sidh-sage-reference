# Copyright (C) 2016 InfoSec Global Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

########################################################################

#PublicKey Generation
def key_gen(m,n,P,Q,l,e,Po,Qo,E):
    R = m*P+n*Q
    Ec = E
    Pi = Po
    Qi = Qo
    for i in range(0,e):
        K = (l^(e-1-i))*R
        phi = EllipticCurveIsogeny(Ec, K)
        Ec = phi.codomain()
        R = phi(R)
        Pi = phi(Pi)
        Qi = phi(Qi)
    return Ec, Pi, Qi

#PK for initiator
def key_gen_init(mA,nA):
    return key_gen(mA,nA,PA,QA,lA,eA,PB,QB,E)

#PK for responder
def key_gen_resp(mB,nB):
    return key_gen(mB,nB,PB,QB,lB,eB,PA,QA,E)

########################################################################

#Key Derivation from Own Private and Other Public information
#Input: own private key, own order, other EC, other aux points
def key_der(m,n,l,e,E,P,Q):
    R = m*P+n*Q
    Ec = E
    for i in range(0,e):
        K = (l^(e-1-i))*R
        phi = EllipticCurveIsogeny(Ec, K)
        Ec = phi.codomain()
        R = phi(R)
    key = Ec.j_invariant()
    return key

#Key Derivation by initiator
def key_der_init(m,n,E,P,Q):
    return key_der(m,n,lA,eA,E,P,Q)

#Key Derivation by responder
def key_der_resp(m,n,E,P,Q):
    return key_der(m,n,lB,eB,E,P,Q)

########################################################################

#Testing and test vectors
def test_KA(mA,nA,mB,nB,E):
    In = key_gen_init(mA,nA)
    Re = key_gen_resp(mB,nB)
    kA = key_der_init(mA,nA,Re[0],Re[1],Re[2])
    kB = key_der_resp(mB,nB,In[0],In[1],In[2])
    if kA == kB:
        return kA
    else:
        return "ERROR in KA"



#Parameters

#Curve
lA = 2
lB = 3
eA = 63
eB = 41
f = 11
p = (lA^eA)*(lB^eB)*f - 1
F.<z> = GF(p^2, modulus=x^2+1)
E = EllipticCurve([F(1), F(0)])

########################################################################

#Points
PAx = F(2374093068336250774107936421407893885897*z + 2524646701852396349308425328218203569693)
PAy = F(1944869260414574206229153243510104781725*z + 1309099413211767078055232768460483417201)
PA = E(PAx,PAy)
QA = E(F(-PAx),F(z*PAy))

PBx = F(1556716033657530876728525059284431761206*z + 1747407329595165241335131647929866065215)
PBy = F(3456956202852028835529419995475915388483*z + 1975912874247458572654720717155755005566)
PB = E(PBx,PBy)
QB = E(F(-PBx),F(z*PBy))

#Example
mA = 2575042839726612324
nA = 8801426132580632841

mB = 4558164392438856871
nB = 20473135767366569910

test_KA(mA,nA,mB,nB,E)
