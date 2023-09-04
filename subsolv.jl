########################################################################################################
### GCMMA-MMA-Julia                                                                                  ### 
###                                                                                                  ###
### This file is part of GCMMA-MMA-Julia. GCMMA-MMA-Julia is licensed under the terms of GNU         ###
### General Public License as published by the Free Software Foundation. For more information and    ###
### the LICENSE file, see <https://github.com/pollinico/TopOpt_jl/blob/main/LICENSE>.                ###
###                                                                                                  ###
### The orginal work is written by Krister Svanberg in MATLAB.                                       ###
### This is the Julia version of the code written by Nicolò Pollini.                                 ###
### version 18-05-2023                                                                               ###
########################################################################################################

#-------------------------------------------------------------
#
#    Copyright (C) 2006 Krister Svanberg
#
#    This file, subsolv.m, is part of GCMMA-MMA-code.
#    
#    GCMMA-MMA-code is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License as 
#    published by the Free Software Foundation; either version 3 of 
#    the License, or (at your option) any later version.
#    
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    (file COPYING) along with this file.  If not, see 
#    <http://www.gnu.org/licenses/>.
#    
#    You should have received a file README along with this file,
#    containing contact information.  If not, see
#    <http://www.smoptit.se/> or e-mail mmainfo@smoptit.se or krille@math.kth.se.
#
#    Version Dec 2006.
#
#
function subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d)
    #
    # This function subsolv solves the MMA subproblem:
    #         
    # minimize   SUM[ p0j/(uppj-xj) + q0j/(xj-lowj) ] + a0*z +
    #          + SUM[ ci*yi + 0.5*di*(yi)^2 ],
    #
    # subject to SUM[ pij/(uppj-xj) + qij/(xj-lowj) ] - ai*z - yi <= bi,
    #            alfaj <=  xj <=  betaj,  yi >= 0,  z >= 0.
    #        
    # Input:  m, n, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d.
    # Output: xmma,ymma,zmma, slack variables and Lagrange multiplers.
    #
    een = ones(n)
    eem = ones(m)
    epsi = 1
    epsvecn = epsi*een
    epsvecm = epsi*eem
    x = 0.5*(alfa+beta)
    y = copy(eem)
    z = 1
    lam = copy(eem)
    xsi = een./(x-alfa)
    xsi = max.(xsi,een)
    eta = een./(beta-x)
    eta = max.(eta,een)
    mu  = max.(eem,0.5*c)
    zet = 1
    s = copy(eem)
    itera = 0
    while epsi > epsimin
        epsvecn = epsi*een
        epsvecm = epsi*eem
        ux1 = upp-x
        xl1 = x-low
        ux2 = ux1.*ux1
        xl2 = xl1.*xl1
        uxinv1 = een./ux1
        xlinv1 = een./xl1
        plam = p0 + P'*lam 
        qlam = q0 + Q'*lam 
        gvec = P*uxinv1 + Q*xlinv1
        dpsidx = plam./ux2 - qlam./xl2 
        rex = dpsidx - xsi + eta
        rey = c + d.*y - mu - lam
        rez = a0 - zet - a'*lam
        relam = gvec - a*z - y + s - b
        rexsi = xsi.*(x-alfa) - epsvecn
        reeta = eta.*(beta-x) - epsvecn
        remu = mu.*y - epsvecm
        rezet = zet*z - epsi
        res = lam.*s - epsvecm
        residu1 = [rex' rey' rez]'
        residu2 = [relam' rexsi' reeta' remu' rezet res']'
        residu = [residu1' residu2']'
        residunorm = sqrt(residu'*residu)
        residumax = maximum(abs.(residu))
        ittt = 0;
        while (residumax > 0.9*epsi) & (ittt < 200)
            ittt = ittt + 1
            itera = itera + 1
            ux1 = upp-x
            xl1 = x-low
            ux2 = ux1.*ux1
            xl2 = xl1.*xl1
            ux3 = ux1.*ux2
            xl3 = xl1.*xl2
            uxinv1 = een./ux1
            xlinv1 = een./xl1
            uxinv2 = een./ux2
            xlinv2 = een./xl2
            plam = p0 + P'*lam
            qlam = q0 + Q'*lam
            gvec = P*uxinv1 + Q*xlinv1
            GG = P*spdiagm(n,n,uxinv2) - Q*spdiagm(n,n,xlinv2)
            dpsidx = plam./ux2 - qlam./xl2 
            delx = dpsidx - epsvecn./(x-alfa) + epsvecn./(beta-x)
            dely = c + d.*y - lam - epsvecm./y
            delz = a0 - a'*lam - epsi/z
            dellam = gvec - a*z - y - b + epsvecm./lam
            diagx = plam./ux3 + qlam./xl3
            diagx = 2*diagx + xsi./(x-alfa) + eta./(beta-x)
            diagxinv = een./diagx
            diagy = d + mu./y
            diagyinv = eem./diagy
            diaglam = s./lam
            diaglamyi = diaglam+diagyinv
            if m < n
                blam = dellam + dely./diagy - GG*(delx./diagx)
                bb = [blam' delz]'
                Alam = spdiagm(m,m,diaglamyi) + GG*spdiagm(n,n,diagxinv)*GG'
                AA = [Alam     a;
                        a'    -zet/z ]
                solut = AA\bb
                dlam = solut[1:m]
                dz = solut[m+1]
                dx = -delx./diagx - (GG'*dlam)./diagx
            else
                diaglamyiinv = eem./diaglamyi
                dellamyi = dellam + dely./diagy
                Axx = spdiagm(n,n,diagx) + GG'*spdiagm(m,m,diaglamyiinv)*GG
                azz = zet/z + a'*(a./diaglamyi)
                axz = -GG'*(a./diaglamyi)
                bx = delx + GG'*(dellamyi./diaglamyi)
                bz  = delz - a'*(dellamyi./diaglamyi)
                AA = [Axx   axz;
                        axz'  azz ]
                bb = [-bx' -bz]'
                solut = AA\bb
                dx  = solut[1:n]
                dz = solut[n+1]
                dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi
            end
        #
            dy = -dely./diagy + dlam./diagy
            dxsi = -xsi + epsvecn./(x-alfa) - (xsi.*dx)./(x-alfa)
            deta = -eta + epsvecn./(beta-x) + (eta.*dx)./(beta-x)
            dmu  = -mu + epsvecm./y - (mu.*dy)./y
            dzet = -zet + epsi/z - zet*dz/z
            ds   = -s + epsvecm./lam - (s.*dlam)./lam
            xx  = [ y'  z  lam'  xsi'  eta'  mu'  zet  s']'
            dxx = [dy' dz dlam' dxsi' deta' dmu' dzet ds']'
        #    
            stepxx = -1.01*dxx./xx
            stmxx  = maximum(stepxx)
            stepalfa = -1.01*dx./(x-alfa)
            stmalfa = maximum(stepalfa)
            stepbeta = 1.01*dx./(beta-x)
            stmbeta = maximum(stepbeta)
            stmalbe  = max(stmalfa,stmbeta)
            stmalbexx = max(stmalbe,stmxx)
            stminv = max(stmalbexx,1)
            steg = 1/stminv
        #
            xold   =  copy(x)
            yold   =  copy(y)
            zold   =  copy(z)
            lamold =  copy(lam)
            xsiold =  copy(xsi)
            etaold =  copy(eta)
            muold  =  copy(mu)
            zetold =  copy(zet)
            sold   =  copy(s)
        #
            itto = 0
            resinew = 2*residunorm
            while (resinew > residunorm) & (itto < 50)
                itto = itto+1
                x   =   xold + steg*dx
                y   =   yold + steg*dy
                z   =   zold + steg*dz
                lam = lamold + steg*dlam
                xsi = xsiold + steg*dxsi
                eta = etaold + steg*deta
                mu  = muold  + steg*dmu
                zet = zetold + steg*dzet
                s   = sold + steg*ds
                ux1 = upp-x
                xl1 = x-low
                ux2 = ux1.*ux1
                xl2 = xl1.*xl1
                uxinv1 = een./ux1
                xlinv1 = een./xl1
                plam = p0 + P'*lam
                qlam = q0 + Q'*lam
                gvec = P*uxinv1 + Q*xlinv1
                dpsidx = plam./ux2 - qlam./xl2
                rex = dpsidx - xsi + eta
                rey = c + d.*y - mu - lam
                rez = a0 - zet - a'*lam
                relam = gvec - a*z - y + s - b
                rexsi = xsi.*(x-alfa) - epsvecn
                reeta = eta.*(beta-x) - epsvecn
                remu = mu.*y - epsvecm
                rezet = zet*z - epsi
                res = lam.*s - epsvecm
                residu1 = [rex' rey' rez]'
                residu2 = [relam' rexsi' reeta' remu' rezet res']'
                residu = [residu1' residu2']'
                resinew = sqrt(residu'*residu)
                steg = steg/2
            end
        residunorm = copy(resinew)
        residumax = maximum(abs.(residu))
        steg = 2 * steg
        end
        if ittt > 198
            println("epsi: ", epsi)
            println("ittt: ", ittt)
        end
        epsi = 0.1*epsi
    end
    xmma   =  copy(x)
    ymma   =  copy(y)
    zmma   =  copy(z)
    lamma  =  copy(lam)
    xsimma =  copy(xsi)
    etamma =  copy(eta)
    mumma  =  copy(mu)
    zetmma =  copy(zet)
    smma   =  copy(s)
    #-------------------------------------------------------------
    return xmma,ymma,zmma,lamma,xsimma,etamma,mumma,zetmma,smma
end
