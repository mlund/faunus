c  ********************************************************************
c  ********************************************************************
c  **
c  ** MONTE CARLO SIMULATION OF MULTICOMPONENT ISOTROPIC             **
c  ** IONIC SOLUTION                                                 **
c  **                                                                **   
c  ********************************************************************
c  ********************************************************************


      PROGRAM BULK   

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INCLUDE 'bulk.inc'  
      DIMENSION uula(2),uuta(2),suu(25),tijd(10)
     * ,mtot(mxspec),macc(mxspec),menrj(mxspec),mhcrj(mxspec) 
      DATA uuta / 2*0. / 
      DATA suu / 25*0. /        
      DATA mtot,macc,menrj,mhcrj / 40*0  /         
      DATA kn1,kn2,key,kntot / 4*0 /
      DATA tijd / 10*0. /                  
   43 FORMAT('MONTE CARLO SIMULATION OF MULTICOMPONENT ISOTROPIC    '
     *      ,'IONIC SYSTEM  ',/,/,'Written by P. Bolhuis, B. Jonsson, '
     *      ,'T. Akesson, February 1992',/)
   44 FORMAT('************************************************************
     ************************',/)

      iclk = 0 
      lclk = 0
 
      kkk=12
      lll=13 
      mmm=15
      jjj=6

      WRITE(jjj,44)
      WRITE(jjj,43)
      WRITE(jjj,44)
c         
c        
c  ********************************************************************
c  ** physical constants and conversion factors  
c  ********************************************************************

      pi = acos(-1.)                            
      epsx=8.85418782d-12                      
      ech = 1.60219d-19                       
      avno = 6.0223d23                       
      bk=8.314                              

c  ********************************************************************
c  ** sample units are angstrom and elementary charge  
c  ** kilojoules are used for a mole   
c  **                                 
c  ** MEANING OF THE INPUT VARIABLES                              **
c  **
c  ** ink =  0   from file   
c  ** ink =  1   from slump                         
c  ** nspec      number of ionic components
c  ** hion       array contains the input info on the ions
c  **            (1=nr,2=number,3=radius,4=charge,5=displacement)
c  ** ny(1,2,3)  number of MC loops 
c  ** dtemp      temperature
C  ** eps        relative dielectric constant
C  ** box        box size                                    
C  ** nwint      a widom/collision calc is performed every nwint cycle
C  ** nwins      number of widom particles inserted each calc.
c  ** nfix       where is the chem. pot. measured
c  **            0 = averaged over the cell
c  ** 
c  ********************************************************************
c  ********************************************************************
c                                                          
      OPEN(lll,file='bulk.conf',form='formatted')      
      OPEN(mmm,file='bulk.inp',form='formatted')      
      OPEN(kkk,file='bulk.macro',form='formatted')
c                 
      READ(mmm,*)
      READ(mmm,*)nspec   
      READ(mmm,*)       
      READ(mmm,*)      
      READ(mmm,*)(hion(l,1),hion(l,2),hion(l,3),hion(l,4), hion(l,5)
     *            ,l=1,nspec)
      READ(mmm,*)  
      READ(mmm,*) 
      READ(mmm,*)ny1,ny2,ny3 
      READ(mmm,*)
      READ(mmm,*) 
      READ(mmm,*)ink,islu
      islu = -islu
      READ(mmm,*) 
      READ(mmm,*)
      READ(mmm,*)box
      READ(mmm,*)  
      READ(mmm,*) 
      READ(mmm,*)dtemp,eps
      READ(mmm,*)   
      READ(mmm,*)  
      READ(mmm,*)nwins,nwint,nfix
      IF (nfix.GT.2)  nfix=2           
      IF (islu.EQ.0)  islu=7             
      IF (nwint.EQ.0) nwint=2*ny1         
c                                   


c  ********************************************************************
c  ** BOXSIZE VARIABLES                                              **
c  ********************************************************************
      box2 =box/2.
      box2i=1./box2

c  ********************************************************************
c  **    SET UP HARD CORE AND CHARGE VECTORS                         **
c  **                                                                **
c  **    hc2v contains in the first column the particle type given   **
c  **    by hion(i,1), whereas in the following columns the squared  **
c  **    hard core distances between the particles are given.        **
c  **    chv gives for every particle its charge                     **
c  **    npart is the total number of particles                      **
c  ********************************************************************


      npart =0
      DO 100 i= 1,nspec
        num   = INT(hion(i,2))
        DO 101 j= 1,num
          npart =npart + 1 
          ispc(npart) = i
          DO 102 k= 1,nspec
             hc2v(npart,k) = (hion(i,3) + hion(k,3)) ** 2
102       CONTINUE
          chv(npart) = hion(i,4)
          dp(npart)  = hion(i,5)
101     CONTINUE
100   CONTINUE  
      DO 103 i= 1,nspec
        caver(i)=(hion(i,2))/(box**3)
103   CONTINUE
c      DO 82 i=1,npart
c   82 ei(i)=0.     
      nkonf=ny1*ny2*ny3
      IF(ny3.GT.25) ny3=25
      t = dtemp  
      abeta = - 1. / (bk*dtemp)

c  ********************************************************************
c  ** READ CONFIGURATION FROM FILE                                   **
c  ********************************************************************
     
      IF(ink.EQ.0) THEN
        READ(lll,*) (x6(l),y6(l),z6(l),l=1,npart)
        WRITE(jjj,*) 'Configuration read from file' 
      ENDIF

      IF (ink.EQ.1) CALL slump 
      CALL collision
      CALL widom

c  ********************************************************************
c  ** WRITE INPUT DATA                                               **
c  ********************************************************************

      nwtot = nwins*(INT(ny1/nwint)*ny2*ny3)
      WRITE(jjj,800) npart
      WRITE(jjj,801)(l,INT(hion(l,2)),hion(l,3),hion(l,4),
     *               hion(l,5),caver(l)*1.d27/avno,l=1,nspec)
      WRITE(jjj,802) eps,dtemp,box,ny1,ny2,ny3,nkonf,
     *               DBLE(nkonf/npart),nwtot
  800 FORMAT (/,'VARIABLES FROM INPUT FILE',/,
     */,' number of particles        =',i10,/,
     */,' species   number   radius    charge   dp   average conc',/)
  801 FORMAT (i6,i10,f10.1,f10.1,f10.2,f10.6)
  802 FORMAT (/,' dielectric constant         =',f10.1,/
     * ,' temperature                 =',f10.1,/
     * ,' box size (AA)               =',f10.1,/
     * ,' number of configurations    =',i6,' *',i6,' *',i6,' =',i8,/
     * ,' number of conf. per part    =',f10.1,/
     * ,' widom test per species      =',i10,/)

      ecf=ech**2*1.d10*avno/(eps*epsx*4.*pi)



c  ********************************************************************
c  ** ENERGY EVALUTATION INITIAL CONFIGURATION                       **
c  ********************************************************************

      CALL liv
      ytot1  = xww1*ecf*1.d-3
      WRITE(jjj,'(a)')'Initial configuration'
      WRITE(jjj,'(a,e12.5)')'Coulomb energy (kJ/mol)  =',ytot1


c  ********************************************************************
c  ** START OF MONTE CARLO LOOP                                      **
c  **                                                                **
c  ** Loop over ny1*ny2*ny3 particle moves.                          **
c  ** il    = current particle                                       **
c  ** ei    = new coulomb energy of particle il                      **
c  ** esa   = all coulomb interaction energies of the particle       **
c  ** xww1  = total energy (also ulaa(1))                            **
c  ** dp(1..10)  displacement parameter for every ion            
c  **                                                                **
c  ********************************************************************

      DO 68 my3=1,ny3
        uula(1) = 0.
        DO 67 my2=1,ny2
          DO 66 my1=1,ny1 
            il=1+MOD(kntot,npart) 
            ispec = ispc(il) 
            kntot = kntot + 1
            mtot(ispec) = mtot(ispec) +1
            tx6 = x6(il) + dp(il) * (ran2(islu)-0.5)
            ty6 = y6(il) + dp(il) * (ran2(islu)-0.5)
            tz6 = z6(il) + dp(il) * (ran2(islu)-0.5) 
            IF (tx6.GT.box2) tx6=tx6 - box
            IF (tx6.LT.-box2) tx6=tx6 + box
            IF (ty6.GT.box2) ty6=ty6 - box
            IF (ty6.LT.-box2) ty6=ty6 + box
            IF (tz6.GT.box2) tz6=tz6 - box
            IF (tz6.LT.-box2) tz6=tz6 + box
            isos=99 
            CALL qin 
            IF(isos.NE.0) GO TO 63
c            deltae=0.            
c            DO 2 j=1,npart
c              deltae=deltae+ei(j)-esa(j,il)
c2           CONTINUE
c            dzxp=abeta*(deltae*ecf)  
            dzxp=abeta*(du*ecf)  
            IF(dzxp.LT.-80.) GO TO 64
            IF(dzxp.GT.0.)   GO TO 62
            IF(EXP(dzxp).LT.ran2(islu)) GO TO 64
62          macc(ispec)=macc(ispec)+1
c   ****    trial config accepted
            x6(il)=tx6 
            y6(il)=ty6
            z6(il)=tz6
c            DO 3 j=1,npart
c              esa(il,j)=ei(j)
c              esa(j,il)=ei(j) 
c3           CONTINUE  
c            xww1   =xww1+deltae 
            xww1   = xww1+du
            GO TO 65

   63       mhcrj(ispec)= mhcrj(ispec) + 1 
            GO TO 65
   64       menrj(ispec)= menrj(ispec) + 1
   65       uula(1) = uula(1)+ xww1
            IF(MOD(my1,nwint).EQ.0) THEN  
              CALL collision1            
              CALL widom1               
            ENDIF                 
     
   66     CONTINUE 
   67   CONTINUE  

        qww1=xww1 
        qfww1 = 0.d0
        CALL liv
        IF(xww1.NE.0) qww1=(qww1-xww1)/xww1  
        WRITE(kkk,*)
        WRITE(kkk,44)
	WRITE(kkk,73) my3,qww1,qfww1
        IF(abs(qww1).GT.0.001)STOP 

        CALL liv

        IF(nwint.LE.ny1) THEN 
           CALL collision2   
           CALL widom2      
        ENDIF                                                                     
        uula(1)=uula(1)/DBLE(ny1*ny2)
        uuta(1)=uuta(1)+uula(1)
        uula(1)=uula(1)*ecf*1.d-3
        suu(my3)=uula(1)
c  ***********  calculate total acceptance and rejection ********
        key = 0
        kn1 = 0
        kn2 = 0
        DO 61 i=1,nspec
          key = key + macc(i)
          kn1 = kn1 + mhcrj(i)
          kn2 = kn2 + menrj(i)
61      CONTINUE
        WRITE(kkk,810)kntot,key,kn2,kn1
        WRITE(kkk,'(a,e12.5)')'Coulomb energy (kj/mol)  =',uula(1)
   68 CONTINUE

   73 FORMAT(/' MACROSTEP  ', I10,/,/,'Checkparameters : ',3E10.3)
  810 FORMAT(/,' total number of configurations    =',i8,/, 
     * ' accepted configurations           =',i8,/,  
     * ' energy rejected configurations    =',i8,/, 
     * ' hard core rejected configurations =',i8,/)
  811 FORMAT(/,'Divided per species',/,'species      accepted   ',
     * ' energy rej.   hard core rej.',/,10(i6,3(6X,f8.4),/),/)


c
c  ********************************************************************
c  ** END OF MONTE CARLO LOOP                                        **
c  **                                                                **
c  ** Calculation of the average energies and the pressures          **
c  **                                                                **
c  ********************************************************************

      yn3    = 1./FLOAT(ny3)
      uuta(1)= yn3*uuta(1)
      uuta(1)= uuta(1)*ecf*1.d-3
      CALL earth(suu,uuvar,uuta(1),ny3) 
      WRITE(jjj,'(/)')
      WRITE(jjj,44) 
      WRITE(jjj,'(/,a,I3,a,/)') 'FINAL RESULTS AFTER ',ny3,
     *          ' MACROSTEPS'
      WRITE(jjj,810)kntot,key,kn2,kn1
      WRITE(jjj,811) (l,DBLE(macc(l))/mtot(l),DBLE(menrj(l))
     *         /mtot(l),DBLE(mhcrj(l))/mtot(l),l=1,nspec)
      WRITE(jjj,'(a,2e12.5)')'Coulomb energy (kj/mol)  =',uuta(1),uuvar
      WRITE(jjj,*)
      IF(nwint.LE.ny1) THEN 
         CALL collision3   
         CALL widom3      
      ENDIF              
      REWIND lll   
      WRITE(lll,771)(x6(l),y6(l),z6(l),l=1,npart) 
  771 format(5e16.8) 
c 
c  ********************************************************************
c  **                                                                **
c  ** CALCULATION OF PRESSURES FROM THE DENSITIES                    **
c  **                                                                **
c  ********************************************************************

      pid = 0. 
      DO 34 k = 1,nspec  
 34     pid    = pid + caver(k)*1.d27/avno
      pexen  = uuta(1)*1.d30/(3*bk*dtemp*box**3*avno)
      pexenv = uuvar  *1.d30/(3*bk*dtemp*box**3*avno)
      ptot   = pid + pcollav + pexen 
      ptotv = SQRT(pcollv*pcollv + pexenv*pexenv)
      WRITE(jjj,'(/,a,/)') 'TOTAL BULK PRESSURE '
      WRITE(jjj,729) 'Ideal pressure       ',pid,0.0
      WRITE(jjj,729) 'Energy con. <E/3V>   ',pexen,pexenv
      WRITE(jjj,729) 'Collision press      ',pcollav,pcollv
      WRITE(jjj,731)
      WRITE(jjj,729) 'Bulk pressure        ' ,ptot,ptotv
      write(jjj,729) 'uuta',uuta(1),0.0

729   FORMAT(a,F12.4,F10.4)
731   FORMAT(70('_'))

      STOP  
      END  



c  ********************************************************************
c  ********************************************************************
c  **                                                                **
c  ** SUBROUTINE EARTH                                               **
c  **                                                                **
c  ** Calculates the average and deviation of nnn values stored in   **
c  ** array a(1..25)                                                 **
c  **                                                                **
c  ********************************************************************
c  ********************************************************************

      SUBROUTINE earth(a,bbbwww,xb,nnn)   
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION a(25)      
      yak=1./(nnn*(nnn-1))
      b=0.               
      xb = 0.
      DO 2 k=1,nnn   
    2   xb=xb+a(k)    
      xb=xb/nnn    
      DO 1 k=1,nnn 
    1   b=b+(a(k)-xb)*(a(k)-xb)  
      b=SQRT(yak*b)           
      bbbwww = b             
      RETURN                
      END                  
      
c  ********************************************************************
c  ********************************************************************
c  **                                                                **
c  ** SUBROUTINE EARTH2                                              **
c  **                                                                **
c  ** Calculates the average and deviation of nnn values stored in   **
c  ** matrix a(1..25,1..num) and stores them in arrays xb and bbbwww **
c  **                                                                **
c  ********************************************************************
c  ********************************************************************

      SUBROUTINE earth2(a,bbbwww,xb,nnn,num)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (mxspec = 10)
      DIMENSION a(25,mxspec),xb(mxspec),bbbwww(mxspec)
      yak=1./(nnn*(nnn-1))   
      DO 3 i=1,num          
      b=0. 
      xb(i) = 0.                
      DO 2 k=1,nnn     
    2   xb(i)=xb(i)+a(k,i) 
      xb(i)=xb(i)/nnn     
      DO 1 k=1,nnn  
    1   b=b+(a(k,i)-xb(i))*(a(k,i)-xb(i))  
      b=SQRT(yak*b)   
    3 bbbwww(i) = b  
      RETURN        
      END          
      

c  ********************************************************************
c  ********************************************************************
c  **                                                                **
c  ** SUBROUTINE SLUMP                                               **
c  **                                                                **
c  ** Set up of the random initial configuration.                    **
c  **                                                                **
c  ********************************************************************
c  ********************************************************************

      SUBROUTINE slump 
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INCLUDE 'bulk.inc'  
 
      WRITE(jjj,*) 'Random configuration generated'
      nsl=0  
      nslmax=100000
      k12=0
    1 nsl=nsl+1
      IF(nsl.GT.nslmax)THEN 
        WRITE(jjj,*)' too dense system'
        STOP 
      ENDIF
      x6tt=(ran2(islu)-.5)*box 
      y6tt=(ran2(islu)-.5)*box
      z6tt=(ran2(islu)-.5)*box
      ispec = ispc(k12+1)
      DO 2 i=1,k12 
        ddx=x6tt-x6(i)
        ddy=y6tt-y6(i)
        ddz=z6tt-z6(i)
        ddx=ddx-AINT(ddx*box2i)*box  
        ddy=ddy-AINT(ddy*box2i)*box 
        ddz=ddz-AINT(ddz*box2i)*box 
        r2=ddx**2+ddy**2+ddz**2 
        IF(r2.LT.hc2v(i,ispec)) GO TO 1
    2 CONTINUE 
      k12=k12+1 
      x6(k12)=x6tt
      y6(k12)=y6tt
      z6(k12)=z6tt 
      IF(k12.EQ.npart) RETURN 
      GO TO 1
      END 

c  ********************************************************************
c  ********************************************************************
c  **                                                                **
c  ** SUBROUTINE QIN                                                 **
c  **                                                                **
c  ** Essential part of the program. Checks for overlap and calculate**
c  ** the energy change after a MC step.                             **
c  **                                                                **
c  ********************************************************************
c  ********************************************************************

      SUBROUTINE qin
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INCLUDE 'bulk.inc'  

      du = 0.d0
      ttx6 = x6(il)
      tty6 = y6(il)
      ttz6 = z6(il)
      do k = 1,il-1
        ddx= ABS(tx6-x6(k))  
        ddy= ABS(ty6-y6(k)) 
        ddz= ABS(tz6-z6(k))
        IF (ddx.GT.box2) ddx=ddx-box
        IF (ddy.GT.box2) ddy=ddy-box
        IF (ddz.GT.box2) ddz=ddz-box
        rt22=ddx*ddx+ddy*ddy+ddz*ddz  
        IF(rt22.LT.hc2v(k,ispec)) RETURN 
        un = chv(k)/sqrt(rt22)
        ddx= ABS(ttx6-x6(k))  
        ddy= ABS(tty6-y6(k)) 
        ddz= ABS(ttz6-z6(k))
        IF (ddx.GT.box2) ddx=ddx-box
        IF (ddy.GT.box2) ddy=ddy-box
        IF (ddz.GT.box2) ddz=ddz-box
        rt22=ddx*ddx+ddy*ddy+ddz*ddz  
        uo = chv(k)/sqrt(rt22)
        du = un-uo+du
      enddo
      do k = il+1,npart
        ddx= ABS(tx6-x6(k))  
        ddy= ABS(ty6-y6(k)) 
        ddz= ABS(tz6-z6(k))
        IF (ddx.GT.box2) ddx=ddx-box
        IF (ddy.GT.box2) ddy=ddy-box
        IF (ddz.GT.box2) ddz=ddz-box
        rt22=ddx*ddx+ddy*ddy+ddz*ddz  
        IF(rt22.LT.hc2v(k,ispec)) RETURN 
        un = chv(k)/sqrt(rt22)
        ddx= ABS(ttx6-x6(k))  
        ddy= ABS(tty6-y6(k)) 
        ddz= ABS(ttz6-z6(k))
        IF (ddx.GT.box2) ddx=ddx-box
        IF (ddy.GT.box2) ddy=ddy-box
        IF (ddz.GT.box2) ddz=ddz-box
        rt22=ddx*ddx+ddy*ddy+ddz*ddz  
        uo = chv(k)/sqrt(rt22)
        du = un-uo+du
      enddo
      isos=0   
c      chil=chv(il) 
c      DO 12 k=1,npart
c   12 ei(k)=
c      ei(il)=0
      du = chv(il)*du
      RETURN       
      END        

c  ********************************************************************
c  ********************************************************************
c  **                                                                **
c  ** SUBROUTINE LIV                                                 **
c  **                                                                **
c  ** Gives the total coulombic energy and the force by recalculating**
c  ** every particleinteraction. Used to check QIN                   **
c  **                                                                **
c  ********************************************************************
c  ********************************************************************

      SUBROUTINE LIV 
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INCLUDE 'bulk.inc'  
c      DIMENSION rw2(1000),rwi(1000)  
      xww1=0.             
      DO 5 i = 1,npart-1 
      tx6 = x6(i) 
      ty6 = y6(i)
      tz6 = z6(i) 
      ip1 = i+1  
      DO 40 k=ip1,npart 
        ddx=tx6-x6(k)  
        ddy=ty6-y6(k) 
        ddz=tz6-z6(k) 
        ddx=ddx-aint(ddx*box2i)*box  
        ddy=ddy-aint(ddy*box2i)*box 
        ddz=ddz-aint(ddz*box2i)*box 
        rw2(k)=ddx*ddx+ddy*ddy+ddz*ddz   
   40 CONTINUE
     
      DO 4 k=ip1,npart
        uj1=chv(i)*chv(k)/sqrt(rw2(k))
        xww1=xww1+uj1 
c        esa(i,k)=uj1 
    4 CONTINUE      
    5 CONTINUE     
c      DO 6 i=1,npart-1 
c      DO 6 k=i+1,npart
c        esa(k,i)=esa(i,k)  
c    6 CONTINUE
c      DO 7 i=1,npart
c      esa(i,i) =0. 
c    7 CONTINUE 
      RETURN       
      END         


c  ********************************************************************
c  **                                                                **
c  ** SUBROUTINE COLLISION                                           **
c  **                                                                **
c  ** This routine checkes if a particle of type i in one halfcel is **
c  ** able to collide with a patricle of type j in the other. The    **
c  ** contribution to the pressure is calculation by intergration    **i
c  ** At a random position around the real particle of type i a ghost**
c  ** particle of type j is inserted and the potential with every    **
c  ** other particle in the system is calculated. This gives the     **
c  ** collision probability                                          **
c  **                                                                **
c  ********************************************************************
c  ********************************************************************

      SUBROUTINE Collision 
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INCLUDE 'bulk.inc'  
      DIMENSION rel(mxspec)
     * ,scoll(25,mxspec,mxspec),coll(mxspec,mxspec)
      DO 1 i=1, mxspec
        DO 1 k=1,mxspec
          coll(i,k)=0.  
            DO 1 j=1,ny3                                                     
              scoll(j,i,k)=0.                                                
   1  CONTINUE                                                                  
      RETURN  
      
c     --------------------------------------------------                        
 
      ENTRY Collision1  
      DO 110 i=1,nwins    
        num=1                                                                  
        DO 140 isp=1,nspec
          nisp = INT(hion(isp,2))
          DO 121 ksp = 1,nspec          
              nnn  = ran2(islu)*nisp + num
              dis  = hion(isp,3) + hion(ksp,3) 
            aa=(2*ran2(islu)-1)*dis 
            v   = 2 * pi * ran2(islu)
            wz6 = z6(nnn) + aa
            g2  = SQRT(dis * dis - aa*aa)+.00001  
            wx6 = x6(nnn) + g2*COS(v)        
            wy6 = y6(nnn) + g2*SIN(v)       
            wx6 = wx6 - AINT(wx6*box2i)*box
            wy6 = wy6 - AINT(wy6*box2i)*box
            wz6 = wz6 - AINT(wz6*box2i)*box
            urej=0                                                                 
            DO 11 k=1, npart
              ddx=ABS(wx6-x6(k))  
              ddy=ABS(wy6-y6(k)) 
              ddz=ABS(wz6-z6(k))
              IF (ddx.GT.box2) ddx = ddx - box
              IF (ddy.GT.box2) ddy = ddy - box
              IF (ddz.GT.box2) ddz = ddz - box
              rw2(k)=ddx*ddx+ddy*ddy+ddz*ddz   
              IF (rw2(k).LT.hc2v(k,ksp)) GO TO 121     
   11       CONTINUE                                                               
c            IF(urej.EQ.0) THEN 
c              CALL vrsqrt(rw2,rwi,npart)  

              do 9453 k=1,npart
 9453         rwi(k)=1./sqrt(rw2(k))
              wtot2=0  
              DO 151 k=1,npart
                wtot2 = wtot2 + rwi(k)*chv(k)
 151          CONTINUE         
              coll(isp,ksp)=EXP(abeta*wtot2*hion(ksp,4)*ecf)
     *                    + coll(isp,ksp)   
c            ENDIF                        
 121      CONTINUE 
          num  = num + nisp
 140    CONTINUE                                                            
 110  CONTINUE 
      RETURN  
c     -------------------------------------------------  
      ENTRY Collision2     
      nwtot = nwins*INT(ny1/nwint)*ny2
      WRITE (kkk,'(/,/,a)') 'COLLISION PRESSURE MATRIX'
      WRITE (kkk,'(/,a,i6)') 'Total collision trials per species ',nwtot
      WRITE (kkk,2010) (i,i=1,nspec)
      DO 1020 k=1, nspec 
        DO 1010 i= 1, nspec                                                    
          scoll(my3,i,k)=coll(i,k)/nwtot
          coll(i,k)=0  
 1010   CONTINUE 
        WRITE(kkk,2012) (scoll(my3,i,k),i=1,nspec)
 1020 CONTINUE       
      RETURN     
2010  FORMAT('Species          ',10(12X,I6),/)   
2011  FORMAT(/,I4,'    Col. samples',10('            ',I6))
2012  FORMAT('        pressure    ',10E18.6,/)
2013  FORMAT(I4,'        ',10E12.5)
c     ------------------------------------------------------                    
      ENTRY collision3    
      WRITE(jjj,'(/,a,/)') 'COLLISION MATRIX  AND RELATIVE ERROR'          
      WRITE (jjj,2010) (i,i=1,nspec)
      DO 1032 i=1,nspec                                                      
        DO 1030 k=1,nspec            
          collav=0                     
          DO 1031 j=1,ny3             
            cwi(j,i,k) =scoll(j,i,k)
1031        dum(j)=scoll(j,i,k)      
          CALL EARTH(dum,collv,collav,ny3) 
          rel(k)=0                     
          IF(collav.NE.0) rel(k)=collv/collav  
          coll(i,k)=caver(i)*collav
1030    CONTINUE 
        WRITE(jjj,2013) i,(coll(i,k),rel(k),k=1,nspec)
1032  CONTINUE 
      WRITE(jjj,'(/)')  
      RETURN   
      END     

      
c  ********************************************************************
c  ********************************************************************
c  **                                                                **
c  ** SUBROUTINE WIDOM                                               **
c  **                                                                **
c  ** Puts a ghost particle anywhere in the cell and calculates the  **
c  ** energy with the rest of the system. Averaging of the exponent  **
c  ** results in the chemical potential. A rescaling of the ion-     **
c  ** charges is neccessary to assure neutrality. The integration    **
c  ** is performed in entry 2 with a simpson's rule.                 **
c  ** The total chemical potential is divided in to the ideal (chid),**
c  ** hard core (chhc) and the electric contribution (chel). For     **
c  ** each species the chem.pot. is given cellaverage, at the wall,  **
c  ** and in the midplan.
c  **                                                                **
c  ********************************************************************
c  ********************************************************************

      SUBROUTINE widom   
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INCLUDE 'bulk.inc'  
      DIMENSION chel(25,mxspec),chhc(25,mxspec),chex(25,mxspec)
     *,chto(25,mxspec),chexw(25,mxspec),dch1(25,mxspec)
     *,dch2(25,mxspec),chint(mxspec,11),ewnom(mxspec,11)
     *,ewden(mxspec,11), chinta(mxspec,11)
     *,expuw(mxspec), chid(mxspec)
     *,chelav(mxspec),chelv(mxspec),chhcav(mxspec),chhcv(mxspec) 
     *,chexav(mxspec),chexv(mxspec),chtoav(mxspec),chtov(mxspec)  
     *,chexwa(mxspec),chexwv(mxspec)
     *,ihc(mxspec),irej(mxspec),mwcn(25),ihcall(0:5)
      CHARACTER*80 str(mxspec)

      DO 11 i=0,nfix                                                            
        DO 16 j=1,nspec    
          IF(caver(j).NE.0) THEN                                               
            chid(j+nspec*i) = DLOG( caver(j) / avno * 1.d27) 
          ELSE  
            chid(j+nspec*i) = -77  
          ENDIF
 16     CONTINUE
 11   CONTINUE  
      DO 12 j=1,mxspec 
        DO 13 I=1,25   
          chel(i,j)=0 
          chhc(i,j)=0 
          chex(i,j)=0 
          chto(i,j)=0 
 13     CONTINUE   
        chelav(j)=0 
        chhcav(j)=0 
        chexav(j)=0 
        chtoav(j)=0 
        expuw(j)=0 
        ihc(j)=0  
        DO 15 k=1,11 
          ewden(j,k)=0 
          ewnom(j,k)=0 
          chinta(j,k)=0 
 15      CONTINUE 
 12   CONTINUE   
      DO 14 j=1,25    
         mwcn(j)=0 
 14   CONTINUE  
      DO 17 j= 0,5
        ihcall(j) =0
 17   CONTINUE
      ntel =0
 
      str(1) = 'Chemical potentials averaged over the whole cell'
      str(2) = 'Chemical potentials at the wall'
      str(3) = 'Chemical potentials at the midplane'

      RETURN 


      ENTRY widom1

      mwcn(my3)=mwcn(my3) + 1 
      DO  127 mp = 0,nfix
        DO 110 i=1,nwins
          x = box*(ran2(islu)-0.5) 
          y = box*(ran2(islu)-0.5)
          z = 0 
          IF (mp.LE.1) THEN
            z = box*(ran2(islu)-0.5)  
            IF (mp.EQ.1) z = SIGN(box2,z)
          ENDIF
          DO 120 j=1,npart  
            ddx = ABS(x  - x6(j))
            IF (ddx.GT.box2) ddx= ddx - box                   
            ddy = ABS(y  - y6(j)) 
            IF (ddy.GT.box2) ddy= ddy - box                   
            ddz = ABS(z  - z6(j))                   
            IF (ddz.GT.box2) ddz= ddz - box                   
            rw2(j)=ddx*ddx + ddy*ddy + ddz*ddz 
 120      CONTINUE 
                             
          irsum=0
          DO 130 j= 1,nspec
            jm = mp*nspec + j
            irej(jm)=0
            DO 131 k=1,npart    
              IF(rw2(k).LT.hc2v(k,j)) irej(jm)=1
 131        CONTINUE
            irsum = irsum + irej(jm) 
 130      CONTINUE 
          IF(irsum.EQ.nspec) THEN
            ihcall(mp) =ihcall(mp)+1
c            WRITE(*,*) 'irsum eq nspec'
            GOTO 110               
          ENDIF 
                     
c          CALL vrsqrt(rw2,rwi,npart) 

          do 9375 k=1,npart
 9375     rwi(k)=1./sqrt(rw2(k))
          wtot2=0     
          wtot3=0   
          nl =  1 
          DO 150 j=1,nspec 
            nh = INT (hion(j,2))
            wsum=0              
            DO 151 k= nl, nl + nh -1      
              wsum=wsum+rwi(k)   
 151        CONTINUE               
            nl = nl + nh
            wtot2=wtot2 + wsum*hion(j,4) 
            wtot3=wtot3 + wsum
 150      CONTINUE  
c          WRITE(*,*) 'first :nergy',wtot2,wtot3

          uj1 = 0.d0


          wtot2 = wtot2 + uj1
c          WRITE(*,*) 'energy',wtot2,wtot3
          DO 160 j=1,nspec 
            jm = mp*nspec + j
            IF(irej(jm).EQ.1)THEN 
              ihc(jm)=ihc(jm)+1
              GOTO 160 
            ENDIF  
            expuw(jm)=expuw(jm)+EXP(abeta*wtot2*hion(j,4)*ecf)  
            DO 170 k1=0,10  
              k=k1+1  
              ew=hion(j,4)*(wtot2-k1*0.1*hion(j,4)*wtot3/npart) 
              ewla=ew*k1*0.1   
              ewd=EXP(abeta*ecf*ewla)  
              ewden(jm,k)=ewden(jm,k)+ewd   
              ewnom(jm,k)=ewnom(jm,k)-ew*abeta*ecf*ewd  
 170        CONTINUE                      
 160      CONTINUE                     
 110    CONTINUE
127   CONTINUE 
      RETURN                     
                
      ENTRY widom2 
      ntocp = nspec*(nfix +1)               
      DO 1025 i=1,ntocp
        DO 1030 j=1,11  
          IF(ewden(i,j).EQ.0)THEN 
            WRITE(JJJ,*)' WIDOM DENOMINATOR EQUALS ZERO',i,j 
          ELSE 
            chint(i,j)=ewnom(i,j)/ewden(i,j)  
            ewnom(i,j)=0 
            ewden(i,j)=0  
            chinta(i,j)=chinta(i,j)+1./ny3*chint(i,j) 
          ENDIF  
 1030   CONTINUE  
        aint4=chint(i,2)+chint(i,4)+chint(i,6)+chint(i,8)+chint(i,10) 
        aint2=chint(i,3)+chint(i,5)+chint(i,7)+chint(i,9)  
        aint1=chint(i,1)+chint(i,11)   
        chel(my3,i)=1./30.*(aint1+2*aint2+4*aint4) 
 1025 CONTINUE  
      nwtot=mwcn(my3)*nwins
      nwtot = nwins*INT(ny1/nwint)*ny2
      DO 1020 i=1,ntocp 
        ihc(i)=ihc(i)+ihcall(INT((i-1)/nspec))
        chhc(my3,i)=-DLOG(DBLE(nwtot-ihc(i))/nwtot) 
        chexw(my3,i)=-DLOG(expuw(i)/nwtot)  
        chex(my3,i)=chel(my3,i)+chhc(my3,i)  
        chto(my3,i)=chex(my3,i)+chid(i)  
        expuw(i)=0   
        ihc(i)=0  
 1020 CONTINUE 
      DO 1035 j=1,nspec   
        dch1(my3,j)=chto(my3,j)-chex(my3,nspec+j)
        dch2(my3,j)=chto(my3,j)-chex(my3,2*nspec +j)
1035  CONTINUE
      WRITE(kkk,2010)nwtot  
      DO 1036 j=0,nfix
        ihcall(j)=0 
        WRITE(kkk,'(/,a)') str(j+1)
        WRITE(kkk,2015)
        WRITE(kkk,2020)(i,chid(i),chhc(my3,i),chel(my3,i),
     *    chex(my3,i),chto(my3,i),chexw(my3,i),i=j*nspec + 1,
     *    (j+1)*nspec) 
1036  CONTINUE 
      RETURN  

      ENTRY widom3
      nwtot=0
      DO 2100 i=1,ny3   
         nwtot=nwtot+mwcn(i)*nwins  
 2100 CONTINUE 
      
      ntocp = nspec*(nfix +1) 
      i=1
      CALL earth2(chel,chelv,chelav,ny3,ntocp)
      call earth2(chexw,chexwv,chexwa,ny3,ntocp) 
      call earth2(chhc,chhcv,chhcav,ny3,ntocp) 
      call earth2(chex,chexv,chexav,ny3,ntocp)  
      call earth2(chto,chtov,chtoav,ny3,ntocp)  
      WRITE (jjj,'(a,/)') 'CONTACT CORRELATION g(r)'
      WRITE (jjj,2001) (i,i=1,nspec)
      pcollav = 0.
      pcollv  = 0.
      DO 2140 i=1,nspec  
        DO 2150  k=1,nspec                                                     
          DO 2160 j= 1,ny3
            dum(j)=cwi(j,i,k)*EXP(chexw(j,k))
2160      CONTINUE  
          CALL EARTH(dum,cwiv,cwiav,ny3)
          cwi(11,i,k) = cwiv
          cwi(12,i,k) = cwiav
2150    CONTINUE 
        WRITE(jjj,2002) i,(cwi(12,i,k),cwi(11,i,k),k=1,nspec)
2140  CONTINUE  
      WRITE (jjj,'(/,a,/)') 'CONTACT PRESSURE MATRIX'
      WRITE (jjj,2001) (i,i=1,nspec)
      DO 2170 i=1,nspec  
        DO 2180  k=1,nspec                                                     
          DO 2190 j= 1,ny3
            dum(j)=2./3.*pi*caver(i)*cwi(j,i,k)*
     *        EXP(chexw(j,k)+chid(k))*(hion(i,3)+hion(k,3))**3
2190      CONTINUE  
          CALL EARTH(dum,cwiv,cwiav,ny3)
          cwi(11,i,k) = cwiv
          cwi(12,i,k) = cwiav
          pcollav = pcollav + cwiav
          pcollv  = pcollv  + cwiv*cwiv
2180    CONTINUE
        WRITE(jjj,2002) i,(cwi(12,i,k),cwi(11,i,k),k=1,nspec)
2170  CONTINUE  
      pcollv = SQRT(pcollv)
      WRITE(jjj,'(/,a,F12.5,F10.5,/)') 'Total Collision pressure   =        
     *  ' ,pcollav,pcollv



      WRITE(jjj,44)
      WRITE(JJJ,2034) 
      WRITE(jjj,2035)((i-1.)/10.,i=1,11) 
      DO 2031 i=1,ntocp
      WRITE(JJJ,2040)i,(chinta(i,j),j=1,11) 
 2031 CONTINUE                             
      WRITE(jjj,2010)nwtot  
      DO 2101 j=0,nfix
        WRITE(jjj,'(/,a)') str(j+1)
        WRITE(jjj,2015)
        WRITE(jjj,2030)(MOD(i-1,nspec)+1,chid(i),chhcav(i),chhcv(i)
     *    ,chelav(i),chelv(i),chexav(i),chexv(i),chtoav(i),chtov(i)
     *    ,chexwa(i) ,chexwv(i),i=j*nspec + 1,nspec*(j+1))
        WRITE(jjj,2033)( EXP((2*chexav(i-2)+chexav(i-1))/3) )
2101  CONTINUE
        

                        
 2001 FORMAT('Species      ',10(I12,12X))
 2002 FORMAT(I4, '        ',10(F12.5,F10.5))
 2010 FORMAT(/,/,'CHEMICAL POTENTIALS IN UNITS OF KT, MEASURED',I8,
     *' TIMES')
 2015 FORMAT(/,'Type  Ideal     Hard Core       El. Static       Excess'
     *,'         Total         Uncorrected Excess') 
 2020 FORMAT(10(I3,F8.4,5(F9.4,'       ')/))     
 2030 FORMAT(10(I3,F8.4,5(F9.4,F7.4)/))         
 2033 FORMAT('2:1 ',10(F8.4)/)
 2034 FORMAT(/,/,'INTEGRAND IN WIDOM CORRECTION FOR DIFFERENT CHARGE PARAMET
     *AMETERS') 
 2035 FORMAT('TYPE',11('  ',F4.2,'  ')/) 
 2040 FORMAT(10(I3,11F8.3)/) 
   44 FORMAT(/,'************************************************************
     **********************************************')
      RETURN  
      END   
