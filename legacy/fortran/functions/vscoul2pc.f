************************************************************************
*                                                                      *
************************************************************************
      Subroutine GenTab
      Implicit Real*8 (a-h,o-z)
      dimension ix0(1:2),ix1(1:2)
      equivalence (x0,ix0)
      equivalence (x1,ix1)
      Include 'tabdef.inc'
*     Data fInfty/z7ff0000000000000/
      fnc(x)=exp(-sqrt(x))/sqrt(x)
      x0=0.0d0
      x1=0.0d0 
      ix0(2)=ind0
      ix1(2)=ind1
      index=0
      Do 100 i=0,nPol
*        c(i)=fInfty
         c(i)=0.0d0
100   Continue
      c(1) =-fnc(x0)*(1.0d0/sqrt(x0)+1.0d0/x0)
      c(0) = fnc(x0)-c(1)*x0
      index=index+nPol+1
      ErrMax=0.0d0
      Write(*,'(a,2f12.8)') ' Smallest argument/value:',x0,fnc(x0)
      Write(*,'(a,2f12.8)') ' Largest argument/value: ',x1,fnc(x1)
      Write(*,'(a,i12)') ' Polynomial degree:      ',nPol
      Do 200 i=ind0,ind1,ind2
         ix0(2)=i
         ix1(2)=i+ind2
         Call MkFit(x0,x1,c(index),nPol)
         Call MkErr(x0,x1,c(index),nPol,err)
*        Write(*,'(1x,a,2f12.6,d10.3)') ' *** Intervall:   ',x0,x1,err
         ErrMax=Max(ErrMax,err)
         index=index+nPol+1
200   Continue
      index=index-nPol-1
      Do 300 i=0,nPol
         c(i+index)=0.0d0
300   Continue
      index=index+nPol+1
      Write(*,'(a,i12)') ' Size of table:          ',index 
      Write(*,'(a,e10.3)') ' Max error of algorithm:',ErrMax
      Return
      End

************************************************************************
*                                                                      *
************************************************************************
      Subroutine MkErr(x0,x1,Coef,nPol,ErrMax)
      Implicit Real*8 (a-h,o-z)
      Dimension Coef(0:nPol)
      fnc(x)=exp(-sqrt(x))/sqrt(x)
      ErrMax=0.0d0
      Do 300 x=x0,x1,0.01d0*(x1-x0)
         Exact=fnc(x)
         t=0.0d0
         Do 310 k=nPol,0,-1
            t=t*x+Coef(k)
310      Continue
         Fit=t
         Err=Exact-Fit
*        Write(*,'(1x,3f14.8,d10.3)') x,Exact,Fit,Err
         ErrMax=Max(ErrMax,Abs(Err))
300   Continue
*     Write(*,'(1x,a,d12.3)') '     Err:         ',ErrMax
      Return
      End

************************************************************************
*                                                                      *
************************************************************************
      Subroutine MkFit(x0,x1,Coef,nPol)
      Implicit Real*8 (a-h,o-z)
      Parameter (MxPol=5)
      Dimension Coef(0:nPol)
      Dimension Hess(0:MxPol,0:MxPol)
      Dimension Absc(0:MxPol)
      Dimension Fit(0:MxPol),Exact(0:MxPol)
      fnc(x)=exp(-sqrt(x))/sqrt(x)
      pi=4.0d0*atan(1.0d0)
      Do 100 i=0,nPol
         Absc(i)=0.5d0*(x0+x1)-0.5d0*(x1-x0)*Cos(pi*(2*i+1)/(2*nPol+2))
100   Continue
      Do 200 i=0,nPol
         t=1.0d0
         Do 210 j=0,nPol
            Hess(i,j)=t
            t=t*Absc(i)
210      Continue
         Coef(i)=fnc(Absc(i))
200   Continue
      Call Gauss(nPol+1,MxPol+1,Hess,Coef,Coef)
*     Write(*,'(1x,a,5f12.6)') '     Abscissae:   ',(Absc(i),i=0,nPol)
*     Write(*,'(1x,a,5f12.6)') '     Coefficients:',(Coef(i),i=0,nPol)
      Do 300 i=0,nPol
         Exact(i)=fnc(Absc(i))
         t=0.0d0
         Do 310 k=nPol,0,-1
            t=t*Absc(i)+Coef(k)
310      Continue
         Fit(i)=t
300   Continue
*     Write(*,'(1x,a,5f12.6)') '     Exact:       ',(Exact(i),i=0,nPol)
*     Write(*,'(1x,a,5f12.6)') '     Fit:         ',(Fit(i),i=0,nPol)
      Return
      End

************************************************************************
*                                                                      *
************************************************************************
c@Process opt(3) source list 
      Subroutine vscoul(x,y,ix,n)
      Implicit Real*8 (a-h,o-z)
      Dimension x(n),y(n),ix(2,n)
      Include 'tabdef.inc'
      IndMax=iShft(ind1-indx,-nShift)
      IndMin=0
      Do 100 i=1,n
c        index=iShft(ix(1,i)-indx,-nShift)
         index=iShft(ix(2,i)-indx,-nShift)
c        If(ix(1,i).lt.ind0) Then
         If(ix(2,i).lt.ind0) Then
            index=IndMin
c        Else If(ix(1,i).gt.ind1) Then
         Else If(ix(2,i).gt.ind1) Then
            index=IndMax
         End If
         index=Max(Min(index,IndMax),IndMin)
         m=(nPol+1)*index
         p=x(i)
*        t=c(m)+p*(c(m+1)+p*(c(m+2)+p*(c(m+3)+p*c(m+4))))
         t=c(m+4)
         t=p*t+c(m+3)
         t=p*t+c(m+2)
         t=p*t+c(m+1)
         t=p*t+c(m)
         y(i)=t
100   Continue
      Return
      End
