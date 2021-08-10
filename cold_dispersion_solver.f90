!!!!!!!!!!! to solve the cold plasma dispersion relationship!!!!!!!
!!!!!!!!!!! Stix 1992 Waves in plasma: page 8 and 9!!!!!!!!!!
!!!!!!!!!!! by Jingchun Li, Yaoru Qu!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Email: lijc@sustech.edu.cn !!!!!!!!!!!!!!


module parameters
implicit none
real(8),parameter::c=2.99792458d8,epsilon0=8.854187817d-12,qe=1.60217662d-19,mp=1.6726231d-27,b0=1.0d0,pi=3.1415926535898
real(8)::c2,theta,mu0,m_i,m_e,q_i,q_e,n_i,n_e,c2v,mp2me,Oi,Oe
real(8)::Oi2,Oe2,wi,wi2,we2,we4,wi4,wp2,ct2,ww(10,280),w(10,280),wex,wex2,rpol(5,280)
real(8),dimension(280)::kc,k,p0,p2,p4,p6,p8,p10,kc2,kc4
real(8)::pol(11)
integer i,j,z
    end module
 
 
program cold
use parameters
implicit none
interface
    function calw(pol)
    implicit none!                
    real(8)::calw(10),pol(11)
    end function
end interface
real(8)::tem(10),tem2(10),b(11)
integer n,m,L
n=11
m=10!
theta=89.0d0/180.0*pi!theta
c2v=15d0
mp2me=1800

ct2=cos(theta)**2
c2=c*c
m_i=mp
q_i=qe
q_e=-qe
mu0=1/(c2*epsilon0)
n_e=c2v*c2v*epsilon0*b0*b0/mp
n_i=n_e
m_e=m_i/mp2me
Oi=q_i*b0/m_i
Oe=-q_e*b0/m_e
Oi2=Oi**2
Oe2=Oe**2
wi2=n_i*q_i**2/(epsilon0*m_i)
we2=n_e*q_e**2/(epsilon0*m_e)
we4=we2**2
wi4=wi2**2
wi=sqrt(wi2)
wp2=wi2+we2
wex=wp2+Oi*Oe
wex2=wex**2
do i=1,280
kc(i)=10.0**(-12.2+0.02*i)*1d20
end do
kc2=kc*kc
kc4=kc2*kc2

!!Xie
!p0=-kc4*oi2*oe2*wp2*ct2
!p2=kc4*(wp2*(Oe2+Oi2-Oe*Oi)*ct2+Oi*Oe*WEx)+kc2*wp2*Oi*Oe*WEx*(1+ct2)
!!difference is here
!p4=-(kc4*(Oe2+Oi2+wp2)+2*kc2*WEx2+kc2*wp2*(Oe2+Oi2-Oe*Oi)*(1+ct2)+wp2*WEx2)
!p6=kc4+(2*kc2+wp2)*(Oe2+Oi2+2*wp2)+WEx2
!p8=-(2*kc2+oe2+oi2+3*wp2)

!!origin
!p0=-kc4*oi2*oe2*wp2*ct2
!p6=kc4+2*Oi2*we2-3*we4+3*wi4+Oi2*wp2+6*we2*wp2+2*kc2*(Oe2+Oi2+2*wp2)+Oe2*(Oi2-2*we2+3*wp2) !do not affact
!p4=2*Oe2*we4+2*Oe*Oi*we4-Oi2*we4+2*we2*we4-3*Oe2*wi4-3*we2*wi4-wi2*wi4-Oe2*Oi2*wp2-2*Oe2*we2*wp2-2*Oe*Oi*we2*wp2-2*Oi2*we2*wp2&
!&-3*we4*wp2-kc4*(Oe2+Oi2+wp2)-kc2*(-2*we4+2*wi4+4*we2*wp2+Oi2*(3*we2+wp2)+Oe2*(2*oi2-3*we2+4*wp2))&
!&-kc2*(oe2*we2+oi2*wi2)*ct2   ! do not affact
!p2=-2*oe*oi*we2*we4+oe2*we2*wi4+2*oe*oi*we2*wi4+oe2*wi2*wi4+2*oe*oi*we4*wp2+oi2*we4*wp2&
!&+kc4*(oi2*we2+oe2*(oi2+wi2))+kc2*(2*oe*oi*we2*wi2+oi2*we2*(we2+wp2)+oe2*(-we4+2*wi4+oi2*wp2+we2*wp2))&
!&+kc2*(-2*oe*oi*we2*wi2+oi2*we2*wi2+oe2*(-we4+oi2*wp2+we2*wp2)+kc2*(oe2*we2+oi2*wi2))*ct2
!p8=-(2*kc2+oe2+oi2+3*wp2)


!!simplify
p0=-kc4*oi2*oe2*wp2*ct2
p6=kc4+2*kc2*(Oe2+Oi2+2*wp2)
p4=2*Oe2*we4+2*Oe*Oi*we4-Oi2*we4+2*we2*we4-3*Oe2*wi4-Oe2*Oi2*wp2-2*Oe2*we2*wp2-2*Oe*Oi*we2*wp2-2*Oi2*we2*wp2&
&-3*we4*wp2-kc4*(Oe2+Oi2+wp2)-kc2*(-2*we4+2*wi4+4*we2*wp2+Oi2*(3*we2+wp2)+Oe2*(2*oi2-3*we2+4*wp2))&
&-kc2*(oe2*we2+oi2*wi2)*ct2   
p2=kc4*(oi2*we2+oe2*(oi2+wi2))+kc2*(2*oe*oi*we2*wi2+oi2*we2*(we2+wp2)+oe2*(-we4+2*wi4+oi2*wp2+we2*wp2))+&
&kc2*(-2*oe*oi*we2*wi2+oi2*we2*wi2+oe2*(-we4+oi2*wp2+we2*wp2)+kc2*(oe2*we2+oi2*wi2))*ct2
p8=-(2*kc2+oe2+oi2+3*wp2)


rpol(1,:)=p0  ! can delete
rpol(2,:)=p2
rpol(3,:)=p4
rpol(4,:)=p6
rpol(5,:)=p8
do i=1,280
pol=(/1d0,0d0,p8(i),0d0,p6(i),0d0,p4(i),0d0,p2(i),0d0,p0(i)/)
![1,0,p8,0,p6,0,p4,0,p2,0,p1,0] the input
call DSRRT(pol,tem,tem2,n,m,L,b)  !Newton 
    do j=1,10
    WW(j,i)=tem(j)
    end do
end do
W=WW/oi
k=kc/c*sqrt(m_i*c2*epsilon0/n_e/q_e**2)
write(*,*)shape(w)
open(unit=10,file='data.dat')

!write(10,*)pol
!do j=1,5
!write(10,"(280(E28.20E3,'   '))")rpol(j,:)
!end do

write(10,"(280(E28.20E3,'   '))")k
do j=1,10
write(10,"(280(E28.20E3,'   '))")w(j,:)
end do
close(10)
end program cold
    


SUBROUTINE DSRRT(A,XR,XI,N,M,L,B)
real(8)::a(n),xr(m),xi(m),b(n),pp,w,q,p,x,y,x1,y1,g,dx,u,v,pq,dy,g1,t,dd,dc,u1,v1
integer L,k,m,is,n,i
IF (ABS(A(1))+1.0.EQ.1.0) THEN
	L=0
    WRITE(*,5)
    RETURN
END IF
5	FORMAT(1X,'  ERR')
L=1
K=M
IS=0
W=1.0
DO I=1,N
	B(I)=A(I)/A(1)
end do
20	PP=ABS(B(K+1))
IF (PP.LT.1.0E-12) THEN
  XR(K)=0.0
  XI(K)=0.0
  K=K-1
  IF (K.EQ.1) THEN
    XR(K)=-B(2)*W/B(1)
    XI(K)=0.0
    RETURN
  END IF
  GOTO 20
END IF
Q=PP**(1.0/K)
P=Q
W=W*P
DO I=1,K
  B(I+1)=B(I+1)/Q
  Q=Q*P
end do
X=0.0001
X1=X
Y=0.2
Y1=Y
G=1.0d+37
DX=1.0
40	U=B(1)
V=0.0
DO I=1,K
  P=U*X1
  Q=V*Y1
  PQ=(U+V)*(X1+Y1)
  U=P-Q+B(I+1)
  V=PQ-P-Q
end do
G1=U*U+V*V
IF (G1.LT.G) GOTO 105
IF (IS.NE.0) GOTO 80
60	T=T/1.67
X1=X-T*DX
Y1=Y-T*DY
IF (K.GE.50) THEN
  P=SQRT(X1*X1+Y1*Y1)
  Q=EXP(85.0/K)
  IF (P.GE.Q) GOTO 60
END IF
IF (T.GE.1.0E-03) GOTO 40
IF (G.LE.1.0E-18) GOTO 90
65	IS=1
DD=SQRT(DX*DX+DY*DY)
IF (DD.GT.1.0) DD=1.0
DC=6.28/(K+4.5)
70	C=0.0
80	C=C+DC
DX=DD*COS(C)
DY=DD*SIN(C)
X1=X+DX
Y1=Y+DY
IF (C.LE.6.29) GOTO 40
DD=DD/1.67
IF (DD.GT.1.0E-07) GOTO 70
90	IF (ABS(Y).LE.1.0E-06) THEN
	P=-X
    Y=0.0
	Q=0.0
	ELSE
	P=-2.0*X
	Q=X*X+Y*Y
	XR(K)=X*W
	XI(K)=-Y*W
	K=K-1
	END IF
DO I=1,K
  B(I+1)=B(I+1)-B(I)*P
  B(I+2)=B(I+2)-B(I)*Q
end do
XR(K)=X*W
XI(K)=Y*W
K=K-1
IF (K.EQ.1) THEN
  XR(K)=-B(2)*W/B(1)
  XI(K)=0.0
  RETURN
END IF
GOTO 20
105	G=G1
X=X1
Y=Y1
IS=0
IF (G.LE.1.0E-22) GOTO 90
U1=K*B(1)
V1=0.0
DO I=2,K
  P=U1*X
  Q=V1*Y
  PQ=(U1+V1)*(X+Y)
  U1=P-Q+(K-I+1)*B(I)
  V1=PQ-P-Q
end do
P=U1*U1+V1*V1
IF (P.LE.1.0E-20) GOTO 65
DX=(U*U1+V*V1)/P
DY=(U1*V-V1*U)/P
T=1.0+4.0/K
GOTO 60
END
