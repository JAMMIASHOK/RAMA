
MODULE main

CONTAINS

SUBROUTINE bandred(a,d,e)
!
! This subroutine transforms a real symmetric band matrix a,
! of order n and band width iw, to tridiagonal form by an appropriate
! sequence of Jacobi rotations. During the transformation the
! property of the band matrix is maintained. The method yields
! a tridiagonal matrix, the diagonal elements of which are in
! d(n) and off-diagonal elements in e(n).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN OUT)::a(:,:)
 REAL(iwp),INTENT(OUT)::d(0:),e(0:)
 INTEGER::iw,n2,n,k,maxr,irr,ir,kr,j,jm,iugl,j2,l,jl,maxl,i
 REAL(iwp)::g,b,s,c,c2,s2,cs,u,u1,zero=0.0_iwp,one=1.0_iwp,two=2.0_iwp
 n=UBOUND(a,1)
 iw=UBOUND(a,2)-1
 n2=n-2
 IF(n2>=1)THEN
   loop_1: DO k=1,n2
     maxr=iw
     IF(n-k<iw)maxr=n-k
     loop_2: DO irr=2,maxr
       ir=2+maxr-irr
       kr=k+ir
       loop_3: DO j=kr,n,iw
         IF(j/=kr)THEN
           IF(ABS(g)<TINY(g))EXIT loop_3
           jm=j-iw
           b=-a(jm-1,iw+1)/g
           iugl=j-iw
         ELSE
           IF(ABS(a(k,ir+1))<TINY(a(k,ir+1)))EXIT loop_3
           b=-a(k,ir)/a(k,ir+1)
           iugl=k
         END IF
         s=one/SQRT(one+b*b)
         c=b*s
         c2=c*c
         s2=s*s
         cs=c*s
         u=c2*a(j-1,1)-two*cs*a(j-1,2)+s2*a(j,1)
         u1=s2*a(j-1,1)+two*cs*a(j-1,2)+c2*a(j,1)
         a(j-1,2)=cs*(a(j-1,1)-a(j,1))+(c2-s2)*a(j-1,2)
         a(j-1,1)=u
         a(j,1)=u1
         j2=j-2
         DO l=iugl,j2
           jl=j-l
           u=c*a(l,jl)-s*a(l,jl+1)
           a(l,jl+1)=s*a(l,jl)+c*a(l,jl+1)
           a(l,jl)=u
         END DO
         jm=j-iw
         IF(j/=kr)a(jm-1,iw+1)=c*a(jm-1,iw+1)-s*g
         maxl=iw-1
         IF(n-j<iw-1)maxl=n-j
         IF(maxl>0)THEN
           DO l=1,maxl
             u=c*a(j-1,l+2)-s*a(j,l+1)
             a(j,l+1)=s*a(j-1,l+2)+c*a(j,l+1)
             a(j-1,l+2)=u
           END DO
         END IF
         IF(j+iw<=n)THEN
           g=-s*a(j,iw+1)
           a(j,iw+1)=c*a(j,iw+1)
         END IF
       END DO loop_3
     END DO loop_2
   END DO loop_1
 END IF
 e(1)=zero
 d(1:n)=a(1:n,1)
 IF(2<=n)THEN
   DO i=2,n
     e(i)=a(i-1,2)
   END DO
 END IF
RETURN
END SUBROUTINE bandred

FUNCTION bandwidth(g) RESULT(nband)
!
! This function finds the element bandwidth from g.
!
 IMPLICIT NONE     
 INTEGER,INTENT(IN)::g(:) 
 INTEGER::nband
 nband=MAXVAL(g,1,g>0)-MINVAL(g,1,g>0)
RETURN
END FUNCTION bandwidth

SUBROUTINE banmul(kb,loads,ans)
! This subroutine multiplies a symmetrical band kb by the vector loads.
! 
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::kb(:,:),loads(0:)
 REAL(iwp),INTENT(OUT)::ans(0:)
 INTEGER::neq,nband,i,j 
 REAL(iwp)::x,zero=0.0_iwp
 neq=UBOUND(kb,1) 
 nband=UBOUND(kb,2)-1
 DO i=1,neq
   x=zero
   DO j=nband+1,1,-1
     IF(i+j>nband+1)x=x+kb(i,j)*loads(i+j-nband-1)
   END DO
   DO j=nband,1,-1
     IF(i-j<neq-nband)x=x+kb(i-j+nband+1,j)*loads(i-j+nband+1)
   END DO
   ans(i)=x
 END DO
RETURN
END SUBROUTINE banmul
SUBROUTINE bantmul(kb,loads,ans)
!
! This subroutine multiplies an unsymmetrical band kb by the vector loads.
! Could be much improved for vector processors.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::kb(:,:),loads(0:)
 REAL(iwp),INTENT(OUT)::ans(0:)
 INTEGER::i,j,k,l,m,n,iw
 REAL(iwp)::x,zero=0.0_iwp
 n=SIZE(kb,1)
 l=SIZE(kb,2)
 iw=(l-1)/2
 DO i=1,n
   x=zero
   k=iw+2
   DO j=1,l
     k=k-1
     m=i-k+1
     IF(m<=n.AND.m>=1)x=x+kb(i,j)*loads(m)
   END DO
   ans(i)=x
 END DO
RETURN
END SUBROUTINE bantmul 
SUBROUTINE beamdis(loads,nf,ratmax,interp,nels,ell,argv,nlen,ips)
!
! This subroutine produces a PostScript output file "*.dis" displaying
! the deformed 1-D finite element mesh.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::loads(0:),ratmax,ell(:)
 INTEGER,INTENT(IN)::ips,nf(:,:),nels,interp,nlen
 REAL(iwp)::width,height,scale=72,sxy,xo,yo,x,y,dismag,vmax
 REAL(iwp)::xmin,xmax,ymin,ymax,xnow,dmax,zero=0.0_iwp,pt5=0.5_iwp,       &
   opt5=1.5_iwp,fpt5=5.5_iwp,d8=8.0_iwp,ept5=8.5_iwp,d11=11.0_iwp,        &
   one=1.0_iwp,two=2.0_iwp,thr=3.0_iwp,                                   &
   localx,globalx,na,nb,nc,nd,ll,wold,wnew
 INTEGER::i,ii,j,jj,nn,nel,nod
 CHARACTER(LEN=15)::argv
 REAL(iwp),ALLOCATABLE::xcoord(:)
 OPEN(ips,FILE=argv(1:nlen)//'.dis')
!
 nn=nels+1
 ALLOCATE(xcoord(nn))

 xmin=zero
 xmax=zero
 xnow=zero
 ymin=zero
 ymax=zero

 xcoord(1)=xnow

 DO i=2,nn
   xnow=xnow+ell(i-1)
   xcoord(i)=xnow
 END DO

 xmax=xcoord(nn)

 DO i=1,nn
    IF(loads(nf(1,i))<ymin)ymin=loads(nf(1,i))
    IF(loads(nf(1,i))>ymax)ymax=loads(nf(1,i))
 END DO

 width=xmax-xmin
 height=ymax-ymin
 dmax=ratmax*width
 IF(height>width)dmax=ratmax*height
!
 vmax=zero
 DO i=1,nn
     IF(ABS(loads(nf(1,i)))>vmax)vmax=ABS(loads(nf(1,i)))
 END DO
 dismag=dmax/vmax
!
 ymin=zero
 ymax=zero

 DO i=1,nn
   IF(dismag*loads(nf(1,i))<ymin)                            &
     ymin=dismag*loads(nf(1,i))
   IF(dismag*loads(nf(1,i))>ymax)                            &
     ymax=dismag*loads(nf(1,i))
!
   IF(loads(nf(1,i))<ymin)ymin=loads(nf(1,i))
   IF(loads(nf(1,i))>ymax)ymax=loads(nf(1,i))
 END DO
!
 width =xmax-xmin
 height=ymax-ymin
!
!                       allow 1.5" margin minimum on each side of figure
!
!                       portrait mode 
!
 IF(height.GE.d11/ept5*width)THEN
!
!                       height governs the scale
!
   sxy=scale*d8/height
   xo=scale*pt5*(ept5-d8*width/height)
   yo=scale*opt5
 ELSE
!
!                       width governs the scale
!
   sxy=scale*fpt5/width
   xo=scale*opt5
   yo=scale*pt5*(d11-fpt5*height/width)
 END IF
!
!
!                       start PostScript output
!
 WRITE(ips,'(a)')'%!PS-Adobe-1.0'
 WRITE(ips,'(a)')'%%DocumentFonts: none'
 WRITE(ips,'(a)')'%%Pages: 1'
 WRITE(ips,'(a)')'%%EndComments'
 WRITE(ips,'(a)')'/m {moveto} def'
 WRITE(ips,'(a)')'/l {lineto} def'
 WRITE(ips,'(a)')'/s {stroke} def'
 WRITE(ips,'(a)')'/c {closepath} def'
 WRITE(ips,'(a)')'%%EndProlog'
 WRITE(ips,'(a)')'%%Page: 0 1'
 WRITE(ips,'(a)')'gsave'
!
!                       draw the deformed mesh
!
 WRITE(ips,'(2f9.2,a)') xo, yo, ' translate'
 WRITE(ips,'(f9.2,a)') 0.5, ' setlinewidth'

 wnew=loads(nf(1,1))

 DO i=1,nels
    DO j=1,interp

    wold=wnew
    localx=j*ell(i)/interp
    globalx=localx+xcoord(i)
    ll=ell(i)
    na=(one/(ll**3))*(ll**3-thr*ll*localx**2+two*localx**3)
    nb=(one/(ll**2))*(localx*ll**2-two*ll*localx**2+localx**3)
    nc=(one/(ll**3))*(thr*ll*localx**2-two*localx**3)
    nd=(one/(ll**2))*(localx**3-ll*localx**2)
    wnew=na*loads(nf(1,i))+nb*loads(nf(2,i))+nc*loads(nf(1,1+i))+nd*loads(nf(2,i+1))

    x=sxy*((globalx-ell(i)/interp)-xmin)
    y=sxy*(dismag*wold-ymin)
    WRITE(ips,'(2f9.2,a)') x, y,' m'
    x=sxy*(globalx-xmin)
    y=sxy*(dismag*wnew-ymin)
    WRITE(ips,'(2f9.2,a)') x, y,' l'
    WRITE(ips,'(a)')'c s'

    END DO

 END DO
!
 WRITE(ips,'(a)')'grestore'
 WRITE(ips,'(a)')'showpage'
 CLOSE(ips)
!
RETURN
END SUBROUTINE beamdis

SUBROUTINE beam_ge(ge,ell)
!
! This subroutine forms the beam geometric matrix for stability analysis.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::ell
 REAL(iwp),INTENT(OUT)::ge(:,:)
 REAL(iwp)::pt1=0.1_iwp,opt2=1.2_iwp,two=2.0_iwp,d15=15.0_iwp,d30=30.0_iwp
 ge(1,1)=opt2/ell
 ge(1,2)=pt1
 ge(2,1)=pt1
 ge(1,3)=-opt2/ell
 ge(3,1)=-opt2/ell
 ge(1,4)=pt1
 ge(4,1)=pt1
 ge(2,2)=two*ell/d15
 ge(2,3)=-pt1
 ge(3,2)=-pt1
 ge(2,4)=-ell/d30
 ge(4,2)=-ell/d30
 ge(3,3)=opt2/ell
 ge(3,4)=-pt1
 ge(4,3)=-pt1
 ge(4,4)=two*ell/d15
RETURN
END SUBROUTINE beam_ge
SUBROUTINE beam_km(km,ei,ell)
!
! This subroutine forms the stiffness matrix of a
! beam element (bending only).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::ei,ell
 REAL(iwp),INTENT(OUT)::km(:,:)
 REAL(iwp)::two=2.0_iwp,d4=4.0_iwp,d6=6.0_iwp,d12=12.0_iwp
 km(1,1)=d12*ei/(ell*ell*ell) 
 km(3,3)=km(1,1)
 km(1,2)=d6*ei/(ell*ell) 
 km(2,1)=km(1,2) 
 km(1,4)=km(1,2)
 km(4,1)=km(1,4) 
 km(1,3)=-km(1,1) 
 km(3,1)=km(1,3) 
 km(3,4)=-km(1,2)
 km(4,3)=km(3,4) 
 km(2,3)=km(3,4) 
 km(3,2)=km(2,3)
 km(2,2)=d4*ei/ell
 km(4,4)=km(2,2) 
 km(2,4)=two*ei/ell 
 km(4,2)=km(2,4)
RETURN
END SUBROUTINE beam_km                                                      
SUBROUTINE beam_mm(mm,fs,ell)
!
! This subroutine forms the consistent mass matrix of a beam element.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::fs,ell
 REAL(iwp),INTENT(OUT)::mm(:,:)
 REAL(iwp)::fac
 fac=(fs*ell)/420.0_iwp 
 mm(1,1)=156.0_iwp*fac
 mm(3,3)=mm(1,1)
 mm(1,2)=22.0_iwp*ell*fac
 mm(2,1)=mm(1,2)
 mm(3,4)=-mm(1,2)
 mm(4,3)=mm(3,4)
 mm(1,3)=54.0_iwp*fac
 mm(3,1)=mm(1,3)
 mm(1,4)=-13.0_iwp*ell*fac
 mm(4,1)=mm(1,4)
 mm(2,3)=-mm(1,4)
 mm(3,2)=mm(2,3)
 mm(2,2)=4.0_iwp*(ell**2)*fac
 mm(4,4)=mm(2,2)
 mm(2,4)=-3.0_iwp*(ell**2)*fac
 mm(4,2)=mm(2,4)
RETURN
END SUBROUTINE beam_mm  
SUBROUTINE bee8(bee,coord,xi,eta,det)
!
! Analytical version of the bee matrix for an 8-node quad      
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::coord(:,:),xi,eta
 REAL(iwp),INTENT(OUT)::bee(:,:),det
 REAL::x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6,y7,y8,xi2,xi3,et2,et3,   &
   xi2et2,xi1et1,xi2et1,xi1et2,xy12,xy13,xy14,xy15,xy16,xy17,xy18,xy23,   &
   xy24,xy25,xy26,xy27,xy28,xy34,xy35,xy36,xy37,xy38,xy45,xy46,xy47,xy48, &
   xy56,xy57,xy58,xy67,xy68,xy78,xcu,xsq,xon,ecu,esq,eon,x2e2,x1e1,x2e1,  &
   x1e2,cons,bot,top,dn1dx,dn2dx,dn3dx,dn4dx,dn5dx,dn6dx,dn7dx,dn8dx,     &
   dn1dy,dn2dy,dn3dy,dn4dy,dn5dy,dn6dy,dn7dy,dn8dy,zero=0.0_iwp,          &
   one=1.0_iwp,two=2.0_iwp,d3=3.0_iwp,d4=4.0_iwp,d5=5.0_iwp,d6=6.0_iwp,   &
   d8=8.0_iwp,pt125=0.125_iwp
 x1=coord(1,1)
 x2=coord(2,1)
 x3=coord(3,1)
 x4=coord(4,1)
 x5=coord(5,1)
 x6=coord(6,1)
 x7=coord(7,1)
 x8=coord(8,1)
 y1=coord(1,2)
 y2=coord(2,2)
 y3=coord(3,2)
 y4=coord(4,2)
 y5=coord(5,2)
 y6=coord(6,2)
 y7=coord(7,2)
 y8=coord(8,2)
 xi2=xi*xi
 xi3=xi2*xi
 et2=eta*eta
 et3=et2*eta
 xi2et2=xi2*et2
 xi1et1=xi*eta
 xi2et1=xi2*eta
 xi1et2=xi*et2
 xy12=x1*y2-y1*x2
 xy13=x1*y3-y1*x3
 xy14=x1*y4-y1*x4
 xy15=x1*y5-y1*x5
 xy16=x1*y6-y1*x6
 xy17=x1*y7-y1*x7
 xy18=x1*y8-y1*x8
 xy23=x2*y3-y2*x3
 xy24=x2*y4-y2*x4
 xy25=x2*y5-y2*x5
 xy26=x2*y6-y2*x6
 xy27=x2*y7-y2*x7
 xy28=x2*y8-y2*x8
 xy34=x3*y4-y3*x4
 xy35=x3*y5-y3*x5
 xy36=x3*y6-y3*x6
 xy37=x3*y7-y3*x7
 xy38=x3*y8-y3*x8
 xy45=x4*y5-y4*x5
 xy46=x4*y6-y4*x6
 xy47=x4*y7-y4*x7
 xy48=x4*y8-y4*x8
 xy56=x5*y6-y5*x6
 xy57=x5*y7-y5*x7
 xy58=x5*y8-y5*x8
 xy67=x6*y7-y6*x7
 xy68=x6*y8-y6*x8
 xy78=x7*y8-y7*x8
 xcu=-d8*xy48+d4*(-xy14+xy38+xy47+xy58)+two*(xy13+xy15-xy37-xy57)
 xsq=two*(-xy13+xy14-xy17+xy18+xy24-xy28-xy34+xy35-xy38-xy45+xy46+xy47-   &
   xy57+xy58+xy68-xy78)+one*(-xy12+xy16-xy23-xy25+xy27-xy36-xy56-xy67) 
 xon=d8*xy48+two*(xy14-xy18+xy34-xy38-xy45-xy47-xy58-xy78)+xy12-xy16+xy23-&
   xy25+xy27+xy36-xy56-xy67
 ecu=-d8*xy26+d4*(xy16+xy25+xy36+xy27)+two*(-xy15-xy17-xy35-xy37)
 esq=two*(-xy12+xy13-xy16+xy17-xy23+xy24+xy25-xy27-xy28-xy35+xy36+xy46-   &
   xy56+xy57-xy67+xy68)-xy14+xy18-xy34+xy38-xy45-xy47-xy58-xy78
 eon=d8*xy26+two*( xy12-xy16-xy23-xy25-xy27-xy36-xy56+xy67)+xy14-xy18-    &
   xy34+xy38-xy45+xy47-xy58+xy78
 x2e2=d6*(xy24-xy28+xy46+xy68)+d3*(-xy12+xy13-xy14+xy16-xy17+xy18-xy23-   &
   xy25+xy27-xy34+xy35-xy36+xy38-xy45-xy47-xy56+xy57-xy58-xy67-xy78)
 x1e1=d8*(-xy24-xy28+xy46-xy68)+d6*(-xy12+xy18+xy23+xy34-xy45-xy56+xy67+  &
   xy78)+two*(xy14-xy16+xy25+xy27-xy36+xy38-xy47+xy58)
 x2e1=d8*(xy24+xy28+xy46-xy68)+d5*( xy17-xy18-xy34+xy35-xy45+xy78)+d4*    &
   (xy12-xy16-xy23-xy25-xy27-xy36-xy56+xy67)+d3*(-xy14+xy15+xy37-xy38-    &
   xy47+xy58)
 x1e2=d8*(-xy24+xy28+xy46+xy68)+d5*(xy12-xy13+xy23+xy57-xy56-xy67)+d4*    &
   (xy14-xy18+xy34-xy38-xy45-xy47-xy58-xy78)+d3*(-xy15+xy16+xy25-xy27-    &
   xy36+xy37)
 cons= two*(-xy24+xy28-xy46-xy68)
!
 bot=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 bot=bot+xi2et2*x2e2+xi1et1*x1e1
 bot=bot+xi2et1*x2e1+xi1et2*x1e2+cons
!
 det=pt125*bot
!
 xcu=two*y3-d4*y4+two*y5
 xsq=-y2-two*y3+two*y4+y6-two*y7+two*y8
 xon=y2+two*y4 -y6-two*y8
 ecu=-two*y5+d4*y6-two*y7
 esq=-two*y2+two*y3-y4-two*y6+two*y7+y8
 eon=two*y2+y4-two*y6-y8
 x2e2=-d3*y2+d3*y3-d3*y4+d3*y6-d3*y7+d3*y8
 x1e1=-d6*y2+two*y4-two*y6+d6*y8
 x2e1=d4*y2-d3*y4+d3*y5-d4*y6+d5*y7-d5*y8
 x1e2=d5*y2-d5*y3+d4*y4-d3*y5+d3*y6-d4*y8
 cons=zero
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn1dx=top/bot
!
 xcu=zero
 xsq=y1-y3+two*y4-y5+y7-two*y8
 xon=-y1+y3-y5+y7      
 ecu=d4*y5-d8*y6+d4*y7
 esq=two*y1-two*y3+two*y4+two*y5-two*y7-two*y8
 eon=-two*y1-two*y3-two*y5+d8*y6-two*y7      
 x2e2= d3*y1-d3*y3+d6*y4-d3*y5+d3*y7-d6*y8
 x1e1= d6*y1+d6*y3-d8*y4+two*y5+two*y7-d8*y8
 x2e1=-d4*y1-d4*y3+d8*y4-d4*y5-d4*y7+d8*y8
 x1e2=-d5*y1+d5*y3-d8*y4+d3*y5-d3*y7+d8*y8
 cons=-two*y4+two*y8
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn2dx=top/bot
!
 xcu=-two*y1-two*y7+d4*y8
 xsq=two*y1+y2-two*y4+two*y5-y6-two*y8
 xon=-y2+two*y4+y6-two*y8
 ecu=-two*y5+d4*y6-two*y7
 esq=-two*y1+two*y2-y4-two*y5+two*y6+y8
 eon=two*y2-y4-two*y6+y8
 x2e2=-d3*y1+d3*y2-d3*y4+d3*y5-d3*y6+d3*y8
 x1e1=-d6*y2+d6*y4-two*y6+two*y8
 x2e1=+d4*y2-d5*y4+d5*y5-d4*y6+d3*y7-d3*y8
 x1e2=d5*y1-d5*y2+d4*y4-d3*y6+d3*y7-d4*y8
 cons=zero
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn3dx=top/bot
!
 xcu=d4*y1+d4*y7-d8*y8
 xsq=-two*y1-two*y2+two*y3-two*y5+two*y6+two*y7
 xon=-two*y1-two*y3-two*y5-two*y7+d8*y8
 ecu=zero
 esq=y1-two*y2+y3-y5+two*y6-y7
 eon=-y1+y3-y5+y7 
 x2e2=d3*y1-d6*y2+d3*y3-d3*y5+d6*y6-d3*y7  
 x1e1=-two*y1+d8*y2-d6*y3-d6*y5+d8*y6-two*y7 
 x2e1=d3*y1-d8*y2+d5*y3-d5*y5+d8*y6-d3*y7 
 x1e2=-d4*y1+d8*y2-d4*y3-d4*y5+d8*y6-d4*y7  
 cons=two*y2-two*y6
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn4dx=top/bot
!
 xcu=-two*y1-two*y7+d4*y8
 xsq=y2-two*y3+two*y4-y6-two*y7+two*y8
 xon=y2 +two*y4-y6-two*y8
 ecu=two*y1-d4*y2+two*y3
 esq=-two*y2+two*y3+y4-two*y6+two*y7-y8
 eon=two*y2+y4-two*y6-y8
 x2e2=d3*y2-d3*y3+d3*y4-d3*y6+d3*y7-d3*y8 
 x1e1=-two*y2+d6*y4-d6*y6+two*y8
 x2e1=-d3*y1+d4*y2-d5*y3+d5*y4-d4*y6+d3*y8
 x1e2=d3*y1-d3*y2+d4*y4-d5*y6+d5*y7-d4*y8 
 cons=zero
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn5dx=top/bot
!
 xcu=zero     
 xsq=-y1+y3-two*y4+y5-y7+two*y8
 xon=y1-y3+y5-y7 
 ecu=-d4*y1+d8*y2-d4*y3
 esq=two*y1-two*y3-two*y4+two*y5-two*y7+two*y8
 eon=two*y1-d8*y2+two*y3+two*y5+two*y7 
 x2e2=-d3*y1+d3*y3-d6*y4+d3*y5-d3*y7+d6*y8 
 x1e1=two*y1+two*y3-d8*y4+d6*y5+d6*y7-d8*y8
 x2e1=d4*y1+d4*y3-d8*y4+d4*y5+d4*y7-d8*y8
 x1e2=-d3*y1+d3*y3-d8*y4+d5*y5-d5*y7+d8*y8 
 cons=two*y4-two*y8
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn6dx=top/bot
!
 xcu=two*y3-d4*y4+two*y5         
 xsq=two*y1-y2-two*y4+two*y5+y6-two*y8
 xon=-y2+two*y4+y6-two*y8
 ecu=two*y1-d4*y2+two*y3
 esq=-two*y1+two*y2+y4-two*y5+two*y6-y8
 eon=two*y2-y4-two*y6+y8
 x2e2=d3*y1-d3*y2+d3*y4-d3*y5+d3*y6-d3*y8 
 x1e1=-two*y2+two*y4-d6*y6+d6*y8
 x2e1=-d5*y1+d4*y2-d3*y3+d3*y4-d4*y6+d5*y8
 x1e2=d3*y2-d3*y3+d4*y4-d5*y5+d5*y6-d4*y8 
 cons=zero
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn7dx=top/bot
!
 xcu=-d4*y3+d8*y4-d4*y5         
 xsq=-two*y1+two*y2+two*y3-two*y5-two*y6+two*y7 
 xon=two*y1+two*y3-d8*y4+two*y5+two*y7 
 ecu=zero      
 esq=-y1+two*y2-y3+y5-two*y6+y7 
 eon=y1-y3+y5-y7  
 x2e2=-d3*y1+d6*y2-d3*y3+d3*y5-d6*y6+d3*y7  
 x1e1=-d6*y1+d8*y2-two*y3-two*y5+d8*y6-d6*y7  
 x2e1=d5*y1-d8*y2+d3*y3-d3*y5+d8*y6-d5*y7  
 x1e2=d4*y1-d8*y2+d4*y3+d4*y5-d8*y6+d4*y7  
 cons=-two*y2+two*y6
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn8dx=top/bot
!
 y1=-x1
 y2=-x2
 y3=-x3
 y4=-x4
 y5=-x5
 y6=-x6
 y7=-x7
 y8=-x8
!
 xcu=two*y3-d4*y4+two*y5
 xsq=-y2-two*y3+two*y4+y6-two*y7+two*y8
 xon=y2+two*y4-y6-two*y8
 ecu=-two*y5+d4*y6-two*y7
 esq=-two*y2+two*y3-y4-two*y6+two*y7+y8
 eon=two*y2+y4-two*y6-y8
 x2e2=-d3*y2+d3*y3-d3*y4 +d3*y6-d3*y7+d3*y8
 x1e1=-d6*y2+two*y4-two*y6+d6*y8
 x2e1=d4*y2-d3*y4+d3*y5-d4*y6+d5*y7-d5*y8
 x1e2=d5*y2-d5*y3+d4*y4-d3*y5+d3*y6-d4*y8
 cons=zero
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn1dy=top/bot
!
 xcu=zero
 xsq=y1-y3+two*y4-y5+y7-two*y8
 xon=-y1+y3-y5+y7 
 ecu=d4*y5-d8*y6+d4*y7
 esq=two*y1-two*y3+two*y4+two*y5-two*y7-two*y8
 eon=-two*y1-two*y3-two*y5+d8*y6-two*y7 
 x2e2=d3*y1-d3*y3+d6*y4-d3*y5+d3*y7-d6*y8
 x1e1=d6*y1+d6*y3-d8*y4+two*y5+two*y7-d8*y8
 x2e1=-d4*y1-d4*y3+d8*y4-d4*y5-d4*y7+d8*y8
 x1e2=-d5*y1+d5*y3-d8*y4+d3*y5-d3*y7+d8*y8
 cons=-two*y4+two*y8
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn2dy=top/bot
!
 xcu=-two*y1-two*y7+d4*y8
 xsq=two*y1+y2-two*y4+two*y5-y6-two*y8
 xon=-y2+two*y4+y6-two*y8
 ecu=-two*y5+d4*y6-two*y7
 esq=-two*y1+two*y2-y4-two*y5+two*y6+y8
 eon=two*y2-y4-two*y6+y8
 x2e2=-d3*y1+d3*y2-d3*y4+d3*y5-d3*y6+d3*y8
 x1e1=-d6*y2+d6*y4-two*y6+two*y8
 x2e1=d4*y2-d5*y4+d5*y5-d4*y6+d3*y7-d3*y8
 x1e2=d5*y1-d5*y2+d4*y4-d3*y6+d3*y7-d4*y8
 cons=zero
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn3dy=top/bot
!
 xcu=d4*y1+d4*y7-d8*y8
 xsq=-two*y1-two*y2+two*y3-two*y5+two*y6+two*y7
 xon=-two*y1-two*y3-two*y5-two*y7+d8*y8
 ecu=zero
 esq=y1-two*y2+y3-y5+two*y6-y7
 eon=-y1+y3-y5+y7 
 x2e2=d3*y1-d6*y2+d3*y3-d3*y5+d6*y6-d3*y7  
 x1e1=-two*y1+d8*y2-d6*y3-d6*y5+d8*y6-two*y7 
 x2e1=d3*y1-d8*y2+d5*y3-d5*y5+d8*y6-d3*y7 
 x1e2=-d4*y1+d8*y2-d4*y3-d4*y5+d8*y6-d4*y7  
 cons=two*y2-two*y6
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn4dy=top/bot
!
 xcu=-two*y1-two*y7+d4*y8
 xsq=y2-two*y3+two*y4-y6-two*y7+two*y8
 xon=y2+two*y4-y6-two*y8
 ecu=two*y1-d4*y2+two*y3
 esq=-two*y2+two*y3+y4-two*y6+two*y7-y8
 eon=two*y2+y4-two*y6-y8
 x2e2=d3*y2-d3*y3+d3*y4-d3*y6+d3*y7-d3*y8 
 x1e1=-two*y2+d6*y4-d6*y6+two*y8
 x2e1=-d3*y1+d4*y2-d5*y3+d5*y4-d4*y6 +d3*y8
 x1e2=d3*y1-d3*y2+d4*y4-d5*y6+d5*y7-d4*y8 
 cons=zero
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn5dy=top/bot
!
 xcu=zero          
 xsq=-y1+y3-two*y4+y5-y7+two*y8
 xon=y1-y3+y5-y7 
 ecu=-d4*y1+d8*y2-d4*y3
 esq=two*y1-two*y3-two*y4+two*y5-two*y7+two*y8
 eon=two*y1-d8*y2+two*y3+two*y5+two*y7 
 x2e2=-d3*y1+d3*y3-d6*y4+d3*y5-d3*y7+d6*y8 
 x1e1=two*y1+two*y3-d8*y4+d6*y5+d6*y7-d8*y8
 x2e1=d4*y1+d4*y3-d8*y4+d4*y5+d4*y7-d8*y8
 x1e2=-d3*y1+d3*y3-d8*y4+d5*y5-d5*y7+d8*y8 
 cons=two*y4-two*y8
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn6dy=top/bot
!
 xcu=two*y3-d4*y4+two*y5         
 xsq=two*y1-y2-two*y4+two*y5+y6-two*y8
 xon=-y2+two*y4+y6-two*y8
 ecu=two*y1-d4*y2+two*y3
 esq=-two*y1+two*y2+y4-two*y5+two*y6-y8
 eon=two*y2-y4-two*y6+y8
 x2e2=d3*y1-d3*y2+d3*y4-d3*y5+d3*y6-d3*y8 
 x1e1=-two*y2+two*y4-d6*y6+d6*y8
 x2e1=-d5*y1+d4*y2-d3*y3+d3*y4-d4*y6+d5*y8
 x1e2=d3*y2-d3*y3+d4*y4-d5*y5+d5*y6-d4*y8 
 cons=zero
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn7dy=top/bot
!
 xcu=-d4*y3+d8*y4-d4*y5         
 xsq=-two*y1+two*y2+two*y3-two*y5-two*y6+two*y7 
 xon=two*y1+two*y3-d8*y4+two*y5+two*y7 
 ecu=zero      
 esq=-y1+two*y2-y3+y5-two*y6+y7 
 eon=y1-y3+y5-y7  
 x2e2=-d3*y1+d6*y2-d3*y3+d3*y5-d6*y6+d3*y7  
 x1e1=-d6*y1+d8*y2-two*y3-two*y5+d8*y6-d6*y7  
 x2e1=d5*y1-d8*y2+d3*y3-d3*y5+d8*y6-d5*y7  
 x1e2=d4*y1-d8*y2+d4*y3+d4*y5-d8*y6+d4*y7  
 cons=-two*y2+two*y6
!
 top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
 top=top+xi2et2*x2e2+xi1et1*x1e1
 top=top+xi2et1*x2e1+xi1et2*x1e2+cons
!
 dn8dy=top/bot
!
 bee=zero
 bee(1,1)=dn1dx
 bee(1,3)=dn2dx
 bee(1,5)=dn3dx
 bee(1,7)=dn4dx
 bee(1,9)=dn5dx   
 bee(1,11)=dn6dx
 bee(1,13)=dn7dx   
 bee(1,15)=dn8dx
 bee(2,2)=dn1dy
 bee(2,4)=dn2dy
 bee(2,6)=dn3dy
 bee(2,8)=dn4dy
 bee(2,10)=dn5dy
 bee(2,12)=dn6dy
 bee(2,14)=dn7dy
 bee(2,16)=dn8dy
 bee(3,1)=dn1dy
 bee(3,3)=dn2dy
 bee(3,5)=dn3dy
 bee(3,7)=dn4dy
 bee(3,9)=dn5dy
 bee(3,11)=dn6dy
 bee(3,13)=dn7dy
 bee(3,15)=dn8dy
 bee(3,2)=dn1dx
 bee(3,4)=dn2dx   
 bee(3,6)=dn3dx
 bee(3,8)=dn4dx
 bee(3,10)=dn5dx
 bee(3,12)=dn6dx   
 bee(3,14)=dn7dx
 bee(3,16)=dn8dx
 RETURN
END SUBROUTINE bee8
SUBROUTINE beemat(bee,deriv)
!
! This subroutine forms the bee matrix in 2-d (ih=3 or 4) or 3-d (ih=6).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::deriv(:,:)
 REAL(iwp),INTENT(OUT)::bee(:,:)
 INTEGER::k,l,m,n,ih,nod
 REAL::x,y,z
 bee=0.0_iwp
 ih=UBOUND(bee,1)
 nod=UBOUND(deriv,2)
 SELECT CASE (ih)
 CASE(3,4)
   DO m=1,nod
     k=2*m
     l=k-1
     x=deriv(1,m)
     y=deriv(2,m)
     bee(1,l)=x
     bee(3,k)=x
     bee(2,k)=y
     bee(3,l)=y
   END DO
 CASE(6)
   DO m=1,nod
     n=3*m
     k=n-1
     l=k-1
     x=deriv(1,m)
     y=deriv(2,m)
     z=deriv(3,m)
     bee(1,l)=x
     bee(4,k)=x
     bee(6,n)=x
     bee(2,k)=y
     bee(4,l)=y
     bee(5,n)=y
     bee(3,n)=z
     bee(5,k)=z
     bee(6,l)=z
   END DO
 CASE DEFAULT
   WRITE(*,*)'wrong dimension for nst in bee matrix'        
 END SELECT   
RETURN
END SUBROUTINE beemat

SUBROUTINE bisect(d,e,acheps,ifail)
!
! This subroutine finds the eigenvalues of a tridiagonal matrix,
! given with its diagonal elements in the array d(n) and
! its subdiagonal elements in the last n - 1 stores of the
! array e(n), using ql transformations. The eigenvalues are
! overwritten on the diagonal elements in the array d in
! ascending order. The subroutine will fail if any one
! eigenvalue takes more than 30 iterations.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::acheps
 REAL(iwp),INTENT(IN OUT)::d(0:),e(0:)
 INTEGER,INTENT(IN OUT)::ifail
 INTEGER::m,n,i,l,j,i1,m1,ii,aux
 REAL(iwp)::b,f,h,g,p,r,c,s,zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp,         &
   two=2.0_iwp
 n=UBOUND(d,1)
 IF(n/=1)THEN
   DO i=2,n
     e(i-1)=e(i)
   END DO
 END IF
 e(n)=zero
 b=zero
 f=zero
 loop_1: DO l=1,n
   j=0
   h=acheps*(ABS(d(l))+ABS(e(l)))
   IF(b<h)b=h
   ! look for small sub diagonal element
   DO  m=l,n
     IF(ABS(e(m))<=b)EXIT
   END DO
   IF(m/=l)THEN
     loop_2: DO
       IF(j==30)THEN
         ifail=1
         EXIT loop_1
       END IF
       j = j + 1
       ! form shift
       g = d(l)
       h = d(l+1) - g
       IF(ABS(h)<ABS(e(l)))THEN
         p=h*pt5/e(l)
         r=SQRT(p*p+one)
         h=p+r
         IF(p<zero)h=p-r
         d(l)=e(l)/h
       ELSE
         p=two*e(l)/h
         r=SQRT(p*p+one)
         d(l)=e(l)*p/(one+r)
       END IF
       h=g-d(l)
       i1=l+1
       IF(i1<=n)THEN
         DO i=i1,n
           d(i)=d(i)-h
         END DO
       END IF
       f=f+h
       ! ql transformation
       p=d(m)
       c=one
       s=zero
       m1=m-1
       loop_4: DO ii=l,m1
         i=m1-ii+l
         g=c*e(i)
         h=c*p
         IF(ABS(p)>=ABS(e(i)))THEN
           c= e(i)/p
           r=SQRT(c*c+one)
           e(i+1)=s*p*r
           s=c/r
           c=one/r
         ELSE
           c=p/e(i)
           r=SQRT(c*c+one)
           e(i+1)=s*e(i)*r
           s=one/r
           c=c/r
         END IF
         p=c*d(i)-s*g
         d(i+1)=h+s*(c*g+s*d(i))
       END DO loop_4
       e(l)=s*p
       d(l)=c*p
       IF(ABS(e(l))<=b)EXIT loop_2
     END DO loop_2
   END IF
   p=d(l)+f
   ! order eigenvalue
   aux=0
   IF(l/=1)THEN
     loop_3: DO ii=2,l
       i=l-ii+2
       IF(p>=d(i-1))THEN
         aux=1
         EXIT loop_3
       END IF
       d(i)=d(i-1)
     END DO loop_3
   END IF
   IF(aux==0)THEN
     i=1
   END IF
   d(i) = p
   ifail=0
 END DO loop_1
RETURN
END SUBROUTINE bisect
SUBROUTINE bmat_nonaxi(bee,radius,coord,deriv,fun,iflag,lth)
!
! This subroutine forms the strain-displacement matrix for
! axisymmetric solids subjected to non-axisymmetric loading.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::deriv(:,:),fun(:),coord(:,:)
 REAL(iwp),INTENT(OUT)::bee(:,:),radius
 INTEGER,INTENT(IN)::iflag,lth
 INTEGER::nod,k,l,m,n
 nod=UBOUND(deriv,2) 
 bee=0.0_iwp
 radius=SUM(fun*coord(:,1))
 DO m=1,nod
   n=3*m      
   k=n-1      
   l=k-1
   bee(1,l)=deriv(1,m)
   bee(2,k)=deriv(2,m)
   bee(3,l)=fun(m)/radius
   bee(3,n)=iflag*lth*bee(3,l) 
   bee(4,l)=deriv(2,m)
   bee(4,k)=deriv(1,m)
   bee(5,k)=-iflag*lth*fun(m)/radius 
   bee(5,n)=deriv(2,m)
   bee(6,l)=bee(5,k) 
   bee(6,n)=deriv(1,m)-fun(m)/radius
 END DO
RETURN
END SUBROUTINE bmat_nonaxi


SUBROUTINE checon(loads,oldlds,tol,converged)
!
! This subroutine sets converged to .FALSE. if relative change in loads
! and oldlds is greater than tol and updates oldlds.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::loads(0:),tol
 REAL(iwp),INTENT(IN OUT)::oldlds(0:)
 LOGICAL,INTENT(OUT)::converged
 CONVERGED=.TRUE.
 CONVERGED=(MAXVAL(ABS(loads-oldlds))/MAXVAL(ABS(loads))<=tol)
 oldlds=loads
RETURN
END SUBROUTINE checon
 subroutine chobk1(kb,loads)
!Choleski back-substitution
 implicit none
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),intent(in)::kb(:,:)
 REAL(iwp),intent(in out)::loads(0:)
 integer::iw,n,i,j,k,l,m
 REAL(iwp)::x
 n=size(kb,1)
 iw=size(kb,2)-1
 loads(1)=loads(1)/kb(1,iw+1)
 do i=2,n
   x=.0
   k=1
   if(i<=iw+1)k=iw-i+2
   do j=k,iw
     x=x+kb(i,j)*loads(i+j-iw-1)
   end do
   loads(i)=(loads(i)-x)/kb(i,iw+1)
 end do
return
end subroutine chobk1


 subroutine chobk2(kb,loads)
!Choleski back-substitution
 implicit none
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),intent(in)::kb(:,:)
 REAL(iwp),intent(in out)::loads(0:)
 integer::iw,n,i,j,k,l,m
 REAL(iwp)::x
 n=size(kb,1)
 iw=size(kb,2)-1
 loads(n)=loads(n)/kb(n,iw+1)
 do i=n-1,1,-1
   x=0.0
   l=i+iw
   if(i>n-iw)l=n
   m=i+1
   do j=m,l
     x=x+kb(j,iw+i-j+1)*loads(j)
   end do
   loads(i)=(loads(i)-x)/kb(i,iw+1)
 end do
 return
 end subroutine chobk2

subroutine cholin(kb)
! Choleski reduction on kb(l,iw+1) stored as a lower triangle
 implicit none
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),intent(in out)::kb(:,:)
 integer::i,j,k,l,ia,ib,n,iw
 REAL(iwp)::x
 n=ubound(kb,1)
 iw=ubound(kb,2)-1
 do i=1,n
   x=.0
   do j=1,iw
     x=x+kb(i,j)**2
   end do
   kb(i,iw+1)=sqrt(kb(i,iw+1)-x)
   do k=1,iw
     x=.0
     if(i+k<=n)then
       if(k/=iw)then
         do l=iw-k,1,-1
           x=x+kb(i+k,l)*kb(i,l+k)
         end do
       end if
       ia=i+k
       ib=iw-k+1
       kb(ia,ib)=(kb(ia,ib)-x)/kb(i,iw+1)
     end if
   end do
 end do
return
end subroutine cholin  
SUBROUTINE contour(loads,g_coord,g_num,ned,argv,nlen,ips)
!
! This subroutine produces a PostScript output file "*.con" displaying
! a contour map of loads over the finite element mesh.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN):: loads(0:),g_coord(:,:)
 INTEGER,INTENT(IN)::g_num(:,:),ned,ips,nlen
 CHARACTER(*),INTENT(IN)::argv
 REAL(iwp)::xmin,xmax,ymin,ymax,width,height,scale=72,sxy,xo,yo
 REAL(iwp)::pmin,pmax,ratio,x12,y12,x23,y23,x34,y34,x41,y41,x1,y1
 REAL(iwp)::pt5=0.5_iwp,elvn=11.0_iwp,ept5=8.5_iwp,eight=8.0_iwp,         &
   onept5=1.5_iwp,fpt5=5.5_iwp
 LOGICAL::s12,s23,s34,s41,draw
 INTEGER,ALLOCATABLE::corner(:,:)
 REAL(iwp),ALLOCATABLE::cont(:)
 INTEGER::i,j,k,l,nn,nels,nci,nod,ns,i1,i2,j1,j2
!
 OPEN(ips,FILE=argv(1:nlen)//'.con')
!
!                       compute size of mesh
!
 nn=UBOUND(g_coord,2)
 nod= UBOUND(g_num,1)
 nels=UBOUND(g_num,2)
 xmin=g_coord(1,1)
 xmax=g_coord(1,1)
 ymin=g_coord(2,1)
 ymax=g_coord(2,1)
 DO i=2,nn
   IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)      
   IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)      
   IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)      
   IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)      
 END DO
 width =xmax-xmin
 height=ymax-ymin
!
!                       allow 1.5" margin minimum on each side of figure
!
!                       portrait mode 
!
 IF(height.GE.elvn/ept5*width)THEN 
!
!                       height governs the scale
!
   sxy=scale*eight/height
   xo=scale*pt5*(ept5-eight*width/height)
   yo=scale*onept5
 ELSE
!
!                       width governs the scale
!
   sxy=scale*fpt5/width
   xo=scale*onept5
   yo=scale*pt5*(elvn-fpt5*height/width)
 END IF
!
!                       find range of potentials and contour values
!
 nci=ned+1
 pmin=MINVAL(loads(1:))
 pmax=MAXVAL(loads(1:))
 ALLOCATE(cont(nci))
 cont(1)=pmin
 DO i=2,nci
   cont(i)=cont(i-1)+(pmax-pmin)/ned
 END DO
!
!                       start PostScript output
!
 WRITE(ips,'(a)')'%!PS-Adobe-1.0'
 WRITE(ips,'(a)')'%%DocumentFonts: none'
 WRITE(ips,'(a)')'%%Pages: 1'
 WRITE(ips,'(a)')'%%EndComments'
 WRITE(ips,'(a)')'/m {moveto} def'
 WRITE(ips,'(a)')'/l {lineto} def'
 WRITE(ips,'(a)')'/s {stroke} def'
 WRITE(ips,'(a)')'%%EndProlog'
 WRITE(ips,'(a)')'%%Page: 0 1'
 WRITE(ips,'(a)')'gsave'
!
 WRITE(ips,'(2f9.2,a)') xo, yo, ' translate'
 WRITE(ips,'(f9.2,a)') 1.0, ' setlinewidth'
!
!                       draw the mesh outline
!
 IF(nod==3.OR.nod==6.OR.nod==10.OR.nod==15)ns=3
 IF(nod==4.OR.nod==8.OR.nod==9)ns=4
 ALLOCATE(corner(ns,2))
 IF(nod== 3)corner=RESHAPE((/1,2,3,2,3,1/),(/3,2/))
 IF(nod== 6)corner=RESHAPE((/1,3,5,3,5,1/),(/3,2/))
 IF(nod==10)corner=RESHAPE((/1,4,7,4,7,1/),(/3,2/))
 IF(nod==15)corner=RESHAPE((/1,5,9,5,9,1/),(/3,2/))
 IF(nod== 4)corner=RESHAPE((/1,2,3,4,2,3,4,1/),(/4,2/))
 IF(nod== 8)corner=RESHAPE((/1,3,5,7,3,5,7,1/),(/4,2/))
 IF(nod== 9)corner=RESHAPE((/1,3,5,7,3,5,7,1/),(/4,2/))
 DO i=1,nels
   DO j=1,ns
     draw=.TRUE.
     i1=g_num(corner(j,1),i)
     i2=g_num(corner(j,2),i)
     DO k=1,nels
       DO l=1,ns
         j1=g_num(corner(l,1),k)
         j2=g_num(corner(l,2),k)
         IF((i1==j2).AND.(i2==j1))THEN
           draw=.FALSE.
           EXIT
         END IF
       END DO
       IF(.NOT.draw)EXIT
     END DO
     IF(draw)THEN
       x1=sxy*(g_coord(1,i1)-xmin)
       y1=sxy*(g_coord(2,i1)-ymin)
       WRITE(ips,'(2f9.2,a)') x1, y1,' m'
       x1=sxy*(g_coord(1,i2)-xmin)
       y1=sxy*(g_coord(2,i2)-ymin)
       WRITE(ips,'(2f9.2,a)') x1, y1,' l'
       WRITE(ips,'(a)')' s'
     END IF
   END DO
 END DO
!
!                       check intersection of contours with each element
!
 WRITE(ips,'(f9.2,a)') 0.5, ' setlinewidth'
 DO i=2,nci-1
   DO j=1,nels
     s12=.FALSE.
     s23=.FALSE.
     s34=.FALSE.
     s41=.FALSE.
     IF((loads(g_num(1,j))<=cont(i).AND.loads(g_num(2,j))>cont(i)).OR.    &
       (loads(g_num(2,j))<=cont(i).AND.loads(g_num(1,j))>cont(i)))        &
       s12=.TRUE.
     IF((loads(g_num(2,j))<=cont(i).AND.loads(g_num(3,j))>cont(i)).OR.    &
       (loads(g_num(3,j))<=cont(i).AND.loads(g_num(2,j))>cont(i)))        &
       s23=.TRUE.
     IF((loads(g_num(3,j))<=cont(i).AND.loads(g_num(4,j))>cont(i)).OR.    &
       (loads(g_num(4,j))<=cont(i).AND.loads(g_num(3,j))>cont(i)))        &
       s34=.TRUE.
     IF((loads(g_num(4,j))<=cont(i).AND.loads(g_num(1,j))>cont(i)).OR.    &
       (loads(g_num(1,j))<=cont(i).AND.loads(g_num(4,j))>cont(i)))        &
       s41=.TRUE.
     IF(s12)THEN
       ratio=(cont(i)-loads(g_num(1,j)))/                                 &
       (loads(g_num(2,j))-loads(g_num(1,j)))
       x12=sxy*(g_coord(1,g_num(1,j))+                                    &
         ratio*(g_coord(1,g_num(2,j))-g_coord(1,g_num(1,j)))-xmin)
       y12=sxy*(g_coord(2,g_num(1,j))+                                    &
         ratio*(g_coord(2,g_num(2,j))-g_coord(2,g_num(1,j)))-ymin)
     END IF
     IF(s23)THEN
       ratio=(cont(i)-loads(g_num(2,j)))/                                 &
         (loads(g_num(3,j))-loads(g_num(2,j)))
       x23=sxy*(g_coord(1,g_num(2,j))+                                    &
         ratio*(g_coord(1,g_num(3,j))-g_coord(1,g_num(2,j)))-xmin)
       y23=sxy*(g_coord(2,g_num(2,j))+                                    &
         ratio*(g_coord(2,g_num(3,j))-g_coord(2,g_num(2,j)))-ymin)
     END IF
     IF(s34)THEN
       ratio=(cont(i)-loads(g_num(3,j)))/                                 &
         (loads(g_num(4,j))-loads(g_num(3,j)))
       x34=sxy*(g_coord(1,g_num(3,j))+                                    &
         ratio*(g_coord(1,g_num(4,j))-g_coord(1,g_num(3,j)))-xmin)
       y34=sxy*(g_coord(2,g_num(3,j))+                                    &
         ratio*(g_coord(2,g_num(4,j))-g_coord(2,g_num(3,j)))-ymin)
     END IF
     IF(s41)THEN
       ratio=(cont(i)-loads(g_num(4,j)))/                                 &
         (loads(g_num(1,j))-loads(g_num(4,j)))
       x41=sxy*(g_coord(1,g_num(4,j))+                                    &
         ratio*(g_coord(1,g_num(1,j))-g_coord(1,g_num(4,j)))-xmin)
       y41=sxy*(g_coord(2,g_num(4,j))+                                    &
         ratio*(g_coord(2,g_num(1,j))-g_coord(2,g_num(4,j)))-ymin)
     END IF
!
!                       draw contours 
!
     IF(s12)THEN
       WRITE(ips,'(2f9.2,a)')x12,y12,' m'
       IF(s23)WRITE(ips,'(2f9.2,a)')x23,y23,' l s'
       IF(s34)WRITE(ips,'(2f9.2,a)')x34,y34,' l s'
       IF(s41)WRITE(ips,'(2f9.2,a)')x41,y41,' l s'
     END IF
     IF(s23)THEN
       WRITE(ips,'(2f9.2,a)')x23,y23,' m'
       IF(s34)WRITE(ips,'(2f9.2,a)')x34,y34,' l s'
       IF(s41)WRITE(ips,'(2f9.2,a)')x41,y41,' l s'
     END IF
     IF(s34.AND.s41)THEN
       WRITE(ips,'(2f9.2,a)')x34,y34,' m'
       WRITE(ips,'(2f9.2,a)')x41,y41,' l s'
     END IF
   END DO
 END DO
!                       close output file
 WRITE(ips,'(a)')'grestore'
 WRITE(ips,'(a)')'showpage'
 CLOSE(ips)
!
RETURN
END SUBROUTINE contour

SUBROUTINE cross_product(b,c,a)
!
! This subroutine forms the cross product of two vectors, a = b x c
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::b(:),c(:)
 REAL(iwp),INTENT(OUT)::a(:,:)
 INTEGER::ib,ic,i,j
 ib=SIZE(b)
 ic=SIZE(c)
 DO i=1,ib
   DO j=1,ic
     a(i,j)=b(i)*c(j)
   END DO
 END DO
RETURN
END SUBROUTINE cross_product   
SUBROUTINE deemat(dee,e,v)
!
! This subroutine returns the elastic dee matrix for ih=3 (plane strain),
! ih=4 (axisymmetry or plane strain elastoplasticity) or ih=6
! (three dimensions).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::e,v
 REAL(iwp),INTENT(OUT)::dee(:,:)
 REAL(iwp)::v1,v2,c,vv,zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp,two=2.0_iwp
 INTEGER::i,ih
 dee=zero  
 ih=UBOUND(dee,1)
 v1=one-v
 c=e/((one+v)*(one-two*v))
 SELECT CASE(ih)
 CASE(3)
   dee(1,1)=v1*c
   dee(2,2)=v1*c
   dee(1,2)=v*c
   dee(2,1)=v*c
   dee(3,3)=pt5*c*(one-two*v)
 CASE(4)
   dee(1,1)=v1*c
   dee(2,2)=v1*c
   dee(4,4)=v1*c
   dee(3,3)=pt5*c*(one-two*v) 
   dee(1,2)=v*c
   dee(2,1)=v*c
   dee(1,4)=v*c
   dee(4,1)=v*c
   dee(2,4)=v*c
   dee(4,2)=v*c
 CASE(6)
   v2=v/(one-v)
   vv=(one-two*v)/(one-v)*pt5
   DO i=1,3
     dee(i,i)=one
   END DO
   DO i=4,6
     dee(i,i)=vv
   END DO
   dee(1,2)=v2
   dee(2,1)=v2
   dee(1,3)=v2
   dee(3,1)=v2
   dee(2,3)=v2
   dee(3,2)=v2
   dee=dee*e/(two*(one+v)*vv)
 CASE DEFAULT
   WRITE(*,*)'wrong size for dee matrix'
 END SELECT
RETURN
END SUBROUTINE deemat    
FUNCTION determinant(jac)RESULT(det)
!
! This function returns the determinant of a 1x1, 2x2 or 3x3
! Jacobian matrix.
!
 IMPLICIT NONE    
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::jac(:,:)
 REAL(iwp)::det
 INTEGER::it 
 it=UBOUND(jac,1)  
 SELECT CASE(it)
 CASE(1)
   det=1.0_iwp
 CASE(2)
   det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
 CASE(3)
   det=jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3))
   det=det-jac(1,2)*(jac(2,1)*jac(3,3)-jac(3,1)*jac(2,3))
   det=det+jac(1,3)*(jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2))
 CASE DEFAULT
   WRITE(*,*)' wrong dimension for Jacobian matrix'
 END SELECT
RETURN
END FUNCTION determinant
SUBROUTINE dismsh(loads,nf,ratmax,g_coord,g_num,argv,nlen,ips)
!
! This subroutine produces a PostScript output file "*.dis" displaying
! the deformed finite element mesh.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::g_coord(:,:),loads(0:),ratmax
 INTEGER,INTENT(IN)::g_num(:,:),ips,nf(:,:),nlen
 CHARACTER(*),INTENT(IN)::argv
 REAL(iwp)::width,height,scale=72,sxy,xo,yo,x,y,dismag,vmax
 REAL(iwp)::xmin,xmax,ymin,ymax,dmax,zero=0.0_iwp,pt5=0.5_iwp,            &
   opt5=1.5_iwp,fpt5=5.5_iwp,d8=8.0_iwp,ept5=8.5_iwp,d11=11.0_iwp
 INTEGER::i,ii,j,jj,nn,nel,nod
 OPEN(ips,FILE=argv(1:nlen)//'.dis')
!
 nn=UBOUND(nf,2)
 xmin=g_coord(1,1)
 xmax=g_coord(1,1)
 ymin=g_coord(2,1)
 ymax=g_coord(2,1)
 DO i=2,nn
   IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)      
   IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)      
   IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)      
   IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)      
 END DO
 width=xmax-xmin
 height=ymax-ymin
 dmax=ratmax*width
 IF(height>width)dmax=ratmax*height
!
 vmax=zero
 DO i=1,nn
   DO j=1,2
     IF(ABS(loads(nf(j,i)))>vmax)vmax=ABS(loads(nf(j,i)))
   END DO
 END DO
 dismag=dmax/vmax
!
 xmin=g_coord(1,1)
 xmax=g_coord(1,1)
 ymin=g_coord(2,1)
 ymax=g_coord(2,1)
 DO i=1,nn
   IF(g_coord(1,i)+dismag*loads(nf(1,i))<xmin)                            &
     xmin=g_coord(1,i)+dismag*loads(nf(1,i))
   IF(g_coord(1,i)+dismag*loads(nf(1,i))>xmax)                            &
     xmax=g_coord(1,i)+dismag*loads(nf(1,i))      
   IF(g_coord(2,i)+dismag*loads(nf(2,i))<ymin)                            &
     ymin=g_coord(2,i)+dismag*loads(nf(2,i))      
   IF(g_coord(2,i)+dismag*loads(nf(2,i))>ymax)                            &
     ymax=g_coord(2,i)+dismag*loads(nf(2,i))      
!
   IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)
   IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)
   IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)
   IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)
 END DO
!
 width =xmax-xmin
 height=ymax-ymin
!
!                       allow 1.5" margin minimum on each side of figure
!
!                       portrait mode 
!
 IF(height.GE.d11/ept5*width)THEN
!
!                       height governs the scale
!
   sxy=scale*d8/height
   xo=scale*pt5*(ept5-d8*width/height)
   yo=scale*opt5
 ELSE
!
!                       width governs the scale
!
   sxy=scale*fpt5/width
   xo=scale*opt5
   yo=scale*pt5*(d11-fpt5*height/width)
 END IF
!
 nel=UBOUND(g_num,2)
 nod=UBOUND(g_num,1)
!
!                       start PostScript output
!
 WRITE(ips,'(a)')'%!PS-Adobe-1.0'
 WRITE(ips,'(a)')'%%DocumentFonts: none'
 WRITE(ips,'(a)')'%%Pages: 1'
 WRITE(ips,'(a)')'%%EndComments'
 WRITE(ips,'(a)')'/m {moveto} def'
 WRITE(ips,'(a)')'/l {lineto} def'
 WRITE(ips,'(a)')'/s {stroke} def'
 WRITE(ips,'(a)')'/c {closepath} def'
 WRITE(ips,'(a)')'%%EndProlog'
 WRITE(ips,'(a)')'%%Page: 0 1'
 WRITE(ips,'(a)')'gsave'
!
!                       draw the deformed mesh
!
 WRITE(ips,'(2f9.2,a)') xo, yo, ' translate'
 WRITE(ips,'(f9.2,a)') 0.5, ' setlinewidth'
 IF(nod==5)nod=4
 IF(nod==9)nod=8
 IF(nod==10)nod=9
 IF(nod==15)nod=12
 DO i=1,nel
   ii=g_num(1,i)
   IF(ii==0)CYCLE
   x=sxy*(g_coord(1,ii)+dismag*loads(nf(1,ii))-xmin)
   y=sxy*(g_coord(2,ii)+dismag*loads(nf(2,ii))-ymin)
   WRITE(ips,'(2f9.2,a)') x, y,' m'
   DO j=2,nod
     jj=g_num(j,i)
     x=sxy*(g_coord(1,jj)+dismag*loads(nf(1,jj))-xmin)
     y=sxy*(g_coord(2,jj)+dismag*loads(nf(2,jj))-ymin)
     WRITE(ips,'(2f9.2,a)') x, y,' l'
   END DO
   WRITE(ips,'(a)')'c s'
 END DO
!
 WRITE(ips,'(a)')'grestore'
 WRITE(ips,'(a)')'showpage'
 CLOSE(ips)
!
RETURN
END SUBROUTINE dismsh


























SUBROUTINE ecmat(ecm,fun,ndof,nodof)
!
! This subroutine forms the element consistent mass matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::fun(:)
 REAL(iwp),INTENT(OUT)::ecm(:,:)
 INTEGER,INTENT(IN)::nodof,ndof
 INTEGER::nod,i,j
 REAL::nt(ndof,nodof),tn(nodof,ndof),zero=0.0_iwp
 nod=ndof/nodof
 nt=zero
 tn=zero
 DO i=1,nod 
   DO j=1,nodof
     nt((i-1)*nodof+j,j)=fun(i)
     tn(j,(i-1)*nodof+j)=fun(i)
   END DO
 END DO
 ecm=MATMUL(nt,tn)
RETURN
END SUBROUTINE ecmat
SUBROUTINE elmat(area,rho,emm)
!
! This subroutine forms the "analytical" lumped mass matrix for
! quadrilateral 4- or 8-node plane strain elements.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::area,rho  
 REAL(iwp),INTENT(OUT)::emm(:,:)
 REAL(iwp)::zero=0.0_iwp,pt2=0.2_iwp,pt25=0.25_iwp
 INTEGER::i,ndof
 ndof=UBOUND(emm,1)
 emm=zero
 SELECT CASE(ndof)
 CASE(8) 
   DO i=1,8
     emm(i,i)=pt25*area*rho
   END DO
 CASE(16)
   DO i=1,16
     emm(i,i)=pt2*area*rho
   END DO
   DO i=1,13,4
     emm(i,i)=pt25*emm(3,3)
   END DO
   DO i=2,14,4
     emm(i,i)=pt25*emm(3,3)
   END DO                     
 CASE DEFAULT
   WRITE(*,*)"Wrong number of nodes for rectangular element"
 END SELECT
RETURN
END SUBROUTINE elmat   

SUBROUTINE exc_nods(noexe,exele,g_num,totex,ntote,nf)
!
! This subroutine forms the nodes removed in an excavation lift.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN)::noexe,exele(:),g_num(:,:)
 INTEGER,INTENT(IN OUT)::totex(:),ntote,nf(:,:)
 INTEGER::i,jj,k,iel,modex,nodex,ncheck,nels,nod
 nels=UBOUND(g_num,2)
 nod=UBOUND(g_num,1)
 ntote=ntote+noexe
 totex(ntote-noexe+1:ntote)=exele
 DO i=1,noexe
   DO k=1,8
     nodex=0
     ncheck=g_num(k,exele(i))
     DO iel=1,nels
       modex=0
       DO jj=1,ntote
         IF(iel==totex(jj))THEN
           modex=1
           EXIT
         END IF
       END DO
       IF(modex==1)CYCLE
       DO jj=1,nod
         IF(ncheck==g_num(jj,iel))THEN
           nodex=1
           EXIT
         END IF
       END DO
       IF(nodex==1)EXIT
     END DO
     IF(nodex==0)nf(:,ncheck)=0
   END DO
 END DO
RETURN
END SUBROUTINE exc_nods
SUBROUTINE fkdiag(kdiag,g)
!
! This subroutine computes the skyline profile.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN)::g(:)
 INTEGER,INTENT(OUT)::kdiag(:)
 INTEGER::idof,i,iwp1,j,im,k
 idof=SIZE(g)
 DO i=1,idof
   iwp1=1
   IF(g(i)/=0)THEN
     DO j=1,idof
       IF(g(j)/=0)THEN
         im=g(i)-g(j)+1
         IF(im>iwp1)iwp1=im
       END IF
     END DO
     k=g(i)
     IF(iwp1>kdiag(k))kdiag(k)=iwp1
   END IF
 END DO
RETURN
END SUBROUTINE fkdiag
SUBROUTINE fmacat(vmfl,acat)
!
! This subroutine sets up an intermediate matrix acat.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::vmfl(:)
 REAL(iwp),INTENT(OUT)::acat(:,:)
 REAL(iwp)::temp(4,4),zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp,three=3.0_iwp
 INTEGER::i,j
 temp=zero
 temp(1,1)=one
 temp(1,2)=-pt5
 temp(1,4)=-pt5
 temp(2,1)=-pt5
 temp(2,2)=one
 temp(2,4)=-pt5
 temp(3,3)=three
 temp(4,1)=-pt5
 temp(4,2)=-pt5
 temp(4,4)=one
 DO i=1,4
   DO j=1,4
     acat(i,j)=vmfl(i)*vmfl(j)
   END DO
 END DO
 acat=temp-acat
RETURN 
END SUBROUTINE fmacat
SUBROUTINE fmdsig(dee,e,v)
!
! This subroutine returns the elastic dee matrix for ih=3 (plane stress),
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::e,v
 REAL(iwp),INTENT(OUT)::dee(:,:)
 REAL(iwp)::zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp
 INTEGER::ih
 dee=zero  
 ih=UBOUND(dee,1)
 SELECT CASE(ih)
 CASE(3)
   dee=zero
   dee(1,1)=e/(one-v*v)
   dee(2,2)=dee(1,1)
   dee(3,3)=pt5*e/(one+v)
   dee(1,2)=v*dee(1,1)
   dee(2,1)=dee(1,2)
 CASE DEFAULT
   WRITE(*,*)'wrong size for dee matrix'
 END SELECT
RETURN
END SUBROUTINE fmdsig   
SUBROUTINE fmkdke(km,kp,c,ke,kd,theta)
!
! This subroutine builds up the 'coupled' stiffnesses ke and kd from 
! the 'elastic' stiffness km, fluid stiffness kp and coupling matrix c.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::km(:,:),kp(:,:),c(:,:),theta
 REAL(iwp),INTENT(OUT)::ke(:,:),kd(:,:)
 INTEGER::ndof
 REAL::one=1.0_iwp
 ndof=SIZE(km,1)
 ke(:ndof,:ndof)=theta*km
 ke(:ndof,ndof+1:)=theta*c
 ke(ndof+1:,:ndof)=theta*TRANSPOSE(c)
 ke(ndof+1:,ndof+1:)=-theta**2*kp
 kd(:ndof,:ndof)=(theta-one)*km
 kd(:ndof,ndof+1:)=(theta-one)*c
 kd(ndof+1:,:ndof)=ke(ndof+1:,:ndof)
 kd(ndof+1:,ndof+1:)=theta*(one-theta)*kp
RETURN
END SUBROUTINE fmkdke

SUBROUTINE fmplat(d2x,d2y,d2xy,points,aa,bb,i)
!
! This subroutine forms the 2nd derivatives for rectangular
! plate bending elements.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::points(:,:),aa,bb
 REAL(iwp),INTENT(OUT)::d2x(:),d2y(:),d2xy(:)
 INTEGER,INTENT(IN)::i
 REAL(iwp)::x,e,xp1,xp12,xp13,ep1,ep12,ep13,p1,q1,p2,q2,p3,q3,p4,q4,dp1,  &
   dq1,dp2,dq2,dp3,dq3,dp4,dq4,d2p1,d2p2,d2p3,d2p4,d2q1,d2q2,d2q3,d2q4,   &
   pt25=0.25_iwp,pt375=0.375_iwp,pt5=0.5_iwp,pt75=0.75_iwp,one=1.0_iwp,   &
   opt5=1.5_iwp,d3=3.0_iwp
 x=points(i,1)
 e=points(i,2)
 xp1=x+one
 xp12=xp1*xp1
 xp13=xp12*xp1
 ep1=e+one
 ep12=ep1*ep1
 ep13=ep12*ep1
 p1=one-pt75*xp12+pt25*xp13
 q1=one-pt75*ep12+pt25*ep13
 p2=pt5*aa*xp1*(one-xp1+pt25*xp12)
 q2=pt5*bb*ep1*(one-ep1+pt25*ep12)
 p3=pt25*xp12*(d3-xp1)
 q3=pt25*ep12*(d3-ep1)
 p4=pt25*aa*xp12*(pt5*xp1-one)
 q4=pt25*bb*ep12*(pt5*ep1-one)
 dp1=opt5*xp1*(pt5*xp1-one)
 dq1=opt5*ep1*(pt5*ep1-one)
 dp2=aa*(pt5-xp1+pt375*xp12)
 dq2=bb*(pt5-ep1+pt375*ep12)
 dp3=opt5*xp1*(one-pt5*xp1)
 dq3=opt5*ep1*(one-pt5*ep1)
 dp4=pt5*aa*xp1*(pt75*xp1-one)
 dq4=pt5*bb*ep1*(pt75*ep1-one)
 d2p1=opt5*x
 d2p2=pt25*aa*(d3*x-one)
 d2p3=-d2p1
 d2p4=pt25*aa*(d3*x+one)
 d2q1=opt5*e
 d2q2=pt25*bb*(d3*e-one)
 d2q3=-d2q1
 d2q4=pt25*bb*(d3*e+one)
 d2x(1)=d2p1*q1
 d2x(2)=d2p2*q1
 d2x(3)=d2p1*q2
 d2x(4)=d2p2*q2
 d2x(5)=d2p1*q3
 d2x(6)=d2p2*q3
 d2x(7)=d2p1*q4
 d2x(8)=d2p2*q4
 d2x(9)=d2p3*q3
 d2x(10)=d2p4*q3
 d2x(11)=d2p3*q4
 d2x(12)=d2p4*q4
 d2x(13)=d2p3*q1
 d2x(14)=d2p4*q1
 d2x(15)=d2p3*q2
 d2x(16)=d2p4*q2
 d2y(1)=p1*d2q1
 d2y(2)=p2*d2q1
 d2y(3)=p1*d2q2
 d2y(4)=p2*d2q2
 d2y(5)=p1*d2q3
 d2y(6)=p2*d2q3
 d2y(7)=p1*d2q4
 d2y(8)=p2*d2q4
 d2y(9)=p3*d2q3
 d2y(10)=p4*d2q3
 d2y(11)=p3*d2q4
 d2y(12)=p4*d2q4
 d2y(13)=p3*d2q1
 d2y(14)=p4*d2q1
 d2y(15)=p3*d2q2
 d2y(16)=p4*d2q2
 d2xy(1)=dp1*dq1
 d2xy(2)=dp2*dq1
 d2xy(3)=dp1*dq2
 d2xy(4)=dp2*dq2
 d2xy(5)=dp1*dq3
 d2xy(6)=dp2*dq3
 d2xy(7)=dp1*dq4
 d2xy(8)=dp2*dq4
 d2xy(9)=dp3*dq3
 d2xy(10)=dp4*dq3
 d2xy(11)=dp3*dq4
 d2xy(12)=dp4*dq4
 d2xy(13)=dp3*dq1
 d2xy(14)=dp4*dq1
 d2xy(15)=dp3*dq2
 d2xy(16)=dp4*dq2
RETURN
END SUBROUTINE fmplat
SUBROUTINE fmrmat(vmfl,dsbar,dlam,dee,rmat)
!
! This subroutine forms the rmat matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::vmfl(:),dsbar,dlam,dee(:,:)
 REAL(iwp),INTENT(OUT)::rmat(:,:)
 REAL(iwp)::acat(4,4),acatc(4,4),qmat(4,4),temp(4,4),con,zero=0.0_iwp,    &
   pt5=0.5_iwp,one=1.0_iwp,three=3.0_iwp
 INTEGER::i,j  
 temp=zero
 temp(1,1)=one
 temp(1,2)=-pt5
 temp(1,4)=-pt5
 temp(2,1)=-pt5
 temp(2,2)=one
 temp(2,4)=-pt5
 temp(3,3)=three
 temp(4,1)=-pt5
 temp(4,2)=-pt5
 temp(4,4)=one
 DO i=1,4
   DO j=1,4
     acat(i,j)=vmfl(i)*vmfl(j)
   END DO
 END DO
 acat=(temp-acat)/dsbar
 acatc=matmul(dee,acat)
 qmat=acatc*dlam
 DO i=1,4
   qmat(i,i)=qmat(i,i)+one
 END DO
 DO i=1,4
   con=qmat(i,i)
   qmat(i,i)=one
   qmat(i,:)=qmat(i,:)/con
   DO j=1,4
     IF(j/=i)THEN
       con=qmat(j,i)
       qmat(j,i)=zero
       qmat(j,:)=qmat(j,:)-qmat(i,:)*con
     END IF
   END DO
 END DO
 rmat=matmul(qmat,dee)
RETURN
END SUBROUTINE fmrmat
SUBROUTINE formaa(vmfl,rmat,daatd)
!
! This subroutine modifies the dee matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::vmfl(:),rmat(:,:)
 REAL(iwp),INTENT(OUT)::daatd(:,:)
 REAL(iwp)::flowt(1,4),flowa(4,1)   
 flowt(1,:)=vmfl 
 flowa(:,1)=vmfl
 daatd=MATMUL(MATMUL(MATMUL(rmat,flowa),flowt),rmat)
RETURN
END SUBROUTINE formaa
subroutine formkb(kb,km,g)
! lower triangular global stiffness kb stored as kb(n,iw+1)
 implicit none
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),intent(in)::km(:,:)
 REAL(iwp),intent(out)::kb(:,:)
 integer,intent(in)::g(:)
 integer::iw,idof,i,j,icd
 idof=size(km,1)
 iw=size(kb,2)-1
 do i=1,idof
   if(g(i)>0)then 
     do j=1,idof
       if(g(j)>0)then
         icd=g(j)-g(i)+iw+1
         if(icd-iw-1<=0)kb(g(i),icd)=kb(g(i),icd)+km(i,j)
       end if
     end do
   end if
 end do
return
end subroutine formkb
SUBROUTINE formke(km,kp,c,ke,theta)
!
! This subroutine creates the ke matrix for incremental Biot.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::km(:,:),kp(:,:),c(:,:),theta
 REAL(iwp),INTENT(OUT)::ke(:,:)
 INTEGER::ndof
 ndof=UBOUND(km,1)
 ke(:ndof,:ndof)=km 
 ke(:ndof,ndof+1:)=c
 ke(ndof+1:,:ndof)=TRANSPOSE(c) 
 ke(ndof+1:,ndof+1:)=-theta*kp
RETURN  
END SUBROUTINE formke
SUBROUTINE formku(ku,km,g)
!
! This subroutine assembles element matrices into symmetrical
! global matrix (stored as an upper rectangle).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::km(:,:)
 REAL(iwp),INTENT(OUT)::ku(:,:)
 INTEGER,INTENT(IN)::g(:) 
 INTEGER::i,j,icd,ndof
 ndof=UBOUND(km,1)
 DO i=1,ndof
   IF(g(i)/=0)THEN
     DO j=1,ndof
       IF(g(j)/=0)THEN
         icd=g(j)-g(i)+1
         IF(icd>=1)ku(g(i),icd)=ku(g(i),icd)+km(i,j)
       END IF
     END DO
   END IF
 END DO
RETURN
END SUBROUTINE formku
SUBROUTINE formlump(diag,emm,g)
!
! This subroutine forms the lumped global mass matrix as a vector diag.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::emm(:,:)  
 REAL(iwp),INTENT(OUT)::diag(0:)
 INTEGER,INTENT(IN)::g(:)
 INTEGER::i,ndof
 ndof=UBOUND(emm,1)
 DO i=1,ndof
   diag(g(i))=diag(g(i))+emm(i,i) 
 END DO
RETURN
END SUBROUTINE formlump   

SUBROUTINE formm(stress,m1,m2,m3)
!
! This subroutine forms the derivatives of the invariants with respect to
! stress in 2- or 3-d.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::stress(:)
 REAL(iwp),INTENT(OUT)::m1(:,:),m2(:,:),m3(:,:)
 REAL(iwp)::sx,sy,txy,tyz,tzx,sz,dx,dy,dz,sigm,zero=0.0_iwp,one=1.0_iwp,  &
   two=2.0_iwp,three=3.0_iwp,six=6.0_iwp,nine=9.0_iwp
 INTEGER::nst,i,j
 nst=UBOUND(stress,1)
 SELECT CASE(nst)
 CASE(4)
   sx=stress(1)
   sy=stress(2)
   txy=stress(3)
   sz=stress(4)
   dx=(two*sx-sy-sz)/three
   dy=(two*sy-sz-sx)/three
   dz=(two*sz-sx-sy)/three
   sigm=(sx+sy+sz)/three
   m1=zero
   m2=zero
   m3=zero
   m1(1,1:2)=one
   m1(2,1:2)=one
   m1(4,1:2)=one
   m1(1,4)=one
   m1(4,4)=one
   m1(2,4)=one
   m1=m1/nine/sigm
   m2(1,1)=two/three
   m2(2,2)=two/three
   m2(4,4)= two/three
   m2(2,4)=-one/three
   m2(4,2)=-one/three
   m2(1,2)=-one/three
   m2(2,1)=-one/three
   m2(1,4)=-one/three
   m2(4,1)=-one/three
   m2(3,3)=two
   m3(3,3)=-dz
   m3(1:2,3)=txy/three
   m3(3,1:2)=txy/three
   m3(3,4)=-two*txy/three
   m3(4,3)=-two*txy/three
   m3(1,1)=dx/three
   m3(2,4)=dx/three
   m3(4,2)=dx/three
   m3(2,2)=dy/three
   m3(1,4)=dy/three
   m3(4,1)=dy/three
   m3(4,4)=dz/three
   m3(1,2)=dz/three
   m3(2,1)=dz/three
 CASE(6)
   sx=stress(1)
   sy=stress(2)    
   sz=stress(3)
   txy=stress(4)  
   tyz=stress(5) 
   tzx=stress(6)
   sigm=(sx+sy+sz)/three
   dx=sx-sigm  
   dy=sy-sigm 
   dz=sz-sigm
   m1=zero
   m2=zero
   m1(1:3,1:3)=one/(three*sigm)
   DO i=1,3 
     m2(i,i)=two 
     m2(i+3,i+3)=six 
   END DO
   m2(1,2)=-one
   m2(1,3)=-one 
   m2(2,3)=-one
   m3(1,1)=dx
   m3(1,2)=dz 
   m3(1,3)=dy 
   m3(1,4)=txy  
   m3(1,5)=-two*tyz
   m3(1,6)=tzx 
   m3(2,2)=dy 
   m3(2,3)=dx 
   m3(2,4)=txy
   m3(2,5)=tyz 
   m3(2,6)=-two*tzx 
   m3(3,3)=dz
   m3(3,4)=-two*txy
   m3(3,5)=tyz 
   m3(3,6)=tzx
   m3(4,4)=-three*dz 
   m3(4,5)=three*tzx
   m3(4,6)=three*tyz
   m3(5,5)=-three*dx
   m3(5,6)=three*txy 
   m3(6,6)=-three*dy
   DO i=1,6 
     DO j=i+1,6
       m1(j,i)=m1(i,j) 
       m2(j,i)=m2(i,j)   
       m3(j,i)=m3(i,j)
     END DO
   END DO
   m1=m1/three
   m2=m2/three
   m3=m3/three
 CASE DEFAULT
   WRITE(*,*)"nst size not recognised in formm"
 END SELECT
RETURN   
END SUBROUTINE formm

SUBROUTINE formnf(nf)
!
! This subroutine forms the nf matrix.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN OUT)::nf(:,:)
 INTEGER::i,j,m
 m=0
 DO j=1,UBOUND(nf,2)
   DO i=1,UBOUND(nf,1)
     IF(nf(i,j)/=0)THEN
       m=m+1
       nf(i,j)=m
     END IF
   END DO
 END DO
RETURN
END SUBROUTINE formnf 
SUBROUTINE formtb(pb,km,g)
!
! This subroutine assembles an unsymmetrical band matrix pb from
! element constituent matrices km.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::km(:,:)
 INTEGER,INTENT(IN)::g(:)
 REAL(iwp),INTENT(OUT)::pb(:,:)
 INTEGER::i,j,idof,icd,iw
 idof=SIZE(km,1)
 iw=(SIZE(pb,2)-1)/2
 DO i=1,idof
   IF(g(i)/=0)THEN
     DO j=1,idof
       IF(g(j)/=0)THEN
         icd=g(j)-g(i)+iw+1
         pb(g(i),icd)=pb(g(i),icd)+km(i,j)
       END IF
     END DO
   END IF
 END DO
RETURN
END SUBROUTINE formtb
SUBROUTINE formupv(ke,c11,c12,c21,c23,c32)
!
! This subroutine forms the unsymmetrical stiffness matrix
! for the u-p-v version of the Navier-Stokes equations.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::c11(:,:),c21(:,:),c23(:,:),c32(:,:),c12(:,:)
 REAL(iwp),INTENT(OUT)::ke(:,:)  
 INTEGER::nod,nodf,ntot
 nod=UBOUND(c11,1)
 nodf=UBOUND(c21,1) 
 ntot=nod+nodf+nod
 ke(1:nod,1:nod)=c11  
 ke(1:nod,nod+1:nod+nodf)=c12
 ke(nod+1:nod+nodf,1:nod)=c21  
 ke(nod+1:nod+nodf,nod+nodf+1:ntot)=c23
 ke(nod+nodf+1:ntot,nod+1:nod+nodf)=c32 
 ke(nod+nodf+1:ntot,nod+nodf+1:ntot)=c11
RETURN
END SUBROUTINE formupv 
SUBROUTINE form_s(gg,ell,kappa,omega,gamma,s)
!
! This subroutine forms the s vector in bicgstab(l)
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::gg(:,:),kappa
 INTEGER,INTENT(IN)::ell
 REAL(iwp),INTENT(OUT)::omega,gamma(:),s(:)
 REAL(iwp)::HH(ell-1,ell-1),gamma0(ell+1),p(ell-1),q(ell-1),gamma1(ell+1),&
   ngamma0,ngamma1,cosine,zero=0.0_iwp,one=1.0_iwp
 hh=-gg(2:ell,2:ell)
 CALL invert(hh)
 p=matmul(hh,gg(2:ell,1))
 q=matmul(hh,gg(2:ell,ell+1))
 gamma0(1)=one
 gamma0(ell+1)=zero
 gamma0(2:ell)=p
 gamma1(1)=zero
 gamma1(ell+1)=one
 gamma1(2:ell)=q
 ngamma0=DOT_PRODUCT(gamma0,MATMUL(gg,gamma0))
 ngamma1=DOT_PRODUCT(gamma1,MATMUL(gg,gamma1))
 omega=DOT_PRODUCT(gamma0,MATMUL(gg,gamma1))
 cosine=ABS(omega)/SQRT(ABS(ngamma0*ngamma1))
 omega=omega/ngamma1
 IF(cosine<kappa)omega=(kappa/cosine)*omega
 gamma=gamma0-omega*gamma1
 s(1:ell)=gamma(2:ell+1)
 s(ell+1)=zero
RETURN
CONTAINS
SUBROUTINE invert(matrix)
!
! This subroutine inverts a small square matrix onto itself.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN OUT)::matrix(:,:)
 REAL(iwp)::det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
 INTEGER::ndim,i,k
 ndim=UBOUND(matrix,1)
 IF(ndim==2)THEN
   det=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
   j11=matrix(1,1)
   matrix(1,1)=matrix(2,2)
   matrix(2,2)=j11
   matrix(1,2)=-matrix(1,2)
   matrix(2,1)=-matrix(2,1)
   matrix=matrix/det
 ELSE IF(ndim==3)THEN
   det=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3))
   det=det-matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
   det=det+matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2))
   j11=matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
   j21=-matrix(2,1)*matrix(3,3)+matrix(3,1)*matrix(2,3)
   j31=matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
   j12=-matrix(1,2)*matrix(3,3)+matrix(3,2)*matrix(1,3)
   j22=matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
   j32=-matrix(1,1)*matrix(3,2)+matrix(3,1)*matrix(1,2)
   j13=matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
   j23=-matrix(1,1)*matrix(2,3)+matrix(2,1)*matrix(1,3)
   j33=matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
   matrix(1,1)=j11
   matrix(1,2)=j12
   matrix(1,3)=j13
   matrix(2,1)=j21
   matrix(2,2)=j22
   matrix(2,3)=j23
   matrix(3,1)=j31
   matrix(3,2)=j32
   matrix(3,3)=j33
   matrix=matrix/det
 ELSE
   DO k=1,ndim
     con=matrix(k,k)
     matrix(k,k)=1.0_iwp
     matrix(k,:)=matrix(k,:)/con
     DO i=1,ndim
       IF(i/=k)THEN
         con=matrix(i,k)
         matrix(i,k)=0.0_iwp
         matrix(i,:)=matrix(i,:)-matrix(k,:)*con
       END IF
     END DO
   END DO
 END IF
RETURN
END SUBROUTINE invert
END SUBROUTINE form_s










SUBROUTINE fsparv(kv,km,g,kdiag)
!
! This subroutine assembles element matrices into a symmetric skyline
! global matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::g(:),kdiag(:)
 REAL(iwp),INTENT(IN)::km(:,:)
 REAL(iwp),INTENT(OUT)::kv(:) 
 INTEGER::i,idof,k,j,iw,ival
 idof=UBOUND(g,1)
 DO i=1,idof
   k=g(i)
   IF(k/=0)THEN
     DO j=1,idof
       IF(g(j)/=0)THEN
         iw=k-g(j)
         IF(iw>=0)THEN
           ival=kdiag(k)-iw
           kv(ival)=kv(ival)+km(i,j) 
         END IF
       END IF
     END DO
   END IF
 END DO
RETURN
END SUBROUTINE fsparv

SUBROUTINE gauss_band(pb,work)
!
! This subroutine performs gaussian reduction of an unsymmetric
! banded matrix pb. Array work used as working space.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN OUT)::pb(:,:),work(:,:)
 REAL(iwp)::s,zero=0.0_iwp,small=1.e-10_iwp
 INTEGER::n,iwp1,iq,iqp,iwp11,i,j,k,l,ip,k1
 n=UBOUND(pb,1)
 iwp1=(UBOUND(pb,2)-1)/2+1
 iq=2*iwp1-1
 iqp=iwp1
 iwp11=iwp1-1
 DO i=1,iwp11
   DO j=1,iq
     IF(j>=iwp1+i)THEN
       pb(i,j)=zero
       pb(n-i+1,j)=zero
     ELSE 
       pb(i,j)=pb(i,j+iwp1-i)
     END IF
   END DO
 END DO
 DO k=1,n
   l=k+iwp1-1
   IF(l>n)l=n
   ip=0
   s=small
   DO i=k,l
     IF(ABS(pb(i,1))<=s)CYCLE
     s=ABS(pb(i,1))
     ip=i
   END DO
   IF(ip==0)THEN
     WRITE(6,'("singular")')
     EXIT
   END IF
   IF(k==n)EXIT
   work(iwp1,k)=ip
   iqp=iqp-1
   j=iwp1+ip-k
   IF(iqp<j)iqp=j
   IF(j/=iwp1)THEN
     DO j=1,iqp
       s=pb(k,j)
       pb(k,j)=pb(ip,j)
       pb(ip,j)=s
     END DO
   END IF
   k1=k+1
   DO i=k1,l
     s=pb(i,1)/pb(k,1)
     DO j=2,iq
       IF(j>iqp)THEN
         pb(i,j-1)=pb(i,j)
       ELSE 
         pb(i,j-1)=pb(i,j)-s*pb(k,j)
       END IF
     END DO
     pb(i,iq)=zero
     work(i-k,k)=s
   END DO
 END DO
RETURN
END SUBROUTINE gauss_band

SUBROUTINE getname(argv,nlen)
!
! This subroutine the base name of data file.
!
 IMPLICIT NONE
 INTEGER::narg
 INTEGER,INTENT(OUT)::nlen
 INTEGER::lnblnk,iargc
 CHARACTER(*),INTENT(OUT)::argv
 LOGICAL found
 narg=iargc()
 IF(narg.lt.1)THEN
   WRITE(*,*)'Please enter the base name of data file: '
   READ(*,*) argv
  ELSE
   CALL getarg(1,argv)
 ENDIF
 nlen=lnblnk(argv)
 INQUIRE(file=argv(1:nlen)//'.dat',exist=found)
 IF(.not.found)THEN
  WRITE(*,*)'Data file not found: ',argv(1:nlen)//'.dat'
  WRITE(*,*)'Please create or check spelling.'
  STOP
 ENDIF
RETURN
END SUBROUTINE getname
SUBROUTINE glob_to_axial(axial,global,coord)
!
! This subroutine transforms the global end reactions
! into an axial force for rod elements (2- or 3-d).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::global(:),coord(:,:)
 REAL(iwp),INTENT(OUT)::axial
 REAL(iwp)::add,ell,zero=0.0_iwp
 INTEGER::ndim,i
 ndim=UBOUND(coord,2)
 add=zero
 DO i=1,ndim
   add=add+(coord(2,i)-coord(1,i))**2
 END DO
 ell=SQRT(add)
 axial=zero
 DO i=1,ndim
   axial=axial+(coord(2,i)-coord(1,i))/ell*global(ndim+i)
 END DO
RETURN
END SUBROUTINE glob_to_axial                       
SUBROUTINE glob_to_loc(local,global,gamma,coord)
!
! This subroutine transforms the global end reactions and
! moments into the local system (2- or 3-d).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::global(:),gamma,coord(:,:)
 REAL(iwp),INTENT(OUT)::local(:)
 REAL(iwp)::t(12,12),r0(3,3),x1,x2,y1,y2,z1,z2,xl,yl,zl,pi,gamrad,cg,sg,  &
   den,ell,x,sum,zero=0.0_iwp,one=1.0_iwp,d180=180.0_iwp
 INTEGER::i,j,k,ndim
 ndim=UBOUND(coord,2)
 SELECT CASE(ndim)
 CASE(2)
   x1=coord(1,1)
   y1=coord(1,2)
   x2=coord(2,1)
   y2=coord(2,2)
   ell=SQRT((x2-x1)**2+(y2-y1)**2)
   cg=(x2-x1)/ell
   sg=(y2-y1)/ell
   local(1)=cg*global(1)+sg*global(2)
   local(2)=cg*global(2)-sg*global(1)
   local(3)=global(3)
   local(4)=cg*global(4)+sg*global(5)
   local(5)=cg*global(5)-sg*global(4)
   local(6)=global(6)
 CASE(3)
   x1=coord(1,1)
   y1=coord(1,2)
   z1=coord(1,3)
   x2=coord(2,1)
   y2=coord(2,2)
   z2=coord(2,3)
   xl=x2-x1
   yl=y2-y1
   zl=z2-z1
   ell=SQRT(xl*xl+yl*yl+zl*zl)
   t=zero
   pi=ACOS(-one)
   gamrad=gamma*pi/d180
   cg=COS(gamrad)
   sg=SIN(gamrad)
   den=ell*SQRT(xl*xl+zl*zl)
   IF(den/=zero)THEN
     r0(1,1)=xl/ell
     r0(1,2)=yl/ell
     r0(1,3)=zl/ell
     r0(2,1)=(-xl*yl*cg-ell*zl*sg)/den
     r0(2,2)=den*cg/(ell*ell)
     r0(2,3)=(-yl*zl*cg+ell*xl*sg)/den
     r0(3,1)=(xl*yl*sg-ell*zl*cg)/den
     r0(3,2)=-den*sg/(ell*ell)
     r0(3,3)=(yl*zl*sg+ell*xl*cg)/den
   ELSE
     r0(1,1)=zero
     r0(1,3)=zero
     r0(2,2)=zero
     r0(3,2)=zero
     r0(1,2)=one
     r0(2,1)=-cg
     r0(3,3)=cg
     r0(2,3)=sg
     r0(3,1)=sg
   END IF
   DO i=1,3
     DO j=1,3 
       x=r0(i,j)
       DO k=0,9,3
         t(i+k,j+k)=x
       END DO
     END DO
   END DO
   DO i=1,12
     sum=zero
     DO j=1,12
       sum=sum+t(i,j)*global(j)
     END DO
     local(i)=sum
   END DO
 END SELECT
RETURN
END SUBROUTINE glob_to_loc    
SUBROUTINE hinge(coord,holdr,action,react,prop,iel,etype,gamma)
!
! This subroutine forms the end forces and moments to be
! applied to a member if a joint has gone plastic.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::holdr(:,:),coord(:,:),action(:),prop(:,:),gamma(:)
 REAL(iwp),INTENT(OUT)::react(:)
 INTEGER,INTENT(IN)::etype(:),iel
 REAL(iwp)::ell,x1,x2,y1,y2,z1,z2,csch,snch,bm1,bm2,bm3,bm4,s1,s2,s3,s4,  &
   s5,mpy,mpz,mpx,gam,zero=0.0_iwp
 REAL(iwp),ALLOCATABLE::global(:),local(:),total(:)
 INTEGER::ndim,ndof
 ndim=UBOUND(coord,2)
 ndof=UBOUND(action,1)
 ALLOCATE(global(ndof),local(ndof),total(ndof))
 bm1=zero
 bm2=zero
 bm3=zero
 bm4=zero
 total(:)=holdr(:,iel)
 SELECT CASE(ndim)
 CASE(1)
   mpy=prop(2,etype(iel))
   ell=coord(2,1)-coord(1,1)
   s1=total(2)+action(2)
   s2=total(4)+action(4)
   IF(ABS(s1)>mpy)THEN
     if(s1> zero)bm1= mpy-s1
     if(s1<=zero)bm1=-mpy-s1
   END IF
   IF(ABS(s2)>mpy)THEN
     IF(s2> zero)bm2= mpy-s2
     IF(s2<=zero)bm2=-mpy-s2
   END IF
   react(1)= (bm1+bm2)/ell
   react(2)=bm1
   react(3)=-react(1)
   react(4)=bm2
 CASE(2)
   mpy=prop(3,etype(iel))  
   x1=coord(1,1)
   y1=coord(1,2)
   x2=coord(2,1)
   y2=coord(2,2)
   ell=SQRT((y2-y1)**2+(x2-x1)**2)
   csch=(x2-x1)/ell
   snch=(y2-y1)/ell
   s1=total(3)+action(3)
   s2=total(6)+action(6)
   IF(ABS(s1)>mpy)THEN 
     IF(s1> zero)bm1= mpy-s1
     IF(s1<=zero)bm1=-mpy-s1
   END IF
   IF(ABS(s2)>mpy)THEN
     IF(s2> zero)bm2= mpy-s2
     IF(s2<=zero)bm2=-mpy-s2
   END IF
   react(1)=-(bm1+bm2)*snch/ell
   react(2)= (bm1+bm2)*csch/ell
   react(3)=bm1
   react(4)=-react(1)
   react(5)=-react(2)
   react(6)=bm2
 CASE(3)
   gam=gamma(iel)
   mpy=prop(5,etype(iel))
   mpz=prop(6,etype(iel))
   mpx=prop(7,etype(iel))
   x1=coord(1,1)
   y1=coord(1,2)
   z1=coord(1,3)
   x2=coord(2,1)
   y2=coord(2,2)
   z2=coord(2,3)
   ell=SQRT((z2-z1)**2+(y2-y1)**2+(x2-x1)**2)
   global=total+action
   call glob_to_loc(local,global,gam,coord)
   global=zero
   s1=local(5)
   s2=local(11)
   IF(ABS(s1)>mpy)THEN 
     IF(s1> zero)bm1= mpy-s1
     IF(s1<=zero)bm1=-mpy-s1
   END IF
   IF(ABS(s2)>mpy)THEN
     IF(s2> zero)bm2= mpy-s2
     IF(s2<=zero)bm2=-mpy-s2
   END IF
   local( 3)=-(bm1+bm2)/ell
   local( 9)=-local(3)
   local( 5)= bm1
   local(11)= bm2
   s3=local(6)
   s4=local(12)
   IF(ABS(s3)>mpz)THEN
     IF(s3> zero)bm1= mpz-s3
     IF(s3<=zero)bm1=-mpz-s3
   END IF
   IF(ABS(s4)>mpy)THEN
     IF(s4> zero)bm2= mpz-s4
     IF(s4<=zero)bm2=-mpz-s4
   END IF
   local( 2)=(bm3+bm4)/ell
   local( 8)=-local(2)
   local( 6)= bm3
   local(12)= bm4
   s5=local(4)
   IF(ABS(s5)>mpx)THEN
     IF(s5> zero)global(4)= mpx-s5
     IF(s5<=zero)global(4)=-mpx-s5
   END IF
   local(10)=-local(4)
   CALL loc_to_glob(local,react,gam,coord)
 END SELECT
RETURN
CONTAINS
SUBROUTINE loc_to_glob(local,global,gamma,coord)
!
! This subroutine transforms the local end reactions and
! moments into the global system (3-d).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::local(:),gamma,coord(:,:)
 REAL(iwp),INTENT(OUT)::global(:)
 REAL(iwp)::t(12,12),r0(3,3),x1,x2,y1,y2,z1,z2,xl,yl,zl,pi,gamrad,cg,sg,  &
   den,ell,x,sum,zero=0.0_iwp,one=1.0_iwp,d180=180.0_iwp
 INTEGER::i,j,k
 x1=coord(1,1)
 y1=coord(1,2)
 z1=coord(1,3)
 x2=coord(2,1)
 y2=coord(2,2)
 z2=coord(2,3)
 xl=x2-x1
 yl=y2-y1
 zl=z2-z1
 ell=SQRT(xl*xl+yl*yl+zl*zl)
 t=zero
 pi=ACOS(-one)
 gamrad=gamma*pi/d180
 cg=COS(gamrad)
 sg=SIN(gamrad)
 den=ell*SQRT(xl*xl+zl*zl)
 IF(den/=zero)THEN
   r0(1,1)=xl/ell
   r0(2,1)=yl/ell
   r0(3,1)=zl/ell
   r0(1,2)=(-xl*yl*cg-ell*zl*sg)/den
   r0(2,2)=den*cg/(ell*ell)
   r0(3,2)=(-yl*zl*cg+ell*xl*sg)/den
   r0(1,3)=(xl*yl*sg-ell*zl*cg)/den
   r0(2,3)=-den*sg/(ell*ell)
   r0(3,3)=(yl*zl*sg+ell*xl*cg)/den
 ELSE
   r0(1,1)=zero
   r0(3,1)=zero
   r0(2,2)=zero
   r0(2,3)=zero
   r0(2,1)=one
   r0(1,2)=-cg
   r0(3,3)=cg
   r0(3,2)=sg
   r0(1,3)=sg
 END IF
 DO i=1,3
   DO j=1,3 
     x=r0(i,j)
     DO k=0,9,3
       t(i+k,j+k)=x
     END DO
   END DO
 END DO
 DO i=1,12
   sum=zero
   DO j=1,12
     sum=sum+t(i,j)*local(j)
   END DO
   global(i)=sum
 END DO
RETURN
END SUBROUTINE loc_to_glob    
SUBROUTINE glob_to_loc(local,global,gamma,coord)
!
! This subroutine transforms the global end reactions and
! moments into the local system (2- or 3-d).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::global(:),gamma,coord(:,:)
 REAL(iwp),INTENT(OUT)::local(:)
 REAL(iwp)::t(12,12),r0(3,3),x1,x2,y1,y2,z1,z2,xl,yl,zl,pi,gamrad,cg,sg,  &
   den,ell,x,sum,zero=0.0_iwp,one=1.0_iwp,d180=180.0_iwp
 INTEGER::i,j,k,ndim
 ndim=UBOUND(coord,2)
 SELECT CASE(ndim)
 CASE(2)
   x1=coord(1,1)
   y1=coord(1,2)
   x2=coord(2,1)
   y2=coord(2,2)
   ell=SQRT((x2-x1)**2+(y2-y1)**2)
   cg=(x2-x1)/ell
   sg=(y2-y1)/ell
   local(1)=cg*global(1)+sg*global(2)
   local(2)=cg*global(2)-sg*global(1)
   local(3)=global(3)
   local(4)=cg*global(4)+sg*global(5)
   local(5)=cg*global(5)-sg*global(4)
   local(6)=global(6)
 CASE(3)
   x1=coord(1,1)
   y1=coord(1,2)
   z1=coord(1,3)
   x2=coord(2,1)
   y2=coord(2,2)
   z2=coord(2,3)
   xl=x2-x1
   yl=y2-y1
   zl=z2-z1
   ell=SQRT(xl*xl+yl*yl+zl*zl)
   t=zero
   pi=ACOS(-one)
   gamrad=gamma*pi/d180
   cg=COS(gamrad)
   sg=SIN(gamrad)
   den=ell*SQRT(xl*xl+zl*zl)
   IF(den/=zero)THEN
     r0(1,1)=xl/ell
     r0(1,2)=yl/ell
     r0(1,3)=zl/ell
     r0(2,1)=(-xl*yl*cg-ell*zl*sg)/den
     r0(2,2)=den*cg/(ell*ell)
     r0(2,3)=(-yl*zl*cg+ell*xl*sg)/den
     r0(3,1)=(xl*yl*sg-ell*zl*cg)/den
     r0(3,2)=-den*sg/(ell*ell)
     r0(3,3)=(yl*zl*sg+ell*xl*cg)/den
   ELSE
     r0(1,1)=zero
     r0(1,3)=zero
     r0(2,2)=zero
     r0(3,2)=zero
     r0(1,2)=one
     r0(2,1)=-cg
     r0(3,3)=cg
     r0(2,3)=sg
     r0(3,1)=sg
   END IF
   DO i=1,3
     DO j=1,3 
       x=r0(i,j)
       DO k=0,9,3
         t(i+k,j+k)=x
       END DO
     END DO
   END DO
   DO i=1,12
     sum=zero
     DO j=1,12
       sum=sum+t(i,j)*global(j)
     END DO
     local(i)=sum
   END DO
 END SELECT
RETURN
END SUBROUTINE glob_to_loc    
END SUBROUTINE hinge
SUBROUTINE interp(k,dtim,rt,rl,al,nstep)
!
! This subroutine forms the load/time functions by interpolation.
! If dtim is not an exact multiple it stops one short.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::dtim,rt(:),rl(:)
 INTEGER,INTENT(IN)::k,nstep
 REAL(iwp),INTENT(IN OUT)::al(:,:)
 INTEGER::np,i,j
 REAL(iwp)::t,val 
 np=SIZE(rt)
 al(1,k)=rl(1)
 t=rt(1)
 DO j=2,nstep
   t=t+dtim
   DO i=2,np
     IF(t.LE.rt(i))THEN
       val=rl(i-1)+((t-rt(i-1))/(rt(i)-rt(i-1)))*(rl(i)-rl(i-1))
       EXIT
     END IF
   END DO
   al(j,k)=val
 END DO
RETURN
END SUBROUTINE interp
SUBROUTINE invar(stress,sigm,dsbar,theta)
!
! This subroutine forms the stress invariants in 2- or 3-d.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::stress(:)
 REAL(iwp),INTENT(OUT),OPTIONAL::sigm,dsbar,theta
 REAL(iwp)::sx,sy,sz,txy,dx,dy,dz,xj3,sine,s1,s2,s3,s4,s5,s6,ds1,ds2,ds3, &
   d2,d3,sq3,zero=0.0_iwp,small=1.e-10_iwp,one=1.0_iwp,two=2.0_iwp,       &
   three=3.0_iwp,six=6.0_iwp,thpt5=13.5_iwp
 INTEGER::nst 
 nst=UBOUND(stress,1)
 SELECT CASE(nst)
 CASE(4)
   sx=stress(1)
   sy=stress(2)
   txy=stress(3)
   sz=stress(4)
   sigm=(sx+sy+sz)/three
   dsbar=SQRT((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+six*txy**2)/SQRT(two)
   IF(dsbar<small)THEN
     theta=zero
   ELSE
     dx=(two*sx-sy-sz)/three
     dy=(two*sy-sz-sx)/three
     dz=(two*sz-sx-sy)/three
     xj3=dx*dy*dz-dz*txy**2
     sine=-thpt5*xj3/dsbar**3
     IF(sine>=one)sine=one
     IF(sine<-one)sine=-one
     theta=ASIN(sine)/three
   END IF
 CASE(6)
   sq3=SQRT(three)
   s1=stress(1)  
   s2=stress(2)
   s3=stress(3) 
   s4=stress(4)
   s5=stress(5)
   s6=stress(6)
   sigm=(s1+s2+s3)/three
   d2=((s1-s2)**2+(s2-s3)**2+(s3-s1)**2)/six+s4*s4+s5*s5+s6*s6
   ds1=s1-sigm 
   ds2=s2-sigm  
   ds3=s3-sigm
   d3=ds1*ds2*ds3-ds1*s5*s5-ds2*s6*s6-ds3*s4*s4+two*s4*s5*s6
   dsbar=sq3*SQRT(d2)
   IF(dsbar<small)THEN
     theta=zero
   ELSE
     sine=-three*sq3*d3/(two*SQRT(d2)**3)
     IF(sine>=one)sine=one 
     IF(sine<-one)sine=-one 
     theta=ASIN(sine)/three
   END IF
 CASE DEFAULT
   WRITE(*,*)"wrong size for nst in invar"
 END SELECT
RETURN
END SUBROUTINE invar



SUBROUTINE invert(matrix)
!
! This subroutine inverts a small square matrix onto itself.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN OUT)::matrix(:,:)
 REAL(iwp)::det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
 INTEGER::ndim,i,k
 ndim=UBOUND(matrix,1)
 IF(ndim==2)THEN
   det=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
   j11=matrix(1,1)
   matrix(1,1)=matrix(2,2)
   matrix(2,2)=j11
   matrix(1,2)=-matrix(1,2)
   matrix(2,1)=-matrix(2,1)
   matrix=matrix/det
 ELSE IF(ndim==3)THEN
   det=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3))
   det=det-matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
   det=det+matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2))
   j11=matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
   j21=-matrix(2,1)*matrix(3,3)+matrix(3,1)*matrix(2,3)
   j31=matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
   j12=-matrix(1,2)*matrix(3,3)+matrix(3,2)*matrix(1,3)
   j22=matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
   j32=-matrix(1,1)*matrix(3,2)+matrix(3,1)*matrix(1,2)
   j13=matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
   j23=-matrix(1,1)*matrix(2,3)+matrix(2,1)*matrix(1,3)
   j33=matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
   matrix(1,1)=j11
   matrix(1,2)=j12
   matrix(1,3)=j13
   matrix(2,1)=j21
   matrix(2,2)=j22
   matrix(2,3)=j23
   matrix(3,1)=j31
   matrix(3,2)=j32
   matrix(3,3)=j33
   matrix=matrix/det
 ELSE
   DO k=1,ndim
     con=matrix(k,k)
     matrix(k,k)=1.0_iwp
     matrix(k,:)=matrix(k,:)/con
     DO i=1,ndim
       IF(i/=k)THEN
         con=matrix(i,k)
         matrix(i,k)=0.0_iwp
         matrix(i,:)=matrix(i,:)-matrix(k,:)*con
       END IF
     END DO
   END DO
 END IF
RETURN
END SUBROUTINE invert







SUBROUTINE linmul_sky(kv,disps,loads,kdiag)
!
! This subroutine forms the product of symmetric matrix stored as
! a skyline and a vector.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::kv(:),disps(0:)
 REAL(iwp),INTENT(OUT)::loads(0:)
 INTEGER,INTENT(IN)::kdiag(:)
 INTEGER::n,i,j,low,lup,k
 REAL(iwp)::x,zero=0.0_iwp
 n=UBOUND(disps,1)
 DO i=1,n
   x=zero 
   lup=kdiag(i)
   IF(i==1)low=lup
   IF(i/=1)low=kdiag(i-1)+1
   DO j=low,lup
     x=x+kv(j)*disps(i+j-lup) 
   END DO
   loads(i)=x
   IF(i==1)CYCLE   
   lup=lup-1
   DO j=low,lup
     k=i+j-lup-1
     loads(k)=loads(k)+kv(j)*disps(i)        
   END DO
 END DO
RETURN
END SUBROUTINE linmul_sky 
 
INTEGER FUNCTION lnblnk( str )
 CHARACTER*(*),INTENT(IN)::str
 CHARACTER*1 space, tab, null, char
 INTEGER j,i
 DATA space/' '/, tab/'	'/
 null = CHAR(0)
 i = LEN( str )
 DO j = i, 1, -1
  IF(str(j:j).ne.space.and.str(j:j).ne.null.and.str(j:j).ne.tab)THEN
   lnblnk = j
   RETURN
  ENDIF
 ENDDO
 lnblnk = 0
 RETURN
END FUNCTION lnblnk


SUBROUTINE load_function(lf,dtim,al)
!
! This subroutine forms the increment of load at each
! calculation time step.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::lf(:,:),dtim
 REAL(iwp),INTENT(IN OUT)::al(:)
 REAL::time,aold,anew
 INTEGER::nlfp,nincs,i,j
 nincs=SIZE(al)
 nlfp=UBOUND(lf,2)
 aold=lf(2,1)                        
 time=lf(1,1)
 DO i=1,nincs
   time=time+dtim
   DO j=1,nlfp
     IF(time<lf(1,j))THEN
       anew=lf(2,j-1)+                                                    &
         (time-lf(1,j-1))*(lf(2,j)-lf(2,j-1))/(lf(1,j)-lf(1,j-1))
       al(i)=anew-aold
       aold=anew
       EXIT
     END IF
   END DO
 END DO
RETURN
END SUBROUTINE load_function
SUBROUTINE loc_to_glob(local,global,gamma,coord)
!
! This subroutine transforms the local end reactions and
! moments into the global system (3-d).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::local(:),gamma,coord(:,:)
 REAL(iwp),INTENT(OUT)::global(:)
 REAL(iwp)::t(12,12),r0(3,3),x1,x2,y1,y2,z1,z2,xl,yl,zl,pi,gamrad,cg,sg,  &
   den,ell,x,sum,zero=0.0_iwp,one=1.0_iwp,d180=180.0_iwp
 INTEGER::i,j,k
 x1=coord(1,1)
 y1=coord(1,2)
 z1=coord(1,3)
 x2=coord(2,1)
 y2=coord(2,2)
 z2=coord(2,3)
 xl=x2-x1
 yl=y2-y1
 zl=z2-z1
 ell=SQRT(xl*xl+yl*yl+zl*zl)
 t=zero
 pi=ACOS(-one)
 gamrad=gamma*pi/d180
 cg=COS(gamrad)
 sg=SIN(gamrad)
 den=ell*SQRT(xl*xl+zl*zl)
 IF(den/=zero)THEN
   r0(1,1)=xl/ell
   r0(2,1)=yl/ell
   r0(3,1)=zl/ell
   r0(1,2)=(-xl*yl*cg-ell*zl*sg)/den
   r0(2,2)=den*cg/(ell*ell)
   r0(3,2)=(-yl*zl*cg+ell*xl*sg)/den
   r0(1,3)=(xl*yl*sg-ell*zl*cg)/den
   r0(2,3)=-den*sg/(ell*ell)
   r0(3,3)=(yl*zl*sg+ell*xl*cg)/den
 ELSE
   r0(1,1)=zero
   r0(3,1)=zero
   r0(2,2)=zero
   r0(2,3)=zero
   r0(2,1)=one
   r0(1,2)=-cg
   r0(3,3)=cg
   r0(3,2)=sg
   r0(1,3)=sg
 END IF
 DO i=1,3
   DO j=1,3 
     x=r0(i,j)
     DO k=0,9,3
       t(i+k,j+k)=x
     END DO
   END DO
 END DO
 DO i=1,12
   sum=zero
   DO j=1,12
     sum=sum+t(i,j)*local(j)
   END DO
   global(i)=sum
 END DO
RETURN
END SUBROUTINE loc_to_glob    
SUBROUTINE mcdpl(phi,psi,dee,stress,pl)
!
! This subroutine forms the plastic stress/strain matrix
! for a Mohr-Coulomb material (phi,psi in degrees).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::stress(:),dee(:,:),phi,psi  
 REAL(iwp),INTENT(OUT)::pl(:,:)
 REAL(iwp),ALLOCATABLE::dfds(:),dqds(:),ddqds(:),dfdsd(:)
 REAL(iwp)::t1,t2,t3,t4,t5,t6,t8,t10,t12,t13,t14,t15,t16,t17,t18,t19,t20, &
   t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,   &
   t38,t39,t40,t41,t42,t43,t44,t45,t46,t48,t49,t50,t51,t53,t54,t55,t56,   &
   t60,t61,t63,t64,t68,t69,t70,t71,t73,t74,t77,t79,t80,t82,t83,t84,t85,   &
   t86,t89,t92,t93,t94,t97,t98,t101,t103,t106,t110,t111,t113,t122,t129,   &
   t133,t140,t145,t152,t166,t186,t206,pm,pi,phir,snph,snth,sq3,sx,sy,sz,  &
   txy,tyz,tzx,zero=0.0_iwp,pt49=0.49_iwp,one=1.0_iwp,two=2.0_iwp,        &
   d3=3.0_iwp,d4=4.0_iwp,d6=6.0_iwp,d180=180.0_iwp,snps,psir
 REAL(iwp)::denom
 INTEGER::i,j,ih
 ih=SIZE(stress)
 ALLOCATE(dfds(ih),dqds(ih),ddqds(ih),dfdsd(ih))
 pi=ACOS(-one) 
 phir=phi*pi/d180
 snph=SIN(phir)
 psir=psi*pi/d180
 snps=SIN(psir)
 sq3=SQRT(d3)
 SELECT CASE(ih)
 CASE(4)
   sx=stress(1)
   sy=stress(2)
   txy=stress(3) 
   sz=stress(4)   
   t3=one/d3
   t4=(sx+sy+sz)*t3
   t8=sz-t4
   t10=txy**2
   t16=(sx-sy)**2
   t18=(sy-sz)**2
   t20=(sz-sx)**2
   t25=SQRT((t16+t18+t20)/d6+t10)
   t26=t25**2
   t30=d3*sq3*((sx-t4)*(sy-t4)*t8-t8*t10)/two/t26/t25
   IF(t30>one)t30=one
   IF(t30<-one)t30=-one
   t31=ASIN(t30)
   t33=SIN(t31*t3)
   snth=-t33
   IF(ABS(snth).GT.pt49)THEN
     pm=-one
     if(snth.LT.zero)pm=one
     t2=snph/d3
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t14=SQRT((t4+t6+t8)*t10+t12)
     t17=one/t14/sq3
     t18=one/two
     t19=t17*t18
     t21=d3+pm*snph
     t23=two*sy
     t24=two*sz
     t31=two*sx
     dfds(1)=t2+t19*t21*(d4*sx-t23-t24)*t10/two
     dfds(2)=t2+t19*t21*(-t31+d4*sy-t24)*t10/two
     dfds(3)=t17*t18*t21*txy
     dfds(4)=t2+t19*t21*(-t23+d4*sz-t31)*t10/two
     t2=snps/d3
     t21=d3+pm*snps
     dqds(1)=t2+t19*t21*(d4*sx-t23-t24)*t10/two
     dqds(2)=t2+t19*t21*(-t31+d4*sy-t24)*t10/two
     dqds(3)=t17*t18*t21*txy
     dqds(4)=t2+t19*t21*(-t23+d4*sz-t31)*t10/two
   ELSE
     t1=one/d3
     t2=snph*t1
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t13=(t4+t6+t8)*t10+t12
     t14=SQRT(t13)
     t16=d3*sq3
     t18=(sx+sy+sz)*t1
     t19=sx-t18
     t20=sy-t18
     t21=t19*t20
     t22=sz-t18
     t25=t21*t22-t22*t12
     t26=one/two
     t28=t14**2
     t30=one/t28/t14
     t31=t16*t25*t26*t30
     IF(t31>one)t31=one
     IF(t31<-one)t31=-one
     t33=ASIN(t31)
     t34=t33*t1
     t35=COS(t34)
     t36=SIN(t34)
     t38=one/sq3
     t41=one/t14*(t35+t36*snph*t38)
     t43=two*sy
     t44=two*sz
     t46=(d4*sx-t43-t44)*t10
     t49=one-t1
     t53=t19*t1*t22
     t54=t21*t1
     t55=t1*t12
     t60=t16*t25
     t61=t28**2
     t64=t26/t61/t14
     t68=t16*(t49*t20*t22-t53-t54+t55)*t26*t30-d3/two*t60*t64*t46
     t70=d3**2
     t71=sq3**2
     t73=t25**2
     t74=two**2
     t77=t13**2
     t83=SQRT(one-t70*t71*t73/t74/t77/t13)
     t84=one/t83
     t85=t84*t1
     t89=t2*t38
     t94=two*sx
     t97=(-t94+d4*sy-t44)*t10
     t101=t1*t20*t22
     t111=t16*(-t101+t19*t49*t22-t54+t55)*t26*t30-d3/two*t60*t64*t97
     t129=-two*t16*t22*txy*t26*t30-d3*t60*t64*txy
     t140=(-t43+d4*sz-t94)*t10
     t152=t16*(-t101-t53+t21*t49-t49*t12)*t26*t30-d3/two*t60*t64*t140
     dfds(1)=t2+t41*t46/two+t14*(-t36*t68*t85+t35*t68*t84*t89)
     dfds(2)=t2+t41*t97/two+t14*(-t36*t111*t85+t35*t111*t84*t89)
     dfds(3)=t41*txy+t14*(-t36*t129*t85+t35*t129*t84*t89)
     dfds(4)=t2+t41*t140/two+t14*(-t36*t152*t85+t35*t152*t84*t89)
     t2=snps*t1
     t41=one/t14*(t35+t36*snps*t38)
     t89=t2*t38
     dqds(1)=t2+t41*t46/two+t14*(-t36*t68*t85+t35*t68*t84*t89)
     dqds(2)=t2+t41*t97/two+t14*(-t36*t111*t85+t35*t111*t84*t89)
     dqds(3)=t41*txy+t14*(-t36*t129*t85+t35*t129*t84*t89)
     dqds(4)=t2+t41*t140/two+t14*(-t36*t152*t85+t35*t152*t84*t89)
   END IF
 CASE(6)
   sx=stress(1)
   sy=stress(2)
   sz=stress(3)   
   txy=stress(4) 
   tyz=stress(5) 
   tzx=stress(6) 
   t3=one/d3
   t4=(sx+sy+sz)*t3
   t5=sx-t4
   t6=sy-t4
   t8=sz-t4
   t10=tyz**2
   t12=tzx**2
   t14=txy**2
   t23=(sx-sy)**2
   t25=(sy-sz)**2
   t27=(sz-sx)**2
   t32=SQRT((t23+t25+t27)/d6+t14+t10+t12)
   t33=t32**2
   t37=d3*sq3*(t5*t6*t8-t5*t10-t6*t12-t8*t14+two*txy*tyz*tzx)/two/t33/t32
   IF(t37>one)t37=one
   IF(t37<-one)t37=-one
   t38=ASIN(t37)
   t40=SIN(t38*t3)
   snth=-t40
   IF(ABS(snth).GT.pt49)THEN
     pm=-one
     IF(snth.LT.zero)pm=one  
     t2=snph/d3
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t13=tyz**2
     t14=tzx**2
     t16=SQRT((t4+t6+t8)*t10+t12+t13+t14)
     t19=one/t16/sq3
     t20=one/two
     t21=t19*t20
     t23=d3+pm*snph
     t25=two*sy
     t26=two*sz
     t33=two*sx
     t48=t20*t23
     dfds(1)=t2+t21*t23*(d4*sx-t25-t26)*t10/two
     dfds(2)=t2+t21*t23*(-t33+d4*sy-t26)*t10/two
     dfds(3)=t2+t21*t23*(-t25+d4*sz-t33)*t10/two
     dfds(4)=t19*t48*txy
     dfds(5)=t19*t48*tyz
     dfds(6)=t19*t48*tzx
     t2=snps/d3
     t23=d3+pm*snps
     t48=t20*t23
     dqds(1)=t2+t21*t23*(d4*sx-t25-t26)*t10/two
     dqds(2)=t2+t21*t23*(-t33+d4*sy-t26)*t10/two
     dqds(3)=t2+t21*t23*(-t25+d4*sz-t33)*t10/two
     dqds(4)=t19*t48*txy
     dqds(5)=t19*t48*tyz
     dqds(6)=t19*t48*tzx
   ELSE
     t1=one/d3
     t2=snph*t1
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t13=tyz**2
     t14=tzx**2
     t15=(t4+t6+t8)*t10+t12+t13+t14
     t16=SQRT(t15)
     t18=d3*sq3
     t20=(sx+sy+sz)*t1
     t21=sx-t20
     t22=sy-t20
     t23=t21*t22
     t24=sz-t20
     t29=two*txy
     t32=t23*t24-t21*t13-t22*t14-t24*t12+t29*tyz*tzx
     t33=one/two
     t35=t16**2
     t37=one/t35/t16
     t39=t18*t32*t33*t37
     IF(t39>one)t39=one
     IF(t39<-one)t39=-one
     t40=ASIN(t39)
     t41=t40*t1
     t42=COS(t41)
     t43=SIN(t41)
     t45=one/sq3
     t48=one/t16*(t42+t43*snph*t45)
     t50=two*sy
     t51=two*sz
     t53=(d4*sx-t50-t51)*t10
     t56=one-t1
     t60=t21*t1*t24
     t61=t23*t1
     t63=t1*t14
     t64=t1*t12
     t69=t18*t32
     t70=t35**2
     t73=t33/t70/t16
     t77=t18*(t56*t22*t24-t60-t61-t56*t13+t63+t64)*t33*t37-               &
       d3/two*t69*t73*t53
     t79=d3**2
     t80=sq3**2
     t82=t32**2
     t83=two**2
     t86=t15**2
     t92=SQRT(one-t79*t80*t82/t83/t86/t15)
     t93=one/t92
     t94=t93*t1
     t98=t2*t45
     t103=two*sx
     t106=(-t103+d4*sy-t51)*t10
     t110=t1*t22*t24
     t113=t1*t13
     t122=t18*(-t110+t21*t56*t24-t61+t113-t56*t14+t64)*t33*t37-           &
       d3/two*t69*t73*t106
     t133=(-t50+d4*sz-t103)*t10
     t145=t18*(-t110-t60+t23*t56+t113+t63-t56*t12)*t33*t37-               &
       d3/two*t69*t73*t133
     t166=t18*(-two*t24*txy+two*tyz*tzx)*t33*t37-d3*t69*t73*txy
     t186=t18*(-two*t21*tyz+t29*tzx)*t33*t37-d3*t69*t73*tyz
     t206=t18*(-two*t22*tzx+t29*tyz)*t33*t37-d3*t69*t73*tzx
     dfds(1)=t2+t48*t53/two+t16*(-t43*t77*t94+t42*t77*t93*t98)
     dfds(2)=t2+t48*t106/two+t16*(-t43*t122*t94+t42*t122*t93*t98)
     dfds(3)=t2+t48*t133/two+t16*(-t43*t145*t94+t42*t145*t93*t98)
     dfds(4)=t48*txy+t16*(-t43*t166*t94+t42*t166*t93*t98)
     dfds(5)=t48*tyz+t16*(-t43*t186*t94+t42*t186*t93*t98)
     dfds(6)=t48*tzx+t16*(-t43*t206*t94+t42*t206*t93*t98)
     t2=snps*t1
     t48=one/t16*(t42+t43*snps*t45)
     t98=t2*t45
     dqds(1)=t2+t48*t53/two+t16*(-t43*t77*t94+t42*t77*t93*t98)
     dqds(2)=t2+t48*t106/two+t16*(-t43*t122*t94+t42*t122*t93*t98)
     dqds(3)=t2+t48*t133/two+t16*(-t43*t145*t94+t42*t145*t93*t98)
     dqds(4)=t48*txy+t16*(-t43*t166*t94+t42*t166*t93*t98)
     dqds(5)=t48*tyz+t16*(-t43*t186*t94+t42*t186*t93*t98)
     dqds(6)=t48*tzx+t16*(-t43*t206*t94+t42*t206*t93*t98)
   END IF
 END SELECT
 ddqds=MATMUL(dee,dqds)
 dfdsd=MATMUL(dee,dfds) 
 denom=DOT_PRODUCT(dfdsd,dqds)
 DO i=1,ih
   DO j=1,ih
     pl(i,j)=ddqds(i)*dfdsd(j)/denom
   END DO
 END DO
 DEALLOCATE(dfds,dqds,ddqds,dfdsd)
RETURN
END SUBROUTINE mcdpl
SUBROUTINE mesh(g_coord,g_num,argv,nlen,ips)
!
! This subroutine produces a PostScript output file "*.msh" displaying
! the undeformed finite element mesh.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::g_coord(:,:)
 INTEGER,INTENT(IN)::g_num(:,:),ips,nlen
 CHARACTER(*),INTENT(IN)::argv
 REAL(iwp)::xmin,xmax,ymin,ymax,width,height,scale=72,sxy,xo,yo,x,y,      &
   pt5=0.5_iwp,opt5=1.5_iwp,fpt5=5.5_iwp,d8=8.0_iwp,ept5=8.5_iwp,         &
   d11=11.0_iwp
 INTEGER::i,ii,j,jj,nn,nod,nel
 OPEN(ips,FILE=argv(1:nlen)//'.msh')
!
!                       compute size of mesh
!
 nn=UBOUND(g_coord,2)
 xmin=g_coord(1,1)
 xmax=g_coord(1,1)
 ymin=g_coord(2,1)
 ymax=g_coord(2,1)
 DO i=2,nn
   IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)      
   IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)      
   IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)      
   IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)      
 END DO
 width =xmax-xmin
 height=ymax-ymin
!
!                       allow 1.5" margin minimum on each side of figure
!
 IF(height.GE.d11/ept5*width)THEN
!
!                       height governs the scale
!
   sxy=scale*d8/height
   xo=scale*pt5*(ept5-d8*width/height)
   yo=scale*opt5
 ELSE
!
!                       width governs the scale
!
   sxy=scale*fpt5/width
   xo=scale*opt5
   yo=scale*pt5*(d11-fpt5*height/width)
 END IF
!
!                       start PostScript output
!
 WRITE(ips,'(a)')'%!PS-Adobe-1.0'
 WRITE(ips,'(a)')'%%DocumentFonts: none'
 WRITE(ips,'(a)')'%%Pages: 1'
 WRITE(ips,'(a)')'%%EndComments'
 WRITE(ips,'(a)')'/m {moveto} def'
 WRITE(ips,'(a)')'/l {lineto} def'
 WRITE(ips,'(a)')'/s {stroke} def'
 WRITE(ips,'(a)')'/c {closepath} def'
 WRITE(ips,'(a)')'%%EndProlog'
 WRITE(ips,'(a)')'%%Page: 0 1'
 WRITE(ips,'(a)')'gsave'
 WRITE(ips,'(2f9.2,a)') xo, yo, ' translate'
 WRITE(ips,'(f9.2,a)') 0.5, ' setlinewidth'
!
!                       draw the mesh
!
 nod=UBOUND(g_num,1)
 nel=UBOUND(g_num,2)
 IF(nod==5)nod=4
 IF(nod==9)nod=8
 IF(nod==10)nod=9
 IF(nod==15)nod=12
 DO i=1,nel
   ii=g_num(1,i)
   IF(ii==0)CYCLE
   x=sxy*(g_coord(1,ii)-xmin)
   y=sxy*(g_coord(2,ii)-ymin)
   WRITE(ips,'(2f9.2,a)')x,y,' m'
   DO j=2,nod
     jj=g_num(j,i)
     x=sxy*(g_coord(1,jj)-xmin)
     y=sxy*(g_coord(2,jj)-ymin)
     WRITE(ips,'(2f9.2,a)') x, y,' l'
   END DO
   WRITE(ips,'(a)')'c s'
 END DO
!
!                       close output file
!
 WRITE(ips,'(a)')'grestore'
 WRITE(ips,'(a)')'showpage'
 CLOSE(ips)
!
RETURN
END SUBROUTINE mesh
SUBROUTINE mocouf(phi,c,sigm,dsbar,theta,f)
!
! This subroutine calculates the value of the yield function
! for a Mohr-Coulomb material (phi in degrees).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::phi,c,sigm,dsbar,theta   
 REAL(iwp),INTENT(OUT)::f
 REAL(iwp)::phir,snph,csph,csth,snth,one=1.0_iwp,d3=3.0_iwp,d4=4.0_iwp,   &
   d180=180.0_iwp
 phir=phi*d4*ATAN(one)/d180
 snph=SIN(phir) 
 csph=COS(phir) 
 csth=COS(theta)
 snth=SIN(theta)
 f=snph*sigm+dsbar*(csth/SQRT(d3)-snth*snph/d3)-c*csph
RETURN
END SUBROUTINE mocouf
SUBROUTINE mocouq(psi,dsbar,theta,dq1,dq2,dq3)
!
! This subroutine forms the derivatives of a Mohr-Coulomb potential
! function with respect to the three invariants (psi in degrees).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::psi,dsbar,theta
 REAL(iwp),INTENT(OUT)::dq1,dq2,dq3
 REAL(iwp)::psir,snth,snps,sq3,c1,csth,cs3th,tn3th,tnth,zero=0.0_iwp,     &
   pt49=0.49_iwp,pt5=0.5_iwp,one=1.0_iwp,d3=3.0_iwp,d4=4.0_iwp,           &
   d180=180.0_iwp
 psir=psi*d4*ATAN(one)/d180 
 snth=SIN(theta) 
 snps=SIN(psir)
 sq3=SQRT(d3)  
 dq1=snps
 if(ABS(snth).GT.pt49)THEN
   c1=one
   IF(snth.LT.zero)c1=-one
   dq2=(sq3*pt5-c1*snps*pt5/sq3)*sq3*pt5/dsbar 
   dq3=zero
 ELSE
   csth=COS(theta)
   cs3th=COS(d3*theta)
   tn3th=TAN(d3*theta)
   tnth=snth/csth
   dq2=sq3*csth/dsbar*((one+tnth*tn3th)+snps*(tn3th-tnth)/sq3)*pt5
   dq3=pt5*d3*(sq3*snth+snps*csth)/(cs3th*dsbar*dsbar)
 END IF
RETURN
END SUBROUTINE mocouq
FUNCTION norm(x)RESULT(l2n)
!
! THis function calculates the l2 norm of vector x
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::x(:)
 REAL(iwp)::l2n
 l2n=SQRT(SUM(x**2))
RETURN
END FUNCTION norm









SUBROUTINE num_to_g(num,nf,g)
!
! This subroutine finds the g vector from num and nf.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN)::num(:),nf(:,:)  
 INTEGER,INTENT(OUT)::g(:)
 INTEGER::i,k,nod,nodof 
 nod=UBOUND(num,1) 
 nodof=UBOUND(nf,1)
 DO i=1,nod
   k=i*nodof
   g(k-nodof+1:k)=nf(:,num(i))
 END DO
RETURN
END SUBROUTINE num_to_g   
SUBROUTINE pin_jointed(km,ea,coord)
!
! This subroutine forms the stiffness matrix of a
! general rod element (1-, 2- or 3-d).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::ea,coord(:,:)
 REAL(iwp),INTENT(OUT)::km(:,:)
 INTEGER::ndim
 REAL(iwp)::ell,cs,sn,x1,x2,y1,y2,z1,z2,a,b,c,d,e,f,xl,yl,zl,one=1.0_iwp
 ndim=UBOUND(coord,2)
 SELECT CASE(ndim)
 CASE(1)
   ell=coord(2,1)-coord(1,1)
   km(1,1)=one
   km(1,2)=-one
   km(2,1)=-one
   km(2,2)=one  
 CASE(2)
   x1=coord(1,1)
   y1=coord(1,2)
   x2=coord(2,1)
   y2=coord(2,2)
   ell=SQRT((y2-y1)**2+(x2-x1)**2)
   cs=(x2-x1)/ell
   sn=(y2-y1)/ell
   a=cs*cs
   b=sn*sn
   c=cs*sn
   km(1,1)=a
   km(3,3)=a
   km(1,3)=-a
   km(3,1)=-a
   km(2,2)=b
   km(4,4)=b
   km(2,4)=-b
   km(4,2)=-b
   km(1,2)=c
   km(2,1)=c
   km(3,4)=c
   km(4,3)=c
   km(1,4)=-c
   km(4,1)=-c
   km(2,3)=-c
   km(3,2)=-c
 CASE(3)
   x1=coord(1,1)
   y1=coord(1,2)
   z1=coord(1,3)
   x2=coord(2,1)
   y2=coord(2,2)
   z2=coord(2,3)
   xl=x2-x1
   yl=y2-y1
   zl=z2-z1
   ell=SQRT(xl*xl+yl*yl+zl*zl)
   xl=xl/ell
   yl=yl/ell
   zl=zl/ell
   a=xl*xl
   b=yl*yl
   c=zl*zl
   d=xl*yl
   e=yl*zl
   f=zl*xl
   km(1,1)=a
   km(4,4)=a
   km(2,2)=b
   km(5,5)=b
   km(3,3)=c
   km(6,6)=c
   km(1,2)=d
   km(2,1)=d
   km(4,5)=d
   km(5,4)=d
   km(2,3)=e
   km(3,2)=e
   km(5,6)=e
   km(6,5)=e
   km(1,3)=f
   km(3,1)=f
   km(4,6)=f
   km(6,4)=f
   km(1,4)=-a
   km(4,1)=-a
   km(2,5)=-b
   km(5,2)=-b
   km(3,6)=-c
   km(6,3)=-c
   km(1,5)=-d
   km(5,1)=-d
   km(2,4)=-d
   km(4,2)=-d
   km(2,6)=-e
   km(6,2)=-e
   km(3,5)=-e
   km(5,3)=-e
   km(1,6)=-f
   km(6,1)=-f
   km(3,4)=-f
   km(4,3)=-f
 END SELECT
 km=km*ea/ell
RETURN
END SUBROUTINE pin_jointed       
SUBROUTINE rect_km(km,coord,e,v)
!
! This subroutine forms the "analytical" stiffness matrix for
! rectangular 4- or 8-node plane strain elements using nip=4.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::coord(:,:),e,v
 REAL(iwp),INTENT(OUT)::km(:,:)
 INTEGER::nod,i,j
 REAL(iwp)::aa,bb,t2,t3,t4,t5,t6,t7,t8,t10,t11,t15,t17,t18,t19,t20,t22,   &
   t23,t25,t29,t31,t33,t34,t37,t39,t40,t41,t45,t46,t48,t49,t50,t51,t52,   &
   t53,t56,t57,t58,t59,t61,t63,t66,t67,t72,t78,t82,t84,t87,t88,t92,t96,   &
   t97,t101,t102,t104,t107,t108,t110,t113,t117,t123,t126,t127,t134,t141,  &
   t150,d0=0.0_iwp,d1=1.0_iwp,d2=2.0_iwp,d4=4.0_iwp,d5=5.0_iwp,d6=6.0_iwp,&
   d7=7.0_iwp,d8=8.0_iwp,d9=9.0_iwp,d12=12.0_iwp,d16=16.0_iwp,            &
   d17=17.0_iwp,d18=18.0_iwp,d24=24.0_iwp,d72=72.0_iwp
 nod=UBOUND(coord,1)
 SELECT CASE(nod)
 CASE(4)
   aa=coord(3,1)-coord(2,1)
   bb=coord(2,2)-coord(1,2)
   t2=-d1+d2*v
   t3=e*t2
   t4=aa**2
   t6=-d1+v
   t7=d2*e*t6
   t8=bb**2
   t10=t3*t4+t7*t8
   t11=d1/aa
   t15=d1/(d1+v)
   t17=d1/t2
   t18=d1/bb*t15*t17
   t20=t10*t11*t18/d6
   t23=e*t15*t17/d8
   t25=-e*t2*t4
   t31=(t25+e*t6*t8)*t11*t18/d6
   t37=e*(d4*v-d1)*t15*t17/d8
   t40=-t10*t11*t18/d12
   t46=(-t25-d4*e*t6*t8)*t11*t18/d12
   t48=t3*t8
   t49=t7*t4+t48
   t52=t49*t11*t18/d6
   t58=(-d4*e*t6*t4+t48)*t11*t18/d12
   t61=-t49*t11*t18/d12
   t67=(e*t6*t4-t48)*t11*t18/d6
   km(1,1)=t20
   km(1,2)=-t23
   km(1,3)=t31
   km(1,4)=t37
   km(1,5)=t40
   km(1,6)=t23
   km(1,7)=t46
   km(1,8)=-t37
   km(2,2)=t52
   km(2,3)=-t37
   km(2,4)=t58
   km(2,5)=t23
   km(2,6)=t61
   km(2,7)=t37
   km(2,8)=t67
   km(3,3)=t20
   km(3,4)=t23
   km(3,5)=t46
   km(3,6)=t37
   km(3,7)=t40
   km(3,8)=-t23
   km(4,4)=t52
   km(4,5)=-t37
   km(4,6)=t67
   km(4,7)=-t23
   km(4,8)=t61
   km(5,5)=t20
   km(5,6)=-t23
   km(5,7)=t31
   km(5,8)=t37
   km(6,6)=t52
   km(6,7)=-t37
   km(6,8)=t58
   km(7,7)=t20
   km(7,8)=t23
   km(8,8)=t52
 CASE(8)
   aa=coord(5,1)-coord(3,1)
   bb=coord(3,2)-coord(1,2)
   t2=-d1+d2*v
   t3=e*t2
   t4=aa**2
   t5=t3*t4
   t6=-d1+v
   t7=d2*e*t6
   t8=bb**2
   t11=d1/aa
   t15=d1/(d1+v)
   t17=d1/t2
   t18=d1/bb*t15*t17
   t19=d5*(t5+t7*t8)*t11*t18
   t20=t19/d18
   t22=e*t15*t17
   t23=d17/d72*t22
   t25=-d4+d8*v
   t29=-e*t6*t8
   t33=(-e*t25*t4-t29)*t11*t18/d9
   t34=d12*v
   t37=t15*t17
   t39=e*(t34-d1)*t37/d18
   t40=e*t6
   t41=t40*t8
   t45=(t5+t41)*t11*t18/d6
   t46=d4*v
   t50=e*(t46-d1)*t37/d24
   t51=d8*e*t6
   t53=-t5-t51*t8
   t56=t53*t11*t18/d18
   t57=t22/d18
   t58=t19/36
   t59=d7/d72*t22
   t61=e*(-d2+t46)
   t63=-t61*t4-t41
   t66=t63*t11*t18/d9
   t67=d4*e*t6
   t72=(t5+t67*t8)*t11*t18/d12
   t78=(t5-d16*e*t6*t8)*t11*t18/d18
   t82=e*(t34-d5)*t37/d18
   t84=t3*t8
   t87=d5*(t7*t4+t84)*t11*t18
   t88=t87/d18
   t92=-e*t2*t8
   t96=(-d16*e*t6*t4-t92)*t11*t18/d18
   t97=t67*t4
   t101=(t97+t84)*t11*t18/d12
   t102=t40*t4
   t104=-t102-t61*t8
   t107=t104*t11*t18/d9
   t108=t87/36
   t110=-t51*t4-t84
   t113=t110*t11*t18/d18
   t117=(t102+t84)*t11*t18/d6
   t123=(t102-e*t25*t8)*t11*t18/d9
   t126=-d4/d9*t63*t11*t18
   t127=d2/d9*t22
   t134=-d2/d9*t110*t11*t18
   t141=-d2/d9*t53*t11*t18
   t150=-d4/d9*t104*t11*t18
   km(1,1)=t20
   km(1,2)=-t23
   km(1,3)=t33
   km(1,4)=t39
   km(1,5)=t45
   km(1,6)=-t50
   km(1,7)=t56
   km(1,8)=t57
   km(1,9)=t58
   km(1,10)=-t59
   km(1,11)=t66
   km(1,12)=t57
   km(1,13)=t72
   km(1,14)=t50
   km(1,15)=t78
   km(1,16)=-t82
   km(2,2)=t88
   km(2,3)=-t82
   km(2,4)=t96
   km(2,5)=t50
   km(2,6)=t101
   km(2,7)=t57
   km(2,8)=t107
   km(2,9)=-t59
   km(2,10)=t108
   km(2,11)=t57
   km(2,12)=t113
   km(2,13)=-t50
   km(2,14)=t117
   km(2,15)=t39
   km(2,16)=t123
   km(3,3)=t126
   km(3,4)=d0
   km(3,5)=t33
   km(3,6)=t82
   km(3,7)=d0
   km(3,8)=t127
   km(3,9)=t66
   km(3,10)=t57
   km(3,11)=d4/d9*(t5+t29)*t11*t18
   km(3,12)=d0
   km(3,13)=t66
   km(3,14)=-t57
   km(3,15)=d0
   km(3,16)=-t127
   km(4,4)=t134
   km(4,5)=-t39
   km(4,6)=t96
   km(4,7)=t127
   km(4,8)=d0
   km(4,9)=t57
   km(4,10)=t113
   km(4,11)=d0
   km(4,12)=d2/d9*(t97+t92)*t11*t18
   km(4,13)=-t57
   km(4,14)=t113
   km(4,15)=-t127
   km(4,16)=d0
   km(5,5)=t20
   km(5,6)=t23
   km(5,7)=t78
   km(5,8)=t82
   km(5,9)=t72
   km(5,10)=-t50
   km(5,11)=t66
   km(5,12)=-t57
   km(5,13)=t58
   km(5,14)=t59
   km(5,15)=t56
   km(5,16)=-t57
   km(6,6)=t88
   km(6,7)=-t39
   km(6,8)=t123
   km(6,9)=t50
   km(6,10)=t117
   km(6,11)=-t57
   km(6,12)=t113
   km(6,13)=t59
   km(6,14)=t108
   km(6,15)=-t57
   km(6,16)=t107
   km(7,7)=t141
   km(7,8)=d0
   km(7,9)=t78
   km(7,10)=t39
   km(7,11)=d0
   km(7,12)=-t127
   km(7,13)=t56
   km(7,14)=-t57
   km(7,15)=d2/d9*(-t5+d4*e*t6*t8)*t11*t18
   km(7,16)=d0
   km(8,8)=t150
   km(8,9)=-t82
   km(8,10)=t123
   km(8,11)=-t127
   km(8,12)=d0
   km(8,13)=-t57
   km(8,14)=t107
   km(8,15)=d0
   km(8,16)=d4/d9*(-t102-t92)*t11*t18
   km(9,9)=t20
   km(9,10)=-t23
   km(9,11)=t33
   km(9,12)=t39
   km(9,13)=t45
   km(9,14)=-t50
   km(9,15)=t56
   km(9,16)=t57
   km(10,10)=t88
   km(10,11)=-t82
   km(10,12)=t96
   km(10,13)=t50
   km(10,14)=t101
   km(10,15)=t57
   km(10,16)=t107
   km(11,11)=t126
   km(11,12)=d0
   km(11,13)=t33
   km(11,14)=t82
   km(11,15)=d0
   km(11,16)=t127
   km(12,12)=t134
   km(12,13)=-t39
   km(12,14)=t96
   km(12,15)=t127
   km(12,16)=d0
   km(13,13)=t20
   km(13,14)=t23
   km(13,15)=t78
   km(13,16)=t82
   km(14,14)=t88
   km(14,15)=-t39
   km(14,16)=t123
   km(15,15)=t141
   km(15,16)=d0
   km(16,16)=t150
 CASE DEFAULT
   WRITE(*,*)"Wrong number of nodes for rectangular element"
 END SELECT
 DO i=1,nod*2
   DO j=i+1,nod*2
     km(j,i)=km(i,j)
   END DO
 END DO
RETURN
END SUBROUTINE rect_km
SUBROUTINE rigid_jointed(km,prop,gamma,etype,iel,coord) 
!
! This subroutine forms the stiffness matrix of a
! general beam/column element (1-, 2- or 3-d).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::gamma(:),coord(:,:),prop(:,:)
 INTEGER,INTENT(IN)::etype(:),iel
 REAL(iwp),INTENT(OUT)::km(:,:)
 INTEGER::ndim,i,j,k
 REAL(iwp)::ell,x1,x2,y1,y2,z1,z2,c,s,e1,e2,e3,e4,pi,xl,yl,zl,cg,sg,den,  &
   ea,ei,eiy,eiz,gj,a1,a2,a3,a4,a5,a6,a7,a8,sum,gamrad,x,t(12,12),        &
   tt(12,12),cc(12,12),r0(3,3),zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp,      &
   two=2.0_iwp,d4=4.0_iwp,d6=6.0_iwp,d12=12.0_iwp,d180=180.0_iwp
 ndim=UBOUND(coord,2)
 SELECT CASE(ndim)
 CASE(1)
   ei=prop(1,etype(iel))
   ell=coord(2,1)-coord(1,1)
   km(1,1)=d12*ei/(ell*ell*ell) 
   km(3,3)=km(1,1)
   km(1,2)=d6*ei/(ell*ell) 
   km(2,1)=km(1,2) 
   km(1,4)=km(1,2)
   km(4,1)=km(1,4) 
   km(1,3)=-km(1,1) 
   km(3,1)=km(1,3) 
   km(3,4)=-km(1,2)
   km(4,3)=km(3,4) 
   km(2,3)=km(3,4) 
   km(3,2)=km(2,3)
   km(2,2)=d4*ei/ell
   km(4,4)=km(2,2) 
   km(2,4)=two*ei/ell 
   km(4,2)=km(2,4)
 CASE(2)
   ea=prop(1,etype(iel))
   ei=prop(2,etype(iel))
   x1=coord(1,1)
   y1=coord(1,2)
   x2=coord(2,1)
   y2=coord(2,2)
   ell=SQRT((y2-y1)**2+(x2-x1)**2)
   c=(x2-x1)/ell
   s=(y2-y1)/ell
   e1=ea/ell
   e2=d12*ei/(ell*ell*ell)
   e3=ei/ell
   e4=d6*ei/(ell*ell)
   km(1,1)=c*c*e1+s*s*e2
   km(4,4)=km(1,1)
   km(1,2)=s*c*(e1-e2)
   km(2,1)=km(1,2)
   km(4,5)=km(1,2)
   km(5,4)=km(4,5)
   km(1,3)=-s*e4
   km(3,1)=km(1,3)
   km(1,6)=km(1,3)
   km(6,1)=km(1,6)
   km(3,4)=s*e4 
   km(4,3)=km(3,4)
   km(4,6)=km(3,4)
   km(6,4)=km(4,6)
   km(1,4)=-km(1,1) 
   km(4,1)=km(1,4)
   km(1,5)=s*c*(-e1+e2)
   km(5,1)=km(1,5)
   km(2,4)=km(1,5)
   km(4,2)=km(2,4)
   km(2,2)=s*s*e1+c*c*e2
   km(5,5)=km(2,2)
   km(2,5)=-km(2,2)
   km(5,2)=km(2,5)
   km(2,3)=c*e4
   km(3,2)=km(2,3)
   km(2,6)=km(2,3)
   km(6,2)=km(2,6)
   km(3,3)=d4*e3
   km(6,6)=km(3,3)
   km(3,5)=-c*e4
   km(5,3)=km(3,5)
   km(5,6)=km(3,5)
   km(6,5)=km(5,6)
   km(3,6)=two*e3
   km(6,3)=km(3,6)
 CASE(3)
   ea=prop(1,etype(iel))
   eiy=prop(2,etype(iel))
   eiz=prop(3,etype(iel))
   gj=prop(4,etype(iel))
   x1=coord(1,1)
   y1=coord(1,2)
   z1=coord(1,3)
   x2=coord(2,1)
   y2=coord(2,2)
   z2=coord(2,3)
   xl=x2-x1
   yl=y2-y1
   zl=z2-z1
   ell=SQRT(xl*xl+yl*yl+zl*zl)
   km=zero
   t=zero
   tt=zero
   a1=ea/ell
   a2=d12*eiz/(ell*ell*ell)
   a3=d12*eiy/(ell*ell*ell)
   a4=d6*eiz/(ell*ell)
   a5=d6*eiy/(ell*ell)
   a6=d4*eiz/ell
   a7=d4*eiy/ell
   a8=gj/ell
   km(1,1)=a1
   km(7,7)=a1
   km(1,7)=-a1
   km(7,1)=-a1
   km(2,2)=a2
   km(8,8)=a2
   km(2,8)=-a2
   km(8,2)=-a2
   km(3,3)=a3
   km(9,9)=a3
   km(3,9)=-a3
   km(9,3)=-a3
   km(4,4)=a8
   km(10,10)=a8
   km(4,10)=-a8
   km(10,4)=-a8
   km(5,5)=a7
   km(11,11)=a7
   km(5,11)=pt5*a7
   km(11,5)=pt5*a7
   km(6,6)=a6
   km(12,12)=a6
   km(6,12)=pt5*a6
   km(12,6)=pt5*a6
   km(2,6)=a4
   km(6,2)=a4
   km(2,12)=a4
   km(12,2)=a4
   km(6,8)=-a4
   km(8,6)=-a4
   km(8,12)=-a4
   km(12,8)=-a4
   km(5,9)=a5
   km(9,5)=a5
   km(9,11)=a5
   km(11,9)=a5
   km(3,5)=-a5
   km(5,3)=-a5
   km(3,11)=-a5
   km(11,3)=-a5
   pi=ACOS(-one)
   gamrad=gamma(iel)*pi/d180
   cg=COS(gamrad)
   sg=SIN(gamrad)
   den=ell*SQRT(xl*xl+zl*zl)
   IF(den/=zero)THEN
     r0(1,1)=xl/ell
     r0(1,2)=yl/ell
     r0(1,3)=zl/ell
     r0(2,1)=(-xl*yl*cg-ell*zl*sg)/den
     r0(2,2)=den*cg/(ell*ell)
     r0(2,3)=(-yl*zl*cg+ell*xl*sg)/den
     r0(3,1)=(xl*yl*sg-ell*zl*cg)/den
     r0(3,2)=-den*sg/(ell*ell)
     r0(3,3)=(yl*zl*sg+ell*xl*cg)/den
   ELSE
     r0(1,1)=zero
     r0(1,3)=zero
     r0(2,2)=zero
     r0(3,2)=zero
     r0(1,2)=one
     r0(2,1)=-cg
     r0(3,3)=cg
     r0(2,3)=sg
     r0(3,1)=sg
   END IF
     DO i=1,3
       DO j=1,3 
       x=r0(i,j)
       DO k=0,9,3
         t(i+k,j+k)=x
         tt(j+k,i+k)=x
       END DO
     END DO
   END DO
   DO i=1,12
     DO j=1,12
       sum=zero
       DO k=1,12
         sum=sum+km(i,k)*t(k,j)
       END DO
       cc(i,j)=sum
     END DO
   END DO
   DO i=1,12
     DO j=1,12
       sum=zero
       DO k=1,12
         sum=sum+tt(i,k)*cc(k,j)
       END DO
       km(i,j)=sum
     END DO
   END DO
 END SELECT
RETURN
END SUBROUTINE rigid_jointed                
SUBROUTINE rod_km(km,ea,length)
!
! This subroutine forms the stiffness matrix of a 1-d "rod" element.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::ea,length
 REAL(iwp),INTENT(OUT)::km(:,:)
 REAL(iwp)::one=1.0_iwp
 km(1,1)=one
 km(2,2)=one
 km(1,2)=-one
 km(2,1)=-one
 km=km*ea/length
RETURN
END SUBROUTINE rod_km
SUBROUTINE rod_mm(mm,length)
!
! This subroutine forms the consistent mass matrix of a 1-d "rod" element.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),intent(in)::length
 REAL(iwp),intent(out)::mm(:,:)
 REAL(iwp)::one=1.0_iwp,d3=3.0_iwp,d6=6.0_iwp
 mm(1,1)=one/d3
 mm(1,2)=one/d6
 mm(2,1)=one/d6
 mm(2,2)=one/d3
 mm=mm*length
RETURN
END SUBROUTINE rod_mm
SUBROUTINE sample(element,s,wt)
!
! This subroutine returns the local coordinates and weighting coefficients
! of the integrating points.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(OUT)::s(:,:)
 REAL(iwp),INTENT(OUT),OPTIONAL::wt(:)
 CHARACTER(*),INTENT(IN)::element
 INTEGER::nip
 REAL(iwp)::root3,r15,w(3),v(9),b,c
 root3=1.0_iwp/SQRT(3.0_iwp)
 r15=0.2_iwp*SQRT(15.0_iwp)
 nip=UBOUND(s,1)
 w=(/5.0_iwp/9.0_iwp,8.0_iwp/9.0_iwp,5.0_iwp/9.0_iwp/)
 v=(/5.0_iwp/9.0_iwp*w,8.0_iwp/9.0_iwp*w,5.0_iwp/9.0_iwp*w/)
 SELECT CASE(element)
 CASE('line')
   SELECT CASE(nip)
   CASE(1)
     s(1,1)=0.0_iwp
     wt(1) =2.0_iwp
   CASE(2)
     s(1,1)=-0.577350269189626_iwp
     s(2,1)= 0.577350269189626_iwp
     wt(1) = 1.000000000000000_iwp
     wt(2) = 1.000000000000000_iwp
   CASE(3)
     s(1,1)=-0.774596669241484_iwp
     s(2,1)= 0.000000000000000_iwp
     s(3,1)= 0.774596669241484_iwp
     wt(1) = 0.555555555555556_iwp
     wt(2) = 0.888888888888889_iwp
     wt(3) = 0.555555555555556_iwp
   CASE(4)
     s(1,1)=-0.861136311594053_iwp
     s(2,1)=-0.339981043584856_iwp
     s(3,1)= 0.339981043584856_iwp
     s(4,1)= 0.861136311594053_iwp
     wt(1) = 0.347854845137454_iwp
     wt(2) = 0.652145154862546_iwp
     wt(3) = 0.652145154862546_iwp
     wt(4) = 0.347854845137454_iwp
   CASE(5)
     s(1,1)=-0.906179845938664_iwp
     s(2,1)=-0.538469310105683_iwp
     s(3,1)= 0.000000000000000_iwp
     s(4,1)= 0.538469310105683_iwp
     s(5,1)= 0.906179845938664_iwp
     wt(1) = 0.236926885056189_iwp
     wt(2) = 0.478628670499366_iwp
     wt(3) = 0.568888888888889_iwp
     wt(4) = 0.478628670499366_iwp
     wt(5) = 0.236926885056189_iwp
   CASE(6)
     s(1,1)=-0.932469514203152_iwp
     s(2,1)=-0.661209386466265_iwp
     s(3,1)=-0.238619186083197_iwp
     s(4,1)= 0.238619186083197_iwp
     s(5,1)= 0.661209386466265_iwp
     s(6,1)= 0.932469514203152_iwp
     wt(1) = 0.171324492379170_iwp
     wt(2) = 0.360761573048139_iwp
     wt(3) = 0.467913934572691_iwp
     wt(4) = 0.467913934572691_iwp
     wt(5) = 0.360761573048139_iwp
     wt(6) = 0.171324492379170_iwp
   CASE(7)
     s(1,1)=-0.9491079123427585245261897_iwp
     s(2,1)=-0.7415311855993944398638648_iwp
     s(3,1)=-0.4058451513773971669066064_iwp
     s(4,1)= 0.000000000000000_iwp
     s(5,1)= 0.4058451513773971669066064_iwp
     s(6,1)= 0.7415311855993944398638648_iwp
     s(7,1)= 0.9491079123427585245261897_iwp
     wt(1) = 0.1294849661688696932706114_iwp
     wt(2) = 0.2797053914892766679014678_iwp
     wt(3) = 0.3818300505051189449503698_iwp
     wt(4) = 0.4179591836734693877551020_iwp
     wt(5) = 0.3818300505051189449503698_iwp
     wt(6) = 0.2797053914892766679014678_iwp
     wt(7) = 0.1294849661688696932706114_iwp
   CASE(8)
     s(1,1)=-0.9602898564975362316835609_iwp
     s(2,1)=-0.7966664774136267395915539_iwp
     s(3,1)=-0.5255324099163289858177390_iwp
     s(4,1)=-0.1834346424956498049394761_iwp
     s(5,1)= 0.1834346424956498049394761_iwp
     s(6,1)= 0.5255324099163289858177390_iwp
     s(7,1)= 0.7966664774136267395915539_iwp
     s(8,1)= 0.9602898564975362316835609_iwp
     wt(1) = 0.1012285362903762591525314_iwp
     wt(2) = 0.2223810344533744705443560_iwp
     wt(3) = 0.3137066458778872873379622_iwp
     wt(4) = 0.3626837833783619829651504_iwp
     wt(5) = 0.3626837833783619829651504_iwp
     wt(6) = 0.3137066458778872873379622_iwp
     wt(7) = 0.2223810344533744705443560_iwp
     wt(8) = 0.1012285362903762591525314_iwp
   CASE(9)
     s(1,1)=-0.9681602395076260898355762_iwp
     s(2,1)=-0.8360311073266357942994298_iwp    
     s(3,1)=-0.6133714327005903973087020_iwp
     s(4,1)=-0.3242534234038089290385380_iwp    
     s(5,1)= 0.000000000000000_iwp                            
     s(6,1)= 0.3242534234038089290385380_iwp                            
     s(7,1)= 0.6133714327005903973087020_iwp                            
     s(8,1)= 0.8360311073266357942994298_iwp                            
     s(9,1)= 0.9681602395076260898355762_iwp                            
     wt(1) = 0.0812743883615744119718922_iwp                            
     wt(2) = 0.1806481606948574040584720_iwp                            
     wt(3) = 0.2606106964029354623187429_iwp                            
     wt(4) = 0.3123470770400028400686304_iwp                            
     wt(5) = 0.3302393550012597631645251_iwp                            
     wt(6) = 0.3123470770400028400686304_iwp                            
     wt(7) = 0.2606106964029354623187429_iwp                            
     wt(8) = 0.1806481606948574040584720_iwp                            
     wt(9) = 0.0812743883615744119718922_iwp                            
   CASE(10)
     s(1,1)=-0.9739065285171717200779640_iwp            
     s(2,1)=-0.8650633666889845107320967_iwp 
     s(3,1)=-0.6794095682990244062343274_iwp 
     s(4,1)=-0.4333953941292471907992659_iwp 
     s(5,1)=-0.1488743389816312108848260_iwp 
     s(6,1)= 0.1488743389816312108848260_iwp 
     s(7,1)= 0.4333953941292471907992659_iwp 
     s(8,1)= 0.6794095682990244062343274_iwp 
     s(9,1)= 0.8650633666889845107320967_iwp 
    s(10,1)= 0.9739065285171717200779640_iwp 
     wt(1) = 0.0666713443086881375935688_iwp                     
     wt(2) = 0.1494513491505805931457763_iwp                     
     wt(3) = 0.2190863625159820439955349_iwp                     
     wt(4) = 0.2692667193099963550912269_iwp                     
     wt(5) = 0.2955242247147528701738930_iwp                     
     wt(6) = 0.2955242247147528701738930_iwp                      
     wt(7) = 0.2692667193099963550912269_iwp                     
     wt(8) = 0.2190863625159820439955349_iwp                     
     wt(9) = 0.1494513491505805931457763_iwp                     
    wt(10) = 0.0666713443086881375935688_iwp                     
   CASE DEFAULT                              
     WRITE(*,*)"Wrong number of integrating points for a line"
   END SELECT
 CASE('triangle')
 SELECT CASE(nip)
   CASE(1)
     s(1,1)= 0.333333333333333_iwp
     s(1,2)= 0.333333333333333_iwp
     wt(1) = 0.500000000000000_iwp
   CASE(3)
     s(1,1)= 0.500000000000000_iwp
     s(1,2)= 0.500000000000000_iwp
     s(2,1)= 0.500000000000000_iwp
     s(2,2)= 0.000000000000000_iwp
     s(3,1)= 0.000000000000000_iwp
     s(3,2)= 0.500000000000000_iwp
     wt(1:3)=0.333333333333333_iwp
     wt=0.5_iwp*wt
   CASE(4)
     s(1,1)= 0.6_iwp
     s(1,2)= 0.2_iwp
     s(2,1)= 0.2_iwp
     s(2,2)= 0.6_iwp
     s(3,1)= 0.2_iwp
     s(3,2)= 0.2_iwp
     s(4,1)= 0.333333333333333_iwp
     s(4,2)= 0.333333333333333_iwp
     wt(1:3)= 0.520833333333333_iwp
     wt(4)=  -0.5625_iwp
     wt=0.5_iwp*wt
   CASE(6)
     s(1,1)= 0.816847572980459_iwp
     s(1,2)= 0.091576213509771_iwp
     s(2,1)= 0.091576213509771_iwp
     s(2,2)= 0.816847572980459_iwp
     s(3,1)= 0.091576213509771_iwp
     s(3,2)= 0.091576213509771_iwp
     s(4,1)= 0.108103018168070_iwp
     s(4,2)= 0.445948490915965_iwp
     s(5,1)= 0.445948490915965_iwp
     s(5,2)= 0.108103018168070_iwp
     s(6,1)= 0.445948490915965_iwp
     s(6,2)= 0.445948490915965_iwp
     wt(1:3)=0.109951743655322_iwp
     wt(4:6)=0.223381589678011_iwp
     wt=0.5_iwp*wt
   CASE(7)
     s(1,1)= 0.333333333333333_iwp
     s(1,2)= 0.333333333333333_iwp
     s(2,1)= 0.797426985353087_iwp
     s(2,2)= 0.101286507323456_iwp
     s(3,1)= 0.101286507323456_iwp
     s(3,2)= 0.797426985353087_iwp
     s(4,1)= 0.101286507323456_iwp
     s(4,2)= 0.101286507323456_iwp
     s(5,1)= 0.470142064105115_iwp
     s(5,2)= 0.059715871789770_iwp
     s(6,1)= 0.059715871789770_iwp
     s(6,2)= 0.470142064105115_iwp
     s(7,1)= 0.470142064105115_iwp
     s(7,2)= 0.470142064105115_iwp
     wt(1) = 0.225000000000000_iwp
     wt(2:4)=0.125939180544827_iwp
     wt(5:7)=0.132394152788506_iwp
     wt=0.5_iwp*wt
   CASE(12)
     s(1,1)= 0.873821971016996_iwp
     s(1,2)= 0.063089014491502_iwp
     s(2,1)= 0.063089014491502_iwp
     s(2,2)= 0.873821971016996_iwp
     s(3,1)= 0.063089014491502_iwp
     s(3,2)= 0.063089014491502_iwp
     s(4,1)= 0.501426509658179_iwp
     s(4,2)= 0.249286745170910_iwp
     s(5,1)= 0.249286745170910_iwp
     s(5,2)= 0.501426509658179_iwp
     s(6,1)= 0.249286745170910_iwp
     s(6,2)= 0.249286745170910_iwp
     s(7,1) =0.053145049844817_iwp
     s(7,2) =0.310352451033784_iwp
     s(8,1) =0.310352451033784_iwp
     s(8,2) =0.053145049844817_iwp
     s(9,1) =0.053145049844817_iwp
     s(9,2) =0.636502499121398_iwp
     s(10,1)=0.310352451033784_iwp
     s(10,2)=0.636502499121398_iwp
     s(11,1)=0.636502499121398_iwp
     s(11,2)=0.053145049844817_iwp
     s(12,1)=0.636502499121398_iwp
     s(12,2)=0.310352451033784_iwp
     wt(1:3)=0.050844906370207_iwp
     wt(4:6)=0.116786275726379_iwp
     wt(7:12)=0.082851075618374_iwp
     wt=0.5_iwp*wt
   CASE(16)
     s(1,1)=0.333333333333333_iwp
     s(1,2)=0.333333333333333_iwp
     s(2,1)=0.658861384496478_iwp
     s(2,2)=0.170569307751761_iwp
     s(3,1)=0.170569307751761_iwp
     s(3,2)=0.658861384496478_iwp
     s(4,1)=0.170569307751761_iwp
     s(4,2)=0.170569307751761_iwp
     s(5,1)=0.898905543365938_iwp
     s(5,2)=0.050547228317031_iwp
     s(6,1)=0.050547228317031_iwp
     s(6,2)=0.898905543365938_iwp
     s(7,1)=0.050547228317031_iwp
     s(7,2)=0.050547228317031_iwp
     s(8,1)=0.081414823414554_iwp
     s(8,2)=0.459292588292723_iwp
     s(9,1)=0.459292588292723_iwp
     s(9,2)=0.081414823414554_iwp
     s(10,1)=0.459292588292723_iwp
     s(10,2)=0.459292588292723_iwp
     s(11,1)=0.008394777409958_iwp
     s(11,2)=0.263112829634638_iwp
     s(12,1)=0.008394777409958_iwp
     s(12,2)=0.728492392955404_iwp
     s(13,1)=0.263112829634638_iwp
     s(13,2)=0.008394777409958_iwp
     s(14,1)=0.263112829634638_iwp
     s(14,2)=0.728492392955404_iwp
     s(15,1)=0.728492392955404_iwp
     s(15,2)=0.008394777409958_iwp
     s(16,1)=0.728492392955404_iwp
     s(16,2)=0.263112829634638_iwp
     wt(1)=0.144315607677787_iwp
     wt(2:4)=0.103217370534718_iwp
     wt(5:7)=0.032458497623198_iwp
     wt(8:10)=0.095091634267284_iwp
     wt(11:16)=0.027230314174435_iwp
     wt=0.5_iwp*wt
   CASE DEFAULT
     WRITE(*,*)"wrong number of integrating points for a triangle"
   END SELECT
 CASE('quadrilateral')
   SELECT CASE(nip)
   CASE(1)
     s(1,1)=0.0_iwp
     s(1,2)=0.0_iwp
     wt(1)=4.0_iwp
   CASE(4)
     s(1,1)=-root3
     s(1,2)= root3
     s(2,1)= root3
     s(2,2)= root3
     s(3,1)=-root3
     s(3,2)=-root3
     s(4,1)= root3
     s(4,2)=-root3
     wt=1.0_iwp
   CASE(9)
     s(1:7:3,1)=-r15
     s(2:8:3,1)=0.0_iwp
     s(3:9:3,1)=r15
     s(1:3,2)  =r15
     s(4:6,2)  =0.0_iwp
     s(7:9,2)  =-r15
     wt= v
   CASE(16)
     s(1:13:4,1)=-0.861136311594053_iwp
     s(2:14:4,1)=-0.339981043584856_iwp
     s(3:15:4,1)= 0.339981043584856_iwp
     s(4:16:4,1)= 0.861136311594053_iwp
     s(1:4,2)   = 0.861136311594053_iwp
     s(5:8,2)   = 0.339981043584856_iwp
     s(9:12,2)  =-0.339981043584856_iwp
     s(13:16,2) =-0.861136311594053_iwp
     wt(1)      = 0.121002993285602_iwp
     wt(4)      = wt(1)
     wt(13)     = wt(1)
     wt(16)     = wt(1)
     wt(2)      = 0.226851851851852_iwp
     wt(3)      = wt(2)
     wt(5)      = wt(2)
     wt(8)      = wt(2)
     wt(9)      = wt(2)
     wt(12)     = wt(2)
     wt(14)     = wt(2)
     wt(15)     = wt(2)
     wt(6)      = 0.425293303010694_iwp
     wt(7)      = wt(6)
     wt(10)     = wt(6)
     wt(11)     = wt(6)
   CASE(25)
     s(1:21:5,1)= 0.906179845938664_iwp
     s(2:22:5,1)= 0.538469310105683_iwp
     s(3:23:5,1)= 0.0_iwp
     s(4:24:5,1)=-0.538469310105683_iwp
     s(5:25:5,1)=-0.906179845938664_iwp
     s( 1: 5,2) = 0.906179845938664_iwp
     s( 6:10,2) = 0.538469310105683_iwp
     s(11:15,2) = 0.0_iwp
     s(16:20,2) =-0.538469310105683_iwp
     s(21:25,2) =-0.906179845938664_iwp
     wt(1) =0.056134348862429_iwp
     wt(2) =0.113400000000000_iwp
     wt(3) =0.134785072387521_iwp
     wt(4) =0.113400000000000_iwp
     wt(5) =0.056134348862429_iwp
     wt(6) =0.113400000000000_iwp
     wt(7) =0.229085404223991_iwp
     wt(8) =0.272286532550750_iwp
     wt(9) =0.229085404223991_iwp
     wt(10)=0.113400000000000_iwp
     wt(11)=0.134785072387521_iwp
     wt(12)=0.272286532550750_iwp
     wt(13)=0.323634567901235_iwp
     wt(14)=0.272286532550750_iwp
     wt(15)=0.134785072387521_iwp
     wt(16)=0.113400000000000_iwp
     wt(17)=0.229085404223991_iwp
     wt(18)=0.272286532550750_iwp
     wt(19)=0.229085404223991_iwp
     wt(20)=0.113400000000000_iwp
     wt(21)=0.056134348862429_iwp
     wt(22)=0.113400000000000_iwp
     wt(23)=0.134785072387521_iwp
     wt(24)=0.113400000000000_iwp
     wt(25)=0.056134348862429_iwp
   CASE DEFAULT
     WRITE(*,*)"wrong number of integrating points for a quadrilateral"
   END SELECT
 CASE('tetrahedron')
!                       for tetrahedra weights multiplied by 1/6
   SELECT CASE(nip)
   CASE(1)
     s(1,1)=0.25_iwp
     s(1,2)=0.25_iwp
     s(1,3)=0.25_iwp
     wt(1)=1.0_iwp/6.0_iwp
   CASE(4)
     s(1,1)=0.58541020_iwp
     s(1,2)=0.13819660_iwp
     s(1,3)=s(1,2)
     s(2,2)=s(1,1)
     s(2,3)=s(1,2)
     s(2,1)=s(1,2)
     s(3,3)=s(1,1)
     s(3,1)=s(1,2)
     s(3,2)=s(1,2)
     s(4,1)=s(1,2)
     s(4,2)=s(1,2)
     s(4,3)=s(1,2)
     wt(1:4)=0.25_iwp/6.0_iwp
   CASE(5)
     s(1,1)=0.25_iwp
     s(1,2)=0.25_iwp
     s(1,3)=0.25_iwp
     s(2,1)=0.5_iwp
     s(2,2)=1.0_iwp/6.0_iwp
     s(2,3)=s(2,2)
     s(3,2)=0.5_iwp
     s(3,3)=1.0_iwp/6.0_iwp
     s(3,1)=s(3,3)
     s(4,3)=0.5_iwp
     s(4,1)=1.0_iwp/6.0_iwp
     s(4,2)=s(4,1)
     s(5,1)=1.0_iwp/6.0_iwp
     s(5,2)=s(5,1)
     s(5,3)=s(5,1)
     wt(1)=-0.8_iwp
     wt(2)=9.0_iwp/20.0_iwp
     wt(3:5)=wt(2)
     wt=wt/6.0_iwp
   CASE DEFAULT
     WRITE(*,*)"wrong number of integrating points for a tetrahedron"
   END SELECT
 CASE('hexahedron')
   SELECT CASE(nip)
   CASE(1)
     s(1,1:3)=0.0_iwp
     wt(1)=8.0_iwp
   CASE(8)
     s(1,1)= root3
     s(1,2)= root3
     s(1,3)= root3
     s(2,1)= root3
     s(2,2)= root3
     s(2,3)=-root3
     s(3,1)= root3
     s(3,2)=-root3
     s(3,3)= root3
     s(4,1)= root3
     s(4,2)=-root3
     s(4,3)=-root3
     s(5,1)=-root3
     s(5,2)= root3
     s(5,3)= root3
     s(6,1)=-root3
     s(6,2)=-root3
     s(6,3)= root3
     s(7,1)=-root3
     s(7,2)= root3
     s(7,3)=-root3
     s(8,1)=-root3
     s(8,2)=-root3
     s(8,3)=-root3
     wt=1.0_iwp
   CASE(14)
     b=0.795822426_iwp
     c=0.758786911_iwp
     wt(1:6)=0.886426593_iwp
     wt(7:14)=0.335180055_iwp
     s(1,1)=-b
     s(2,1)=b
     s(3,2)=-b
     s(4,2)=b
     s(5,3)=-b
     s(6,3)=b
     s(7:,:)=c
     s(7,1)=-c
     s(7,2)=-c
     s(7,3)=-c
     s(8,2)=-c
     s(8,3)=-c
     s(9,1)=-c
     s(9,3)=-c
     s(10,3)=-c
     s(11,1)=-c
     s(11,2)=-c
     s(12,2)=-c
     s(13,1)=-c
   CASE(15)
     b=1.0_iwp
     c       =0.674199862_iwp
     wt(1)   =1.564444444_iwp
     wt(2:7) =0.355555556_iwp
     wt(8:15)=0.537777778_iwp
     s(2,1)=-b
     s(3,1)=b
     s(4,2)=-b
     s(5,2)=b
     s(6,3)=-b
     s(7,3)=b
     s(8:,:)=c
     s(8,1)=-c
     s(8,2)=-c
     s(8,3)=-c
     s(9,2)=-c
     s(9,3)=-c
     s(10,1)=-c
     s(10,3)=-c
     s(11,3)=-c
     s(12,1)=-c
     s(12,2)=-c
     s(13,2)=-c
     s(14,1)=-c
   CASE(27)
     wt=(/5.0_iwp/9.0_iwp*v,8.0_iwp/9.0_iwp*v,5.0_iwp/9.0_iwp*v/)
     s(1:7:3,1)=-r15
     s(2:8:3,1)=0.0_iwp
     s(3:9:3,1)=r15
     s(1:3,3)=r15
     s(4:6,3)=0.0_iwp
     s(7:9,3)=-r15
     s(1:9,2)=-r15
     s(10:16:3,1)=-r15
     s(11:17:3,1)=0.0_iwp
     s(12:18:3,1)=r15
     s(10:12,3)=r15
     s(13:15,3)=0.0_iwp
     s(16:18,3)=-r15
     s(10:18,2)=0.0_iwp
     s(19:25:3,1)=-r15
     s(20:26:3,1)=0.0_iwp
     s(21:27:3,1)=r15
     s(19:21,3)= r15
     s(22:24,3)=0.0_iwp
     s(25:27,3)=-r15
     s(19:27,2)= r15
   CASE DEFAULT
     WRITE(*,*)"wrong number of integrating points for a hexahedron"
   END SELECT
 CASE DEFAULT
   WRITE(*,*)"not a valid element type"
 END SELECT
RETURN
END SUBROUTINE sample
SUBROUTINE seep4(kp,coord,perm)
!
! This contains the "analytical" conductivity matrix
! for a four node quadrilateral based on nip=4.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::coord(:,:),perm(:,:)
 REAL(iwp),INTENT(OUT)::kp(:,:)
 INTEGER::i,j
 REAL(iwp)::x1,x2,x3,x4,y1,y2,y3,y4,y31,y42,x31,x42,x34,y34,a2,f1,f2,kx,  &
   ky,ht(4),it(4),temp,s1,s2,t1,t2,t3,gdiag,goppo,two=2.0_iwp,d3=3.0_iwp
! 
!                        nodal coordinates
! 
 x1=coord(1,1)
 x2=coord(2,1)
 x3=coord(3,1)
 x4=coord(4,1)
 y1=coord(1,2)
 y2=coord(2,2)
 y3=coord(3,2)
 y4=coord(4,2)
 kx=perm(1,1)
 ky=perm(2,2)
! 
!                       procedure definitions
! 
 y31=(y3-y1)
 y42=(y4-y2)
 x31=(x3-x1)
 x42=(x4-x2)
 x34=x3-x4
 y34=y3-y4
 ht(1)=x4-x1
 ht(2)=x1-x2
 ht(3)=x2-x3
 it(1)=y4-y1
 it(2)=y1-y2
 it(3)=y2-y3
 a2=y42*x31-y31*x42
 f1=2*y42*x34-2*x42*y34-a2
 f2=2*y31*x34-2*x31*y34-a2
 s1=kx*y42**2+ky*x42**2
 s2=kx*y34**2+ky*x34**2
 t1=kx*y42*y34+ky*x42*x34
 t2=kx*y31*y34+ky*x31*x34
 t3=kx*y42*y31+x31*x42*ky
! 
!                       calculate parent terms
! 
 gdiag=(-(two*a2+f1)*s1)/(two*(d3*a2**2-f1**2))
 goppo=-(s1*(two*a2+f2)+two*t1*(a2+f2)+two*a2*s2)/(two*(d3*a2**2-f2**2))
 kp(1,1)=gdiag+goppo
!
 gdiag=((two*a2+f1)*t3-(a2+f1)*t1)/(two*(d3*a2**2-f1**2))
 goppo=((two*a2+f2)*t3+(a2+f2)*t2)/(two*(d3*a2**2-f2**2))
 kp(2,1)=gdiag+goppo
!
 gdiag=(s1*a2)/(two*(d3*a2**2-f1**2))
 goppo=((two*a2+f2)*s1+(a2+f2)*(two*t1-t3)+two*a2*(s2-t2))/               &
   (two*(d3*a2**2-f2**2))
 kp(3,1)=gdiag+goppo
!
!                       use parent terms and transform
! 
 temp=y31
 y31=y42
 y42=-temp
 temp=x31 
 x31=x42
 x42=-temp
 x34=ht(1)
 y34=it(1)
!
 f1=2*y42*x34-2*x42*y34-a2
 f2=2*y31*x34-2*x31*y34-a2
 s1=kx*y42**2+ky*x42**2
 s2=kx*y34**2+ky*x34**2
 t1=kx*y42*y34+ky*x42*x34
 t2=kx*y31*y34+ky*x31*x34
 t3=kx*y42*y31+x31*x42*ky
!
 gdiag=(-(two*a2+f1)*s1)/(two*(d3*a2**2-f1**2))
 goppo=-(s1*(two*a2+f2)+two*t1*(a2+f2)+two*a2*s2)/(two*(d3*a2**2-f2**2))
 kp(2,2)=gdiag+goppo
!
 gdiag=((two*a2+f1)*t3-(a2+f1)*t1)/(two*(d3*a2**2-f1**2))
 goppo=((two*a2+f2)*t3+(a2+f2)*t2)/(two*(d3*a2**2-f2**2))
 kp(3,2)=gdiag+goppo
!
 gdiag=(s1*a2)/(two*(d3*a2**2-f1**2))
 goppo=((two*a2+f2)*s1+(a2+f2)*(two*t1-t3)+two*a2*(s2-t2))/               &
   (two*(d3*a2**2-f2**2))
 kp(4,2)=gdiag+goppo
!
 temp=y31
 y31=y42
 y42=-temp
 temp=x31
 x31=x42
 x42=-temp
 x34=ht(2)
 y34=it(2)
!
 f1=2*y42*x34-2*x42*y34-a2
 f2=2*y31*x34-2*x31*y34-a2
 s1=kx*y42**2+ky*x42**2
 s2=kx*y34**2+ky*x34**2
 t1=kx*y42*y34+ky*x42*x34
 t2=kx*y31*y34+ky*x31*x34
 t3=kx*y42*y31+x31*x42*ky
!
 gdiag=(-(two*a2+f1)*s1)/(two*(d3*a2**2-f1**2))
 goppo=-(s1*(two*a2+f2)+two*t1*(a2+f2)+two*a2*s2)/(two*(d3*a2**2-f2**2))
 kp(3,3)=gdiag+goppo
!
 gdiag=((two*a2+f1)*t3-(a2+f1)*t1)/(two*(d3*a2**2-f1**2))
 goppo=((two*a2+f2)*t3+(a2+f2)*t2)/(two*(d3*a2**2-f2**2))
 kp(4,3)=gdiag+goppo
!
 temp=y31
 y31=y42
 y42=-temp
 temp=x31
 x31=x42
 x42=-temp
 x34=ht(3)
 y34=it(3)
!
 f1=2*y42*x34-2*x42*y34-a2
 f2=2*y31*x34-2*x31*y34-a2
 s1=kx*y42**2+ky*x42**2
 s2=kx*y34**2+ky*x34**2
 t1=kx*y42*y34+ky*x42*x34
 t2=kx*y31*y34+ky*x31*x34
 t3=kx*y42*y31+x31*x42*ky
!
 gdiag=(-(two*a2+f1)*s1)/(two*(d3*a2**2-f1**2))
 goppo=-(s1*(two*a2+f2)+two*t1*(a2+f2)+two*a2*s2)/(two*(d3*a2**2-f2**2))
 kp(4,4)=gdiag+goppo
!
 gdiag=((two*a2+f1)*t3-(a2+f1)*t1)/(two*(d3*a2**2-f1**2))
 goppo=((two*a2+f2)*t3+(a2+f2)*t2)/(two*(d3*a2**2-f2**2))
 kp(4,1)=gdiag+goppo
! 
!                       mirror matrix about main diagonal
! 
 DO i=2,4
   DO j=1,i-1
     kp(j,i)=kp(i,j)
   END DO
 END DO
RETURN
END SUBROUTINE seep4  
SUBROUTINE shape_der(der,points,i)
!
!   This subroutine produces derivatives of shape functions withe respect
!   to local coordinates.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::i
 REAL(iwp),INTENT(IN)::points(:,:)
 REAL(iwp),INTENT(OUT)::der(:,:)
 REAL(iwp)::eta,xi,zeta,xi0,eta0,zeta0,etam,etap,xim,xip,c1,c2,c3 
 REAL(iwp)::t1,t2,t3,t4,t5,t6,t7,t8,t9,x2p1,x2m1,e2p1,e2m1,zetam,zetap
 REAL,PARAMETER::zero=0.0_iwp,pt125=0.125_iwp,pt25=0.25_iwp,pt5=0.5_iwp,  &
   pt75=0.75_iwp,one=1.0_iwp,two=2.0_iwp,d3=3.0_iwp,d4=4.0_iwp,d5=5.0_iwp,&
   d6=6.0_iwp,d8=8.0_iwp,d9=9.0_iwp,d10=10.0_iwp,d11=11.0_iwp,            &
   d12=12.0_iwp,d16=16.0_iwp,d18=18.0_iwp,d27=27.0_iwp,d32=32.0_iwp,      &
   d36=36.0_iwp,d54=54.0_iwp,d64=64.0_iwp,d128=128.0_iwp
 INTEGER::xii(20),etai(20),zetai(20),l,ndim,nod
 ndim=UBOUND(der,1)
 nod= UBOUND(der,2)
 SELECT CASE(ndim)
 CASE(1)   ! one dimensional elements
   xi=points(i,1)
   SELECT CASE(nod)
   CASE(2)
     der(1,1)=-pt5 
     der(1,2)= pt5
   CASE(3)
     t1=-one-xi 
     t2=-xi  
     t3=one-xi
     der(1,1)=-(t3+t2)/two  
     der(1,2)=(t3+t1)    
     der(1,3)=-(t2+t1)/two   
   CASE(4)
     t1=-one-xi 
     t2=-one/d3-xi 
     t3=one/d3-xi 
     t4=one-xi
     der(1,1)=-(t3*t4+t2*t4+t2*t3)*d9/d16     
     der(1,2)=(t3*t4+t1*t4+t1*t3)*d27/d16 
     der(1,3)=-(t2*t4+t1*t4+t1*t2)*d27/d16 
     der(1,4)=(t2*t3+t1*t3+t1*t2)*d9/d16   
   CASE(5)
     t1=-one-xi 
     t2=-pt5-xi 
     t3=-xi 
     t4=pt5-xi 
     t5=one-xi
     der(1,1)=-(t3*t4*t5+t2*t4*t5+t2*t3*t5+t2*t3*t4)*two/d3   
     der(1,2)=(t3*t4*t5+t1*t4*t5+t1*t3*t5+t1*t3*t4)*d8/d3
     der(1,3)=-(t2*t4*t5+t1*t4*t5+t1*t2*t5+t1*t2*t4)*d4 
     der(1,4)=(t2*t3*t5+t1*t3*t5+t1*t2*t5+t1*t2*t3)*d8/d3
     der(1,5)=-(t2*t3*t4+t1*t3*t4+t1*t2*t4+t1*t2*t3)*two/d3
   CASE DEFAULT
     WRITE(*,*)"wrong number of nodes in shape_der"        
   END SELECT
 CASE(2)      ! two dimensional elements
   xi=points(i,1)
   eta=points(i,2) 
   c1=xi 
   c2=eta 
   c3=one-c1-c2
   etam=pt25*(one-eta)
   etap=pt25*(one+eta)
   xim= pt25*(one-xi)
   xip= pt25*(one+xi)
   x2p1=two*xi+one 
   x2m1=two*xi-one 
   e2p1=two*eta+one 
   e2m1=two*eta-one
   SELECT CASE(nod)
   CASE(3)
     der(1,1)=one
     der(1,3)=zero
     der(1,2)=-one
     der(2,1)=zero
     der(2,3)=one
     der(2,2)=-one
   CASE(6) 
     der(1,1)=d4*c1-one 
     der(1,6)=d4*c2
     der(1,5)=zero  
     der(1,4)=-d4*c2
     der(1,3)=-(d4*c3-one)
     der(1,2)=d4*(c3-c1)
     der(2,1)=zero
     der(2,6)=d4*c1 
     der(2,5)=d4*c2-one
     der(2,4)=d4*(c3-c2)
     der(2,3)=-(d4*c3-one)  
     der(2,2)=-d4*c1
   CASE(10)                          
     der(1,1)=(d27*c1**2-d18*c1+two)/two
     der(1,9)=(d9*(d6*c1-one)*c2)/two
     der(1,8)=(d9*(d3*c2-one)*c2)/two
     der(1,7)=zero
     der(1,6)=-(d9*(d3*c2-one)*c2)/two
     der(1,5)= (d9*(d6*c1+d6*c2-d5)*c2)/two
     der(1,4)=-(d27*c1**2+d54*c1*c2-d36*c1+d27*c2**2-d36*c2+d11)/two
     der(1,3)= (d9*(d9*c1**2+d12*c1*c2-d10*c1+d3*c2**2-d5*c2+two))/two
     der(1,2)=-(d9*(d9*c1**2+d6*c1*c2-d8*c1-c2+one))/two
     der(1,10)=-d27*(((c2-one)+c1)+c1)*c2
     der(2,1)=zero
     der(2,9)= (d9*(d3*c1-one)*c1)/two
     der(2,8)= (d9*(d6*c2-one)*c1)/two
     der(2,7)=(d27*c2**2-d18*c2+two)/two
     der(2,6)=-(d9*((c1+c2-one)*(d6*c2-one)+(d3*c2-one)*c2))/two
     der(2,5)= (d9*(d3*c1**2+d12*c1*c2-d5*c1+d9*c2**2-d10*c2+two))/two
     der(2,4)=-(d27*c1**2+d54*c1*c2-d36*c1+d27*c2**2-d36*c2+d11)/two
     der(2,3)= (d9*(d6*c1+d6*c2-d5)*c1)/two
     der(2,2)=-(d9*(d3*c1-one)*c1)/two
     der(2,10)=-d27*(((c2-one)+c1)+c2)*c1
   CASE(15)                          
     t1=c1-pt25  
     t2=c1-pt5 
     t3=c1-pt75   
     t4=c2-pt25
     t5=c2-pt5   
     t6=c2-pt75 
     t7=c3-pt25  
     t8=c3-pt5 
     t9=c3-pt75
     der(1,1)=d32/d3*(t2*t3*(t1+c1)+c1*t1*(t3+t2))
     der(1,12)=d128/d3*c2*(t2*(t1+c1)+c1*t1) 
     der(1,11)=d64*c2*t4*(t1+c1)
     der(1,10)=d128/d3*c2*t4*t5  
     der(1,9)=zero 
     der(1,8)=-d128/d3*c2*t4*t5
     der(1,7)=-d64*c2*t4*(t7+c3) 
     der(1,6)=-d128/d3*c2*(t8*(t7+c3)+c3*t7)
     der(1,5)=-d32/d3*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
     der(1,4)=d128/d3*(c3*t7*t8-c1*(t8*(t7+c3)+c3*t7))
     der(1,3)=d64*(c3*t7*(t1+c1)-c1*t1*(t7+c3))
     der(1,2)=d128/d3*(c3*(t2*(t1+c1)+c1*t1)-c1*t1*t2)
     der(1,13)=d128*c2*(c3*(t1+c1)-c1*t1) 
     der(1,15)=d128*c2*t4*(c3-c1)
     der(1,14)=d128*c2*(c3*t7-c1*(t7+c3))
     der(2,1)=zero 
     der(2,12)=d128/d3*c1*t1*t2
     der(2,11)=d64*c1*t1*(t4+c2)
     der(2,10)=d128/d3*c1*(t5*(t4+c2)+c2*t4)
     der(2,9)=d32/d3*(t5*t6*(t4+c2)+c2*t4*(t6+t5))
     der(2,8)=d128/d3*((c3*(t5*(t4+c2)+c2*t4))-c2*t4*t5)
     der(2,7)=d64*(c3*t7*(t4+c2)-c2*t4*(t7+c3))
     der(2,6)=d128/d3*(c3*t7*t8-c2*(t8*(t7+c3)+c3*t7))
     der(2,5)=-d32/d3*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
     der(2,4)=-d128/d3*c1*(t8*(t7+c3)+c3*t7)
     der(2,3)=-d64*c1*t1*(t7+c3)  
     der(2,2)=-d128/d3*c1*t1*t2
     der(2,13)=d128*c1*t1*(c3-c2)
     der(2,15)=d128*c1*(c3*(t4+c2)-c2*t4)
     der(2,14)=d128*c1*(c3*t7-c2*(c3+t7))        
   CASE (4)                                                              
     der(1,1)=-etam
     der(1,2)=-etap
     der(1,3)=etap
     der(1,4)=etam
     der(2,1)=-xim
     der(2,2)=xim
     der(2,3)=xip
     der(2,4)=-xip
   CASE(5)
     der(1,1)=-etam+pt5*xi*(one-eta**2)
     der(1,2)=-etap+pt5*xi*(one-eta**2)
     der(1,3)=etap+pt5*xi*(one-eta**2)
     der(1,4)=etam+pt5*xi*(one-eta**2)
     der(1,5)=-two*xi*(one-eta**2)
     der(2,1)=-xim+pt5*eta*(one-xi**2)
     der(2,2)=xim+pt5*eta*(one-xi**2)
     der(2,3)=xip+pt5*eta*(one-xi**2)
     der(2,4)=-xip+pt5*eta*(one-xi**2)
     der(2,5)=-two*eta*(one-xi**2)
   CASE(8)
     der(1,1)=etam*(two*xi+eta)
     der(1,2)=-d8*etam*etap
     der(1,3)=etap*(two*xi-eta)
     der(1,4)=-d4*etap*xi
     der(1,5)=etap*(two*xi+eta)
     der(1,6)=d8*etap*etam
     der(1,7)=etam*(two*xi-eta)
     der(1,8)=-d4*etam*xi
     der(2,1)=xim*(xi+two*eta)
     der(2,2)=-d4*xim*eta
     der(2,3)=xim*(two*eta-xi)
     der(2,4)=d8*xim*xip
     der(2,5)=xip*(xi+two*eta)
     der(2,6)=-d4*xip*eta
     der(2,7)=xip*(two*eta-xi)
     der(2,8)=-d8*xim*xip   
   CASE(9)
     etam=eta-one
     etap=eta+one
     xim=xi-one
     xip=xi+one
     der(1,1)=pt25*x2m1*eta*etam  
     der(1,2)=-pt5*x2m1*etap*etam
     der(1,3)=pt25*x2m1*eta*etap  
     der(1,4)=-xi*eta*etap
     der(1,5)=pt25*x2p1*eta*etap  
     der(1,6)=-pt5*x2p1*etap*etam
     der(1,7)=pt25*x2p1*eta*etam  
     der(1,8)=-xi*eta*etam
     der(1,9)=two*xi*etap*etam    
     der(2,1)=pt25*xi*xim*e2m1
     der(2,2)=-xi*xim*eta        
     der(2,3)=pt25*xi*xim*e2p1
     der(2,4)=-pt5*xip*xim*e2p1   
     der(2,5)=pt25*xi*xip*e2p1
     der(2,6)=-xi*xip*eta        
     der(2,7)=pt25*xi*xip*e2m1
     der(2,8)=-pt5*xip*xim*e2m1   
     der(2,9)=two*xip*xim*eta
   CASE DEFAULT
     WRITE(*,*)"wrong number of nodes in shape_der"        
   END SELECT
 CASE(3)  ! d3 dimensional elements
   xi=points(i,1)
   eta=points(i,2)
   zeta=points(i,3)
   etam=one-eta 
   xim=one-xi
   zetam=one-zeta
   etap=eta+one 
   xip=xi+one 
   zetap=zeta+one
   SELECT CASE(nod)
   CASE(4)
     der(1:3,1:4)=zero
     der(1,1)=one
     der(2,2)=one  
     der(3,3)=one
     der(1,4)=-one 
     der(2,4)=-one 
     der(3,4)=-one  
   CASE(8)
     der(1,1)=-pt125*etam*zetam    
     der(1,2)=-pt125*etam*zetap
     der(1,3)= pt125*etam*zetap     
     der(1,4)= pt125*etam*zetam
     der(1,5)=-pt125*etap*zetam    
     der(1,6)=-pt125*etap*zetap
     der(1,7)= pt125*etap*zetap     
     der(1,8)= pt125*etap*zetam
     der(2,1)=-pt125*xim*zetam     
     der(2,2)=-pt125*xim*zetap
     der(2,3)=-pt125*xip*zetap     
     der(2,4)=-pt125*xip*zetam
     der(2,5)= pt125*xim*zetam      
     der(2,6)= pt125*xim*zetap
     der(2,7)= pt125*xip*zetap      
     der(2,8)= pt125*xip*zetam
     der(3,1)=-pt125*xim*etam      
     der(3,2)= pt125*xim*etam
     der(3,3)= pt125*xip*etam       
     der(3,4)=-pt125*xip*etam
     der(3,5)=-pt125*xim*etap      
     der(3,6)= pt125*xim*etap
     der(3,7)= pt125*xip*etap       
     der(3,8)=-pt125*xip*etap  
   CASE(14) ! type 6 element
     der(1,1)= (two*xi*eta+two*xi*zeta+d4*xi+eta*zeta+eta+zeta)*          &
       (eta-one)*(zeta-one)/d8
     der(1,2)=-(two*xi*eta-two*xi*zeta+d4*xi-eta*zeta+eta-zeta)*          &
       (eta-one)*(zeta+one)/d8
     der(1,3)=-(two*xi*eta-two*xi*zeta+d4*xi+eta*zeta-eta+zeta)*          &
       (eta-one)*(zeta+one)/d8
     der(1,4)= (two*xi*eta+two*xi*zeta+d4*xi-eta*zeta-eta-zeta)*          &
       (eta-one)*(zeta-one)/d8
     der(1,5)= -(eta-one)*(zeta+one)*(zeta-one)*xi 
     der(1,6)=-(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
     der(1,7)=  (eta+one)*(eta-one)*(zeta+one)*xi
     der(1,8)= (eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
     der(1,9)= -(eta+one)*(eta-one)*(zeta-one)*xi  
     der(1,10)= (two*xi*eta-two*xi*zeta-d4*xi+eta*zeta+eta-zeta)*         &
       (eta+one)*(zeta-one)/d8
     der(1,11)=-(two*xi*eta+two*xi*zeta-d4*xi-eta*zeta+eta+zeta)*         &
       (eta+one)*(zeta+one)/d8
     der(1,12)=-(two*xi*eta+two*xi*zeta-d4*xi+eta*zeta-eta-zeta)*         &
       (eta+one)*(zeta+one)/d8
     der(1,13)= (two*xi*eta-two*xi*zeta-d4*xi-eta*zeta-eta+zeta)*         &
       (eta+one)*(zeta-one)/d8
     der(1,14)=  (eta+one)*(zeta+one)*(zeta-one)*xi
     der(2,1)= (two*xi*eta+xi*zeta+xi+two*eta*zeta+d4*eta+zeta)*          &
       (xi-one)*(zeta-one)/d8                                  
     der(2,2)=-(two*xi*eta-xi*zeta+xi-two*eta*zeta+d4*eta-zeta)*          &
       (xi-one)*(zeta+one)/d8
     der(2,3)=-(two*xi*eta-xi*zeta+xi+two*eta*zeta-d4*eta+zeta)*          &
       (xi+one)*(zeta+one)/d8
     der(2,4)= (two*xi*eta+xi*zeta+xi-two*eta*zeta-d4*eta-zeta)*          &
       (xi+one)*(zeta-one)/d8
     der(2,5)=-(xi+one)*(xi-one)*(zeta+one)*(zeta-one)/two
     der(2,6)= -(xi-one)*(zeta+one)*(zeta-one)*eta
     der(2,7)=  (xi+one)*(xi-one)*(zeta+one)*eta
     der(2,8)=  (xi+one)*(zeta+one)*(zeta-one)*eta
     der(2,9)= -(xi+one)*(xi-one)*(zeta-one)*eta
     der(2,10)= (two*xi*eta-xi*zeta-xi+two*eta*zeta+d4*eta-zeta)*         &
       (xi-one)*(zeta-one)/d8
     der(2,11)=-(two*xi*eta+xi*zeta-xi-two*eta*zeta+d4*eta+zeta)*         &
       (xi-one)*(zeta+one)/d8
     der(2,12)=-(two*xi*eta+xi*zeta-xi+two*eta*zeta-d4*eta-zeta)*         &
       (xi+one)*(zeta+one)/d8
     der(2,13)= (two*xi*eta-xi*zeta-xi-two*eta*zeta-d4*eta+zeta)*         &
       (xi+one)*(zeta-one)/d8
     der(2,14)= (xi+one)*(xi-one)*(zeta+one)*(zeta-one)/two
     der(3,1)= (xi*eta+two*xi*zeta+xi+two*eta*zeta+eta+d4*zeta)*          &
       (xi-one)*(eta-one)/d8
     der(3,2)=-(xi*eta-two*xi*zeta+xi-two*eta*zeta+eta-d4*zeta)*          &
       (xi-one)*(eta-one)/d8
     der(3,3)=-(xi*eta-two*xi*zeta+xi+two*eta*zeta-eta+d4*zeta)*          &
       (xi+one)*(eta-one)/d8
     der(3,4)= (xi*eta+two*xi*zeta+xi-two*eta*zeta-eta-d4*zeta)*          &
       (xi+one)*(eta-one)/d8
     der(3,5)= -(xi+one)*(xi-one)*(eta-one)*zeta  
     der(3,6)= -(xi-one)*(eta+one)*(eta-one)*zeta  
     der(3,7)= (xi+one)*(xi-one)*(eta+one)*(eta-one)/two
     der(3,8)=  (xi+one)*(eta+one)*(eta-one)*zeta
     der(3,9)=-(xi+one)*(xi-one)*(eta+one)*(eta-one)/two
     der(3,10)= (xi*eta-two*xi*zeta-xi+two*eta*zeta+eta-d4*zeta)*         &
       (xi-one)*(eta+one)/d8
     der(3,11)=-(xi*eta+two*xi*zeta-xi-two*eta*zeta+eta+d4*zeta)*         &
       (xi-one)*(eta+one)/d8
     der(3,12)=-(xi*eta+two*xi*zeta-xi+two*eta*zeta-eta-d4*zeta)*         &
       (xi+one)*(eta+one)/d8
     der(3,13)= (xi*eta-two*xi*zeta-xi-two*eta*zeta-eta+d4*zeta)*         &
       (xi+one)*(eta+one)/d8
     der(3,14)=  (xi+one)*(xi-one)*(eta+one)*zeta
   CASE(20)
     xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
     etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
     zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
     DO l=1,20
       xi0=xi*xii(l)
       eta0=eta*etai(l)
       zeta0=zeta*zetai(l)
       IF(l==4.OR.l==8.OR.l==16.OR.l==20)THEN
         der(1,l)=-pt5*xi*(one+eta0)*(one+zeta0)
         der(2,l)=pt25*etai(l)*(one-xi*xi)*(one+zeta0)
         der(3,l)=pt25*zetai(l)*(one-xi*xi)*(one+eta0)
       ELSE IF(l>=9.AND.l<=12)THEN
         der(1,l)=pt25*xii(l)*(one-eta*eta)*(one+zeta0)
         der(2,l)=-pt5*eta*(one+xi0)*(one+zeta0)
         der(3,l)=pt25*zetai(l)*(one+xi0)*(one-eta*eta)
       ELSE IF(l==2.OR.l==6.OR.l==14.OR.l==18) THEN
         der(1,l)=pt25*xii(l)*(one+eta0)*(one-zeta*zeta)
         der(2,l)=pt25*etai(l)*(one+xi0)*(one-zeta*zeta)
         der(3,l)=-pt5*zeta*(one+xi0)*(one+eta0)
       ELSE
         der(1,l)=pt125*xii(l)*(one+eta0)*(one+zeta0)*                    &
           (two*xi0+eta0+zeta0-one)
         der(2,l)=pt125*etai(l)*(one+xi0)*(one+zeta0)*                    &
           (xi0+two*eta0+zeta0-one)
         der(3,l)=pt125*zetai(l)*(one+xi0)*(one+eta0)*                    &
           (xi0+eta0+two*zeta0-one)
       END IF
     END DO 
   CASE DEFAULT
     WRITE(*,*)"wrong number of nodes in shape_der"        
   END SELECT
 CASE DEFAULT
   WRITE(*,*)"wrong number of dimensions in shape_der"
 END SELECT
RETURN
END SUBROUTINE shape_der

SUBROUTINE shape_fun(fun,points,i)
!
!   This subroutine computes the values of the shape functions.
!   to local coordinates
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(in)::i
 REAL(iwp),INTENT(IN)::points(:,:)
 REAL(iwp),INTENT(OUT)::fun(:)
 REAL(iwp)::eta,xi,etam,etap,xim,xip,zetam,zetap,c1,c2,c3     
 REAL(iwp)::t1,t2,t3,t4,t5,t6,t7,t8,t9
 REAL(iwp)::zeta,xi0,eta0,zeta0
 INTEGER::xii(20),etai(20),zetai(20),l,ndim,nod
 REAL,PARAMETER::pt125=0.125_iwp,pt25=0.25_iwp,pt5=0.5_iwp,pt75=0.75_iwp, &
   one=1.0_iwp,two=2.0_iwp,d3=3.0_iwp,d4=4.0_iwp,d8=8.0_iwp,d9=9.0_iwp,   &
   d16=16.0_iwp,d27=27.0_iwp,d32=32.0_iwp,d64=64.0_iwp,d128=128.0_iwp
 ndim=UBOUND(points,2)
 nod=UBOUND(fun,1)  
 SELECT CASE(ndim)
 CASE(1) ! one dimensional case
   xi=points(i,1)
   SELECT CASE(nod)
   CASE(2)
     t1=-one-xi 
     t2= one-xi
     fun(1)=t2/two 
     fun(2)=-t1/two
   CASE(3)
     t1=-one-xi 
     t2=-xi 
     t3=one-xi
     fun(1)=t2*t3/two 
     fun(2)=-t1*t3 
     fun(3)=t1*t2/two
   CASE(4)
     t1=-one-xi 
     t2=-one/d3-xi 
     t3=one/d3-xi 
     t4=one-xi
     fun(1)=t2*t3*t4*d9/d16  
     fun(2)=-t1*t3*t4*d27/d16
     fun(3)=t1*t2*t4*d27/d16 
     fun(4)=-t1*t2*t3*d9/d16
   CASE(5)
     t1=-one -xi 
     t2=-pt5-xi 
     t3=-xi 
     t4=pt5-xi 
     t5=one-xi
     fun(1)=t2*t3*t4*t5*two/d3 
     fun(2)=-t1*t3*t4*t5*d8/d3
     fun(3)=t1*t2*t4*t5*d4 
     fun(4)=-t1*t2*t3*t5*d8/d3
     fun(5)=t1*t2*t3*t4*two/d3
   CASE DEFAULT
     WRITE(*,*)"wrong number of nodes in shape_fun"
   END SELECT
 CASE(2) ! two dimensional case
   c1=points(i,1)
   c2=points(i,2)
   c3=one-c1-c2 
   xi=points(i,1)
   eta=points(i,2)
   etam=pt25*(one-eta)
   etap=pt25*(one+eta)
   xim=pt25*(one-xi)
   xip=pt25*(one+xi)
   SELECT CASE(nod)
   CASE(3)
     fun = (/c1,c3,c2/)  
   CASE(6)
     fun(1)=(two*c1-one)*c1 
     fun(2)=d4*c3*c1
     fun(3)=(two*c3-one)*c3 
     fun(4)=d4*c2*c3      
     fun(5)=(two*c2-one)*c2
     fun(6)=d4*c1*c2 
   CASE(10)
     fun(1)= ((d3*c1-one)*(d3*c1-two)*c1)/two
     fun(2)= -(d9*(d3*c1-one)*(c1+c2-one)*c1)/two
     fun(3)=  (d9*(d3*c1+d3*c2-two)*(c1+c2-one)*c1)/two
     fun(4)=-((d3*c1+d3*c2-one)*(d3*c1+d3*c2-two)*(c1+c2-one))/two    
     fun(5)=  (d9*(d3*c1+d3*c2-two)*(c1+c2-one)*c2)/two
     fun(6)= -(d9*(c1+c2-one)*(d3*c2-one)*c2)/two
     fun(7)= ((d3*c2-one)*(d3*c2-two)*c2)/two
     fun(8)=  (d9*(d3*c2-one)*c1*c2)/two
     fun(9)=  (d9*(d3*c1-one)*c1*c2)/two
     fun(10)=-d27*((c2-one)+c1)*c1*c2
   CASE(15)
     t1=c1-pt25  
     t2=c1-pt5 
     t3=c1-pt75   
     t4=c2-pt25
     t5=c2-pt5   
     t6=c2-pt75 
     t7=c3-pt25  
     t8=c3-pt5 
     t9=c3-pt75
     fun(1)=d32/d3*c1*t1*t2*t3   
     fun(2)=d128/d3*c3*c1*t1*t2
     fun(3)=d64*c3*c1*t1*t7      
     fun(4)=d128/d3*c3*c1*t7*t8
     fun(5)=d32/d3*c3*t7*t8*t9   
     fun(6)=d128/d3*c2*c3*t7*t8
     fun(7)=d64*c2*c3*t4*t7      
     fun(8)=d128/d3*c2*c3*t4*t5
     fun(9)=d32/d3*c2*t4*t5*t6   
     fun(10)=d128/d3*c1*c2*t4*t5
     fun(11)=d64*c1*c2*t1*t4     
     fun(12)=d128/d3*c1*c2*t1*t2
     fun(13)=d128*c1*c2*t1*c3    
     fun(15)=d128*c1*c2*c3*t4
     fun(14)=d128*c1*c2*c3*t7      
   CASE(4)
     fun=(/d4*xim*etam,d4*xim*etap,d4*xip*etap,d4*xip*etam/)
   CASE(5)
     fun=(/d4*xim*etam-pt25*(one-xi**2)*(one-eta**2),			  &
           d4*xim*etap-pt25*(one-xi**2)*(one-eta**2),			  &
           d4*xip*etap-pt25*(one-xi**2)*(one-eta**2),			  &
           d4*xip*etam-pt25*(one-xi**2)*(one-eta**2),			  &
           (one-xi**2)*(one-eta**2)/)
   CASE(8)
     fun=(/d4*etam*xim*(-xi-eta-one),d32*etam*xim*etap,                   &
           d4*etap*xim*(-xi+eta-one),d32*xim*xip*etap,                    &
           d4*etap*xip*(xi+eta-one), d32*etap*xip*etam,                   &
           d4*xip*etam*(xi-eta-one), d32*xim*xip*etam/)
   CASE(9)
     etam=eta-one
     etap=eta+one
     xim=xi-one
     xip=xi+one
     fun=(/pt25*xi*xim*eta*etam,-pt5*xi*xim*etap*etam,                    &
           pt25*xi*xim*eta*etap,-pt5*xip*xim*eta*etap,                    &
           pt25*xi*xip*eta*etap,-pt5*xi*xip*etap*etam,                    &
           pt25*xi*xip*eta*etam,-pt5*xip*xim*eta*etam,                    &
           xip*xim*etap*etam/)
   CASE DEFAULT
     WRITE(*,*)"wrong number of nodes in shape_fun"
   END SELECT
 CASE(3) ! d3 dimensional case
   xi=points(i,1)
   eta=points(i,2)
   zeta=points(i,3)
   etam=one-eta 
   xim=one-xi  
   zetam=one-zeta
   etap=eta+one 
   xip=xi+one   
   zetap=zeta+one
   SELECT CASE(nod)
   CASE(4)
     fun(1)=xi   
     fun(2)=eta 
     fun(3)=zeta 
     fun(4)=one-fun(1)-fun(2)-fun(3)
   CASE(8)
     fun=(/pt125*xim*etam*zetam,pt125*xim*etam*zetap,                     &
           pt125*xip*etam*zetap,pt125*xip*etam*zetam,                     &
           pt125*xim*etap*zetam,pt125*xim*etap*zetap,                     &
           pt125*xip*etap*zetap,pt125*xip*etap*zetam/)
   CASE(14) ! type 6 element
     fun(1) = (xi*eta+xi*zeta+two*xi+eta*zeta+two*eta+two*zeta+two)*      &
       (xi-one)*(eta-one)*(zeta-one)/d8
     fun(2) =-(xi*eta-xi*zeta+two*xi-eta*zeta+two*eta-two*zeta+two)*      &
       (xi-one)*(eta-one)*(zeta+one)/d8
     fun(3) =-(xi*eta-xi*zeta+two*xi+eta*zeta-two*eta+two*zeta-two)*      &
       (xi+one)*(eta-one)*(zeta+one)/d8
     fun(4) = (xi*eta+xi*zeta+two*xi-eta*zeta-two*eta-two*zeta-two)*      &
       (xi+one)*(eta-one)*(zeta-one)/d8
     fun(5) =-(xi+one)*(xi-one)*(eta-one)*(zeta+one)*(zeta-one)/two
     fun(6) =-(xi-one)*(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
     fun(7) = (xi+one)*(xi-one)*(eta+one)*(eta-one)*(zeta+one)/two
     fun(8) = (xi+one)*(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
     fun(9) =-(xi+one)*(xi-one)*(eta+one)*(eta-one)*(zeta-one)/two
     fun(10)= (xi*eta-xi*zeta-two*xi+eta*zeta+two*eta-two*zeta-two)*      &
       (xi-one)*(eta+one)*(zeta-one)/d8
     fun(11)=-(xi*eta+xi*zeta-two*xi-eta*zeta+two*eta+two*zeta-two)*      &
       (xi-one)*(eta+one)*(zeta+one)/d8
     fun(12)=-(xi*eta+xi*zeta-two*xi+eta*zeta-two*eta-two*zeta+two)*      &
       (xi+one)*(eta+one)*(zeta+one)/d8
     fun(13)= (xi*eta-xi*zeta-two*xi-eta*zeta-two*eta+two*zeta+two)*      &
       (xi+one)*(eta+one)*(zeta-one)/d8
     fun(14)= (xi+one)*(xi-one)*(eta+one)*(zeta+one)*(zeta-one)/two
   CASE(20)
     xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
     etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
     zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
     DO l=1,20
       xi0=xi*xii(l)
       eta0=eta*etai(l)
       zeta0=zeta*zetai(l)
       IF(l==4.OR.l==8.OR.l==16.OR.l==20)THEN
         fun(l)=pt25*(one-xi*xi)*(one+eta0)*(one+zeta0)
       ELSE IF(l>=9.AND.l<=12)THEN
         fun(l)=pt25*(one+xi0)*(one-eta*eta)*(one+zeta0)
       ELSE IF(l==2.OR.l==6.OR.l==14.OR.l==18)THEN
         fun(l)=pt25*(one+xi0)*(one+eta0)*(one-zeta*zeta)
       ELSE
         fun(l)=pt125*(one+xi0)*(one+eta0)*(one+zeta0)*(xi0+eta0+zeta0-2)
       END IF
     END DO
   CASE DEFAULT
     WRITE(*,*)"wrong number of nodes in shape_fun"
   END SELECT
 CASE DEFAULT
   WRITE(*,*)"wrong number of dimensions in shape_fun"
 END SELECT
RETURN
END SUBROUTINE shape_fun 
SUBROUTINE solve_band(pb,work,loads)
!
! This subroutine performs Gaussian forward and back-substitution
! on the reduced unsymmetric band matrix pb.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::pb(:,:),work(:,:)
 REAL(iwp),INTENT(OUT)::loads(0:)
 INTEGER::iwp1,n,n1,i,iv,l,iq,iv1
 REAL(iwp)::s,pt5=0.5_iwp
 iwp1=(UBOUND(pb,2)-1)/2+1
 n=UBOUND(pb,1)
 iq=2*iwp1-1  
 n1=n-1
 DO iv=1,n1
   i=INT(work(iwp1,iv)+pt5)
   IF(i/=iv)THEN
     s=loads(iv)
     loads(iv)=loads(i)
     loads(i)=s
   END IF
   l=iv+iwp1-1
   IF(l>n)l=n
   iv1=iv+1
   DO i=iv1,l
     loads(i)=loads(i)-work(i-iv,iv)*loads(iv)
   END DO
 END DO
 loads(n)=loads(n)/pb(n,1)
 iv=n-1
 DO WHILE(iv/=0)
   s=loads(iv)
   l=iq
   IF(iv+l-1>n)l=n-iv+1
   DO i=2,l
     s=s-pb(iv,i)*loads(iv+i-1)
     loads(iv)=s/pb(iv,1)
   END DO
 iv=iv-1
 END DO
RETURN
END SUBROUTINE solve_band
SUBROUTINE spabac(kv,loads,kdiag)
!
! This subroutine performs Cholesky forward and back-substitution
! on a symmetric skyline global matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::kv(:)
 REAL(iwp),INTENT(IN OUT)::loads(0:)
 INTEGER,INTENT(IN)::kdiag(:)
 INTEGER::n,i,ki,l,m,j,it,k
 REAL(iwp)::x
 n=UBOUND(kdiag,1)
 loads(1)=loads(1)/kv(1)
 DO i=2,n
   ki=kdiag(i)-i
   l=kdiag(i-1)-ki+1 
   x=loads(i)
   IF(l/=i)THEN
     m=i-1
     DO j=l,m 
       x=x-kv(ki+j)*loads(j)
     END DO
   END IF
   loads(i)=x/kv(ki+i)
 END DO
 DO it=2,n
   i=n+2-it
   ki=kdiag(i)-i
   x=loads(i)/kv(ki+i)
   loads(i)=x
   l=kdiag(i-1)-ki+1
   IF(l/=i)THEN
     m=i-1
     DO k=l,m
       loads(k)=loads(k)-x*kv(ki+k)
     END DO
   END IF
 END DO
 loads(1)=loads(1)/kv(1)
RETURN
END SUBROUTINE spabac               
SUBROUTINE spabac_gauss(kv,loads,kdiag)
!
! This subroutine performs Gaussian forwrad and back-substitution on a
! skyline matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::kv(:)
 REAL(iwp),INTENT(IN OUT)::loads(0:)
 INTEGER,INTENT(IN)::kdiag(:)
 REAL(iwp)::num,den,fac,asum,zero=0.0_iwp
 INTEGER::i,j,l,n,ii,jj,l1,l2
 n=UBOUND(kdiag,1)
 DO j=1,n-1
   den=kv(kdiag(j))
   ii=0
   DO i=j+1,n
     ii=ii+1
     l=kdiag(i)-ii
     IF(l-kdiag(i-1)>zero)THEN
       num=kv(l)
       fac=num/den
       loads(i)=loads(i)-fac*loads(j)
     END IF
   END DO
 END DO
 loads(n)=loads(n)/kv(kdiag(n))
 DO i=n-1,1,-1
   jj=0
   asum=zero
   DO j=i+1,n
     jj=jj+1
     l1=kdiag(i+jj)-jj
     l2=kdiag(i+jj-1)
     IF(l1-l2>zero)asum=asum+kv(l1)*loads(j)
   END DO
   loads(i)=(loads(i)-asum)/kv(kdiag(i))
 END DO
RETURN
END SUBROUTINE spabac_gauss  















SUBROUTINE sparin(kv,kdiag)
!
! This subroutine performs Cholesky factorisation on a symmetric
! skyline global matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN OUT)::kv(:)
 INTEGER,INTENT(IN)::kdiag(:)
 INTEGER::n,i,ki,l,kj,j,ll,m,k
 REAL(iwp)::x
 n=UBOUND(kdiag,1)  
 kv(1)=SQRT(kv(1))
 DO i=2,n
   ki=kdiag(i)-i
   l=kdiag(i-1)-ki+1
   DO j=l,i
     x=kv(ki+j)
     kj=kdiag(j)-j
     IF(j/=1)THEN
       ll=kdiag(j-1)-kj+1
       ll=max(l,ll)
       IF(ll/=j)THEN
         m=j-1
         DO k=ll,m 
           x=x-kv(ki+k)*kv(kj+k) 
         END DO
       END IF
     END IF
     kv(ki+j)=x/kv(kj+j)
   END DO
   kv(ki+i)=SQRT(x)
 END DO
RETURN
END SUBROUTINE sparin
SUBROUTINE sparin_gauss(kv,kdiag)
!
! This subroutine performs Gaussian factorisation of a skyline matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::kdiag(:)
 REAL(iwp),INTENT(OUT)::kv(:)
 REAL(iwp)::num,den,fac,zero=0.0_iwp
 INTEGER::n,ii,i,j,k,l,kk,l1,l2,l3
 n=UBOUND(kdiag,1)
 DO j=1,n-1
   den=kv(kdiag(j))
   ii=0                 
   DO i=j+1,n
     ii=ii+1
     l=kdiag(i)-ii
     IF(l-kdiag(i-1)>zero)THEN
       num=kv(l)
       fac=num/den
       kk=-1
       DO k=i,n
         kk=kk+1
         l1=kdiag(i+kk)-kk
         l2=l1-ii
         l3=kdiag(i+kk-1)
         IF(l2-l3>zero)kv(l1)=kv(l1)-fac*kv(l2)
       END DO 
     END IF
   END DO
 END DO
RETURN
END SUBROUTINE sparin_gauss
SUBROUTINE stability(kv,gv,kdiag,tol,limit,iters,evec,eval)
!
! This subroutine computes the smallest eigenvalue in a beam
! stability analysis.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::limit,kdiag(:)
 INTEGER,INTENT(OUT)::iters
 INTEGER::neq
 REAL(iwp),INTENT(IN OUT)::kv(:),gv(:),tol,eval
 REAL(iwp),INTENT(OUT)::evec(:)
 REAL(iwp)::big,zero=0.0_iwp,one=1.0_iwp
 LOGICAL::converged
 REAL(iwp),ALLOCATABLE::x0(:),x1(:)
 neq=UBOUND(kdiag,1)
 ALLOCATE(x0(0:neq),x1(0:neq))
 CALL sparin(kv,kdiag)
 iters=0
 x0=zero
 x0(1)=1.0_iwp
 DO
   iters=iters+1
   CALL linmul_sky(gv,x0,x1,kdiag)
   CALL spabac(kv,x1,kdiag)  
   big=MAXVAL(x1(1:))
   IF(ABS(MINVAL(x1(1:)))>big)big=MINVAL(x1(1:))
   x1=x1/big
   converged=(MAXVAL(ABS(x1(1:)-x0(1:)))/MAXVAL(ABS(x1(1:)))<tol)
   x0=x1
   IF(converged.OR.iters==limit)EXIT
 END DO
 x1(1:)=x1(1:)/SQRT(SUM(x1(1:)**2))
 evec=x1
 eval=one/big
RETURN
CONTAINS
SUBROUTINE sparin(kv,kdiag)
!
! This subroutine performs Cholesky factorisation on a symmetric
! skyline global matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN OUT)::kv(:)
 INTEGER,INTENT(IN)::kdiag(:)
 INTEGER::n,i,ki,l,kj,j,ll,m,k
 REAL(iwp)::x
 n=UBOUND(kdiag,1)  
 kv(1)=SQRT(kv(1))
 DO i=2,n
   ki=kdiag(i)-i
   l=kdiag(i-1)-ki+1
   DO j=l,i
     x=kv(ki+j)
     kj=kdiag(j)-j
     IF(j/=1)THEN
       ll=kdiag(j-1)-kj+1
       ll=max(l,ll)
       IF(ll/=j)THEN
         m=j-1
         DO k=ll,m 
           x=x-kv(ki+k)*kv(kj+k) 
         END DO
       END IF
     END IF
     kv(ki+j)=x/kv(kj+j)
   END DO
   kv(ki+i)=SQRT(x)
 END DO
RETURN
END SUBROUTINE sparin
SUBROUTINE spabac(kv,loads,kdiag)
!
! This subroutine performs Cholesky forward and back-substitution
! on a symmetric skyline global matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::kv(:)
 REAL(iwp),INTENT(IN OUT)::loads(0:)
 INTEGER,INTENT(IN)::kdiag(:)
 INTEGER::n,i,ki,l,m,j,it,k
 REAL(iwp)::x
 n=UBOUND(kdiag,1)
 loads(1)=loads(1)/kv(1)
 DO i=2,n
   ki=kdiag(i)-i
   l=kdiag(i-1)-ki+1 
   x=loads(i)
   IF(l/=i)THEN
     m=i-1
     DO j=l,m 
       x=x-kv(ki+j)*loads(j)
     END DO
   END IF
   loads(i)=x/kv(ki+i)
 END DO
 DO it=2,n
   i=n+2-it
   ki=kdiag(i)-i
   x=loads(i)/kv(ki+i)
   loads(i)=x
   l=kdiag(i-1)-ki+1
   IF(l/=i)THEN
     m=i-1
     DO k=l,m
       loads(k)=loads(k)-x*kv(ki+k)
     END DO
   END IF
 END DO
 loads(1)=loads(1)/kv(1)
RETURN
END SUBROUTINE spabac               
SUBROUTINE linmul_sky(kv,disps,loads,kdiag)
!
! This subroutine forms the product of symmetric matrix stored as
! a skyline and a vector.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::kv(:),disps(0:)
 REAL(iwp),INTENT(OUT)::loads(0:)
 INTEGER,INTENT(IN)::kdiag(:)
 INTEGER::n,i,j,low,lup,k
 REAL(iwp)::x,zero=0.0_iwp
 n=UBOUND(disps,1)
 DO i=1,n
   x=zero 
   lup=kdiag(i)
   IF(i==1)low=lup
   IF(i/=1)low=kdiag(i-1)+1
   DO j=low,lup
     x=x+kv(j)*disps(i+j-lup) 
   END DO
   loads(i)=x
   IF(i==1)CYCLE   
   lup=lup-1
   DO j=low,lup
     k=i+j-lup-1
     loads(k)=loads(k)+kv(j)*disps(i)        
   END DO
 END DO
RETURN
END SUBROUTINE linmul_sky  
END SUBROUTINE stability
SUBROUTINE stiff4(km,coord,e,v)
!
! This subroutine generates the "analytical" stiffness matrix
! for a four node quadrilateral in plane strain based on nip=4.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::coord(:,:),e,v
 REAL(iwp),INTENT(OUT)::km(:,:)
 REAL(iwp)::e1,e2,g,x1,x2,x3,x4,y1,y2,y3,y4,a2,a2st3,alph,beta,c11,c21,   &
   c31,c41,c51,c61,c71,c81,c22,c32,c42,c52,c62,c72,c82,c33,c43,c53,c63,   &
   c73,c83,c44,c54,c64,c74,c84,c55,c65,c75,c85,c66,c76,c86,c77,c87,c88,   &
   f1,f2,s1,s2,s3,s4,t1,t2,t3,t4,pt5=0.5_iwp,one=1.0_iwp,two=2.0_iwp,     &
   d3=3.0_iwp,zero=0.0_iwp
 INTEGER::i,j
 e1=e*(one-v)/(one+v)/(1-two*v)
 e2=v*e1/(one-v)
 g=e/two/(one+v)      
!
 x1=coord(1,1)
 x2=coord(2,1)
 x3=coord(3,1)
 x4=coord(4,1)
 y1=coord(1,2)
 y2=coord(2,2)
 y3=coord(3,2)
 y4=coord(4,2)
!
 a2=(x4-x2)*(y3-y1)-(x3-x1)*(y4-y2)
 a2st3=d3*a2*a2
!
 CALL groupa(x1,x2,x3,x4,y1,y2,y3,y4,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c11=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(1,1)=c11
!
 CALL groupa(x2,x3,x4,x1,y2,y3,y4,y1,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c33=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(3,3)=c33
!
 CALL groupa(x3,x4,x1,x2,y3,y4,y1,y2,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c55=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(5,5)=c55
!
 CALL groupa(x4,x1,x2,x3,y4,y1,y2,y3,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c77=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(7,7)=c77
!
 CALL groupa(y4,y3,y2,y1,x4,x3,x2,x1,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c88=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(8,8)=c88
!
 CALL groupa(y1,y4,y3,y2,x1,x4,x3,x2,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c22=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(2,2)=c22
!
 CALL groupa(y2,y1,y4,y3,x2,x1,x4,x3,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c44=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(4,4)=c44
!
 CALL groupa(y3,y2,y1,y4,x3,x2,x1,x4,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c66=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(6,6)=c66
! 
 CALL groupb(x1,x2,x3,x4,y1,y2,y3,y4,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c21=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(2,1)=c21
!
 CALL groupb(x2,x3,x4,x1,y2,y3,y4,y1,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c43=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(4,3)=c43
!
 CALL groupb(x3,x4,x1,x2,y3,y4,y1,y2,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c65=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(6,5)=c65
!
 CALL groupb(x4,x1,x2,x3,y4,y1,y2,y3,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c87=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(8,7)=c87
!
 CALL groupc(x1,x2,x3,x4,y1,y2,y3,y4,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c31=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(3,1)=c31
!
 CALL groupc(x2,x3,x4,x1,y2,y3,y4,y1,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c53=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(5,3)=c53
!
 CALL groupc(x3,x4,x1,x2,y3,y4,y1,y2,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c75=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(7,5)=c75
!
 CALL groupc(x4,x1,x2,x3,y4,y1,y2,y3,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c71=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(7,1)=c71
! 
 CALL groupc(y4,y3,y2,y1,x4,x3,x2,x1,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c86=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(8,6)=c86
!
 CALL groupc(y1,y4,y3,y2,x1,x4,x3,x2,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c82=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(8,2)=c82
!
 CALL groupc(y2,y1,y4,y3,x2,x1,x4,x3,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c42=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(4,2)=c42
!
 CALL groupc(y3,y2,y1,y4,x3,x2,x1,x4,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c64=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(6,4)=c64
! 
 CALL groupd(x1,x2,x3,x4,y1,y2,y3,y4,s1,s2,s3,s4,t1,t2,t3,t4)
 CALL f1f2(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c41=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(4,1)=c41
!
 CALL groupd(x2,x3,x4,x1,y2,y3,y4,y1,s1,s2,s3,s4,t1,t2,t3,t4)
 CALL f1f2(x2,x3,x4,x1,y2,y3,y4,y1,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c63=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(6,3)=c63
!
 CALL groupd(x3,x4,x1,x2,y3,y4,y1,y2,s1,s2,s3,s4,t1,t2,t3,t4)
 CALL f1f2(x3,x4,x1,x2,y3,y4,y1,y2,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c85=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(8,5)=c85
!
 CALL groupd(x4,x1,x2,x3,y4,y1,y2,y3,s1,s2,s3,s4,t1,t2,t3,t4)
 CALL f1f2(x4,x1,x2,x3,y4,y1,y2,y3,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c72=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(7,2)=c72
!
 CALL groupd(y1,y2,y3,y4,x1,x2,x3,x4,s1,s2,s3,s4,t1,t2,t3,t4)
 CALL f1f2(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c32=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(3,2)=c32
!
 CALL groupd(y2,y3,y4,y1,x2,x3,x4,x1,s1,s2,s3,s4,t1,t2,t3,t4)
 CALL f1f2(x2,x3,x4,x1,y2,y3,y4,y1,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c54=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(5,4)=c54
!
 CALL groupd(y3,y4,y1,y2,x3,x4,x1,x2,s1,s2,s3,s4,t1,t2,t3,t4)
 CALL f1f2(x3,x4,x1,x2,y3,y4,y1,y2,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c76=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(7,6)=c76
!
 CALL groupd(y4,y1,y2,y3,x4,x1,x2,x3,s1,s2,s3,s4,t1,t2,t3,t4)
 CALL f1f2(x4,x1,x2,x3,y4,y1,y2,y3,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c81=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(8,1)=c81
!
 CALL groupe(x1,x2,x3,x4,y1,y2,y3,y4,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c51=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(5,1)=c51
!
 CALL groupe(x2,x3,x4,x1,y2,y3,y4,y1,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c73=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(7,3)=c73
!
 CALL groupe(y2,y1,y4,y3,x2,x1,x4,x3,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c84=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(8,4)=c84
!
 CALL groupe(y3,y2,y1,y4,x3,x2,x1,x4,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c62=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(6,2)=c62
!
 CALL groupf(x1,x2,x3,x4,y1,y2,y3,y4,s1,s2,s3,s4,t1,t2,t3,t4)
 CALL f1f2(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c61=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(6,1)=c61
!
 CALL groupf(x2,x3,x4,x1,y2,y3,y4,y1,s1,s2,s3,s4,t1,t2,t3,t4)
 CALL f1f2(x2,x3,x4,x1,y2,y3,y4,y1,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c83=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(8,3)=c83
!
 CALL groupf(y1,y2,y3,y4,x1,x2,x3,x4,s1,s2,s3,s4,t1,t2,t3,t4)
 CALL f1f2(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c52=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(5,2)=c52
!
 CALL groupf(y2,y3,y4,y1,x2,x3,x4,x1,s1,s2,s3,s4,t1,t2,t3,t4)
 CALL f1f2(x2,x3,x4,x1,y2,y3,y4,y1,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c74=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*pt5
 km(7,4)=c74
!
 DO i=1,8
   DO j=i+1,8 
     km(i,j)=km(j,i) 
   END DO 
 END DO
!
RETURN 
CONTAINS 
SUBROUTINE groupa(x1,x2,x3,x4,y1,y2,y3,y4,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::x1,x2,x3,x4,y1,y2,y3,y4
 REAL(iwp),INTENT(OUT)::s1,s2,s3,s4,t1,t2,t3,t4,f1,f2
 s1=two*(y4-y2)**2
 s2=two*(x4-x2)**2
 s3=-s1/two
 s4=-s2/two
 t1=(y2-y3)**2+(y3-y4)**2+(y4-y2)**2
 t2=(x2-x3)**2+(x3-x4)**2+(x4-x2)**2
 t3=(y4-y3)**2-(y3-y2)**2
 t4=(x4-x3)**2-(x3-x2)**2
 f1=(x1+x3)*(y4-y2)-(y1+y3)*(x4-x2)-two*(x2*y4-x4*y2)
 f2=(y2+y4)*(x3-x1)-(x2+x4)*(y3-y1)-two*(x3*y1-x1*y3)
RETURN
END SUBROUTINE groupa
SUBROUTINE groupb(x1,x2,x3,x4,y1,y2,y3,y4,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::x1,x2,x3,x4,y1,y2,y3,y4
 REAL(iwp),INTENT(OUT)::s1,s2,s3,s4,t1,t2,t3,t4,f1,f2
 s1=two*(x2-x4)*(y4-y2)
 s2=s1
 s3=-s1/2
 s4=s3
 t1=x2*(y4-two*y2+y3)+x3*(y2-two*y3+y4)+x4*(y2-two*y4+y3)
 t2=t1
 t3=x2*(y2-y3)+x3*(y4-y2)+x4*(y3-y4)
 t4=t3
 f1=(x1+x3)*(y4-y2)-(y1+y3)*(x4-x2)-two*(x2*y4-x4*y2)
 f2=(y2+y4)*(x3-x1)-(x2+x4)*(y3-y1)-two*(x3*y1-x1*y3)
RETURN
END SUBROUTINE groupb
SUBROUTINE groupc(x1,x2,x3,x4,y1,y2,y3,y4,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::x1,x2,x3,x4,y1,y2,y3,y4
 REAL(iwp),INTENT(OUT)::s1,s2,s3,s4,t1,t2,t3,t4,f1,f2
 s1=(y4-y2)*(two*y1-y3-y4)
 s2=(x4-x2)*(two*x1-x3-x4)
 s3=(y4-y2)*(y4-y1)
 s4=(x4-x2)*(x4-x1)
 t1=(y3-y1)*(two*y2-y3-y4)
 t2=(x3-x1)*(two*x2-x3-x4)
 t3=(y3-y1)*(y3-y2)
 t4=(x3-x1)*(x3-x2)
 f1=(x1+x3)*(y4-y2)-(y1+y3)*(x4-x2)-two*(x2*y4-x4*y2)
 f2=(y2+y4)*(x3-x1)-(x2+x4)*(y3-y1)-two*(x3*y1-x1*y3)
RETURN
END SUBROUTINE groupc
SUBROUTINE groupd(x1,x2,x3,x4,y1,y2,y3,y4,s1,s2,s3,s4,t1,t2,t3,t4)
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::x1,x2,x3,x4,y1,y2,y3,y4
 REAL(iwp),INTENT(OUT)::s1,s2,s3,s4,t1,t2,t3,t4
 s1=(x3-x1)*(y4-y2)+(x4-x1)*(y4-y2)
 s2=(y3-y1)*(x4-x2)+(y4-y1)*(x4-x2)
 s3=(x4-x1)*(y2-y4)
 s4=(y4-y1)*(x2-x4)
 t1=(x3-x1)*(y4-y2)+(x3-x1)*(y3-y2)
 t2=(y3-y1)*(x4-x2)+(y3-y1)*(x3-x2)
 t3=(x3-x1)*(y2-y3)
 t4=(y3-y1)*(x2-x3)
RETURN
END SUBROUTINE groupd
SUBROUTINE groupe(x1,x2,x3,x4,y1,y2,y3,y4,s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::x1,x2,x3,x4,y1,y2,y3,y4
 REAL(iwp),INTENT(OUT)::s1,s2,s3,s4,t1,t2,t3,t4,f1,f2
 s1=-(y4-y2)**2
 s2=-(x4-x2)**2
 s3=zero
 s4=zero
 t1=(y3+y1)*(y4+y2)-two*(y4-y2)**2-two*(y1*y3+y2*y4)
 t2=(x3+x1)*(x4+x2)-two*(x4-x2)**2-two*(x1*x3+x2*x4)
 t3=(y4-y2)*(y1-y2+y3-y4)
 t4=(x4-x2)*(x1-x2+x3-x4)
 f1=(x1+x3)*(y4-y2)-(y1+y3)*(x4-x2)-two*(x2*y4-x4*y2)
 f2=(y2+y4)*(x3-x1)-(x2+x4)*(y3-y1)-two*(x3*y1-x1*y3)
RETURN
END SUBROUTINE groupe
SUBROUTINE groupf(x1,x2,x3,x4,y1,y2,y3,y4,s1,s2,s3,s4,t1,t2,t3,t4)
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::x1,x2,x3,x4,y1,y2,y3,y4
 REAL(iwp),INTENT(OUT)::s1,s2,s3,s4,t1,t2,t3,t4
 s1=(x4-x2)*(y4-y2)
 s2=s1
 s3=zero
 s4=zero
 t1=(x4-x2)*(y4-y2)+(x2-x1)*(y2-y3)+(x4-x1)*(y4-y3)
 t2=(y4-y2)*(x4-x2)+(y2-y1)*(x2-x3)+(y4-y1)*(x4-x3)
 t3=(x2-x1)*(y3-y2)+(x4-x1)*(y4-y3)
 t4=(y2-y1)*(x3-x2)+(y4-y1)*(x4-x3)
RETURN
END SUBROUTINE groupf
SUBROUTINE f1f2(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2)
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::x1,x2,x3,x4,y1,y2,y3,y4
 REAL(iwp),INTENT(OUT)::f1,f2
 f1=(x1+x3)*(y4-y2)-(y1+y3)*(x4-x2)-two*(x2*y4-x4*y2)
 f2=(y2+y4)*(x3-x1)-(x2+x4)*(y3-y1)-two*(x3*y1-x1*y3)
RETURN
END SUBROUTINE f1f2
END SUBROUTINE stiff4

SUBROUTINE vecmsh(loads,nf,ratmax,cutoff,g_coord,g_num,argv,nlen,ips)
!
! This subroutine produces a PostScript output file "*.vec" displaying
! the nodal displacement vectors.
!
IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::g_coord(:,:),loads(0:),ratmax,cutoff
 INTEGER,INTENT(IN)::g_num(:,:),ips,nf(:,:),nlen
 REAL(iwp)::width,height,scale=72,sxy,xo,yo,x1,y1,x2,y2,dismag,           &
   zero=0.0_iwp,pt5=0.5_iwp,opt5=1.5_iwp,fpt5=5.5_iwp,d8=8.0_iwp,         &
   ept5=8.5_iwp,d11=11.0_iwp,xmin,xmax,ymin,ymax,dmax,vlen,vmax
 INTEGER::i,j,k,l,nn,nels,nod,ns,i1,i2,j1,j2
 INTEGER,ALLOCATABLE::corner(:,:)
 CHARACTER(*),INTENT(IN)::argv
 LOGICAL::draw
!                       formats
 OPEN(ips,FILE=argv(1:nlen)//'.vec')
!                       open output file and compute scale factors
 nn=UBOUND(nf,2)
!
 xmin=g_coord(1,1)
 xmax=g_coord(1,1)
 ymin=g_coord(2,1)
 ymax=g_coord(2,1)
 DO i=2,nn
   IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)      
   IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)      
   IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)      
   IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)      
 END DO
 width =xmax-xmin
 height=ymax-ymin
 dmax=ratmax*width
 IF(height>width)dmax=ratmax*height
!
 vmax=zero
 DO i=1,nn
   DO j=1,2
     IF(ABS(loads(nf(j,i)))>vmax)vmax=ABS(loads(nf(j,i)))
   END DO
 END DO
 dismag=dmax/vmax
!
 xmin=g_coord(1,1)
 xmax=g_coord(1,1)
 ymin=g_coord(2,1)
 ymax=g_coord(2,1)
!
 DO i=1,nn
   IF(g_coord(1,i)+dismag*loads(nf(1,i)) < xmin)                          &
     xmin=g_coord(1,i)+dismag*loads(nf(1,i))      
   IF(g_coord(1,i)+dismag*loads(nf(1,i)) > xmax)                          &
     xmax=g_coord(1,i)+dismag*loads(nf(1,i))      
   IF(g_coord(2,i)+dismag*loads(nf(2,i)) < ymin)                          &
     ymin=g_coord(2,i)+dismag*loads(nf(2,i))      
   IF(g_coord(2,i)+dismag*loads(nf(2,i)) > ymax)                          &
     ymax=g_coord(2,i)+dismag*loads(nf(2,i))      
!
   IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)
   IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)
   IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)
   IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)
 END DO
!
 width=xmax-xmin
 height=ymax-ymin
!
!                       allow 1.5" margin minimum on each side of figure
!
!                       Portrait mode 
!
 IF(height.GE.d11/ept5*width)THEN
!
!                       height governs the scale
!
   sxy=scale*d8/height
   xo=scale*pt5*(ept5-d8*width/height)
   yo=scale*opt5
 ELSE
!
!                       width governs the scale
!
   sxy=scale*fpt5/width
   xo=scale*opt5
   yo=scale*pt5*(d11-fpt5*height/width)
 END IF
!
 WRITE(ips,'(a)')'%!PS-Adobe-1.0'
 WRITE(ips,'(a)')'%%DocumentFonts: none'
 WRITE(ips,'(a)')'%%Pages: 1'
 WRITE(ips,'(a)')'%%EndComments'
 WRITE(ips,'(a)')'/m {moveto} def'
 WRITE(ips,'(a)')'/l {lineto} def'
 WRITE(ips,'(a)')'/s {stroke} def'
 WRITE(ips,'(a)')'/c {closepath} def'
 WRITE(ips,'(a)')'/edef {exch def} bind def'
 WRITE(ips,'(a)')                                                         &
   '/arrow {/@to_y edef /@to_x edef /@from_y edef /@from_x edef'
 WRITE(ips,'(a)')'/@dx @to_x @from_x sub def /@dy @to_y @from_y sub def'
 WRITE(ips,'(a)')'/@length @dx @dx mul @dy @dy mul add sqrt def'
 WRITE(ips,'(a)')'/@angle @dy @dx atan def'
 WRITE(ips,'(a)')'gsave @from_x @from_y translate @angle rotate'
 WRITE(ips,'(a)')                                                         &
   '0 0 moveto @length 0 lineto currentpoint stroke newpath moveto'
 WRITE(ips,'(a)')'-4 -2 rlineto @length 0 moveto'
 WRITE(ips,'(a)')'-4  2 rlineto stroke grestore'
 WRITE(ips,'(a)')'} def'
 WRITE(ips,'(a)')'/*sf {'
 WRITE(ips,'(a)')'exch findfont exch'
 WRITE(ips,'(a)')                                                         &
   'dup type /arraytype eq {makefont}{scalefont} ifelse setfont'
 WRITE(ips,'(a)')'} bind def'
 WRITE(ips,'(a)')'/languagelevel where'
 WRITE(ips,'(a)')'{pop languagelevel} {1} ifelse'
 WRITE(ips,'(a)')'2 lt { % ifelse'
 WRITE(ips,'(a)')'/sf /*sf load def'
 WRITE(ips,'(a)')'} { % else'
 WRITE(ips,'(a)')'/sf /selectfont load def'
 WRITE(ips,'(a)')'} ifelse'
 WRITE(ips,'(a)')'%%EndProlog'
 WRITE(ips,'(a)')'%%Page: 0 1'
 WRITE(ips,'(a)')'gsave'
!
 WRITE(ips,'(2f9.2,a)')xo,yo, ' translate'
 WRITE(ips,'(f9.2,a)')0.5,' setlinewidth'
!
!                       draw the displacement vectors
!
 vmax=zero
 DO i=1,nn
   vlen=loads(nf(1,i))**2+loads(nf(2,i))**2
   IF(vlen>vmax)vmax=vlen
 END DO
 vmax=SQRT(vmax)*cutoff      
 DO i=1,nn
   vlen=SQRT(loads(nf(1,i))**2+loads(nf(2,i))**2)
   x1=sxy*(g_coord(1,i)-xmin)
   y1=sxy*(g_coord(2,i)-ymin)
   x2=sxy*(g_coord(1,i)+dismag*loads(nf(1,i))-xmin)
   y2=sxy*(g_coord(2,i)+dismag*loads(nf(2,i))-ymin)
   IF(vlen>vmax)THEN
     WRITE(ips,'(2f9.2,a,2f9.2,a)') x1, y1,' ', x2, y2, ' arrow'
     WRITE(ips,'(a)') 's'
   END IF
 END DO
!
!                       draw the mesh border                          
!
 nels=UBOUND(g_num,2)
 nod=UBOUND(g_num,1)
 IF(nod==3.OR.nod==6.OR.nod==10.OR.nod==15)ns=3
 IF(nod==4.OR.nod==5.OR.nod==8.OR.nod==9)ns=4
 ALLOCATE(corner(ns,2))
 IF(nod== 3)corner=RESHAPE((/1,2,3,2,3,1/),(/3,2/))
 IF(nod== 6)corner=RESHAPE((/1,3,5,3,5,1/),(/3,2/))
 IF(nod==10)corner=RESHAPE((/1,4,7,4,7,1/),(/3,2/))
 IF(nod==15)corner=RESHAPE((/1,5,9,5,9,1/),(/3,2/))
 IF(nod== 4)corner=RESHAPE((/1,2,3,4,2,3,4,1/),(/4,2/))
 IF(nod== 5)corner=RESHAPE((/1,2,3,4,2,3,4,1/),(/4,2/))
 IF(nod== 8)corner=RESHAPE((/1,3,5,7,3,5,7,1/),(/4,2/))
 IF(nod== 9)corner=RESHAPE((/1,3,5,7,3,5,7,1/),(/4,2/))
 DO i=1,nels
   DO j=1,ns
     draw=.TRUE.
     i1=g_num(corner(j,1),i)
     i2=g_num(corner(j,2),i)
     DO k=1,nels
       DO l=1,ns
         j1=g_num(corner(l,1),k)
         j2=g_num(corner(l,2),k)
         IF((i1==j2).AND.(i2==j1))THEN
           draw=.FALSE.
           EXIT
         END IF
       END DO
       IF(.NOT.draw)EXIT
     END DO
     IF(draw)THEN
       x1=sxy*(g_coord(1,i1)-xmin)
       y1=sxy*(g_coord(2,i1)-ymin)
       WRITE(ips,'(2f9.2,a)')x1, y1,' m'
       x1=sxy*(g_coord(1,i2)-xmin)
       y1=sxy*(g_coord(2,i2)-ymin)
       WRITE(ips,'(2f9.2,a)')x1, y1,' l'
       WRITE(ips,'(a)')' s'
     END IF
   END DO
 END DO
!                       close output file?
 WRITE(ips,'(a)')'grestore'
 WRITE(ips,'(a)')'showpage'
 CLOSE(ips)
RETURN
END SUBROUTINE vecmsh
SUBROUTINE vmdpl(dee,stress,pl)
!
! This subroutine forms the plastic stress/strain matrix
! for a von-Mises material.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::stress(:),dee(:,:)  
 REAL(iwp),INTENT(OUT)::pl(:,:)
 REAL(iwp),ALLOCATABLE::dfds(:),ddfds(:)
 REAL(iwp)::t2,t6,t10,t14,t16,t17,t18,t19,t21,t22,t23,t25,t26,t30,sq2,sx, &
   sy,sz,txy,tyz,tzx,one=1.0_iwp,two=2.0_iwp,d4=4.0_iwp,d6=6.0_iwp
 REAL(iwp)::denom
 INTEGER::i,j,ih
 ih=SIZE(stress)
 ALLOCATE(dfds(ih),ddfds(ih))
 sq2=SQRT(two)
 SELECT CASE(ih)
 CASE(4)
   sx=stress(1)
   sy=stress(2)
   txy=stress(3) 
   sz=stress(4)   
   t2=sx**2
   t6=sy**2
   t10=sz**2
   t14=txy**2
   t17=SQRT(two*t2-two*sx*sy+two*t6-two*sy*sz+two*t10-two*sz*sx+d6*t14)
   t19=one/sq2/t17
   t21=two*sy
   t22=two*sz
   t26=two*sx
   dfds(1)=t19*(d4*sx-t21-t22)/two
   dfds(2)=t19*(-t26+d4*sy-t22)/two
   dfds(3)=t19*d6*txy
   dfds(4)=t19*(-t21+d4*sz-t26)/two
 CASE(6)
   sx=stress(1)
   sy=stress(2)
   sz=stress(3)   
   txy=stress(4) 
   tyz=stress(5) 
   tzx=stress(6) 
   t2=sx**2
   t6=sy**2
   t10=sz**2
   t14=txy**2
   t16=tyz**2
   t18=tzx**2
   t21=SQRT(two*t2-two*sx*sy+two*t6-two*sy*sz+two*t10-                    &
     two*sz*sx+d6*t14+d6*t16+d6*t18)
   t23=one/sq2/t21
   t25=two*sy
   t26=two*sz
   t30=two*sx
   dfds(1)=t23*(d4*sx-t25-t26)/two
   dfds(2)=t23*(-t30+d4*sy-t26)/two
   dfds(3)=t23*(-t25+d4*sz-t30)/two
   dfds(4)=t23*d6*txy
   dfds(5)=t23*d6*tyz
   dfds(6)=t23*d6*tzx
 END SELECT
 ddfds=MATMUL(dee,dfds) 
 denom=DOT_PRODUCT(ddfds,dfds)
 DO i=1,ih
   DO j=1,ih
     pl(i,j)=ddfds(i)*ddfds(j)/denom
   END DO
 END DO
 DEALLOCATE(dfds,ddfds)
RETURN
END SUBROUTINE vmdpl
SUBROUTINE vmflow(stress,dsbar,vmfl)
!
! This subroutine forms the von-Mises flow vector.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::stress(:),dsbar
 REAL(iwp),INTENT(OUT)::vmfl(:)
 REAL(iwp)::sigm,onept5=1.5_iwp,two=2.0_iwp,d3=3.0_iwp
 sigm=(stress(1)+stress(2)+stress(4))/d3
 vmfl(1)=stress(1)-sigm
 vmfl(2)=stress(2)-sigm
 vmfl(3)=stress(3)*two 
 vmfl(4)=stress(4)-sigm
 vmfl=vmfl*onept5/dsbar
RETURN
END SUBROUTINE vmflow

END MODULE main
