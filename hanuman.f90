MODULE Hanuman

CONTAINS

SUBROUTINE num_eq(eqn,strt)
    !
! This subroutine forms the eqn matrix.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN OUT)::eqn(:,:)
 integer,intent(in)::strt
 INTEGER::i,j,m
 m=strt
  DO j=1,UBOUND(eqn,2)
   DO i=1,UBOUND(eqn,1)
       m=m+1
       eqn(i,j)=m
       
     
   END DO
 END DO
RETURN
END SUBROUTINE num_eq

SUBROUTINE end_to_est(e_nd_num,eqn,e_st_v)
    !
    ! This subroutine finds the element steering vector from element node number and equation matrix.
    !
     IMPLICIT NONE
     INTEGER,INTENT(IN)::e_nd_num(:),eqn(:,:)  
     INTEGER,INTENT(INOUT)::e_st_v(:)
     INTEGER::i,k,nod,nodof 
     nod=UBOUND(e_nd_num,1) 
     nodof=UBOUND(eqn,1)
     DO i=1,nod
       k=i*nodof
       e_st_v(k-nodof+1:k)=eqn(:,e_nd_num(i))
     END DO
    RETURN
    END SUBROUTINE end_to_est
    
    SUBROUTINE truss(kele,ea,e_coord)
        !
        ! This subroutine forms the stiffness matrix of a
        ! general rod element (1-, 2- or 3-d).
        !
         IMPLICIT NONE
         INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
         REAL(iwp),INTENT(IN)::ea,e_coord(:,:)
         REAL(iwp),INTENT(OUT)::kele(:,:)
         INTEGER::ndim
         REAL(iwp)::ell,cs,sn,x1,x2,y1,y2,z1,z2,a,b,c,d,e,f,xl,yl,zl,one=1.0_iwp
         ndim=UBOUND(e_coord,2)
           x1=e_coord(1,1)
           y1=e_coord(1,2)
           x2=e_coord(2,1)
           y2=e_coord(2,2)
           ell=SQRT((y2-y1)**2+(x2-x1)**2)
           cs=(x2-x1)/ell
           sn=(y2-y1)/ell
           a=cs*cs
           b=sn*sn
           c=cs*sn
           kele(1,1)=a
           kele(3,3)=a
           kele(1,3)=-a
           kele(3,1)=-a
           kele(2,2)=b
           kele(4,4)=b
           kele(2,4)=-b
           kele(4,2)=-b
           kele(1,2)=c
           kele(2,1)=c
           kele(3,4)=c
           kele(4,3)=c
           kele(1,4)=-c
           kele(4,1)=-c
           kele(2,3)=-c
           kele(3,2)=-c
         kele=kele*ea/ell
        RETURN
    END SUBROUTINE truss


    SUBROUTINE assemble(kv,kele,e_st_v)
        !
        ! maps elements stiffness to corresponding rows and columns in global stiffness matrix
        !
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
        REAL(iwp),INTENT(IN)::kele(:,:)
        INTEGER,INTENT(IN)::e_st_v(:)
        REAL(iwp),INTENT(OUT)::kv(:,:)
        INTEGER::i,j,siz
        siz=ubound(e_st_v,1)

        DO i=1,siz
            DO j=1,siz
                kv(e_st_v(i),e_st_v(j))=kele(i,j)


            END do
        END do

        RETURN

    END SUBROUTINE assemble

    
    SUBROUTINE submat(neq,count,res_eq,ksub,kstruct,full_load,load_sub)

        !
        ! forms submatrix of stiffness and loads by eleminating displaccement boundary conditions from global stiffness matrix and full_load matrix
        !

        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
        REAL(iwp),INTENT(IN)::kstruct(:,:),full_load(:)
        INTEGER,INTENT(IN)::res_eq(:)
        REAL(iwp),INTENT(OUT)::ksub(:,:),load_sub(:)
        REAL(iwp),ALLOCATABLE::temp(:)
        INTEGER::i,j,k,count,neq,p,count1,count2,t=1,r,q=1,count3
        allocate(temp((neq-count)*(neq-count)))
        do i=1,neq
            do j=1,neq
                count1=0
                count2=0
                do k=1,count
                    if(i==res_eq(k))then
                        count1=count1+1
                    endif
                end do
                do k=1,count
                    if(j==res_eq(k))then
                        count2=count2+1
                    endif
                end do
                    if(count1==0.and.count2==0)then
                        temp(t)=kstruct(i,j)
                        
                        
                        t=t+1
                        
                    endif
            end do 
            
            
        end do
        
        t=1
        r=neq-count
        
        do i=1,r
            do j=1,r
                ksub(i,j)=temp(t)
                t=t+1
            end do
        end do
        
        do i=1,neq
            count3=0
            do k=1,count
             if(i==res_eq(k))then
                count3=count3+1
             endif
            end do
            if(count3==0)then
            load_sub(q)=full_load(i)
            q=q+1
            endif
        end do
        
        RETURN
    end subroutine submat

    SUBROUTINE mesh_size(element,nodp,nelsp,nnp,nxe,nye,nze)
        !
        !  This subroutine returns the number of elements (nels) and the number
        !  of nodes (nn) in a 2-d geometry-created mesh.
        !
         IMPLICIT NONE
         CHARACTER(LEN=15),INTENT(IN)::element
         INTEGER,INTENT(IN)::nodp,nxe,nye
         INTEGER,INTENT(IN),OPTIONAL::nze
         INTEGER,INTENT(OUT)::nelsp,nnp

         nelsp=nxe*nye
       IF(nodp==4)nnp=(nxe+1)*(nye+1)
       IF(nodp==5)nnp=(nxe+1)*(nye+1)+nxe*nye
       IF(nodp==8)nnp=(2*nxe+1)*(nye+1)+(nxe+1)*nye
       IF(nodp==9)nnp=(2*nxe+1)*(2*nye+1)

       RETURN
    END SUBROUTINE mesh_size


    SUBROUTINE geom_rect(element,iel,x_coords,y_coords,e_coordp,e_nd_nump,dir)
        !
        ! This subroutine forms the coordinates and connectivity for a
        ! rectangular mesh of rt. angled triangular elements (3, 6, 10 or 15-node)
        ! or quadrilateral elements (4, 8 or 9-node) counting in the
        ! x- or y-dir. 
        !
         IMPLICIT NONE
         INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
         REAL(iwp),INTENT(IN)::x_coords(:),y_coords(:)
         REAL(iwp),INTENT(OUT)::e_coordp(:,:)
         CHARACTER(LEN=15),INTENT(IN)::element
         CHARACTER(LEN=1),INTENT(IN)::dir
         INTEGER,INTENT(IN)::iel
         INTEGER,INTENT(OUT)::e_nd_nump(:)
         INTEGER::ip,iq,jel,fac1,nod,nxe,nye
         REAL(iwp)::pt5=0.5_iwp,two=2.0_iwp,d3=3.0_iwp 
         nxe=UBOUND(x_coords,1)-1
         nod=UBOUND(e_nd_nump,1)
         nye=UBOUND(y_coords,1)-1
         iq=(iel-1)/nxe+1
         ip=iel-(iq-1)*nxe
         IF(dir=='x'.OR.dir=='r')THEN
            e_nd_nump(1)=iq*(nxe+1)+ip		        		
            e_nd_nump(2)=(iq-1)*(nxe+1)+ip				
            e_nd_nump(3)=e_nd_nump(2)+1					
            e_nd_nump(4)=e_nd_nump(1)+1
         ENDIF
         e_coordp(1:2,1)=x_coords(ip)
         e_coordp(3:4,1)=x_coords(ip+1)
         e_coordp(1,2)=y_coords(iq+1)
         e_coordp(2:3,2)=y_coords(iq)
         e_coordp(4,2)=e_coordp(1,2)
         RETURN
        END SUBROUTINE geom_rect

    SUBROUTINE sample(s,wt)
            !
            ! This subroutine returns the local coordinates and weighting coefficients
            ! of the integrating points.
            !
             IMPLICIT NONE
             INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
             REAL(iwp),INTENT(OUT)::s(:,:)
             REAL(iwp),INTENT(OUT)::wt(:)
             !CHARACTER(*),INTENT(IN)::element
             INTEGER::nip
             REAL(iwp)::root3,r15,w(3),v(9),b,c
             root3=1.0_iwp/SQRT(3.0_iwp)
             r15=0.2_iwp*SQRT(15.0_iwp)
             nip=UBOUND(s,1)
             w=(/5.0_iwp/9.0_iwp,8.0_iwp/9.0_iwp,5.0_iwp/9.0_iwp/)
             v=(/5.0_iwp/9.0_iwp*w,8.0_iwp/9.0_iwp*w,5.0_iwp/9.0_iwp*w/)
            s(1,1)=-root3
            s(1,2)= root3
            s(2,1)= root3
            s(2,2)= root3
            s(3,1)=-root3
            s(3,2)=-root3
            s(4,1)= root3
            s(4,2)=-root3
            wt=1.0_iwp
              
            RETURN
    END SUBROUTINE sample
    SUBROUTINE cmat(dee,e,v)
            !
            ! This subroutine returns the elastic constitutive matrix for ih=3 (plane stress),
            ! 
            !
            
             IMPLICIT NONE
             INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
             REAL(iwp),INTENT(IN)::e,v
             REAL(iwp),INTENT(OUT)::dee(:,:)
             REAL(iwp)::v1,v2,c,vv,zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp,two=2.0_iwp
             INTEGER::i,ih
             dee=zero  
             ih=UBOUND(dee,1)
             print*,"reached sub"
             c=e/(one-(v**2))
             dee(1,1)=1
             dee(2,2)=1
             dee(1,2)=v*c
             dee(2,1)=v*c
             dee(3,3)=pt5*c*(one-v)
            RETURN
    END SUBROUTINE cmat

        
    SUBROUTINE shape_fun(fun,points,i)
            !
            !   This subroutine computes the values of the shape functions.
            !   to local coordinates
            !
             IMPLICIT NONE
             INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
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
             c1=points(i,1)
             c2=points(i,2)
             c3=one-c1-c2 
             xi=points(i,1)
             eta=points(i,2)
             etam=pt25*(one-eta)
             etap=pt25*(one+eta)
             xim=pt25*(one-xi)
             xip=pt25*(one+xi)
             fun=(/d4*xim*etam,d4*xim*etap,d4*xip*etap,d4*xip*etam/)
             
             RETURN
    END SUBROUTINE shape_fun

    SUBROUTINE shape_der(der,points,i)
        !
        !   This subroutine produces derivatives of shape functions withe respect
        !   to local coordinates.
        !
         IMPLICIT NONE
         INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
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
         t1=-one-xi 
        t2=-one/d3-xi 
        t3=one/d3-xi 
        t4=one-xi
        der(1,1)=-(t3*t4+t2*t4+t2*t3)*d9/d16     
        der(1,2)=(t3*t4+t1*t4+t1*t3)*d27/d16 
        der(1,3)=-(t2*t4+t1*t4+t1*t2)*d27/d16 
        der(1,4)=(t2*t3+t1*t3+t1*t2)*d9/d16 

        RETURN
    END SUBROUTINE shape_der

    SUBROUTINE beemat(bee,deriv)
        !
        ! This subroutine forms the bee matrix in 2-d (ih=3 or 4) or 3-d (ih=6).
        !
         IMPLICIT NONE
         INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
         REAL(iwp),INTENT(IN)::deriv(:,:)
         REAL(iwp),INTENT(OUT)::bee(:,:)
         INTEGER::k,l,m,n,ih,nod
         REAL::x,y,z
         bee=0.0_iwp
         ih=UBOUND(bee,1)
         nod=UBOUND(deriv,2)
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
        RETURN
        END SUBROUTINE 
        FUNCTION determinant(jac)RESULT(det)
            !
            ! This function returns the determinant of a 1x1, 2x2 or 3x3
            ! Jacobian matrix.
            !
             IMPLICIT NONE    
             INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
             REAL(iwp),INTENT(IN)::jac(:,:)
             REAL(iwp)::det
             INTEGER::it 
             it=UBOUND(jac,1) 
             det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
             RETURN
        END FUNCTION determinant

        SUBROUTINE invert(matrix)
            !
            ! This subroutine inverts a small square matrix onto itself.
            !
             IMPLICIT NONE
             INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
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
             endif
                RETURN
        END SUBROUTINE invert


        SUBROUTINE beam(km,prop,etype,iel,coord) 
            !
            ! This subroutine forms the stiffness matrix of a
            ! general beam/column element (1-, 2- or 3-d).
            !
             IMPLICIT NONE
             INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
             REAL(iwp),INTENT(IN)::coord(:,:),prop(:,:)
             INTEGER,INTENT(IN)::etype(:),iel
             REAL(iwp),INTENT(OUT)::km(:,:)
             INTEGER::ndim,i,j,k
             REAL(iwp)::ell,x1,x2,y1,y2,z1,z2,c,s,e1,e2,e3,e4,pi,xl,yl,zl,cg,sg,den,  &
               ea,ei,eiy,eiz,gj,a1,a2,a3,a4,a5,a6,a7,a8,sum,gamrad,x,t(12,12),        &
               tt(12,12),cc(12,12),r0(3,3),zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp,      &
               two=2.0_iwp,d4=4.0_iwp,d6=6.0_iwp,d12=12.0_iwp,d180=180.0_iwp
             ndim=UBOUND(coord,2)
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
   RETURN
        END SUBROUTINE beam

END MODULE Hanuman