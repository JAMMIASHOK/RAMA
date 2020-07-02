MODULE Hanuman

CONTAINS

SUBROUTINE num_eq(eqn)
    !
! This subroutine forms the eqn matrix.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN OUT)::eqn(:,:)
 INTEGER::i,j,m
 m=0
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
    ! This subroutine finds the g vector from num and nf.
    !
     IMPLICIT NONE
     INTEGER,INTENT(IN)::e_nd_num(:),eqn(:,:)  
     INTEGER,INTENT(OUT)::e_st_v(:)
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
        INTEGER::i,j,siz=4
        

        DO i=1,siz
            DO j=1,siz
                kv(e_st_v(i),e_st_v(j))=kele(i,j)


            END do
        END do

        RETURN

    END SUBROUTINE assemble

    
    SUBROUTINE submat(neq,count,res_eq,ksub,kstruct)

        !
        ! forms submatrix by eleminating displaccement boundary conditions from global stiffness matrix
        !

        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
        REAL(iwp),INTENT(IN)::kstruct(:,:)
        INTEGER,INTENT(IN)::res_eq(:)
        REAL(iwp),INTENT(OUT)::ksub(:,:)
        REAL(iwp),ALLOCATABLE::temp(:)
        INTEGER::i,j,k,count,neq,p,count1,count2,t=1,r
        allocate(temp(neq-count))
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

        RETURN
    end subroutine submat




    

END MODULE Hanuman