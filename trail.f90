program trail
    implicit none
    integer,allocatable::a(:,:)
    integer::i,j

allocate(a(2,2))
do i=1,2
    do j=1,2
     a(i,j)=i+j
    end do
end do
print*,a
deallocate(a)

allocate(a(3,3))
do i=1,3
    do j=1,3
     a(i,j)=i+j
    end do
end do
print*,a


end program trail



!___________________________________________________________________________!
!                                                                           !
!                           PLATE COMPUTATIONS                              !
!___________________________________________________________________________!

allocate(eqnb(kele(ndofp,ndofp)),e_coord(nodp,ndim),g_coord(ndim,nnp),    &
e_disp(ndofp),actionp(ndofp),g_nd_num(nodp,nelsp),e_nd_num(nodp),e_st_v(ndofp),g_e_st_m(ndofp,nelsp),&
petype(nelsp),pprop(nprops,pnp_types),x_coords(nxe+1),y_coords(nye+1),points(nip,ndim),shp_fun(nodp),jac(ndim,ndim),&
loc_sp_der(ndim,nodp),sp_der(ndim,nodp),str_disp_m(nst,ndofp))
READ(10,*)pprop
 petype=1
 IF(pnp_types>1)read(10,*)petype
 READ(10,*)x_coords,y_coords

!-----------------------loop the elements to find global arrays sizes-----
 
 elements_1: DO iel=1,nelsp
   CALL geom_rect(element,iel,x_coords,y_coords,e_coord,e_nd_num,dir)
   CALL num_to_g(num,nf,g)
   CALL fkdiag(kdiag,g)
   g_num(:,iel)=num
   g_coord(:,num)=TRANSPOSE(e_coord)
   g_g(:,iel)=g
 END DO elements_1
 


























!----------------------boundary conditions----------------------
  bndry=eqn
  READ(10,*)nr,(k,bndry(:,k),i=1,nr)  !boundary conditions
    do i=1,ndim
      do j=1,nn
        IF(bndry(i,j)==0)count=count+1
    end do
  end do
  ALLOCATE(ksub(neq-count,neq-count),load_sub(neq-count))
  ALLOCATE(res_eq(count))
  ksub=0
  count=0
  do j=1,nn
    do i=1,nodof
      IF(bndry(i,j)==0) THEN
      count=count+1
      res_eq(count)=eqn(i,j)
      ENDIF
    end do
  end do
  print*,res_eq
 !--------------loads----------------------
  load=0
  read(10,*)loaded_nodes,(k,load(:,k),i=1,loaded_nodes)
  print*,load
  allocate(load_eq(loaded_nodes),load_val(loaded_nodes))
  
  do j=1,nn
    do i=1,nodof
      
      if(load(i,j)/=0)then
        load_eq(load_count)=eq_num
        load_val(load_count)=load(i,j)
        load_count=load_count+1
      endif
      eq_num=eq_num+1
    end do
  end do
  
  full_load=0
  do i=1,load_count-1
    full_load(load_eq(i))=load_val(i)
  end do
  print*,full_load
 
 CALL submat(neq,count,res_eq,ksub,kstruct,full_load,load_sub)

 
 write(unit=11,fmt=255)ksub
 255 format(9f12.3)
! write(11,*)" "
 write(unit=11,fmt=265)load_sub
 265 format(9f12.3)
!----------------------------SOLUTION PROCEDURE----------------------------------




























 SUBROUTINE mesh_size(element,nodp,nelsp,nnp,nxe,nye,nze)
    !
    !  This subroutine returns the number of elements (nels) and the number
    !  of nodes (nn) in a 2-d geometry-created mesh.
    !
     IMPLICIT NONE
     CHARACTER(LEN=15),INTENT(IN)::element
     INTEGER,INTENT(IN)::nod,nxe,nye
     INTEGER,INTENT(IN),OPTIONAL::nze
     INTEGER,INTENT(OUT)::nels,nn
     IF(element=="triangle")THEN
       nels=nxe*nye*2
       IF(nod==3)nn=(nxe+1)*(nye+1)
       IF(nod==6)nn=(2*nxe+1)*(2*nye+1)
       IF(nod==10)nn=(3*nxe+1)*(3*nye+1)
       IF(nod==15)nn=(4*nxe+1)*(4*nye+1)
     ELSE IF(element=="quadrilateral")THEN
       nels=nxe*nye
       IF(nod==4)nn=(nxe+1)*(nye+1)
       IF(nod==5)nn=(nxe+1)*(nye+1)+nxe*nye
       IF(nod==8)nn=(2*nxe+1)*(nye+1)+(nxe+1)*nye
       IF(nod==9)nn=(2*nxe+1)*(2*nye+1)
     ELSE IF(element=="hexahedron")THEN
       nels=nxe*nye*nze
       IF(nod==8)nn=(nxe+1)*(nye+1)*(nze+1)
       IF(nod==14)nn=4*nxe*nye*nze+2*(nxe*nye+nye*nze+nze*nxe)+nxe+nye+nze+1
       IF(nod==20)nn=((2*nxe+1)*(nze+1)+(nxe+1)*nze)*(nye+1)+                 &
         (nxe+1)*(nze+1)*nye
     END IF
    RETURN
    END SUBROUTINE mesh_size

    SUBROUTINE geom_rect(element,iel,x_coords,y_coords,coord,num,dir)
        !
        ! This subroutine forms the coordinates and connectivity for a
        ! rectangular mesh of rt. angled triangular elements (3, 6, 10 or 15-node)
        ! or quadrilateral elements (4, 8 or 9-node) counting in the
        ! x- or y-dir. 
        !
         IMPLICIT NONE
         INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
         REAL(iwp),INTENT(IN)::x_coords(:),y_coords(:)
         REAL(iwp),INTENT(OUT)::coord(:,:)
         CHARACTER(LEN=15),INTENT(IN)::element
         CHARACTER(LEN=1),INTENT(IN)::dir
         INTEGER,INTENT(IN)::iel
         INTEGER,INTENT(OUT)::num(:)
         INTEGER::ip,iq,jel,fac1,nod,nxe,nye
         REAL(iwp)::pt5=0.5_iwp,two=2.0_iwp,d3=3.0_iwp 
         nxe=UBOUND(x_coords,1)-1
         nod=UBOUND(num,1)
         IF(element=='triangle')THEN
           nye=(UBOUND(y_coords,1)-1)*2
           IF(dir=='x'.OR.dir=='r')THEN
             jel=2*nxe*((iel-1)/(2*nxe))
             ip=(iel-jel+1)/2
             iq=2*((iel-1)/(2*nxe)+1)-1+((iel/2)*2)/iel
           ELSE  
             jel=(iel-1)/nye
             ip=jel+1
             iq=iel-nye*jel
           END IF
           SELECT CASE(nod)
           CASE(3)
             IF(MOD(iq,2)/=0)THEN
               IF(dir=='x'.OR.dir=='r')THEN
                 num(1)=(nxe+1)*(iq-1)/2+ip
                 num(2)=num(1)+1              
                 num(3)=(nxe+1)*(iq+1)/2+ip
               ELSE
                 num(1)=(ip-1)*(nye+2)/2+(iq+1)/2
                 num(2)=num(1)+(nye+2)/2
                 num(3)=num(1)+1
               END IF
        !
               coord(1,1)=x_coords(ip)
               coord(1,2)=y_coords((iq+1)/2)
               coord(2,1)=x_coords(ip+1)   
               coord(2,2)=y_coords((iq+1)/2)
               coord(3,1)=x_coords(ip)   
               coord(3,2)=y_coords((iq+3)/2)
             ELSE
               IF(dir=='x'.OR.dir=='r')THEN
                 num(1)=(nxe+1)*iq/2+ip+1     
                 num(2)=num(1)-1               
                 num(3)=(nxe+1)*(iq-2)/2+ip+1
               ELSE
                 num(1)=ip*(nye+2)/2+(iq+2)/2
                 num(2)=(ip-1)*(nye+2)/2+(iq+1)/2+1
                 num(3)=num(1)-1
               END IF
        !
               coord(1,1)=x_coords(ip+1)
               coord(1,2)=y_coords((iq+2)/2)
               coord(2,1)=x_coords(ip)   
               coord(2,2)=y_coords((iq+2)/2)
               coord(3,1)=x_coords(ip+1) 
               coord(3,2)=y_coords(iq/2)
             END IF
           CASE(6)
             IF(MOD(iq,2)/=0)THEN
               IF(dir=='x'.OR.dir=='r')THEN
                 num(1)=(iq-1)*(2*nxe+1)+2*ip-1
                 num(2)=num(1)+1 
                 num(3)=num(1)+2 
                 num(4)=(iq-1)*(2*nxe+1)+2*nxe+2*ip+1
                 num(5)=(iq+1)*(2*nxe+1)+2*ip-1
                 num(6)=num(4)-1 
               ELSE
                 num(1)=2*(nye+1)*(ip-1)+iq
                 num(2)=2*(nye+1)*(ip-1)+nye+1+iq
                 num(3)=2*(nye+1)*ip+iq
                 num(4)=num(2)+1
                 num(5)=num(1)+2 
                 num(6)=num(1)+1
               END IF
        !
               coord(1,1)=x_coords(ip)
               coord(1,2)=y_coords((iq+1)/2)
               coord(3,1)=x_coords(ip+1)   
               coord(3,2)=y_coords((iq+1)/2)
               coord(5,1)=x_coords(ip)   
               coord(5,2)=y_coords((iq+3)/2)
             ELSE
               IF(dir=='x'.OR.dir=='r')THEN
                 num(1)=iq*(2*nxe+1)+2*ip+1
                 num(2)=num(1)-1 
                 num(3)=num(1)-2 
                 num(4)=(iq-2)*(2*nxe+1)+2*nxe+2*ip+1
                 num(5)=(iq-2)*(2*nxe+1)+2*ip+1
                 num(6)=num(4)+1 
               ELSE 
                 num(1)=2*(nye+1)*ip+iq+1 
                 num(2)=2*(nye+1)*(ip-1)+nye+iq+2
                 num(3)=2*(nye+1)*(ip-1)+iq+1
                 num(4)=num(2)-1 
                 num(5)=num(1)-2
                 num(6)=num(1)-1
               END IF
        !
               coord(1,1)=x_coords(ip+1)
               coord(1,2)=y_coords((iq+2)/2)
               coord(3,1)=x_coords(ip)   
               coord(3,2)=y_coords((iq+2)/2)
               coord(5,1)=x_coords(ip+1) 
               coord(5,2)=y_coords(iq/2)
             END IF
             coord(2,:)=pt5*(coord(1,:)+coord(3,:))
             coord(4,:)=pt5*(coord(3,:)+coord(5,:))
             coord(6,:)=pt5*(coord(5,:)+coord(1,:))
           CASE(10)
             IF(MOD(iq,2)/=0)THEN
               IF(dir=='x'.OR.dir=='r')THEN
                 num(1)=(iq-1)/2*(3*nxe+1)*3+3*ip-2
                 num(2)=num(1)+1
                 num(3)=num(1)+2
                 num(4)=num(1)+3
                 num(5)=(iq-1)/2*(3*nxe+1)*3+3*nxe+1+3*ip
                 num(6)=(iq-1)/2*(3*nxe+1)*3+6*nxe+2+3*ip-1
                 num(7)=(iq-1)/2*(3*nxe+1)*3+9*nxe+3+3*ip-2
                 num(8)=num(6)-1
                 num(9)=num(5)-2
                 num(10)=num(9)+1
               ELSE
                 num(1)=(9*(nye-2)/2+12)*(ip-1)+3*(iq-1)/2+1
                 num(2)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)/2+4+3*(iq-1)/2+1
                 num(3)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)+8+3*(iq-1)/2+1
                 num(4)=(9*(nye-2)/2+12)*(ip-1)+9*(nye-2)/2+12+3*(iq-1)/2+1
                 num(5)=num(3)+1 
                 num(6)=num(2)+2
                 num(7)=num(1)+3
                 num(8)=num(1)+2
                 num(9)=num(1)+1
                 num(10)=num(2)+1
               END IF
        !
               coord(1,1)=x_coords(ip)
               coord(2,1)=x_coords(ip)+(x_coords(ip+1)-x_coords(ip))/d3
               coord(3,1)=x_coords(ip)+two*(x_coords(ip+1)-x_coords(ip))/d3
               coord(4,1)=x_coords(ip+1)
               coord(4,2)=y_coords((iq+1)/2)
               coord(5,2)=y_coords((iq+1)/2)+                                     &
                 (y_coords((iq+3)/2)-y_coords((iq+1)/2))/d3
               coord(6,2)=y_coords((iq+1)/2)+                                     &
                 two*(y_coords((iq+3)/2)-y_coords((iq+1)/2))/d3
               coord(7,2)=y_coords((iq+3)/2)
             ELSE
               IF(dir=='x'.OR.dir=='r')THEN
                 num(1)=(iq-2)/2*(3*nxe+1)*3+9*nxe+3+3*ip+1
                 num(2)=num(1)-1
                 num(3)=num(1)-2
                 num(4)=num(1)-3
                 num(5)=(iq-2)/2*(3*nxe+1)*3+6*nxe+2+3*ip-1
                 num(6)=(iq-2)/2*(3*nxe+1)*3+3*nxe+1+3*ip
                 num(7)=(iq-2)/2*(3*nxe+1)*3+3*ip+1
                 num(8)=num(6)+1
                 num(9)=num(5)+2
                 num(10)=num(9)-1
               ELSE
                 num(1)=(9*(nye-2)/2+12)*(ip-1)+9*(nye-2)/2+12+3*iq/2+1
                 num(2)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)+8+3*iq/2+1
                 num(3)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)/2+4+3*iq/2+1
                 num(4)=(9*(nye-2)/2+12)*(ip-1)+3*iq/2+1
                 num(5)=num(3)-1
                 num(6)=num(2)-2
                 num(7)=num(1)-3
                 num(8)=num(1)-2
                 num(9)=num(1)-1
                 num(10)=num(2)-1
               END IF
        !
               coord(1,1)=x_coords(ip+1)
               coord(2,1)=x_coords(ip+1)-(x_coords(ip+1)-x_coords(ip))/d3
               coord(3,1)=x_coords(ip+1)-two*(x_coords(ip+1)-x_coords(ip))/d3
               coord(4,1)=x_coords(ip)
               coord(4,2)=y_coords((iq+2)/2)
               coord(5,2)=y_coords((iq+2)/2)-(y_coords((iq+2)/2)-y_coords(iq/2))/d3
               coord(6,2)=y_coords((iq+2)/2)-                                     &
                 two*(y_coords((iq+2)/2)-y_coords(iq/2))/d3
               coord(7,2) =y_coords(iq/2)
             END IF
             coord(5,1)=coord(3,1)
             coord(6,1)=coord(2,1)
             coord(7,1)=coord(1,1)
             coord(8,1)=coord(1,1)
             coord(9,1)=coord(1,1)
             coord(10,1)=coord(2,1)
             coord(1,2)=coord(4,2)
             coord(2,2)=coord(4,2)
             coord(3,2)=coord(4,2)
             coord(8,2)=coord(6,2)
             coord(9,2)=coord(5,2)
             coord(10,2)=coord(5,2)
           CASE(15)
             IF(MOD(iq,2)/=0)THEN
               IF(dir=='x'.OR.dir=='r')THEN
               fac1=4*(4*nxe+1)*(iq-1)/2
                 num(1)=fac1+4*ip-3
                 num(2)=num(1)+1
                 num(3)=num(1)+2
                 num(4)=num(1)+3
                 num(5)=num(1)+4
                 num(6)=fac1+ 4*nxe+1+4*ip
                 num(7)=fac1+ 8*nxe+1+4*ip
                 num(8)=fac1+12*nxe+1+4*ip
                 num(9)=fac1+16*nxe+1+4*ip
                 num(10)=num(8)-1
                 num(11)=num(7)-2
                 num(12)=num(6)-3
                 num(13)=num(12)+1
                 num(14)=num(12)+2
                 num(15)=num(11)+1
               ELSE
                 fac1=4*(2*nye+1)*(ip-1)+2*iq-1 
                 num(1)=fac1
                 num(2)=fac1+2*nye+1
                 num(3)=fac1+4*nye+2 
                 num(4)=fac1+6*nye+3 
                 num(5)=fac1+8*nye+4
                 num(6)=fac1+6*nye+4 
                 num(7)=fac1+4*nye+4 
                 num(8)=fac1+2*nye+4
                 num(9)=fac1+4 
                 num(10)=fac1+3 
                 num(11)=fac1+2 
                 num(12)=fac1+1
                 num(13)=fac1+2*nye+2 
                 num(14)=fac1+4*nye+3
                 num(15)=fac1+2*nye+3  
               END IF
        !
               coord(1,1)=x_coords(ip)
               coord(1,2)=y_coords((iq+1)/2)
               coord(5,1)=x_coords(ip+1)   
               coord(5,2)=y_coords((iq+1)/2)
               coord(9,1)=x_coords(ip)   
               coord(9,2)=y_coords((iq+3)/2)
             ELSE
               IF(dir=='x'.OR.dir=='r')THEN
                 fac1=4*(4*nxe+1)*(iq-2)/2
                 num(1)=fac1+16*nxe+5+4*ip
                 num(2)=num(1)-1
                 num(3)=num(1)-2
                 num(4)=num(1)-3
                 num(5)=num(1)-4
                 num(6)=fac1+12*nxe+1+4*ip
                 num(7)=fac1+8*nxe+1+4*ip
                 num(8)=fac1+4*nxe+1+4*ip
                 num(9)=fac1+4*ip+1
                 num(10)=num(8)+1
                 num(11)=num(7)+2
                 num(12)=num(6)+3
                 num(13)=num(12)-1
                 num(14)=num(12)-2
                 num(15)=num(11)-1
               ELSE
                 fac1=4*(2*nye+1)*(ip-1)+2*iq+8*nye+5 
                 num(1)=fac1 
                 num(2)=fac1-2*nye-1
                 num(3)=fac1-4*nye-2 
                 num(4)=fac1-6*nye-3 
                 num(5)=fac1-8*nye-4
                 num(6)=fac1-6*nye-4  
                 num(7)=fac1-4*nye-4 
                 num(8)=fac1-2*nye-4
                 num(9)=fac1-4
                 num(10)=fac1-3 
                 num(11)=fac1-2 
                 num(12)=fac1-1
                 num(13)=fac1-2*nye-2  
                 num(14)=fac1-4*nye-3
                 num(15)=fac1-2*nye-3 
               END IF
        !
               coord(1,1)=x_coords(ip+1)
               coord(1,2)=y_coords((iq+2)/2)
               coord(5,1)=x_coords(ip)   
               coord(5,2)=y_coords((iq+2)/2)
               coord(9,1)=x_coords(ip+1) 
               coord(9,2)=y_coords(iq/2)
             END IF
             coord(3,:)=pt5*(coord(1,:)+coord(5,:))
             coord(7,:)=pt5*(coord(5,:)+coord(9,:))
             coord(11,:)=pt5*(coord(9,:)+coord(1,:))
             coord(2,:)=pt5*(coord(1,:)+coord(3,:))
             coord(4,:)=pt5*(coord(3,:)+coord(5,:))
             coord(6,:)=pt5*(coord(5,:)+coord(7,:))
             coord(8,:)=pt5*(coord(7,:)+coord(9,:))
             coord(10,:)=pt5*(coord(9,:)+coord(11,:))
             coord(12,:)=pt5*(coord(11,:)+coord(1,:))
             coord(15,:)=pt5*(coord(7,:)+coord(11,:))
             coord(14,:)=pt5*(coord(3,:)+coord(7,:))
             coord(13,:)=pt5*(coord(2,:)+coord(15,:))
           CASE DEFAULT
             WRITE(11,'(a)')"Wrong number of nodes for triangular element"
             STOP
           END SELECT
         ELSE
           nye=UBOUND(y_coords,1)-1
           IF(dir=='x'.OR.dir=='r')THEN
             iq=(iel-1)/nxe+1
             ip=iel-(iq-1)*nxe
           ELSE
             ip=(iel-1)/nye+1
             iq=iel-(ip-1)*nye
           END IF
           SELECT CASE(nod)
           CASE(4)
             IF(dir=='x'.OR.dir=='r')THEN
               e_nd_num(1)=iq*(nxe+1)+ip		        		
               e_nd_num(2)=(iq-1)*(nxe+1)+ip				
               e_nd_num(3)=e_nd_num(2)+1					
               e_nd_num(4)=e_nd_num(1)+1					
             ELSE
               num(1)=(ip-1)*(nye+1)+iq+1
               num(2)=num(1)-1
               num(3)=ip*(nye+1)+iq
               num(4)=num(3)+1
             END IF
        !
             coord(1:2,1)=x_coords(ip)
             coord(3:4,1)=x_coords(ip+1)
             coord(1,2)=y_coords(iq+1)
             coord(2:3,2)=y_coords(iq)
             coord(4,2)=coord(1,2)
           CASE(5)
             IF(dir=='x'.OR.dir=='r')THEN
               num(1)=iq*(2*nxe+1)+ip
               num(2)=(iq-1)*(2*nxe+1)+ip
               num(3)=num(2)+1
               num(4)=num(1)+1
               num(5)=iq*(2*nxe+1)+ip-nxe
             ELSE
               num(1)=(ip-1)*(2*nye+1)+iq+1
               num(2)=num(1)-1
               num(3)=ip*(2*nye+1)+iq
               num(4)=num(3)+1
               num(5)=ip*(2*nye+1)+iq-nye
             END IF
        !
             coord(1:2,1)=x_coords(ip)
             coord(3:4,1)=x_coords(ip+1)
             coord(1,2)=y_coords(iq+1)
             coord(2:3,2)=y_coords(iq)
             coord(4,2)=coord(1,2)
             coord(5,:)=0.25_iwp*(coord(1,:)+coord(2,:)+coord(3,:)+coord(4,:))
           CASE(8)
             IF(dir=='x'.OR.dir=='r')THEN
               num(1)=iq*(3*nxe+2)+2*ip-1                 
               num(2)=iq*(3*nxe+2)+ip-nxe-1		  
               num(3)=(iq-1)*(3*nxe+2)+2*ip-1		   
               num(4)=num(3)+1
               num(5)=num(4)+1
               num(6)=num(2)+1
               num(7)=num(1)+2
               num(8)=num(1)+1
             ELSE
               num(1)=(ip-1)*(3*nye+2)+2*iq+1
               num(2)=num(1)-1
               num(3)=num(1)-2
               num(4)=(ip-1)*(3*nye+2)+2*nye+iq+1
               num(5)=ip*(3*nye+2)+2*iq-1
               num(6)=num(5)+1
               num(7)=num(5)+2
               num(8)=num(4)+1
             END IF
        !
             coord(1:3,1)=x_coords(ip)
             coord(5:7,1)=x_coords(ip+1)
             coord(4,1)=pt5*(coord(3,1)+coord(5,1))
             coord(8,1)=pt5*(coord(7,1)+coord(1,1))
             coord(1,2)=y_coords(iq+1)
             coord(7:8,2)=y_coords(iq+1)
             coord(3:5,2)=y_coords(iq)
             coord(2,2)=pt5*(coord(1,2)+coord(3,2))
             coord(6,2)=pt5*(coord(5,2)+coord(7,2))
           CASE(9)
             IF(dir=='x'.OR.dir=='r')THEN
               num(1)=iq*(4*nxe+2)+2*ip-1
               num(2)=iq*(4*nxe+2)+2*ip-nxe-4
               num(3)= (iq-1)*(4*nxe+2)+2*ip-1
               num(4)=num(3)+1
               num(5)=num(4)+1
               num(6)=num(2)+2
               num(7)=num(1)+2
               num(8)=num(1)+1
               num(9)=num(2)+1
             ELSE
               num(1)=(ip-1)*2*(2*nye+1)+2*iq+1
               num(2)=num(1)-1
               num(3)=num(1)-2
               num(4)=(ip-1)*2*(2*nye+1)+2*nye+2*iq
               num(5)=ip*2*(2*nye+1)+2*iq-1
               num(6)=num(5)+1
               num(7)=num(5)+2
               num(8)=num(4)+2
               num(9)=num(4)+1
             END IF
        !
             coord(1:3,1)=x_coords(ip)
             coord(5:7,1)=x_coords(ip+1)
             coord(4,1)=pt5*(coord(3,1)+coord(5,1))
             coord(8,1)=pt5*(coord(7,1)+coord(1,1))
             coord(1,2)=y_coords(iq+1)
             coord(7:8,2)=y_coords(iq+1)
             coord(3:5,2)=y_coords(iq)
             coord(2,2)=pt5*(coord(1,2)+coord(3,2))
             coord(6,2)=pt5*(coord(5,2)+coord(7,2))
             coord(9,:)=pt5*(coord(4,:)+coord(8,:))
           CASE DEFAULT
             WRITE(11,'(a)')"Wrong number of nodes for quadrilateral element"
             STOP
           END SELECT
         END IF
        RETURN
        END SUBROUTINE geom_rect