!-----------------------------------------------------!
!          Racking And Monotonic Analysis             !
!                 of CFS wall panlels                 !
!_____________________________________________________!

PROGRAM RAMA
    USE Hanuman
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
 INTEGER::fixed_freedoms,i,iel,k,loaded_nodes,ndim,ndof,nels,neq,nod=2, &
   nodof,nn,nprops=1,np_types,nr,nlen,count=0,j
 REAL(iwp)::axial,penalty=1.0e20_iwp,zero=0.0_iwp
!-------------------Dynamic Arrays---------------!
 INTEGER,ALLOCATABLE::etype(:),e_st_v(:),g_e_st_m(:,:),g_nd_num(:,:),eqn(:,:), &
   node(:),e_nd_num(:),sense(:),bndry(:,:),res_eq(:)
 REAL(iwp),ALLOCATABLE::action(:),e_coord(:,:),e_dis(:),g_coord(:,:),kele(:,:), &
   kstruct(:,:),k_fill(:,:),loads(:),prop(:,:),value(:),kv(:,:),ksub(:,:)
 CHARACTER(LEN=15)::filename
!-----------------PROGRAM STARTS HERE-------------!
 PRINT*,"ENTER THE INPUT FILENAME"
 READ*,filename
 OPEN(UNIT=10,FILE=filename)
 OPEN(UNIT=11,FILE='RESULT.dat',STATUS='NEW')
 READ(10,*)nels,nn,ndim,np_types
 nodof=ndim
 ndof=nod*nodof
 ALLOCATE(eqn(nodof,nn),bndry(nodof,nn),kele(ndof,ndof),e_coord(nod,ndim),g_coord(ndim,nn),    &
   e_dis(ndof),action(ndof),g_nd_num(nod,nels),e_nd_num(nod),e_st_v(ndof),g_e_st_m(ndof,nels),&
   etype(nels),prop(nprops,np_types))
 READ(10,*)prop
 etype=1
 IF(np_types>1)READ(10,*)etype
 READ(10,*)g_coord
 
 READ(10,*)g_nd_num 
 
 eqn=1
 CALL num_eq(eqn)
 
 neq=MAXVAL(eqn)
 ALLOCATE(kstruct(neq,neq),kv(neq,neq))
!-----------------------Loop over elements to form global steering matrix--------------!
 elements_1: DO iel=1,nels
   e_nd_num=g_nd_num(:,iel)
   CALL end_to_est(e_nd_num,eqn,e_st_v)
   g_e_st_m(:,iel)=e_st_v   
  END DO elements_1

!-----------------------global stiffness matrix assembly------------------
  kstruct=0
  elements_2: DO iel=1,nels
    e_nd_num=g_nd_num(:,iel)
    e_coord=TRANSPOSE(g_coord(:,e_nd_num))
    
    245 format(1x,2f6.2)
    CALL truss(kele,prop(1,etype(iel)),e_coord)
    kv=0
    e_st_v=g_e_st_m(:,iel)
    CALL assemble(kv,kele,e_st_v)
    kstruct=kv+kstruct
  END DO elements_2
  
!----------------------Load and boundary conditions----------------------
  bndry=eqn
  READ(10,*)nr,(k,bndry(:,k),i=1,nr)  

  do i=1,ndim
    do j=1,nn
      IF(bndry(i,j)==0)count=count+1
    end do
  end do
  ALLOCATE(ksub(neq-count,neq-count))
  ALLOCATE(res_eq(count))
  count=0
  do i=1,ndim
    do j=1,nn
      IF(bndry(i,j)==0) THEN
      count=count+1
      res_eq(count)=eqn(i,j)
      ENDIF
    end do
  end do
  
  print*,res_eq

 
 

END PROGRAM RAMA

