!-----------------------------------------------------!
!          Racking And Monotonic Analysis             !
!                 of CFS wall panlels                 !
!_____________________________________________________!

PROGRAM RAMA
    USE Hanuman
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)
 INTEGER::fixed_freedoms,i,iel,k,loaded_nodes,ndim,ndof,nels,neq,nod=2, &
   nodof,nn,nprops=1,np_types,nr,nlen,count=0,j,eq_num=1,load_count=1
 REAL(iwp)::axial,penalty=1.0e20_iwp,zero=0.0_iwp
!-------------------Dynamic Arrays---------------!
 INTEGER,ALLOCATABLE::etype(:),e_st_v(:),g_e_st_m(:,:),g_nd_num(:,:),eqn(:,:), &
   node(:),e_nd_num(:),sense(:),bndry(:,:),res_eq(:),load(:,:),load_eq(:)
 REAL(iwp),ALLOCATABLE::action(:),e_coord(:,:),e_dis(:),g_coord(:,:),kele(:,:), &
   kstruct(:,:),k_fill(:,:),loads(:),prop(:,:),value(:),kv(:,:),ksub(:,:),load_val(:),full_load(:),load_sub(:)
 CHARACTER(LEN=15)::filename
!-----------------PROGRAM STARTS HERE-------------!
 PRINT*,"ENTER THE INPUT FILENAME"
 READ*,filename
 OPEN(UNIT=10,FILE=filename)
 OPEN(UNIT=11,FILE='RESULT.dat',STATUS='NEW')
 READ(10,*)nels,nn,ndim,np_types
 nodof=ndim
 ndof=nod*nodof
 ALLOCATE(eqn(nodof,nn),bndry(nodof,nn),load(nodof,nn),kele(ndof,ndof),e_coord(nod,ndim),g_coord(ndim,nn),    &
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
 ALLOCATE(kstruct(neq,neq),kv(neq,neq),full_load(neq))

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
    CALL truss(kele,prop(1,etype(iel)),e_coord)
    kv=0
    e_st_v=g_e_st_m(:,iel)
    CALL assemble(kv,kele,e_st_v)
    kstruct=kv+kstruct
  END DO elements_2
  write(unit=11,fmt=245)kstruct
  245 format(12f12.3)
 ! write(11,*)" "
  
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
  print*,load_eq
  print*,load_val
  full_load=0
  do i=1,load_count-1
    full_load(load_eq(i))=load_val(i)
  end do
  print*,full_load
 
 CALL submat(neq,count,res_eq,ksub,kstruct,full_load,load_sub)
 print*,"resched here"
 
 write(unit=11,fmt=255)ksub
 255 format(9f12.3)
! write(11,*)" "
 write(unit=11,fmt=265)load_sub
 265 format(9f12.3)
!----------------------------SOLUTION PROCEDURE----------------------------------




END PROGRAM RAMA

