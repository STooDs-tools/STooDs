module STOODS_formulas
!~**********************************************************************
!~* Purpose: Load and handle all formulas defined in STOODS model.
!~*          This is an interface between TextFile_model and STOODS_types
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon, Univerity of Adelaide
!~**********************************************************************
!~* Created: 08/07/2019
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~*     1. XXX
!~**********************************************************************
!~* Quick description of public procedures:
!~*    1. XXX
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use STOODS_types ! STOODS objects and methods

implicit none
Private
Public::loadFormulas,getParentPar,getHyperPar

type::indices
    integer(mik),allocatable::indx(:)
end type indices

type::formulas
    integer(mik)::nIN=undefIN ! number of inputs (covariates, processes dim-distance and dim-covariates in stoods model)
    integer(mik)::nPAR=undefIN ! number of parameters (parameters and hyper-parameters in stoods model)
    integer(mik)::nFORM=undefIN ! total number of formulas
    character(250),allocatable::nameIN(:) ! size nIN, name of inputs
    character(250),allocatable::namePAR(:) ! size nPAR, name of parameters
    character(fmaxlength),allocatable::f(:) ! size nFORM, all formulas
    type(indices),allocatable::ifParentPar(:) ! size model%nVar. ifParentPar(k)%indx contains indices of formulas to be applied for retrieving parent parameters of the kth variable
    type(indices),allocatable::ifHyperPar(:) ! size model%nPro. ifHyperPar(k)%indx contains indices of formulas to be applied for retrieving hyperparameters of the kth process
    integer(mik),allocatable::iPar(:) ! size model%nPar, indices of parameters in namePAR
    integer(mik),allocatable::iCov(:) ! size model%nCov, indices for parameters in nameIN
    integer(mik),allocatable::iPro(:) ! size model%nPro, indices for processes in nameIN
    type(indices),allocatable::iHyperPar(:) ! size model%nPro. iHyperPar(k)%indx contains indices of hyperparameters of process k in namePAR
    type(indices),allocatable::iHyperCov(:) ! size model%nDim. iHyperCov(k)%indx contains indices of covariates of dimension k in nameIN
    integer(mik),allocatable::iDist(:) ! size model%nDim, indices for distance in nameIN
end type formulas
type(formulas)::FF

character(250),parameter::emptyFormulaString='-6666'

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine loadFormulas(model,err,mess)
!^**********************************************************************
!^* Purpose: load all formulas contained in model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 08/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.model, stoods model object
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use TextFile_model,only:TxtMdl_define
type(stoodsModel),intent(in)::model
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='loadFormulas'
integer(mik)::i,j,k
character(250)::parse(2)

err=0;mess=''
!==============================================================
! Compute sizes and allocate
FF%nPAR=model%nPar+sum(model%process(:)%hyper%nPar)
FF%nIN=model%nCov+model%nPro+model%nDim+sum(model%dims(:)%k)
FF%nFORM=sum(model%nParentPar)+hyperdistNformula*model%nPro
if(allocated(FF%nameIN)) deallocate(FF%nameIN);allocate(FF%nameIN(FF%nIN))
if(allocated(FF%namePAR)) deallocate(FF%namePAR);allocate(FF%namePAR(FF%nPAR))
if(allocated(FF%f)) deallocate(FF%f);allocate(FF%f(FF%nFORM))
if(allocated(FF%ifParentPar)) deallocate(FF%ifParentPar);allocate(FF%ifParentPar(model%nVar))
if(allocated(FF%ifHyperPar)) deallocate(FF%ifHyperPar);allocate(FF%ifHyperPar(model%nPro))
if(allocated(FF%iHyperPar)) deallocate(FF%iHyperPar);allocate(FF%iHyperPar(model%nPro))
if(allocated(FF%iHyperCov)) deallocate(FF%iHyperCov);allocate(FF%iHyperCov(model%nDim))
if(allocated(FF%iPar)) deallocate(FF%iPar);allocate(FF%iPar(model%nPar))
if(allocated(FF%iCov)) deallocate(FF%iCov);allocate(FF%iCov(model%nCov))
if(allocated(FF%iPro)) deallocate(FF%iPro);allocate(FF%iPro(model%nPro))
if(allocated(FF%iDist)) deallocate(FF%iDist);allocate(FF%iDist(model%nDim))

!==============================================================
! populate namePAR
k=0
! parameters
do j=1,model%nPar
    FF%namePAR(k+j) = model%par(j)%name
    FF%iPar(j)=k+j
enddo
k=k+model%nPar
! processes
do i=1,model%nPro
    if(allocated(FF%iHyperPar(i)%indx)) deallocate(FF%iHyperPar(i)%indx)
    allocate(FF%iHyperPar(i)%indx(model%process(i)%hyper%nPar))
    do j=1,model%process(i)%hyper%nPar ! hyper-parameters
        FF%namePAR(k+j)=trim(model%process(i)%hyper%par(j)%name)
        FF%iHyperPar(i)%indx(j)=k+j
    enddo
    k=k+model%process(i)%hyper%nPar
enddo

!==============================================================
! populate nameIN
k=0
! covariates
do j=1,model%nCov
    FF%nameIN(k+j) = model%covName(j)
    FF%iCov(j)=k+j
enddo
k=k+model%nCov
! processes
do j=1,model%nPro
    FF%nameIN(k+j) = model%process(j)%name
    FF%iPro(j)=k+j
enddo
k=k+model%nPro
! dimensions: distance and covariates
do i=1,model%nDim
    if(allocated(FF%iHyperCov(i)%indx)) deallocate(FF%iHyperCov(i)%indx)
    allocate(FF%iHyperCov(i)%indx(model%dims(i)%k))
    FF%nameIN(k+1)=trim(distance_str)//trim(model%dims(i)%name)
    FF%iDist(i)=k+1
    k=k+1
    do j=1,model%dims(i)%k ! covariates of the dimension
        FF%nameIN(k+j)=trim(model%dims(i)%Zname(j))
        FF%iHyperCov(i)%indx(j)=k+j
    enddo
    k=k+model%dims(i)%k
enddo

! check for duplicates
if(isAnyDuplicate( (/FF%nameIN,FF%namePAR/) )) then
    parse(1)=trim(procname);parse(2)='variable'
    call warningDuplicate(strings=(/FF%nameIN,FF%namePAR/),parse=parse)
    err=1;mess=trim(procname)//':duplicate variable names';return
endif

! populate formulas
k=0
if(sum(model%nParentPar)>0) FF%f( (k+1) : (k+sum(model%nParentPar)) ) = model%f(:) ! formulas for parent parameters
k=k+sum(model%nParentPar)
do i=1,model%nPro ! mean/covariance formulas for hyperdist of each process
    FF%f(k+1)=model%process(i)%hyper%fMu
    FF%f(k+2)=model%process(i)%hyper%fSigma
    ! change empty formulas to avoid a bug when loading through call to TxtMdl_define()
    if(len_trim(FF%f(k+1))==0) FF%f(k+1)=emptyFormulaString
    if(len_trim(FF%f(k+2))==0) FF%f(k+2)=emptyFormulaString
    k=k+2
enddo

! call Text File model to load formulas
call TxtMdl_define(nameIN=FF%nameIN,namePAR=FF%namePAR,&
                   formulas=FF%f,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

! Compute formula indices
k=0
do i=1,model%nVar
    if(allocated(FF%ifParentPar(i)%indx)) deallocate(FF%ifParentPar(i)%indx)
    allocate(FF%ifParentPar(i)%indx(model%nParentPar(i)))
    do j=1,model%nParentPar(i)
        FF%ifParentPar(i)%indx(j)=k+j
    enddo
    k=k+model%nParentPar(i)
enddo

do i=1,model%nPro
    if(allocated(FF%ifHyperPar(i)%indx)) deallocate(FF%ifHyperPar(i)%indx)
    allocate(FF%ifHyperPar(i)%indx(hyperdistNformula))
    do j=1,hyperdistNformula
        FF%ifHyperPar(i)%indx(j)=k+j
    enddo
    k=k+hyperdistNformula
enddo

end subroutine loadFormulas
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getParentPar(model,i,parentPar,feas,err,mess)
!^**********************************************************************
!^* Purpose: get parent parameters associated with the ith data y(i)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 08/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.model, stoods model object
!^*    2.i, data index
!^* OUT
!^*    1.parentPar, parent parameters
!^*    2.feas, is feasible?
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
use TextFile_model,only:TxtMdl_Apply
type(stoodsModel),intent(in)::model
integer(mik),intent(in)::i
real(mrk),intent(out)::parentPar(:)
logical,intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='getParentPar'
real(mrk)::IN(1,FF%nIN),theta(FF%nPAR),OUT(1,size(parentPar))
integer(mik)::j,k,p,ivar
logical::vfeas(1)

err=0;mess='';feas=.true.;parentPar=undefRN
IN=undefRN;theta=undefRN
theta(FF%iPar)=model%par(:)%val ! parent parameters
IN(1,FF%iCov)=model%dat%x(i,:) ! covariates
do j=1,model%nPro
    k=model%proDimIndx(j) ! This process is associated with the kth dimension
    p=model%dat%idims(i,k) ! y(i) is associated with the pth value of this dimension (e.g. the pth site or time step)
    if(p <= 0) then ! this process is not active for data y(i)
        IN(1,FF%iPro(j))=undefRN
    else
        IN(1,FF%iPro(j))=model%process(j)%val(p) ! process value
    endif
enddo

! apply formula
ivar=model%dat%ivar(i) ! y(i) data is for variable ivar
do j=1,model%nParentPar(ivar) ! all parent parameters for this variable
    k=FF%ifParentPar(ivar)%indx(j) ! the kth formula has to be applied for this parent par
    call TxtMdl_Apply(IN=IN,theta=theta,whichOUT=k,OUT=OUT,&
                      feas=vfeas,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    feas=vfeas(1)
    if(.not.feas) return
    parentPar(j)=OUT(1,1)
enddo
end subroutine getParentPar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getHyperPar(model,i,feas,err,mess)
!^**********************************************************************
!^* Purpose: get hyperparameters associated with a process
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 10/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* INOUT
!^*    1.model, stoods model object, to be updated
!^* IN
!^*    2.i, process index
!^* OUT
!^*    1.feas, is feasible?
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use TextFile_model,only:TxtMdl_Apply
use MultivariateDistribution_tools,only:initPar,NNGP_getNN
use utilities_dmsl_kit,only:number_string
! arguments
type(stoodsModel),intent(inout)::model
integer(mik),intent(in)::i
logical,intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='getHyperPar'
type(stoodsProcess)::pro
type(stoodsDimension)::dim
real(mrk)::IN(1,FF%nIN),theta(FF%nPAR),OUT(1,1)
integer(mik)::j,k,m,idim
logical::vfeas(1)
character(250)::hyperdist

err=0;mess='';feas=.true.
IN=undefRN;theta=undefRN
pro=model%process(i) ! process
idim=model%proDimIndx(i) ! This process is associated with the (idim)th dimension
dim=model%dims(idim) ! associated dimension
hyperdist=model%process(i)%hyper%name ! hyperdist name

!==============================================================
! initialise parameter object
call initPar(DistId=hyperdist,p=dim%n,&
             par=model%process(i)%hyper%mdpar,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

!==============================================================
! scalar parameters
if(pro%hyper%nSpar>0) then
    if(pro%hyper%nSpar /= size(model%process(i)%hyper%mdpar%spar)) then
        mess=trim(procname)//': nSpar ['//trim(number_string(pro%hyper%nSpar))//']'//&
             ' is incorrect for hyper-distribution ['//trim(hyperdist)//']'
        err=1;return
    endif
    model%process(i)%hyper%mdpar%spar=pro%hyper%par(pro%hyper%iSpar)%val
endif

!==============================================================
! vector parameter (mean)
if(len_trim(pro%hyper%fmu)>0) then
    theta(FF%iHyperPar(i)%indx)= pro%hyper%par(:)%val ! hyperparameters
    k=FF%ifHyperPar(i)%indx(1) ! the kth formula has to be applied for vector parameter
    do j=1,dim%n
        IN(1,FF%iHyperCov(idim)%indx)=dim%z(j,:) ! hyper-covariates
        call TxtMdl_Apply(IN=IN,theta=theta,whichOUT=k,OUT=OUT,&
                          feas=vfeas,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        feas=vfeas(1)
        if(.not.feas) return
        model%process(i)%hyper%mdpar%vpar(j)=OUT(1,1)
    enddo
endif

!==============================================================
! matrix parameter (covariance)
if(len_trim(pro%hyper%fsigma)>0) then
    k=FF%ifHyperPar(i)%indx(2) ! the kth formula has to be applied for vector parameter
    do j=1,dim%n
        IN(1,FF%iHyperCov(idim)%indx)=dim%z(j,:) ! hyper-covariates
        do m=1,j
            IN(1,FF%iDist(idim))=dim%D(j,m) ! distance
            call TxtMdl_Apply(IN=IN,theta=theta,whichOUT=k,OUT=OUT,&
                              feas=vfeas,err=err,mess=mess)
            if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
            feas=vfeas(1)
            if(.not.feas) return
            model%process(i)%hyper%mdpar%mpar(j,m)=OUT(1,1)
            model%process(i)%hyper%mdpar%mpar(m,j)=OUT(1,1)
        enddo
    enddo
    model%process(i)%hyper%mdpar%isCov=.true.
    model%process(i)%hyper%mdpar%saveInvert=.true.
endif

!==============================================================
! NNGP: get nearest parents
! 2DO: this is not efficient, since nearest parents should never change
!      monitor time spent here, and if needed :
!      (i) do more efficient sorting
!      (ii) re-organize stuff to avoid re-assigning nearest parents (using SAMI?)
if(trim(hyperdist)=='NNGP') then
    call NNGP_getNN(k=nint(model%process(i)%hyper%mdpar%spar(1)),&
                    sigma=model%process(i)%hyper%mdpar%mpar,&
                    Pa=model%process(i)%hyper%mdpar%i1,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
endif


end subroutine getHyperPar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module STOODS_formulas
