module STOODS_types
!~**********************************************************************
!~* Purpose: Definition of types for models in Space, Time Or Other
!~*          DimensionS (STOODS) & methods to manipulate them
!~*          (specify from config files, fold & unfold, etc.)
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon, University of Adelaide
!~**********************************************************************
!~* Created: 24/06/2019
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~*     1. XXX
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. stoodsModel, object that contains everything (data, model, priors, etc.)
!~*     2. Config_Model_Read, read STOODS model from config file
!~*     3. fold (vector to model) and unfold (model to vector)
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use MultivariateDistribution_tools,only:mDistPar ! parameter object for a multivariate distribution

implicit none
Private
Public::Config_Model_Read,fold,unfold,unfoldNames,applyConstraints,&
        isAnyDuplicate,warningDuplicate

integer(mik),parameter,public:: fmaxlength=2500 ! maximum length for stoods' model formulae
integer(mik),parameter,public:: hyperdistNformula=2 ! number of formulas for hyperdist (2: mean and covariance matrix)
character(1),parameter,public::distance_str='D' ! string prefix for distance (e.g. Dspace for distance in 'space' dimension)
character(250),parameter,public::FIX_str='FIX' ! string for a fixed (non-inferred) parameter
logical, parameter::printWarning=.true.  ! .false. to disable all console warnings
real(mrk),parameter::jump_priorRatio=0.01_mrk ! initial jump is jump_priorRatio * 0.5 * prior interquantile
real(mrk),parameter::jump_valRatio=0.01_mrk   ! if prior interquantile can't be computed, initial jump is jump_valRatio * parameter value
real(mrk),parameter::jump_default=0.1_mrk    ! if nothing works, initial jump is jump_default
real(mrk),parameter::jump_prob=0.84_mrk      ! probability used to define interquantile
character(250),parameter::initFile_val='overrideVal.init'      ! name of the file to override initial values (useful to continue a MCMC run)
character(250),parameter::initFile_jump='overrideJump.init'    ! name of the file to override initial jumps (useful to continue a MCMC run)

! Parameter
type,public::stoodsParam
    character(250)::name=''                     ! name of the parameter
    real(mrk)::val=undefRN                      ! value
    character(250)::priorDist=''                ! prior distribution
    real(mrk), allocatable::priorPar(:)         ! prior parameters
    real(mrk)::jump=undefRN                     ! MCMC: jump size
    real(mrk)::move=undefRN                     ! MCMC: move rate (or count)
end type stoodsParam

! Distance
type:: distance
    character(250)::funk='' ! distance function
    integer(mik)::nPar=undefIN ! number of distance parameters
    type(stoodsParam), allocatable::par(:) ! size nPar: distance parameters
end type distance

! Dimension
type,public:: stoodsDimension
    character(250)::name='' ! 'space', 'time', 'duration', etc.
    character(250)::file='' ! file where the dataset is stored
    integer(mik)::p=undefIN ! coordinates are defined in R^p (1 for time, 2 for lon-lat space, 3 for lon-lat-elevation space, etc.)
    integer(mik)::n=undefIN ! number of points
    real(mrk), allocatable::coord(:,:) ! n*p: coordinates
    type(distance)::dist ! properties of distance calculation
    real(mrk), allocatable::D(:,:) ! n*n: pairwise distance matrix
    integer(mik)::k=undefIN ! number of covariates
    character(250),allocatable::zName(:) ! size k: name of the covariates
    real(mrk), allocatable::z(:,:) ! n*k: covariates
end type stoodsDimension

! Hyper-distribution
type::hyperDist
    character(250)::name='' ! name of the hyper-distribution
    integer(mik)::nPar=undefIN ! total number of parameters needed to evaluate the hyper-pdf
    type(stoodsParam),allocatable::par(:) ! size nPar: all parameters
    integer(mik)::nSpar=undefIN ! number of scalar parameters of the hyper-pdf (see MultivariateDistribution_tools)
    integer(mik),allocatable::iSpar(:) ! size nSpar: indices of scalar parameters in par(:)
    character(fmaxlength)::fMu='' ! formula to retrieve the mean vector (vpar in MultivariateDistribution_tools) from parameters and covariates
    character(fmaxlength)::fSigma='' ! formula to retrieve the pairwise dependence matrix (mpar in MultivariateDistribution_tools) from parameters, covariates and pairwise distances
    type(mDistPar)::mdpar ! 'mDistPar' object in MultivariateDistribution_tools
end type hyperDist

! Stochastic process
type,public:: stoodsProcess
    character(250)::name=''         ! e.g. 'location parameter in space', 'hidden climate index', etc.
    character(250)::dimName=''      ! Nname of the dimension the process is associated with.
    real(mrk)::initVal=undefRN      ! initial value
    character(250)::initCovName=''  ! dimension-covariate used to set initial values
    real(mrk), allocatable::val(:)  ! size dim%n: values
    real(mrk), allocatable::jump(:) ! size dim%n: MCMC: jump size
    real(mrk), allocatable::move(:) ! size dim%n: MCMC: move rate (or count)
    integer(mik)::constraint=undefIN ! constraint index. 0 = none, 1 = centred, 2 = centred & scaled.
    type(hyperDist)::hyper ! finite-dimensional hyperdistribution of the process
end type stoodsProcess

! stoods Data
type::stoodsData
    character(250)::name='' ! name of the dataset
    character(250)::file='' ! file where the dataset is stored
    integer(mik)::n=undefIN ! number of values
    ! target variable (predictand)
    real(mrk),allocatable::y(:) ! size n: values taken by the target variable
    ! multivariate case - each variable may have a different distribution
    integer(mik),allocatable::ivar(:) ! size n: variable index
    ! covariates (predictors)
    real(mrk),allocatable::x(:,:) ! size n*nCov: values taken by the covariates
    ! dimension indices
    integer(mik),allocatable::idims(:,:) ! size n*nDim: dimension indices
    ! censoring
    integer(mik),allocatable::cType(:) ! size n: censoring type; <0:less than, >0:more than, 0:none or interval, see below.
    real(mrk),allocatable::cWidth(:) ! size n: half-width of censoring interval. Example: (y=2.1,cWidth=0.3) <=> y in [1.8;2.4]
end type stoodsData

! stoods model
type,public::stoodsModel
    character(250)::name='' ! name of the model
    character(250)::workspace='' ! workspace where configuration files are stored
    ! data layer
    integer(mik)::nVar=undefIN ! number of variables
    character(250),allocatable::varName(:) ! size nVar: variable names
    character(250),allocatable::parentDist(:) ! size nVar: parent distribution of data y (e.g. Bernoulli, Gaussian, GEV)
    integer(mik),allocatable::nParentPar(:)   ! size nVar: number of parent parameters (e.g. 1 for Bernoulli, 2 for Gaussian, 3 for GEV)
    character(250),allocatable::parentParName(:) ! size sum(nDpar): name of each parent parameter (e.g. location, scale, shape)
    type(stoodsData)::dat
    ! covariates (predictors)
    integer(mik)::nCov=undefIN ! number of covariates
    character(250),allocatable::covName(:) ! size nCov: name of the covariates
    ! parameters
    integer(mik)::nPar=undefIN ! number of parameters
    type(stoodsParam),allocatable::par(:) ! size nPar: parameters
    ! dimensions
    integer(mik)::nDim=undefIN ! number of dimensions
    type(stoodsDimension),allocatable::dims(:) ! size nDim: dimensions
    ! processes
    integer(mik)::nPro=undefIN ! number of processes
    type(stoodsProcess),allocatable::process(:) ! size nPro: processes
    character(250),allocatable::proDim(:) ! size nPro: name of the dimension associated with each process
    integer(mik),allocatable::proDimIndx(:)  ! size nPro: index of the dimension associated with each process
    ! formulas
    character(fmaxlength),allocatable::f(:) ! size sum(nParentPar): formula to retrieve each parent parameters from parameters, covariates and processes
    ! misc.
    integer(mik)::nVec=undefIN ! number of elements when the model is unfolded into a vector
end type stoodsModel

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function getUnfoldedSize(model)
!^**********************************************************************
!^* Purpose: Compute numbers of parameters when model is unfolded into vector
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 03/07/2019
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
!^*    1.getUnfoldedSize, total number of parameters when model is unfolded into vector
!^**********************************************************************
type(stoodsModel),intent(in)::model
integer(mik)::getUnfoldedSize
! locals
character(250),parameter::procname='getUnfoldedSize'
integer(mik)::i,j,nVec

nVec=0
! Parameters
nVec=nVec+model%nPar
! Processes
do i=1,model%nPro
    ! process values
    nVec=nVec+size(model%process(i)%val)
    ! hyperparameters
    nVec=nVec+model%process(i)%hyper%nPar
enddo
! Dimensions
do i=1,model%nDim
    nVec=nVec+model%dims(i)%dist%nPar
enddo
getUnfoldedSize=nVec

end function getUnfoldedSize
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine fold(vector,model,err,mess,what)
!^**********************************************************************
!^* Purpose: fold a vector into a model object
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 08/07/2019
!^**********************************************************************
!^* Comments: model should already be specified (arrays allocated etc.)
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.vector of numerical values
!^*    2.[OPTIONAL] what, 'val' [default] or 'jump'
!^* INOUT
!^*    1.model, stoods model object
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
real(mrk),intent(in)::vector(:)
type(stoodsModel),intent(inout)::model
character(*),optional,intent(in)::what
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='fold'
integer(mik)::k,i,n
character(4)::w

err=0;mess=''
if(present(what)) then;w=what;else;w='val ';endif

k=0
! Parameters
n=model%nPar
if(n>0) then
    select case(trim(w))
    case('val')
      model%par(:)%val = vector( (k+1) : (k+n) )
    case('jump')
      model%par(:)%jump = vector( (k+1) : (k+n) )
    end select
endif
k=k+n
! Processes
do i=1,model%nPro
    ! process values
    n=size(model%process(i)%val)
    if(n>0) then
        select case(trim(w))
        case('val')
            model%process(i)%val = vector( (k+1) : (k+n) )
        case('jump')
            model%process(i)%jump = vector( (k+1) : (k+n) )
        end select
    endif
    k=k+n
    ! hyperparameters
    n=model%process(i)%hyper%nPar
    if(n>0) then
        select case(trim(w))
        case('val')
            model%process(i)%hyper%par(:)%val = vector( (k+1) : (k+n) )
        case('jump')
            model%process(i)%hyper%par(:)%jump = vector( (k+1) : (k+n) )
        end select
    endif
    k=k+n
enddo
! Dimensions
do i=1,model%nDim
    n=model%dims(i)%dist%nPar
    if(n>0) then
        select case(trim(w))
        case('val')
            model%dims(i)%dist%par(:)%val = vector( (k+1) : (k+n) )
       case('jump')
            model%dims(i)%dist%par(:)%jump = vector( (k+1) : (k+n) )
        end select
    endif
    k=k+n
enddo

end subroutine fold
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine unfold(model,vector,err,mess,what)
!^**********************************************************************
!^* Purpose: unfold a model object into a vector
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 03/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.model, stoods model object
!^*    2.[OPTIONAL] what, 'val' [default], 'jump' or 'move'
!^* OUT
!^*    1.vector, model unfolded into vector
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
type(stoodsModel),intent(in)::model
character(*),optional,intent(in)::what
real(mrk),intent(out)::vector(model%nVec)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='unfold'
integer(mik)::k,i,n
character(4)::w

err=0;mess='';vector=undefRN
if(present(what)) then;w=what;else;w='val ';endif
k=0
! Parameters
n=model%nPar
if(n>0) then
    select case(trim(w))
    case('val')
        vector( (k+1) : (k+n) ) = model%par(:)%val
    case('move')
        vector( (k+1) : (k+n) ) = model%par(:)%move
    case('jump')
        vector( (k+1) : (k+n) ) = model%par(:)%jump
    end select
endif
k=k+n
! Processes
do i=1,model%nPro
    ! process values
    n=size(model%process(i)%val)
    if(n>0) then
        select case(trim(w))
        case('val')
            vector( (k+1) : (k+n) ) = model%process(i)%val
        case('move')
            vector( (k+1) : (k+n) ) = model%process(i)%move
        case('jump')
            vector( (k+1) : (k+n) ) = model%process(i)%jump
        end select
    endif
    k=k+n
    ! hyperparameters
    n=model%process(i)%hyper%nPar
    if(n>0) then
        select case(trim(w))
        case('val')
            vector( (k+1) : (k+n) ) = model%process(i)%hyper%par(:)%val
        case('move')
            vector( (k+1) : (k+n) ) = model%process(i)%hyper%par(:)%move
        case('jump')
            vector( (k+1) : (k+n) ) = model%process(i)%hyper%par(:)%jump
        end select
    endif
    k=k+n
enddo
! Dimensions
do i=1,model%nDim
    n=model%dims(i)%dist%nPar
    if(n>0) then
        select case(trim(w))
        case('val')
            vector( (k+1) : (k+n) ) = model%dims(i)%dist%par(:)%val
        case('move')
            vector( (k+1) : (k+n) ) = model%dims(i)%dist%par(:)%move
        case('jump')
            vector( (k+1) : (k+n) ) = model%dims(i)%dist%par(:)%jump
        end select
    endif
    k=k+n
enddo

end subroutine unfold
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine unfoldNames(model,vector,err,mess)
!^**********************************************************************
!^* Purpose: unfold all names in a model object into a vector
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 12/07/2019
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
!^*    1.vector, names unfolded into vector
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:number_string
type(stoodsModel),intent(in)::model
character(*),intent(out)::vector(model%nVec)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='unfoldNames'
integer(mik)::k,i,j,n

err=0;mess='';vector=''
k=0
! Parameters
n=model%nPar
if(n>0) vector( (k+1) : (k+n) ) = model%par(:)%name
k=k+n
! Processes
do i=1,model%nPro
    ! process values
    n=size(model%process(i)%val)
    do j=1,n
        vector(k+j) = trim(model%process(i)%name)//'_'//trim(number_string(j))
    enddo
    k=k+n
    ! hyperparameters
    n=model%process(i)%hyper%nPar
    if(n>0) vector( (k+1) : (k+n) ) = model%process(i)%hyper%par(:)%name
    k=k+n
enddo
! Dimensions
do i=1,model%nDim
    n=model%dims(i)%dist%nPar
    if(n>0) vector( (k+1) : (k+n) ) = model%dims(i)%dist%par(:)%name
    k=k+n
enddo

end subroutine unfoldNames
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Config_Model_Read(workspace,file,model,err,mess)
!^**********************************************************************
!^* Purpose: Read config file for model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 28/06/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: improve initialisation of jump size for processes
!^**********************************************************************
!^* IN
!^*    1.workspace
!^*    2.file
!^* OUT
!^*    1.model, STOODS model object
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::workspace,file
type(stoodsModel),intent(out)::model
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_Model_Read'
integer(mik)::unt,i,indx,nVec,nInf
character(250)::cfile,name,parse(2)
character(250),allocatable::cfiles(:),namesALL(:)
real(mrk),allocatable::vec(:)
logical:: feas,found

err=0;mess=''
!==============================================================
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(workspace)//trim(file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
!==============================================================
! basics
model%workspace=workspace
read(unt,*,iostat=err) model%name
if(err/=0) then;mess=trim(procname)//':read error';return;endif
!==============================================================
! Variables
read(unt,*,iostat=err) model%nVar
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(model%nVar<=0) then;err=1;mess=trim(procname)//':nVar<=0';return;endif
if(allocated(model%varName)) deallocate(model%varName);allocate(model%varName(model%nVar))
read(unt,*,iostat=err) model%varName
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! Check names are all different, warning otherwise
if(isAnyDuplicate(model%varName))then
    parse(1)=trim(procname);parse(2)='variable'
    call warningDuplicate(strings=model%varName,parse=parse)
endif
if(allocated(model%parentDist)) deallocate (model%parentDist);allocate(model%parentDist(model%nVar))
read(unt,*,iostat=err) model%parentDist
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(allocated(model%nParentPar)) deallocate(model%nParentPar);allocate(model%nParentPar(model%nVar))
read(unt,*,iostat=err) model%nParentPar
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(allocated(model%parentParName)) deallocate(model%parentParName);allocate(model%parentParName(sum(model%nParentPar)))
read(unt,*,iostat=err) model%parentParName
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! Check names are all different, warning otherwise
if(isAnyDuplicate(model%parentParName))then
    parse(1)=trim(procname);parse(2)='parent parameter'
    call warningDuplicate(strings=model%parentParName,parse=parse)
endif
!==============================================================
! Covariates
read(unt,*,iostat=err) model%nCov
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(model%nCov<0) then;err=1;mess=trim(procname)//':nCov<0';return;endif
if(allocated(model%covName)) deallocate(model%covName);allocate(model%covName(model%nCov))
if(model%nCov>0) then
    read(unt,*,iostat=err) model%covName
else
    read(unt,*,iostat=err) ! model%covName
endif
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! Check names are all different, fatal error otherwise
if(isAnyDuplicate(model%covName))then
    parse(1)=trim(procname);parse(2)='covariate'
    call warningDuplicate(strings=model%covName,parse=parse)
    err=1;mess=trim(procname)//':duplicate covariate names';return
endif
!==============================================================
! Parameters
read(unt,*,iostat=err) model%nPar
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(model%nPar<0) then;err=1;mess=trim(procname)//':nPar<0';return;endif
if(allocated(model%par)) deallocate(model%par);allocate(model%par(model%nPar))
if(model%nPar>0) then
    read(unt,*,iostat=err) model%par(:)%name
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
    ! Check names are all different, fatal error otherwise
    if(isAnyDuplicate(model%par(:)%name))then
        parse(1)=trim(procname);parse(2)='parameter'
        call warningDuplicate(strings=model%par(:)%name,parse=parse)
        err=1;mess=trim(procname)//':duplicate parameter names';return
    endif
    read(unt,*,iostat=err) cfile
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
    call Config_Parameter_Read(trim(workspace)//trim(cfile),model,err,mess)
    if(err/=0) then;mess=trim(procname)//':read parameters error:'//trim(mess);return;endif
else
    read(unt,*,iostat=err) ! model%par(:)%name
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
    read(unt,*,iostat=err) ! cfile
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
endif
!==============================================================
! Dimensions
read(unt,*,iostat=err) model%nDim
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(model%nDim<0) then;err=1;mess=trim(procname)//':nDim<0';return;endif
if(allocated(model%dims)) deallocate(model%dims);allocate(model%dims(model%nDim))
if(model%nDim>0) then
    if(allocated(cfiles)) deallocate(cfiles);allocate(cfiles(model%nDim))
    read(unt,*,iostat=err) cfiles
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
    do i=1,model%nDim
        call Config_Dimension_Read(file=trim(workspace)//trim(cfiles(i)),&
                                   workspace=trim(workspace),&
                                   d=model%dims(i),err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':read dimension error:'//trim(mess);return;endif
    enddo
    ! Check names are all different, fatal error otherwise
    if(isAnyDuplicate(model%dims(:)%name))then
        parse(1)=trim(procname);parse(2)='dimension'
        call warningDuplicate(strings=model%dims(:)%name,parse=parse)
        err=1;mess=trim(procname)//':duplicate dimension names';return
    endif
else
    read(unt,*,iostat=err) ! cfiles
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
endif
!==============================================================
! Processes
read(unt,*,iostat=err) model%nPro
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(model%nPro<0) then;err=1;mess=trim(procname)//':nPro<0';return;endif
if(allocated(model%process)) deallocate(model%process);allocate(model%process(model%nPro))
if(allocated(model%proDim)) deallocate(model%proDim);allocate(model%proDim(model%nPro))
if(allocated(model%proDimIndx)) deallocate(model%proDimIndx);allocate(model%proDimIndx(model%nPro))
if(model%nPro>0) then
    read(unt,*,iostat=err) model%process(:)%name
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
    ! Check names are all different, fatal error otherwise
    if(isAnyDuplicate(model%process(:)%name))then
        parse(1)=trim(procname);parse(2)='process'
        call warningDuplicate(strings=model%process(:)%name,parse=parse)
        err=1;mess=trim(procname)//':duplicate process names';return
    endif
    if(allocated(cfiles)) deallocate(cfiles);allocate(cfiles(model%nPro))
    read(unt,*,iostat=err) cfiles
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
    do i=1,model%nPro
        ! read process configuration file
        name=model%process(i)%name
        call Config_Process_Read(file=trim(workspace)//trim(cfiles(i)),name=name,&
                                   p=model%process(i),err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':read process error:'//trim(mess);return;endif
        ! associate dimension
        model%proDim(i)=model%process(i)%dimName
        indx=findInList(txt=model%proDim(i),list=model%dims(:)%name)
        if(indx<=0) then
            err=1
            mess=trim(procname)//': process ['//trim(model%process(i)%name)//']:'//&
                 ' unrecognized dimension ['//trim(model%process(i)%dimName)//']'
            return
        else
            model%proDimIndx(i)=indx
        endif
        ! Initialize values
        if(allocated(model%process(i)%val)) deallocate(model%process(i)%val)
        allocate(model%process(i)%val( model%dims(indx)%n ))
        ! look for initCovName in available dimension-covariates
        indx=findInList(txt=model%process(i)%initCovName,list=model%dims(model%proDimIndx(i))%zname)
        if(indx>0) then ! dimension-covariate exists, use this as initial values
            model%process(i)%val=model%dims(model%proDimIndx(i))%z(:,indx)
        else ! interpret string as an initial value
            read(model%process(i)%initCovName,*,iostat=err) model%process(i)%initVal
            if(err/=0) then
                mess=trim(procname)//': process ['//trim(model%process(i)%name)//']:'//&
                     ' unable to interpret starting value'
                return
            endif
            model%process(i)%val=model%process(i)%initVal
        endif
        call applyConstraints(model%process(i),feas)
        if(.not.feas) then
            err=1
            mess=trim(procname)//': process ['//trim(model%process(i)%name)//']:'//&
                 ' unfeasible constraint'
            return
        endif
        ! initialize jump size and move rates
        indx=model%proDimIndx(i)
        if(allocated(model%process(i)%move)) deallocate(model%process(i)%move)
        allocate(model%process(i)%move( model%dims(indx)%n ))
        model%process(i)%move=0._mrk
        if(allocated(model%process(i)%jump)) deallocate(model%process(i)%jump)
        allocate(model%process(i)%jump( model%dims(indx)%n ))
        model%process(i)%jump=jump_default ! 2DO: can be improved!
    enddo
else
    read(unt,*,iostat=err) ! model%process(:)%name
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
    read(unt,*,iostat=err) ! cfiles
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
endif
!==============================================================
! Formulas
if(allocated(model%f)) deallocate(model%f);allocate(model%f(sum(model%nParentPar)))
do i=1,sum(model%nParentPar)
    read(unt,*,iostat=err) model%f(i)
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
enddo
!==============================================================
! Data
read(unt,*,iostat=err) cfile
if(err/=0) then;mess=trim(procname)//':read error';return;endif
call Config_Data_Read(file=trim(workspace)//trim(cfile),model=model,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':read data error:'//trim(mess);return;endif
!==============================================================
close(unt)
if(allocated(cfiles)) deallocate(cfiles)
!==============================================================
! Compute sizes
model%nVec=getUnfoldedSize(model)
!==============================================================
! read files for overriding default vals/jumps (if they exist)
allocate(vec(model%nVec))
! initial values
cfile=trim(workspace)//trim(initFile_val)
inquire(file=cfile,exist=found)
if(found) then
    open(unit=unt,file=cfile, status='old', iostat=err)
    if(err/=0) then;mess=trim(procname)//':open error, file: '//trim(cfile);return;endif
    read(unt,*,iostat=err) vec
    if(err/=0) then;mess=trim(procname)//':read error, file: '//trim(cfile);return;endif
    call fold(vec,model,err,mess,'val')
    if(err/=0) then;mess=trim(procname)//trim(mess);return;endif
endif
! initial jumps
cfile=trim(workspace)//trim(initFile_jump)
inquire(file=cfile,exist=found)
if(found) then
    open(unit=unt,file=cfile, status='old', iostat=err)
    if(err/=0) then;mess=trim(procname)//':open error, file: '//trim(cfile);return;endif
    read(unt,*,iostat=err) vec
    if(err/=0) then;mess=trim(procname)//':read error, file: '//trim(cfile);return;endif
    call fold(vec,model,err,mess,'jump')
    if(err/=0) then;mess=trim(procname)//trim(mess);return;endif
endif
deallocate(vec)
!==============================================================
! Further "no-duplicate" check
if(allocated(namesALL)) deallocate(namesALL);allocate(namesALL(model%nCov+model%nPro+model%nPar))
namesALL=(/model%covName(:),model%process(:)%name,model%par(:)%name/)
if(isAnyDuplicate(namesALL))then
    parse(1)=trim(procname);parse(2)='variable'
    call warningDuplicate(strings=namesALL,parse=parse)
    err=1;mess=trim(procname)//':duplicate names';return
endif
end subroutine Config_Model_Read
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Config_Data_Read(file,model,err,mess)
!^**********************************************************************
!^* Purpose: Read config file for data
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 28/06/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file
!^* INOUT
!^*    1.model, STOODS model object
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
type(stoodsModel),intent(inout)::model
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_Data_Read'
integer(mik)::unt,nh,nr,nc,yc,ivarc,cTypec,cWidthc,i
integer(mik)::xc(model%nCov),idimc(model%nDim)
real(mrk),allocatable::foo(:)

err=0;mess=''
!==============================================================
! Read config
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
read(unt,*,iostat=err) model%dat%name
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) model%dat%file
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) nh
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) nr
if(err/=0) then;mess=trim(procname)//':read error';return;endif
model%dat%n=nr
if(allocated(model%dat%y)) deallocate(model%dat%y);allocate(model%dat%y(nr))
if(allocated(model%dat%ivar)) deallocate(model%dat%ivar);allocate(model%dat%ivar(nr))
if(allocated(model%dat%x)) deallocate(model%dat%x);allocate(model%dat%x(nr,model%nCov))
if(allocated(model%dat%idims)) deallocate(model%dat%idims);allocate(model%dat%idims(nr,model%nDim))
if(allocated(model%dat%cType)) deallocate(model%dat%cType);allocate(model%dat%cType(nr))
if(allocated(model%dat%cWidth)) deallocate(model%dat%cWidth);allocate(model%dat%cWidth(nr))
read(unt,*,iostat=err) nc
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) yc
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(model%nVar>1) then
    read(unt,*,iostat=err) ivarc
else
    read(unt,*,iostat=err)
endif
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(model%nCov>0) then
    read(unt,*,iostat=err) xc
else
    read(unt,*,iostat=err)
endif
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(model%nDim>0) then
    read(unt,*,iostat=err) idimc
else
    read(unt,*,iostat=err)
endif
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) cTypec
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) cWidthc
if(err/=0) then;mess=trim(procname)//':read error';return;endif
close(unt)
!==============================================================
! Read Data
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(model%dat%file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open data error';return;endif
do i=1,nh
    read(unt,*,iostat=err)
    if(err/=0) then;mess=trim(procname)//':read data error [headers]';return;endif
enddo
if(allocated(foo)) deallocate(foo);allocate(foo(nc))
do i=1,nr
    read(unt,*,iostat=err) foo
    if(err/=0) then;mess=trim(procname)//':read data error';return;endif
    model%dat%y(i)=foo(yc)
    if(model%nVar>1) then
        model%dat%ivar(i)=foo(ivarc)
    else
        model%dat%ivar(i)=1
    endif
    if(model%nCov>0) model%dat%x(i,:)=foo(xc)
    if(model%nDim>0) model%dat%idims(i,:)=foo(idimc)
    if(cTypec==0) then
        model%dat%cType(i)=0
    else
        model%dat%cType(i)=foo(cTypec)
    endif
    if(cWidthc==0) then
        model%dat%cWidth(i)=0._mrk
    else
        model%dat%cWidth(i)=foo(cWidthc)
    endif
enddo
close(unt)
if(allocated(foo)) deallocate(foo)
!==============================================================

end subroutine Config_Data_Read
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Config_Parameter_Read(file,model,err,mess)
!^**********************************************************************
!^* Purpose: Read config file for parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 01/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file
!^* INOUT
!^*    1.model, STOODS model object
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
use Distribution_tools, only:GetParNumber
character(*), intent(in)::file
type(stoodsModel),intent(inout)::model
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_Parameter_Read'
integer(mik)::unt,i,np
character(250)::name

err=0;mess=''
!==============================================================
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
!==============================================================
do i=1,model%nPar
    name=model%par(i)%name
    call Config_OneParameter_Read(unt=unt,parname=name,par=model%par(i),err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
enddo
!==============================================================
close(unt)

end subroutine Config_Parameter_Read
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Config_OneParameter_Read(unt,parname,par,err,mess)
!^**********************************************************************
!^* Purpose: Read config file for A SINGLE parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 01/07/2019
!^**********************************************************************
!^* Comments: file already opened by calling sub
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.unt, unit of already opened file
!^*    2.[parname], par name if already existing (optional)
!^* INOUT
!^*    1.par, STOODS parameter object
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use Distribution_tools, only:GetParNumber
integer(mik),intent(in)::unt
character(*),intent(in),optional::parname
type(stoodsParam),intent(inout)::par
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_OneParameter_Read'
integer(mik)::np

err=0;mess=''
read(unt,*,iostat=err) par%name
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! warning if name differs from the already existing parname
if(present(parname)) then
    if(trim(parname)/=trim(par%name)) then
        call warningInconsistent(par%name,parname,(/procname,''/))
        par%name=parname
    endif
endif
read(unt,*,iostat=err) par%val
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) par%priorDist
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(trim(par%priorDist)==FIX_str) then
    np=0
else
    call GetParNumber(DistID=trim(par%priorDist), npar=np, err=err, mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
endif
if(allocated(par%priorPar)) deallocate(par%priorPar);allocate(par%priorPar(np))
if(np>0) then
    read(unt,*,iostat=err) par%priorPar
else
    read(unt,*,iostat=err) ! par%priorPar
endif
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! Jump size and move rate
par%jump=setJump(par%priorDist,par%priorPar,par%val)
par%move=0._mrk

end subroutine Config_OneParameter_Read
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Config_Dimension_Read(file,workspace,d,err,mess)
!^**********************************************************************
!^* Purpose: Read config file for one dimension
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 01/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, config file
!^*    2.workspace
!^* INOUT
!^*    1.d, STOODS dimension object
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
use Geodesy_tools,only:GetDistance
character(*),intent(in)::file,workspace
type(stoodsDimension),intent(inout)::d
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_Dimension_Read'
integer(mik)::unt,i,nh,nr,nc,p,k
integer(mik),allocatable::icoord(:),icov(:)
real(mrk),allocatable::foo(:)

err=0;mess=''
!==============================================================
! Read config
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
read(unt,*,iostat=err) d%name
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) d%file
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) nh
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) nr
if(err/=0) then;mess=trim(procname)//':read error';return;endif
d%n=nr
read(unt,*,iostat=err) nc
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) p
if(err/=0) then;mess=trim(procname)//':read error';return;endif
d%p=p
if(allocated(d%coord)) deallocate(d%coord);allocate(d%coord(nr,p))
if(allocated(icoord)) deallocate(icoord);allocate(icoord(p))
read(unt,*,iostat=err) icoord
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) k
if(err/=0) then;mess=trim(procname)//':read error';return;endif
d%k=k
if(allocated(d%zName)) deallocate(d%zName);allocate(d%zName(k))
if(allocated(icov)) deallocate(icov);allocate(icov(k))
if(allocated(d%z)) deallocate(d%z);allocate(d%z(nr,k))
read(unt,*,iostat=err) d%zName
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) icov
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) d%dist%funk
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) d%dist%nPar
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(allocated(d%dist%par)) deallocate(d%dist%par);allocate(d%dist%par(d%dist%nPar))
do i=1,d%dist%nPar
    call Config_OneParameter_Read(unt=unt,par=d%dist%par(i),err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
enddo
close(unt)
!==============================================================
! Read dimension data
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(d%file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open data error';return;endif
do i=1,nh
    read(unt,*,iostat=err)
    if(err/=0) then;mess=trim(procname)//':read data error [headers]';return;endif
enddo
if(allocated(foo)) deallocate(foo);allocate(foo(nc))
do i=1,nr
    read(unt,*,iostat=err) foo
    if(err/=0) then;mess=trim(procname)//':read data error';return;endif
    d%coord(i,:)=foo(icoord)
    if(d%k>0) d%z(i,:)=foo(icov)
enddo
close(unt)
if(allocated(foo)) deallocate(foo)
!==============================================================
! Compute distances
if(allocated(d%D)) deallocate(d%D);allocate(d%D(d%n,d%n))
call GetDistance(formula=d%dist%funk,pts=d%coord,&
                 par=d%dist%par(:)%val,cov=d%z,&
                 D=d%D,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
!==============================================================

end subroutine Config_Dimension_Read
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Config_Process_Read(file,name,p,err,mess)
!^**********************************************************************
!^* Purpose: Read config file for one process
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 02/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, config file
!^*    2.[name], already-existing name (optional)
!^* INOUT
!^*    1.p, STOODS process object
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*),intent(in)::file
character(*),intent(in),optional::name
type(stoodsProcess),intent(inout)::p
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_Process_Read'
integer(mik)::unt,i
character(250)::parname

err=0;mess=''
!==============================================================
! Open
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
!==============================================================
! Process name
read(unt,*,iostat=err) p%name
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! warning if name differs from the already existing parname
if(present(name)) then
    if(trim(name)/=trim(p%name)) then
        call warningInconsistent(p%name,name,(/procname,''/))
        p%name=name
    endif
endif
!==============================================================
! Associated dimension
read(unt,*,iostat=err) p%dimName
if(err/=0) then;mess=trim(procname)//':read error';return;endif
!==============================================================
! Constrained process?
read(unt,*,iostat=err) p%constraint
if(err/=0) then;mess=trim(procname)//':read error';return;endif
!==============================================================
! Initial value
read(unt,*,iostat=err) p%initCovName
if(err/=0) then;mess=trim(procname)//':read error';return;endif
!==============================================================
! Hyper-distribution
read(unt,*,iostat=err) p%hyper%name
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! All parameters needed to evaluate the hyperdist
read(unt,*,iostat=err) p%hyper%nPar
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(allocated(p%hyper%par)) deallocate(p%hyper%par);allocate(p%hyper%par(p%hyper%nPar))
if(p%hyper%nPar>0) then
    read(unt,*,iostat=err) p%hyper%par(:)%name
else
    read(unt,*,iostat=err) ! p%hyper%par(:)%name
endif
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! Scalar parameters
read(unt,*,iostat=err) p%hyper%nSpar
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(allocated(p%hyper%iSpar)) deallocate(p%hyper%iSpar);allocate(p%hyper%iSpar(p%hyper%nSpar))
if(p%hyper%nSpar>0) then
    read(unt,*,iostat=err) p%hyper%iSpar
else
    read(unt,*,iostat=err) ! p%hyper%iSpar
endif
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(any(p%hyper%iSpar<=0).or.any(p%hyper%iSpar>p%hyper%nPar)) then
    mess=trim(procname)//':iSpar:invalid index';return
endif
! vector (typically mean) and matrix (typically covariance) parameters
read(unt,*,iostat=err) p%hyper%fmu
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) p%hyper%fsigma
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! Priors
do i=1,p%hyper%nPar
    parname=p%hyper%par(i)%name
    call Config_OneParameter_Read(unt=unt,parname=parname,par=p%hyper%par(i),err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
enddo
close(unt)
!==============================================================

end subroutine Config_Process_Read
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function isAnyDuplicate(strings)
!^**********************************************************************
!^* Purpose: look for duplicates in a string vector
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 01/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.strings, string vector
!^* OUT
!^*    1.logical, is any string duplicate?
!^**********************************************************************
character(*), intent(in)::strings(:)
logical::isAnyDuplicate
! locals
integer(mik)::i,j,n

isAnyDuplicate=.false.
n=size(strings)
if(n>=2) then
    do i=2,n
        do j=1,(i-1)
            if(trim(strings(i))==trim(strings(j))) then
                isAnyDuplicate=.true.;return
            endif
        enddo
    enddo
endif

end function isAnyDuplicate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine warningDuplicate(strings,parse)
!^**********************************************************************
!^* Purpose: console warning for duplicate names
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 01/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.strings (vector): string vector where duplicate has ben detected
!^*    2.parse (vector): strings to be passed in warning message
!^**********************************************************************
character(*), intent(in)::strings(:),parse(:)
! locals
integer(mik)::i

if(.not. printWarning) return

write(*,*) '/!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ '
write(*,*) 'WARNING from subroutine: '//trim(parse(1))
write(*,*) 'Duplicate '//trim(parse(2))//' names in list below.'
do i=1,size(strings)
    write(*,*) trim(strings(i))
enddo
write(*,*) '/!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ '

end subroutine warningDuplicate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine warningInconsistent(name1,name2,parse)
!^**********************************************************************
!^* Purpose: console warning for inconsistent names
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 01/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.name1
!^*    1.name2
!^*    2.parse (vector): strings to be passed in warning message
!^**********************************************************************
character(*), intent(in)::name1,name2,parse(:)
! locals
integer(mik)::i

if(.not. printWarning) return

write(*,*) '/!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ '
write(*,*) 'WARNING from subroutine: '//trim(parse(1))
write(*,*) 'Conflicting parameter names.'
write(*,*) 'In current configuration file '//trim(parse(2))//':'
write(*,*) '==========> name = '//trim(name1)
write(*,*) 'In model configuration file:'
write(*,*) '==========> name = '//trim(name2)
write(*,*) 'Only the latter will be used.'
write(*,*) '/!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ '

end subroutine warningInconsistent
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function findInList(txt,list)
!^**********************************************************************
!^* Purpose: find first occurrence of text in a list of strings and return index
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 02/07/2019
!^**********************************************************************
!^* Comments: return 0 if txt not found in list
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.txt, text to b found
!^*    2.list, list of strings
!^* OUT
!^*    1.integer, indx of first occurrence of txt in list, 0 if not found
!^**********************************************************************
character(*), intent(in)::txt,list(:)
integer(mik)::findInList
! locals
integer(mik)::i,n

findInList=0
n=size(list)
do i=1,n
    if(trim(txt)==trim(list(i))) then
        findInList=i;return
    endif
enddo

end function findInList
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setJump(priorDist,priorPar,val)
!^**********************************************************************
!^* Purpose: set initial jump size
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 12/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.priorDist, string: prior distribution
!^*    2.priorPar, real vector: prior parameters
!^*    3.val, real: parameter initial value
!^* OUT
!^*    1.real: initial jump
!^**********************************************************************
use Distribution_tools, only:getQuantile
character(*), intent(in)::priorDist
real(mrk),intent(in)::priorPar(:),val
real(mrk)::setJump
! locals
character(250),parameter::procname='setJump'
real(mrk)::q1,q2
integer(mik)::err
character(250)::mess
logical::feas

setJump=undefRN
!==============================================================
! try computing interquantile range
call GetQuantile(DistId=trim(priorDist),p=jump_prob,par=priorPar,q=q2,feas=feas,err=err,mess=mess)
if(feas .and. err==0) then
    call GetQuantile(DistId=trim(priorDist),p=1._mrk-jump_prob,par=priorPar,q=q1,feas=feas,err=err,mess=mess)
endif
!==============================================================
! assign jump
if(feas .and. err==0) then
    setJump=jump_priorRatio*0.5_mrk*(q2-q1)
else
    if(val==0) then
        setJump = jump_default
    else
        setJump = jump_valRatio * abs(val)
    endif
endif
!==============================================================

end function setJump
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine applyConstraints(proc,feas)
!^**********************************************************************
!^* Purpose: Apply identifiability constraints
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 26/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* INOUT
!^*    1.proc: process to be constrained
!^* OUT
!^*    1.feas: feasible?
!^**********************************************************************
type(stoodsProcess),intent(inout)::proc
logical,intent(out)::feas
! locals
integer(mik)::n
real(mrk)::u,v,delta

feas=.true.
if(proc%constraint==0) return ! nothing to do
if(proc%constraint>2) then ! should never happen
    proc%val=undefRN;return
endif
n=size(proc%val)
if(proc%constraint==1) then
    proc%val(n)=-1._mrk*sum(proc%val(1:n-1))
    return
endif
if(proc%constraint==2) then
    u=sum(proc%val(1:(n-2)))
    v=dot_product(proc%val(1:(n-2)),proc%val(1:(n-2)))
    delta=2._mrk*n-2._mrk*v-u**2
    if(delta<0._mrk) then;feas=.false.;return;endif
    proc%val(n-1)=-0.5_mrk*(u+sqrt(delta))
    proc%val(n)=-0.5_mrk*(u-sqrt(delta))
endif
end subroutine applyConstraints
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module STOODS_types
