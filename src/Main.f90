program main

use kinds_dmsl_kit
use STOODS_types,only:stoodsModel,Config_Model_Read,unfold
use STOODS_formulas,only:loadFormulas
use STOODS_inference,only:getLogLkh,getLogPrior,getLogH,getLogPost
use STOODS_mcmc,only:OATsampler,MCMCTunings,Config_MCMC_Read
use uniran1_dmsl_mod,only: seed_uniran

implicit none

character(250),parameter::version='0.1.1 March 2021'
character(250),parameter::wkFile='Config.txt'
character(250),parameter::modelFile='model.config'
character(250),parameter::mcmcFile='mcmc.config'
type(stoodsModel)::model
integer(mik)::err,i,narg,seed
character(250)::mess
character(2500)::wk,workspace,arg
type(MCMCTunings)::mcmc

err=0;mess='ok'

! command line arguments
i=1;wk=''
narg=command_argument_count()
do while (i<=narg)
     call get_command_argument(i,arg)
     select case (arg)
     case ('-wk', '--workspace')
        i=i+1
        if(i<=narg) then
            call get_command_argument(i,arg)
            wk=trim(arg);i=i+1
        else
            call fatalExit('-wk requires a path to a directory')
        endif
     case ('-sd', '--seed')
        i=i+1
        if(i<=narg) then
            call get_command_argument(i,arg)
            read(arg,*,iostat=err) seed
            if(err==0) then
                call seed_uniran(put=seed)
                i=i+1
            else
                call fatalExit('-sd requires the seed as an integer number')
            endif
        else
            call fatalExit('-sd requires the seed as an integer number')
        endif
     case ('-rd', '--random')
        call seed_uniran(CPUtime=.true.)
        i=i+1
    case ('-v', '--version')
        write(*,*) 'version: ', trim(version)
        STOP
     case ('-h', '--help')
        call print_help()
        STOP
     case default
        write(*,*) 'Unrecognized command-line option: ', trim(arg)
        call print_help()
        call fatalExit('')
     end select
  end do

! Get Workspace
if (wk=='') then ! workspace not passed through command line, try reading it from wkFile
    open(unit=1,file=wkFile)
    read(1,*) workspace
    close(1)
else
    workspace=trim(wk)
endif

call Config_Model_Read(workspace=trim(workspace),file=modelFile,model=model,err=err,mess=mess)
if(err/=0) call fatalExit(mess)

call loadFormulas(model,err,mess)
if(err/=0) call fatalExit(mess)

call Config_MCMC_Read(workspace=trim(workspace),file=mcmcFile,mcmc=mcmc,err=err,mess=mess)
if(err/=0) call fatalExit(mess)

call OATsampler(model,mcmc,err,mess)
if(err/=0) call fatalExit(mess)

write(*,*) 'Finito!'
!read(*,*)

!*********!
contains
!*********!

subroutine fatalExit(mess)
character(*), intent(in)::mess
write(*,*) '+++++++FATAL+++++++++'
write(*,*) trim(mess)
read(*,*)
STOP
end subroutine fatalExit

subroutine print_help()
    write(*,'(a)') 'usage: STooDs [OPTIONS]'
    write(*,'(a)') 'available options:'
    write(*,'(a)') '  -wk path, --workspace path: set path to workspace'
    write(*,'(a)') '  -sd k, --seed k:            set seed to k (k should be an integer)'
    write(*,'(a)') '  -rd, --random:              randomize seed (=> non-reproducible MCMC runs)'
    write(*,'(a)') '  -v, --version:              print version information and exit'
    write(*,'(a)') '  -h, --help:                 print help and exit'
end subroutine print_help

end program main
