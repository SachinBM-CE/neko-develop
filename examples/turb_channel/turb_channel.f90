! Martin Karp 13/3-2023
! updated initial condition Philipp Schlatter 09/07/2024
module user

  use neko
  
  !> TorchFort =======================================================
  use bc, only: bc_t
  use wall_model_bc, only: wall_model_bc_t
  use rlwm, only: rlwm_t 
  use math
  use comm, only : pe_rank, pe_size, NEKO_COMM
  use torchfort
  use mpi_f08, only : MPI_Gatherv, MPI_DOUBLE_PRECISION, MPI_IN_PLACE
  !===================================================================
  
  implicit none

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%mesh_setup => user_mesh_scale
    user%compute => usercheck  ! TorchFort
  end subroutine user_setup

!==================================================================================================================================
!> To access previous state and action
!==================================================================================================================================  
  subroutine usercheck(time)
    type(time_state_t), intent(in) :: time
    
    ! Extract time information from time object
    real(kind=rp) :: t
    integer :: tstep
    
    ! Local variables
    integer :: i, j, k, e, n, i_wm, res, ierr
    class(bc_t), pointer :: bc
    type(wall_model_bc_t), pointer :: wall_bc
    real(kind=rp):: reward_sum
	
	t = time%t          ! Current time
    tstep = time%tstep  ! Current timestep
    
	do k = 1, neko_user_access%case%fluid%bcs_vel%size()
	
		bc => neko_user_access%case%fluid%bcs_vel%get(k)
    
        select type (bc)
    
            type is (wall_model_bc_t)
            wall_bc => bc
    
            select type(this => wall_bc%wall_model)
       
                type is (rlwm_t)
            
                if (allocated(this%state)) then
				
					! Gather current state from all ranks
					call MPI_Gatherv(this%state, 2*this%n_nodes, MPI_DOUBLE_PRECISION, &
									 this%global_state, 2*this%recvcounts, 2*this%displs, MPI_DOUBLE_PRECISION, &
									 0, NEKO_COMM, ierr)
					! if (pe_rank == 0) then
						! print *, "shape(this%global_state)", shape(this%global_state)
					! end if
							
					! Gather current action from all ranks  
					call MPI_Gatherv(this%action, this%n_nodes, MPI_DOUBLE_PRECISION, &
									 this%global_action, this%recvcounts, this%displs, MPI_DOUBLE_PRECISION, &
									 0, NEKO_COMM, ierr)
					! if (pe_rank == 0) then
						! print *, "shape(this%global_action)", shape(this%global_action)
					! end if
							  
					! Gather current reward from all ranks  
					call MPI_Gatherv(this%reward%x, this%n_nodes, MPI_DOUBLE_PRECISION, &
									 this%global_reward, this%recvcounts, this%displs, MPI_DOUBLE_PRECISION, &
									 0, NEKO_COMM, ierr)
					! if (pe_rank == 0) then
						! print *, "shape(this%global_reward)", shape(this%global_reward)
					! end if
							   
					! Gather current terminal from all ranks  		
					call MPI_Gatherv(this%terminal%x, this%n_nodes, MPI_DOUBLE_PRECISION, &
									 this%global_terminal, this%recvcounts, this%displs, MPI_DOUBLE_PRECISION, &
									 0, NEKO_COMM, ierr)
					! if (pe_rank == 0) then
						! print *, "shape(this%global_terminal)", shape(this%global_terminal)
					! end if

					! Only rank 0 processes global arrays
					if (pe_rank == 0) then
				  
						! Copy global arrays to their "old" and "older" versions
						if (allocated(this%global_state_older)) then
							call copy(this%global_state_older, this%global_state_old, size(this%global_state_old))
						end if
				  
						if (allocated(this%global_state_old)) then
							call copy(this%global_state_old, this%global_state, size(this%global_state))
						end if
				  
						if (allocated(this%global_action_older)) then
							call copy(this%global_action_older, this%global_action_old, size(this%global_action_old))
						end if
				  
						if (allocated(this%global_action_old)) then
							call copy(this%global_action_old, this%global_action, size(this%global_action))
						end if
				  
						! ! Print some global state information for debugging
						! do i = 1, min(5, this%total_agents)  ! Print first 5 global agents
						   ! if (i == 1) then
							  ! write(*, '(A6, A20, A20, A20, A20, A20, A20)') &
							  ! 'i', 'global_state', 'global_state_old', 'global_state_older', &
							  ! 'global_action', 'global_action_old', 'global_action_older' 
						   ! end if
						 
						   ! if (allocated(this%global_state_old) .and. allocated(this%global_state_older) .and. &
							 ! allocated(this%global_action_old) .and. allocated(this%global_action_older)) then
							 ! write(*, '(I6, ES20.5, ES20.5, ES20.5, ES20.5, ES20.5, ES20.5)') &
							 ! i, this%global_state(1,i), this%global_state_old(1,i), this%global_state_older(1,i), &
							 ! this%global_action(1,i), this%global_action_old(1,i), this%global_action_older(1,i)
						   ! else
							 ! write(*, '(I6, ES20.5, A19, A19, ES20.5, A19, A19)') &
							 ! i, this%global_state(1,i), 'not_allocated', 'not_allocated', &
							 ! this%global_action(1,i), 'not_allocated', 'not_allocated'
						   ! end if
						! end do
				  
					end if
					
					if ((mod(tstep, this%tsteps_rl) == 0) .or. (abs(t - time%end_time) .le. 1.0e-3_rp)) then
						reward_sum = glsum(this%total_reward%x(:), this%n_nodes)
						this%episode = this%episode + 1
						! ::: WANDB LOGGING :::
						res = torchfort_rl_off_policy_wandb_log(this%tf_key, "reward_sum", this%episode, reward_sum)
						! res = torchfort_rl_off_policy_wandb_log(this%tf_key, "actor_loss", this%episode, this%p_loss_val)
						! res = torchfort_rl_off_policy_wandb_log(this%tf_key, "critic_loss", this%episode, this%q_loss_val)
						! res = torchfort_rl_off_policy_wandb_log(this%tf_key, "rew_out", this%episode, this%reward_out(1))
						! res = torchfort_rl_off_policy_wandb_log(this%tf_key, "tau_new", this%episode, this%tau_new_l(1))
					end if
			  
					do i = 1, this%n_nodes

						if ((mod(tstep, this%tsteps_rl) == 0) .or. (abs(t - time%end_time) .le. 1.0e-3_rp)) then
							this%terminal%x(i) = 1.0_rp
							this%total_reward%x(i) = 0.0_rp
						else
							this%terminal%x(i) = 0.0_rp
						end if

					end do
			  
					call copy(this%terminal_older%x, this%terminal_old%x, size(this%terminal_old%x))
					call copy(this%terminal_old%x, this%terminal%x, size(this%terminal%x))
          
				end if
			end select
        end select
    end do
  end subroutine usercheck

  ! Rescale mesh, we create a mesh with some refinement close to the wall.
  ! initial mesh: 0..4, -1..1, 0..1.5
  ! mesh size (4*pi,2*delta,4/3*pi)
  ! New mesh can easily be genreated with genmeshbox
  ! OBS refinement is not smooth and the constant values are a bit ad hoc.
  ! Stats converge close to reference DNS
  subroutine user_mesh_scale(msh, time)
    type(mesh_t), intent(inout) :: msh
    type(time_state_t), intent(in) :: time
    integer :: i, p, nvert

    real(kind=rp) :: d, y, viscous_layer, visc_el_h, el_h
    real(kind=rp) :: center_el_h, dist_from_wall
    integer :: el_in_visc_lay, el_in_y
    real(kind=rp) :: llx, llz

    ! target mesh size
    llx = 4.*pi
    llz = 4./3.*pi

    ! rescale mesh
    el_in_y = 18
    el_in_visc_lay = 2
    viscous_layer = 0.0888889
    el_h = 2.0_rp/el_in_y
    visc_el_h = viscous_layer/el_in_visc_lay
    center_el_h = (1.0_rp-viscous_layer)/(el_in_y/2-el_in_visc_lay)

    nvert = size(msh%points)
    do i = 1, nvert
       msh%points(i)%x(1) = llx/4.*msh%points(i)%x(1)
       y = msh%points(i)%x(2)
       if ((1-abs(y)) .le. (el_in_visc_lay*el_h)) then
          dist_from_wall = (1-abs(y))/el_h*visc_el_h
       else
          dist_from_wall = viscous_layer + (1-abs(y)- &
               el_in_visc_lay*el_h)/el_h*center_el_h
       end if
       if (y .gt. 0) msh%points(i)%x(2) = 1.0_rp - dist_from_wall
       if (y .lt. 0) msh%points(i)%x(2) = -1.0_rp + dist_from_wall
       msh%points(i)%x(3) = 2./3.*llz*msh%points(i)%x(3)
    end do

  end subroutine user_mesh_scale

  ! User defined initial condition
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    real(kind=rp) :: uvw(3), x, y, z
    type (field_t), pointer :: u, v, w
    integer :: i

    if (scheme_name .eq. 'fluid') then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       do i = 1, u%size()
          x = u%dof%x(i,1,1,1)
          y = u%dof%y(i,1,1,1)
          z = u%dof%z(i,1,1,1)

          uvw = channel_ic(x, y, z)

          u%x(i,1,1,1) = uvw(1)
          v%x(i,1,1,1) = uvw(2)
          w%x(i,1,1,1) = uvw(3)
       end do
    end if
  end subroutine initial_conditions

  ! Kind of brute force with rather large initial disturbances
  function channel_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: uvw(3)
    real(kind=rp) :: ux, uy, uz, eps, Re_tau, yp, Re_b, alpha, beta
    real(kind=rp) :: C, k, kx, kz, eps1, ran

    real(kind=rp) :: llx, llz

    llx = 4.*pi
    llz = 4./3.*pi

    Re_tau = 180
    C = 5.17
    k = 0.41
    Re_b = 2800

    yp = (1-y)*Re_tau
    if (y .lt. 0) yp = (1+y)*Re_tau

    ! Reichardt function
    ux = 1/k*log(1.0+k*yp) + (C - (1.0/k)*log(k)) * &
         (1.0 - exp(-yp/11.0) - yp/11*exp(-yp/3.0))
    ux = ux * Re_tau/Re_b

    ! actually, sometimes one may not use the turbulent profile, but
    ! rather the parabolic lamianr one
    ! ux = 1.5*(1-y**2)

    ! add perturbations to trigger turbulence
    ! base flow
    uvw(1) = ux
    uvw(2) = 0
    uvw(3) = 0

    ! first, large scale perturbation
    eps = 0.05
    kx = 3
    kz = 4
    alpha = kx * 2*PI/llx
    beta = kz * 2*PI/llz
    uvw(1) = uvw(1) + eps*beta * sin(alpha*x)*cos(beta*z)
    uvw(2) = uvw(2) + eps * sin(alpha*x)*sin(beta*z)
    uvw(3) = uvw(3) -eps*alpha * cos(alpha*x)*sin(beta*z)

    ! second, small scale perturbation
    eps = 0.005
    kx = 17
    kz = 13
    alpha = kx * 2*PI/llx
    beta = kz * 2*PI/llz
    uvw(1) = uvw(1) + eps*beta * sin(alpha*x)*cos(beta*z)
    uvw(2) = uvw(2) + eps * sin(alpha*x)*sin(beta*z)
    uvw(3) = uvw(3) -eps*alpha * cos(alpha*x)*sin(beta*z)

    ! finally, random perturbations only in y
    eps1 = 0.001
    ran = sin(-20*x*z+y**3*tan(x*z**2)+100*z*y-20*sin(x*y*z)**5)
    uvw(2) = uvw(2) + eps1*ran

  end function channel_ic

end module user
