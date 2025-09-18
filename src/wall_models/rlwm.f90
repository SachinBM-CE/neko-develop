! Copyright (c) 2024, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!
!> Implements `rlwm_t`.
module rlwm
  use field, only: field_t
  use num_types, only : rp
  use json_module, only : json_file
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use wall_model, only : wall_model_t
  use field_registry, only : neko_field_registry
  use json_utils, only : json_get_or_default
  use rlwm_cpu, only : rlwm_compute_cpu
  ! use rlwm_device, only : rlwm_compute_device
  use field_math, only: field_invcol3
  use vector, only : vector_t
  use math, only: masked_gather_copy_0
  use device_math, only: device_masked_gather_copy_0
  use scratch_registry, only : neko_scratch_registry
  
  ! ****** TorchFort ******
  use torchfort
  use comm, only : pe_rank, pe_size, NEKO_COMM
  use iso_c_binding
  use math
  use tensor
  use num_types, only : sp
  use mpi_f08, only : MPI_INTEGER, MPI_SUCCESS, MPI_SUM, MPI_Allgather, MPI_Allreduce
  use json_utils, only : json_get
  use utils, only : neko_error
  use operators, only : grad, dudxyz

  implicit none
  private

  !> Wall model based on rlwm's law of the wall.
  !! Reference: http://dx.doi.org/10.1115/1.3641728
  type, public, extends(wall_model_t) :: rlwm_t
     !> The von Karman coefficient.
     real(kind=rp) :: kappa = 0.41_rp
     !> The log-law intercept.
     real(kind=rp) :: B = 5.2_rp
     !> The kinematic viscosity.
     type(vector_t) :: nu
	 
	 !> TorchFort =================================================================================================================
	 !> JSON INPUTS 
	 character(len=256) :: tf_key, yaml_path, log_dir
	 integer :: model_device, rb_device, start_train_tstep, tsteps_rl, n_epochs
	 real(kind=rp) :: tau_true
	 !> Vectors
	 type(vector_t) :: ui_l, vi_l, wi_l, normu_l, magu_l, vg_l, utau_l, tau_old_l, tau_new_l, & 
					   l_star, u_plus, g_plus, h_plus, slope, intercept, &
					   error_new, error_old, rel_error, &
					   reward, total_reward, reward_out, base_reward, bonus_reward, &
					   terminal, terminal_old, terminal_older
	 !> 4D Arrays 
	 real(kind=rp), dimension(:,:,:,:), allocatable :: dudy	
	 !> 2D Arrays
	 real(kind=rp), dimension(:,:), allocatable :: state, action
	 !> MPI
	 integer :: total_agents, episode=0
	 integer, dimension(:), allocatable :: recvcounts, displs, global_recvcounts, global_displs
	 real(kind=rp), dimension(:), allocatable :: global_reward, global_terminal
	 real(kind=rp), dimension(:,:), allocatable :: global_state, global_state_old, global_state_older, & 
												   global_action, global_action_old, global_action_older
	 !> Fields 
	 type(field_t), pointer :: reward_field => null(), slope_field => null(), intercept_field => null()
	 !> Single precision loss values
	 real(kind=sp) :: p_loss_val, q_loss_val
	 !=============================================================================================================================

   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => rlwm_init
     !> Partial constructor from JSON, meant to work as the first stage of
     !! initialization before the `finalize` call.
     procedure, pass(this) :: partial_init => rlwm_partial_init
     !> Finalize the construction using the mask and facet arrays of the bc.
     procedure, pass(this) :: finalize => rlwm_finalize
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => rlwm_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => rlwm_free
     !> Compute the kinematic viscosity at the wall.
     procedure, pass(this) :: compute_nu => rlwm_compute_nu
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute => rlwm_compute
  end type rlwm_t

contains
  !> Constructor from JSON.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param json A dictionary with parameters.
  subroutine rlwm_init(this, scheme_name, coef, msk, facet, h_index, json)
    class(rlwm_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, B
	
	! ****** TorchFort ****** 
	integer :: res, ierr, i
	! ***********************
	
	call json_get_or_default(json, "kappa", kappa, 0.41_rp)
    call json_get_or_default(json, "B", B, 5.2_rp)
	
	call this%init_from_components(scheme_name, coef, msk, facet, h_index, kappa, B)

  end subroutine rlwm_init

  !> Constructor from JSON.
  !! @param coef SEM coefficients.
  !! @param json A dictionary with parameters.
  subroutine rlwm_partial_init(this, coef, json)
    class(rlwm_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json
	
	! ****** TorchFort ****** 
	integer :: res, ierr, i
	real(kind=rp) :: tmp_real
	character(len=:), allocatable :: tmp_string
	! ***********************

    call this%partial_init_base(coef, json)
    call json_get_or_default(json, "kappa", this%kappa, 0.41_rp)
    call json_get_or_default(json, "B", this%B, 5.2_rp)

    print *, "rlwm_partial_init called"

	call json_get(json, "tf_key", tmp_string)
	this%tf_key = trim(tmp_string)

	call json_get(json, "yaml_path", tmp_string)
	this%yaml_path = trim(tmp_string)

	call json_get(json, "log_dir", tmp_string)
	this%log_dir = trim(tmp_string)

	call json_get_or_default(json, "model_device", tmp_real, -1.0_rp)
	this%model_device = int(tmp_real)

	call json_get_or_default(json, "rb_device", tmp_real, -1.0_rp)
	this%rb_device = int(tmp_real)

	call json_get_or_default(json, "tau_true", tmp_real, 0.002_rp)
	this%tau_true = tmp_real
	
	call json_get_or_default(json, "start_train_tstep", tmp_real, 1000.0_rp)
	this%start_train_tstep = int(tmp_real)
	
	call json_get_or_default(json, "tsteps_rl", tmp_real, 100.0_rp)
	this%tsteps_rl = int(tmp_real)
	
    call json_get_or_default(json, "n_epochs", tmp_real, 10.0_rp)
	this%n_epochs = int(tmp_real)
	
	res = torchfort_set_manual_seed(123)
	if (res /= TORCHFORT_RESULT_SUCCESS) stop
	print *, "Result of set_manual_seed : ", res
	print *
	
	res = torchfort_rl_off_policy_create_system(this%tf_key, this%yaml_path, this%model_device, this%rb_device)
	if (res /= TORCHFORT_RESULT_SUCCESS) stop
	print *, "Result of create_system : ", res
	print *
	
	! res = torchfort_rl_off_policy_create_distributed_system(this%tf_key, & 
	! this%yaml_path, NEKO_COMM, this%model_device, this%rb_device)
	! if (res /= TORCHFORT_RESULT_SUCCESS) stop
	! print *, "Result of create_distributed_system : ", res
	
	! print *, "===> pe_size :", pe_size	
	! print *, "this%n_nodes : ", this%n_nodes, "from pe_rank: ", pe_rank
	! print *, "this%msk(0) : ", this%msk(0), "from pe_rank: ", pe_rank
	
	allocate(this%dudy(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))

  end subroutine rlwm_partial_init

  !> Finalize the construction using the mask and facet arrays of the bc.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  subroutine rlwm_finalize(this, msk, facet)
    class(rlwm_t), intent(inout) :: this
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
	
	! ****** TorchFort ****** 
	integer :: res, ierr, i
	! ***********************

    call this%finalize_base(msk, facet)
    call this%nu%init(this%n_nodes)
	
	!> Vectors
	call this%ui_l%init(this%n_nodes)
	call this%vi_l%init(this%n_nodes)
	call this%wi_l%init(this%n_nodes)
	call this%normu_l%init(this%n_nodes)
	call this%magu_l%init(this%n_nodes)
	call this%vg_l%init(this%n_nodes)
	call this%utau_l%init(this%n_nodes)
	call this%tau_old_l%init(this%n_nodes)
	call this%tau_new_l%init(this%n_nodes)
	
	call this%l_star%init(this%n_nodes)
	call this%u_plus%init(this%n_nodes)
	call this%g_plus%init(this%n_nodes)
	call this%h_plus%init(this%n_nodes)
	call this%slope%init(this%n_nodes)
	call this%intercept%init(this%n_nodes)
	
    if (this%n_nodes > 0) then
        allocate(this%state(2, this%n_nodes), this%action(1, this%n_nodes))
    else
        allocate(this%state(2, 1), this%action(1, 1))
    end if
	
	call this%error_new%init(this%n_nodes)
	call this%error_old%init(this%n_nodes)
	call this%rel_error%init(this%n_nodes)
	
	call this%reward%init(this%n_nodes)
	call this%total_reward%init(this%n_nodes)
	call this%reward_out%init(this%n_nodes)
	call this%base_reward%init(this%n_nodes)
	call this%bonus_reward%init(this%n_nodes)
	
	call this%terminal%init(this%n_nodes)
	call this%terminal_old%init(this%n_nodes)
	call this%terminal_older%init(this%n_nodes)
	
	!> Collecting receive counts, displacements & total agents ====================================================================
	
	! allocate(this%recvcounts(pe_size), this%displs(pe_size))
	! allocate(this%global_recvcounts(pe_size), this%global_displs(pe_size))
	allocate(this%recvcounts(0:(pe_size-1)), this%displs(0:(pe_size-1)))
	allocate(this%global_recvcounts(0:(pe_size-1)), this%global_displs(0:(pe_size-1)))
	
	! Gather the number of agents (this%n_nodes) from all ranks onto all ranks
    call MPI_Allgather(this%n_nodes, 1, MPI_INTEGER, this%recvcounts, 1, MPI_INTEGER, NEKO_COMM, ierr)
    if (ierr /= MPI_SUCCESS) then
        call neko_error("MPI_Allgather failed in rlwm_finalize")
    end if
    print *, "recvcounts = ", this%recvcounts
	
    ! Calculate the displacements for MPI_Gatherv (all ranks need this)
    ! this%displs(1) = 0
    ! do i = 2, pe_size
        ! this%displs(i) = this%displs(i-1) + this%recvcounts(i-1)
        ! if (pe_rank == 0) then
            ! print *, "i = ", i, "displs(i) = ", this%displs(i), "displs(i-1) = ", this%displs(i-1), &
                     ! "recvcounts(i-1)", this%recvcounts(i-1)
        ! end if
    ! end do
	! do i = 1, pe_size
        ! this%global_recvcounts(i) = this%recvcounts(i) * 2
        ! this%global_displs(i) = this%displs(i) * 2
    ! end do
    this%displs(0) = 0
    do i = 1, (pe_size - 1)
       this%displs(i) = this%displs(i-1) + this%recvcounts(i-1)
       if (pe_rank == 0) then
           print *, "i = ", i, "    displs(i) = ", this%displs(i), "    displs(i-1) = ", this%displs(i-1), &
                     "    recvcounts(i-1)", this%recvcounts(i-1)
       end if
    end do
    do i = 0, pe_size-1
        this%global_recvcounts(i) = this%recvcounts(i) * 2
        this%global_displs(i) = this%displs(i) * 2
    end do
	
	! Get the total number of agents across all ranks
    call MPI_Allreduce(this%n_nodes, this%total_agents, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    if (ierr /= MPI_SUCCESS) then
        call neko_error("MPI_Allreduce failed in rlwm_finalize")
    end if
    print *, ">>>> total_agents = ", this%total_agents
	
	!==============================================================================================================================
	
    ! Allocate global arrays
    if (pe_rank == 0) then
        allocate(this%global_state(2,this%total_agents))
		allocate(this%global_state_old(2,this%total_agents))
		allocate(this%global_state_older(2,this%total_agents))
        allocate(this%global_action(1,this%total_agents))
		allocate(this%global_action_old(1,this%total_agents))
		allocate(this%global_action_older(1,this%total_agents))
		allocate(this%global_reward(this%total_agents))
		allocate(this%global_terminal(this%total_agents))
	else
		! For non-root processes, these can be unallocated or size 1
		allocate(this%global_state(1,1))
		allocate(this%global_state_old(1,1))
		allocate(this%global_state_older(1,1))
		allocate(this%global_action(1,1))
		allocate(this%global_action_old(1,1))
		allocate(this%global_action_older(1,1))
		allocate(this%global_reward(1))
		allocate(this%global_terminal(1))
    end if
	
    call neko_field_registry%add_field(this%dof, "reward", ignore_existing = .true.)
    this%reward_field => neko_field_registry%get_field("reward")
	
    call neko_field_registry%add_field(this%dof, "slope", ignore_existing = .true.)
    this%slope_field => neko_field_registry%get_field("slope")
	
    call neko_field_registry%add_field(this%dof, "intercept", ignore_existing = .true.)
    this%intercept_field => neko_field_registry%get_field("intercept")
	
	print *, "**********************************"
	print *, "|||| Initialization Completed ||||"
	print *, "**********************************"
	print *	
	
  end subroutine rlwm_finalize

  !> Constructor from components.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param kappa The von Karman coefficient.
  !! @param B The log-law intercept.
  subroutine rlwm_init_from_components(this, scheme_name, coef, msk, &
       facet, h_index, kappa, B)
    class(rlwm_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    real(kind=rp), intent(in) :: kappa
    real(kind=rp), intent(in) :: B

    call this%free()
    call this%init_base(scheme_name, coef, msk, facet, h_index)

    this%kappa = kappa
    this%B = B

    call this%nu%init(this%n_nodes)
  end subroutine rlwm_init_from_components

  !> Compute the kinematic viscosity vector.
  subroutine rlwm_compute_nu(this)
    class(rlwm_t), intent(inout) :: this
    type(field_t), pointer :: temp
    integer :: idx

    call neko_scratch_registry%request_field(temp, idx)
    call field_invcol3(temp, this%mu, this%rho)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_masked_gather_copy_0(this%nu%x_d, temp%x_d, this%msk_d, &
            temp%size(), this%nu%size())
    else
       call masked_gather_copy_0(this%nu%x, temp%x, this%msk, temp%size(), &
            this%nu%size())
    end if

    call neko_scratch_registry%relinquish_field(idx)
  end subroutine rlwm_compute_nu

  !> Destructor for the rlwm_t (base) class.
  subroutine rlwm_free(this)
    class(rlwm_t), intent(inout) :: this

    call this%free_base()
	
	if (allocated(this%dudy)) deallocate(this%dudy)
	if (allocated(this%state)) deallocate(this%state)
	if (allocated(this%action)) deallocate(this%action)	
	
	!> Vectors
	call this%ui_l%free()
	call this%vi_l%free()
	call this%wi_l%free()
	call this%normu_l%free()
	call this%magu_l%free()
	call this%vg_l%free()
	call this%utau_l%free()
	call this%tau_old_l%free()
	call this%tau_new_l%free()
	
	call this%l_star%free()
	call this%u_plus%free()
	call this%g_plus%free()
	call this%h_plus%free()
	call this%slope%free()
	call this%intercept%free()
	
	call this%error_new%free()
	call this%error_old%free()
	call this%rel_error%free()
	
	call this%reward%free()
	call this%total_reward%free()
	call this%reward_out%free()
	call this%base_reward%free()
	call this%bonus_reward%free()
	
	if (allocated(this%recvcounts)) deallocate(this%recvcounts)
	if (allocated(this%displs)) deallocate(this%displs)
	if (allocated(this%global_recvcounts)) deallocate(this%global_recvcounts)
	if (allocated(this%global_displs)) deallocate(this%global_displs)	
	
	if (allocated(this%global_state)) deallocate(this%global_state)
	if (allocated(this%global_state_old)) deallocate(this%global_state_old)
	if (allocated(this%global_state_older)) deallocate(this%global_state_older)
	if (allocated(this%global_action)) deallocate(this%global_action)
	if (allocated(this%global_action_old)) deallocate(this%global_action_old)
	if (allocated(this%global_action_older)) deallocate(this%global_action_older)
	if (allocated(this%global_reward)) deallocate(this%global_reward)
	if (allocated(this%global_terminal)) deallocate(this%global_terminal)
	
	nullify(this%reward_field)

  end subroutine rlwm_free

  !> Compute the wall shear stress.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine rlwm_compute(this, t, tstep)
    class(rlwm_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u
    type(field_t), pointer :: v
    type(field_t), pointer :: w
    integer :: i
    real(kind=rp) :: ui, vi, wi, magu, utau, normu, guess

    call this%compute_nu()

    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")
	
	! Gradient Tensor
	call dudxyz(this%dudy, u%x, this%coef%drdy, this%coef%dsdy, this%coef%dtdy, this%coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call rlwm_compute_device(u%x_d, v%x_d, w%x_d, this%ind_r_d, &
            this%ind_s_d, this%ind_t_d, this%ind_e_d, &
            this%n_x%x_d, this%n_y%x_d, this%n_z%x_d, &
            this%nu%x_d, this%h%x_d, &
            this%tau_x%x_d, this%tau_y%x_d, this%tau_z%x_d, &
            this%n_nodes, u%Xh%lx, &
            this%kappa, this%B, tstep)
    else
       call rlwm_compute_cpu(u%x, v%x, w%x, &
            this%ind_r, this%ind_s, this%ind_t, this%ind_e, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%nu%x, this%h%x, &
            this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, u%Xh%lx, u%msh%nelv, &
            this%kappa, this%B, tstep, & 
			this%tf_key, this%yaml_path, this%log_dir, &
			this%model_device, this%rb_device, this%start_train_tstep, this%tsteps_rl, this%n_epochs, this%tau_true, &
			this%ui_l%x, this%vi_l%x, this%wi_l%x, this%normu_l%x, this%magu_l%x, this%vg_l%x, this%utau_l%x, &
			this%tau_old_l%x, this%tau_new_l%x, &
			this%l_star%x, this%u_plus%x, this%g_plus%x, this%h_plus%x, this%slope%x, this%intercept%x, this%dudy, &
			this%error_new%x, this%error_old%x, this%rel_error%x, &
			this%reward%x, this%total_reward%x, this%reward_out%x, this%base_reward%x, this%bonus_reward%x, &
			this%terminal%x, this%terminal_old%x, this%terminal_older%x, &
			this%recvcounts, this%displs, this%total_agents, this%state, this%action, this%global_state, this%global_action, &
			this%episode, this%global_state_older, this%global_action_older, this%global_reward, this%global_terminal, &
			this%p_loss_val, this%q_loss_val, &
			this%msk, this%reward_field, this%slope_field, this%intercept_field, &
			this%global_recvcounts, this%global_displs)
    end if

  end subroutine rlwm_compute
end module rlwm
