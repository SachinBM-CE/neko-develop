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
  
  ! TorchFort
  use torchfort
  use operators, only : grad, dudxyz
  use comm, only : pe_rank, pe_size, NEKO_COMM
  use iso_c_binding
  use math
  use tensor
  use num_types, only : sp
  use mpi_f08
  use json_utils, only : json_get
  use utils, only : neko_error

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
	 
	 ! TorchFort
	 integer :: model_device, rb_device, total_agents
	 integer :: episode=0, start_train_tstep, tsteps_rl, n_epochs
	 character(len=256) :: tf_key
	 character(len=256) :: yaml_path 
	 character(len=256) :: log_dir 
	 real(kind=sp) :: p_loss_val, q_loss_val
	 real(kind=rp) :: tau_true
	 real(kind=rp), dimension(:,:,:,:), allocatable :: dudy
	 real(kind=rp), dimension(:), allocatable :: ui_l, vi_l, wi_l, normu_l, magu_l, h_l
	 real(kind=rp), dimension(:), allocatable :: vg_l, utau_l, tau_old_l, tau_new_l
	 real(kind=rp), dimension(:), allocatable :: l_star, u_plus, g_plus, h_plus
	 real(kind=rp), dimension(:,:), allocatable :: state, action, state_transposed, action_transposed
	 real(kind=rp), dimension(:,:), allocatable :: state_old, state_older, action_old, action_older
	 real(kind=rp), dimension(:), allocatable :: terminal, terminal_old, terminal_older, reward, total_reward, reward_out
	 real(kind=rp), dimension(:), allocatable :: error_new_l, error_old_l, base_reward_l, rel_error_l, bonus_l
	 integer, dimension(:), allocatable :: recvcounts, displs
	 real(kind=rp), dimension(:,:), allocatable :: global_state, global_action
	 real(kind=rp), dimension(:,:), allocatable :: global_state_old, global_state_older, global_action_old, global_action_older
	 real(kind=rp), dimension(:), allocatable :: global_reward, global_terminal
	 real(kind=rp), dimension(3) :: running_mean, running_var
     real(kind=rp) :: running_count
	 !> Pointers to the state and reward fields in the registry
     type(field_t), pointer :: reward_field => null()
	 
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => rlwm_init
     !> Partial constructor from JSON, meant to work as the first stage of
     !! initialization before the `finalize` call.
     procedure, pass(this) :: partial_init => rlwm_partial_init
     !> Finalize the construction using the mask and facet arrays of the bc.
     procedure, pass(this) :: finalize => rlwm_finalize
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          rlwm_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => rlwm_free
     !> Compute the kinematic viscosity at the wall.
     procedure, pass(this) :: compute_nu => rlwm_compute_nu
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute => rlwm_compute
	 !> TorchFort
	 ! procedure, private, pass(this) :: calculate_reward, normalize_state, print_outputs, print_rewards, print_states, print_debug_info
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
	real(kind=rp) :: tmp_real
	character(len=:), allocatable :: tmp_string
	! ***********************

    print *, "RLWM_init called"
	
	call json_get_or_default(json, "kappa", kappa, 0.41_rp)
    call json_get_or_default(json, "B", B, 5.2_rp)

	call json_get(json, "tf_key", tmp_string)
	this%tf_key = tmp_string

	call json_get(json, "yaml_path", tmp_string)
	this%yaml_path = tmp_string

	call json_get(json, "log_dir", tmp_string)
	this%log_dir = tmp_string

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

    call this%init_from_components(scheme_name, coef, msk, facet, h_index, &
         kappa, B)
		 
	res = torchfort_rl_off_policy_create_system(this%tf_key, this%yaml_path, this%model_device, this%rb_device)
	if (res /= TORCHFORT_RESULT_SUCCESS) stop
	print *, "Result of create_system : ", res
	print *
	
	! res = torchfort_rl_off_policy_create_distributed_system(this%tf_key, & 
	! this%yaml_path, NEKO_COMM, this%model_device, this%rb_device)
	! if (res /= TORCHFORT_RESULT_SUCCESS) stop
	! print *, "Result of create_distributed_system : ", res
	
	allocate(this%dudy(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
	
	allocate(this%ui_l(this%n_nodes), this%vi_l(this%n_nodes), this%wi_l(this%n_nodes), this%h_l(this%n_nodes))
	allocate(this%normu_l(this%n_nodes), this%magu_l(this%n_nodes), this%vg_l(this%n_nodes), this%utau_l(this%n_nodes))
	allocate(this%tau_old_l(this%n_nodes), this%tau_new_l(this%n_nodes))
	allocate(this%l_star(this%n_nodes), this%u_plus(this%n_nodes), this%g_plus(this%n_nodes), this%h_plus(this%n_nodes))
	allocate(this%state(3, this%n_nodes), this%action(1, this%n_nodes))
	allocate(this%state_transposed(3, this%n_nodes), this%action_transposed(1, this%n_nodes))
	
	allocate(this%state_old(3, this%n_nodes), this%state_older(3, this%n_nodes))
	allocate(this%action_old(1, this%n_nodes), this%action_older(1, this%n_nodes))
	allocate(this%terminal(this%n_nodes), this%terminal_old(this%n_nodes), this%terminal_older(this%n_nodes))
	allocate(this%reward(this%n_nodes), this%total_reward(this%n_nodes), this%reward_out(this%n_nodes))
	
	allocate(this%error_new_l(this%n_nodes), this%error_old_l(this%n_nodes))
	allocate(this%base_reward_l(this%n_nodes), this%rel_error_l(this%n_nodes), this%bonus_l(this%n_nodes))
	
	allocate(this%recvcounts(pe_size), this%displs(pe_size))
	
	print *, "===> pe_size :", pe_size	
	print *, "this%n_nodes : ", this%n_nodes, "from pe_rank: ", pe_rank
	
	! Gather the number of agents (this%n_nodes) from all ranks onto all ranks
    call mpi_allgather(this%n_nodes, 1, MPI_INTEGER, this%recvcounts, 1, MPI_INTEGER, NEKO_COMM, ierr)
    if (ierr /= MPI_SUCCESS) then
        call neko_error("MPI_ALLGATHER failed in RLWM_init")
    end if
    print *, "recvcounts = ", this%recvcounts
	
    ! Calculate the displacements for MPI_Gatherv (all ranks need this)
    this%displs(1) = 0
    do i = 2, pe_size
        this%displs(i) = this%displs(i-1) + this%recvcounts(i-1)
        if (pe_rank == 0) then
            print *, "i = ", i, "displs(i) = ", this%displs(i), "displs(i-1) = ", this%displs(i-1), &
                     "recvcounts(i-1)", this%recvcounts(i-1)
        end if
    end do
	
	! Get the total number of agents across all ranks
    call mpi_allreduce(this%n_nodes, this%total_agents, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    if (ierr /= MPI_SUCCESS) then
        call neko_error("MPI_ALLREDUCE failed in RLWM_init")
    end if
    print *, ">>>> total_agents = ", this%total_agents
	
    ! Allocate global arrays
    if (pe_rank == 0) then
        allocate(this%global_state(3,this%total_agents))
        allocate(this%global_action(1,this%total_agents))
		allocate(this%global_state_old(3,this%total_agents))
		allocate(this%global_state_older(3,this%total_agents))
		allocate(this%global_action_old(1,this%total_agents))
		allocate(this%global_action_older(1,this%total_agents))
		allocate(this%global_reward(this%total_agents))
		allocate(this%global_terminal(this%total_agents))
    else
        ! Allocate dummy arrays on non-root ranks (or you could use nullify/not allocate)
        allocate(this%global_state(1, 1))
        allocate(this%global_action(1, 1))
    end if
	
	! Initialize running statistics
    this%running_mean = 0.0_rp
    this%running_var = 1.0_rp
    this%running_count = 1e-8_rp
	
    call neko_field_registry%add_field(this%dof, "reward", ignore_existing = .true.)
    this%reward_field => neko_field_registry%get_field("reward")
	
	print *, "**********************************"
	print *, "|||| Initialization Completed ||||"
	print *, "**********************************"
	print *
	
  end subroutine rlwm_init

  !> Constructor from JSON.
  !! @param coef SEM coefficients.
  !! @param json A dictionary with parameters.
  subroutine rlwm_partial_init(this, coef, json)
    class(rlwm_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json

    call this%partial_init_base(coef, json)
    call json_get_or_default(json, "kappa", this%kappa, 0.41_rp)
    call json_get_or_default(json, "B", this%B, 5.2_rp)

  end subroutine rlwm_partial_init

  !> Finalize the construction using the mask and facet arrays of the bc.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  subroutine rlwm_finalize(this, msk, facet)
    class(rlwm_t), intent(inout) :: this
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)

    call this%finalize_base(msk, facet)
    call this%nu%init(this%n_nodes)
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
	
	if (allocated(this%ui_l)) deallocate(this%ui_l)
	if (allocated(this%vi_l)) deallocate(this%vi_l)
	if (allocated(this%wi_l)) deallocate(this%wi_l)
	if (allocated(this%h_l)) deallocate(this%h_l)
	if (allocated(this%normu_l)) deallocate(this%normu_l)
	if (allocated(this%magu_l)) deallocate(this%magu_l)
	
	if (allocated(this%vg_l)) deallocate(this%vg_l)
	if (allocated(this%utau_l)) deallocate(this%utau_l)
	if (allocated(this%tau_old_l)) deallocate(this%tau_old_l)
	if (allocated(this%tau_new_l)) deallocate(this%tau_new_l)
	
    if (allocated(this%l_star)) deallocate(this%l_star)
    if (allocated(this%u_plus)) deallocate(this%u_plus)
    if (allocated(this%g_plus)) deallocate(this%g_plus)
    if (allocated(this%h_plus)) deallocate(this%h_plus)
	
	if (allocated(this%state)) deallocate(this%state)
	if (allocated(this%action)) deallocate(this%action)
    if (allocated(this%state_transposed)) deallocate(this%state_transposed)
    if (allocated(this%action_transposed)) deallocate(this%action_transposed)
	
	if (allocated(this%state_old)) deallocate(this%state_old)
	if (allocated(this%state_older)) deallocate(this%state_older)
	if (allocated(this%action_old)) deallocate(this%action_old)
	if (allocated(this%action_older)) deallocate(this%action_older)
	
	if (allocated(this%terminal)) deallocate(this%terminal)
	if (allocated(this%terminal_old)) deallocate(this%terminal_old)
	if (allocated(this%terminal_older)) deallocate(this%terminal_older)	
	if (allocated(this%reward)) deallocate(this%reward)
    if (allocated(this%total_reward)) deallocate(this%total_reward)
	if (allocated(this%reward_out)) deallocate(this%reward_out)
	
	if (allocated(this%error_new_l)) deallocate(this%error_new_l)
    if (allocated(this%error_old_l)) deallocate(this%error_old_l)
    if (allocated(this%base_reward_l)) deallocate(this%base_reward_l)
    if (allocated(this%rel_error_l)) deallocate(this%rel_error_l)
    if (allocated(this%bonus_l)) deallocate(this%bonus_l)
	
	if (allocated(this%recvcounts)) deallocate(this%recvcounts)
	if (allocated(this%displs)) deallocate(this%displs)
	if (allocated(this%global_state)) deallocate(this%global_state)
	if (allocated(this%global_action)) deallocate(this%global_action)
	if (allocated(this%global_state_old)) deallocate(this%global_state_old)
	if (allocated(this%global_action_old)) deallocate(this%global_action_old)
	if (allocated(this%global_state_older)) deallocate(this%global_state_older)
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
            this%kappa, this%B, tstep, t, this%coef, &
			this%start_train_tstep, this%tsteps_rl, this%n_epochs, &
			this%tau_true, this%tf_key, &
			this%dudy, this%ui_l, this%vi_l, this%wi_l, this%normu_l, &
			this%magu_l, this%h_l, this%vg_l, this%utau_l, &
			this%tau_old_l, this%tau_new_l, this%l_star, &
			this%u_plus, this%g_plus, this%h_plus, &
			this%state, this%action, this%reward, &
			this%error_new_l, this%error_old_l, this%base_reward_l, &
			this%rel_error_l, this%bonus_l, this%total_reward, &
			this%reward_out, this%running_mean, this%running_var, &
			this%running_count, this%total_agents, this%reward_field, &
			this%recvcounts, this%displs, this%global_state, &
			this%global_action, this%global_state_older, &
			this%global_action_older, this%global_reward, &
			this%global_terminal, this%p_loss_val, this%q_loss_val)
    end if

  end subroutine rlwm_compute
end module rlwm
