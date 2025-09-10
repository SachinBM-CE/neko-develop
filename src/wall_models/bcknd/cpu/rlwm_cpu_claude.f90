! ! Copyright (c) 2025, The Neko Authors
! ! All rights reserved.
! !
! ! Redistribution and use in source and binary forms, with or without
! ! modification, are permitted provided that the following conditions
! ! are met:
! !
! !   * Redistributions of source code must retain the above copyright
! !     notice, this list of conditions and the following disclaimer.
! !
! !   * Redistributions in binary form must reproduce the above
! !     copyright notice, this list of conditions and the following
! !     disclaimer in the documentation and/or other materials provided
! !     with the distribution.
! !
! !   * Neither the name of the authors nor the names of its
! !     contributors may be used to endorse or promote products derived
! !     from this software without specific prior written permission.
! !
! ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! ! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! ! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! ! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! ! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! ! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! ! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! ! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! ! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! ! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! ! POSSIBILITY OF SUCH DAMAGE.
! !
! !> Implements the CPU kernel for the `rlwm_t` type.
! module rlwm_cpu
  ! use num_types, only : rp
  ! use logger, only : neko_log, NEKO_LOG_DEBUG, LOG_SIZE
  ! implicit none
  ! private

  ! public :: rlwm_compute_cpu

! contains
  ! !> Compute the wall shear stress on cpu using rlwm's model.
  ! !! @param t The time value.
  ! !! @param tstep The current time-step.
  ! subroutine rlwm_compute_cpu(u, v, w, ind_r, ind_s, ind_t, ind_e, &
       ! n_x, n_y, n_z, nu, h, tau_x, tau_y, tau_z, n_nodes, lx, nelv, &
       ! kappa, B, tstep)
    ! integer, intent(in) :: n_nodes, lx, nelv, tstep
    ! real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w
    ! integer, intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    ! real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h, nu
    ! real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    ! real(kind=rp), intent(in) :: kappa, B
    ! integer :: i
    ! real(kind=rp) :: ui, vi, wi, magu, utau, normu, guess

    ! do i=1, n_nodes
       ! ! Sample the velocity
       ! ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       ! vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       ! wi = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))

       ! ! Project on tangential direction
       ! normu = ui * n_x(i) + vi * n_y(i) + wi * n_z(i)

       ! ui = ui - normu * n_x(i)
       ! vi = vi - normu * n_y(i)
       ! wi = wi - normu * n_z(i)

       ! magu = sqrt(ui**2 + vi**2 + wi**2)

       ! ! Get initial guess for Newton solver
       ! if (tstep .eq. 1) then
          ! guess = sqrt(magu * nu(i) / h(i))
       ! else
          ! guess = tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2
          ! guess = sqrt(sqrt(guess))
       ! end if

       ! utau = solve_cpu(magu, h(i), guess, nu(i), kappa, B)

       ! ! Distribute according to the velocity vector
       ! tau_x(i) = -utau**2 * ui / magu
       ! tau_y(i) = -utau**2 * vi / magu
       ! tau_z(i) = -utau**2 * wi / magu
    ! end do

  ! end subroutine rlwm_compute_cpu

  ! !> Newton solver for the algebraic equation defined by the law on cpu.
  ! !! @param u The velocity value.
  ! !! @param y The wall-normal distance.
  ! !! @param guess Initial guess.
  ! !! @param nu The molecular kinematic viscosity.
  ! !! @param kappa The von Karman constant.
  ! !! @param B The log-law intercept.
  ! function solve_cpu(u, y, guess, nu, kappa, B) result(utau)
    ! real(kind=rp), intent(in) :: u
    ! real(kind=rp), intent(in) :: y
    ! real(kind=rp), intent(in) :: guess
    ! real(kind=rp), intent(in) :: nu, kappa, B
    ! real(kind=rp) :: yp, up, utau
    ! real(kind=rp) :: error, f, df, old
    ! integer :: niter, k, maxiter
    ! character(len=LOG_SIZE) :: log_msg

    ! utau = guess

    ! maxiter = 100

    ! do k=1, maxiter
       ! up = u / utau
       ! yp = y * utau / nu
       ! niter = k
       ! old = utau

       ! ! Evaluate function and its derivative
       ! f = (up + exp(-kappa*B)* &
            ! (exp(kappa*up) - 1.0_rp - kappa*up - 0.5_rp*(kappa*up)**2 - &
            ! 1.0_rp/6*(kappa*up)**3) - yp)

       ! df = (-y / nu - u/utau**2 - kappa*up/utau*exp(-kappa*B) * &
            ! (exp(kappa*up) - 1 - kappa*up - 0.5*(kappa*up)**2))

       ! ! Update solution
       ! utau = utau - f / df

       ! error = abs((old - utau)/old)

       ! if (error < 1e-3) then
          ! exit
       ! endif

    ! enddo

    ! if (niter .eq. maxiter) then
       ! write(log_msg, *) "Newton not converged", error, f, utau, old, guess
       ! call neko_log%message(log_msg, NEKO_LOG_DEBUG)
    ! end if
  ! end function solve_cpu
! end module rlwm_cpu

! Copyright (c) 2025, The Neko Authors
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
!> Implements the CPU kernel for the `rlwm_t` type with RL integration.
module rlwm_cpu
  use num_types, only : rp, sp
  use logger, only : neko_log, NEKO_LOG_DEBUG, LOG_SIZE
  use field, only : field_t
  use operators, only : dudxyz
  use comm, only : pe_rank, pe_size, NEKO_COMM
  use math
  use mpi_f08
  use utils, only : neko_error
  
  ! TorchFort
  use torchfort
  use coefs, only : coef_t
  
  implicit none
  private

  public :: rlwm_compute_cpu

contains
  !> Compute the wall shear stress on cpu using rlwm's RL model.
  subroutine rlwm_compute_cpu(u, v, w, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, nu, h, tau_x, tau_y, tau_z, n_nodes, lx, nelv, &
       kappa, B, tstep, t, coef, &
       start_train_tstep, tsteps_rl, n_epochs, tau_true, tf_key, &
       dudy, ui_l, vi_l, wi_l, normu_l, magu_l, h_l, vg_l, utau_l, &
       tau_old_l, tau_new_l, l_star, u_plus, g_plus, h_plus, &
       state, action, reward, error_new_l, error_old_l, base_reward_l, &
       rel_error_l, bonus_l, total_reward, reward_out, &
       running_mean, running_var, running_count, total_agents, &
       reward_field, recvcounts, displs, global_state, global_action, &
       global_state_older, global_action_older, global_reward, &
       global_terminal, p_loss_val, q_loss_val)
    
    integer, intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w
    integer, intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h, nu
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    real(kind=rp), intent(in) :: kappa, B, t, tau_true
    integer, intent(in) :: start_train_tstep, tsteps_rl, n_epochs, total_agents
    character(len=256), intent(in) :: tf_key
    type(coef_t), intent(in) :: coef
    
    ! RL-specific arrays
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: dudy
    real(kind=rp), dimension(n_nodes), intent(inout) :: ui_l, vi_l, wi_l, normu_l, magu_l, h_l
    real(kind=rp), dimension(n_nodes), intent(inout) :: vg_l, utau_l, tau_old_l, tau_new_l
    real(kind=rp), dimension(n_nodes), intent(inout) :: l_star, u_plus, g_plus, h_plus
    real(kind=rp), dimension(3, n_nodes), intent(inout) :: state
    real(kind=rp), dimension(1, n_nodes), intent(inout) :: action
    real(kind=rp), dimension(n_nodes), intent(inout) :: reward, total_reward, reward_out
    real(kind=rp), dimension(n_nodes), intent(inout) :: error_new_l, error_old_l, base_reward_l
    real(kind=rp), dimension(n_nodes), intent(inout) :: rel_error_l, bonus_l
    real(kind=rp), dimension(3), intent(inout) :: running_mean, running_var
    real(kind=rp), intent(inout) :: running_count
    type(field_t), pointer, intent(in) :: reward_field
    
    ! MPI arrays
    integer, dimension(pe_size), intent(in) :: recvcounts, displs
    real(kind=rp), dimension(:,:), intent(inout) :: global_state, global_action
    real(kind=rp), dimension(:,:), intent(inout) :: global_state_older, global_action_older
    real(kind=rp), dimension(:), intent(inout) :: global_reward, global_terminal
    real(kind=sp), intent(inout) :: p_loss_val, q_loss_val
    
    ! Local variables
    integer :: i, res, epoch, ierr
    real(kind=rp) :: ui, vi, wi, magu, utau, normu
	logical :: is_ready = .false.

    if (tstep == 1) then
      print *, ""
      print *, "======================================================"
      print *, "  Verifying RLWM parameters from JSON file:"
      print *, "======================================================"
      print '(" 1> tf_key            : ", A)', tf_key
      print '(" 2> tau_true          : ", F0.6)', tau_true  
      print '(" 3> start_train_tstep : ", I0)', start_train_tstep
      print '(" 4> tsteps_rl         : ", I0)', tsteps_rl
      print '(" 5> n_epochs          : ", I0)', n_epochs
      print *, "======================================================"
      print *, ""
    end if

    ! Gradient Tensor
    call dudxyz(dudy, u, coef%drdy, coef%dsdy, coef%dtdy, coef)

    do i=1, n_nodes
      ! Sample the velocity
      ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
      vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
      wi = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))

      ! Project on tangential direction
      normu = ui * n_x(i) + vi * n_y(i) + wi * n_z(i)
      ui = ui - normu * n_x(i)
      vi = vi - normu * n_y(i)
      wi = wi - normu * n_z(i)

      ! Magnitude of Velocity
      magu = sqrt(ui**2 + vi**2 + wi**2)
      
      ! Arrays to use in the do loop for action
      ui_l(i) = ui
      vi_l(i) = vi
      wi_l(i) = wi
      normu_l(i) = normu
      magu_l(i) = magu
      
      !=====================================
      ! Get initial guess for Newton solver
      !=====================================
      if (tstep .eq. 1) then
        ! Guess velocity
        vg_l(i) = sqrt(magu_l(i) * nu(i) / h(i))
        ! Friction velocity
        utau_l(i) = solve_cpu(magu_l(i), h(i), vg_l(i), nu(i), kappa, B)
        h_l(i) = h(i)
        tau_old_l(i) = utau_l(i)**2
        ! Distribute according to the velocity vector
        tau_x(i) = -utau_l(i)**2 * ui_l(i) / magu_l(i)
        tau_y(i) = -utau_l(i)**2 * vi_l(i) / magu_l(i)
        tau_z(i) = -utau_l(i)**2 * wi_l(i) / magu_l(i)
          
      !=====================================
      ! Spalding for collecting statistics
      !=====================================
      else if (tstep .le. start_train_tstep) then
        ! Magnitude of Shear Stress
        tau_old_l(i) = sqrt(tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2)
        ! Guess Velocity
        vg_l(i) = sqrt(tau_old_l(i))
        ! Friction velocity
        utau_l(i) = solve_cpu(magu_l(i), h(i), vg_l(i), nu(i), kappa, B)
        h_l(i) = h(i)
        tau_new_l(i) = utau_l(i)**2
        ! Calculation of Baseline Reward (using Spalding)
        call calculate_reward_cpu(i, tau_true, tau_new_l(i), tau_old_l(i), &
                                  error_new_l(i), error_old_l(i), base_reward_l(i), &
                                  rel_error_l(i), bonus_l(i), reward(i), total_reward(i))
        ! Copy reward value to the field
        reward_field%x(ind_r(i), ind_s(i), ind_t(i), ind_e(i)) = reward(i)
        ! Distribute according to the velocity vector
        tau_x(i) = -utau_l(i)**2 * ui_l(i) / magu_l(i)
        tau_y(i) = -utau_l(i)**2 * vi_l(i) / magu_l(i)
        tau_z(i) = -utau_l(i)**2 * wi_l(i) / magu_l(i)
         
      else
        ! Magnitude of Shear Stress
        tau_old_l(i) = sqrt(tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2)  
        utau_l(i) = sqrt(tau_old_l(i))
        vg_l(i) = utau_l(i)
      end if
      
      ! Normalization w.r.t. viscous scales (nu, utau)
      l_star(i) = nu(i) / (utau_l(i) + 1e-6)
      u_plus(i) = magu_l(i) / (utau_l(i) + 1e-6)
      g_plus(i) = dudy(ind_r(i), ind_s(i), ind_t(i), ind_e(i)) &
                  / ( (utau_l(i) + 1e-6) / l_star(i))
      h_plus(i) = h(i) / (l_star(i) + 1e-6)

      ! Changing to normalized states
      state(1,i) = h_plus(i) * g_plus(i) - (log(h_plus(i))/kappa)
      state(2,i) = u_plus(i) - (log(h_plus(i))/kappa)
      state(3,i) = h_plus(i)

    end do

    ! Normalize the state w.r.t. mean & std
    call normalize_state_cpu(state, n_nodes, total_agents, &
                             running_mean, running_var, running_count)

    ! --- MPI Gatherv ----------------------------------------------------------------------------
    call mpi_gatherv(state, 3*n_nodes, MPI_DOUBLE_PRECISION, &
                     global_state, 3*recvcounts, 3*displs, MPI_DOUBLE_PRECISION, &
                     0, NEKO_COMM, ierr)

    ! Get actions
    if (pe_rank .eq. 0) then
      res = torchfort_rl_off_policy_predict_explore(tf_key, global_state, global_action)
      if (res /= TORCHFORT_RESULT_SUCCESS) stop
    end if

    ! --- MPI Scatterv ----------------------------------------------------------------------------
    call mpi_scatterv(global_action, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                     action, n_nodes, MPI_DOUBLE_PRECISION, &
                     0, NEKO_COMM, ierr)

    ! Replay Buffer & Train (When using RL)
    if (tstep .gt. start_train_tstep) then
      
      do i=1, n_nodes
        tau_new_l(i) = tau_old_l(i) * action(1,i)
        utau_l(i) = sqrt(tau_new_l(i))
        ! Distribute according to the velocity vector
        tau_x(i) = -utau_l(i)**2 * ui_l(i) / magu_l(i)
        tau_y(i) = -utau_l(i)**2 * vi_l(i) / magu_l(i)
        tau_z(i) = -utau_l(i)**2 * wi_l(i) / magu_l(i)
        
        call calculate_reward_cpu(i, tau_true, tau_new_l(i), tau_old_l(i), &
                                  error_new_l(i), error_old_l(i), base_reward_l(i), &
                                  rel_error_l(i), bonus_l(i), reward(i), total_reward(i))
      end do

      if (pe_rank .eq. 0) then
        res = torchfort_rl_off_policy_update_replay_buffer(tf_key, global_state_older, &
        global_action_older, global_state, global_reward, global_terminal)
        if (res /= TORCHFORT_RESULT_SUCCESS) stop
        print *, "result of update_replay_buffer [Multi]: ", res
      end if

      do epoch = 1, n_epochs
        if (mod(tstep, tsteps_rl) == 0 .and. tstep > 3) then
          res = torchfort_rl_off_policy_is_ready(tf_key, is_ready)
          if (res /= TORCHFORT_RESULT_SUCCESS) stop
          if (is_ready) then
            res = torchfort_rl_off_policy_train_step(tf_key, p_loss_val, q_loss_val)
            print *, "p_loss_val = ", p_loss_val
            print *, "q_loss_val = ", q_loss_val
            print *, "result of train_step: ", res
          end if
        end if
      end do

      res = torchfort_rl_off_policy_evaluate(tf_key, global_state_older, global_action_older, reward_out)
      res = torchfort_rl_off_policy_save_checkpoint(tf_key, "checkpoint_dir")
    end if

  end subroutine rlwm_compute_cpu

  !> Calculate reward for agent i
  subroutine calculate_reward_cpu(i, tau_true, tau_new, tau_old, &
                                  error_new, error_old, base_reward, &
                                  rel_error, bonus, reward, total_reward)
    integer, intent(in) :: i
    real(kind=rp), intent(in) :: tau_true, tau_new, tau_old
    real(kind=rp), intent(out) :: error_new, error_old, base_reward
    real(kind=rp), intent(out) :: rel_error, bonus, reward
    real(kind=rp), intent(inout) :: total_reward
    
    error_new = abs(tau_true - tau_new)
    error_old = abs(tau_true - tau_old)
    base_reward = (error_new - error_old) / tau_true
    
    rel_error = error_new / tau_true
    if (rel_error < 0.01_rp) then
      bonus = 1.0_rp - rel_error
    else
      bonus = 0.0_rp
    end if
    
    reward = base_reward + bonus
    total_reward = total_reward + reward
  end subroutine calculate_reward_cpu

  !> Normalize state using running statistics
  subroutine normalize_state_cpu(state, n_nodes, total_agents, &
                                 running_mean, running_var, running_count)
    integer, intent(in) :: n_nodes, total_agents
    real(kind=rp), dimension(3, n_nodes), intent(inout) :: state
    real(kind=rp), dimension(3), intent(inout) :: running_mean, running_var
    real(kind=rp), intent(inout) :: running_count
    
    integer :: feat_idx, ierr
    real(kind=rp) :: batch_count, total_count
    real(kind=rp), dimension(3) :: local_batch_mean, local_batch_var
    real(kind=rp), dimension(3) :: global_batch_mean, global_batch_var
    real(kind=rp), dimension(3) :: delta, new_mean, m_a, m_b, m_2, new_var

    batch_count = real(n_nodes, kind=rp)
    if (batch_count == 0.0_rp) return

    do feat_idx = 1, 3
      local_batch_mean(feat_idx) = sum(state(feat_idx, :)) / batch_count
      local_batch_var(feat_idx) = sum((state(feat_idx, :) - local_batch_mean(feat_idx))**2) / batch_count
    end do

    call mpi_allreduce(local_batch_mean * batch_count, global_batch_mean, 3, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    global_batch_mean = global_batch_mean / total_agents

    do feat_idx = 1, 3
      local_batch_var(feat_idx) = sum((state(feat_idx, :) - global_batch_mean(feat_idx))**2)
    end do
    call mpi_allreduce(local_batch_var, global_batch_var, 3, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    global_batch_var = global_batch_var / total_agents

    delta = global_batch_mean - running_mean
    total_count = running_count + total_agents

    new_mean = running_mean + delta * total_agents / total_count
    m_a = running_var * running_count
    m_b = global_batch_var * total_agents
    m_2 = m_a + m_b + (delta**2) * running_count * total_agents / total_count
    new_var = m_2 / total_count

    running_mean = new_mean
    running_var = new_var
    running_count = total_count

    do feat_idx = 1, 3
      call cadd(state(feat_idx, :), -running_mean(feat_idx), n_nodes)
      call cmult(state(feat_idx, :), 1.0_rp / (sqrt(running_var(feat_idx)) + 1e-8), n_nodes)
    end do

  end subroutine normalize_state_cpu

  !> Newton solver for the algebraic equation defined by the law on cpu.
  function solve_cpu(u, y, guess, nu, kappa, B) result(utau)
    real(kind=rp), intent(in) :: u, y, guess, nu, kappa, B
    real(kind=rp) :: yp, up, utau, error, f, df, old
    integer :: niter, k, maxiter
    character(len=LOG_SIZE) :: log_msg

    utau = guess
    maxiter = 100

    do k=1, maxiter
      up = u / utau
      yp = y * utau / nu
      niter = k
      old = utau

      f = (up + exp(-kappa*B)* &
           (exp(kappa*up) - 1.0_rp - kappa*up - 0.5_rp*(kappa*up)**2 - &
           1.0_rp/6*(kappa*up)**3) - yp)

      df = (-y / nu - u/utau**2 - kappa*up/utau*exp(-kappa*B) * &
           (exp(kappa*up) - 1 - kappa*up - 0.5*(kappa*up)**2))

      utau = utau - f / df
      error = abs((old - utau)/old)

      if (error < 1e-3) then
        exit
      endif
    enddo

    if (niter .eq. maxiter) then
      write(log_msg, *) "Newton not converged", error, f, utau, old, guess
      call neko_log%message(log_msg, NEKO_LOG_DEBUG)
    end if
  end function solve_cpu

end module rlwm_cpu