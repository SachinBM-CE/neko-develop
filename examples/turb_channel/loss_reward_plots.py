import re
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

def plot_rlwm_metrics(log_data, log_file_path):
    """
    Parses rlwm model log data to plot losses, rewards, and tau values
    and saves the plots to the directory of the log file.

    Args:
        log_data (str): A string containing the log file content.
        log_file_path (str): The path to the log file.
    """
    # Dictionaries to store parsed data from the 'rlwm' model
    rlwm_actor_loss = {}
    rlwm_critic_loss = {}
    rlwm_rew_out = {}
    rlwm_tau_new = {}
    rlwm_reward_sum = {}

    # Read and parse the log data line by line
    for line in log_data.splitlines():
        # Parse actor loss (model: actor with train_loss)
        if "model: actor" in line and "train_loss:" in line:
            match_actor_loss = re.search(r"step: (\d+), train_loss: ([\deE\+\-\.]+)", line)
            if match_actor_loss:
                step, loss = int(match_actor_loss.group(1)), float(match_actor_loss.group(2))
                rlwm_actor_loss[step] = loss
        
        # Parse critic loss (model: critic_0 with train_loss)
        elif "model: critic_0" in line and "train_loss:" in line:
            match_critic_loss = re.search(r"step: (\d+), train_loss: ([\deE\+\-\.]+)", line)
            if match_critic_loss:
                step, loss = int(match_critic_loss.group(1)), float(match_critic_loss.group(2))
                rlwm_critic_loss[step] = loss
        
        # Parse rew_out (if it exists in your logs)
        elif "model: rlwm" in line and "rew_out:" in line:
            match_rew_out = re.search(r"step: (\d+), rew_out: ([\deE\.\+\-]+)", line)
            if match_rew_out:
                step, rew = int(match_rew_out.group(1)), float(match_rew_out.group(2))
                rlwm_rew_out[step] = rew
        
        # Parse tau_new (if it exists in your logs)
        elif "model: rlwm" in line and "tau_new:" in line:
            match_tau_new = re.search(r"step: (\d+), tau_new: ([\deE\.\+\-]+)", line)
            if match_tau_new:
                step, tau = int(match_tau_new.group(1)), float(match_tau_new.group(2))
                rlwm_tau_new[step] = tau
        
        # Parse reward_sum (model: rlwm)
        elif "model: rlwm" in line and "reward_sum:" in line:
            match_reward_sum = re.search(r"step: (\d+), reward_sum: ([\deE\.\+\-]+)", line)
            if match_reward_sum:
                step, reward_sum = int(match_reward_sum.group(1)), float(match_reward_sum.group(2))
                rlwm_reward_sum[step] = reward_sum / 12800 # / (12800*100)
    
    # Debug: Print what we found
    print(f"Actor loss entries found: {len(rlwm_actor_loss)}")
    print(f"Critic loss entries found: {len(rlwm_critic_loss)}")
    print(f"Reward entries found: {len(rlwm_reward_sum)}")
    print(f"Rew_out entries found: {len(rlwm_rew_out)}")
    print(f"Tau_new entries found: {len(rlwm_tau_new)}")
    
    # Get the directory path from the log file path
    log_dir = os.path.dirname(log_file_path)
    
    # --- Plot 1: Actor vs. Critic Loss for rlwm model ---
    common_loss_steps = sorted(list(set(rlwm_actor_loss.keys()) & set(rlwm_critic_loss.keys())))
    if common_loss_steps:
        steps = common_loss_steps
        actor_loss = [rlwm_actor_loss[step] for step in common_loss_steps]
        critic_loss = [rlwm_critic_loss[step] for step in common_loss_steps]

        fig, ax1 = plt.subplots(figsize=(12, 6))
        ax1.set_xlabel('Time Step')
        ax1.set_ylabel('Actor Loss', color='tab:red')
        ax1.plot(steps, actor_loss, color='tab:red', label='Actor Loss', markersize=4)
        ax1.tick_params(axis='y', labelcolor='tab:red')
        ax1.grid(True, alpha=0.3)
        
        ax2 = ax1.twinx()
        ax2.set_ylabel('Critic Loss', color='tab:blue')
        ax2.plot(steps, critic_loss, color='tab:blue', label='Critic Loss', markersize=4)
        ax2.tick_params(axis='y', labelcolor='tab:blue')

        plt.title('Actor Loss vs. Critic Loss Over Time')
        fig.tight_layout()
        
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(100))
            
        output_filename = os.path.join(log_dir, "rlwm_actor_vs_critic_loss.png")
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Plot successfully saved to {os.path.abspath(output_filename)}")
        plt.close(fig)
    else:
        print("No common steps found between actor and critic loss data")

    # --- Plot 2: Reward vs. Tau for rlwm model ---
    common_reward_steps = sorted(list(set(rlwm_rew_out.keys()) & set(rlwm_tau_new.keys())))
    if common_reward_steps:
        steps = common_reward_steps
        rew_out = [rlwm_rew_out[step] for step in common_reward_steps]
        tau_new = [rlwm_tau_new[step] for step in common_reward_steps]

        fig, ax1 = plt.subplots(figsize=(12, 6))
        ax1.set_xlabel('Time Step')
        ax1.set_ylabel('rew_out', color='tab:green')
        ax1.plot(steps, rew_out, color='tab:green', label='rew_out', markersize=4)
        ax1.tick_params(axis='y', labelcolor='tab:green')
        ax1.grid(True, alpha=0.3)

        ax2 = ax1.twinx()
        ax2.set_ylabel('tau_new', color='tab:orange')
        ax2.plot(steps, tau_new, color='tab:orange', label='tau_new', markersize=4)
        ax2.tick_params(axis='y', labelcolor='tab:orange')

        plt.title('Reward vs. Tau Over Time')
        fig.tight_layout()
        output_filename = os.path.join(log_dir, "rlwm_rewards_and_tau.png")
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Plot successfully saved to {os.path.abspath(output_filename)}")
        plt.close(fig)
    else:
        print("No rew_out or tau_new data found in logs")
        
    # --- Plot 3: mean reward per episode for a single agent ---
    if rlwm_reward_sum:
        steps = sorted(rlwm_reward_sum.keys())
        reward_sums = [rlwm_reward_sum[step] for step in steps]

        fig, ax = plt.subplots(figsize=(12, 6))
        ax.set_xlabel('Episode')
        ax.set_ylabel('Mean Reward of a single agent')
        ax.plot(steps, reward_sums, color='tab:purple', label='Mean Reward', markersize=6)
        ax.grid(True, alpha=0.3)
        ax.legend()

        plt.title('Mean reward per episode for a single agent')
        fig.tight_layout()
        output_filename = os.path.join(log_dir, "rlwm_cumulative_reward_sum.png")
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Plot successfully saved to {os.path.abspath(output_filename)}")
        plt.close(fig)
    else:
        print("No reward_sum data found in logs")
    

if __name__ == "__main__":
    log_file_path = "/home/sachinbm/cpu-neko-tf/neko-develop/examples/turb_channel/01_logdir/torchfort.log"
    try:
        with open(log_file_path, "r") as f:
            log_data_string = f.read()
        plot_rlwm_metrics(log_data_string, log_file_path)
    except FileNotFoundError:
        print(f"Error: The file '{log_file_path}' was not found.")