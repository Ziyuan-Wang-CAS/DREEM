import argparse
import os
import shutil

def copy_control_stats(exp_date, exp_name, control_date, control_name):

    base_dir = "/sibcb1/hanshuolab1/wangziyuan/DREEM/"
    base_exp_folder = f"results_{exp_date}_{exp_name}"
    base_control_folder = f"results_{control_date}_{control_name}"

    for item in os.listdir(base_dir):  # Iterate through items in base directory
        if base_exp_folder == item and os.path.isdir(os.path.join(base_dir, item)):
            exp_base_path = os.path.join(base_dir, item)

            # Search for analysis region within the base experimental folder
            for sub_item in os.listdir(exp_base_path):
                if base_exp_folder in sub_item and os.path.isdir(os.path.join(exp_base_path, sub_item)):
                    exp_path = os.path.join(exp_base_path, sub_item)
                    analysis_region = sub_item.replace(base_exp_folder + '_', '')  # Extract analysis region

                    # Search for control group folder within the base directory
                    control_path = None
                    for control_item in os.listdir(base_dir):
                        if base_control_folder == control_item and os.path.isdir(os.path.join(base_dir, control_item)):
                            control_base_path = os.path.join(base_dir, control_item)
                            # 
                            for control_sub_item in os.listdir(control_base_path):
                                if base_control_folder in control_sub_item and os.path.isdir(os.path.join(control_base_path, control_sub_item)) and analysis_region in control_sub_item:
                                    control_path = os.path.join(control_base_path, control_sub_item)
                                    break
                            if control_path:
                                break
                    if not control_path:
                        continue

                    control_stats_file = os.path.join(control_path, "control_group_stats.json")
                    stats_file_destination = os.path.join(exp_path, "control_group_stats.json")

                    if not os.path.exists(control_stats_file):
                        print(f"[error] Control stats file not found: {control_stats_file}")
                        continue
                    if not os.path.exists(exp_path):
                        print(f"[error] Experimental folder not found: {exp_path}")
                        continue

                    try:
                        shutil.copy2(control_stats_file, stats_file_destination)  # Copy control stats file to experimental folder
                        rel_control_stats_file = os.path.relpath(control_stats_file, base_dir)
                        rel_stats_file_destination = os.path.relpath(stats_file_destination, base_dir)
                        print(f"[INFO] Copied: {rel_control_stats_file} -> {rel_stats_file_destination}\n")
                    except Exception as e:
                        print(f"[error] Failed to copy {rel_control_stats_file} to {rel_stats_file_destination}: {e}\n")
            break

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Copies control_group_stats.json to multiple experimental groups.")
    parser.add_argument("exp_date", help="Analysis date of experimental group (YYYYMMDD)")
    parser.add_argument("exp_name", help="Data name of experimental group")
    parser.add_argument("control_date", help="Analysis date of control group (YYYYMMDD)")
    parser.add_argument("control_name", help="Data name of control group")

    args = parser.parse_args()
    copy_control_stats(args.exp_date, args.exp_name, args.control_date, args.control_name)