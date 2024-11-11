import os

def read_traj_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return lines

def write_traj_chunk(chunk, index, output_dir):
    file_name = os.path.join(output_dir, f'chunk_{index}.traj')
    with open(file_name, 'w') as file:
        file.writelines(chunk)

def divide_traj_into_blocks(lines, num_blocks=5):
    block_size = len(lines) // num_blocks
    blocks = [lines[i*block_size:(i+1)*block_size] for i in range(num_blocks)]
    return blocks

def process_traj_files(input_dir):
    for file_name in os.listdir(input_dir):
        traj_path = os.path.join(input_dir, file_name)
        traj_lines = read_traj_file(traj_path)
        traj_blocks = divide_traj_into_blocks(traj_lines, num_blocks = 5)

        cut_file_name = file_name[:-13]
        output_dir = os.path.join(input_dir, 'output', cut_file_name)
        os.makedirs(output_dir, exist_ok = True)

        for i, block in enumerate(traj_blocks, 1):
            write_traj_chunk(block, i, output_dir)

def main():
    process_traj_files('/home/daniel/Downloads/trajectory_data/CRALBP/sampling/R234_lp')

if __name__ == "__main__":
    main()
