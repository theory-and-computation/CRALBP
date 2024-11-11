import os

def cutFile(inFile, outFile, num_lines):
    lines = []
    if os.path.isdir(inFile):
        return
    with open(inFile, 'r') as readobj:
        lines = readobj.readlines()
        if len(lines) > num_lines:
            lines = lines[num_lines:]

    with open(outFile, 'w+') as writeobj:
        lines = writeobj.writelines(lines)


def processTrajFiles(inDir, outDir):
    for file_name in os.listdir(inDir):
        inTraj_path = os.path.join(inDir, file_name)
        outTraj_path = os.path.join(outDir, file_name)
        cutFile(inTraj_path, outTraj_path, 250000)


def main():
    processTrajFiles('/home/daniel/Downloads/trajectory_data/CRALBP/sampling/lp_bound/updated_windows', '/home/daniel/Downloads/trajectory_data/CRALBP/sampling/lp_bound/updated_windows_cut')

if __name__ == "__main__":
    main()
