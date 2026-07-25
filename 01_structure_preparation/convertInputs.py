"""
refactored: getGaussianInput.py ver 1.0 + batch conversions of both DA and DAD to GJF + runBatchConvertLog

Dependency: Open Babel's command prompt scripts 
- Open Babel 3.1 https://anaconda.org/conda-forge/openbabel; remember to set path to Anaconda 

Two processes here to define as arguments in CLI
1. SMILES --> GJF and TXT for structural optimizations
2. LOG --> GJF and TXT for TD DFT

"""
import pandas as pd 
import os 
import sys
from pathlib import Path

def read_smiles_simple(path):
    smiles = []
    with open(path, 'r', encoding='utf-8') as f:
        next(f)  # 跳过表头
        for lineno, line in enumerate(f, start=2):
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            smiles.append(parts[1])
            print(parts[1])
    return smiles

def read_smiles_pandas(path):
    # 自动按空白分隔，要求第一行为表头包含 Index 和 SMILES（或第二列为 SMILES）
    df = pd.read_csv(path, sep=r'\s+', engine='python')
    # 如果有名为 SMILES 的列就直接用它；否则使用第二列
    if 'SMILES' in df.columns:
        return df['SMILES'].astype(str).tolist()
    else:
        return df.iloc[:, 1].astype(str).tolist()
# 1. Convert SMILES to gjf files with 3d coordinates first 
def convert_smiles_to_gjf(smiles_list, molecule_name_list,calc_dir):
    for i in range(len(smiles_list)):
        smi = smiles_list[i]
        filename = molecule_name_list[i]
        print("Currently converting SMILES {} to file {}.gjf in {}".format(smi, filename, calc_dir))
        #控制在calc文件夹下执行
        # sys_text = 'obabel -:"{}" -O {}.gjf --gen3d --ff uff --minimize --steps 500'.format(smi, Path.joinpath(calc_dir,filename))
        sys_text = 'obabel -:"{}" -O {}.gjf --gen3d'.format(smi, Path.joinpath(calc_dir,filename))

        # os.chdir(calc_dir)
        os.system(sys_text)
        # os.chdir(os.pardir)
    print("Initial conversion to GJFs DONE!\nTotal number of molecule files in this batch: {}\n\n".format(len(smiles_list)))
    return  

# 1. Convert log files to .gjf again using Open Babel, but now for TD-DFT calculations 
def convert_log_to_gjf_td(log_list,calc_dir):
    """
    Given log list in the directory, convert them to _td.gjf files 
    Then also return the molecule_name_list to be used for later steps
    """
    molecule_name_list = []
    for i in range(len(log_list)):
        filename = log_list[i][:-4] # DA_24877.log to DA_24877 is filename
        newfilename = filename + '_td' # append _td to filename
        molecule_name_list.append(newfilename)
        print("Currently converting {}.log to {}.gjf".format(filename, newfilename))
        sys_text = 'obabel "{}".log -O {}.gjf'.format(Path(calc_dir).joinpath(filename), Path(calc_dir).joinpath(newfilename)) 
        os.system(sys_text)
    print("Initial conversion to TD GJFs DONE!\nTotal number of molecule files in this batch: {}\n\n".format(len(log_list)))
    return molecule_name_list

# 2. Edit gjf files to be runnable and create txt job file for every molecule 
# refactored version for both structural opt and td dft gjf files
def edit_gjf_file(filename, process = None,calc_dir=None):
    """
    process
    - SO for structural opt
    - TD for TD-DFT 
    - TDlink for TD-DFT with Guess=Read (no molecular structure needed)
    """
    # 对于TDlink模式，不需要读取现有的gjf文件，直接创建新的
    if process == 'TDlink':
        # 直接创建新的gjf文件，不需要分子结构（使用Guess=Read）
        pass
    else:
        # 对于SO和TD模式，需要读取现有的gjf文件并处理
        with open('{}.gjf'.format(Path(calc_dir).joinpath(filename)), 'r') as f:
            lines = f.readlines()
        with open('{}.gjf'.format(Path(calc_dir).joinpath(filename)), 'w') as f:
            for line in lines[4:]: # only take the lines from line 4 onwards and remove the rest on top
                f.write(line)

    # lines to be added to the top of every gjf file for the job to run subsequently
    if process == 'SO':
        gjf_lines_to_add = ["%chk={}.chk".format(Path(calc_dir).joinpath(filename)), 
                            "%mem=80000MB",  # 使用大部分内存
                            "%nproc=56",     # 使用全部CPU核心
                            "#p b3lyp/6-31g(d,p) opt freq scrf=(solvent=water)", 
                            "",
                            "Title Card Required"] 
    # #p b3lyp/6-31g(d) opt
    elif process == 'TD':
        gjf_lines_to_add = ["%chk={}.chk".format(Path(calc_dir).joinpath(filename)),
                            "%mem=80000MB",  # 使用大部分内存
                            "%nproc=56",     # 使用全部CPU核心
                            "#p b3lyp/6-31g(d,p) td=(nstates=6, 50-50) em=GD3BJ scrf=(solvent=water)",
                            "",
                            "Title Card Required"]
    elif process == 'TDlink':
        gjf_lines_to_add = ["%oldchk={}.chk".format(Path(calc_dir).joinpath(filename[:-3])),
                            "%chk={}.chk".format(Path(calc_dir).joinpath(filename)),
                            "%mem=80000MB",  # 使用大部分内存
                            "%nproc=56",     # 使用全部CPU核心
                            "#p b3lyp/6-31g(d,p) td=(nstates=6, 50-50) opt Guess=Read em=GD3BJ scrf=(solvent=water)",
                            "",
                            "Title Card Required"]
                            
    else:
        print("No such method for GJF defined.")
        return
        
    def prepend_multiple_lines(file_name, list_of_lines, process_mode=None): 
        # to add lines to the top of every gjf file 
        dummy_file = file_name + '.gjf'
        
        if process_mode == 'TDlink':
            # 对于TDlink模式，从原始gjf文件读取电荷和自旋多重度
            original_filename = filename[:-3]  # 去掉_td后缀，得到原始文件名
            original_gjf_path = Path(calc_dir).joinpath(original_filename + '.gjf')
            
            charge_spin = "0 1"  # 默认值
            if original_gjf_path.exists():
                try:
                    with open(original_gjf_path, 'r') as f:
                        lines = f.readlines()
                        # 查找电荷和自旋多重度行（通常在分子结构之前）
                        for line in lines:
                            line = line.strip()
                            if line and not line.startswith('%') and not line.startswith('#') and not line.startswith('Title'):
                                # 检查是否是电荷和自旋多重度行（两个数字）
                                parts = line.split()
                                if len(parts) == 2 and parts[0].replace('-', '').replace('+', '').isdigit() and parts[1].isdigit():
                                    charge_spin = line
                                    break
                except Exception as e:
                    print(f"[WARN] Could not read charge/spin from {original_gjf_path}: {e}")
                    print(f"[WARN] Using default charge/spin: {charge_spin}")
            else:
                print(f"[WARN] Original gjf file not found: {original_gjf_path}")
                print(f"[WARN] Using default charge/spin: {charge_spin}")
            
            # 创建完整的gjf文件
            with open(dummy_file, 'w') as write_obj:
                for line in list_of_lines:
                    write_obj.write(line + '\n')
                # 添加空行和电荷自旋多重度
                write_obj.write('\n')
                write_obj.write(f'{charge_spin}\n')
                write_obj.write('\n')     # 空行结束
        else:
            # 对于其他模式，读取现有文件并添加头部
            with open(file_name, 'r') as read_obj, open(dummy_file, 'w') as write_obj:
                for line in list_of_lines:
                    write_obj.write(line + '\n')
                for line in read_obj:
                    write_obj.write(line)
            os.remove(file_name)
        
        os.rename(dummy_file, file_name)
        return

    prepend_multiple_lines('{}.gjf'.format(Path(calc_dir).joinpath(filename)), gjf_lines_to_add, process)
    return

# 3. Create a job sh file for SLURM 
def create_job_sh(filename,calc_dir):
    # SBATCH header for SLURM
    sbatch_lines = [
        '#!/bin/bash\n',
        '#SBATCH -J gaussian_{}\n'.format(filename),
        '#SBATCH -p cnall\n',
        '#SBATCH -N 1\n',
        '#SBATCH -o {}.out.%j\n'.format(filename),
        '#SBATCH -e {}.err.%j\n'.format(filename),
        '#SBATCH --no-requeue\n',
        '#SBATCH --ntasks-per-node=1\n',   # 每个节点运行1个作业
        '#SBATCH --cpus-per-task=56\n',    # 每个作业使用全部56个CPU核心
        '\n',
        f'cd {calc_dir}\n',
        '\n',
        'g16 <{}.gjf > {}.log\n'.format(filename, filename),
        '\n',
        '##--- END HERE ---\n'
    ]
    
    #filename在calc_dir下
    filename1 = Path(calc_dir).joinpath(filename)
    out_filename = filename1.with_suffix('.sh')
    with open(out_filename, 'w') as f:
        for line in sbatch_lines:
            f.write(line)
    return

# Together: combining both gjf and txt steps
# given molecule name for the generated molecule, get completed gjf and txt files in directory
def generate_gjf_sh_batch(molecule_name_list, process = None,calc_dir=None):
    if calc_dir is None:
        print("No directory provided. Please provide a calc directory.")
        return
    for filename in molecule_name_list:
        create_job_sh(filename,calc_dir)
        if process in ["SO", "TD", "TDlink"]:
            edit_gjf_file(filename, process,calc_dir)
        else:
            print("No such method for GJF conversion.")
            return
        print('Both gjf and txt files created for: {} !'.format(filename))
    return 


#################################################

# command line executions 
"""
python convertGaussianInputs.py X 
- X = 1: Convert SMILES from CSV in the directory to GJF and TXT files
- X = 2: Convert all LOG files in the directory to td.GJF and td.TXT files

"""
def main():
    if len(sys.argv) < 3:
        print(
            "Usage: python convertInputs.py <mode> <csv_path> [output_dir]\n"
            "Example: python convertInputs.py 1 data/raw/pic16.csv work/pic16\n"
            "Modes: 1=SO, 2=TD, 3=TDlink"
        )
        return

    process = int(sys.argv[1])
    csv_filename = sys.argv[2]
    #判断两个参数，若错误给示例
    if process not in [1,2,3]:
        print("Please provide a valid process. Example: python convertInputs.py 1 data.csv")
        print("Process options: 1=SO, 2=TD, 3=TDlink")
        return
    # 仓库根目录。脚本位于 01_structure_preparation/ 下。
    WORK_DIR = Path(__file__).resolve().parents[1]

    # 优先使用命令行给出的仓库相对路径；为兼容旧用法，再尝试 data/<文件名>。
    csv_arg = Path(csv_filename)
    csv_candidates = (
        [csv_arg]
        if csv_arg.is_absolute()
        else [WORK_DIR / csv_arg, WORK_DIR / "data" / csv_arg]
    )
    csv_path_file = next((path for path in csv_candidates if path.exists()), None)
    if csv_path_file is None:
        print("CSV file not found. Checked:")
        for path in csv_candidates:
            print("  {}".format(path))
        return

    # 可选第三参数允许把结构直接写入 work/<batch>；未提供时保持旧行为。
    if len(sys.argv) >= 4:
        output_arg = Path(sys.argv[3])
        calc_dir = output_arg if output_arg.is_absolute() else WORK_DIR / output_arg
    else:
        calc_dir = WORK_DIR / csv_path_file.stem

    calc_dir.mkdir(parents=True, exist_ok=True)
    #说明处理阶段
    print("process: ",process)
    if process == 1:
        print("Get files for Structural Optimization: Converting SMILES in CSV to GJF and TXT files...")
        sample_df = pd.read_csv(csv_path_file)
        print("CSV selected: {}".format(csv_path_file))
        #我可能没有MOLNAME
        if 'MOLNAME' in sample_df.columns:
            names = list(sample_df['MOLNAME'].astype(str))
        else:
            # 没有 MOLNAME 就自动编号：mol_0001, mol_0002, ...
            names = [f"mol_{i:04d}" for i in range(1, len( sample_df['SMILES']) + 1)]
        #创建一个文件夹，用csv的文件名作为后续计算目录
        
        convert_smiles_to_gjf(list(sample_df['SMILES']), names,calc_dir)
        generate_gjf_sh_batch(names, process = "SO",calc_dir=calc_dir)
        print("DONE! All GJF and TXT files are generated in this directory.")

    elif process == 2:
        print("Get files for TD DFT: Converting all LOG in directory to td.GJF and td.TXT files...")
        #需要再calc_dir下来listdir
        dir_log_list = []
        for file in os.listdir(calc_dir):
            if file.endswith(".log") and not file.endswith(".td.log"):
                #只要文件名，不要路径
                dir_log_list.append(file)
        print("Total LOG files in this directory to be converted: ", len(dir_log_list))
        molecule_name_list = convert_log_to_gjf_td(dir_log_list,calc_dir)
        generate_gjf_sh_batch(molecule_name_list, process = "TD",calc_dir=calc_dir)
        print("DONE! All td.GJF and td.TXT files are generated in this directory.")
    elif process == 3:
        print("Get files for TD DFT with Guess=Read: Creating td.GJF and td.SH files using existing CHK files...")
        #需要再calc_dir下来listdir，找到所有.log文件（排除.td.log）
        dir_log_list = []
        for file in os.listdir(calc_dir):
            if file.endswith(".log") and not file.endswith(".td.log"):
                #只要文件名，不要路径
                dir_log_list.append(file)
        print("Total LOG files in this directory to be converted: ", len(dir_log_list))
        
        # 直接生成molecule_name_list，不需要用obabel转换
        molecule_name_list = []
        for log_file in dir_log_list:
            filename = log_file[:-4]  # 去掉.log后缀
            newfilename = filename + '_td'  # 添加_td后缀
            molecule_name_list.append(newfilename)
            print("Creating TD-DFT files for: {} (using Guess=Read from {}.chk)".format(newfilename, filename))
        
        # 直接生成gjf和sh文件，不需要obabel转换
        generate_gjf_sh_batch(molecule_name_list, process = "TDlink",calc_dir=calc_dir)
        print("DONE! All td.GJF and td.SH files are generated in this directory.")

    else:
        print("There is no such method. Please recheck script.")

if __name__ == "__main__":
    main()
