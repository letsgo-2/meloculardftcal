#!/bin/bash
set -euo pipefail

#################### 用户可配置区域 ####################

# 是否生成 第一步 opt/freq 的 sbatch 提交脚本：1=生成，0=不生成
GENERATE_SBATCH="${GENERATE_SBATCH:-1}"

# 是否生成 第二步 TD-DFT 的 inp 和 sbatch：1=生成，0=不生成
GENERATE_TD="${GENERATE_TD:-1}"

# sbatch 资源配置（两步统一用这个）
SBATCH_PARTITION="${SBATCH_PARTITION:-amd_512q}"
SBATCH_NODES="${SBATCH_NODES:-1}"
SBATCH_NTASKS="${SBATCH_NTASKS:-32}"          # 建议和 inp 中 %pal nprocs 对齐

# 作业运行的工作目录（提交时会 cd 到这里，再进入 OUTPUT_DIR）
WORK_DIR="${WORK_DIR:-/path/to/dft-workflow}"

# ORCA 二进制路径
ORCA_BIN="${ORCA_BIN:-/path/to/orca/orca}"
MODULE_INIT="${MODULE_INIT:-/path/to/modules/module.sh}"
ORCA_MODULE="${ORCA_MODULE:-orca/6_0_1_linux_x86-64_shared_openmpi416_avx2}"

########################################################

if [[ $# -lt 1 ]]; then
    echo "用法: $0 GJF_DIR [OUTPUT_DIR]" >&2
    exit 1
fi

GJF_DIR=${1%/}                         # 去掉末尾的 /
OUTPUT_DIR=${2:-"${GJF_DIR}inp"}      # 默认输出目录：<输入目录名>inp

mkdir -p "$OUTPUT_DIR"

# 允许没有匹配文件时不报错
shopt -s nullglob

# 同时匹配 .gjf 和 .GJF
for gjf_file in "$GJF_DIR"/*.gjf "$GJF_DIR"/*.GJF; do
    [[ -e "$gjf_file" ]] || continue

    base_name=$(basename "$gjf_file")
    base_name=${base_name%.*}         # 去掉扩展名，如 pic16_0001
    # 如果文件名前缀是 mol，改用目录名作为前缀（如目录 pic3 -> pic3_0001）
    dir_prefix=$(basename "$GJF_DIR")
    if [[ "$base_name" == mol_* ]]; then
        base_name="${dir_prefix}_${base_name#mol_}"
    elif [[ "$base_name" == mol* ]]; then
        base_name="${dir_prefix}_${base_name#mol}"
    fi

    inp_file="$OUTPUT_DIR/${base_name}.inp"

    # 从 gjf 中自动找到第一行“电荷 自旋多重度”
    charge_spin=$(awk '
        NF==2 && $1 ~ /^-?[0-9]+$/ && $2 ~ /^[0-9]+$/ {
            print $1, $2
            exit
        }
    ' "$gjf_file")

    # 坐标开始行 = 上面那行的下一行
    start_coord_line=$(awk '
        NF==2 && $1 ~ /^-?[0-9]+$/ && $2 ~ /^[0-9]+$/ {
            print NR+1
            exit
        }
    ' "$gjf_file")

    if [[ -z "${charge_spin:-}" || -z "${start_coord_line:-}" ]]; then
        echo "警告: 无法在 $gjf_file 中找到 charge/multiplicity 行，跳过。" >&2
        continue
    fi

#     ######## 第一步：生成 ORCA 的 opt/freq .inp 文件 ########
#     {
#         echo "! B3LYP def2-SVP Opt Freq RIJCOSX"
#         echo
#         echo "%maxcore 3000"
#         echo "%pal nprocs ${SBATCH_NTASKS} end"
#         echo "%base \"$base_name\""
#         echo
#         echo
#         echo "* xyz $charge_spin"
#         # 从坐标开始的行复制到文件末尾
#         sed -n "${start_coord_line},\$p" "$gjf_file"
#         echo
#         echo "*"
#     } > "$inp_file"

#     echo "Converted: $(basename "$gjf_file") -> $(basename "$inp_file")"

#     ######## 第一步：生成 sbatch 提交脚本（如果开关开启） ########
#     if [[ "$GENERATE_SBATCH" == "1" ]]; then
#         sbatch_file="$OUTPUT_DIR/${base_name}.sh"

#         cat > "$sbatch_file" <<EOF
# #!/bin/bash
# #SBATCH -J ${base_name}
# #SBATCH -p ${SBATCH_PARTITION}
# #SBATCH -N ${SBATCH_NODES}
# #SBATCH -n ${SBATCH_NTASKS}
# #SBATCH -o ${base_name}.out.%j
# #SBATCH -e ${base_name}.err.%j
# export OMP_NUM_THREADS=1
# source /path/to/modules/module.sh
# module load orca/6_0_1_linux_x86-64_shared_openmpi416_avx2

# cd ${WORK_DIR}/${OUTPUT_DIR}

# ${ORCA_BIN}   ${base_name}.inp > ${base_name}.log
# EOF

#         chmod +x "$sbatch_file"
#         echo "  -> SBATCH script: $(basename "$sbatch_file")"
#     fi

    ######## 第二步：生成 TD-DFT 输入和提交脚本 ########
    if [[ "$GENERATE_TD" == "1" ]]; then
        td_prefix="${base_name}td"            # 如 pic16_0001td
        td_inp_file="$OUTPUT_DIR/${td_prefix}.inp"
        td_sbatch_file="$OUTPUT_DIR/${td_prefix}.sh"
        td_job_name="${base_name}_td"         # 如 pic16_0001_td

        # 生成 TD 的 .inp
        {
            echo "! wB97X-D3 ZORA ZORA-def2-TZVP SARC/J RIJCOSX RI-SOMF(1X) TightSCF CPCM(Water)"
            echo "%tddft"
            echo "    NRoots   10"
            echo "    Triplets true"
            echo "    dosoc true"
            echo "end"
            echo "%basis"
            echo "  NewGTO I \"SARC-ZORA-TZVP\" end"
            echo "end"
            echo "%pal"
            echo "  nprocs ${SBATCH_NTASKS}"
            echo "end"
            echo "%maxcore 3000"
            echo "%base \"${td_prefix}\""
            echo
            echo "%scf"
            echo "  Guess MORead"
            echo "end"
            echo
            echo "%moinp \"${base_name}.gbw\""
            echo
            echo "* xyzfile $charge_spin  ${base_name}.xyz"
            echo
        } > "$td_inp_file"

        echo "  -> TD-DFT input: $(basename "$td_inp_file")"

        # 生成 TD 的 sbatch 提交脚本
        cat > "$td_sbatch_file" <<EOF
#!/bin/bash
#SBATCH -J ${td_job_name}
#SBATCH -p ${SBATCH_PARTITION}
#SBATCH -N ${SBATCH_NODES}
#SBATCH -n ${SBATCH_NTASKS}
#SBATCH -o ${td_prefix}.out.%j
#SBATCH -e ${td_prefix}.err.%j
export OMP_NUM_THREADS=1
source ${MODULE_INIT}
module load ${ORCA_MODULE}

cd ${WORK_DIR}/${OUTPUT_DIR}

${ORCA_BIN}  ${td_prefix}.inp > ${td_prefix}.log
EOF

        chmod +x "$td_sbatch_file"
        echo "  -> TD SBATCH script: $(basename "$td_sbatch_file")"
    fi

done

echo "Conversion complete! Files saved in $OUTPUT_DIR/"
