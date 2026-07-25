#!/usr/bin/env python3
"""
从指定的CSV文件读取SMILES，生成分子结构图片
图片命名格式：pic16_0001.png

用法:
    # 使用RDKit（默认），输出到 data/image/pic16/
    python generate_molecule_images.py data/pic16.csv
    
    # 使用RDKit，指定输出目录
    python generate_molecule_images.py data/pic16.csv --tool rdkit --output custom_dir
    
    # 使用Open Babel，输出到 data/image/pic16ob/（自动加ob后缀）
    python generate_molecule_images.py data/pic16.csv --tool openbabel
    
    # 使用Open Babel，指定输出目录
    python generate_molecule_images.py data/pic16.csv --tool openbabel --output custom_dir
    
    # 测试模式：只生成指定序号的图片（用于调试）
    python generate_molecule_images.py data/pic16.csv --tool openbabel --test 1
    python generate_molecule_images.py data/pic16.csv --tool rdkit --test 5

输出目录规则:
    - RDKit（默认）: data/image/pic16/
    - Open Babel（未指定输出）: data/image/pic16ob/
    - 指定输出目录: 使用指定的目录/pic16/
"""

import os
import sys
import argparse
import pandas as pd
import subprocess
from pathlib import Path

# 尝试导入RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("警告: RDKit未安装，将无法使用rdkit工具")

def check_obabel_command():
    """检查obabel命令行工具是否可用"""
    try:
        # obabel --version 会输出版本信息到stdout，但返回非零退出码
        # 所以检查stdout中是否有版本信息即可
        result = subprocess.run(['obabel', '--version'], 
                              capture_output=True, 
                              text=True, 
                              timeout=5)
        # 检查stdout或stderr中是否有"Open Babel"版本信息
        if 'Open Babel' in result.stdout or 'Open Babel' in result.stderr:
            return True
        return False
    except (FileNotFoundError, subprocess.TimeoutExpired):
        # 如果命令不存在，会抛出FileNotFoundError
        return False
    except Exception:
        # 其他异常也认为不可用
        return False

def generate_molecule_image_rdkit(smiles, output_path, size=(400, 400)):
    """
    使用RDKit从SMILES生成分子结构图片
    
    Args:
        smiles: SMILES字符串
        output_path: 输出图片路径
        size: 图片尺寸 (width, height)
    """
    if not RDKIT_AVAILABLE:
        return False
    
    try:
        # 从SMILES创建分子对象
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"警告: 无法解析SMILES: {smiles}")
            return False
        
        # 生成2D结构图片
        img = Draw.MolToImage(mol, size=size)
        img.save(output_path)
        return True
    except Exception as e:
        print(f"错误: RDKit生成图片失败 {output_path}: {e}")
        return False

def generate_molecule_image_openbabel(smiles, output_path, size=(400, 400)):
    """
    使用Open Babel命令行工具从SMILES生成分子结构图片
    
    Args:
        smiles: SMILES字符串
        output_path: 输出图片路径
        size: 图片尺寸 (width, height)
    """
    try:
        # 使用obabel命令行工具
        # -:SMILES 格式输入SMILES，-O 指定输出文件，--gen2d生成2D坐标，-xp设置图片宽度
        cmd = [
            'obabel',
            '-:' + smiles,  # 输入SMILES（使用-:格式）
            '-O', str(output_path),  # 输出文件
            '--gen2d',  # 生成2D坐标
            '-xp', str(size[0]),  # 图片宽度（像素）
        ]
        
        # 打印调试信息（仅在测试模式或详细模式下）
        if '--test' in sys.argv or '-v' in sys.argv:
            print(f"执行命令: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, 
                              capture_output=True, 
                              text=True, 
                              timeout=10)
        
        # 打印命令输出用于调试（仅在测试模式）
        if '--test' in sys.argv:
            if result.stdout:
                print(f"obabel stdout: {result.stdout}")
            if result.stderr:
                print(f"obabel stderr: {result.stderr}")
            print(f"返回码: {result.returncode}")
        
        if result.returncode == 0 and os.path.exists(output_path):
            # 检查文件大小，确保不是空文件（至少应该大于1KB才算有效）
            file_size = os.path.getsize(output_path)
            if '--test' in sys.argv:
                print(f"生成的文件大小: {file_size} 字节")
            
            # PNG文件至少应该有几百字节，如果太小可能是错误信息
            if file_size > 1000:  # 至少1KB
                return True
            else:
                if '--test' in sys.argv:
                    print(f"警告: 生成的图片文件太小（{file_size}字节），可能是错误")
                    # 尝试读取文件内容看看是什么
                    try:
                        with open(output_path, 'rb') as f:
                            content = f.read(200)
                            print(f"文件内容预览: {content[:100]}")
                    except:
                        pass
                return False
        else:
            if '--test' in sys.argv:
                print(f"错误: obabel命令执行失败")
                if result.stderr:
                    print(f"错误信息: {result.stderr}")
            return False
    except FileNotFoundError:
        print("错误: 找不到obabel命令，请确保Open Babel已安装并在PATH中")
        return False
    except Exception as e:
        if '--test' in sys.argv:
            print(f"错误: Open Babel命令行生成图片失败 {output_path}: {e}")
            import traceback
            traceback.print_exc()
        return False

def generate_molecule_image(smiles, output_path, tool='rdkit', size=(400, 400)):
    """
    从SMILES生成分子结构图片（根据工具选择）
    
    Args:
        smiles: SMILES字符串
        output_path: 输出图片路径
        tool: 使用的工具 ('rdkit' 或 'openbabel')
        size: 图片尺寸 (width, height)
    """
    if tool == 'rdkit':
        return generate_molecule_image_rdkit(smiles, output_path, size)
    elif tool == 'openbabel':
        # 只使用命令行工具
        if check_obabel_command():
            return generate_molecule_image_openbabel(smiles, output_path, size)
        else:
            print("错误: Open Babel命令行工具不可用，请安装: apt-get install openbabel")
            return False
    else:
        print(f"错误: 未知的工具: {tool}")
        return False

def process_csv_file(csv_file_path, output_base_dir=None, tool='rdkit'):
    """
    处理指定的CSV文件
    
    Args:
        csv_file_path: CSV文件路径
        output_base_dir: 输出基础目录，如果为None则使用data/image
        tool: 使用的工具 ('rdkit' 或 'openbabel')
    """
    # 检查文件是否存在
    csv_path = Path(csv_file_path)
    if not csv_path.exists():
        print(f"错误: 文件不存在: {csv_file_path}")
        return
    
    # 获取文件名（不含扩展名），如 pic16
    file_base = csv_path.stem
    
    # 确定输出目录
    if output_base_dir is None:
        # 如果CSV文件在data目录下，使用data/image/文件名
        if csv_path.parent.name == "data":
            base_output = csv_path.parent / "image" / file_base
        else:
            # 否则在CSV文件所在目录创建image/文件名
            base_output = csv_path.parent / "image" / file_base
        
        # 如果使用Open Babel且没有指定输出目录，在默认目录后加"ob"
        if tool == 'openbabel':
            output_dir = Path(str(base_output) + 'ob')
        else:
            output_dir = base_output
    else:
        output_dir = Path(output_base_dir) / file_base
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"处理文件: {csv_path.name}")
    print(f"使用工具: {tool}")
    print(f"输出目录: {output_dir}")
    
    # 检查工具可用性
    if tool == 'rdkit' and not RDKIT_AVAILABLE:
        print("错误: RDKit不可用，请安装: pip install rdkit")
        return
    elif tool == 'openbabel':
        if not check_obabel_command():
            print("错误: Open Babel命令行工具不可用，请安装: apt-get install openbabel")
            return
    
    try:
        # 读取CSV文件
        df = pd.read_csv(csv_path)
        
        # 检查是否有SMILES列
        if 'SMILES' not in df.columns:
            print(f"错误: {csv_path.name} 中没有SMILES列")
            return
        
        total_images = 0
        failed_images = 0
        
        # 处理每一行的SMILES
        for idx, row in df.iterrows():
            smiles = row['SMILES']
            
            if pd.isna(smiles) or not smiles:
                continue
            
            # 生成图片文件名，格式：pic16_0001.png
            image_filename = f"{file_base}_{idx+1:04d}.png"
            output_path = output_dir / image_filename
            
            # 生成图片
            if generate_molecule_image(smiles, str(output_path), tool=tool):
                total_images += 1
                if (idx + 1) % 50 == 0:
                    print(f"  已处理 {idx + 1} 个分子...")
            else:
                failed_images += 1
        
        print(f"\n完成!")
        print(f"成功生成 {total_images} 张图片")
        if failed_images > 0:
            print(f"失败 {failed_images} 张图片")
        print(f"图片保存在: {output_dir}")
                    
    except Exception as e:
        print(f"错误: 处理文件 {csv_path.name} 时出错: {e}")
        import traceback
        traceback.print_exc()

def test_single_image(csv_file_path, index, tool='rdkit', output_base_dir=None):
    """
    测试生成单个指定序号的图片
    
    Args:
        csv_file_path: CSV文件路径
        index: 要生成的图片序号（从1开始）
        tool: 使用的工具 ('rdkit' 或 'openbabel')
        output_base_dir: 输出基础目录
    """
    # 检查文件是否存在
    csv_path = Path(csv_file_path)
    if not csv_path.exists():
        print(f"错误: 文件不存在: {csv_file_path}")
        return
    
    # 获取文件名（不含扩展名），如 pic16
    file_base = csv_path.stem
    
    # 确定输出目录
    if output_base_dir is None:
        if csv_path.parent.name == "data":
            base_output = csv_path.parent / "image" / file_base
        else:
            base_output = csv_path.parent / "image" / file_base
        
        if tool == 'openbabel':
            output_dir = Path(str(base_output) + 'ob')
        else:
            output_dir = base_output
    else:
        output_dir = Path(output_base_dir) / file_base
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"测试模式: 生成第 {index} 个分子的图片")
    print(f"处理文件: {csv_path.name}")
    print(f"使用工具: {tool}")
    print(f"输出目录: {output_dir}")
    
    # 检查工具可用性
    if tool == 'rdkit' and not RDKIT_AVAILABLE:
        print("错误: RDKit不可用，请安装: pip install rdkit")
        return
    elif tool == 'openbabel':
        if not check_obabel_command():
            print("错误: Open Babel命令行工具不可用，请安装: apt-get install openbabel")
            return
    
    try:
        # 读取CSV文件
        df = pd.read_csv(csv_path)
        
        # 检查是否有SMILES列
        if 'SMILES' not in df.columns:
            print(f"错误: {csv_path.name} 中没有SMILES列")
            return
        
        # 检查序号是否有效（转换为0-based索引）
        if index < 1 or index > len(df):
            print(f"错误: 序号 {index} 超出范围（1-{len(df)}）")
            return
        
        idx = index - 1  # 转换为0-based索引
        row = df.iloc[idx]
        smiles = row['SMILES']
        
        if pd.isna(smiles) or not smiles:
            print(f"错误: 第 {index} 行的SMILES为空")
            return
        
        print(f"\nSMILES: {smiles}")
        
        # 生成图片文件名
        image_filename = f"{file_base}_{index:04d}.png"
        output_path = output_dir / image_filename
        
        print(f"输出文件: {output_path}")
        print("-" * 50)
        
        # 生成图片
        if generate_molecule_image(smiles, str(output_path), tool=tool):
            print(f"\n✓ 成功生成图片: {output_path}")
            if os.path.exists(output_path):
                file_size = os.path.getsize(output_path)
                print(f"  文件大小: {file_size} 字节")
        else:
            print(f"\n✗ 生成图片失败")
                    
    except Exception as e:
        print(f"错误: 处理文件 {csv_path.name} 时出错: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='从CSV文件读取SMILES并生成分子结构图片')
    parser.add_argument('csv_file', help='CSV文件路径，例如: data/pic16.csv')
    parser.add_argument('--output', '-o', default=None, help='输出基础目录 (默认: data/image)')
    parser.add_argument('--tool', '-t', 
                       choices=['rdkit', 'openbabel'], 
                       default='rdkit',
                       help='使用的工具: rdkit 或 openbabel (默认: rdkit)')
    parser.add_argument('--test', type=int, metavar='INDEX',
                       help='测试模式: 只生成指定序号的图片（从1开始），例如: --test 1')
    
    args = parser.parse_args()
    
    if args.test:
        # 测试模式：只生成指定序号的图片
        test_single_image(args.csv_file, args.test, args.tool, args.output)
    else:
        # 正常模式：生成所有图片
        process_csv_file(args.csv_file, args.output, args.tool)

