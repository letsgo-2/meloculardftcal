# meloculardftcal Agent Interface

[中文接口](AGENT_INTERFACE_CN.md) | [English interface](AGENT_INTERFACE_EN.md)

本文件是语言入口。智能体应根据工作语言读取完整接口：

- 中文：`AGENT_INTERFACE_CN.md`
- English: `AGENT_INTERFACE_EN.md`

机器可读注册表和只读查询入口：

```bash
python3 workflow.py show
python3 workflow.py stages
python3 workflow.py parameters
python3 workflow.py paths --check
```

所有项目内路径均相对于仓库根目录
`https://github.com/letsgo-2/meloculardftcal`。
