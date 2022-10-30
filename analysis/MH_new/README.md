# 分析流程及脚本介绍

## 分析流程：
1. PDF的产生：
- PDF产生路径：/junofs/users/miaoyu/supernova/wenlj/simulation/examples
在submit路径下的python文件用于批量产生脚本，对于较长的PDF一般分段产生。

- PDF文件路径：/junofs/users/miaoyu/supernova/production/PDFs/10kpc
其中Garching和Burrows分别来自于lhl独立模拟代码和SNEWPY产生子产生的inital flux，放到lhl独立模拟包后的一维时间谱；

2. Data的产生：
- Data产生脚本：/junofs/users/miaoyu/supernova/production/PDFs/10kpc/jobs
使用python脚本产生批量文件，可以指定所用的模型，时间窗口，，距离和文件号码（两种质量顺序和不同探测道不同阈值都一起产生了），每个文件产生1000个事例，如需更大样本则产生多个作业。

- Data存放路径：
/junofs/users/miaoyu/supernova/production/Data
不同距离的两类模型Garching/Burrows，都是基于10kpc的PDF经过scale之后产生的。

3. 分析脚本：

- 路径： /junofs/users/miaoyu/supernova/analysis/MH_new

- 脚本：

        1). 对于每个通道 channel_analyser.py，构造函数里需要传递所有具体的参数；

        2). 联合分析脚本 combining_analysis.py，有多个不同的args可以在运行时往里动态传递参数。
- 批量分析：jobs目录里的py脚本用于批量产生作业和提交

- 结果： results/路径下以csv格式存储。


4. 结果检查：
- 路径：/junofs/users/miaoyu/supernova/analysis/MH_new
- validate.py用于检查一些基础的分布，输出到pdf文件进行浏览；
- visualizer.ipynb进行高级绘图；

最后更新于2022.10.30