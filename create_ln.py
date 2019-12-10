import os
import sys,re
dir_path = "/home/zhluo/Project/lung_cancer/data"
sDir = "/home/zhluo/Project/lung_cancer/fastq_ln"
files = os.listdir(dir_path)

for one_file in files:
    regx = re.compile("11100_10495_(.*)_H3KYTBGXC_(.*)_R(.*).fastq.gz")
    result = regx.search(one_file)
    if not result:
        continue
    number = result.groups()[1]
    number_2 = result.groups()[2]
    cmd = "ln -s %s %s" %(os.path.join(dir_path, one_file), os.path.join("/home/zhluo/Project/lung_cancer/fastq_ln", number + "_R" + number_2 + ".fastq.gz"))
    os.system(cmd)