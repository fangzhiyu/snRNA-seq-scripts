# -- coding: utf-8 -
import os
import time

def traverse_directory(directory):
    # 获取当前文件夹下的所有文件和子文件夹
    items = os.listdir(directory)
    for item in items:
        item_path = os.path.join(directory, item)
        if os.path.isfile(item_path):
            # 处理文件
            handle_file(item_path)
        elif os.path.isdir(item_path):
            # 处理子文件夹
            traverse_directory(item_path)

def handle_file(file_path):
    # 处理文件，这里可以记录文件路径和其他信息
    print("Processing file:", file_path)
    # 记录文件路径和时间戳等信息
    record_file_info(file_path)

def record_file_info(file_path):
    # 记录文件信息到文件中，可以保存为.txt或者.rds文件
    timestamp =  time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getctime(file_path)))
    with open("file_info.txt", "w") as f:
        f.write(f"{file_path}\t{timestamp}\n")

# 通过rds类来记录读入的rds文件的路径和时间戳
class rds:
    def __init__(self, path):
        self.check_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))
        self.dir = path
        self.file_path = []
    
    def expand(self):
        #在self.dir下找到所有的rds文件
        for root, dirs, files in os.walk(self.dir):
            for file in files:
                if file.endswith(".rds"):
                    self.file_path.append(os.path.join(root, file))

    def keep_newest(self):
        #找到最新的rds文件 每一组rds文件用.来分隔， 最后一个是后缀名，倒数第二个是六位数字的时间戳，再往前的部分则是样品名字。合并样品名字，保留最新的，保存它的绝对地址
        #例如 SNI.240322.rds 和 SNI.240321.rds 保留 SNI.240322.rds
        sample = {}
        self.before = len(self.file_path)
        for file in self.file_path:
            print(file)
            sample_name_pool = [parts for parts in file.split('/')[-1].split('.') if parts!="rds" and not (parts.isdigit() and len(parts)==6)]
            sample_name = ".".join(sample_name_pool)
            print("\t",sample_name)
            # print(sample_name)
            #时间为文件生成的时间
            time_stamp = os.path.getctime(file)
            if sample_name in sample:
                if time_stamp > sample[sample_name]:
                    sample[sample_name] = time_stamp
            else:
                sample[sample_name] = time_stamp
        for file in self.file_path:
            sample_name_pool = [parts for parts in file.split('/')[-1].split('.') if parts!="rds" and not (parts.isdigit() and len(parts)==6)]
            sample_name = ".".join(sample_name_pool)
            time_stamp = os.path.getctime(file)
            if time_stamp != sample[sample_name]:
                self.file_path.remove(file)
        self.after = len(self.file_path)
        print(f"Before: {self.before}, After: {self.after}")

    def output(self):
        #将最新的rds文件保存到一个新的文件夹里
        new_dir = self.dir + '_newest'
        #如果已经有这个文件夹，那就不用创建
        if not os.path.exists(self.dir + '_newest'):
            os.makedirs(new_dir)
        for file in self.file_path:
            os.system(f"cp {file} {new_dir}")

    def automate(self):
        self.expand()
        self.keep_newest()
        self.output()
        


class TreeNode:
    def __init__(self, name):
        self.name = name
        self.children = []

def build_tree_from_file(file_path):
    root = TreeNode("/")
    with open(file_path, 'r') as file:
        for line in file:
            path = line.strip()  # 去除每行末尾的换行符和空格
            add_to_tree(root, path.split('/'))
    return root

def add_to_tree(node, path_parts):
    if len(path_parts) == 1:
        node.children.append(TreeNode(path_parts[0]))
    else:
        for child in node.children:
            if child.name == path_parts[0]:
                add_to_tree(child, path_parts[1:])
                return
        new_node = TreeNode(path_parts[0])
        node.children.append(new_node)
        add_to_tree(new_node, path_parts[1:])

def print_tree(node, indent='', last=True):
    print(indent, end='')
    if last:
        print('└─', end='')
        indent += '  '
    else:
        print('├─', end='')
        indent += '│ '
    print(node.name)
    for i, child in enumerate(node.children):
        print_tree(child, indent, i == len(node.children) - 1)





if __name__ == "__main__":
    # 设置根文件夹路径
    root_directory = "/home/zhanglab02/2_filterred"
    # 遍历根文件夹
    traverse_directory(root_directory)
    # 示例：从文件中读取目录信息并构建树状结构
    file_path = '/home/zhanglab02/scripts/snRNA-seq-scripts/file_info.txt'
    root = build_tree_from_file(file_path)
    print_tree(root)
    #自动搜索所有的rds文件，并且将最新的rds文件保存到一个新的文件夹里
    search = rds(root_directory)
    search.automate()


