from functools import reduce
from itertools import chain
import os, pathlib, re

pattern = re.compile(r"terminals(\d+).csv")
counter = 1
for (path, folders, files) in chain(os.walk("SolidObstacles"),os.walk("SoftObstacles")):
    for file in files:
        finds = pattern.findall(file)
        if finds != []:
            for run in range(1,4):
                print(f"make.target.{counter}:\n\tcargo run --release -- {path/pathlib.Path('terminals'+finds[0])}.csv {path/pathlib.Path('obstacles'+finds[0])}.csv > {pathlib.Path('experiments')/path/pathlib.Path(f'Instance{finds[0]}Run{run}.csv')}")
                counter += 1
all = [f"make.target.{counter} " for counter in range(1, counter)]
print("experiments: " + reduce(lambda x,y: x+y, all,""))