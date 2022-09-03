from functools import reduce
from itertools import chain
import os, pathlib, re

repetitions = 3

targets = []
pattern = re.compile(r"terminals(\d+).csv")
counter = 1
for (path, folders, files) in chain(os.walk("SolidObstacles"),os.walk("SoftObstacles")):
    for file in files:
        finds = pattern.findall(file)
        if finds != []:
            for run in range(1,repetitions+1):
                print(f"make.target.{counter}:\n\tcargo run --release -- {path/pathlib.Path('terminals'+finds[0])}.csv {path/pathlib.Path('obstacles'+finds[0])}.csv {run} > {pathlib.Path('experiments')/path/pathlib.Path(f'Instance{finds[0]}Run{run}.csv')}")
                counter += 1
            ins = [f"make.target.{counter-i} " for i in range(1, repetitions+1)]
            affix = "Solid" if "Solid" in path else "Soft"
            print(f"{affix}Instance{finds[0]}: {reduce(lambda x,y: x+y, ins,'')}")
            targets.append(f"{affix}Instance{finds[0]} ")
targets.sort(key=lambda s:int(re.findall(r"(\d+)",s)[0]))
print("experiments: " + reduce(lambda x,y: x+y, targets,""))