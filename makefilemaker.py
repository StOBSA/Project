from functools import reduce
from itertools import chain
import os, pathlib, re

repetitions = 30

targets = []
sizing = []
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
            name = f"{affix}Instance{finds[0]}" if not "Sizing" in path else f"{affix}Sizing{finds[0]}"
            print(f"{name}: {reduce(lambda x,y: x+y, ins,'')}")
            targets.append(f"{name} ")
targets.sort(key=lambda s:int(re.findall(r"(\d+)",s)[0]))
print("experiments: " + reduce(lambda x,y: x+y, targets,""))