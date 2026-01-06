import math
import sys


def mean(lst):
    return sum(lst) / len(lst)


def std(lst):
    m = mean(lst)
    var = mean([(el - m) ** 2 for el in lst])
    return math.sqrt(var)


with open(sys.argv[1]) as f:
    lines = f.readlines()

start_list = []
for i, line in enumerate(lines):
    if '   Step' in line:
        start_list.append(i + 1)

if not start_list:
    raise ValueError("Begin of MD trajectory not found")

stop_list = []
for i, line in enumerate(lines):
    if 'Loop time' in line:
        stop_list.append(i + 1)

if not stop_list:
    raise ValueError("End of MD trajectory not found")
print("start_list=", start_list)
print("stop_list=", stop_list)

i = start_list[-1]
j = stop_list[-1] - 1
data_lines = lines[i:j]

data_lines = list(map(lambda s: s.strip().split(), data_lines))

data = [[float(e) for e in line] for line in data_lines]
try:
    n_atoms = int(lines[j].strip().split()[-2])
    print("Number of atoms: ", n_atoms)
except Exception as e:
    print("Error: ", e)
en_data = [l[1] for l in data]
std_enrgy_deviation = std(en_data) / abs(mean(en_data))

print("Std. of NVE energy deviation: ", std_enrgy_deviation)
print("Std. of NVE energy deviation per atom: ", std_enrgy_deviation / n_atoms)
if std_enrgy_deviation > 4e-07:
    # print("ERROR: Std. of NVE energy deviation is too large.")
    raise ValueError("Std. of NVE energy deviation is too large.")
    sys.exit(1)
