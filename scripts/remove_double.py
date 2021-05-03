import itertools

a = ''.join(i for i, _ in itertools.groupby("modeling"))
print(list(itertools.groupby("modeling")))