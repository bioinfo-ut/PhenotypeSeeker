mylits = [i for i in range(10000)]
print(mylits)

print([(mylits[i:i+1024], int(i/1024)) for i in range(0, len(mylits), 1024)])

# print(mylits[9990:10500])