import sys
import math





litst = list(range(10))
print(litst)
print([str(a) + ("b" if a < 5 else "") for a in litst])





# print(sys.version)
# print(sys.executable)

# myl = list(range(21))
# print(myl)

# num_samps = len(myl)

# print()
# print()

# j = 0
# print(num_samps//2**(j+1))
# print(num_samps//2**(j))
# if (num_samps//2**(j)) % 2 == 0:
# 	print("even")
# else:
# 	print("not even")
# myl2 = [myl[i: i + 2**(j+2) if (num_samps <= i + 1.5*(2**(j+1)) and (num_samps//2**(j)) % 2 ) else i + 2**(j+1) : 2**j] for i in range(0, num_samps, 2**(j+1)) if i + 0.5*(2**(j+1)) < num_samps]
# print(myl2)
# myl2 = [myl[i: i + 2**(j+2) if num_samps <= i + 1.5*(2**(j+1)) else i + 2**(j+1) : 2**j] for i in range(0, num_samps, 2**(j+1)) if i + 0.5*(2**(j+1)) < num_samps]
# print(myl2)
# print()

# j = 1
# print(num_samps//2**(j+1))
# print(num_samps//2**(j))
# if (num_samps//2**(j)) % 2 == 0:
# 	print("even")
# else:
# 	print("not even")
# myl2 = [myl[i: i + 2**(j+2) if (num_samps <= i + 1.5*(2**(j+1)) and (num_samps//2**(j)) % 2 ) else i + 2**(j+1) : 2**j] for i in range(0, num_samps, 2**(j+1)) if i + 0.5*(2**(j+1)) < num_samps]
# print(myl2)
# myl2 = [myl[i: i + 2**(j+2) if num_samps <= i + 1.5*(2**(j+1)) else i + 2**(j+1) : 2**j] for i in range(0, num_samps, 2**(j+1)) if i + 0.5*(2**(j+1)) < num_samps]
# print(myl2)
# print()

# j = 2
# print(num_samps//2**(j+1))
# print(num_samps//2**(j))
# if (num_samps//2**(j)) % 2 == 0:
# 	print("even")
# else:
# 	print("not even")
# myl2 = [myl[i: i + 2**(j+2) if (num_samps <= i + 1.5*(2**(j+1)) and (num_samps//2**(j)) % 2 ) else i + 2**(j+1) : 2**j] for i in range(0, num_samps, 2**(j+1)) if i + 0.5*(2**(j+1)) < num_samps]
# print(myl2)
# myl2 = [myl[i: i + 2**(j+2) if num_samps <= i + 1.5*(2**(j+1)) else i + 2**(j+1) : 2**j] for i in range(0, num_samps, 2**(j+1)) if i + 0.5*(2**(j+1)) < num_samps]
# print(myl2)
# print()

# j = 3
# print(num_samps//2**(j+1))
# print(num_samps//2**(j))
# if (num_samps//2**(j)) % 2 == 0:
# 	print("even")
# else:
# 	print("not even")
# myl2 = [myl[i: i + 2**(j+2) if (num_samps <= i + 1.5*(2**(j+1)) and (num_samps//2**(j)) % 2 ) else i + 2**(j+1) : 2**j] for i in range(0, num_samps, 2**(j+1)) if i + 0.5*(2**(j+1)) < num_samps]
# print(myl2)
# myl2 = [myl[i: i + 2**(j+2) if num_samps <= i + 1.5*(2**(j+1)) else i + 2**(j+1) : 2**j] for i in range(0, num_samps, 2**(j+1)) if i + 0.5*(2**(j+1)) < num_samps]
# print(myl2)
# print()