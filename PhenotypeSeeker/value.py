#!/usr/bin/python

from multiprocess import Value, Array, Manager
from multiprocessing.sharedctypes import RawArray
import os
import time
from ctypes import c_char, c_char_p

if __name__ == '__main__':

	l = Manager().list()
	n = Manager().Namespace()

	protsess = os.fork()
	if protsess == 0:
		print("child")
		mystr = "h".encode()
		l.append("mystriiiing")
		n.bla = "boo"
		print(l)
		print(id(l))
	else:
		print("paarent")
		time.sleep(2)
		print(l)
		print(id(l))
		print(n.bla)