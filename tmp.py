import multiprocessing as mpp
#import bindex
import earthquake_site.py
#

class First(object):
  def __init__(self):
    super(First, self).__init__()
    print("first")

class Second(object):
  def __init__(self):
    super(Second, self).__init__()
    print("second")

class Third(First, Second, mpp.Process):
  def __init__(self):
    super(Third, self).__init__()
    print("that's it")
