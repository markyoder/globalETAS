import sys
import os
import glob

#
#print('mypath: {}'.format(__package__))
#
#print('path?: .{}'.format(os.path.split(sys.modules[__name__].__file__)[0]))
#print('__file__?'.format(__file__))
local_path = os.path.split(sys.modules[__name__].__file__)[0]
excluded_py_files=['__init__.py']
#
if not local_path in sys.path:
	sys.path.append(local_path)
#
__all__ = [os.path.splitext(os.path.split(g)[-1])[0] for g in glob.glob(os.path.join(local_path,'*.py')) if not os.path.split(g)[-1] in excluded_py_files]


